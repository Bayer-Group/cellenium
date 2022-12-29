import scanpy as sc
import numpy as np
import h5ad_preparation as prep


#downloads an H5AD file and returns the file handle
def download_h5ad():
    # pancreas data from here: https://www.nature.com/articles/s41592-021-01336-8
    url = "https://figshare.com/ndownloader/files/24539828"
    adata = sc.read("../../scratch/pancreas_original.h5ad", backup_url=url)
    return adata

#add cellenium stuff to hormonize metadata
def add_cellenium_settings(adata):
    prep.cellenium_settings(
        adata,
        main_sample_attributes=['celltype']
        # TODO add NCIT, MeSH, tax_id, title, description, pubmed_id/link
    )
    
#check if an integer but not simply using the type
def isinteger(x):
    return np.equal(np.mod(x, 1), 0)    
        
#basic umap calculation, doesn't take batch effect into account
def calculate_umap(adata):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, layer="norm_log_expression", inplace=True)
    pca = sc.pp.scale(adata.layers["norm_log_expression"][:,adata.var.highly_variable])
    adata.obsm["X_pca"] = sc.pp.pca(pca)
    sc.pp.neighbors(adata, n_pcs=20, copy=False)
    adata = sc.tl.umap(adata, copy=True)
    return adata

    

#download file
adata = download_h5ad()


# TODO make sure this is a sparse matrix. We're creating a 2.5GB file out of 316MB input.
if not scipy.sparse.issparse(adata.X):
    adata.layers["counts"] = scipy.sparse.csr_matrix(adata.X)
else:
    adata.layers["counts"] = adata.X


# TODO log scale the data (anyway, 'counts' data doesn't look like counts...) 
expression_vals = scipy.sparse.find(adata.layers["counts"])[2]
if all(isinteger(expression_vals)):
    adata.layers["norm_log_expression"] = adata.layers["counts"]
    adata = sc.pp.normalize_total(adata, target_sum=1e4, layer="norm_log_expression", inplace=True)
    sc.pp.log1p(adata, layer="norm_log_expression")
else:
    if np.max(expression_vals) > 1000:
        adata.layers["norm_log_expression"] = adata.layers["counts"]
        sc.pp.log1p(adata, layer="norm_log_expression")
    else:
        adata.layers["norm_log_expression"] = adata.layers["counts"]


        
# TODO we need a study with UMAP projection or calculate one here        
if adata.obsm is None:
    adata = calculate_umap(adata)
else:
    if not "X_umap" in adata.obsm:
        adata = calculate_umap(adata)

        
# harmonize cellenium metadata        
add_cellenium_settings(adata)


# for testing, subsample the dataset
query = np.array([s in ["celseq"] for s in adata.obs.tech])
sc.pp.highly_variable_genes(adata, n_top_genes=200, batch_key="tech", subset=True)
adata = adata[query].copy()
adata = adata[:, adata.var_names].copy()


prep.calculate_differentially_expressed_genes(adata, ['celltype'], 'norm_log_expression')
adata.write('../scratch/pancreas_subset.h5ad')
