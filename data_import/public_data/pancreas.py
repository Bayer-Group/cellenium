import scanpy as sc
import numpy as np
import os
import scipy
#import h5ad_preparation as prep


# downloads an H5AD file and returns the file handle
def download_h5ad():
    # pancreas data from here: https://www.nature.com/articles/s41592-021-01336-8
    url = "https://figshare.com/ndownloader/files/24539828"
    adata = sc.read("../../scratch/pancreas_original.h5ad", backup_url=url)
    return adata

# downloads an H5AD file and returns the file handle
def get_h5ad_from_url(url, basedir, filename):
    datadir = basedir + "/" + filename + ".h5ad"
    adata = sc.read(datadir, backup_url=url)
    return adata

# downloads a dataset from the sfaira zoo and returns file handle, see here for a list of available datasets: https://theislab.github.io/sfaira-portal/Datasets
def get_sfaira_h5ad(sfaira_id, basedir):
    
    datadir = os.path.join(basedir, 'sfaira/data/')
    metadir = os.path.join(basedir, 'sfaira/meta/')
    cachedir = os.path.join(basedir, 'sfaira/cache/')

    ds = sfaira.data.Universe(data_path=datadir, meta_path=metadir, cache_path=cachedir)
    ds.subset(key="id", values=[sfaira_id])
    ds.download(verbose=1)
    ds.load(verbose=1)
    ds.datasets[sfaira_id].streamline_metadata(schema="sfaira")  # convert the metadata annotation to the sfaira standard
    adata = ds.datasets[sfaira_id].adata  # get the anndata object
    return adata

# check if an integer but not simply using the type
def isinteger(x):
    return np.equal(np.mod(x, 1), 0)    

# checks if the matrix at .X is sparse, if not make it so. Stores matrix in specified layer
def make_sparse(adata, layer = "counts"):
    if not scipy.sparse.issparse(adata.X):
        adata.layers[layer] = scipy.sparse.csr_matrix(adata.X)
    else:
        adata.layers[layer] = adata.X
    return adata

# make sure final data is in log space with normalized expression, layer to check is layer_from, layer to normalize to is layer_to
def make_norm_expression(adata, layer_from = "counts", layer_to = "norm_log_expression"):
    expression_vals = scipy.sparse.find(adata.layers[layer_from])[2]
    if all(isinteger(expression_vals)):
        adata.layers[layer_to] = adata.layers[layer_from]
        adata = sc.pp.normalize_total(adata, target_sum=1e4, layer=layer_to, inplace=True)
        sc.pp.log1p(adata, layer=layer_to)
    else:
        if np.max(expression_vals) > 1000:
            adata.layers[layer_to] = adata.layers[layer_from]
            sc.pp.log1p(adata, layer=layer_to)
        else:
            adata.layers[layer_to] = adata.layers[layer_from]
    return adata

# basic umap calculation, doesn't take batch effect into account
def calculate_umap(adata, layer = "norm_log_expression"):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, layer=layer, inplace=True)
    pca = sc.pp.scale(adata.layers[layer][:,adata.var.highly_variable])
    adata.obsm["X_pca"] = sc.pp.pca(pca)
    sc.pp.neighbors(adata, n_pcs=20, copy=False)
    adata = sc.tl.umap(adata, copy=True)
    return adata

# check if there is a umap or calculate one here        
def make_umap(adata):
    if adata.obsm is None:
        adata = calculate_umap(adata)
    else:
        if not "X_umap" in adata.obsm:
            adata = calculate_umap(adata)
    return adata

        
# add cellenium stuff to hormonize metadata
def add_cellenium_settings(adata):
    prep.cellenium_settings(
        adata,
        main_sample_attributes=['celltype']
        # TODO add NCIT, MeSH, tax_id, title, description, pubmed_id/link
    )
    
        


# pancreas url example
adata = get_h5ad_from_url("https://figshare.com/ndownloader/files/24539828", "cellenium/scratch", "pancreas_atlas")
adata = make_sparse(adata)
adata = make_norm_expression(adata)
adata = make_umap(adata)

# lung sfaira example
#basedir = ""
#sfaira_id = "homosapiens_lungparenchyma_2019_10x3v2_madissoon_001_10.1186/s13059-019-1906-x"
#adata = get_sfaira_h5ad(sfaira_id, "cellenium/scratch")
#adata = make_sparse(adata)
#adata = make_norm_expression(adata)
#adata = make_umap(adata)



# harmonize cellenium metadata        
add_cellenium_settings(adata)


# for testing, subsample the dataset
query = np.array([s in ["celseq"] for s in adata.obs.tech])
sc.pp.highly_variable_genes(adata, n_top_genes=200, batch_key="tech", subset=True)
adata = adata[query].copy()
adata = adata[:, adata.var_names].copy()


prep.calculate_differentially_expressed_genes(adata, ['celltype'], 'norm_log_expression')
adata.write('../scratch/pancreas_subset.h5ad')
