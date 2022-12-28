import scanpy as sc
import numpy as np
import h5ad_preparation as prep

def add_cellenium_settings(adata):
    prep.cellenium_settings(
        adata,
        main_sample_attributes=['celltype']
    )
    

def download_prepare_h5ad():
    # pancreas data from here: https://www.nature.com/articles/s41592-021-01336-8
    url = "https://figshare.com/ndownloader/files/24539828"
    adata = sc.read("../scratch/pancreas_original.h5ad", backup_url=url)
    
    # TODO make sure this is a sparse matrix. We're creating a 2.5GB file out of 316MB input.
    # TODO log scale the data (anyway, 'counts' data doesn't look like counts...)    
    
    add_cellenium_settings(adata)
    prep.calculate_differentially_expressed_genes(adata, ['celltype'], 'counts')
    adata.write('../scratch/pancreas.h5ad')

    query = np.array([s in ["celseq"] for s in adata.obs.tech])
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=200,
        batch_key="tech",
        subset=True
    )
    adata = adata[query].copy()
    adata = adata[:, adata.var_names].copy()
    add_cellenium_settings(adata)
    prep.calculate_differentially_expressed_genes(adata, ['celltype'], 'counts')
    adata.write('../scratch/pancreas_subset.h5ad')
