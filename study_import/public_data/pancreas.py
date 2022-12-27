import scanpy as sc
import numpy as np


def download_prepare_h5ad():
    # pancreas data from here: https://www.nature.com/articles/s41592-021-01336-8
    url = "https://figshare.com/ndownloader/files/24539828"
    adata = sc.read("pancreas.h5ad", backup_url=url)
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
    adata.write('../scratch/pancreas_subset.h5ad')
