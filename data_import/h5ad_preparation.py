from typing import List
import tqdm
import scanpy as sc
import pandas as pd
import numpy as np
import os
import scipy
from anndata import AnnData
import logging
from pathlib import Path


logging.basicConfig(format='%(asctime)s.%(msecs)03d %(process)d %(levelname)s %(name)s:%(lineno)d %(message)s',
                    datefmt='%Y%m%d-%H%M%S', level=logging.INFO)

# default to the .gitignore-d scratch directory in this repo
basedir = Path(__file__).parent.parent.joinpath('scratch').resolve()
logging.info('local study files stored in: %s', basedir)


# downloads an H5AD file and returns the file handle
def get_h5ad_from_url(url: str, filename: str) -> AnnData:
    localfile = basedir.joinpath(f"{filename}.h5ad")
    adata = sc.read(localfile, backup_url=url)
    return adata

# downloads a dataset from the sfaira zoo and returns file handle, see here for a list of available datasets: https://theislab.github.io/sfaira-portal/Datasets
def get_sfaira_h5ad(sfaira_id: str) -> AnnData:
    # sfaira depends on tensorflow, we haven't included it in environment.yml, install when needed
    import sfaira
    datadir = basedir.joinpath('sfaira/data/')
    metadir = basedir.joinpath('sfaira/meta/')
    cachedir = basedir.joinpath('sfaira/cache/')

    ds = sfaira.data.Universe(data_path=datadir, meta_path=metadir, cache_path=cachedir)
    ds.subset(key="id", values=[sfaira_id])
    ds.download(verbose=1)
    ds.load(verbose=1)
    ds.datasets[sfaira_id].streamline_metadata(schema="sfaira")  # convert the metadata annotation to the sfaira standard
    adata = ds.datasets[sfaira_id].adata  # get the anndata object
    return adata

# checks if the matrix at .X is sparse, if not make it so. Stores matrix in specified layer
def make_sparse(adata: AnnData, layer = "counts"):
    if not scipy.sparse.issparse(adata.X):
        adata.layers[layer] = scipy.sparse.csr_matrix(adata.X)
    else:
        adata.layers[layer] = adata.X
    return adata

# check if an integer for a layer in anndata but not simply using the type
def isinteger(adata: AnnData, layer: str) -> bool:
    expression_vals = scipy.sparse.find(adata.layers[layer])[2]
    return all(np.equal(np.mod(expression_vals, 1), 0))

# make sure final data is in log space with normalized expression, layer to check is layer_from, layer to normalize to is layer_to
def make_norm_expression(adata: AnnData, layer_from = "counts", layer_to = "norm_log_expression"):
    if isinteger(adata, layer_from):
        adata.layers[layer_to] = adata.layers[layer_from]
        adata = sc.pp.normalize_total(adata, target_sum=1e4, layer=layer_to, inplace=True)
        sc.pp.log1p(adata, layer=layer_to)
    else:
        if np.max(scipy.sparse.find(adata.layers[layer_from])[2]) > 1000:
            adata.layers[layer_to] = adata.layers[layer_from]
            sc.pp.log1p(adata, layer=layer_to)
        else:
            adata.layers[layer_to] = adata.layers[layer_from]
    return adata

# basic umap calculation, doesn't take batch effect into account
def calculate_umap(adata: AnnData, layer = "norm_log_expression"):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, layer=layer, inplace=True)
    pca = sc.pp.scale(adata.layers[layer][:,adata.var.highly_variable])
    adata.obsm["X_pca"] = sc.pp.pca(pca)
    sc.pp.neighbors(adata, n_pcs=20, copy=False)
    adata = sc.tl.umap(adata, copy=True)
    return adata

# check if there is a umap or calculate one here        
def add_umap(adata: AnnData, layer = "norm_log_expression"):
    if adata.obsm is None:
        adata = calculate_umap(adata)
    else:
        if not "X_umap" in adata.obsm:
            adata = calculate_umap(adata)
    return adata
     
# add cellenium stuff to hormonize metadata
def add_cellenium_settings(adata: AnnData, main_attributes: List[str]):
    cellenium_settings(adata, main_sample_attributes=main_attributes)
    # TODO add NCIT, MeSH, tax_id, title, description, pubmed_id/link
    return adata
    
# add differential expression table to the anndata object
def add_differential_expression_tables(adata: AnnData, attributes: List[str], layer: str):
    diff_exp = calculate_differentially_expressed_genes(adata, attributes, layer)
    d = adata.uns.get('cellenium', {})
    adata.uns['cellenium'] = d
    d['differentially_expressed_genes'] = diff_exp
    return adata


def cellenium_settings(
        adata: AnnData,
        main_sample_attributes: List[str]
):
    d = adata.uns.get('cellenium', {})
    adata.uns['cellenium'] = d

    assert isinstance(main_sample_attributes, list)
    for a in main_sample_attributes:
        if a not in adata.obs.columns:
            raise Exception(f"main_sample_attributes: {a} not in observations dataframe")
    d['main_sample_attributes'] = main_sample_attributes


def calculate_differentially_expressed_genes(
        adata: AnnData,
        diffexp_attributes: List[str],
        layer: str,
        ngenes=100,
        diff_exp_min_group_expr=0.1,
        diff_exp_min_group_fc=0.5,
        diff_exp_max_notgroup_expr=1
):
    result_dataframes = []
    for diffexp_attribute in tqdm.tqdm(diffexp_attributes, desc='diff.exp. genes'):
        valid_attribute_group_check = (adata.obs[diffexp_attribute].value_counts() > 1)
        attr_values = valid_attribute_group_check.index[valid_attribute_group_check].tolist()

        sc.tl.rank_genes_groups(adata, diffexp_attribute, groups=attr_values, method='wilcoxon', layer=layer,
                                use_raw=False,
                                n_genes=ngenes)
        sc.tl.filter_rank_genes_groups(adata, min_in_group_fraction=diff_exp_min_group_expr,
                                       min_fold_change=diff_exp_min_group_fc,
                                       max_out_group_fraction=diff_exp_max_notgroup_expr, use_raw=False,
                                       key="rank_genes_groups", key_added="rank_genes_groups_filtered")
        for attr_value in attr_values:
            df = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups_filtered", group=attr_value)
            # remove filtered elements
            df = df[~df["names"].isnull()]
            df['ref_attr_value'] = attr_value
            df['cmp_attr_value'] = '_OTHERS_'
            df['attribute_name'] = diffexp_attribute
            result_dataframes.append(df)
    adata.uns.pop('rank_genes_groups', None)
    adata.uns.pop('rank_genes_groups_filtered', None)
    result_dataframe = pd.concat(result_dataframes, axis=0).reset_index(drop=True)
    logging.info("calculate_differentially_expressed_genes: found a list of genes for these attributes: %s",
                 result_dataframe['attribute_name'].unique().tolist())

    d = adata.uns.get('cellenium', {})
    adata.uns['cellenium'] = d
    d['differentially_expressed_genes'] = result_dataframe.copy()
    return result_dataframe
