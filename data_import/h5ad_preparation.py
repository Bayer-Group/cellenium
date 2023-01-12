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
import Density_Sampling.density_sampling as density_sampling

logging.basicConfig(format='%(asctime)s.%(msecs)03d %(process)d %(levelname)s %(name)s:%(lineno)d %(message)s',
                    datefmt='%Y%m%d-%H%M%S', level=logging.INFO)

# default to the .gitignore-d scratch directory in this repo
basedir = Path(__file__).parent.parent.joinpath('scratch').resolve()
logging.info('local study files stored in: %s', basedir)


# downloads an H5AD file and returns the file handle
def get_h5ad_from_url(url: str, filename: str) -> AnnData:
    if not filename.endswith("original"):
        filename = filename + "_original"
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
    ds.datasets[sfaira_id].streamline_metadata(
        schema="sfaira")  # convert the metadata annotation to the sfaira standard
    adata = ds.datasets[sfaira_id].adata  # get the anndata object
    return adata


def jupyter_h5ad_overview(adata: AnnData):
    from IPython.display import display, HTML, display_pretty
    pd.set_option("display.max_columns", 100)

    def _header(h):
        display(HTML(f'<h2>{h}</h2>'))

    def _df(title, df):
        _header(title)
        display(HTML(df._repr_html_()))

    def _matrix(title, m):
        _header(title)
        display(m.shape)
        display_pretty(m)

    _df('obs', adata.obs)
    _df('var', adata.var)
    _matrix('X', adata.X)
    if adata.raw is not None:
        _matrix('raw', adata.raw)
    for layer in adata.layers:
        _matrix(f'layer "{layer}"', adata.layers[layer])
    _header('uns')
    display_pretty(adata.uns)


def remove_raw_and_layers(adata: AnnData):
    adata.raw = None
    for layer in list(adata.layers.keys()):
        adata.layers.pop(layer)


def _get_X_or_layer(adata: AnnData, layer):
    return adata.layers[layer] if layer else adata.X


def _set_X_or_layer(adata: AnnData, layer, m):
    if layer:
        adata.layers[layer] = m
    else:
        adata.X = m


# checks if the matrix at .X (or layer) is sparse, if not make it so
def make_sparse(adata: AnnData, layer=None):
    if not scipy.sparse.issparse(_get_X_or_layer(adata, layer)):
        _set_X_or_layer(adata, layer, scipy.sparse.csr_matrix(_get_X_or_layer(adata, layer)))
        logging.info('make_sparse: conversion to sparse matrix done')


def filter_outliers(adata: AnnData):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)


# check if an integer for a layer in anndata but not simply using the type
def isinteger(adata: AnnData, layer=None) -> bool:
    expression_vals = scipy.sparse.find(_get_X_or_layer(adata, layer))[2]
    return all(np.equal(np.mod(expression_vals, 1), 0))


# make sure final data is in log space with normalized expression, use .X (or layer)
def make_norm_expression(adata: AnnData, layer=None):
    if isinteger(adata, layer):
        sc.pp.normalize_total(adata, target_sum=1e4, layer=layer, inplace=True)
        sc.pp.log1p(adata, layer=layer)
        logging.info('make_norm_expression: integer values detected - applied normalize_total and log')
    else:
        if np.max(scipy.sparse.find(_get_X_or_layer(adata, layer))[2]) > 1000:
            sc.pp.log1p(adata, layer=layer)
            logging.info('make_norm_expression: high values detected - applied log')
        else:
            logging.info('make_norm_expression: no transformations necessary')


def adata_subset_for_testing(adata, cells_filter_attribute, cells_filter_values, n_top_genes):
    query = np.array([s in cells_filter_values for s in adata.obs[cells_filter_attribute]])
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=cells_filter_attribute, subset=True)
    adata_sub = adata[query].copy()
    adata_sub = adata_sub[:, adata.var_names].copy()
    return adata_sub


# basic umap calculation, doesn't take batch effect into account
def calculate_umap(adata: AnnData, layer=None):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, layer=layer, inplace=True)
    pca = sc.pp.scale(_get_X_or_layer(adata, layer)[:, adata.var.highly_variable])
    adata.obsm["X_pca"] = sc.pp.pca(pca)
    sc.pp.neighbors(adata, n_pcs=20, copy=False)
    sc.tl.umap(adata)


# check if there is a umap or calculate one here
def add_umap(adata: AnnData, layer=None):
    if adata.obsm is None:
        calculate_umap(adata, layer)
    else:
        if not "X_umap" in adata.obsm:
            calculate_umap(adata)


def density_sample_umap(adata: AnnData, desired_samples=50000):
    if len(adata.var) > desired_samples:
        sampled_indices = density_sampling.density_sampling(adata.obsm['X_umap'], metric='euclidean',
                                                            desired_samples=desired_samples)
        _cellenium_uns_dictionary(adata)['umap_density_sampled_indices'] = sampled_indices


# add cellenium stuff to harmonize metadata
def add_cellenium_settings(adata: AnnData, main_attributes: List[str]):
    cellenium_settings(adata, main_sample_attributes=main_attributes)
    # TODO add NCIT, MeSH, tax_id, title, description, pubmed_id/link
    return adata


def _cellenium_uns_dictionary(adata: AnnData) -> dict:
    d = adata.uns.get('cellenium', {})
    adata.uns['cellenium'] = d
    return d


# add differential expression table to the anndata object
# TODO detect automatically which attributes to do differential expression for by fuzzy matching (but optional, I'd say...)
def add_differential_expression_tables(adata: AnnData, attributes: List[str], layer: str):
    diff_exp = calculate_differentially_expressed_genes(adata, attributes, layer)
    _cellenium_uns_dictionary(adata)['differentially_expressed_genes'] = diff_exp
    return adata


def set_cellenium_metadata(
        adata: AnnData,
        title: str,
        description: str,
        taxonomy_id: int,
        ncit_tissue_ids: List[str],
        mesh_disease_ids: List[str],
        X_pseudolayer_name: str,
        main_sample_attributes: List[str]
):
    # lets keep this stable for the jupyter way of h5ad generation

    d = _cellenium_uns_dictionary(adata)

    assert isinstance(main_sample_attributes, list)
    for a in main_sample_attributes:
        if a not in adata.obs.columns:
            raise Exception(f"main_sample_attributes: {a} not in observations dataframe")
    d['main_sample_attributes'] = main_sample_attributes

    assert title is not None
    d['title'] = title
    d['description'] = description
    assert isinstance(taxonomy_id, int)
    d['taxonomy_id'] = taxonomy_id
    assert isinstance(ncit_tissue_ids, list)
    assert len(ncit_tissue_ids) > 0
    d['ncit_tissue_ids'] = ncit_tissue_ids
    assert isinstance(mesh_disease_ids, list)
    d['mesh_disease_ids'] = mesh_disease_ids
    assert X_pseudolayer_name is not None
    d['X_pseudolayer_name'] = X_pseudolayer_name


# cellenium meta data
def cellenium_settings(
        adata: AnnData,
        title: str,
        description: str,
        taxonomy_id: str,
        ncit_tissue_id: List[str],
        mesh_disease_id: List[str],
        pubmed_id: str,
        data_source: str
):
    d = _cellenium_uns_dictionary(adata)

    # assert isinstance(main_sample_attributes, list)
    # for a in main_sample_attributes:
    #    if a not in adata.obs.columns:
    #        raise Exception(f"main_sample_attributes: {a} not in observations dataframe")

    # d['main_sample_attributes'] = main_sample_attributes
    d['title'] = title
    d['description'] = description
    d['taxonomy_id'] = taxonomy_id
    d['ncit_tissue_id'] = ncit_tissue_id
    d['mesh_disease_id'] = mesh_disease_id
    d['pubmed_id'] = pubmed_id
    d['data_source'] = data_source


# add cellenium meta data to AnnData.uns,
def add_cellenium_settings(adata: AnnData,
                           title: str,
                           description: str,
                           taxonomy_id: str,
                           ncit_tissue_id: List[str],
                           mesh_disease_id: List[str],
                           pubmed_id: str,
                           data_source: str
                           ):
    ncit_tissue_id = [x.strip() for x in ncit_tissue_id.split(',')]
    mesh_disease_id = [x.strip() for x in mesh_disease_id.split(',')]
    cellenium_settings(adata, title, description, taxonomy_id, ncit_tissue_id, mesh_disease_id, pubmed_id, data_source)
    # TODO: add fuzzy matching to detect main cell type attributes
    # TODO: add fuzzy matching of cell types to cell ontology
    return adata


# calculate differentially expressed genes using rank_genes_groups from scanpy
def calculate_differentially_expressed_genes(
        adata: AnnData,
        diffexp_attributes: List[str],
        ngenes=100,
        diff_exp_min_group_expr=0.1,
        diff_exp_min_group_fc=0.5,
        diff_exp_max_notgroup_expr=1
):
    result_dataframes = []
    for diffexp_attribute in tqdm.tqdm(diffexp_attributes, desc='diff.exp. genes'):
        valid_attribute_group_check = (adata.obs[diffexp_attribute].value_counts() > 1)
        attr_values = valid_attribute_group_check.index[valid_attribute_group_check].tolist()

        sc.tl.rank_genes_groups(adata, diffexp_attribute, groups=attr_values, method='wilcoxon',
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

    _cellenium_uns_dictionary(adata)['differentially_expressed_genes'] = result_dataframe.copy()
    return result_dataframe
