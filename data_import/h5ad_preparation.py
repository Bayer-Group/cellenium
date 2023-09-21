import logging
import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import cello
import Density_Sampling.density_sampling as density_sampling
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import tqdm
from anndata import AnnData
from muon import MuData
from muon import atac as ac

logging.basicConfig(
    format="%(asctime)s.%(msecs)03d %(process)d %(levelname)s %(name)s:%(lineno)d %(message)s",
    datefmt="%Y%m%d-%H%M%S",
    level=logging.INFO,
)

# default to the .gitignore-d scratch directory in this repo
basedir = Path(__file__).parent.parent.joinpath("scratch").resolve()
logging.info("local study files stored in: %s", basedir)


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

    datadir = basedir.joinpath("sfaira/data/")
    metadir = basedir.joinpath("sfaira/meta/")
    cachedir = basedir.joinpath("sfaira/cache/")

    ds = sfaira.data.Universe(data_path=datadir, meta_path=metadir, cache_path=cachedir)
    ds.subset(key="id", values=[sfaira_id])
    ds.download(verbose=1)
    ds.load(verbose=1)
    ds.datasets[sfaira_id].streamline_metadata(schema="sfaira")  # convert the metadata annotation to the sfaira standard
    adata = ds.datasets[sfaira_id].adata  # get the anndata object
    return adata


def jupyter_h5ad_overview(adata: AnnData):
    from IPython.display import HTML, display, display_pretty

    pd.set_option("display.max_columns", 100)

    def _header(h):
        display(HTML(f"<h2>{h}</h2>"))

    def _df(title, df):
        _header(title)
        display(HTML(df._repr_html_()))

    def _matrix(title, m):
        _header(title)
        display(m.shape)
        display_pretty(m)

    _df("obs", adata.obs)
    _df("var", adata.var)
    _matrix("X", adata.X)
    if adata.raw is not None:
        _matrix("raw", adata.raw)
    for layer in adata.layers:
        _matrix(f'layer "{layer}"', adata.layers[layer])
    _header("uns")
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


def validate_gene_ids(adata: AnnData, taxonomy_id: int):
    # cellenium recognizes HGNC Symbols and Ensembl Gene IDs in the annotation index.
    # Simple validation which checks some housekeeping genes.
    require_ids = {
        9541: ["ATF1", "ENSMFAG00000000109", "AK1", "ENSMFAG00000032155", "REL", "ENSMFAG00000038696"],
        9606: [
            "ATF1",
            "ENSG00000123268",
            "ADA",
            "ENSG00000196839",
            "AK1",
            "ENSG00000106992",
            "REL",
            "ENSG00000162924",
            "SLA",
            "ENSG00000155926",
        ],
        10090: [
            "Atf1",
            "ENSMUSG00000023027",
            "Ada",
            "ENSMUSG00000017697",
            "Ak1",
            "ENSMUSG00000026817",
            "Rel",
            "ENSMUSG00000020275",
            "Sla",
            "ENSMUSG00000022372",
        ],
        10116: ["Atf1", "ENSRNOG00000061088", "Rel", "ENSRNOG00000054437", "Sla", "ENSRNOG00000056714"],
    }
    ids = require_ids[taxonomy_id]
    for id in ids:
        if id in adata.var.index:
            return True
    raise AssertionError(f"None of {ids} where found in adata.var")


# checks if the matrix at .X (or layer) is sparse, if not make it so
def make_sparse(adata: AnnData, layer=None):
    if not scipy.sparse.issparse(_get_X_or_layer(adata, layer)):
        _set_X_or_layer(adata, layer, scipy.sparse.csr_matrix(_get_X_or_layer(adata, layer)))
        logging.info("make_sparse: conversion to sparse matrix done")


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
        logging.info("make_norm_expression: integer values detected - applied normalize_total and log")
    else:
        if np.max(scipy.sparse.find(_get_X_or_layer(adata, layer))[2]) > 1000:
            sc.pp.log1p(adata, layer=layer)
            logging.info("make_norm_expression: high values detected - applied log")
        else:
            logging.info("make_norm_expression: no transformations necessary")


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
        if "X_umap" not in adata.obsm:
            calculate_umap(adata)


def density_sample_umap(adata: AnnData, desired_samples=50000):
    if len(adata.obs) > desired_samples:
        sampled_indices = density_sampling.density_sampling(adata.obsm["X_umap"], metric="euclidean", desired_samples=desired_samples)
        _cellenium_uns_dictionary(adata)["umap_density_sampled_indices"] = sampled_indices


def _cellenium_uns_dictionary(adata: AnnData) -> dict:
    d = adata.uns.get("cellenium", {})
    adata.uns["cellenium"] = d
    return d


# add differential expression table to the anndata object
# TODO detect automatically which attributes to do differential expression for by fuzzy matching (but optional, I'd say...)
def add_differential_expression_tables(adata: AnnData, attributes: List[str], layer: str):
    diff_exp = calculate_differentially_expressed_genes(adata, attributes, layer)
    _cellenium_uns_dictionary(adata)["differentially_expressed_genes"] = diff_exp
    return adata


def set_cellenium_metadata(
    data: AnnData | MuData,
    title: str,
    description: str,
    taxonomy_id: int,
    ncit_tissue_ids: List[str],
    mesh_disease_ids: List[str],
    X_pseudolayer_name: str,
    main_sample_attributes: List[str] | Dict[str, List[str]],
    secondary_sample_attributes: Optional[List[str]] = None,
    import_projections: Optional[List[str]] = None,
    initial_reader_permissions: Optional[List[str]] = None,
    initial_admin_permissions: Optional[List[str]] = None,
    modalities: Optional[List[Dict]] = None,
    pubmed_id: Optional[str] = None,
    geo_accession: Optional[str] = None,
    any_metadata: Optional[Dict[str, str | List[str]]] = None,
):
    if any_metadata is None:
        any_metadata = {}
    if import_projections is None:
        import_projections = ["umap"]
    if secondary_sample_attributes is None:
        secondary_sample_attributes = []

    def _check_cell_annotation(data, attribute: str):
        if attribute not in data.obs.columns:
            raise Exception(f"attribute {attribute} not in observations dataframe")
        count_unique_values = len(data.obs[a].unique())
        if count_unique_values > 100:
            raise Exception(f"attribute {attribute} has {count_unique_values} unique annotations, 100 is the maximum")

    d = _cellenium_uns_dictionary(data)

    assert isinstance(main_sample_attributes, (list, dict))

    if not modalities:
        for a in main_sample_attributes:
            _check_cell_annotation(data, a)
    else:
        for modality in modalities:
            for a in main_sample_attributes[modality]:
                _check_cell_annotation(data.mod[modality], a)

    d["main_sample_attributes"] = main_sample_attributes

    assert title is not None
    d["title"] = title
    d["description"] = description
    assert isinstance(taxonomy_id, int)
    d["taxonomy_id"] = taxonomy_id
    if ncit_tissue_ids is None:
        ncit_tissue_ids = []
    assert isinstance(ncit_tissue_ids, list)
    d["ncit_tissue_ids"] = ncit_tissue_ids
    if mesh_disease_ids is None:
        mesh_disease_ids = []
    assert isinstance(mesh_disease_ids, list)
    d["mesh_disease_ids"] = mesh_disease_ids
    assert X_pseudolayer_name is not None
    d["X_pseudolayer_name"] = X_pseudolayer_name

    for a in secondary_sample_attributes:
        _check_cell_annotation(data, a)
        if a in main_sample_attributes:
            raise Exception(f"secondary_sample_attributes: {a} is also listed in main_sample_attributes, overlap not allowed")
    d["secondary_sample_attributes"] = secondary_sample_attributes

    if not modalities:
        for p in import_projections:
            assert data.obsm[f"X_{p}"] is not None
        d["import_projections"] = import_projections
    else:
        collect = defaultdict(list)
        for modality in modalities:
            for p in import_projections:
                assert data.mod[modality].obsm[f"X_{p}"] is not None
                collect[modality].append(p)
        d["import_projections"] = dict(collect)

        for modality in modalities:
            if "cellenium" in data.mod[modality].uns:
                data.mod[modality].uns["cellenium"]["main_sample_attributes"] = d["main_sample_attributes"][modality]
            else:
                d.mod[modality].uns["cellenium"] = {"main_sample_attributes": d["main_sample_attributes"][modality]}

    d["initial_reader_permissions"] = initial_reader_permissions
    d["initial_admin_permissions"] = initial_admin_permissions
    d["modalities"] = modalities

    if pubmed_id:
        any_metadata["Pubmed ID"] = str(pubmed_id)
    if geo_accession:
        any_metadata["GEO Accession"] = str(geo_accession)
    d["metadata"] = any_metadata


# calculate differentially expressed genes using rank_genes_groups from scanpy
def calculate_differentially_expressed_genes(
    adata: AnnData,
    diffexp_attributes: List[str],
    ngenes=100,
    diff_exp_min_group_expr=0.1,
    diff_exp_min_group_fc=0.5,
    diff_exp_max_notgroup_expr=1,
):
    result_dataframes = []
    for diffexp_attribute in tqdm.tqdm(diffexp_attributes, desc="diff.exp. genes"):
        valid_attribute_group_check = adata.obs[diffexp_attribute].value_counts() > 1
        attr_values = valid_attribute_group_check.index[valid_attribute_group_check].tolist()

        sc.tl.rank_genes_groups(
            adata,
            diffexp_attribute,
            groups=attr_values,
            method="wilcoxon",
            use_raw=False,
            n_genes=ngenes,
        )
        sc.tl.filter_rank_genes_groups(
            adata,
            min_in_group_fraction=diff_exp_min_group_expr,
            min_fold_change=diff_exp_min_group_fc,
            max_out_group_fraction=diff_exp_max_notgroup_expr,
            use_raw=False,
            key="rank_genes_groups",
            key_added="rank_genes_groups_filtered",
        )
        for attr_value in attr_values:
            df = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups_filtered", group=attr_value)
            # remove filtered elements
            df = df[~df["names"].isnull()]
            df["ref_attr_value"] = attr_value
            df["cmp_attr_value"] = "_OTHERS_"
            df["attribute_name"] = diffexp_attribute
            result_dataframes.append(df)
    adata.uns.pop("rank_genes_groups", None)
    adata.uns.pop("rank_genes_groups_filtered", None)
    result_dataframe = pd.concat(result_dataframes, axis=0).reset_index(drop=True)
    logging.info(
        "calculate_differentially_expressed_genes: found a list of genes for these attributes: %s",
        result_dataframe["attribute_name"].unique().tolist(),
    )

    _cellenium_uns_dictionary(adata)["differentially_expressed_genes"] = result_dataframe.copy()
    return result_dataframe


def calculate_differential_peaks(adata: AnnData, diffexp_attributes: List[str], npeaks=100):
    result_dataframes = []
    for diffexp_attribute in tqdm.tqdm(diffexp_attributes, desc="differential peaks"):
        valid_attribute_group_check = adata.obs[diffexp_attribute].value_counts() > 1
        attr_values = valid_attribute_group_check.index[valid_attribute_group_check].tolist()

        ac.tl.rank_peaks_groups(
            adata,
            diffexp_attribute,
            groups=attr_values,
            method="wilcoxon",
            use_raw=False,
            n_genes=npeaks,
        )

        for attr_value in attr_values:
            df = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups", group=attr_value)
            # remove filtered elements
            df = df[~df["names"].isnull()]
            df["ref_attr_value"] = attr_value
            df["cmp_attr_value"] = "_OTHERS_"
            df["attribute_name"] = diffexp_attribute
            result_dataframes.append(df)
    adata.uns.pop("rank_genes_groups", None)
    result_dataframe = pd.concat(result_dataframes, axis=0).reset_index(drop=True)
    logging.info(
        "calculate_differentially_expressed_genes: found a list of genes for these attributes: %s",
        result_dataframe["attribute_name"].unique().tolist(),
    )

    _cellenium_uns_dictionary(adata)["differentially_expressed_genes"] = result_dataframe.copy()
    return result_dataframe


def cello_classify_celltypes(adata: AnnData, cello_clustering_attribute: str):
    if adata.uns["cellenium"]["taxonomy_id"] != 9606:
        logging.info("skipping CellO classification, taxonomy_id is not human")
        return
    resource_dir = basedir.joinpath("cello_resources")
    os.makedirs(resource_dir, exist_ok=True)
    # Mahmoud: CellO makes mistakes sometimes due to ribosomal protein genes, so would be good to filter them out before the CellO call
    remove_ribo = adata.var_names.str.startswith(("RPS", "RPL"))
    adata_cello = adata[:, ~remove_ribo].copy()
    cello.scanpy_cello(
        adata_cello,
        clust_key=cello_clustering_attribute,
        rsrc_loc=resource_dir,
        term_ids=True,
    )
    adata.obs["CellO_celltype"] = adata.obs.join(adata_cello.obs["Most specific cell type"])["Most specific cell type"]

    updated_sample_attributes = ["CellO_celltype"]
    updated_sample_attributes.extend(adata.uns["cellenium"]["main_sample_attributes"])
    adata.uns["cellenium"]["main_sample_attributes"] = updated_sample_attributes
