# H5AD study input format

Cellenium uses AnnData (https://anndata.readthedocs.io/) files for study data import. Its `study_import.py`
command line tool has assumptions / requirements on the contents of h5ad files. The jupyter notebooks in
the `data_import/public_data` directory start with publicly available h5ad files and make modifications
to arrive at the right output structure.

We have developed `h5ad_preparation.py` as a small library of reusable functions for aligning h5ad files
for cellenium usage, and for verifying that requirements are met before the actual data import is attempted. So
many of the details below are automatically covered if you use `h5ad_preparation.py`.

The h5ad requirements are outlined as follows.

## Gene annotation dataframe (`var`)

Only the dataframe index is parsed, the actual columns are ignored.

Index content can be Ensembl Gene IDs, HGNC Symbols or Entrez Gene IDs. Human, mouse and rat genes
are supported currently.

## Cell annotation dataframe (`obs`)

All selected annotation groups (the columns in `obs`) are imported
(see `main_sample_attributes` and `secondary_sample_attributes` below).

Cellenium is designed for categorical cell annotations only (e.g. each distinct annotation value gets a color
assigned in the UI).

## Expression values

The matrix `X` must contain the primary data layer. As scRNA data is typically sparse, this should be a sparse
matrix, i.e. 0 values are left out.

Further data layers can optionally be placed in the `layer` object. `raw` is ignored by cellenium.

## Study metadata

The dictionary `uns['cellenium']` has multiple keys to provide metadata about the study:

required:

* `taxonomy_id`: NCBI ID for species, e.g. 9606 for human
* `title`: study display title
* `ncit_tissue_id`: NCIT ontology tissue IDs (list)
* `mesh_disease_ids`: MeSH disease IDs (list), include `HEALTHY` to identify that the study contains healthy subjects
* `X_pseudolayer_name`: describing the data found in the anonymous `X` matrix, e.g. `normalized`
* `main_sample_attributes`: list of cell annotation column names to include in the cellenium import

optional:

* `description`: study description text
* `secondary_sample_attributes`: further list of cell annotation column names to include in the cellenium import,
  but less prominent in the UI

## UMAP

A UMAP has to be calculated and stored in `obsm` (the `scanpy.tl.umap` function will do this).

## density (based on UMAP) subsampled cells

For studies with many cells (over 50k), a subset of cells is plotted in UMAP plots. The subset is
calculated upfront and stored as a list of indexes in `uns['cellenium']['umap_density_sampled_indices']`.


## Differentially expressed genes

For selected annotation groups, the differential expression of genes for each "cluster" (i.e. each
annotation value vs. all other values in the annotation group, e.g. one specific cell type vs. all
others) is stored in `uns['cellenium']['differentially_expressed_genes']` as a dataframe.

