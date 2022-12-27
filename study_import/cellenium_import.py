import pandas as pd
from anndata import AnnData
from psycopg2 import extras
from sqlalchemy import text
import numpy as np
import scanpy as sc
import scipy.sparse as sparse

from cellenium.study_import.postgres_utils import engine, import_df


def import_study_omics(study_id: int, adata: AnnData):
    omics_df = pd.read_sql(
        "select omics_id, ensembl_gene_id, entrez_gene_ids, hgnc_symbols from omics where tax_id=9606 and omics_type='gene'",
        engine, index_col='omics_id')
    match_dfs = []
    for col in ['ensembl_gene_id', 'entrez_gene_ids', 'hgnc_symbols']:
        match_df = pd.DataFrame(omics_df[[col]])
        match_df.rename(columns={col: 'match_id'}, inplace=True)
        match_df['omics_id'] = match_df.index
        match_df = match_df.explode('match_id', ignore_index=True)
        match_df.drop_duplicates('match_id', inplace=True)
        match_df.set_index('match_id', inplace=True)
        match_dfs.append(match_df)
    match_df = pd.concat(match_dfs)

    adata_genes_df = adata.var.copy()
    adata_genes_df = adata_genes_df.reset_index(names='_var_key')
    adata_genes_df = adata_genes_df.reset_index(names='h5ad_var_index')
    adata_genes_df = adata_genes_df.merge(match_df, how='inner', left_on='_var_key', right_index=True)
    adata_genes_df = adata_genes_df[['h5ad_var_index', 'omics_id']]
    adata_genes_df['study_id'] = study_id
    import_df(adata_genes_df, 'study_omics')
    return adata_genes_df


def import_study_sample(study_id: int, adata: AnnData):
    adata_samples_df = adata.obs.copy()
    adata_samples_df = adata_samples_df.reset_index(names='h5ad_obs_key')
    adata_samples_df = adata_samples_df.reset_index(names='h5ad_obs_index')
    adata_samples_df['study_sample_id'] = range(1, len(adata_samples_df) + 1)
    adata_samples_df = adata_samples_df[['study_sample_id', 'h5ad_obs_index']]
    adata_samples_df['study_id'] = study_id
    adata_samples_df['display_subsampling'] = True
    adata_samples_df
    import_df(adata_samples_df, 'study_sample')
    return adata_samples_df


def import_study_layer_expression(study_id: int, layer_name: str, adata_genes_df, adata_samples_df, adata: AnnData):
    with engine.connect() as connection:
        r = connection.execute(text("""INSERT INTO study_layer (study_id, layer)
                                VALUES (:study_id, :layer)
                                RETURNING study_layer_id"""), {
            'study_id': study_id,
            'layer': layer_name
        })
        study_layer_id = r.fetchone()[0]

        connection.execute(text(f"""create table expression_{study_layer_id}
                    partition of expression
                    (
                        study_layer_id,
                        omics_id,
                        study_sample_ids,
                        values
                        )
                    for values in ( {study_layer_id} )"""))

        sparse_X = sparse.csc_matrix(adata.layers[layer_name])

        map_h5ad_var_index_to_omics_index = np.zeros(shape=[sparse_X.shape[1]], dtype=np.uint32)
        for i, row in adata_genes_df.iterrows():
            map_h5ad_var_index_to_omics_index[row['h5ad_var_index']] = row['omics_id']
        map_h5ad_obs_index_to_studysample_index = np.zeros(shape=[sparse_X.shape[0]], dtype=np.uint32)
        for i, row in adata_samples_df.iterrows():
            map_h5ad_obs_index_to_studysample_index[row['h5ad_obs_index']] = row['study_sample_id']

        for gene_i in range(0, sparse_X.shape[1]):
            csc_gene_data = sparse.find(sparse_X.T[gene_i])
            adata_cell_indexes = csc_gene_data[1]
            adata_values = csc_gene_data[2]

            # omics_ids = map_h5ad_var_index_to_omics_index[adata_gene_indexes]
            omics_id = map_h5ad_var_index_to_omics_index[gene_i]
            if omics_id > 0:
                studysample_ids = map_h5ad_obs_index_to_studysample_index[adata_cell_indexes]

                connection.execute(text("""INSERT INTO expression (study_layer_id, omics_id, study_sample_ids, values)
                            VALUES (:study_layer_id, :omics_id, :study_sample_ids, :values)"""), {
                    'study_layer_id': study_layer_id,
                    'omics_id': omics_id,
                    'study_sample_ids': studysample_ids.tolist(),
                    'values': adata_values.tolist()
                })


def import_study(study_name: str, adata: AnnData) -> int:
    with engine.connect() as connection:
        r = connection.execute(text("""INSERT INTO study (study_name)
            VALUES (:study_name)
            RETURNING study_id"""), {
            'study_name': study_name
        })
        study_id = r.fetchone()[0]

        adata_genes_df = import_study_omics(study_id, adata)
        adata_samples_df = import_study_sample(study_id, adata)

        for layer_name in adata.layers.keys():
            import_study_layer_expression(study_id, layer_name, adata_genes_df, adata_samples_df, adata)

        return study_id


if __name__ == "__main__":
    adata = sc.read_h5ad('../scratch/pancreas_subset.h5ad')
    import_study('pancreas', adata)