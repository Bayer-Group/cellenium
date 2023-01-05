from typing import List

import pandas as pd
from anndata import AnnData
from psycopg2 import extras
from sqlalchemy import text
import numpy as np
import scanpy as sc
import scipy.sparse as sparse
import tqdm
import logging

from huge_palette import huge_palette
from postgres_utils import engine, import_df

logging.basicConfig(format='%(asctime)s.%(msecs)03d %(process)d %(levelname)s %(name)s:%(lineno)d %(message)s',
                    datefmt='%Y%m%d-%H%M%S', level=logging.INFO)


def import_study_omics(study_id: int, adata: AnnData):
    logging.info('importing gene definitions of study')
    omics_df = pd.read_sql(
        "select omics_id, ensembl_gene_id, entrez_gene_ids, hgnc_symbols from omics_all where tax_id=9606 and omics_type='gene'",
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
    adata_genes_df = adata_genes_df.reset_index(names='h5ad_var_key')
    adata_genes_df = adata_genes_df.reset_index(names='h5ad_var_index')
    adata_genes_df = adata_genes_df.merge(match_df, how='inner', left_on='h5ad_var_key', right_index=True)
    adata_genes_df['study_id'] = study_id
    import_df(adata_genes_df[['h5ad_var_index', 'omics_id', 'study_id']], 'study_omics')
    return adata_genes_df[['h5ad_var_index', 'h5ad_var_key', 'omics_id']]


def import_study_sample(study_id: int, adata: AnnData):
    logging.info('importing sample definitions')
    adata_samples_df = adata.obs.copy()
    adata_samples_df = adata_samples_df.reset_index(names='h5ad_obs_key')
    adata_samples_df = adata_samples_df.reset_index(names='h5ad_obs_index')
    adata_samples_df['study_sample_id'] = range(1, len(adata_samples_df) + 1)
    adata_samples_df = adata_samples_df[['study_sample_id', 'h5ad_obs_index']]
    adata_samples_df['study_id'] = study_id
    import_df(adata_samples_df, 'study_sample')

    if 'X_umap' in adata.obsm:
        projection_df = pd.DataFrame({
            'study_id': study_id,
            'study_sample_id': adata_samples_df.study_sample_id,
            'projection_type': 'umap',
            'projection': adata.obsm['X_umap'].tolist(),
            'display_subsampling': True
        })
        import_df(projection_df, 'study_sample_projection')

    return adata_samples_df


def get_annotation_definition_df(h5ad_columns: List[str]):
    annotation_definition_df = pd.read_sql("""select a.annotation_group_id, a.h5ad_column, av.annotation_value_id, av.h5ad_value
            from annotation_group a
            join annotation_value av on av.annotation_group_id = a.annotation_group_id
            where a.h5ad_column = any( %(h5ad_columns)s )""", engine,
                                           params={'h5ad_columns': h5ad_columns})
    return annotation_definition_df


def import_study_sample_annotation(study_id: int, adata_samples_df, adata: AnnData):
    logging.info('importing sample annotations')
    import_sample_annotations = adata.uns['cellenium']['main_sample_attributes'].tolist()
    import_sample_annotations.extend(adata.uns['cellenium'].get('advanced_sample_attributes', []))

    with engine.connect() as connection:
        for annotation_col in import_sample_annotations:
            r = connection.execute(
                text("""SELECT annotation_group_id,
                    (select count(1) from annotation_value v where v.annotation_group_id=annotation_group.annotation_group_id) cnt_values
                    FROM annotation_group WHERE h5ad_column=:h5ad_column"""), {
                    'h5ad_column': annotation_col
                }).fetchone()
            if r is None:
                r = connection.execute(text("""INSERT INTO annotation_group (h5ad_column, display_group)
                            VALUES (:h5ad_column, :h5ad_column)
                            RETURNING annotation_group_id, 0"""), {
                    'h5ad_column': annotation_col
                }).fetchone()
            annotation_group_id = r[0]
            color_index = r[1]

            connection.execute(text("""INSERT INTO study_annotation_group_ui (study_id, annotation_group_id, is_primary, ordering, differential_expression_calculated)
                                                                    VALUES (:study_id, :annotation_group_id, True, :ordering, False)"""),
                               {
                                   'study_id': study_id,
                                   'annotation_group_id': annotation_group_id,
                                   'ordering': import_sample_annotations.index(annotation_col)
                               })

            values = adata.obs[annotation_col].unique().tolist()
            for value in values:
                r = connection.execute(text(
                    "SELECT annotation_value_id FROM annotation_value WHERE annotation_group_id=:annotation_group_id AND h5ad_value=:h5ad_value"),
                    {
                        'annotation_group_id': annotation_group_id,
                        'h5ad_value': value
                    }).fetchone()
                if r is None:
                    connection.execute(text("""INSERT INTO annotation_value (annotation_group_id, h5ad_value, display_value, color)
                                            VALUES (:annotation_group_id, :h5ad_value, :h5ad_value, :color)"""), {
                        'annotation_group_id': annotation_group_id,
                        'h5ad_value': value,
                        'color': huge_palette[color_index]
                    })
                    color_index += 1

    annotation_definition_df = get_annotation_definition_df(import_sample_annotations)

    with engine.connect() as connection:
        adata_sample_annotations = adata.obs.copy()
        adata_sample_annotations = adata_sample_annotations.reset_index()
        adata_sample_annotations = adata_sample_annotations.merge(adata_samples_df,
                                                                  left_index=True, right_on='h5ad_obs_index')
        for h5ad_column in import_sample_annotations:
            h5ad_one_annotation_df = adata_sample_annotations[[h5ad_column, 'study_sample_id']].copy()
            one_annotation_definition_df = annotation_definition_df[annotation_definition_df.h5ad_column == h5ad_column]
            annotation_df = h5ad_one_annotation_df.merge(one_annotation_definition_df,
                                                         left_on=h5ad_column, right_on='h5ad_value')
            annotation_df = annotation_df[['study_sample_id', 'annotation_value_id']].copy()
            annotation_df['study_id'] = study_id
            annotation_df = annotation_df.groupby(['study_id', 'annotation_value_id'])['study_sample_id'].apply(
                list).reset_index().rename(columns={'study_sample_id': 'study_sample_ids'})
            import_df(annotation_df, 'study_sample_annotation')


def import_study_layer_expression(study_id: int, layer_name: str, adata_genes_df, adata_samples_df, adata: AnnData):
    logging.info('importing expression matrix')
    with engine.connect() as connection:
        r = connection.execute(text("""INSERT INTO study_layer (study_id, layer, omics_type)
                                VALUES (:study_id, :layer, 'gene')
                                RETURNING study_layer_id"""), {
            'study_id': study_id,
            'layer': layer_name
        })
        study_layer_id = r.fetchone()[0]

        connection.execute(text("call add_studylayer_partition(:study_layer_id)"),
                           {'study_layer_id': study_layer_id})

        sparse_X = sparse.csc_matrix(adata.layers[layer_name])

        map_h5ad_var_index_to_omics_index = np.zeros(shape=[sparse_X.shape[1]], dtype=np.uint32)
        for i, row in adata_genes_df.iterrows():
            map_h5ad_var_index_to_omics_index[row['h5ad_var_index']] = row['omics_id']
        map_h5ad_obs_index_to_studysample_index = np.zeros(shape=[sparse_X.shape[0]], dtype=np.uint32)
        for i, row in adata_samples_df.iterrows():
            map_h5ad_obs_index_to_studysample_index[row['h5ad_obs_index']] = row['study_sample_id']

        for gene_i in tqdm.tqdm(range(0, sparse_X.shape[1]), desc=f'import expression matrix "{layer_name}"'):
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


def import_differential_expression(study_id: int, adata_genes_df, adata: AnnData):
    if 'differentially_expressed_genes' not in adata.uns['cellenium']:
        return
    logging.info('importing differentially expressed genes')
    df = adata.uns['cellenium']['differentially_expressed_genes']
    df = df.merge(adata_genes_df, left_on='names', right_on='h5ad_var_key')
    annotation_definition_df = get_annotation_definition_df(df['attribute_name'].unique().tolist())
    df = df.merge(annotation_definition_df, left_on=['attribute_name', 'ref_attr_value'],
                  right_on=['h5ad_column', 'h5ad_value'])
    df['study_id'] = study_id
    df.rename(
        columns={'pvals': 'pvalue', 'pvals_adj': 'pvalue_adj', 'scores': 'score', 'logfoldchanges': 'log2_foldchange'},
        inplace=True)
    import_df(df[['study_id', 'omics_id', 'annotation_value_id', 'pvalue', 'pvalue_adj', 'score', 'log2_foldchange']],
              'differential_expression')
    with engine.connect() as connection:
        connection.execute(text("""UPDATE study_annotation_group_ui SET differential_expression_calculated=True
                                    WHERE study_id = :study_id and annotation_group_id = any (:annotation_group_ids)"""),
                           {
                               'study_id': study_id,
                               'annotation_group_ids': df['annotation_group_id'].unique().tolist()
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
        import_study_sample_annotation(study_id, adata_samples_df, adata)
        import_differential_expression(study_id, adata_genes_df, adata)

        for layer_name in adata.layers.keys():
            import_study_layer_expression(study_id, layer_name, adata_genes_df, adata_samples_df, adata)

        return study_id


if __name__ == "__main__":
    adata = sc.read_h5ad('../scratch/pancreas_atlas_subset.h5ad')
    import_study('pancreas_atlas_subset', adata)
    # adata = sc.read_h5ad('../scratch/heart_failure_reichart2022_subset.h5ad')
    # import_study('heart_failure_reichart2022_subset', adata)
    # adata = sc.read_h5ad('../scratch/heart_failure_reichart2022.h5ad')
    # import_study('heart_failure_reichart2022', adata)
