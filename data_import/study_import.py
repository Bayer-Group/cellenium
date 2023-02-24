import argparse
import json
import logging
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from scanpy.pl._tools.scatterplots import _get_palette
import scipy.sparse as sparse
import tqdm
from anndata import AnnData
from psycopg2.extras import Json
from sqlalchemy import text
import h5ad_open

from postgres_utils import engine, import_df, NumpyEncoder

logging.basicConfig(format='%(asctime)s.%(msecs)03d %(process)d %(levelname)s %(name)s:%(lineno)d %(message)s',
                    datefmt='%Y%m%d-%H%M%S', level=logging.INFO)


def import_study_omics(study_id: int, adata: AnnData):
    logging.info('importing gene definitions of study')
    omics_df = pd.read_sql(
        "select omics_id, ensembl_gene_id, entrez_gene_ids, hgnc_symbols from omics_all where tax_id=%(tax_id)s and omics_type='gene'",
        engine,
        params={'tax_id': int(adata.uns['cellenium']['taxonomy_id'])},
        index_col='omics_id')
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
    adata_genes_df.drop_duplicates('omics_id', inplace=True)
    adata_genes_df['study_id'] = study_id
    import_df(adata_genes_df[['h5ad_var_index', 'omics_id', 'study_id']], 'study_omics')
    return adata_genes_df[['h5ad_var_index', 'h5ad_var_key', 'omics_id']]


def import_projection(adata, adata_samples_df, study_id, key):
    projection_df = pd.DataFrame({
        'study_id': study_id,
        'study_sample_id': adata_samples_df.study_sample_id,
        'projection_type': key,
        'projection': adata.obsm[f'X_{key}'][:, 0:2].tolist()
    })
    if f'{key}_density_sampled_indices' in adata.uns['cellenium']:
        projection_df['display_subsampling'] = False
        projection_df.loc[adata.uns['cellenium'][f'{key}_density_sampled_indices'], 'display_subsampling'] = True
    else:
        projection_df['display_subsampling'] = True
    import_df(projection_df, 'study_sample_projection')


def _projection_list(adata: AnnData):
    return adata.uns['cellenium'].get('import_projections', np.array(['umap'])).tolist()


def import_study_sample(study_id: int, adata: AnnData):
    logging.info('importing sample definitions')
    adata_samples_df = adata.obs.copy()
    adata_samples_df = adata_samples_df.reset_index(names='h5ad_obs_key')
    adata_samples_df = adata_samples_df.reset_index(names='h5ad_obs_index')
    adata_samples_df['study_sample_id'] = range(1, len(adata_samples_df) + 1)
    adata_samples_df = adata_samples_df[['study_sample_id', 'h5ad_obs_index']]
    adata_samples_df['study_id'] = study_id
    import_df(adata_samples_df, 'study_sample')
    with engine.connect() as connection:
        connection.execute(text("UPDATE study SET cell_count=:cell_count WHERE study_id=:study_id"), {
            'study_id': study_id,
            'cell_count': len(adata_samples_df)
        })
    for projection in _projection_list(adata):
        import_projection(adata, adata_samples_df, study_id, projection)

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
    secondary_sample_attributes = []
    if adata.uns['cellenium'].get('secondary_sample_attributes') is not None:
        secondary_sample_attributes = adata.uns['cellenium']['secondary_sample_attributes'].tolist()
        import_sample_annotations.extend(secondary_sample_attributes)

    with engine.connect() as connection:
        for annotation_col in import_sample_annotations:
            r = connection.execute(
                text("""SELECT annotation_group_id
                    FROM annotation_group WHERE h5ad_column=:h5ad_column"""), {
                    'h5ad_column': annotation_col
                }).fetchone()
            if r is None:
                r = connection.execute(text("""INSERT INTO annotation_group (h5ad_column, display_group)
                            VALUES (:h5ad_column, :h5ad_column_display)
                            RETURNING annotation_group_id"""), {
                    'h5ad_column': annotation_col,
                    'h5ad_column_display': annotation_col.replace('_', ' ')
                }).fetchone()
            annotation_group_id = r[0]

            connection.execute(text("""INSERT INTO study_annotation_group_ui (study_id, annotation_group_id, is_primary, ordering, differential_expression_calculated)
                                                                    VALUES (:study_id, :annotation_group_id, :is_primary, :ordering, False)"""),
                               {
                                   'study_id': study_id,
                                   'annotation_group_id': annotation_group_id,
                                   'is_primary': annotation_col not in secondary_sample_attributes,
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
                    connection.execute(text("""INSERT INTO annotation_value (annotation_group_id, h5ad_value, display_value)
                                            VALUES (:annotation_group_id, :h5ad_value, :h5ad_value_display)"""),
                                       {
                                           'annotation_group_id': annotation_group_id,
                                           'h5ad_value': value,
                                           'h5ad_value_display': value.replace('_', ' ')
                                       })

    annotation_definition_df = get_annotation_definition_df(import_sample_annotations)

    with engine.connect() as connection:
        adata_sample_annotations = adata.obs.copy()
        adata_sample_annotations = adata_sample_annotations.reset_index()
        adata_sample_annotations = adata_sample_annotations.merge(adata_samples_df,
                                                                  left_index=True, right_on='h5ad_obs_index')
        for h5ad_column in import_sample_annotations:
            palette = _get_palette(adata, h5ad_column)

            h5ad_one_annotation_df = adata_sample_annotations[[h5ad_column, 'study_sample_id']].copy()
            one_annotation_definition_df = annotation_definition_df[annotation_definition_df.h5ad_column == h5ad_column]
            annotation_df = h5ad_one_annotation_df.merge(one_annotation_definition_df,
                                                         left_on=h5ad_column, right_on='h5ad_value')
            annotation_df['color'] = annotation_df.apply(lambda row: palette[row.h5ad_value], axis=1)
            annotation_df = annotation_df[['study_sample_id', 'annotation_value_id', 'color']].copy()
            annotation_df['study_id'] = study_id
            annotation_df = annotation_df.groupby(['study_id', 'annotation_value_id', 'color'])[
                'study_sample_id'].apply(
                list).reset_index().rename(columns={'study_sample_id': 'study_sample_ids'})

            import_df(annotation_df, 'study_sample_annotation')


def import_study_layer_expression(study_id: int, layer_name: str, adata_genes_df, adata_samples_df, adata: AnnData):
    if layer_name is None:
        layer_name = adata.uns['cellenium']['X_pseudolayer_name']
        X = adata.X
    else:
        X = adata.layers[layer_name]
    logging.info(f'importing expression matrix {layer_name}')

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
        connection.connection.commit()

        sparse_X = sparse.csc_matrix(X)

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


def import_study(filename: str, analyze_database: bool) -> int:
    adata = h5ad_open.h5ad_read(filename)
    stored_filename = filename
    if stored_filename.startswith('scratch'):
        # filename inside scratch (scratch will be /h5ad_store in postgres docker)
        stored_filename = Path(filename).relative_to("scratch").as_posix()

    def _config_optional_list(key: str):
        if adata.uns['cellenium'].get(key) is not None:
            return adata.uns['cellenium'][key].tolist()
        return None

    with engine.connect() as connection:
        r = connection.execute(text("""INSERT INTO study (filename, study_name, description, tissue_ncit_ids, disease_mesh_ids, organism_tax_id,
               projections, reader_permissions, admin_permissions, legacy_config)
            VALUES (:filename, :study_name, :description, :tissue_ncit_ids, :disease_mesh_ids, :organism_tax_id,
               :projections, :reader_permissions, :admin_permissions, :legacy_config
            )
            RETURNING study_id"""), {
            'filename': stored_filename,
            'study_name': adata.uns['cellenium']['title'],
            'description': adata.uns['cellenium']['description'],
            'tissue_ncit_ids': adata.uns['cellenium']['ncit_tissue_ids'].tolist(),
            'disease_mesh_ids': adata.uns['cellenium']['mesh_disease_ids'].tolist(),
            'organism_tax_id': adata.uns['cellenium']['taxonomy_id'],
            'projections': _projection_list(adata),
            'reader_permissions': _config_optional_list('initial_reader_permissions'),
            'admin_permissions': _config_optional_list('initial_admin_permissions'),
            'legacy_config': Json(adata.uns['cellenium'].get('legacy_config'),
                                  dumps=lambda data: json.dumps(data, cls=NumpyEncoder))
        })
        study_id = r.fetchone()[0]
        logging.info("importing %s as study_id %s", filename, study_id)

        adata_genes_df = import_study_omics(study_id, adata)
        adata_samples_df = import_study_sample(study_id, adata)
        import_study_sample_annotation(study_id, adata_samples_df, adata)
        import_differential_expression(study_id, adata_genes_df, adata)

        import_study_layer_expression(study_id, None, adata_genes_df, adata_samples_df, adata)
        for layer_name in adata.layers.keys():
            import_study_layer_expression(study_id, layer_name, adata_genes_df, adata_samples_df, adata)

        connection.execute(text("UPDATE study SET visible=True WHERE study_id=:study_id"), {'study_id': study_id})
        logging.info("updating postgres statistics...")
        if analyze_database:
            connection.execute(text("call _analyze_schema()"))
        return study_id


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="cellenium study import tool")
    parser.add_argument('filename', help='h5ad file created for cellenium (e.g. from public_data/*.ipynb notebooks)',
                        type=str)
    parser.add_argument('--analyze-database', action='store_true')
    args = parser.parse_args()
    import_study(args.filename, args.analyze_database)
    logging.info('done')
