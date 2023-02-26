import argparse
import io
import json
import logging
from pathlib import Path
from typing import List, Dict

import mudata
import numpy as np
import pandas as pd
import scanpy as sc
import sys
import tqdm
from anndata import AnnData
from muon import MuData
from psycopg2.extras import Json
from scanpy.pl._tools.scatterplots import _get_palette
from sqlalchemy import text

from postgres_utils import engine, import_df, NumpyEncoder, list_to_pgarray

logging.basicConfig(format='%(asctime)s.%(msecs)03d %(process)d %(levelname)s %(name)s:%(lineno)d %(message)s',
                    datefmt='%Y%m%d-%H%M%S', level=logging.INFO)


def get_all_regions(tax_id):
    existing_region_df = pd.read_sql(
        "select omics_id, display_symbol as region from omics_base where tax_id=%(tax_id)s and omics_type='region'",
        engine,
        params={'tax_id': tax_id},
        index_col='omics_id')

    cols = existing_region_df.columns.tolist()
    tmp = existing_region_df.reset_index()
    match_existing_region_df = tmp.melt(value_vars=cols, id_vars=['omics_id'], value_name='match_id')[
        ['omics_id', 'match_id']] \
        .drop_duplicates() \
        .set_index('match_id')
    return match_existing_region_df


def get_all_regions_ids_already_in(tax_id):
    existing_region_ids_df = pd.read_sql(
        "select region_id from omics_region_gene",
        engine)
    return existing_region_ids_df.region_id.tolist()


def import_region_base_and_study(study_id: int, data: AnnData, metadata: Dict):
    # for the moment just put all in as the genomic ranges changes from study to study
    # future we could intersect with transcription factor binding annotation and use
    # those coordinates
    logging.info('importing genomic ranges')

    match_existing_region_df = get_all_regions(metadata['taxonomy_id'])

    omics_df = pd.read_sql(
        "select omics_id, ensembl_gene_id, entrez_gene_ids, hgnc_symbols from omics_all where tax_id=%(tax_id)s and omics_type='gene'",
        engine,
        params={'tax_id': int(metadata['taxonomy_id'])},
        index_col='omics_id')
    cols = omics_df.columns.tolist()
    tmp = omics_df.explode('entrez_gene_ids').explode('hgnc_symbols').reset_index()
    match_df = tmp.melt(value_vars=cols, id_vars=['omics_id'], value_name='match_id')[['omics_id', 'match_id']] \
        .drop_duplicates() \
        .set_index('match_id')

    df_region = data.uns['atac']['peak_annotation'].reset_index()[['peak', 'gene_name']].rename(
        columns={'peak': 'region'})

    df_region = df_region.merge(match_df, left_on='gene_name', right_index=True, how='left')
    df_region.omics_id = df_region.omics_id.astype('Int64')

    df_region[['chromosome', 'start_position', 'end_position']] = df_region.region.str.extract('(.+):(\d+)-(\d+)',
                                                                                               expand=True)
    df_region['omics_type'] = 'region'
    df_region['tax_id'] = int(metadata['taxonomy_id'])
    df_region[['display_name']] = df_region[['region']]
    df_region[['display_symbol']] = df_region[['region']]

    # get h5ad_var_index
    h5ad_index = data.var.reset_index(names='h5ad_var_key').reset_index(names='h5ad_var_index')[
        ['h5ad_var_index', 'h5ad_var_key']]

    # insert into base
    omics_base_insert = df_region[['omics_type', 'tax_id', 'display_name', 'display_symbol']].drop_duplicates()
    # but not the ones which are already inserted
    regions_not_in_df = omics_base_insert.loc[~omics_base_insert.display_name.isin(match_existing_region_df.index), :]
    import_df(regions_not_in_df, 'omics_base')

    match_existing_region_df = get_all_regions(int(metadata['taxonomy_id']))
    omics_base_insert = omics_base_insert.merge(match_existing_region_df, left_on='display_name', right_index=True)

    # insert into omics_region
    omics_region_insert = df_region[
        ['chromosome', 'start_position', 'end_position', 'region']].drop_duplicates().reset_index(drop=True) \
        .merge(omics_base_insert[['display_name', 'omics_id']], left_on='region', right_on='display_name') \
        .rename(columns={'omics_id': 'region_id'}).drop_duplicates().drop('display_name', axis=1)
    # but not the ones which are already inserted
    omics_regions_not_in_df = omics_region_insert.loc[omics_region_insert.region.isin(regions_not_in_df.display_name),
                              :]
    import_df(omics_regions_not_in_df, 'omics_region')

    # insert into omics_region_gene
    omics_region_gene_insert = df_region[['omics_id', 'region']].dropna().drop_duplicates().reset_index(drop=True) \
        .rename(columns={'omics_id': 'gene_id'}) \
        .merge(omics_base_insert.set_index('display_symbol'), left_on='region', right_index=True)[
        ['gene_id', 'omics_id']] \
        .rename(columns={'omics_id': 'region_id'}) \
        .drop_duplicates()
    # but not the ones in
    all_region_ids = get_all_regions_ids_already_in(metadata['taxonomy_id'])
    omics_region_gene_not_in_df = omics_region_gene_insert.loc[~omics_region_gene_insert.region_id.isin(all_region_ids),
                                  :]
    import_df(omics_region_gene_not_in_df, 'omics_region_gene')

    # insert into study_omics
    study_omics_insert = \
        omics_base_insert.merge(h5ad_index.set_index('h5ad_var_key'), left_on='display_name', right_index=True)[
            ['omics_id', 'h5ad_var_index', 'display_name']] \
            .drop_duplicates().rename(columns={'display_name': 'h5ad_var_key'})
    study_omics_insert['study_id'] = study_id
    import_df(study_omics_insert.drop('h5ad_var_key', axis=1), 'study_omics')
    return study_omics_insert[['h5ad_var_index', 'h5ad_var_key', 'omics_id']]


def import_study_omics_genes(study_id: int, data: AnnData, metadata: Dict):
    logging.info('importing gene definitions of study')
    omics_df = pd.read_sql(
        "select omics_id, ensembl_gene_id, entrez_gene_ids, hgnc_symbols from omics_all where tax_id=%(tax_id)s and omics_type='gene'",
        engine,
        params={'tax_id': int(metadata['taxonomy_id'])},
        index_col='omics_id')

    # generate the mapping gene identifier to omics_id from database
    # could also be done more pandas-like
    #
    # cols = omics_df.columns.tolist()
    # tmp = omics_df.explode('entrez_gene_ids').explode('hgnc_symbols').reset_index()
    # match_df = tmp.melt(value_vars=cols, id_vars=['omics_id'], value_name='match_id')[['omics_id', 'match_id']] \
    #    .drop_duplicates() \
    #    .set_index('match_id')
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

    # now generate the dataframe to be inserted
    data_genes_df = data.var.copy()
    data_genes_df = data_genes_df.reset_index(names='h5ad_var_key')
    data_genes_df = data_genes_df.reset_index(names='h5ad_var_index')
    data_genes_df = data_genes_df.merge(match_df, how='inner', left_on='h5ad_var_key', right_index=True)
    data_genes_df.drop_duplicates('omics_id', inplace=True)
    data_genes_df['study_id'] = study_id
    import_df(data_genes_df[['h5ad_var_index', 'omics_id', 'study_id']], 'study_omics')
    return data_genes_df[['h5ad_var_index', 'h5ad_var_key', 'omics_id']]

def import_study_protein_antibody_tag(study_id: int, data: AnnData, metadata: Dict):
    logging.info('importing protein/antibody definitions of study')
    match_df = pd.read_sql(
        "select omics_id, display_symbol from omics_all where tax_id=%(tax_id)s and omics_type='protein_antibody_tag'",
        engine,
        params={'tax_id': int(metadata['taxonomy_id'])},
        index_col='display_symbol')

    # now generate the dataframe to be inserted
    data_genes_df = data.var.copy()
    data_genes_df = data_genes_df.reset_index(names='h5ad_var_key')
    data_genes_df = data_genes_df.reset_index(names='h5ad_var_index')
    data_genes_df = data_genes_df.merge(match_df, how='inner', left_on='h5ad_var_key', right_index=True)
    data_genes_df.drop_duplicates('omics_id', inplace=True)
    data_genes_df['study_id'] = study_id
    import_df(data_genes_df[['h5ad_var_index', 'omics_id', 'study_id']], 'study_omics')
    return data_genes_df[['h5ad_var_index', 'h5ad_var_key', 'omics_id']]


def import_projection(data, data_samples_df, study_id, key, modality=None):
    study_sample_ids = data.obs.merge(data_samples_df, left_index=True, right_on='h5ad_obs_key')[
        'study_sample_id'].tolist()
    projection_df = pd.DataFrame({
        'study_id': study_id,
        'study_sample_id': study_sample_ids,
        'modality': modality,
        'projection_type': key,
        'projection': data.obsm[f'X_{key}'][:, 0:2].tolist()
    })
    if f'{key}_density_sampled_indices' in data.uns['cellenium']:
        projection_df['display_subsampling'] = False
        projection_df.loc[data.uns['cellenium'][f'{key}_density_sampled_indices'], 'display_subsampling'] = True
    else:
        projection_df['display_subsampling'] = True
    import_df(projection_df, 'study_sample_projection')


def _projection_list(data: AnnData | MuData, filetype='h5ad'):
    if filetype == 'h5ad':
        return data.uns['cellenium'].get('import_projections', np.array(['umap'])).tolist()
    else:
        tmp = data.uns['cellenium'].get('import_projections')
        retlist = []
        for k in tmp.keys():
            retlist.extend([f'{k}:{proj}' for proj in tmp[k]])
        return retlist


def import_study_sample(study_id: int, data: AnnData | MuData, file_extension):
    logging.info('importing sample definitions')
    data_samples_df = data.obs.copy()
    data_samples_df = data_samples_df.reset_index(names='h5ad_obs_key')
    data_samples_df = data_samples_df.reset_index(names='h5ad_obs_index')
    data_samples_df['study_sample_id'] = range(1, len(data_samples_df) + 1)
    data_samples_df = data_samples_df[['study_sample_id', 'h5ad_obs_index', 'h5ad_obs_key']]
    data_samples_df['study_id'] = study_id
    import_df(data_samples_df[['study_sample_id', 'h5ad_obs_index', 'study_id']], 'study_sample')
    with engine.connect() as connection:
        connection.execute(text("UPDATE study SET cell_count=:cell_count WHERE study_id=:study_id"), {
            'study_id': study_id,
            'cell_count': len(data_samples_df)
        })
    if file_extension == 'h5ad':
        for projection in _projection_list(data, file_extension):
            import_projection(data, data_samples_df, study_id, projection)
    else:
        for modality, projections in data.uns['cellenium']['import_projections'].items():
            for projection in projections:
                import_projection(data.mod[modality], data_samples_df, study_id, projection, modality)

    return data_samples_df


def get_annotation_definition_df(h5ad_columns: List[str], modality=None):
    annotation_definition_df = pd.read_sql("""select a.annotation_group_id, a.h5ad_column, av.annotation_value_id, av.h5ad_value, a.modality
            from annotation_group a
            join annotation_value av on av.annotation_group_id = a.annotation_group_id
            where a.h5ad_column = any( %(h5ad_columns)s ) and modality = %(modality)s""", engine,
                                           params={'h5ad_columns': h5ad_columns,
                                                   'modality': modality if modality else ''})

    return annotation_definition_df


def import_study_sample_annotation(study_id: int, data_samples_df, data: AnnData | MuData, modality=None):
    logging.info('importing sample annotations')
    import_sample_annotations = data.uns['cellenium']['main_sample_attributes'].tolist()
    import_sample_annotations.extend(data.uns['cellenium'].get('advanced_sample_attributes', []))
    secondary_sample_attributes = []
    if data.uns['cellenium'].get('secondary_sample_attributes') is not None:
        secondary_sample_attributes = data.uns['cellenium']['secondary_sample_attributes'].tolist()
        import_sample_annotations.extend(secondary_sample_attributes)

    with engine.connect() as connection:
        for annotation_col in import_sample_annotations:
            annotation_col_clean = annotation_col.replace('_', ' ')

            r = connection.execute(
                text("""SELECT annotation_group_id
                    FROM annotation_group WHERE h5ad_column=:h5ad_column AND modality=:modality"""), {
                    'h5ad_column': annotation_col,
                    'modality': modality
                }).fetchone()
            if r is None:
                r = connection.execute(text("""INSERT INTO annotation_group (h5ad_column, display_group, modality)
                            VALUES (:h5ad_column, :h5ad_column_display, :modality)
                            RETURNING annotation_group_id"""), {
                    'h5ad_column': f'{annotation_col}',
                    'h5ad_column_display': annotation_col_clean,
                    'modality': modality

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

            values = data.obs[annotation_col].unique().tolist()
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

    annotation_definition_df = get_annotation_definition_df(import_sample_annotations, modality)

    with engine.connect() as connection:
        data_sample_annotations = data.obs.copy()
        data_sample_annotations = data_sample_annotations.reset_index()
        data_sample_annotations = data_sample_annotations.merge(data_samples_df,
                                                                left_index=True, right_on='h5ad_obs_index')
        for h5ad_column in import_sample_annotations:
            palette = _get_palette(data, h5ad_column)

            h5ad_one_annotation_df = data_sample_annotations[[h5ad_column, 'study_sample_id']].copy()
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


def import_study_layer_expression(study_id: int, layer_name: int, data_genes_df, data_samples_df,
                                  data: AnnData, metadata, omics_type, study_layer_id: int):
    if layer_name is None:
        layer_name = metadata['X_pseudolayer_name']
        X = data.X
    else:
        X = data.layers[layer_name]
    logging.info(f'importing expression matrix {layer_name} {omics_type}')

    with engine.connect() as connection:
        df_expr = generate_dense_expression_df_for_import(data, data_samples_df, data_genes_df, study_layer_id)
        df_expr = df_expr[['study_layer_id', 'omics_id', 'study_sample_ids', 'values']]

        # write df to string
        f = io.StringIO()
        df_expr.to_csv(f, index=False, header=False, sep="|")
        f.seek(0)

        # send to db
        cursor = connection.connection.cursor()
        logging.info('write expression to DB')
        cursor.copy_from(f, 'expression', columns=['study_layer_id', 'omics_id', 'study_sample_ids', 'values'], sep='|')
        connection.connection.commit()
        cursor.close()


'''
        sparse_X = sparse.csc_matrix(X)

        map_h5ad_var_index_to_omics_index = np.zeros(shape=[sparse_X.shape[1]], dtype=np.uint32)
        for i, row in data_genes_df.iterrows():
            map_h5ad_var_index_to_omics_index[row['h5ad_var_index']] = row['omics_id']

        map_h5ad_obs_index_to_studysample_index = np.zeros(shape=[sparse_X.shape[0]], dtype=np.uint32)
        for i, row in data_samples_df.iterrows():
            map_h5ad_obs_index_to_studysample_index[row['h5ad_obs_index']] = row['study_sample_id']

        for gene_i in tqdm.tqdm(range(0, sparse_X.shape[1]), desc=f'import expression matrix "{layer_name}"'):
            csc_gene_data = sparse.find(sparse_X.T[gene_i])
            data_cell_indexes = csc_gene_data[1]
            data_values = csc_gene_data[2]

            # omics_ids = map_h5ad_var_index_to_omics_index[data_gene_indexes]
            omics_id = map_h5ad_var_index_to_omics_index[gene_i]
            if omics_id > 0:
                studysample_ids = map_h5ad_obs_index_to_studysample_index[data_cell_indexes]
                connection.execute(text("""INSERT INTO expression (study_layer_id, omics_id, study_sample_ids, values)
                            VALUES (:study_layer_id, :omics_id, :study_sample_ids, :values)"""), {
                    'study_layer_id': study_layer_id,
                    'omics_id': omics_id,
                    'study_sample_ids': studysample_ids.tolist(),
                    'values': data_values.tolist()
                })
'''


def import_differential_expression(study_id: int, data_genes_df, data: AnnData | MuData, modality=None):
    if 'differentially_expressed_genes' not in data.uns['cellenium']:
        return
    logging.info('importing differentially expressed genes')
    df = data.uns['cellenium']['differentially_expressed_genes']
    df = df.merge(data_genes_df, left_on='names', right_on='h5ad_var_key')

    annotation_definition_df = get_annotation_definition_df(df['attribute_name'].unique().tolist(), modality)
    df = df.merge(annotation_definition_df, left_on=['attribute_name', 'ref_attr_value'],
                  right_on=['h5ad_column', 'h5ad_value'])
    df['study_id'] = study_id
    df.logfoldchanges = df.logfoldchanges.fillna(-1)  # TODO: this is as muon puts logfoldchange to NaN
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


def generate_dense_expression_df_for_import(adata: AnnData, samples_df: pd.DataFrame, data_genes_df: pd.DataFrame,
                                            study_layer_id: int):
    # replace columns with omics_ids
    df = adata.to_df()
    df = df.T.join(data_genes_df.set_index('h5ad_var_key')[['omics_id']]).dropna()
    df.omics_id = df.omics_id.astype(int)
    df = df.set_index('omics_id').T
    # replace the index with study_sample_ids
    df = df.join(samples_df.set_index('h5ad_obs_key')[['study_sample_id']]).set_index('study_sample_id')
    # generate the to be imported dataframe
    collect = []
    for omics_id in tqdm.tqdm(df.columns):
        tmp = df.loc[(df.loc[:, omics_id] > 0), omics_id]
        study_sample_ids = tmp.index.astype(str).tolist()
        values = tmp.astype(str).tolist()
        if len(study_sample_ids) > 0:
            collect.append(
                {"study_sample_ids": list_to_pgarray(study_sample_ids), "values": list_to_pgarray(values),
                 'omics_id': omics_id})
    df_ret = pd.DataFrame(collect)
    df_ret['study_layer_id'] = study_layer_id
    return df_ret


def import_study(filename: str, analyze_database: bool) -> int:
    """
    TODO enable S3 access as we had before:

    adata = h5ad_open.h5ad_read(filename)
    stored_filename = filename
    if stored_filename.startswith('scratch'):
        # filename inside scratch (scratch will be /h5ad_store in postgres docker)
        stored_filename = Path(filename).relative_to("scratch").as_posix()
    """

    file_extension = Path(filename).suffix
    file_extension = file_extension[1:] if file_extension.startswith('.') else file_extension
    if file_extension == 'h5ad':
        data = sc.read_h5ad(filename)
    else:
        data = mudata.read_h5mu(filename)

    def _config_optional_list(key: str):
        if data.uns['cellenium'].get(key) is not None:
            return data.uns['cellenium'][key].tolist()
        return None

    def _generate_study_layer(study_id, layer_name, omics_type=None):
        with engine.connect() as connection:
            r = connection.execute(text("""INSERT INTO study_layer (study_id, layer, omics_type)
                                    VALUES (:study_id, :layer, :omics_type)
                                    RETURNING study_layer_id"""), {
                'study_id': study_id,
                'layer': layer_name,
                'omics_type': omics_type
            })
            study_layer_id = r.fetchone()[0]

            connection.execute(text("call add_studylayer_partition(:study_layer_id)"),
                               {'study_layer_id': study_layer_id})
            connection.connection.commit()

        return study_layer_id

    with engine.connect() as connection:
        r = connection.execute(text("""INSERT INTO study (filename, study_name, description, tissue_ncit_ids, disease_mesh_ids, organism_tax_id,
               projections, reader_permissions, admin_permissions, legacy_config)
            VALUES (:filename, :study_name, :description, :tissue_ncit_ids, :disease_mesh_ids, :organism_tax_id,
               :projections, :reader_permissions, :admin_permissions, :legacy_config
            )
            RETURNING study_id"""), {
            # TODO 'filename': stored_filename,
            'filename': Path(filename).relative_to("scratch").as_posix(),
            # filename inside scratch (scratch will be /h5ad_store in postgres docker)
            'study_name': data.uns['cellenium']['title'],
            'description': data.uns['cellenium']['description'],
            'tissue_ncit_ids': data.uns['cellenium']['ncit_tissue_ids'].tolist(),
            'disease_mesh_ids': data.uns['cellenium']['mesh_disease_ids'].tolist(),
            'organism_tax_id': data.uns['cellenium']['taxonomy_id'],
            'projections': _projection_list(data, file_extension),
            'reader_permissions': _config_optional_list('initial_reader_permissions'),
            'admin_permissions': _config_optional_list('initial_admin_permissions'),
            'legacy_config': Json(data.uns['cellenium'].get('legacy_config'),
                                  dumps=lambda data: json.dumps(data, cls=NumpyEncoder))
        })
        study_id = r.fetchone()[0]
        logging.info("importing %s as study_id %s", filename, study_id)
        if file_extension == 'h5mu':
            modalities = data.uns['cellenium']['modalities']
        else:
            modalities = {'rna': 'gene'}

        data_samples_df = import_study_sample(study_id, data, file_extension)
        for modality in modalities.keys():
            if file_extension == 'h5mu':
                cur_data = data.mod[modality]
            else:
                cur_data = data
            import_study_sample_annotation(study_id, data_samples_df, cur_data, modality)

        for modality in modalities.items():
            data_type = modality[1]  # the data_type
            if file_extension == 'h5mu':
                cur_data = data.mod[modality[0]]
            else:
                cur_data = data
            meta_data = data.uns['cellenium']
            if (data_type == 'gene'):
                data_genes_df = import_study_omics_genes(study_id, cur_data, meta_data)
                import_differential_expression(study_id, data_genes_df, cur_data, modality[0])
            elif (data_type == 'region'):
                data_region_df = import_region_base_and_study(study_id, cur_data,
                                                              meta_data)  # since regions from study to study change we import those which are not yet in the database
                import_differential_expression(study_id, data_region_df, cur_data, modality[0])
            elif (data_type == 'protein_antibody_tag'):
                data_protein_df = import_study_protein_antibody_tag(study_id, cur_data, meta_data)
                import_differential_expression(study_id, data_protein_df, cur_data, modality[0])


        study_layer_id = _generate_study_layer(study_id, 'umap', 'gene')
        for modality in modalities.items():
            omics_type = modality[1]  # the data_type
            if file_extension == 'h5mu':
                cur_data = data.mod[modality[0]]
                cur_data_samples_df = data_samples_df
            else:
                cur_data = data
                cur_data_samples_df = data_samples_df.loc[data_samples_df.h5ad_obs_key.isin(cur_data.obs.index), :]

            meta_data = data.uns['cellenium']

            if omics_type == 'gene':
                import_study_layer_expression(study_id, None, data_genes_df, cur_data_samples_df, cur_data, meta_data,
                                              omics_type, study_layer_id)
            # TODO: layers, e.g. imputed matrices not supported yet
            #                for layer_name in cur_data.layers.keys():
            #                    import_study_layer_expression(study_id, layer_name, data_genes_df, cur_data_samples_df, cur_data,
            #                                                  meta_data, omics_type)
            if omics_type == 'region':
                import_study_layer_expression(study_id, None, data_region_df, cur_data_samples_df, cur_data, meta_data,
                                              omics_type, study_layer_id)
            # TODO: layers not supported yet
            #                for layer_name in cur_data.layers.keys():
            #                    import_study_layer_expression(study_id, layer_name, data_region_df, data_samples_df, cur_data,
            #                                                  meta_data, omics_type)
            if omics_type == 'protein_antibody_tag':
                import_study_layer_expression(study_id, None, data_protein_df, cur_data_samples_df, cur_data, meta_data,
                                              omics_type, study_layer_id)

        connection.execute(text("UPDATE study SET visible=True WHERE study_id=:study_id"), {'study_id': study_id})
        logging.info("updating postgres statistics...")
        if analyze_database:
            connection.execute(text("call _analyze_schema()"))
        return study_id


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="cellenium study import tool")
    parser.add_argument('filename', help='h5ad/h5mu file created for cellenium (e.g. using a jupyter lab notebook).',
                        type=str)
    parser.add_argument('--analyze-database', help='analyses the database schema after insert of study',
                        action='store_true')
    args = parser.parse_args()
    import_study(args.filename, args.analyze_database)
    logging.info('done')
