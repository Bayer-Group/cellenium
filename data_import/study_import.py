import argparse
import io
import json
import logging
import pathlib
from upath import UPath
from typing import List, Dict

import numpy as np
import pandas as pd
import scipy.sparse as sparse
import tqdm
from anndata import AnnData
from muon import MuData
from psycopg2.extras import Json
from scanpy.plotting._tools.scatterplots import _get_palette
from sqlalchemy import text
import os
import h5ad_open
from postgres_utils import (
    import_df,
    NumpyEncoder,
    get_aws_db_engine,
    get_local_db_engine,
)

# true if this code is running in AWS Batch instead of locally (e.g. make target). In AWS Batch mode, the user has
# already registered a study in the web UI and the import code here is triggered by the S3 file upload.
IS_AWS_DEPLOYMENT = os.environ.get("AWS") is not None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s.%(msecs)03d %(process)d %(levelname)s %(name)s:%(lineno)d %(message)s",
    datefmt="%Y%m%d-%H%M%S",
    handlers=[
        logging.FileHandler("./study_import.log"),
        logging.StreamHandler()
    ] if IS_AWS_DEPLOYMENT else [logging.StreamHandler()]
)

if IS_AWS_DEPLOYMENT:
    engine = get_aws_db_engine()
else:
    engine = get_local_db_engine(os.environ.get("ENGINE_URL"))


def get_all_regions(tax_id):
    existing_region_df = pd.read_sql(
        "select omics_id, display_symbol as region from omics_base where tax_id=%(tax_id)s and omics_type='region'",
        engine,
        params={"tax_id": tax_id},
        index_col="omics_id",
    )

    cols = existing_region_df.columns.tolist()
    tmp = existing_region_df.reset_index()
    match_existing_region_df = (
        tmp.melt(value_vars=cols, id_vars=["omics_id"], value_name="match_id")[
            ["omics_id", "match_id"]
        ]
        .drop_duplicates()
        .set_index("match_id")
    )
    return match_existing_region_df


def get_all_regions_ids_already_in(tax_id):
    existing_region_ids_df = pd.read_sql(
        "select region_id from omics_region_gene", engine
    )
    return existing_region_ids_df.region_id.tolist()


def import_region_base_and_study(study_id: int, data: AnnData, metadata: Dict):
    # for the moment just put all in as the genomic ranges changes from study to study
    # future we could intersect with transcription factor binding annotation and use
    # those coordinates
    logging.info("importing genomic ranges")

    match_existing_region_df = get_all_regions(metadata["taxonomy_id"])

    omics_df = pd.read_sql(
        "select omics_id, ensembl_gene_id, entrez_gene_ids, hgnc_symbols from omics_all where tax_id=%(tax_id)s and omics_type='gene'",
        engine,
        params={"tax_id": int(metadata["taxonomy_id"])},
        index_col="omics_id",
    )
    cols = omics_df.columns.tolist()
    tmp = omics_df.explode("entrez_gene_ids").explode("hgnc_symbols").reset_index()
    match_df = (
        tmp.melt(value_vars=cols, id_vars=["omics_id"], value_name="match_id")[
            ["omics_id", "match_id"]
        ]
        .drop_duplicates()
        .set_index("match_id")
    )

    df_region = (
        data.uns["atac"]["peak_annotation"]
        .reset_index()[["peak", "gene_name"]]
        .rename(columns={"peak": "region"})
    )

    df_region = df_region.merge(
        match_df, left_on="gene_name", right_index=True, how="left"
    )
    df_region.omics_id = df_region.omics_id.astype("Int64")

    df_region[
        ["chromosome", "start_position", "end_position"]
    ] = df_region.region.str.extract(r"(.+):(\d+)-(\d+)", expand=True)
    df_region["omics_type"] = "region"
    df_region["tax_id"] = int(metadata["taxonomy_id"])
    df_region[["display_name"]] = df_region[["region"]]
    df_region[["display_symbol"]] = df_region[["region"]]

    # get h5ad_var_index
    h5ad_index = data.var.reset_index(names="h5ad_var_key").reset_index(
        names="h5ad_var_index"
    )[["h5ad_var_index", "h5ad_var_key"]]

    # insert into base
    omics_base_insert = df_region[
        ["omics_type", "tax_id", "display_name", "display_symbol"]
    ].drop_duplicates()
    # but not the ones which are already inserted
    regions_not_in_df = omics_base_insert.loc[
                        ~omics_base_insert.display_name.isin(match_existing_region_df.index), :
                        ]
    import_df(regions_not_in_df, "omics_base", engine)

    match_existing_region_df = get_all_regions(int(metadata["taxonomy_id"]))
    omics_base_insert = omics_base_insert.merge(
        match_existing_region_df, left_on="display_name", right_index=True
    )

    # insert into omics_region
    omics_region_insert = (
        df_region[["chromosome", "start_position", "end_position", "region"]]
        .drop_duplicates()
        .reset_index(drop=True)
        .merge(
            omics_base_insert[["display_name", "omics_id"]],
            left_on="region",
            right_on="display_name",
        )
        .rename(columns={"omics_id": "region_id"})
        .drop_duplicates()
        .drop("display_name", axis=1)
    )
    # but not the ones which are already inserted
    omics_regions_not_in_df = omics_region_insert.loc[
                              omics_region_insert.region.isin(regions_not_in_df.display_name), :
                              ]
    import_df(omics_regions_not_in_df, "omics_region", engine)

    # insert into omics_region_gene
    omics_region_gene_insert = (
        df_region[["omics_id", "region"]]
        .dropna()
        .drop_duplicates()
        .reset_index(drop=True)
        .rename(columns={"omics_id": "gene_id"})
        .merge(
            omics_base_insert.set_index("display_symbol"),
            left_on="region",
            right_index=True,
        )[["gene_id", "omics_id"]]
        .rename(columns={"omics_id": "region_id"})
        .drop_duplicates()
    )
    # but not the ones in
    all_region_ids = get_all_regions_ids_already_in(metadata["taxonomy_id"])
    omics_region_gene_not_in_df = omics_region_gene_insert.loc[
                                  ~omics_region_gene_insert.region_id.isin(all_region_ids), :
                                  ]
    import_df(omics_region_gene_not_in_df, "omics_region_gene", engine)

    # insert into study_omics
    study_omics_insert = (
        omics_base_insert.merge(
            h5ad_index.set_index("h5ad_var_key"),
            left_on="display_name",
            right_index=True,
        )[["omics_id", "h5ad_var_index", "display_name"]]
        .drop_duplicates()
        .rename(columns={"display_name": "h5ad_var_key"})
    )
    study_omics_insert["study_id"] = study_id
    import_df(study_omics_insert.drop("h5ad_var_key", axis=1), "study_omics", engine)
    return study_omics_insert[["h5ad_var_index", "h5ad_var_key", "omics_id"]]


def import_study_omics_genes(study_id: int, data: AnnData, metadata: Dict):
    logging.info("importing gene definitions of study")
    omics_df = pd.read_sql(
        "select omics_id, ensembl_gene_id, entrez_gene_ids, hgnc_symbols from omics_all where tax_id=%(tax_id)s and omics_type='gene'",
        engine,
        params={"tax_id": int(metadata["taxonomy_id"])},
        index_col="omics_id",
    )

    # generate the mapping gene identifier to omics_id from database
    # could also be done more pandas-like
    #
    # cols = omics_df.columns.tolist()
    # tmp = omics_df.explode('entrez_gene_ids').explode('hgnc_symbols').reset_index()
    # match_df = tmp.melt(value_vars=cols, id_vars=['omics_id'], value_name='match_id')[['omics_id', 'match_id']] \
    #    .drop_duplicates() \
    #    .set_index('match_id')
    match_dfs = []

    for col in ["ensembl_gene_id", "entrez_gene_ids", "hgnc_symbols"]:
        match_df = pd.DataFrame(omics_df[[col]])
        match_df.rename(columns={col: "match_id"}, inplace=True)
        match_df["omics_id"] = match_df.index
        match_df = match_df.explode("match_id", ignore_index=True)
        match_df.drop_duplicates("match_id", inplace=True)
        match_df.set_index("match_id", inplace=True)
        match_dfs.append(match_df)
    match_df = pd.concat(match_dfs)

    # now generate the dataframe to be inserted
    data_genes_df = data.var.copy()
    data_genes_df = data_genes_df.reset_index(names="h5ad_var_key")
    data_genes_df = data_genes_df.reset_index(names="h5ad_var_index")
    data_genes_df = data_genes_df.merge(
        match_df, how="inner", left_on="h5ad_var_key", right_index=True
    )
    data_genes_df.drop_duplicates("omics_id", inplace=True)
    data_genes_df["study_id"] = study_id
    import_df(
        data_genes_df[["h5ad_var_index", "omics_id", "study_id"]], "study_omics", engine
    )
    return data_genes_df[["h5ad_var_index", "h5ad_var_key", "omics_id"]]


def import_study_protein_antibody_tag(study_id: int, data: AnnData, metadata: Dict):
    logging.info("importing protein/antibody definitions of study")
    match_df = pd.read_sql(
        "select omics_id, display_symbol from omics_all where tax_id=%(tax_id)s and omics_type='protein_antibody_tag'",
        engine,
        params={"tax_id": int(metadata["taxonomy_id"])},
        index_col="display_symbol",
    )

    # now generate the dataframe to be inserted
    data_genes_df = data.var.copy()
    data_genes_df = data_genes_df.reset_index(names="h5ad_var_key")
    data_genes_df = data_genes_df.reset_index(names="h5ad_var_index")
    data_genes_df = data_genes_df.merge(
        match_df, how="inner", left_on="h5ad_var_key", right_index=True
    )
    data_genes_df.drop_duplicates("omics_id", inplace=True)
    data_genes_df["study_id"] = study_id
    import_df(
        data_genes_df[["h5ad_var_index", "omics_id", "study_id"]], "study_omics", engine
    )
    return data_genes_df[["h5ad_var_index", "h5ad_var_key", "omics_id"]]


def import_projection(data, data_samples_df, study_id, key, modality=None):
    study_sample_ids = data.obs.merge(
        data_samples_df, left_index=True, right_on="h5ad_obs_key"
    )["study_sample_id"].tolist()
    projection_df = pd.DataFrame(
        {
            "study_id": study_id,
            "study_sample_id": study_sample_ids,
            "modality": modality,
            "projection_type": key,
            "projection": data.obsm[f"X_{key}"][:, 0:2].tolist(),
        }
    )
    if f"{key}_density_sampled_indices" in data.uns["cellenium"]:
        projection_df["display_subsampling"] = False
        projection_df.loc[
            data.uns["cellenium"][f"{key}_density_sampled_indices"],
            "display_subsampling",
        ] = True
    else:
        projection_df["display_subsampling"] = True
    import_df(projection_df, "study_sample_projection", engine)


def _projection_list(data: AnnData | MuData, filetype="h5ad"):
    if isinstance(data, AnnData):
        return (
            data.uns["cellenium"].get("import_projections", np.array(["umap"])).tolist()
        )
    else:
        tmp = data.uns["cellenium"].get("import_projections")
        retlist = []
        for k in tmp.keys():
            retlist.extend([f"{k}:{proj}" for proj in tmp[k]])
        return retlist


def import_study_sample(study_id: int, data: AnnData | MuData):
    logging.info("importing sample definitions")
    data_samples_df = data.obs.copy()
    data_samples_df = data_samples_df.reset_index(names="h5ad_obs_key")
    data_samples_df = data_samples_df.reset_index(names="h5ad_obs_index")
    data_samples_df["study_sample_id"] = range(1, len(data_samples_df) + 1)
    data_samples_df = data_samples_df[
        ["study_sample_id", "h5ad_obs_index", "h5ad_obs_key"]
    ]
    data_samples_df["study_id"] = study_id
    import_df(
        data_samples_df[
            ["study_sample_id", "h5ad_obs_index", "h5ad_obs_key", "study_id"]
        ],
        "study_sample",
        engine,
    )
    with engine.connect() as connection:
        connection.execute(
            text("UPDATE study SET cell_count=:cell_count WHERE study_id=:study_id"),
            {"study_id": study_id, "cell_count": len(data_samples_df)},
        )
    if isinstance(data, AnnData):
        for projection in _projection_list(data):
            import_projection(data, data_samples_df, study_id, projection)
    else:
        for modality, projections in data.uns["cellenium"][
            "import_projections"
        ].items():
            for projection in projections:
                import_projection(
                    data.mod[modality], data_samples_df, study_id, projection, modality
                )

    return data_samples_df


def get_annotation_definition_df(h5ad_columns: List[str]):
    annotation_definition_df = pd.read_sql(
        """select a.annotation_group_id, a.h5ad_column, av.annotation_value_id, av.h5ad_value
            from annotation_group a
            join annotation_value av on av.annotation_group_id = a.annotation_group_id
            where a.h5ad_column = any( %(h5ad_columns)s )""",
        engine,
        params={"h5ad_columns": h5ad_columns},
    )

    return annotation_definition_df


def import_study_sample_annotation(study_id: int, data_samples_df, data: AnnData):
    logging.info("importing sample annotations")
    import_sample_annotations = data.uns["cellenium"]["main_sample_attributes"].tolist()
    import_sample_annotations.extend(
        data.uns["cellenium"].get("advanced_sample_attributes", [])
    )
    secondary_sample_attributes = []
    if data.uns["cellenium"].get("secondary_sample_attributes") is not None:
        secondary_sample_attributes = data.uns["cellenium"][
            "secondary_sample_attributes"
        ].tolist()
        import_sample_annotations.extend(secondary_sample_attributes)

    with engine.connect() as connection:
        for annotation_col in import_sample_annotations:
            annotation_col_clean = annotation_col.replace("_", " ")

            r = connection.execute(
                text(
                    """SELECT annotation_group_id
                    FROM annotation_group WHERE h5ad_column=:h5ad_column"""
                ),
                {"h5ad_column": annotation_col},
            ).fetchone()
            if r is None:
                r = connection.execute(
                    text(
                        """INSERT INTO annotation_group (h5ad_column, display_group)
                            VALUES (:h5ad_column, :h5ad_column_display)
                            RETURNING annotation_group_id"""
                    ),
                    {
                        "h5ad_column": f"{annotation_col}",
                        "h5ad_column_display": annotation_col_clean,
                    },
                ).fetchone()
            annotation_group_id = r[0]
            connection.execute(
                text(
                    """INSERT INTO study_annotation_group_ui (study_id, annotation_group_id, is_primary, ordering, differential_expression_calculated)
                                                                    VALUES (:study_id, :annotation_group_id, :is_primary, :ordering, False)"""
                ),
                {
                    "study_id": study_id,
                    "annotation_group_id": annotation_group_id,
                    "is_primary": annotation_col not in secondary_sample_attributes,
                    "ordering": import_sample_annotations.index(annotation_col),
                },
            )

            values = data.obs[annotation_col].unique().tolist()
            for value in values:
                r = connection.execute(
                    text(
                        "SELECT annotation_value_id FROM annotation_value WHERE annotation_group_id=:annotation_group_id AND h5ad_value=:h5ad_value"
                    ),
                    {"annotation_group_id": annotation_group_id, "h5ad_value": value},
                ).fetchone()
                if r is None:
                    connection.execute(
                        text(
                            """INSERT INTO annotation_value (annotation_group_id, h5ad_value, display_value)
                                            VALUES (:annotation_group_id, :h5ad_value, :h5ad_value_display)"""
                        ),
                        {
                            "annotation_group_id": annotation_group_id,
                            "h5ad_value": value,
                            "h5ad_value_display": value.replace("_", " "),
                        },
                    )

    annotation_definition_df = get_annotation_definition_df(import_sample_annotations)

    data_sample_annotations = data.obs.copy()
    data_sample_annotations = data_sample_annotations.merge(
        data_samples_df, left_index=True, right_on="h5ad_obs_key"
    )
    for h5ad_column in import_sample_annotations:
        palette = _get_palette(data, h5ad_column)

        h5ad_one_annotation_df = data_sample_annotations[
            [h5ad_column, "study_sample_id"]
        ].copy()
        one_annotation_definition_df = annotation_definition_df[
            annotation_definition_df.h5ad_column == h5ad_column
            ]
        annotation_df = h5ad_one_annotation_df.merge(
            one_annotation_definition_df, left_on=h5ad_column, right_on="h5ad_value"
        )

        annotation_df["color"] = annotation_df.apply(
            lambda row: palette[row.h5ad_value], axis=1
        )
        annotation_df = annotation_df[
            ["study_sample_id", "annotation_value_id", "color"]
        ].copy()
        annotation_df["study_id"] = study_id
        annotation_df = (
            annotation_df.groupby(["study_id", "annotation_value_id", "color"])[
                "study_sample_id"
            ]
            .apply(list)
            .reset_index()
            .rename(columns={"study_sample_id": "study_sample_ids"})
        )

        import_df(annotation_df, "study_sample_annotation", engine)


def import_study_layer_expression(
        study_id: int,
        layer_name: str | None,
        data_genes_df,
        data_samples_df,
        data: AnnData,
        metadata,
        omics_type,
        study_layer_id: int,
):
    if layer_name is None:
        layer_name = metadata["X_pseudolayer_name"]
        X = data.X
    else:
        X = data.layers[layer_name]
    logging.info(f"importing expression matrix {layer_name} {omics_type}")

    with engine.connect() as connection:
        sparse_X = sparse.csc_matrix(X)

        map_h5ad_var_index_to_omics_index = np.zeros(
            shape=[sparse_X.shape[1]], dtype=np.uint32
        )
        for i, row in data_genes_df.iterrows():
            map_h5ad_var_index_to_omics_index[row["h5ad_var_index"]] = row["omics_id"]

        samples_of_modality_df = data.obs.join(
            data_samples_df.set_index("h5ad_obs_key")[["study_sample_id"]]
        ).reset_index()
        map_h5ad_obs_index_to_studysample_index = np.zeros(
            shape=[sparse_X.shape[0]], dtype=np.uint32
        )
        for i, row in samples_of_modality_df.iterrows():
            map_h5ad_obs_index_to_studysample_index[i] = row["study_sample_id"]

        buffer = io.StringIO()
        buffered_gene_cnt = 0
        cursor = connection.connection.cursor()
        for gene_i in tqdm.tqdm(
                range(0, sparse_X.shape[1]), desc=f'import expression matrix "{layer_name}"'
        ):
            omics_id = map_h5ad_var_index_to_omics_index[gene_i]
            if omics_id > 0:
                csc_gene_data = sparse.find(sparse_X.T[gene_i])
                data_cell_indexes = csc_gene_data[1]
                data_values = csc_gene_data[2]
                if data_values.size > 0:
                    studysample_ids = map_h5ad_obs_index_to_studysample_index[
                        data_cell_indexes
                    ]
                    buffer.write(str(study_layer_id) + "|" + str(omics_id) + "|{")
                    # the [None, :] reshaping is just a workaround to get np.savetxt into writing separators at all
                    np.savetxt(
                        buffer,
                        studysample_ids[None, :],
                        fmt="%i",
                        delimiter=",",
                        newline=" ",
                    )
                    buffer.write("}|{")
                    np.savetxt(
                        buffer,
                        data_values[None, :],
                        fmt="%.8e",
                        delimiter=",",
                        newline=" ",
                    )
                    buffer.write("}\n")
                    buffered_gene_cnt += 1
            if buffered_gene_cnt == 500 or (
                    gene_i == sparse_X.shape[1] - 1 and buffered_gene_cnt > 0
            ):
                buffer.seek(0)
                cursor.copy_from(
                    buffer,
                    "expression",
                    columns=[
                        "study_layer_id",
                        "omics_id",
                        "study_sample_ids",
                        "values",
                    ],
                    sep="|",
                )
                buffer.close()
                buffer = io.StringIO()
                buffered_gene_cnt = 0
        connection.connection.commit()
        buffer.close()
        cursor.close()


def import_differential_expression(study_id: int, data_genes_df, data: AnnData):
    if "differentially_expressed_genes" not in data.uns["cellenium"]:
        return
    logging.info("importing differentially expressed genes")
    df = data.uns["cellenium"]["differentially_expressed_genes"]
    df = df.merge(data_genes_df, left_on="names", right_on="h5ad_var_key")

    annotation_definition_df = get_annotation_definition_df(
        df["attribute_name"].unique().tolist()
    )
    df = df.merge(
        annotation_definition_df,
        left_on=["attribute_name", "ref_attr_value"],
        right_on=["h5ad_column", "h5ad_value"],
    )
    df["study_id"] = study_id
    df.logfoldchanges = df.logfoldchanges.fillna(
        -1
    )  # TODO: this is as muon puts logfoldchange to NaN
    df.rename(
        columns={
            "pvals": "pvalue",
            "pvals_adj": "pvalue_adj",
            "scores": "score",
            "logfoldchanges": "log2_foldchange",
        },
        inplace=True,
    )

    import_df(
        df[
            [
                "study_id",
                "omics_id",
                "annotation_value_id",
                "pvalue",
                "pvalue_adj",
                "score",
                "log2_foldchange",
            ]
        ],
        "differential_expression",
        engine,
    )
    with engine.connect() as connection:
        connection.execute(
            text(
                """UPDATE study_annotation_group_ui SET differential_expression_calculated=True
                                    WHERE study_id = :study_id and annotation_group_id = any (:annotation_group_ids)"""
            ),
            {
                "study_id": study_id,
                "annotation_group_ids": df["annotation_group_id"].unique().tolist(),
            },
        )


def get_or_create_study_id(stored_filename: UPath) -> int:
    with engine.connect() as connection:
        if IS_AWS_DEPLOYMENT:
            r = connection.execute(
                text(
                    """
                    SELECT study_id FROM study WHERE import_file = :filename
                    """
                ),
                dict(filename=str(stored_filename))
            )
            study_row = r.fetchone()
            if study_row is None:
                raise Exception(f"Study with filename {stored_filename} not found")
            study_id = study_row[0]
            return study_id
        else:
            r = connection.execute(
                text(
                    """INSERT INTO study (filename, study_name)
                VALUES (:filename, :filename)
                RETURNING study_id"""
                ),
                dict(filename=str(stored_filename)),
            )
            study_id = r.fetchone()[0]
            return study_id



def update_study_from_file(study_id: int, data: AnnData):
    def _config_optional_list(key: str):
        if data.uns["cellenium"].get(key) is not None:
            return data.uns["cellenium"][key].tolist()
        return None

    study_data = {
        "study_id": study_id,
        "study_name": data.uns["cellenium"]["title"],
        "description": data.uns["cellenium"]["description"],
        "tissue_ncit_ids": data.uns["cellenium"]["ncit_tissue_ids"].tolist(),
        "disease_mesh_ids": data.uns["cellenium"]["mesh_disease_ids"].tolist(),
        "organism_tax_id": data.uns["cellenium"]["taxonomy_id"],
        "projections": _projection_list(data),
        "reader_permissions": _config_optional_list(
            "initial_reader_permissions"
        ),
        "admin_permissions": _config_optional_list("initial_admin_permissions"),
        "legacy_config": Json(
            data.uns["cellenium"].get("legacy_config"),
            dumps=lambda data: json.dumps(data, cls=NumpyEncoder),
        ),
        "import_finished": False,
    }

    with engine.connect() as connection:
        connection.execute(
            text(""" 
                    UPDATE study SET 
                    study_name = :study_name,
                    description = :description,
                    tissue_ncit_ids = :tissue_ncit_ids,
                    disease_mesh_ids = :disease_mesh_ids,
                    organism_tax_id = :organism_tax_id,
                    projections = :projections,
                    /* extend existing arrays to make sure the uploader does not lose permissions herself */
                    reader_permissions = coalesce(reader_permissions, ARRAY[]::text[]) || coalesce(:reader_permissions, ARRAY[]::text[]),
                    admin_permissions = coalesce(admin_permissions, ARRAY[]::text[]) || coalesce(:admin_permissions, ARRAY[]::text[]),
                    legacy_config = :legacy_config,
                    import_finished = :import_finished,
                    import_started = true
                    WHERE study_id = :study_id
                """),
            study_data
        )


def import_study_safe(data: AnnData, study_id: int, filename: str, analyze_database: bool) -> int:
    def _generate_study_layer(study_id, layer_name, omics_type=None):
        with engine.connect() as connection:
            r = connection.execute(
                text(
                    """INSERT INTO study_layer (study_id, layer, omics_type)
                                    VALUES (:study_id, :layer, :omics_type)
                                    RETURNING study_layer_id"""
                ),
                {"study_id": study_id, "layer": layer_name, "omics_type": omics_type},
            )
            study_layer_id = r.fetchone()[0]

            connection.execute(
                text("call add_studylayer_partition(:study_layer_id)"),
                {"study_layer_id": study_layer_id},
            )
            connection.connection.commit()

        return study_layer_id

    logging.info("importing %s as study_id %s", filename, study_id)
    if isinstance(data, MuData):
        modalities = data.uns["cellenium"]["modalities"]
    else:
        modalities = {"rna": "gene"}

    data_samples_df = import_study_sample(study_id, data)
    for modality in modalities.keys():
        if isinstance(data, MuData):
            cur_data = data.mod[modality]
        else:
            cur_data = data
        import_study_sample_annotation(study_id, data_samples_df, cur_data)

    for modality in modalities.items():
        data_type = modality[1]  # the data_type
        if isinstance(data, MuData):
            cur_data = data.mod[modality[0]]
        else:
            cur_data = data
        meta_data = data.uns["cellenium"]
        if data_type == "gene":
            data_genes_df = import_study_omics_genes(study_id, cur_data, meta_data)
            import_differential_expression(study_id, data_genes_df, cur_data)
        elif data_type == "region":
            # since regions from study to study change we import those which are not yet in the database
            data_region_df = import_region_base_and_study(
                study_id, cur_data, meta_data
            )
            import_differential_expression(study_id, data_region_df, cur_data)
        elif data_type == "protein_antibody_tag":
            data_protein_df = import_study_protein_antibody_tag(
                study_id, cur_data, meta_data
            )
            import_differential_expression(study_id, data_protein_df, cur_data)

    study_layer_id = _generate_study_layer(study_id, "default-layer", None)
    for modality in modalities.items():
        omics_type = modality[1]  # the data_type
        if isinstance(data, MuData):
            cur_data = data.mod[modality[0]]
            cur_data_samples_df = data_samples_df
        else:
            cur_data = data
            cur_data_samples_df = data_samples_df.loc[
                                  data_samples_df.h5ad_obs_key.isin(cur_data.obs.index), :
                                  ]

        meta_data = data.uns["cellenium"]
        if omics_type == "gene":
            import_study_layer_expression(
                study_id,
                None,
                data_genes_df,
                cur_data_samples_df,
                cur_data,
                meta_data,
                omics_type,
                study_layer_id,
            )
            # import more layers if they are present in the dataset (gene expression only
            # see https://github.com/Bayer-Group/cellenium/issues/24
            for layer_name in cur_data.layers.keys():
                further_study_layer_id = _generate_study_layer(
                    study_id, layer_name, "gene"
                )
                import_study_layer_expression(
                    study_id,
                    layer_name,
                    data_genes_df,
                    cur_data_samples_df,
                    cur_data,
                    meta_data,
                    omics_type,
                    further_study_layer_id,
                )

        if omics_type == "region":
            import_study_layer_expression(
                study_id,
                None,
                data_region_df,
                cur_data_samples_df,
                cur_data,
                meta_data,
                omics_type,
                study_layer_id,
            )
        if omics_type == "protein_antibody_tag":
            import_study_layer_expression(
                study_id,
                None,
                data_protein_df,
                cur_data_samples_df,
                cur_data,
                meta_data,
                omics_type,
                study_layer_id,
            )

    with engine.connect() as connection:
        connection.execute(
            text("UPDATE study SET visible=True, import_failed=False, import_finished=True WHERE study_id=:study_id"),
            {"study_id": study_id},
        )
        if analyze_database:
            logging.info("updating postgres statistics...")
            connection.execute(text("call _analyze_schema()"))
    return study_id


def import_study(filename: str, analyze_database: bool) -> int:
    """
    TODO enable S3 access as we had before:

    adata = h5ad_open.h5ad_read(filename)
    stored_filename = filename
    if stored_filename.startswith('scratch'):
        # filename inside scratch (scratch will be /h5ad_store in postgres docker)
        stored_filename = Path(filename).relative_to("scratch").as_posix()
    """

    logging.info(f"importing study from file {filename}")
    data = h5ad_open.h5ad_h5mu_read(filename)
    stored_filename = UPath(filename)
    if filename.startswith("scratch"):
        # filename inside scratch (scratch will be /h5ad_store in postgres docker)
        stored_filename = UPath(filename).relative_to("scratch").as_posix()

    study_id = get_or_create_study_id(stored_filename)
    try:
        update_study_from_file(study_id, data)
        import_study_safe(data, study_id, filename, analyze_database)
    except Exception as e:
        logging.exception("import failed for file %s", filename)
        if IS_AWS_DEPLOYMENT:
            with engine.connect() as connection:
                log = str()
                if pathlib.Path("./study_import.log").exists():
                    with open("./study_import.log", "r") as f:
                        log = f.read()
                connection.execute(
                    text(
                        "UPDATE study SET import_failed=True, import_finished=True, import_log=:log WHERE study_id=:study_id"
                    ),
                    {"study_id": study_id, "log": log},
                )
        exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="cellenium study import tool")
    parser.add_argument(
        "filename",
        help="h5ad/h5mu file created for cellenium (e.g. using a jupyter lab notebook).",
        type=str,
    )
    parser.add_argument(
        "--analyze-database",
        help="analyses the database schema after insert of study",
        action="store_true",
    )
    args = parser.parse_args()
    import_study(args.filename, args.analyze_database)
    logging.info("done")
