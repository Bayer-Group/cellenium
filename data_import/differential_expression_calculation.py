from typing import List

import h5ad_open
import h5ad_preparation as prep
import numpy as np
import pandas as pd
from common import engine
from postgres_utils import import_df
from sqlalchemy import text


def read_h5ad(study_id: int):
    with engine.connect() as connection:
        r = connection.execute(
            text("select coalesce(filename, import_file) filename from study where study_id = :study_id"),
            {"study_id": str(study_id)},
        )
        study_row = r.fetchone()
        filename = study_row.filename
    if "/" not in filename:
        # no S3 file and no real path, assume file was imported with demo Makefile
        filename = f"../scratch/{filename}"
    return h5ad_open.h5ad_h5mu_read(filename)


def get_annotation_df(study_id: int, annotation_group_id: int):
    annotation_df = pd.read_sql(
        """
        SELECT ssa.annotation_value_id, ss.h5ad_obs_index
        FROM study_sample_annotation ssa
                 join annotation_value av on av.annotation_value_id = ssa.annotation_value_id
                 cross join UNNEST(study_sample_ids) as sample_id
                 join study_sample ss on ss.study_id = ssa.study_id and ss.study_sample_id = sample_id
        WHERE ssa.study_id = %(study_id)s and av.annotation_group_id= %(annotation_group_id)s""",
        engine,
        params={"study_id": study_id, "annotation_group_id": annotation_group_id},
        index_col="h5ad_obs_index",
    )
    return annotation_df


def save_annotation_group_diffexp_status(study_id: int, annotation_group_ids: List[int]):
    with engine.connect() as connection:
        # mark annotation group to have differentially expressed genes; detect if any pairwise calculations were done
        connection.execute(
            text(
                """UPDATE study_annotation_group_ui SET differential_expression_calculated=True
                   WHERE study_id = :study_id and annotation_group_id = any (:annotation_group_ids);
                   UPDATE study_annotation_group_ui set pairwise_differential_expression_calculated=false
                   where study_id = :study_id;
                   UPDATE study_annotation_group_ui set pairwise_differential_expression_calculated=true
                   where study_id = :study_id and annotation_group_id in (
                        select distinct av.annotation_group_id
                        from differential_expression de
                        join annotation_value av on de.other_annotation_value_id = av.annotation_value_id
                        where study_id = :study_id);
                   """
            ),
            {
                "study_id": study_id,
                "annotation_group_ids": annotation_group_ids,
            },
        )


def differential_expression_calculation(study_id: int, annotation_group_ids: List[int]):
    adata = read_h5ad(study_id)

    genes_df = pd.read_sql(
        """select omics_id, h5ad_var_index
                              from study_omics where study_id= %(study_id)s""",
        engine,
        params={"study_id": study_id},
        index_col="h5ad_var_index",
    )
    genes_index_omics_id_df = adata.var.reset_index().join(genes_df)[["index", "omics_id"]]

    for annotation_group_id in annotation_group_ids:
        # The annotation group we're considering here can be one of the study sample annotations, or it has been
        # dynamically added as a user annotation in cellenium. In any case, we just add another sample annotation
        # column, with the categories being the sample annotation IDs.
        annotation_df = get_annotation_df(study_id, annotation_group_id)
        adata.obs = (
            adata.obs.reset_index()
            .join(annotation_df)[["annotation_value_id"]]
            .fillna("0")
            .astype({"annotation_value_id": "int32"})
            .astype("str")
        )
        prep.calculate_differentially_expressed_genes(adata, ["annotation_value_id"])
        diffexp_df = adata.uns["cellenium"]["differentially_expressed_genes"]
        diffexp_df = diffexp_df.dropna(subset="ref_attr_value")
        diffexp_df.cmp_attr_value = diffexp_df.cmp_attr_value.replace("_OTHERS_", "0")
        df = diffexp_df.merge(genes_index_omics_id_df, left_on="names", right_on="index")
        df = df.dropna(subset="omics_id").astype({"omics_id": "int32", "ref_attr_value": "int32", "cmp_attr_value": "int32"})
        df.rename(
            columns={
                "ref_attr_value": "annotation_value_id",
                "cmp_attr_value": "other_annotation_value_id",
                "pvals": "pvalue",
                "pvals_adj": "pvalue_adj",
                "scores": "score",
                "logfoldchanges": "log2_foldchange",
            },
            inplace=True,
        )
        df.drop_duplicates(subset=["omics_id", "annotation_value_id", "other_annotation_value_id"], inplace=True)
        df.replace(np.inf, 1000, inplace=True)
        df = df[df.annotation_value_id != 0].copy()
        df.other_annotation_value_id.replace(0, None, inplace=True)
        df["study_id"] = study_id
        df.drop(columns=["names", "attribute_name", "index"], inplace=True)

        with engine.connect() as connection:
            connection.execute(
                text(
                    """delete from differential_expression
                       where study_id = :study_id and annotation_value_id = any (:annotation_value_ids)"""
                ),
                {"study_id": study_id, "annotation_value_ids": df["annotation_value_id"].unique().tolist()},
            )
        import_df(df, "differential_expression", engine)
        adata.obs.drop(columns=["annotation_value_id"], inplace=True)

    save_annotation_group_diffexp_status(study_id, annotation_group_ids)
