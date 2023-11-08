-- these functions can be developed in postgres_python_stats_devenv.ipynb
CREATE OR REPLACE FUNCTION expression_ttest(p_study_id int, p_study_layer_id int, p_omics_id int,
                                            p_annotation_group_id int,
                                            p_exclude_annotation_value_ids int[], p_secondary_annotation_group_id int,
                                            p_sample1_annotation_value_id int,
                                            p_sample1_second_annotation_value_id int,
                                            p_sample2_annotation_value_id int,
                                            p_sample2_second_annotation_value_id int)
    RETURNS text AS
$$

from typing import List, Optional
import pandas as pd
import numpy as np
import scipy.stats as stats


def sql_query(query):
    # postgres data retrieval with consistent output, both in the jupyter development
    # environment (plpy is not available) and at runtime inside a plpython3u stored procedure
    try:
        import plpy
    except:
        from postgres_utils import get_local_db_engine
        from sqlalchemy import text
        with get_local_db_engine(None).connect() as connection:
            r = connection.execute(text(query))
            return [row._mapping for row in r.fetchall()]
    r = plpy.execute(query)
    return [row for row in r]


def get_expression_values(study_id: int, study_layer_id: int, omics_id: int, annotation_group_id: int,
                          exclude_annotation_value_ids: List[int],
                          secondary_annotation_group_id: Optional[int] = None):
    expression_records = sql_query(f"""
        select annotation_value_id, second_annotation_value_id, values
        from expression_by_two_annotations({study_id}, {study_layer_id}, ARRAY[{omics_id}],
            {annotation_group_id}, {annotation_group_id if secondary_annotation_group_id is None else secondary_annotation_group_id},
            ARRAY{exclude_annotation_value_ids}::int[],
            false)
        """)
    df = pd.DataFrame(expression_records)
    return df


def expression_ttest(study_id: int, study_layer_id: int, omics_id: int, annotation_group_id: int,
                     exclude_annotation_value_ids: List[int],
                     secondary_annotation_group_id: Optional[int],
                     sample1_annotation_value_id: int,
                     sample1_second_annotation_value_id: Optional[int],
                     sample2_annotation_value_id: int,
                     sample2_second_annotation_value_id: Optional[int]):
    df = get_expression_values(study_id, study_layer_id, omics_id, annotation_group_id, exclude_annotation_value_ids,
                               secondary_annotation_group_id)

    def get_group_values(annotation_value_id: int, second_annotation_value_id: Optional[int]):
        mask = df.annotation_value_id == annotation_value_id
        if secondary_annotation_group_id is not None:
            if second_annotation_value_id is not None:
                mask = mask & (df.second_annotation_value_id == second_annotation_value_id)
            else:
                mask = mask & pd.isna(df.second_annotation_value_id)
        group = df[mask]['values'].tolist()
        if len(group) == 0:
            return None
        return group[0]

    sample1 = get_group_values(sample1_annotation_value_id, sample1_second_annotation_value_id)
    sample2 = get_group_values(sample2_annotation_value_id, sample2_second_annotation_value_id)
    if sample1 is None or sample2 is None:
        return "no expression values in one of the groups"
    pvalue = stats.ttest_ind(sample1, sample2, equal_var=False).pvalue
    if pd.isna(pvalue):
        return "NaN"
    return f"{pvalue:.2e}"


return expression_ttest(p_study_id, p_study_layer_id, p_omics_id,
                        p_annotation_group_id,
                        p_exclude_annotation_value_ids, p_secondary_annotation_group_id,
                        p_sample1_annotation_value_id,
                        p_sample1_second_annotation_value_id,
                        p_sample2_annotation_value_id,
                        p_sample2_second_annotation_value_id)
$$ LANGUAGE plpython3u
    IMMUTABLE
    SECURITY INVOKER
    PARALLEL SAFE;

