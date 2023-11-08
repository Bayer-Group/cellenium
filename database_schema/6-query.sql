DROP FUNCTION IF EXISTS annotation_value_coocurrence(study_id int, annotation_group_id_1 int, annotation_value_id_2 int);
CREATE OR REPLACE FUNCTION annotation_value_coocurrence(study_id int, annotation_group_id_1 int, annotation_group_id_2 int)
    RETURNS TABLE
            (
                annotation_value_id_1 INT,
                annotation_value_id_2 INT,
                occurrence            INT
            )
AS
$$

import pandas as pd
import numpy as np


def sql_query(query):
    # postgres data retrieval with consistent output, both in the jupyter development
    # environment (plpy is not available) and at runtime inside a plpython3u stored procedure
    try:
        import plpy
    except:
        from postgres_utils import engine
        from sqlalchemy import text
        with engine.connect() as connection:
            r = connection.execute(text(query))
            return [row._mapping for row in r.fetchall()]
    r = plpy.execute(query)
    return [row for row in r]


data = sql_query(f"""
        SELECT tmp.annotation_value_id, tmp.study_sample_id
          FROM (
            SELECT ssa.annotation_value_id, UNNEST(ssa.study_sample_ids) as study_sample_id
              FROM study_sample_annotation ssa
              JOIN annotation_value av
                ON ssa.annotation_value_id = av.annotation_value_id
             WHERE av.annotation_group_id in ({annotation_group_id_1},{annotation_group_id_2})
               AND ssa.study_id = {study_id}) tmp;
    """)
df = pd.DataFrame(data)
df['count'] = 1
df = df.pivot(index='annotation_value_id', columns='study_sample_id', values='count').fillna(0).astype(int)
df_cooc = df.dot(df.T)
df_cooc.columns.name = df_cooc.index.name = None
out = df_cooc.stack().reset_index()
out.columns = ['annotation_value_id_1', 'annotation_value_id_2', 'occurrence']
out = out.loc[out['occurrence'] != 0, :].astype(int).reset_index(drop=True)
return out.to_records(index=False).tolist()
$$ LANGUAGE plpython3u
    IMMUTABLE
    SECURITY INVOKER
    PARALLEL SAFE;


create or replace function study_definition_update()
    returns boolean
    language plpgsql
    security definer
    volatile
as
$$
begin
    refresh materialized view study_overview_ontology;
    refresh materialized view study_omics_transposed;
    refresh materialized view study_sample_annotation_subsampling;
    refresh materialized view study_sample_projection_subsampling_transposed;
    return true;
end;
$$;