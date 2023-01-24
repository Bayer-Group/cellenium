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
    SECURITY DEFINER
    PARALLEL SAFE;



drop function if exists ont_codes_info;
create function ont_codes_info(p_ontology text, p_ont_codes text[])
    returns table
            (
                labels     text[],
                parent_ids text[]
            )
    LANGUAGE sql
    STABLE
as
$$
select array_agg(distinct c.label) labels,
       array_agg(parents.ont_code) parent_ids
from ontology o
         join concept c on o.ontid = c.ontid
         cross join concept_all_parents(c) parents
where o.name = p_ontology
  and c.ont_code = any (p_ont_codes)
$$;


drop view if exists study_overview cascade;
create view study_overview
as
select s.study_id,
       s.study_name,
       s.description,
       s.cell_count
from study s
where s.visible = True;




CREATE VIEW study_overview_ontology
AS
SELECT s.study_id,
       'NCIT'            ontology,
       s.tissue_ncit_ids ont_codes,
       ont.labels,
       ont.parent_ids
FROM study s
         cross join ont_codes_info('NCIT', s.tissue_ncit_ids) ont
union all
SELECT s.study_id,
       'MeSH'             ontology,
       s.disease_mesh_ids ont_codes,
       ont.labels,
       ont.parent_ids
FROM study s
         cross join ont_codes_info('MeSH', s.disease_mesh_ids) ont;
comment on view study_overview_ontology is E'@foreignKey (study_id) references study_overview (study_id)|@fieldName study|@foreignFieldName studyOntology';


drop view if exists _all_used_ontology_ids cascade;
create view _all_used_ontology_ids
as
select ontology, i ont_code
from study_overview_ontology
         cross join unnest(ont_codes) i
union all
select ontology, i ont_code
from study_overview_ontology
         cross join unnest(parent_ids) i;


create view tree_ontology
as
with ont_code_lists as (select ontology, array_agg(ont_code) ont_codes
                        from _all_used_ontology_ids
                        group by ontology)
select ont_code_lists.ontology, l.*
from ont_code_lists,
     concept_hierarchy_minimum_trees_parents_lists(ont_code_lists.ontology,
                                                   ont_code_lists.ont_codes) l;


-- TODO materialized view, with refresh no study import
create view annotation_value_combination_sample_count
as
with sample_annotationvalues as (select study_id,
                                        sample_id,
                                        array_agg(ssa.annotation_value_id order by ssa.annotation_value_id) annotation_value_combination
                                 from study_sample_annotation ssa
                                          cross join lateral unnest(ssa.study_sample_ids) sample_id
                                 group by study_id, sample_id)
select study_id, annotation_value_combination, count(1)
from sample_annotationvalues
group by study_id, annotation_value_combination;


create view study_annotation_frontend_group
as
select gui.study_id,
       gui.annotation_group_id,
       gui.is_primary,
       gui.ordering,
       gui.differential_expression_calculated,
       g.display_group
from study_annotation_group_ui gui
         join annotation_group g on gui.annotation_group_id = g.annotation_group_id;
comment on view study_annotation_frontend_group is
    E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName annotationGroups';

create view study_annotation_frontend_value
as
select ssa.study_id,
       v.annotation_group_id,
       ssa.annotation_value_id,
       v.display_value,
       v.color,
       cardinality(ssa.study_sample_ids) sample_count
from study_sample_annotation ssa
         join annotation_value v on ssa.annotation_value_id = v.annotation_value_id;
comment on view study_annotation_frontend_value is
    E'@foreignKey (study_id, annotation_group_id) references study_annotation_frontend_group (study_id, annotation_group_id)|@fieldName group|@foreignFieldName annotationValues';
