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
select min(labels), min(parent_ids)
from (
         -- in case there are no parents, still find the label
         select array_agg(distinct c.label) labels,
                null::text[]                parent_ids
         from ontology o
                  join concept c on o.ontid = c.ontid
         where o.name = p_ontology
           and c.ont_code = any (p_ont_codes)
         union all
         select array_agg(distinct c.label) labels,
                array_agg(parents.ont_code) parent_ids
         from ontology o
                  join concept c on o.ontid = c.ontid
                  cross join concept_all_parents(c) parents
         where o.name = p_ontology
           and c.ont_code = any (p_ont_codes)) x
$$;


drop view if exists study_overview cascade;
create view study_overview
            with (security_invoker = true) -- use current user's permission when querying study table
as
select s.study_id,
       s.study_name,
       s.description,
       s.external_website,
       s.cell_count,
       (select min(sl.study_layer_id)
        from study_layer sl
        where sl.study_id = s.study_id) default_study_layer_id,
       s.metadata
from study s
where s.visible = True;
grant select on study_overview to postgraphile;

create materialized view study_overview_ontology
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
         cross join ont_codes_info('MeSH', s.disease_mesh_ids) ont
union all
SELECT s.study_id,
       'taxonomy'                        ontology,
       ARRAY [s.organism_tax_id]::text[] ont_codes,
       ont.labels,
       ont.parent_ids
FROM study s
         cross join ont_codes_info('taxonomy', ARRAY [s.organism_tax_id]::text[]) ont
union all
select study_cell_ontology_ids.study_id,
       'CO'                          ontology,
       study_cell_ontology_ids.ont_codes,
       ont.labels,
       array_agg(distinct parent_id) parent_ids
from (select ssa.study_id, array_agg(distinct c.ont_code) ont_codes
      from study_sample_annotation ssa
               join annotation_value av on ssa.annotation_value_id = av.annotation_value_id and
                                           av.annotation_group_id in (select annotation_group_id
                                                                      from annotation_group
                                                                      where display_group ilike '%cell%type%')
               join concept c on lower(c.label) = lower(av.display_value) and
                                 c.ontid = (select ontid from ontology where name = 'CO')
      group by ssa.study_id) study_cell_ontology_ids
         cross join ont_codes_info('CO', study_cell_ontology_ids.ont_codes) ont
         cross join unnest(ont.parent_ids) parent_id
group by study_cell_ontology_ids.study_id,
         study_cell_ontology_ids.ont_codes,
         ont.labels;
comment on materialized view study_overview_ontology is E'@foreignKey (study_id) references study_overview (study_id)|@fieldName study|@foreignFieldName studyOntology';
create index study_overview_ontology_1 on study_overview_ontology (study_id);
grant select on study_overview_ontology to postgraphile;


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
grant select on tree_ontology to postgraphile;

create view study_annotation_frontend_group
    with (security_invoker = true)
as
select gui.study_id,
       gui.annotation_group_id,
       gui.is_primary,
       gui.ordering,
       gui.differential_expression_calculated,
       g.display_group,
       ug.created_by_user,
       ug.created_by_user = current_user_email() current_user_is_owner,
       ug.private_to_user
from study_annotation_group_ui gui
         join annotation_group g on gui.annotation_group_id = g.annotation_group_id
         left join user_annotation_group ug on ug.saved_as_annotation_group_id = g.annotation_group_id
where (ug.saved_as_annotation_group_id is null
    or ug.private_to_user = False
    or ug.created_by_user = current_user_email()
          );
comment on view study_annotation_frontend_group is
    E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName annotationGroups';
grant select on study_annotation_frontend_group to postgraphile;


create view study_annotation_frontend_value
as
select ssa.study_id,
       v.annotation_group_id,
       ssa.annotation_value_id,
       v.display_value,
       ssa.color,
       cardinality(ssa.study_sample_ids) sample_count
from study_sample_annotation ssa
         join annotation_value v on ssa.annotation_value_id = v.annotation_value_id;
comment on view study_annotation_frontend_value is
    E'@foreignKey (study_id, annotation_group_id) references study_annotation_frontend_group (study_id, annotation_group_id)|@fieldName group|@foreignFieldName annotationValues';
grant select on study_annotation_frontend_value to postgraphile;


drop view if exists study_admin_details cascade;
create view study_admin_details
as
select s.study_id,
       s.study_name,
       s.filename,
       s.description,
       s.external_website,
       s.cell_count,
       s.visible,
       s.reader_permissions,
       s.admin_permissions,
       s.disease_mesh_ids,
       s.tissue_ncit_ids,
       s.import_started,
       s.import_failed,
       s.import_finished,
       (case when s.import_log is not null then True else False end) as has_import_log,
       case when sv.study_id is not null then True else False end       "reader_permission_granted",
       case when sa.study_id is not null then True else False end       "admin_permission_granted"
from study s
         left join study_visible_currentuser sv on sv.study_id = s.study_id
         left join study_administrable_currentuser sa on sa.study_id = s.study_id;
grant select on study_admin_details to postgraphile;

drop view if exists study_import_log cascade;
create view study_import_log
as
select s.import_file,
       s.import_log,
       s.study_id
from study s
         left join study_visible_currentuser sv on sv.study_id = s.study_id
         left join study_administrable_currentuser sa on sa.study_id = s.study_id;
grant select on study_import_log to postgraphile;

create or replace function study_definition_update()
    returns boolean
    language plpgsql
    security definer
    volatile
as
$$
begin
    refresh materialized view study_overview_ontology;
    return true;
end;
$$;
