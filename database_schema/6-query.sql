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
       s.description
from study s;


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
