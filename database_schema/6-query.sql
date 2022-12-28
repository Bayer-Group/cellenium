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
       s.tissue_ncit_ids,
       ont_tissue.labels      tissue_labels,
       ont_tissue.parent_ids  tissue_parent_ids,
       s.disease_mesh_ids,
       ont_disease.labels     disease_labels,
       ont_disease.parent_ids disease_parent_ids
from study s
         cross join ont_codes_info('NCIT', s.tissue_ncit_ids) ont_tissue
         cross join ont_codes_info('MeSH', s.disease_mesh_ids) ont_disease;

drop view if exists _all_used_ontology_ids cascade;
create view _all_used_ontology_ids
as
select (select array_agg(distinct i)
        from (select i
              from study_overview
                       cross join unnest(tissue_ncit_ids) i
              union
              select i
              from study_overview
                       cross join unnest(tissue_parent_ids) i) x)  "all_tissue_ids",
       (select array_agg(distinct i)
        from (select i
              from study_overview
                       cross join unnest(disease_mesh_ids) i
              union
              select i
              from study_overview
                       cross join unnest(disease_parent_ids) i) x) "all_disease_ids";


create view tree_tissues
as
select *
from concept_hierarchy_minimum_trees_parents_lists('NCIT',
                                                   (select all_tissue_ids from _all_used_ontology_ids));

create view tree_diseases
as
select *
from concept_hierarchy_minimum_trees_parents_lists('MeSH',
                                                   (select all_disease_ids from _all_used_ontology_ids));
