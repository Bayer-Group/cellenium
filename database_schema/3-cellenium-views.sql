-- decide which studies can be analyzed by the
--  current user:
--    all studies which are open to everybody
--  ("reader_permissions is null")
--    and all studies which list of allowed groups
--  has at least one group in common with the
--  current user's groups
CREATE OR REPLACE VIEW study_visible_currentuser
-- see
-- https://www.postgresql.org/docs/current/rules-pr
-- ivileges.html
WITH ( security_barrier
) AS
SELECT
	s.study_id
FROM
	study s
WHERE
	s.reader_permissions IS NULL
	OR s.reader_permissions = ARRAY[]::text[]
	OR s.reader_permissions && current_user_groups ();

GRANT SELECT ON study_visible_currentuser TO postgraphile;

CREATE OR REPLACE VIEW study_administrable_currentuser AS
SELECT
	s.study_id
FROM
	study s
WHERE
	s.admin_permissions IS NULL
	OR s.reader_permissions = ARRAY[]::text[]
	OR s.admin_permissions && current_user_groups ();

GRANT SELECT ON study_administrable_currentuser TO postgraphile;


DROP VIEW IF EXISTS omics_all CASCADE;
CREATE VIEW omics_all AS
SELECT
	b.omics_id,
	b.omics_type,
	b.tax_id,
	b.display_symbol,
	b.display_name,
	og.ensembl_gene_id,
	og.entrez_gene_ids,
	og.hgnc_symbols,
	ogr.region,
	array_remove(array_agg(otfg.gene_id) || array_agg(opatg.gene_id) ||
	array_agg(ogrg.gene_id), NULL) AS linked_genes
FROM
	omics_base b
	LEFT JOIN omics_gene og ON b.omics_id = og.gene_id
	LEFT JOIN omics_region ogr ON b.omics_id = ogr.region_id
	LEFT JOIN omics_region_gene ogrg ON b.omics_id = ogrg.region_id
	LEFT JOIN omics_protein_antibody_tag_gene opatg ON b.omics_id = opatg.protein_antibody_tag_id
	LEFT JOIN omics_transcription_factor_gene otfg ON b.omics_id = otfg.transcription_factor_id
GROUP BY
	og.ensembl_gene_id,
	og.entrez_gene_ids,
	og.hgnc_symbols,
	b.omics_id,
	b.omics_type,
	b.tax_id,
	b.display_symbol,
	b.display_name,
	ogr.region;

GRANT SELECT ON omics_all TO postgraphile;

CREATE OR REPLACE VIEW differential_expression_v WITH ( security_invoker = TRUE
) AS
SELECT
	de.*,
	ob.omics_type,
	ob.display_symbol,
	ob.display_name,
	oa.linked_genes
FROM
	differential_expression de
	JOIN omics_base ob ON de.omics_id = ob.omics_id
	JOIN omics_all oa ON de.omics_id = oa.omics_id;

GRANT SELECT ON differential_expression_v TO postgraphile;
--
-- -- contains all samples which appear in at least
-- -- one projection
-- CREATE OR REPLACE VIEW study_sample_annotation_subsampling AS
-- SELECT
-- 	ssa.study_id,
-- 	ssa.annotation_value_id,
-- 	array_agg(DISTINCT ssp.study_sample_id) study_sample_ids
-- FROM
-- 	study_sample_annotation ssa
-- 	CROSS JOIN unnest(ssa.study_sample_ids) sample_id
-- 	JOIN study_sample_projection ssp ON ssp.study_id = ssa.study_id
-- 		AND ssp.study_sample_id = sample_id
-- WHERE
-- 	ssp.display_subsampling = TRUE
-- GROUP BY
-- 	ssa.study_id,
-- 	ssa.annotation_value_id;
--
-- COMMENT ON VIEW study_sample_annotation_subsampling IS E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName studySampleAnnotationSubsampling';
--
-- GRANT SELECT ON study_sample_annotation_subsampling TO postgraphile;

-- NEW NEW NEW
--  contains all samples which appear in at least
-- one projection


-- DROP VIEW IF EXISTS study_sample_annotation_subsampling;
DROP MATERIALIZED VIEW IF EXISTS study_sample_annotation_subsampling;
CREATE MATERIALIZED VIEW study_sample_annotation_subsampling AS
SELECT
	ssa.study_id,
	ssa.annotation_value_id,
	array_agg(DISTINCT ssp.study_sample_id) study_sample_ids
FROM
	study_sample_annotation ssa
	CROSS JOIN unnest(ssa.study_sample_ids) sample_id
	JOIN study_sample_projection ssp ON ssp.study_id = ssa.study_id
		AND ssp.study_sample_id = sample_id
WHERE
	ssp.display_subsampling = TRUE
GROUP BY
	ssa.study_id,
	ssa.annotation_value_id;

-- COMMENT ON MATERIALIZED VIEW study_sample_annotation_subsampling IS E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName studySampleAnnotationSubsampling';
COMMENT ON MATERIALIZED VIEW study_sample_annotation_subsampling IS NULL;
GRANT SELECT ON study_sample_annotation_subsampling TO postgraphile;
CREATE INDEX study_sample_annotation_subsampling_idx ON study_sample_annotation_subsampling (study_id);

-- CREATE OR REPLACE VIEW study_sample_projection_subsampling_transposed AS
-- SELECT
-- 	study_id,
-- 	projection_type,
-- 	modality,
-- 	array_agg(study_sample_id ORDER BY study_sample_id) study_sample_id,
-- 	array_agg(projection ORDER BY study_sample_id) projection
-- FROM
-- 	study_sample_projection
-- WHERE
-- 	display_subsampling = TRUE
-- GROUP BY
-- 	study_id,
-- 	projection_type,
-- 	modality;
--
-- COMMENT ON VIEW study_sample_projection_subsampling_transposed IS E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName studySampleProjectionSubsamplingTransposed';
-- GRANT SELECT ON study_sample_projection_subsampling_transposed TO postgraphile;

-- DROP VIEW IF EXISTS study_sample_projection_subsampling_transposed;
DROP MATERIALIZED VIEW IF EXISTS study_sample_projection_subsampling_transposed;
CREATE MATERIALIZED VIEW study_sample_projection_subsampling_transposed AS
SELECT
	study_id,
	projection_type,
	modality,
	array_agg(study_sample_id ORDER BY study_sample_id) study_sample_id,
	array_agg(projection ORDER BY study_sample_id) projection
FROM
	study_sample_projection
WHERE
	display_subsampling = TRUE
GROUP BY
	study_id,
	projection_type,
	modality;

-- COMMENT ON MATERIALIZED VIEW study_sample_projection_subsampling_transposed IS E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName studySampleProjectionSubsamplingTransposed';
COMMENT ON MATERIALIZED VIEW study_sample_projection_subsampling_transposed IS NULL;
GRANT SELECT ON study_sample_projection_subsampling_transposed TO postgraphile;
CREATE INDEX study_sample_projection_subsampling_transposed_idx ON study_sample_projection_subsampling_transposed (study_id);

-- CREATE VIEW study_omics_transposed AS
-- SELECT
-- 	study_id,
-- 	array_agg(ob.omics_id ORDER BY ob.omics_id) omics_id,
-- 	array_agg(ob.omics_type ORDER BY ob.omics_id) omics_type,
-- 	array_agg(ob.display_symbol ORDER BY ob.omics_id) display_symbol,
-- 	array_agg(ob.display_name ORDER BY ob.omics_id) display_name
-- FROM
-- 	study_omics
-- 	JOIN omics_base ob ON study_omics.omics_id = ob.omics_id
-- GROUP BY
-- 	study_id;
--
-- COMMENT ON VIEW study_omics_transposed IS E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName studyOmicsTransposed';
-- GRANT SELECT ON study_omics_transposed TO postgraphile;


DROP VIEW IF EXISTS study_omics_transposed;
DROP MATERIALIZED VIEW IF EXISTS study_omics_transposed;
CREATE MATERIALIZED VIEW study_omics_transposed AS
SELECT
	study_id,
	array_agg(ob.omics_id ORDER BY ob.omics_id) omics_id,
	array_agg(ob.omics_type ORDER BY ob.omics_id) omics_type,
	array_agg(ob.display_symbol ORDER BY ob.omics_id) display_symbol,
	array_agg(ob.display_name ORDER BY ob.omics_id) display_name
FROM
	study_omics
	JOIN omics_base ob ON study_omics.omics_id = ob.omics_id
GROUP BY
	study_id;

-- COMMENT ON MATERIALIZED VIEW study_omics_transposed IS E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName studyOmicsTransposed';
COMMENT ON MATERIALIZED VIEW study_omics_transposed IS NULL;
GRANT SELECT ON study_omics_transposed TO postgraphile;
CREATE INDEX study_omics_transposed_idx ON study_omics_transposed (study_id);


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

DROP MATERIALIZED VIEW IF EXISTS study_overview_ontology;
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


-- create or replace view study_annotation_frontend_group AS
-- select gui.study_id,
--        gui.annotation_group_id,
--        gui.is_primary,
--        gui.ordering,
--        gui.differential_expression_calculated,
--        g.display_group,
--        ug.created_by_user,
--        ug.created_by_user = current_user_email() current_user_is_owner,
--        ug.private_to_user
-- from study_annotation_group_ui gui
--          join annotation_group g on gui.annotation_group_id = g.annotation_group_id
--          left join user_annotation_group ug on ug.saved_as_annotation_group_id = g.annotation_group_id
-- where (ug.saved_as_annotation_group_id is null
--     or ug.private_to_user = False
--     or ug.created_by_user = current_user_email()
--           );
-- comment on view study_annotation_frontend_group is
--     E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName annotationGroups';
-- grant select on study_annotation_frontend_group to postgraphile;

DROP VIEW IF EXISTS study_annotation_frontend_group;
create view study_annotation_frontend_group AS
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
-- comment on materialized view study_annotation_frontend_group is E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName annotationGroups';
comment on view study_annotation_frontend_group is NULL;
grant select on study_annotation_frontend_group to postgraphile;

DROP VIEW IF EXISTS study_annotation_frontend_value;
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


