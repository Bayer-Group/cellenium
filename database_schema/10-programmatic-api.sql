DROP MATERIALIZED VIEW IF EXISTS studies_bulk_rna;
CREATE MATERIALIZED VIEW studies_bulk_rna AS
SELECT study_id,
       study_name,
       description,
       external_website,
       import_finished,
       import_failed,
       cell_count,
       tissue_ncit_ids,
       disease_mesh_ids,
       cell_ontology_ids,
       organism_tax_id,
       array(SELECT layer from study_layer WHERE study.study_id = study_layer.study_id) as layers
FROM study
WHERE legacy_config ->> 'studytype' = 'RNA-seq';


DROP MATERIALIZED VIEW IF EXISTS studies_single_cell;
CREATE MATERIALIZED VIEW studies_single_cell AS
SELECT study_id,
       study_name,
       description,
       external_website,
       import_finished,
       import_failed,
       cell_count,
       tissue_ncit_ids,
       disease_mesh_ids,
       cell_ontology_ids,
       organism_tax_id,
       array(SELECT layer from study_layer WHERE study.study_id = study_layer.study_id) as layers
FROM study
WHERE legacy_config ->> 'studytype' = 'scRNA';

