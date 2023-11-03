drop function if exists ensure_text_array cascade;
create function ensure_text_array(a text[]) returns text[]
as
$$
select coalesce(a, ARRAY []::text[])
$$ language sql immutable;



drop view if exists api_study_annotation_overview cascade;
create view api_study_annotation_overview
as
select ssa.study_id,
       ag.h5ad_column                                              annotation_group_column,
       ag.display_group                                            annotation_group_display,
       array_agg(av.h5ad_value order by av.annotation_value_id)    annotation_values,
       array_agg(av.display_value order by av.annotation_value_id) annotation_display_values
from annotation_group ag
         join annotation_value av on ag.annotation_group_id = av.annotation_group_id
         join study_sample_annotation ssa on ssa.annotation_value_id = av.annotation_value_id
group by ssa.study_id, ag.h5ad_column, ag.display_group;
grant select on api_study_annotation_overview to postgraphile;

drop view if exists api_studies_base cascade;
create view api_studies_base
    with (security_invoker = true)
as
SELECT study_id,
       study_name,
       coalesce(legacy_config ->> 'studytype', 'scRNA')                                    study_type,
       description,
       external_website,
       cell_count,
       (select array_agg(asao)
        from api_study_annotation_overview asao
        where asao.study_id = study.study_id)                                              annotations,
       ensure_text_array(tissue_ncit_ids)                                                  tissue_ncit_ids,
       ensure_text_array((select so.labels
                          from study_overview_ontology so
                          where so.ontology = 'NCIT'
                            and so.study_id = study.study_id))                             tissue_ncit_labels,
       ensure_text_array(disease_mesh_ids)                                                 disease_mesh_ids,
       ensure_text_array((select so.labels
                          from study_overview_ontology so
                          where so.ontology = 'MeSH'
                            and so.study_id = study.study_id))                             disease_mesh_labels,
       ensure_text_array(cell_ontology_ids)                                                cell_ontology_ids,
       ensure_text_array((select so.labels
                          from study_overview_ontology so
                          where so.ontology = 'CO'
                            and so.study_id = study.study_id))                             cell_ontology_labels,
       organism_tax_id,
       (select so.labels[1]
        from study_overview_ontology so
        where so.ontology = 'taxonomy'
          and so.study_id = study.study_id)                                                organism_label,
       array(SELECT layer from study_layer WHERE study.study_id = study_layer.study_id) as layers
FROM study
where study.visible = True;
grant select on api_studies_base to postgraphile;

create view api_studies_single_cell
    with (security_invoker = true)
as
SELECT *
from api_studies_base
where study_type = 'scRNA';
grant select on api_studies_single_cell to postgraphile;

create view api_studies_bulk_rna
    with (security_invoker = true)
as
SELECT *
from api_studies_base
where study_type = 'RNA-seq';
grant select on api_studies_single_cell to postgraphile;

/*
query StudyApiExample {
  apiStudiesSingleCellsList {
    studyId
    studyName
    description
    cellCount
    layers
    annotations {
      annotationGroupColumn
      annotationGroupDisplay
      annotationValues
      annotationDisplayValues
    }
    tissueNcitIds
    tissueNcitLabels
    diseaseMeshIds
    diseaseMeshLabels
    organismLabel
    organismTaxId
    externalWebsite
    cellOntologyIds
    cellOntologyLabels
  }
}
 */


--diff_genes = cellenium.differentially_expressed_genes(study_id = study_id, annotation = attribute, annotation_value = annotation_value)

drop view if exists api_omics cascade;
create view api_omics
as
select omics_id,
       omics_type,
       tax_id,
       display_symbol,
       display_name,
       ensembl_gene_id,
       entrez_gene_ids,
       hgnc_symbols,
       region,
       linked_genes
from omics_all;
grant select on api_omics to postgraphile;

drop view if exists api_differential_expression;
create view api_differential_expression
    with (security_invoker = true)
as
select de.study_id,
       s.study_name,
       s.organism_tax_id                                                        study_tax_id,
       de.omics_id                                                              cellenium_internal_omics_id,
       -- Typically genes are measured, not regions or proteins. Simple shortcut for the typical case:
       omics_gene.ensembl_gene_id,
       omics_gene.entrez_gene_ids,
       omics_gene.hgnc_symbols,
       -- All multi omics details can be retrieved in a nested json object:
       (select api_omics from api_omics where api_omics.omics_id = de.omics_id) omics_details,
       ag.h5ad_column                                                           annotation_group_column,
       ag.display_group                                                         annotation_group_display,
       av.h5ad_value                                                            annotation_value,
       av.display_value                                                         annotation_display_value,
       de.log2_foldchange,
       de.pvalue_adj,
       de.pvalue,
       de.score
from differential_expression de
         join study s on de.study_id = s.study_id
         join annotation_value av on de.annotation_value_id = av.annotation_value_id
         join annotation_group ag on av.annotation_group_id = ag.annotation_group_id
         left join omics_gene on de.omics_id = omics_gene.gene_id;
grant select on api_differential_expression to postgraphile;

/*
query DiffExpApiExample {
  apiDifferentialExpressionsList(
    filter: {studyId: {equalTo: 3}, annotationGroupDisplay: {equalToInsensitive: "celltype"}}
    orderBy: SCORE_DESC
    first: 10
  ) {
    ensemblGeneId
    hgncSymbols
    pvalueAdj
    score
    annotationDisplayValue
  }
}
 */

drop type if exists api_expression_by_annotation cascade;
create type api_expression_by_annotation as
(
    study_layer_id           int,
    study_id                 int,
    layer                    text,
    omics_id                 int,
    ensembl_gene_id          text,
    entrez_gene_ids          text[],
    hgnc_symbols             text[],
    annotation_value_id      int,
    annotation_display_value text,
    values                   real[],
    boxplot_params           boxplot_values,
    median                   real,
    q3                       real,
    mean                     real,
    value_count              int,
    expr_samples_fraction    real
);


drop function if exists api_expression_by_annotation;
create function api_expression_by_annotation(p_study_ids int[],
                                             p_annotation_group_column text,
                                             p_ensembl_gene_ids text[] default null,
                                             p_hgnc_symbols text[] default null,
                                             p_entrez_gene_ids text[] default null,
                                             p_omics_ids int[] default null,
                                             p_layer_name text default null
)
    returns setof api_expression_by_annotation
    language sql
    stable
as
$$
select expression_by_annotation.study_layer_id,
       study_layer.study_id,
       study_layer.layer,
       expression_by_annotation.omics_id,
       omics_gene.ensembl_gene_id,
       omics_gene.entrez_gene_ids,
       omics_gene.hgnc_symbols,
       expression_by_annotation.annotation_value_id,
       expression_by_annotation.annotation_display_value,
       expression_by_annotation.values,
       expression_by_annotation.boxplot_params,
       expression_by_annotation.median,
       expression_by_annotation.q3,
       expression_by_annotation.mean,
       expression_by_annotation.value_count,
       expression_by_annotation.expr_samples_fraction
from expression_by_annotation(
             (select ARRAY [min(study_layer_id)]
              from study_layer
              where study_id = any (p_study_ids)
                and (p_layer_name is null or layer = p_layer_name)
              group by study_id),
             (select array_agg(gene_id)
              from omics_gene
              where gene_id = any (p_omics_ids)
                 or hgnc_symbols && p_hgnc_symbols
                 or omics_gene.entrez_gene_ids && p_entrez_gene_ids
                 or ensembl_gene_id = any (p_ensembl_gene_ids)),
             (select annotation_group_id from annotation_group where h5ad_column = p_annotation_group_column),
             ARRAY []::int[])
         join study_layer on expression_by_annotation.study_layer_id = study_layer.study_layer_id
         left join omics_gene on expression_by_annotation.omics_id = omics_gene.gene_id
    ;
$$;

/*
select study_layer_id, omics_id, annotation_value_id, annotation_display_value, q3, expr_samples_fraction
from api_expression_by_annotation(ARRAY [3], 'celltype', p_hgnc_symbols := ARRAY ['ALB','ENG']);

query ExpressionDataAggregatedApiExample {
  apiExpressionByAnnotationList(
    pStudyIds: [3]
    pHgncSymbols: ["ALB", "ENG"]
    pAnnotationGroupColumn: "celltype"
  ) {
    hgncSymbols
    annotationDisplayValue
    median
    exprSamplesFraction
  }
}
 */

DROP VIEW IF EXISTS api_study_h5_download;
DROP FUNCTION IF EXISTS create_study_h5ad_presigned_url;
CREATE OR REPLACE FUNCTION create_study_h5ad_presigned_url(IN file_url text)
    RETURNS text
    LANGUAGE plpython3u
    VOLATILE
AS $$
    import boto3
    from botocore.client import Config
    import os

    if 'S3_BUCKET' not in os.environ:
        raise Exception("S3_BUCKET environment variable not set")
    if 'AWS_REGION' not in os.environ:
        raise Exception("AWS_REGION environment variable not set")

    kwargs = dict()
    if all([key in os.environ for key in ["AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN"]]):
        kwargs = dict(aws_access_key_id=os.environ['AWS_ACCESS_KEY_ID'],
                      aws_secret_access_key=os.environ['AWS_SECRET_ACCESS_KEY'],
                      aws_session_token=os.environ['AWS_SESSION_TOKEN'])
    if file_url is None or not "s3://" in file_url:
     return file_url
    bucket = os.environ["S3_BUCKET"]
    s3_client = boto3.client("s3", config=Config(signature_version='s3v4'), region_name=os.environ['AWS_REGION'], **kwargs)
    return s3_client.generate_presigned_url('get_object',
                                    Params={'Bucket': bucket,
                                            'Key': file_url.replace(f"s3://{bucket}/", "")},
                                    ExpiresIn=60*60*24)
    $$;
revoke all on function create_study_h5ad_presigned_url from postgraphile;
revoke all on function create_study_h5ad_presigned_url from public;

CREATE VIEW api_study_h5_download with (security_invoker = true) AS
SELECT create_study_h5ad_presigned_url(study.import_file) as presigned_url, study_id FROM study where study.visible = True;
grant select on api_study_h5_download to postgraphile;