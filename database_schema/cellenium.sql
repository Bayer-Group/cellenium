create extension if not exists plpython3u;
drop table if exists ontology, concept, concept_synonym, concept_hierarchy;
drop table if exists sample_annotation, sample_annotation_value cascade;
drop table if exists study, annotation, annotation_value, study_sample_annotation_ui, study_omics,
    differential_expression, study_sample, study_sample_annotation, study_layer, expression;
drop table if exists omics, omics_region_gene;

create table ontology
(
    ontid    int primary key,
    name     text,
    umls_sab text
    --cid_range_lower int, -- unused, could be used to filter concept based on ontology, to speed things up
    --cid_range_upper int
);

create table concept
(
    cid                  serial primary key,
    ontid                int references ontology,
    ont_code             text,
    label                text,
    label_tsvector       tsvector,
    umls_disease         boolean,
    umls_semantic_types  text[],
    umls_semantic_groups text[]
);
create unique index concept_i1 on concept (ontid, ont_code);
create index concept_i2 on concept (lower(label), ontid, cid);

create table concept_synonym
(
    cid              int not null references concept,
    synonym          text,
    synonym_tsvector tsvector
);
create index concept_synonym_i1 on concept_synonym (cid);


create table concept_hierarchy
(
    cid        int references concept,
    parent_cid int references concept
);
-- alternative: store the full parent-path(s) for each cid using the ltree data type, see e.g.
-- https://hoverbear.org/blog/postgresql-hierarchical-structures/
create unique index concept_hierarchy_i1 on concept_hierarchy (cid, parent_cid);
create unique index concept_hierarchy_i2 on concept_hierarchy (parent_cid, cid);


CREATE TABLE study
(
    study_id                serial primary key,
    study_name              text not null,
    description             text,
    cluster_color_map       jsonb,

    --attributes              text[],

    cluster_hulls           jsonb,
    plot_coords             jsonb,
    h5adfile_modified_date  timestamptz,
    import_status           text,
    import_status_updated   timestamptz,
    attribute_value_freq    jsonb,
    cell_count              int,
--     qc_status               text,
--     qc_result               jsonb,

    -- for subsampling:
    projection_cell_coords  jsonb,
    projection_cell_indices jsonb
);

-- gene, region etc
CREATE TABLE omics
(
    omics_id         serial primary key,
    omics_type       text  not null,

    tax_id           int   not null,
    display_symbol   text  not null,
    display_name     text,

    -- if omics_type=='gene'
    ensembl_gene_id  text,
    entrez_gene_ids  text[],
    hgnc_symbols     text[],

    -- references omics_element, used if omics_type=='transcription_factor' or 'region' or 'CITE-seq'
    linked_genes     int[] not null default array []:: int[],

    -- if omics_type=='CITE-seq' (protein antibody tag)
    antibody_symbol  text,
    antibody_name    text,
--     gene_symbols         text[],
--     ensembl_gene_ids     text[],

    -- if omics_type=='transcription_factor' (promoter)
    jaspar_matrix_id text,

    -- if omics_type=='region'
    build            text,
    region_chr       text,
    region_start     int,
    region_end       int
);
create unique index omics_element_uq_gene on omics (tax_id, ensembl_gene_id);
create unique index omics_element_uq_antibody on omics (tax_id, antibody_symbol);
create unique index omics_element_uq_promoter on omics (tax_id, jaspar_matrix_id);
create unique index omics_element_uq_region on omics (tax_id, build, region_chr, region_start, region_end);

-- TODO remove dummy data once we have an initial import of genes
insert into omics(omics_type, tax_id, display_symbol, display_name, ensembl_gene_id, entrez_gene_ids, hgnc_symbols)
values ('gene', 9606, 'ACP5', 'acid phosphatase 5, tartrate resistant', 'ENSG00000102575', ARRAY ['54'],
        ARRAY ['ACP5']);
insert into omics(omics_type, tax_id, display_symbol, display_name, ensembl_gene_id, entrez_gene_ids, hgnc_symbols)
values ('gene', 9606, 'VCAN', 'versican', 'ENSG00000038427', ARRAY ['1462'], ARRAY ['VCAN']);



CREATE TABLE omics_region_gene
(
    omics_id        int not null,
    constraint fk_omics_element_region_index
        FOREIGN KEY (omics_id)
            REFERENCES omics (omics_id) ON DELETE CASCADE,
    gene            text,
    ensembl_gene_id text,
    evidence        text,
    evidence_score  real,
    evidence_source text
);
create unique index omics_region_uq on omics_region_gene (omics_id, gene, evidence, evidence_source);

-- e.g. in annotation category 'cell ontology name'
CREATE TABLE annotation
(
    annotation_id serial primary key,
    h5ad_column   text not null,
    display_group text not null
);
create unique index annotation_1 on annotation (h5ad_column);

-- e.g. in annotation category value 'lymphocyte'
CREATE TABLE annotation_value
(
    annotation_value_id serial primary key,
    annotation_id       int  not null,
    constraint fk_annotation
        FOREIGN KEY (annotation_id)
            REFERENCES annotation (annotation_id),

    h5ad_value          text not null,
    display_value       text not null,
    color               text
);
create unique index annotation_value_1 on annotation_value (annotation_id, h5ad_value);

CREATE TABLE study_sample_annotation_ui
(
    study_id                           int     not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,

    annotation_id                      int     not null,
    constraint fk_annotation
        FOREIGN KEY (annotation_id)
            REFERENCES annotation (annotation_id),

    is_primary                         boolean not null,
    ordering                           int     not null,
    differential_expression_calculated boolean not null
);

CREATE TABLE study_sample
(
    study_id            int     not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,

    study_sample_id     int     not null,
    constraint pk_study_sample primary key (study_id, study_sample_id),

    h5ad_obs_index      int     not null,
    display_subsampling boolean not null
);
--create unique index study_sample_i1 on study_sample (study_id, study_sample_id);

CREATE TABLE study_sample_annotation
(
    study_id            int not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,

    study_sample_id     int not null,
    constraint fk_study_sample
        FOREIGN KEY (study_id, study_sample_id)
            REFERENCES study_sample (study_id, study_sample_id) ON DELETE CASCADE,


    annotation_value_id int not null,
    constraint fk_sample_annotation_value
        FOREIGN KEY (annotation_value_id)
            REFERENCES annotation_value (annotation_value_id)
);

CREATE TABLE study_omics
(
    study_id       int not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,

    omics_id       int not null,
    constraint fk_omics_element_region_index
        FOREIGN KEY (omics_id)
            REFERENCES omics (omics_id) ON DELETE CASCADE,

    -- indexing the h5ad .uns['protein_X'] matrix in this study
    h5ad_var_index int not null,
    -- TODO add another h5ad_col_index for second h5ad file (ATAC-seq)?

    -- region as seen in the actual study data before 'fuzzy' region matching with bedtools (expect same build, chromosome)
    region_start   int,
    region_end     int
);

CREATE TABLE differential_expression
(
    study_id            int not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,
    omics_id            int,
    constraint fk_omics_element
        FOREIGN KEY (omics_id)
            REFERENCES omics (omics_id),

    /* can add this, but its redundant
    annotation_id       int not null,
    constraint fk_annotation
        FOREIGN KEY (annotation_id)
            REFERENCES annotation (annotation_id),
     */

    -- differential expression of this group (sample's annotation_value_id) vs. all other groups
    annotation_value_id int not null,
    constraint fk_sample_annotation_value
        FOREIGN KEY (annotation_value_id)
            REFERENCES annotation_value (annotation_value_id),

    pvalue              float,
    pvalue_adj          float,
    score               float,
    log2_foldchange     float
);
create unique index differential_expression_i1 on differential_expression (study_id, annotation_value_id, omics_id);

CREATE TABLE study_layer
(
    study_layer_id serial primary key,
    study_id       int  not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,
    layer          text not null
);
create unique index study_layer_ui1 on study_layer (study_id, layer);


CREATE TABLE expression
(
    study_layer_id   int       not null,
    constraint fk_study_layer_id
        FOREIGN KEY (study_layer_id)
            REFERENCES study_layer (study_layer_id) ON DELETE CASCADE,
    omics_id         int       not null,
    constraint fk_omics_element_region_index
        FOREIGN KEY (omics_id)
            REFERENCES omics (omics_id) ON DELETE CASCADE,

    study_sample_ids integer[] not null,
    values           real[]    not null

) partition by list (study_layer_id);
