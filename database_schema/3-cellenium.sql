CREATE TABLE study
(
    filename           text,
    study_id           serial primary key,
    study_name         text not null,
    description        text,
    external_website   text,
    tissue_ncit_ids    text[],
    disease_mesh_ids   text[],
    cell_ontology_ids  text[],
    organism_tax_id    text,
    --organism_label          text,
    --ontology_ids_with_parents text[]   for search

    cell_count         int,
    projections        text[],
    visible            boolean default False,
    reader_permissions text[],
    admin_permissions  text[],
    legacy_config      jsonb
);


drop table if exists omics_base cascade;
CREATE TYPE omics_type AS ENUM ('gene', 'protein_antibody_tag', 'transcription_factor', 'region');

CREATE TABLE omics_base
(
    omics_id       serial primary key,
    omics_type     omics_type not null,
    tax_id         int        not null,
    display_symbol text       not null,
    display_name   text
);

CREATE TABLE omics_gene
(
    gene_id         int  not null references omics_base primary key,
    ensembl_gene_id text not null,
    entrez_gene_ids text[],
    hgnc_symbols    text[]
);
create unique index omics_gene_1 on omics_gene (ensembl_gene_id);

CREATE TABLE omics_genome_range
(
    genomerange_id not null references omics_base primary key,
    coordinate text NOT NULL
)
create unique index omics_genome_range on omics_genome_range (coordinate);

CREATE TABLE omics_genome_range_gene
(
    genomerange_id int not null references omics_genome_range,
    gene_id                 int not null references omics_gene
);
create unique index omics_genome_range_gene_1 on omics_genome_range_gene (genomerange_id, gene_id);


-- cite-seq
CREATE TABLE omics_protein_antibody_tag
(
    protein_antibody_tag_id int  not null references omics_base primary key,
    tax_id                  int  not null,
    antibody_symbol         text not null
);
create unique index omics_protein_antibody_tag_1 on omics_protein_antibody_tag (tax_id, antibody_symbol);


CREATE TABLE omics_protein_antibody_tag_gene
(
    protein_antibody_tag_id int not null references omics_protein_antibody_tag,
    gene_id                 int not null references omics_gene
);
create unique index omics_protein_antibody_tag_gene_1 on omics_protein_antibody_tag_gene (protein_antibody_tag_id, gene_id);

CREATE TABLE omics_transcription_factor
(
    omics_id         int  not null references omics_base primary key,
    jaspar_matrix_id text not null
);
create unique index omics_transcription_factor_1 on omics_transcription_factor (jaspar_matrix_id);


CREATE TABLE omics_transcription_factor_gene
(
    transcription_factor_id int not null references omics_transcription_factor,
    gene_id                 int not null references omics_gene
);
create unique index omics_transcription_factor_gene_1 on omics_transcription_factor_gene (transcription_factor_id, gene_id);

-- insert into omics_base (omics_id, omics_type, tax_id, display_symbol) values (1000000, 'gene', 9606, 'ACP5');
-- insert into omics_gene (gene_id, ensembl_gene_id, hgnc_symbols) values (1000000, 'ENSG00000102575', ARRAY ['ACP5']);


create view omics_all as
select b.omics_id,
       b.omics_type,
       b.tax_id,
       b.display_symbol,
       b.display_name,
       og.ensembl_gene_id,
       og.entrez_gene_ids,
       og.hgnc_symbols,
       coalesce(
               array_agg(opatg.gene_id),
               array_agg(otfg.gene_id)
           ) linked_genes
from omics_base b
         left join omics_gene og on b.omics_id = og.gene_id
         left join omics_protein_antibody_tag_gene opatg on b.omics_id = opatg.protein_antibody_tag_id
         left join omics_transcription_factor_gene otfg on b.omics_id = otfg.transcription_factor_id
group by b.omics_id, b.omics_type, b.tax_id, b.display_symbol, b.display_name,
         og.ensembl_gene_id, og.entrez_gene_ids, og.hgnc_symbols;

/*
-- TODO add omics_region... tables, same style

    build            text,
    region_chr       text,
    region_start     int,
    region_end       int
);
--create unique index omics_element_uq_region on omics (tax_id, build, region_chr, region_start, region_end);
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

 */

-- e.g. an annotation category, like 'cell ontology name'
CREATE TABLE annotation_group
(
    annotation_group_id serial primary key,
    h5ad_column         text not null,
    display_group       text not null
);
create unique index annotation_group_1 on annotation_group (h5ad_column);

-- e.g. an annotation category value, like 'lymphocyte'
CREATE TABLE annotation_value
(
    annotation_value_id serial primary key,
    annotation_group_id int  not null references annotation_group,

    h5ad_value          text not null,
    display_value       text not null
);
create unique index annotation_value_1 on annotation_value (annotation_group_id, h5ad_value);
create index annotation_value_2 on annotation_value (annotation_group_id) include (annotation_value_id);
create index annotation_value_3 on annotation_value (annotation_value_id) include (display_value, annotation_group_id);

CREATE TABLE study_annotation_group_ui
(
    study_id                           int     not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,

    annotation_group_id                int     not null references annotation_group,

    is_primary                         boolean not null,
    ordering                           int     not null,
    differential_expression_calculated boolean not null
);

CREATE TABLE study_sample
(
    study_id        int not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,

    study_sample_id int not null,
    constraint pk_study_sample primary key (study_id, study_sample_id),

    h5ad_obs_index  int not null
);
--create unique index study_sample_i1 on study_sample (study_id, study_sample_id);


CREATE TABLE study_sample_projection
(
    study_id            int     not null,
    study_sample_id     int     not null,
    constraint fk_study_sample
        FOREIGN KEY (study_id, study_sample_id)
            REFERENCES study_sample (study_id, study_sample_id) ON DELETE CASCADE,
    projection_type     text    not null,
    projection          real[]  not null,
    -- subsampling reduces overlapping points in a projection
    display_subsampling boolean not null
);

CREATE VIEW study_sample_projection_subsampling_transposed
as
select study_id,
       projection_type,
       array_agg(study_sample_id order by study_sample_id) study_sample_id,
       array_agg(projection order by study_sample_id)      projection
from study_sample_projection
where display_subsampling = True
group by study_id, projection_type;
comment on view study_sample_projection_subsampling_transposed is
    E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName studySampleProjectionSubsamplingTransposed';

CREATE TABLE study_sample_annotation
(
    study_id            int   not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,

    annotation_value_id int   not null,
    constraint fk_sample_annotation_value
        FOREIGN KEY (annotation_value_id)
            REFERENCES annotation_value (annotation_value_id),

    -- the samples that are annotated with that value, e.g. that specific cell type
    study_sample_ids    int[] not null,

    color               text
);
create unique index study_sample_annotation_1 on study_sample_annotation (study_id, annotation_value_id);

CREATE VIEW study_sample_annotation_subsampling
as
select ssa.study_id,
       ssa.annotation_value_id,
       array_agg(ssp.study_sample_id) study_sample_ids
from study_sample_annotation ssa
         cross join unnest(ssa.study_sample_ids) sample_id
         join study_sample_projection ssp on ssp.study_id = ssa.study_id and ssp.study_sample_id = sample_id
where ssp.display_subsampling = True
group by ssa.study_id, ssa.annotation_value_id;
comment on view study_sample_annotation_subsampling is
    E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName studySampleAnnotationSubsampling';

CREATE TABLE study_omics
(
    study_id       int not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,

    omics_id       int not null references omics_base,

    -- indexing the h5ad .uns['protein_X'] matrix in this study
    h5ad_var_index int not null,
    -- TODO add another h5ad_col_index for second h5ad file (ATAC-seq)? Or better use h5ad format to combine atac-seq into same h5ad file

    -- region as seen in the actual study data before 'fuzzy' region matching with bedtools (expect same build, chromosome)
    region_start   int,
    region_end     int
);
create unique index study_omics_i1 on study_omics (study_id, omics_id);

CREATE VIEW study_omics_transposed
as
select study_id,
       array_agg(ob.omics_id order by ob.omics_id)       omics_id,
       array_agg(ob.omics_type order by ob.omics_id)     omics_type,
       array_agg(ob.display_symbol order by ob.omics_id) display_symbol,
       array_agg(ob.display_name order by ob.omics_id)   display_name
from study_omics
         join omics_base ob on study_omics.omics_id = ob.omics_id
group by study_id;
comment on view study_omics_transposed is
    E'@foreignKey (study_id) references study (study_id)|@fieldName study|@foreignFieldName studyOmicsTransposed';


CREATE TABLE differential_expression
(
    study_id            int not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,
    omics_id            int not null references omics_base,

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

CREATE VIEW differential_expression_v
    with (security_invoker = true)
AS
SELECT de.*, ob.display_symbol, ob.display_name
FROM differential_expression de
         JOIN omics_base ob on de.omics_id = ob.omics_id;

CREATE TABLE study_layer
(
    study_layer_id serial primary key,
    study_id       int        not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,
    omics_type     omics_type not null,
    layer          text       not null

);
create unique index study_layer_ui1 on study_layer (study_id, layer);


CREATE TABLE expression
(
    study_layer_id   int       not null,
    constraint fk_study_layer_id
        FOREIGN KEY (study_layer_id)
            REFERENCES study_layer (study_layer_id) ON DELETE CASCADE,
    omics_id         int       not null references omics_base,

    -- for sparse data, references study_sample.study_sample_id
    study_sample_ids integer[] not null,
    values           real[]    not null

) partition by list (study_layer_id);


CREATE OR REPLACE PROCEDURE add_studylayer_partition(study_layer_id int)
    LANGUAGE plpgsql
AS
$$
BEGIN
    EXECUTE format(
            'create table expression_%1$s
            partition of expression
            (
                study_layer_id,
                omics_id,
                study_sample_ids,
                values
                )
            for values in ( %1$s );
            comment on table expression_%1$s is ''@omit'';
            create unique index  expression_%1$s_omics_uq on expression_%1$s(omics_id);
            ', study_layer_id);
END ;
$$;
