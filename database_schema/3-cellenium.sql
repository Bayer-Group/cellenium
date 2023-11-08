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
    cell_count         int,
    projections        text[],
    visible            boolean default False,
    import_started     boolean default False,
    import_finished    boolean default False,
    import_failed      boolean default False,
    import_log         text,
    import_file        text, --- s3 file path
    reader_permissions text[],
    admin_permissions  text[],
    metadata           jsonb,
    legacy_config      jsonb
);



ALTER TABLE study
    ENABLE ROW LEVEL SECURITY;

CREATE POLICY study_policy ON study FOR SELECT TO postgraphile
    USING (
        study_id in (select study_id
                     from study_visible_currentuser)
    );
CREATE POLICY study_update_policy ON study FOR UPDATE TO postgraphile
    USING (
        study_id in (select study_id
                     from study_administrable_currentuser)
    );
CREATE POLICY study_delete_policy ON study FOR DELETE TO postgraphile
    USING (
        study_id in (select study_id
                     from study_administrable_currentuser)
    );

DROP POLICY IF EXISTS study_insert_policy ON study;
CREATE POLICY study_insert_policy ON study FOR INSERT TO postgraphile
    WITH CHECK (
    true
    );

grant insert, select, update, delete on study to postgraphile;
grant usage, select ON sequence study_study_id_seq TO postgraphile;

CREATE OR REPLACE FUNCTION create_study_for_current_user(in study_name text)
    RETURNS bool
    LANGUAGE plpgsql
AS
$$
BEGIN
    INSERT INTO study (study_name, reader_permissions, admin_permissions)
    VALUES (study_name, null, null);
    RETURN true;
END;
$$;



drop table if exists omics_base cascade;
CREATE TYPE omics_type AS ENUM ('gene', 'protein_antibody_tag', 'transcription_factor', 'region');

CREATE TABLE omics_base
(
    omics_id                serial primary key,
    omics_type              omics_type not null,
    tax_id                  int        not null,
    display_symbol          text       not null,
    display_name            text,
    display_name_tsvector   tsvector generated always as ( to_tsvector('english', display_name) ) stored,
    display_symbol_tsvector tsvector generated always as ( to_tsvector('english', display_symbol) ) stored
);
grant select on omics_base to postgraphile;
CREATE INDEX omics_base_type_idx ON omics_base (omics_type);


CREATE TABLE omics_gene
(
    gene_id         int  not null references omics_base primary key,
    ensembl_gene_id text not null,
    entrez_gene_ids text[],
    hgnc_symbols    text[]
);
grant select on omics_gene to postgraphile;
CREATE UNIQUE INDEX omics_gene_1 on omics_gene (ensembl_gene_id);

DROP TABLE IF EXISTS omics_region CASCADE;
CREATE TABLE omics_region
(
    region_id      int  not null references omics_base primary key,
    chromosome     text NOT NULL,
    start_position int  NOT NULL,
    end_position   int  NOT NULL,
    region         text NOT NULL
);
grant select on omics_region to postgraphile;

DROP INDEX IF EXISTS omics_region_1;
CREATE UNIQUE INDEx omics_region_1 on omics_region (region);

DROP TABLE IF EXISTS omics_region_gene;
CREATE TABLE omics_region_gene
(
    region_id int not null references omics_region,
    gene_id   int not null references omics_gene
);
grant select on omics_region_gene to postgraphile;
DROP INDEX IF EXISTS omics_region_gene_1;
CREATE UNIQUE INDEX omics_region_gene_1 on omics_region_gene (region_id, gene_id);

COMMENT ON CONSTRAINT "omics_region_region_id_fkey" ON "public"."omics_region" IS E'@fieldName omics_region_gene_region';
alter table omics_region_gene
    alter constraint omics_region_gene_gene_id_fkey DEFERRABLE INITIALLY IMMEDIATE;


-- insert into omics_base (omics_id,omics_type,tax_id,display_symbol,display_name) values (100000, 'region', 9606, 'chr1:120-125', 'chr1:120-125');
-- insert into omics_region (region_id, region) values (100000, 'chr1:120-125');
-- insert into omics_region_gene (region_id, gene_id) values (100000, 1);
-- insert into omics_region_gene (region_id, gene_id) values (100000, 2);
-- insert into omics_region_gene (region_id, gene_id) values (100000, 3);


-- cite-seq
CREATE TABLE omics_protein_antibody_tag
(
    protein_antibody_tag_id int  not null references omics_base primary key,
    tax_id                  int  not null,
    antibody_symbol         text not null
);
grant select on omics_protein_antibody_tag to postgraphile;
create unique index omics_protein_antibody_tag_1 on omics_protein_antibody_tag (tax_id, antibody_symbol);


CREATE TABLE omics_protein_antibody_tag_gene
(
    protein_antibody_tag_id int not null references omics_protein_antibody_tag,
    gene_id                 int not null references omics_gene
);
alter table omics_protein_antibody_tag_gene
    alter constraint omics_protein_antibody_tag_gene_gene_id_fkey DEFERRABLE INITIALLY IMMEDIATE;
grant select on omics_protein_antibody_tag_gene to postgraphile;
create unique index omics_protein_antibody_tag_gene_1 on omics_protein_antibody_tag_gene (protein_antibody_tag_id, gene_id);

CREATE TABLE omics_transcription_factor
(
    omics_id         int  not null references omics_base primary key,
    jaspar_matrix_id text not null
);
grant select on omics_transcription_factor to postgraphile;
create unique index omics_transcription_factor_1 on omics_transcription_factor (jaspar_matrix_id);


CREATE TABLE omics_transcription_factor_gene
(
    transcription_factor_id int not null references omics_transcription_factor,
    gene_id                 int not null references omics_gene
);
alter table omics_transcription_factor_gene
    alter constraint omics_transcription_factor_gene_gene_id_fkey DEFERRABLE INITIALLY IMMEDIATE;
grant select on omics_transcription_factor_gene to postgraphile;
create unique index omics_transcription_factor_gene_1 on omics_transcription_factor_gene (transcription_factor_id, gene_id);

-- insert into omics_base (omics_id, omics_type, tax_id, display_symbol) values (1000000, 'gene', 9606, 'ACP5');
-- insert into omics_gene (gene_id, ensembl_gene_id, hgnc_symbols) values (1000000, 'ENSG00000102575', ARRAY ['ACP5']);


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
grant select on annotation_group to postgraphile;
create unique index annotation_group_1 on annotation_group (h5ad_column);

-- e.g. an annotation category value, like 'lymphocyte'
CREATE TABLE annotation_value
(
    annotation_value_id serial primary key,
    annotation_group_id int  not null references annotation_group,

    h5ad_value          text not null,
    display_value       text not null
);
grant select on annotation_value to postgraphile;
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
grant select on study_annotation_group_ui to postgraphile;
create unique index study_annotation_group_ui_1 on study_annotation_group_ui (study_id, annotation_group_id);

CREATE TABLE study_sample
(
    study_id        int  not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,

    study_sample_id int  not null,
    constraint pk_study_sample primary key (study_id, study_sample_id),

    h5ad_obs_index  int  not null,
    h5ad_obs_key    text not null
);
grant select on study_sample to postgraphile;


CREATE TABLE study_sample_projection
(
    study_id            int     not null,
    study_sample_id     int     not null,
    constraint fk_study_sample
        FOREIGN KEY (study_id, study_sample_id)
            REFERENCES study_sample (study_id, study_sample_id),
    projection_type     text    not null,
    modality            text,
    projection          real[]  not null,
    -- subsampling reduces overlapping points in a projection
    display_subsampling boolean not null
);
grant select on study_sample_projection to postgraphile;
create index study_sample_projection_1 on study_sample_projection (study_id, projection_type, display_subsampling) include (study_sample_id, projection) where (display_subsampling);
create index study_sample_projection_2 on study_sample_projection (study_id, projection_type);
create unique index study_sample_projection_3 on study_sample_projection (study_id, study_sample_id, projection_type, modality);



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
grant select on study_sample_annotation to postgraphile;
create unique index study_sample_annotation_1 on study_sample_annotation (study_id, annotation_value_id);


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
    -- TODO: AS --> should actually be fine as rna and atac have different omics_id
    -- region as seen in the actual study data before 'fuzzy' region matching with bedtools (expect same build, chromosome)
    region_start   int,
    region_end     int
);
grant select on study_omics to postgraphile;
create unique index study_omics_i1 on study_omics (study_id, omics_id);



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

ALTER TABLE differential_expression
    ENABLE ROW LEVEL SECURITY;
CREATE POLICY differential_expression_policy ON differential_expression FOR SELECT TO postgraphile
    USING (
        study_id in (select study_id
                     from study_visible_currentuser)
    );
grant select on differential_expression to postgraphile;


CREATE TABLE study_layer
(
    study_layer_id serial primary key,
    study_id       int  not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,
    omics_type     omics_type,
    layer          text not null
);
grant select on study_layer to postgraphile;
create unique index study_layer_ui1 on study_layer (study_id, layer, omics_type);


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

ALTER TABLE expression
    ENABLE ROW LEVEL SECURITY;
CREATE POLICY expression_policy ON expression FOR SELECT TO postgraphile
    USING (
        study_layer_id in (select sl.study_layer_id
                           from study_visible_currentuser
                                    join study_layer sl on study_visible_currentuser.study_id = sl.study_id)
    );
grant select on expression to postgraphile;


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
END;
$$;

create or replace function delete_all_study_data(p_study_id in int)
    returns boolean
    language plpgsql
    security definer
    volatile
as
$$
declare
    delete_study_layer_id int;
begin

    for delete_study_layer_id in (SELECT study_layer_id FROM study_layer WHERE study_id = p_study_id)
        loop
            EXECUTE 'drop table if exists expression_' || delete_study_layer_id;
        end loop;
    delete from differential_expression where study_id = p_study_id;
    delete from study_omics where study_id = p_study_id;
    delete from study_layer where study_id = p_study_id;
    delete from study_sample_annotation where study_id = p_study_id;
    delete from study_annotation_group_ui where study_id = p_study_id;
    delete from study_sample_projection where study_id = p_study_id;
    delete from study_sample where study_id = p_study_id;
    delete
    from annotation_value
    where annotation_group_id in
          (select saved_as_annotation_group_id from user_annotation_group where study_id = p_study_id);
    delete from user_annotation_group where study_id = p_study_id;
    delete
    from annotation_group
    where annotation_group_id in
          (select saved_as_annotation_group_id from user_annotation_group where study_id = p_study_id);
    delete from study where study_id = p_study_id;
    return true;
end;
$$;

create or replace procedure reset_study(p_study_id in int)
    language plpgsql
as
$$
declare
    keep_study_name        text;
    keep_filename          text;
    keep_import_file       text;
    keep_admin_permissions text[];
begin
    select study_name, filename, import_file, admin_permissions
    into keep_study_name, keep_filename, keep_import_file, keep_admin_permissions
    from study
    where study_id = p_study_id;
    perform delete_all_study_data(p_study_id);
    insert into study (study_id, study_name, filename, import_file, admin_permissions)
    values (p_study_id, keep_study_name, keep_filename, keep_import_file, keep_admin_permissions);
end;
$$;


drop type if exists omics_autocomplete_result cascade;
create type omics_autocomplete_result as
(
    omics_id        integer[],
    omics_type      omics_type[],
    display_symbol  text,
    label_highlight text
);

create aggregate product(integer)
   (stype=float4,sfunc=float4mul,initcond=1);

drop function if exists omics_autocomplete(search_query text, tax_id_filter integer, omics_type_filter omics_type);
create function omics_autocomplete(search_query text, tax_id_filter integer, omics_type_filter omics_type)
    returns setof omics_autocomplete_result
    language sql
    stable
as
$$
with all_concepts_terms as (select ob.display_symbol, array_agg(DISTINCT ob.omics_id) omics_id, array_agg(DISTINCT ob.omics_type) omics_type, array_agg(DISTINCT ob.display_symbol_tsvector) display_symbol_tsvector
                            from omics_base ob
                            where ob.tax_id=tax_id_filter and ob.omics_type=omics_type_filter
                            group by ob.display_symbol
                            ),
     as_tsquery as (
         -- users may double-quote words to find them exactly, otherwise we assume a prefix search
         -- it is also possible to double-quote a phrase of multiple words
         select --to_tsquery('english','cdk | cdk:*') q
                to_tsquery('english', string_agg(
                        case
                            when right(split, 1) = '"' then replace(replace(split, '"', ''), ' ', ' <-> ')
                            else split || ' | ' || split || ':*' end,
                        ' & ')) q
         -- see https://stackoverflow.com/questions/4780728/regex-split-string-preserving-quotes
         from regexp_split_to_table(lower(regexp_replace(search_query, '[^\w|\-|\"]', ' ', 'g')),
                                    E'(?<=^[^"]*(?:"[^"]*\"[^"]*)*) (?=(?:[^"]*"[^"]*")*[^"]*$)') split),
     search_results as (select as_tsquery.q,
                               c.omics_id,
                               c.omics_type,
                               ts_rank(c.display_symbol_tsvector[1], as_tsquery.q, 2) *
                               case
                                   when lower(c.display_symbol) =
                                        lower(regexp_replace(search_query, '[^\w|\-|\"]', ' ', 'g'))
                                       then 2 /*still need to boost the exact match...*/
                                   else 1
                                   end                                     "rank",
                               c.display_symbol,
                               ts_headline(c.display_symbol, as_tsquery.q) label_highlight
                        from all_concepts_terms c,
                             as_tsquery
                        where c.display_symbol_tsvector[1] @@ as_tsquery.q)
select sro.omics_id,
       sro.omics_type,
       sro.display_symbol,
       sro.label_highlight
from search_results sro
order by sro.rank desc
$$;
