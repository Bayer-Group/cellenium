create table ontology
(
    ontid int primary key,
    name  text
);

create table concept
(
    cid            serial primary key,
    ontid          int references ontology,
    ont_code       text,
    label          text,
    label_tsvector tsvector generated always as ( to_tsvector('english', label) ) stored
);
create unique index concept_i1 on concept (ontid, ont_code);
create index concept_i2 on concept (lower(label), ontid, cid);

create table concept_synonym
(
    cid              int not null references concept,
    synonym          text,
    synonym_tsvector tsvector generated always as ( to_tsvector('english', synonym) ) stored
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
