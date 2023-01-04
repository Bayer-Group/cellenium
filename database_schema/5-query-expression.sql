drop type if exists expression_by_omics cascade;
create type expression_by_omics as
(
    omics_id         int,
    study_sample_ids integer[],
    values           real[]
);

drop function if exists expression_by_omics_ids;
create function expression_by_omics_ids(p_study_layer_id int, p_omics_ids int[], p_subsampling boolean)
    returns setof expression_by_omics
    language plpgsql
    immutable
    leakproof
    parallel safe
as
$$
begin
    -- we could use a study-level flag that determines if display_subsampling is False for any sample,
    -- and use the else branch in this case
    if p_subsampling then
        return query select e.omics_id,
                            array_agg(sampleid_v order by val_i) study_sample_ids,
                            array_agg(val_value order by val_i)  values
                     from expression e
                              join study_layer sl on sl.study_layer_id = e.study_layer_id
                              cross join unnest(e.values) with ordinality as val(val_value, val_i)
                              cross join unnest(e.study_sample_ids) with ordinality as sampleid(sampleid_v, sampleid_i)
                              join study_sample s on s.study_id = sl.study_id and s.study_sample_id = val_i
                     where val_i = sampleid_i
                       and e.study_layer_id = p_study_layer_id
                       and e.omics_id = any (p_omics_ids)
                       and s.display_subsampling = True
                     group by e.study_layer_id, e.omics_id;
    else
        return query select e.omics_id,
                            e.study_sample_ids,
                            e.values
                     from expression e
                     where e.study_layer_id = p_study_layer_id
                       and e.omics_id = any (p_omics_ids);
    end if;

end;
$$;


/*
select * from expression_by_omics_ids(1, array [1], False);
select * from expression_by_omics_ids(1, array [1], True);
*/


DROP AGGREGATE IF EXISTS boxplot(real) CASCADE;
DROP FUNCTION IF EXISTS _final_boxplot(a real[]);
DROP TYPE IF EXISTS boxplot_values;

CREATE TYPE boxplot_values AS
(
    q1_whisker real,
    q1         real,
    median     real,
    q3         real,
    q3_whisker real,
    outliers   real[],
    n          integer
);

CREATE OR REPLACE FUNCTION _final_boxplot(a real[])
    RETURNS boxplot_values AS
$$
with data as (select e.e
              from unnest(a) as e
              order by e.e),
     quartiles as (select count(1)                                          as n,
                          percentile_cont(0.25) within group (order by e.e) as q1,
                          percentile_cont(0.5) within group (order by e.e)  as median,
                          percentile_cont(0.75) within group (order by e.e) as q3
                   from data AS e),
     upper_whisker as (select max(e.e) whisker
                       from data e,
                            quartiles q
                       where e between q.q3 and q.q3 + 1.5 * (q.q3 - q.q1)),
     lower_whisker as (select min(e.e) whisker
                       from data e,
                            quartiles q
                       where e between q.q1 - 1.5 * (q.q3 - q.q1) and q.q1)
select coalesce(l.whisker::real, q.q1::real),
       q.q1::real,
       q.median::real,
       q.q3::real,
       coalesce(u.whisker::real, q.q3::real),
       outliers.outlierlist::real[],
       q.n::integer
from quartiles q,
     lower_whisker l,
     upper_whisker u,
     (select coalesce(array_agg(e.e), ARRAY []::real[]) outlierlist
      from data e,
           upper_whisker u,
           lower_whisker l
      where e.e > u.whisker
         or e.e < l.whisker) outliers
$$
    LANGUAGE sql IMMUTABLE;

CREATE AGGREGATE boxplot(real) (
    SFUNC = array_append,
    STYPE = real[],
    INITCOND = '{}',
    FINALFUNC = _final_boxplot
    );


-- box plot calculated in the database:
drop view if exists expression_by_annotation_boxplot;
create view expression_by_annotation_boxplot
as
with sample_annotation as (select ssa.study_id, av.annotation_group_id, ssa.annotation_value_id, study_sample_id
                           from study_sample_annotation ssa
                                    join annotation_value av on ssa.annotation_value_id = av.annotation_value_id
                                    cross join unnest(ssa.study_sample_ids) as study_sample_id)
select e.study_layer_id,
       e.omics_id,
       sa.annotation_group_id,
       sa.annotation_value_id,
       array_agg(val_value) values,
       boxplot(val_value)
from expression e
         join study_layer sl on sl.study_layer_id = e.study_layer_id
         cross join unnest(e.values) with ordinality as val(val_value, val_i)
         cross join unnest(e.study_sample_ids) with ordinality as sampleid(sampleid_v, sampleid_i)
         join sample_annotation sa on sa.study_id = sl.study_id and sa.study_sample_id = sampleid_v
where val_i = sampleid_i
group by e.study_layer_id, e.omics_id, sa.annotation_group_id, sa.annotation_value_id;

-- select * from expression_by_annotation_boxplot where study_layer_id = 1 and omics_id = 116 and annotation_group_id = 1;
