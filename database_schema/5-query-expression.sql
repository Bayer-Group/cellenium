DROP FUNCTION IF EXISTS get_correlated_genes(study_id int, omics_id int);
CREATE OR REPLACE FUNCTION get_correlated_genes(study_id int, omics_id int)
    RETURNS table
        (
            omics_id INT,
            pearson FLOAT
        )
AS
$$
import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed, cpu_count
import numpy as np
import scipy
from numba import njit

@njit
def pearson_corr(v1,v2):
    return np.corrcoef(v1,v2)

def compute_correlation(df, genes, goi):
    collect = []
    for gene in genes:
        tmp = df[[goi,gene]]
        tmp_both = tmp.loc[~(tmp==0).all(axis=1)].to_numpy()
        try:
            r,p = scipy.stats.pearsonr(tmp_both[:,0],tmp_both[:,1])
            if abs(r)>=0.2:
                collect.append({'gene': gene, 'r':r})
        except:
            pass
    return pd.DataFrame(collect)
fn = '../scratch/blood_covid.h5ad'
adata = sc.read(fn)
df = adata.to_df()
chunks = np.array_split(df.columns.drop(goi),8)
result = Parallel(n_jobs=cpu_count())(delayed(compute_correlation)(df = df,genes=chunk, goi=goi) for chunk in chunks)
pd.concat(result).sort_values('r', ascending = False).reset_index(drop = True)
$$ LANGUAGE plpython3u
    IMMUTABLE
    SECURITY DEFINER
    PARALLEL SAFE;


drop type if exists expression_by_omics cascade;
create type expression_by_omics as
(
    omics_id         int,
    study_sample_ids integer[],
    values           real[]
);

drop function if exists expression_by_omics_ids;
create function expression_by_omics_ids(p_study_layer_id int, p_omics_ids int[],
                                        p_subsampling_projection projection_type)
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
    if p_subsampling_projection is not null then
        return query select e.omics_id,
                            array_agg(sample_id order by sample_id) study_sample_ids,
                            array_agg(value order by sample_id)     values
                     from expression e
                              join study_layer sl on sl.study_layer_id = e.study_layer_id
                              cross join unnest(e.study_sample_ids, e.values) as x(sample_id, value)
                              join study_sample_projection sp
                                   on sp.study_id = sl.study_id and sp.study_sample_id = sample_id and
                                      sp.projection_type = p_subsampling_projection and sp.display_subsampling = True
                     where e.study_layer_id = p_study_layer_id
                       and e.omics_id = any (p_omics_ids)
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
select * from expression_by_omics_ids(1, array [1,116], null);
select * from expression_by_omics_ids(1, array [1,116], 'umap');
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
drop view if exists expression_by_annotation cascade;
create view expression_by_annotation
as
with sample_annotation as (select ssa.study_id, av.annotation_group_id, ssa.annotation_value_id, study_sample_id
                           from study_sample_annotation ssa
                                    join annotation_value av on ssa.annotation_value_id = av.annotation_value_id
                                    cross join unnest(ssa.study_sample_ids) as study_sample_id)
select sl.study_id,
       e.study_layer_id,
       e.omics_id,
       sa.annotation_group_id,
       sa.annotation_value_id,
       array_agg(value)                                         values,
       boxplot(value)                                           boxplot_params,
       percentile_cont(0.75) within group (order by value) as   q3,
       count(1) :: real /
       (select cardinality(ssa.study_sample_ids) sample_count
        from study_sample_annotation ssa
        where sl.study_id = ssa.study_id
          and ssa.annotation_value_id = sa.annotation_value_id) expr_cells_fraction
from expression e
         join study_layer sl on sl.study_layer_id = e.study_layer_id
         cross join unnest(e.study_sample_ids, e.values) as x(sample_id, value)
    --          cross join unnest(e.values) with ordinality as val(val_value, val_i)
--          cross join unnest(e.study_sample_ids) with ordinality as sampleid(sampleid_v, sampleid_i)
         join sample_annotation sa on sa.study_id = sl.study_id and sa.study_sample_id = sample_id
group by sl.study_id, e.study_layer_id, e.omics_id, sa.annotation_group_id, sa.annotation_value_id;

-- select * from expression_by_annotation where study_layer_id = 1 and omics_id = 116 and annotation_group_id = 1;

drop view if exists expression_by_celltype;
create view expression_by_celltype
as
select study_id,
       study_layer_id,
       omics_id,
       e.annotation_group_id,
       e.annotation_value_id,
       av.display_value celltype,
       q3,
       expr_cells_fraction
from expression_by_annotation e
         join annotation_value av on e.annotation_value_id = av.annotation_value_id
where e.annotation_group_id = (select annotation_group_id from annotation_group where h5ad_column = 'CellO_celltype');

-- select * from expression_by_celltype where omics_id =8356;
