DROP FUNCTION IF EXISTS get_correlated_genes(study_id int, omics_id int);
CREATE OR REPLACE FUNCTION get_correlated_genes(study_id int, omics_id int)
    RETURNS table
            (
                omics_id       INT,
                display_symbol TEXT,
                display_name   TEXT,
                r              FLOAT
            )
AS
$$

import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed, cpu_count
import numpy as np
import scipy
# from numba import njit
from pathlib import Path
from smart_open import open
import io


def sql_query(query):
    # postgres data retrieval with consistent output, both in the jupyter development
    # environment (plpy is not available) and at runtime inside a plpython3u stored procedure
    try:
        import plpy
    except:
        from postgres_utils import engine
        from sqlalchemy import text
        with engine.connect() as connection:
            r = connection.execute(text(query))
            return [row._mapping for row in r.fetchall()]
    r = plpy.execute(query)
    return [row for row in r]


def h5ad_read(filename):
    if filename.startswith('s3:'):
        # AnnData's read_h5ad passes the "filename" parameter to h5py.File, which supports file-like objects in
        # addition to filename strings. It is able to read an AnnData file directly from S3 using the python
        # file-like object abstraction smart_open provides, however it seeks a lot and that causes read performance
        # to drop significantly. So we're copying the h5ad file into an in-memory file and read from there.
        s3_file_like_obj = open(filename, 'rb')
        memory_file_like_obj = io.BytesIO(s3_file_like_obj.read())
        s3_file_like_obj.close()
        adata = sc.read_h5ad(memory_file_like_obj)
        memory_file_like_obj.close()
        return adata
    else:
        return sc.read_h5ad(filename)


# numba is disabled, the function crashes if the correlation result is empty, e.g.
# for KRAS in the "COVID-19" study:
# @njit
def pearson_corr(m):
    return np.corrcoef(m)[0, 1]


def compute_correlation(m, genes, goi):
    collect = []
    for gene in genes:
        r = pearson_corr(m[[goi, gene]])
        if r >= 0.2:
            collect.append({'h5ad_var_index': gene, 'r': r})
    return pd.DataFrame(collect)


data = sql_query(f"""
    SELECT s.filename,so.h5ad_var_index FROM study s
    JOIN study_omics so
      ON s.study_id=so.study_id
    WHERE omics_id={omics_id} AND s.study_id={study_id}
    """)
goi = data[0].get('h5ad_var_index')
fn = data[0].get('filename')
if not fn.startswith('s3:'):
    fn = Path('/h5ad_store') / Path(fn).name
adata = h5ad_read(fn)
df = adata.X.T.todense()

chunks = np.array_split([_ for _ in range(0, adata.n_vars) if _ != goi], 8)
result = Parallel(n_jobs=cpu_count(), backend="threading")(
    delayed(compute_correlation)(m=df, genes=chunk, goi=goi) for chunk in chunks)
result = pd.concat(result)
if len(result) == 0:
    return []

geneids = tuple(result.h5ad_var_index.tolist())
ret = sql_query(f'''
    SELECT so.omics_id, ob.display_symbol, ob.display_name, so.h5ad_var_index FROM study_omics so
      JOIN omics_base ob
        ON ob.omics_id = so.omics_id
     WHERE  so.study_id = {study_id}
       AND so.h5ad_var_index in {geneids}
''')
out = pd.DataFrame(ret).merge(result, on='h5ad_var_index').drop('h5ad_var_index', axis=1)
out = out[['omics_id', 'display_symbol', 'display_name', 'r']]

return out.sort_values('r', ascending=False).to_records(index=False)
$$ LANGUAGE plpython3u
    IMMUTABLE
    SECURITY INVOKER -- secured by study_policy, as we're retrieving the h5ad filename from the study table
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
                                        p_subsampling_projection text)
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
    if p_subsampling_projection is not null and p_subsampling_projection != '' then
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


-- box plot / dot plot calculated in the database:

drop type if exists expression_by_annotation cascade;
create type expression_by_annotation as
(
    study_layer_id           int,
    omics_id                 int,
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

drop function if exists expression_by_annotation;
create function expression_by_annotation(p_study_layer_ids int[], p_omics_ids int[],
                                         p_annotation_group_id int, p_exclude_annotation_value_ids int[])
    returns setof expression_by_annotation
    language sql
    stable
as
$$
with expr as (select e.study_layer_id, sample_id, e.omics_id, value
              from expression e
                       cross join unnest(e.study_sample_ids, e.values) as x(sample_id, value)
              where omics_id = any (p_omics_ids)
                and study_layer_id = any (p_study_layer_ids)
                -- ordering both CTEs by sample_id is important for the join on "expr.sample_id = annot.sample_id" below
              order by 1, 2, 3),
     annot as (select sl.study_layer_id,
                      sample_id,
                      ssa.annotation_value_id,
                      count(1) over (partition by sl.study_layer_id, ssa.annotation_value_id ) sample_count
               from study_sample_annotation ssa
                        cross join unnest(ssa.study_sample_ids) sample_id
                        join study_layer sl on ssa.study_id = sl.study_id
               where sl.study_layer_id = any (p_study_layer_ids)
                 and ssa.annotation_value_id in (select annotation_value_id
                                                 from annotation_value
                                                 where annotation_group_id = p_annotation_group_id)
                 and sample_id not in (select exclude_sample_id
                                       from study_sample_annotation exclude_ssa
                                                cross join unnest(exclude_ssa.study_sample_ids) exclude_sample_id
                                       where exclude_ssa.annotation_value_id = any (p_exclude_annotation_value_ids))
               order by 1, 2)
select expr.study_layer_id,
       expr.omics_id,
       annot.annotation_value_id,
       (select av.display_value
        from annotation_value av
        where av.annotation_value_id = annot.annotation_value_id) annotation_display_value,
       array_agg(value)                                           values,
       boxplot(value)                                             boxplot_params,
       percentile_cont(0.5) within group (order by value)  as     median,
       percentile_cont(0.75) within group (order by value) as     q3,
       avg(value)                                          as     "mean",
       count(1)                                                   value_count,
       count(1) :: real / annot.sample_count                      expr_samples_fraction
from expr
         join annot on expr.sample_id = annot.sample_id and expr.study_layer_id = annot.study_layer_id
group by expr.study_layer_id, expr.omics_id, annot.annotation_value_id, annot.sample_count
$$;

/*
select study_layer_id, omics_id, annotation_value_id, annotation_display_value, q3, expr_samples_fraction from expression_by_annotation(ARRAY[3],
    ARRAY[10382, 11906], 3, ARRAY[]::int[]);
*/
