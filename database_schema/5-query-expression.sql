DROP FUNCTION IF EXISTS get_correlated_genes(study_id int, omics_id int);
CREATE OR REPLACE FUNCTION get_correlated_genes(study_id int, omics_id int)
    RETURNS table
        (
            omics_id INT,
            display_symbol TEXT,
            display_name TEXT,
            r FLOAT
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

# numba is disabled, the function crashes if the correlation result is empty, e.g.
# for KRAS in the "COVID-19" study:
#@njit
def pearson_corr(m):
    return np.corrcoef(m)[0,1]

def compute_correlation(m, genes, goi):
    collect = []
    for gene in genes:
        r = pearson_corr(m[[goi,gene]])
        if r>=0.2:
            collect.append({'h5ad_var_index': gene, 'r':r})
    return pd.DataFrame(collect)


data = sql_query(f"""
    SELECT s.filename,so.h5ad_var_index FROM study s
    JOIN study_omics so
      ON s.study_id=so.study_id
    WHERE omics_id={omics_id} AND s.study_id={study_id}
    """)
fn = Path('/h5ad_store') / Path(data[0].get('filename')).name
goi = data[0].get('h5ad_var_index')

adata = sc.read(fn)
df = adata.X.T.todense()

chunks = np.array_split([_ for _ in range(0,adata.n_vars) if _!=goi],8)
result = Parallel(n_jobs=cpu_count(), backend="threading")(delayed(compute_correlation)(m = df, genes=chunk, goi=goi) for chunk in chunks)
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
out = pd.DataFrame(ret).merge(result, on = 'h5ad_var_index').drop('h5ad_var_index', axis =1)
out = out[['omics_id','display_symbol','display_name','r']]

return out.sort_values('r', ascending = False).to_records(index=False)
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


-- box plot / dot plot calculated in the database:
drop view if exists expression_by_annotation cascade;
create view expression_by_annotation
    with (security_invoker = true)
as
select sl.study_id,
       e.study_layer_id,
       e.omics_id,
       av.annotation_group_id,
       av.annotation_value_id,
       av.display_value                                       annotation_display_value,
       array_agg(value)                                       values,
       boxplot(value)                                         boxplot_params,
       percentile_cont(0.75) within group (order by value) as q3,
       count(1) :: real / cardinality(ssa.study_sample_ids)   expr_cells_fraction
from expression e
         join study_layer sl on sl.study_layer_id = e.study_layer_id
         cross join unnest(e.study_sample_ids, e.values) as x(sample_id, value)
         cross join annotation_value av
         join study_sample_annotation ssa
              on ssa.study_id = sl.study_id and ssa.annotation_value_id = av.annotation_value_id
where x.sample_id = ANY (ssa.study_sample_ids)
group by sl.study_id, e.study_layer_id, e.omics_id, av.annotation_group_id, av.annotation_value_id, av.display_value,
         cardinality(ssa.study_sample_ids);

/*
select omics_id, annotation_display_value, q3, expr_cells_fraction
from expression_by_annotation
where study_id = 5
  and study_layer_id=5 and annotation_group_id = 12
  and omics_id in (10382, 11906);
*/
