CREATE TABLE user_annotation_group
(
    study_id                          int         not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,

    -- the user annotation group is persisted here in annotation_group:
    saved_as_annotation_group_id      int         not null references annotation_group,

    -- user preference regarding this annotation
    calculate_differential_expression boolean     not null,
    created_by_user                   text,
    creation_timestamp                timestamptz not null,
    private_to_user                   boolean     not null
);
create unique index user_annotation_group_1 on user_annotation_group (study_id, saved_as_annotation_group_id);
create index user_annotation_group_2 on user_annotation_group (saved_as_annotation_group_id) include (study_id);
grant select on user_annotation_group to postgraphile;

drop function if exists user_annotation_define;
create function user_annotation_define(p_study_id int,
                                       p_annotation_group_name text,
                                       p_selected_sample_ids text,
                                       p_unexpressed_samples_omics_ids int[]
)
    returns integer -- created annotation_group_id
    language plpgsql
    security definer -- need to invoke pg_background_launch, which we don't grant to the API user
    volatile
as
$$
declare
    selected_sample_ids_array   int[];
    unexpressed_sample_ids      int[];
    created_annotation_group_id int;
    created_annotation_value_id int;
    worker_pid                  int;
    study_allowed               int;
begin
    select count(1) into study_allowed from study_visible_currentuser where study_id = p_study_id;
    if study_allowed = 0 then
        raise 'no permission to create user annotation for this study';
    end if;

    -- we can't use an array parameter directly, as postgraphile seems to create a statement with more than 65k placeholders (instead of an array placeholder)
    select string_to_array(p_selected_sample_ids, ',') :: int[] into selected_sample_ids_array;

    -- In addition to specifying sample IDs for annotation, it is possible to select all samples which have no expression
    -- in either of the genes in p_unexpressed_samples_omics_ids.
    if p_unexpressed_samples_omics_ids is not null then
        select array_agg(study_sample_id)
        into unexpressed_sample_ids
        from (select study_sample_id
              from study_sample ss
              where ss.study_id = p_study_id
              EXCEPT
              select sample_id study_sample_id
              from expression e
                       cross join unnest(e.study_sample_ids) sample_id
              where e.study_layer_id = (select default_study_layer_id from study_overview where study_id = p_study_id)
                and e.omics_id = any (p_unexpressed_samples_omics_ids)) all_minus_expressed_sample_ids;
    else
        unexpressed_sample_ids := ARRAY [] :: int[];
    end if;

    insert into annotation_group (h5ad_column, display_group)
    values (p_study_id || '_' || p_annotation_group_name, p_annotation_group_name)
    returning annotation_group_id into created_annotation_group_id;

    insert into user_annotation_group (study_id, saved_as_annotation_group_id, calculate_differential_expression,
                                       created_by_user, creation_timestamp, private_to_user)
    values (p_study_id, created_annotation_group_id, false, current_user_email(), now(), true);

    insert into annotation_value (annotation_group_id, h5ad_value, display_value)
    values (created_annotation_group_id, p_annotation_group_name, p_annotation_group_name)
    returning annotation_value_id into created_annotation_value_id;

    insert into study_annotation_group_ui(study_id, annotation_group_id, is_primary, ordering,
                                          differential_expression_calculated)
    values (p_study_id, created_annotation_group_id, false, 0, false);

    insert into study_sample_annotation (study_id, annotation_value_id, study_sample_ids, color)
    values (p_study_id, created_annotation_value_id, selected_sample_ids_array || unexpressed_sample_ids, '#ff0000');

    SELECT pg_background_launch(format('call user_annotation_diffexp_job(%s, %s)', p_study_id,
                                       created_annotation_group_id))
    INTO worker_pid;
    RAISE WARNING 'user_annotation_diffexp_job(%, %) execution started -> pid %', p_study_id, created_annotation_group_id, worker_pid;
    PERFORM pg_sleep(0.5);
    PERFORM pg_background_detach(worker_pid);

    return created_annotation_group_id;
end
$$;

-- select user_annotation_define(3, 'testgroup', ARRAY[1,2,3,4,5,6]);

drop procedure if exists user_annotation_diffexp_job;
create or replace procedure user_annotation_diffexp_job(study_id int,
                                                        annotation_group_id int)
    LANGUAGE plpython3u
AS
$$

from typing import List
import pandas as pd
from smart_open import open
import scanpy as sc
from anndata import AnnData
import io


def sql_query(query, fetch_results=True):
    # postgres data retrieval with consistent output, both in the jupyter development
    # environment (plpy is not available) and at runtime inside a plpython3u stored procedure
    try:
        import plpy
    except:
        from postgres_utils import engine
        from sqlalchemy import text
        with engine.connect() as connection:
            r = connection.execute(text(query))
            if fetch_results:
                return [row._mapping for row in r.fetchall()]
            else:
                return
    r = plpy.execute(query)
    if fetch_results:
        return [row for row in r]
    else:
        return


def _read_h5ad(filename):
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


def read_h5ad(study_id: int):
    filename = sql_query(f"select coalesce(filename, import_file) filename from study where study_id = {study_id}")[0]['filename']
    if not filename.startswith('s3:'):
        try:
            import plpy
            # running in postgres docker image
            filename = f"/h5ad_store/{filename}"
        except:
            # running in devenv jupyter notebook
            filename = f"../scratch/{filename}"
    return _read_h5ad(filename)


# def add_custom_annotation(adata: AnnData, study_id:int, annotation_group_id:int):

def get_annotation_df(study_id: int, annotation_group_id: int):
    annotation_data = sql_query(f"""
       SELECT ssa.annotation_value_id, ss.h5ad_obs_index
     FROM study_sample_annotation ssa
         join annotation_value av on av.annotation_value_id = ssa.annotation_value_id
     cross join UNNEST(study_sample_ids) as sample_id
     join study_sample ss on ss.study_id = ssa.study_id and ss.study_sample_id = sample_id
            WHERE ssa.study_id = {study_id} and av.annotation_group_id= {annotation_group_id}
        """)
    annotation_df = pd.DataFrame(annotation_data)
    annotation_df.set_index('h5ad_obs_index', inplace=True)
    return annotation_df


def add_sample_annotation(adata: AnnData, annotation_df: pd.DataFrame):
    adata.obs = adata.obs.reset_index().join(annotation_df).fillna(0).astype({'annotation_value_id': 'int32'}).astype(
        {'annotation_value_id': 'str'})


diffexp_attribute = 'annotation_value_id'
diff_exp_min_group_expr = 0.1
diff_exp_min_group_fc = 0.5
diff_exp_max_notgroup_expr = 1
ngenes = 100


def find_valid_attribute_values(adata: AnnData):
    valid_attribute_group_check = (adata.obs[diffexp_attribute].value_counts() > 1)
    attr_values = valid_attribute_group_check.index[valid_attribute_group_check].tolist()
    attr_values = [a for a in attr_values if a != '0']
    return attr_values


def calculate_diff_exp(adata: AnnData, attr_values: List[str]):
    sc.tl.rank_genes_groups(adata, diffexp_attribute, groups=attr_values, method='wilcoxon',
                            use_raw=False, n_genes=ngenes)
    sc.tl.filter_rank_genes_groups(adata, min_in_group_fraction=diff_exp_min_group_expr,
                                   min_fold_change=diff_exp_min_group_fc,
                                   max_out_group_fraction=diff_exp_max_notgroup_expr, use_raw=False,
                                   key="rank_genes_groups", key_added="rank_genes_groups_filtered")
    result_dataframes = []
    for attr_value in attr_values:
        diffexp_df = sc.get.rank_genes_groups_df(adata, key="rank_genes_groups_filtered", group=attr_value)
        diffexp_df = diffexp_df[~diffexp_df["names"].isnull()]
        diffexp_df['ref_attr_value'] = attr_value
        result_dataframes.append(diffexp_df)
    diffexp_df = pd.concat(result_dataframes, axis=0).reset_index(drop=True)
    return diffexp_df


def add_omics_ids_for_names(study_id: int, diffexp_df: pd.DataFrame):
    genes_df = pd.DataFrame(sql_query(f"""select omics_id, h5ad_var_index
        from study_omics where study_id={study_id}""")).set_index('h5ad_var_index')
    genes_index_omics_id_df = adata.var.reset_index().join(genes_df)[['index', 'omics_id']]
    diffexp_db_df = diffexp_df.merge(genes_index_omics_id_df, left_on='names', right_on='index')
    diffexp_db_df = diffexp_db_df.dropna(subset='omics_id').astype({'omics_id': 'int32', 'ref_attr_value': 'int32'})
    return diffexp_db_df


def save_differential_expression(study_id: int, annotation_group_id: int, diffexp_db_df: pd.DataFrame):
    for row in diffexp_db_df.to_dict(orient='records'):
        sql_query(f"""insert into differential_expression (study_id, omics_id, annotation_value_id, pvalue, pvalue_adj, score, log2_foldchange)
                      values ({study_id}, {row['omics_id']}, {row['ref_attr_value']}, {row['pvals']}, {row['pvals_adj']}, {row['scores']}, {row['logfoldchanges']} );""",
                  fetch_results=False)
    sql_query(f"""UPDATE study_annotation_group_ui SET differential_expression_calculated=True
                  WHERE study_id = {study_id} and annotation_group_id = {annotation_group_id};""", fetch_results=False)


adata = read_h5ad(study_id)
annotation_df = get_annotation_df(study_id, annotation_group_id)
add_sample_annotation(adata, annotation_df)
attr_values = find_valid_attribute_values(adata)
diffexp_df = calculate_diff_exp(adata, attr_values)
diffexp_db_df = add_omics_ids_for_names(study_id, diffexp_df)
save_differential_expression(study_id, annotation_group_id, diffexp_db_df)
$$;

--call user_annotation_diffexp_job(3, 32);


drop function if exists user_annotation_edit;
create function user_annotation_edit(p_study_id int,
                                     p_annotation_group_id int,
                                     p_private_to_user boolean
)
    returns boolean
    language plpgsql
    security definer
    volatile
as
$$
declare
    annotation_allowed int;
begin
    select count(1)
    into annotation_allowed
    from user_annotation_group
    where study_id = p_study_id
      and saved_as_annotation_group_id = p_annotation_group_id
      and created_by_user = current_user_email();
    if annotation_allowed = 0 then
        raise 'no permission to edit user annotation';
    end if;
    update user_annotation_group
    set private_to_user = p_private_to_user
    where study_id = p_study_id
      and saved_as_annotation_group_id = p_annotation_group_id;
    return true;
end;
$$;

drop function if exists user_annotation_delete;
create function user_annotation_delete(p_study_id int,
                                       p_annotation_group_id int
)
    returns boolean
    language plpgsql
    security definer
    volatile
as
$$
declare
    annotation_allowed int;
begin
    select count(1)
    into annotation_allowed
    from user_annotation_group
    where study_id = p_study_id
      and saved_as_annotation_group_id = p_annotation_group_id
      and created_by_user = current_user_email();
    if annotation_allowed = 0 then
        raise 'no permission to delete user annotation';
    end if;

    delete
    from differential_expression
    where study_id = p_study_id
      and annotation_value_id in
          (select annotation_value_id from annotation_value where annotation_group_id = p_annotation_group_id);
    delete
    from postgres.public.study_sample_annotation
    where study_id = p_study_id
      and annotation_value_id in
          (select annotation_value_id from annotation_value where annotation_group_id = p_annotation_group_id);
    delete from annotation_value where annotation_group_id = p_annotation_group_id;
    delete from study_annotation_group_ui where annotation_group_id = p_annotation_group_id;
    delete from user_annotation_group where saved_as_annotation_group_id = p_annotation_group_id;
    delete from annotation_group where annotation_group_id = p_annotation_group_id;
    return true;
end;
$$;

