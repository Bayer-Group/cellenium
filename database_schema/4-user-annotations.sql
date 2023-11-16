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
    job_id                      text;
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

    select aws_batch_submit_cellenium_cli_job('user-annotation-diffexp',
                                              ARRAY ['differential-expression-calculation', p_study_id::text, created_annotation_group_id::text])
    INTO job_id;
    RAISE WARNING 'user_annotation_define for study % launched batch job %', p_study_id, job_id;

    return created_annotation_group_id;
end
$$;

-- select user_annotation_define(3, 'testgroup', ARRAY[1,2,3,4,5,6]);

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

