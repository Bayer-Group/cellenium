CREATE TABLE user_annotation_group
(
    study_id                          int     not null,
    constraint fk_study_id
        FOREIGN KEY (study_id)
            REFERENCES study (study_id) ON DELETE CASCADE,

    -- the user annotation group is persisted here in annotation_group:
    saved_as_annotation_group_id      int     not null references annotation_group,

    -- user preference regarding this annotation
    calculate_differential_expression boolean not null
);

drop function if exists user_annotation_define;
create function user_annotation_define(p_study_id int,
                                       p_annotation_group_name text,
                                       p_selected_sample_ids int[]
)
    returns integer -- created annotation_group_id
    language plpgsql
    volatile
as
$$
declare
    created_annotation_group_id int;
    created_annotation_value_id int;
begin
    insert into annotation_group (h5ad_column, display_group)
    values (p_study_id || '_' || p_annotation_group_name, p_annotation_group_name)
    returning annotation_group_id into created_annotation_group_id;

    insert into user_annotation_group (study_id, saved_as_annotation_group_id, calculate_differential_expression)
    values (p_study_id, created_annotation_group_id, false);

    insert into annotation_value (annotation_group_id, h5ad_value, display_value, color)
    values (created_annotation_group_id, p_annotation_group_name, p_annotation_group_name, '#ff0000')
    returning annotation_value_id into created_annotation_value_id;

    insert into study_annotation_group_ui(study_id, annotation_group_id, is_primary, ordering,
                                          differential_expression_calculated)
    values (p_study_id, created_annotation_group_id, false, 0, false);

    insert into study_sample_annotation (study_id, annotation_value_id, study_sample_ids)
    values (p_study_id, created_annotation_value_id, p_selected_sample_ids);

    return created_annotation_group_id;
end
$$;

-- select user_annotation_define(3, 'testgroup', ARRAY[1,2,3,4,5,6]);


