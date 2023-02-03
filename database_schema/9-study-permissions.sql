-- If study-level permissions shall be enabled, this function needs to be implemented according to the
-- organization's user group repository, OIDC token format and content etc.
-- The dummy implementation here assumes that the Authentication header contains a list of user groups.
-- This is nice for testing in Postgraphile's GraphiQL webapp, i.e. you can mock a current client's user
-- groups and see how that changes the API responses.
-- In a real use case, e.g. the cloud provider's application load balancer has already authenticated the user,
-- checking that the user is a valid user of your organization and allowed to access the application in general,
-- and retrieved an OIDC token which can contain claims. The current OIDC token is passed by the load balancer
-- to the backend, on each request, using special headers (e.g. x-amzn-oidc-accesstoken in the case of
-- AWS ALB). Our postgraphile_server.js takes care of passing the token header to the SQL context, where
-- it can be retrieved by the current_setting function.
-- The current_user_groups function will be called a lot. In case its response times are slowing cellenium down, consider
-- caching the function result per token input, e.g. using the special postgres SD dictionary variable. You
-- must ensure that no responses are served from the cache beyond the token's lifetime.

CREATE OR REPLACE FUNCTION current_user_groups()
    RETURNS text[]
    LANGUAGE plpython3u
    IMMUTABLE
AS
$$

import plpy


def get_token():
    r = plpy.execute("SELECT current_setting('postgraphile.auth_header_value', TRUE)::VARCHAR token")
    return [row for row in r][0]['token']


def get_groups():
    token = get_token()
    if not token:
        return []
    return token.split(';')


return get_groups()
$$;

DO
$$
    BEGIN
        IF NOT EXISTS(SELECT
                      FROM pg_roles
                      WHERE rolname = 'postgraphile') THEN
            create role postgraphile WITH LOGIN PASSWORD 'postgraphile';
        END IF;
    END
$$;



-- decide which studies can be analyzed by the current user:
-- all studies which are open to everybody ("reader_permissions is null")
-- and all studies which list of allowed groups has at least one group in common with the current user's groups
create view study_visible_currentuser
as
select s.study_id
from study s
where s.reader_permissions is null
   or s.reader_permissions && current_user_groups();


-- make the schema accessible completely to the postgraphile user, and restrict with policies below
grant usage on schema public to postgraphile;
grant all on all tables in schema public to postgraphile;
grant all on all sequences in schema public to postgraphile;
revoke all on study from postgraphile;
grant select on study to postgraphile;

-- protect the core study-specific tables - gene annotations, cell annotation names, layer names don't have to be protected

ALTER TABLE study
    ENABLE ROW LEVEL SECURITY;
CREATE POLICY study_policy ON study FOR SELECT TO postgraphile
    USING (
        study_id in (select study_id
                     from study_visible_currentuser)
    );

ALTER TABLE differential_expression
    ENABLE ROW LEVEL SECURITY;
CREATE POLICY differential_expression_policy ON differential_expression FOR SELECT TO postgraphile
    USING (
        study_id in (select study_id
                     from study_visible_currentuser)
    );

ALTER TABLE expression
    ENABLE ROW LEVEL SECURITY;
CREATE POLICY expression_policy ON expression FOR SELECT TO postgraphile
    USING (
        study_layer_id in (select sl.study_layer_id
                           from study_visible_currentuser
                                    join study_layer sl on study_visible_currentuser.study_id = sl.study_id)
    );

ALTER TABLE user_annotation_group
    ENABLE ROW LEVEL SECURITY;
CREATE POLICY user_annotation_group_policy ON user_annotation_group FOR SELECT TO postgraphile
    USING (
        study_id in (select study_id
                     from study_visible_currentuser)
    );


/*

For testing:
Connect to database with postgraphile user.
Manually set user groups in study table.

Mock groups of the current user (for the dummy implementation of current_user_groups()):
select set_config('postgraphile.auth_header_value','Group1;Group2',false);

Then select e.g. rows from the study table, expression table, or generate plots. The study table
will return less rows depending on the permission definition in each study record.

 */
