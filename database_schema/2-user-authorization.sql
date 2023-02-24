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


/*

For testing:
Connect to database with postgraphile user.
Manually set user groups in study table.

Mock groups of the current user (for the dummy implementation of current_user_groups()):
select set_config('postgraphile.auth_header_value','Group1;Group2',false);

Then select e.g. rows from the study table, expression table, or generate plots. The study table
will return less rows depending on the permission definition in each study record.

 */
