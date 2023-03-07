# Study-level permissions

In the demo setup, cellenium allows all users to see and administer all studies.

In the `study` postgres table, the arrays `reader_permissions` and `admin_permissions` must be defined per study
to enable study-level permissions.

These are implemented with postgres row level security, restricting what the `postgraphile` database user
can see. All tables and views are explicitly granted to the `postgraphile` database user, some tables
have row level security defined to restrict their rows on what is visible to the current user (with the
`study_visible_currentuser` view).

Postgres intersects the user's groups with those required by a study, and if at least one group overlaps
(or the study has no groups defined at all), access is granted.

The current user's groups are provided through the `current_user_groups()` database function which
can make use of an HTTP request token that is present in GraphQL API requests (via a SQL session variable
that postgraphile sets before executing its SQL to serve the current request). It can evaluate a JWT
or access the organization's group database in any way.

Cellenium does not validate the header in any way, it must be deployed in an environment which makes
sure that only authenticated and (on application level) authorized users have access to it.

If cellenium is running behind an AWS ALB that handles OIDC authentication, the header variable
`x-amzn-oidc-accesstoken` is available. Similarly, AWS API Gateway can be set up to validate JWT tokens
and pass valid API requests, together with the token, to a backend.

For demo purposes, see the comments in `docker-compose.yml` and `current_user_groups()` to just
supply a list of fake groups in an Authentication header in Postgraphile's GraphiQL web UI.
