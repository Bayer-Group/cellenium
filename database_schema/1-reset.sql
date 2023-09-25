drop schema if exists public cascade;
create schema public;
create extension if not exists plpython3u;
create extension if not exists pg_background;


-- create postgraphile user (if not exists)
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

grant usage on schema public to postgraphile;


create or replace procedure _analyze_schema() as
$$
DECLARE
    tab RECORD;
BEGIN
    for tab in (select t.relname::varchar AS table_name
                FROM pg_class t
                         JOIN pg_namespace n ON n.oid = t.relnamespace
                WHERE t.relkind = 'r'
                  and n.nspname::varchar = 'public'
                  -- the expression table partitions are analyzed after each study import
                  and t.relname not like 'expression_%'
                order by 1)
        LOOP
            RAISE NOTICE 'ANALYZE %.%', 'public', tab.table_name;
            EXECUTE format('ANALYZE %I.%I', 'public', tab.table_name);
        end loop;
end
$$ LANGUAGE plpgsql;
