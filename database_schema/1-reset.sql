drop schema if exists public cascade;
create schema public;
create extension if not exists plpython3u;
create extension if not exists pg_background;

create procedure _analyze_schema() as
$$
DECLARE
    tab RECORD;
BEGIN
    for tab in (select t.relname::varchar AS table_name
                FROM pg_class t
                         JOIN pg_namespace n ON n.oid = t.relnamespace
                WHERE t.relkind = 'r'
                  and n.nspname::varchar = 'public'
                order by 1)
        LOOP
            RAISE NOTICE 'ANALYZE %.%', 'public', tab.table_name;
            EXECUTE format('ANALYZE %I.%I', 'public', tab.table_name);
        end loop;
end
$$ LANGUAGE plpgsql;
