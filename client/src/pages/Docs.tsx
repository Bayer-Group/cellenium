import React, { ComponentType, useEffect, useMemo, useState } from 'react';
import { Center, createStyles, Group, Loader, ScrollArea, Stack, Text } from '@mantine/core';
import { NavLink, useNavigate, useParams } from 'react-router-dom';
import '../assets/docs_styles.css';

const MDComponents: Record<
  string,
  () => Promise<{
    default: ComponentType<never>;
    attributes: { title: string };
  }>
> = import.meta.glob('../docs/*.md') as never;
const DocsKeys = Object.keys(MDComponents).map((k) => k.replace('../docs/', '').replace('.md', '').toLowerCase());
const DocsMd = Object.keys(MDComponents).map((k) => React.lazy(() => MDComponents[k]() as never));

const HEADER_HEIGHT = 60;

const useStyles = createStyles(() => ({
  inner: {
    height: HEADER_HEIGHT,
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'start',
    padding: '2rem',
  },
  sidebar: {
    borderRight: '1px solid gainsboro',
    height: '100%',
    minWidth: '15rem',
    maxWidth: '20rem',
    overflow: 'hidden',
  },
  current: {
    backgroundColor: 'whitesmoke',
  },
  wrapper: {
    marginTop: '0.25rem',
    padding: '0.25rem',
    position: 'relative',
    display: 'block',
    fontWeight: 'bold',
    borderRadius: '5px',
    '&:hover': {
      backgroundColor: 'whitesmoke',
    },
    textDecoration: 'none',
    cursor: 'pointer',
  },
}));

export function Docs() {
  const navigate = useNavigate();
  const params = useParams();
  const { classes } = useStyles();
  const [loading, setLoading] = useState(true);
  const [attributes, setAttributes] = useState<{ title: string; name: string }[] | null>(null);

  useEffect(() => {
    (async () => {
      setLoading(true);
      const attribs = await Promise.all(Object.keys(MDComponents).map((k) => MDComponents[k]())).then((r) =>
        r.map((m, i) => ({ ...m.attributes, name: DocsKeys[i] })),
      );
      setAttributes(() => attribs);
      setLoading(false);
    })();
  }, []);

  useEffect(() => {
    if (!params.pageId || !DocsKeys.includes(params.pageId)) {
      navigate(`/docs/${DocsKeys[0]}`);
    }
  }, [navigate, params.pageId]);

  const DocsComponent = useMemo(() => {
    if (attributes) {
      if (params.pageId && DocsKeys.includes(params.pageId)) {
        return DocsMd[DocsKeys.indexOf(params.pageId)];
      }
    }
    return undefined;
  }, [attributes, params.pageId]);

  return (
    <>
      {loading && (
        <Center h="100%" w="100%">
          <Loader color="blue" />
        </Center>
      )}
      {!loading && (
        <Group align="start" h="100%" noWrap w="100%" spacing={0}>
          <Stack p="1rem" className={classes.sidebar}>
            <ScrollArea h="100%" maw="20rem">
              {attributes
                ? attributes.map((a, i) => (
                    <NavLink
                      to={`/docs/${a.name}`}
                      key={`${a.title}-attribute`}
                      className={`${classes.wrapper} ${params.pageId && i === DocsKeys.indexOf(params.pageId) ? classes.current : ''}`}
                    >
                      <Text color="dimmed" truncate="end" key={a.title} maw="15rem">
                        {a.title}
                      </Text>
                    </NavLink>
                  ))
                : null}
            </ScrollArea>
          </Stack>
          <Stack pl="md" h="100%" w="100%">
            <ScrollArea h="100%" w="100%">
              <Stack w="100%" pb="5rem" p="1rem">
                {DocsComponent && <DocsComponent />}
              </Stack>
            </ScrollArea>
          </Stack>
        </Group>
      )}
    </>
  );
}
