import React, { ComponentType, useEffect, useMemo, useState } from 'react';
import { NavLink } from 'react-router-dom';
import { createStyles, Group, Header, ScrollArea, Stack, Text, Title } from '@mantine/core';
import ProjPlotIcon from '../assets/images/logo.svg';

const MDComponents: Record<
  string,
  () => Promise<{
    default: ComponentType<never>;
    attributes: { title: string };
  }>
> = import.meta.glob('../docs/*.md') as never;
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
  },
}));

export function Docs() {
  const { classes } = useStyles();
  const [selectedPage, setSelectedPage] = useState(0);
  const [attributes, setAttributes] = useState<{ title: string }[] | null>(null);

  useEffect(() => {
    (async () => {
      const attribs = await Promise.all(Object.keys(MDComponents).map((k) => MDComponents[k]())).then((r) => r.map((m) => m.attributes));
      setAttributes(() => attribs);
    })();
  }, []);

  const DocsComponent = useMemo(() => {
    if (DocsMd.length > 0 && selectedPage < DocsMd.length) {
      return DocsMd[selectedPage];
    }
    return undefined;
  }, [selectedPage]);

  return (
    <Stack h="100%" w="100vw" spacing={0} pos="relative">
      <Header height={HEADER_HEIGHT} zIndex={1000} className={classes.inner} w="100%">
        <NavLink to="/" style={{ textDecoration: 'none', color: 'black' }}>
          <Group spacing={5}>
            <img src={ProjPlotIcon} alt="proj plot icon" />
            <Title>cellenium</Title>
          </Group>
        </NavLink>
      </Header>
      <Group align="start" h="100%">
        <Stack p="1rem" className={classes.sidebar}>
          <ScrollArea h="100%" maw="20rem">
            {attributes
              ? attributes.map((a, i) => (
                  <span className={`${classes.wrapper} ${i === selectedPage ? classes.current : ''}`} key={`${a.title}-attribute`}>
                    <Text color="dimmed" truncate="end" key={a.title} onClick={() => setSelectedPage(i)} maw="15rem">
                      {a.title}
                    </Text>
                  </span>
                ))
              : null}
          </ScrollArea>
        </Stack>
        <Stack p="1rem" pl="2rem" h="100%">
          <ScrollArea h="100%">{DocsComponent && <DocsComponent />}</ScrollArea>
        </Stack>
      </Group>
    </Stack>
  );
}
