import React, { useMemo } from 'react';
import { createStyles, Group, Navbar, Stack, Title, Tooltip, UnstyledButton } from '@mantine/core';
import { useRecoilState } from 'recoil';
import { NavLink } from 'react-router-dom';

import CelleniumLogo from '../../assets/images/logo.svg';
import { pageState } from '../../atoms';
import { CellTypeMarkerIcon, CoExpressionAnalysisIcon, CompareAnnotationsIcon, ExpressionAnalysisIcon, UserAnnotationIcon } from '../../assets/icons/Icons';

const useStyles = createStyles((theme) => ({
  wrapper: {
    display: 'flex',
    flexDirection: 'row',
    width: '100%',
    height: '100%',
  },

  aside: {
    backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.white,
    borderRight: `1px solid ${theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.colors.gray[3]}`,
  },

  main: {
    height: '100%',
    flexGrow: 1,
    backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[6] : theme.colors.gray[0],
  },

  mainLink: {
    width: '100%',
    aspectRatio: '1 / 1',
    padding: theme.spacing.xs,
    borderRadius: theme.radius.md,
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    color: theme.colorScheme === 'dark' ? theme.colors.dark[0] : theme.colors.gray[7],
    userSelect: 'none',
    '&:hover': {
      backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[5] : theme.colors.gray[0],
    },
  },

  mainLinkActive: {
    '&, &:hover': {
      backgroundColor: theme.fn.variant({
        variant: 'light',
        color: theme.colors.gray[4],
      }).background,
      color: theme.fn.variant({ variant: 'light', color: theme.colors.gray[4] }).color,
    },
  },

  title: {
    fontFamily: `Exo-bold, ${theme.fontFamily}`,
    fontSize: '1.8rem',
  },

  titleWrapper: {
    borderBottom: `1px solid ${theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.colors.gray[3]}`,
    backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.white,
  },

  logo: {
    height: 60,
    borderBottom: `1px solid ${theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.colors.gray[3]}`,
  },

  logoLink: {
    textDecoration: 'none',
    color: 'black',
    userSelect: 'none',
    width: '100%',
  },

  childrenStack: {
    height: '100%',
    overflowY: 'auto',
  },
}));

const viewLinks = [
  {
    icon: <CellTypeMarkerIcon />,
    label: 'Cell type marker analysis',
    link: 'CellMarkerAnalysis',
  },
  {
    icon: <ExpressionAnalysisIcon />,
    label: 'Expression analysis',
    link: 'ExpressionAnalysis',
  },
  {
    icon: <CoExpressionAnalysisIcon />,
    label: 'Co-Expression analysis',
    link: 'CoexpressionAnalysis',
  },
  {
    // eslint-disable-next-line react/jsx-no-undef
    icon: <UserAnnotationIcon />,
    label: 'Interactive cell type annotation',
    link: 'CelltypeDiscovery',
  },
  {
    icon: <CompareAnnotationsIcon />,
    label: 'Compare annotations',
    link: 'AnnotationComparison',
  },
];

export function LeftSidePanel({ children }: { children?: React.ReactNode }) {
  const { classes, cx } = useStyles();
  const [page, setPage] = useRecoilState(pageState);

  const mainLinks = useMemo(
    () =>
      viewLinks.map((link) => (
        <Tooltip
          label={link.label}
          position="right"
          withArrow
          transitionProps={{
            duration: 0,
          }}
          key={link.label}
        >
          <UnstyledButton
            onClick={() => setPage(link.link)}
            className={cx(classes.mainLink, {
              [classes.mainLinkActive]: link.link === page,
            })}
          >
            {link.icon}
          </UnstyledButton>
        </Tooltip>
      )),
    [classes, cx, page, setPage],
  );

  return (
    <Navbar zIndex={150} h="100%" miw={300} w={300} top={0} left={0} bottom={0}>
      <Navbar.Section grow className={classes.wrapper}>
        <Stack spacing={0} w={60} align="center" className={classes.aside}>
          <NavLink to="/" className={classes.logoLink}>
            <Group className={classes.logo} align="center" spacing={0} p={0} position="center" w={60}>
              <img src={CelleniumLogo} alt="cellenium logo" />
            </Group>
          </NavLink>
          <Stack pt="sm" pb="sm" className={cx(classes.childrenStack, ['no-scrollbar'])}>
            {mainLinks}
          </Stack>
        </Stack>
        <Stack spacing={0} className={classes.main}>
          <NavLink to="/" className={classes.logoLink}>
            <Group align="center" h={60} w="100%" bg="white" className={classes.titleWrapper} pl="sm">
              <Title className={classes.title}>cellenium</Title>
            </Group>
          </NavLink>
          <Stack p="sm" spacing="md" className={cx(classes.childrenStack, ['no-scrollbar'])}>
            {children}
          </Stack>
        </Stack>
      </Navbar.Section>
    </Navbar>
  );
}
