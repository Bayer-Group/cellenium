import React from 'react';
import { createStyles, Group, Navbar, Stack, Title, Tooltip, UnstyledButton } from '@mantine/core';
import { useRecoilState } from 'recoil';
import { NavLink } from 'react-router-dom';

import CelleniumLogo from '../../assets/images/logo.svg';
import { pageState } from '../../atoms';
import { CellTypeMarkerIcon, CoExpressionAnalysisIcon, CompareAnnotationsIcon, ExpressionAnalysisIcon, UserAnnotationIcon } from '../../assets/icons/Icons';

const useStyles = createStyles((theme) => ({
  wrapper: {
    display: 'flex',
  },
  aside: {
    flex: '0 0 60px',
    backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.white,
    display: 'flex',
    flexDirection: 'column',
    alignItems: 'center',
    borderRight: `1px solid ${theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.colors.gray[3]}`,
  },

  main: {
    flex: 1,
    backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[6] : theme.colors.gray[0],
  },

  navigation: {
    padding: 10,
  },
  mainLink: {
    width: 44,
    height: 44,
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
    boxSizing: 'border-box',
    fontFamily: `Exo-bold, ${theme.fontFamily}`,
    fontSize: '1.8rem',
    marginBottom: theme.spacing.md,
    backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.white,
    paddingTop: 10,
    paddingLeft: 15,
    height: 60,
    borderBottom: `1px solid ${theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.colors.gray[3]}`,
  },

  logo: {
    boxSizing: 'border-box',
    width: '100%',
    display: 'flex',
    justifyContent: 'center',
    height: 60,
    paddingTop: theme.spacing.md,
    marginBottom: theme.spacing.md + 4,

    borderBottom: `1px solid ${theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.colors.gray[3]}`,
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

type Props = {
  children?: React.ReactNode;
};

export function LeftSidePanel({ children }: Props) {
  const { classes, cx } = useStyles();
  const [page, setPage] = useRecoilState(pageState);

  const mainLinks = viewLinks.map((link) => (
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
  ));

  return (
    <Navbar style={{ height: '100vh', minWidth: 300 }} width={{ sm: 300 }}>
      <Navbar.Section grow className={classes.wrapper}>
        <div className={classes.aside}>
          <div className={classes.logo}>
            <NavLink to="/" style={{ textDecoration: 'none', color: 'black', userSelect: 'none' }}>
              <img src={CelleniumLogo} alt="cellenium logo" />
            </NavLink>
          </div>
          <Stack mt="md">{mainLinks}</Stack>
        </div>
        <div className={classes.main}>
          <NavLink to="/" style={{ textDecoration: 'none', color: 'black', userSelect: 'none' }}>
            <Title order={4} className={classes.title}>
              cellenium
            </Title>
          </NavLink>
          <Group className={classes.navigation} spacing="md">
            {children}
          </Group>
        </div>
      </Navbar.Section>
    </Navbar>
  );
}
