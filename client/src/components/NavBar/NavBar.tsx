import ProjPlotIcon from '../../images/logo.svg';
import { Burger, Container, createStyles, Group, Header, Paper, Stack, Title, Transition } from '@mantine/core';
import { useDisclosure } from '@mantine/hooks';
import { NavLink } from 'react-router-dom';

const HEADER_HEIGHT = 60;

const useStyles = createStyles((theme) => ({
  inner: {
    height: HEADER_HEIGHT,
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'space-between',
  },

  burger: {
    [theme.fn.largerThan('md')]: {
      display: 'none',
    },
  },

  links: {
    paddingTop: '2.5rem',
    height: HEADER_HEIGHT,
    display: 'flex',
    flexDirection: 'column',
    justifyContent: 'space-between',
    [theme.fn.smallerThan('md')]: {
      display: 'none',
    },
  },

  mainLinks: {
    marginRight: -theme.spacing.sm,
  },

  mainLink: {
    textTransform: 'uppercase',
    textDecoration: 'none',
    fontSize: 13,
    color: theme.colorScheme === 'dark' ? theme.colors.dark[1] : theme.colors.gray[6],
    padding: `2px ${theme.spacing.sm}px`,
    fontWeight: 700,
    borderLeft: '2px solid transparent',
    transition: 'border-color 100ms ease, color 100ms ease',

    '&:hover': {
      color: theme.colorScheme === 'dark' ? theme.white : theme.black,
      textDecoration: 'none',
    },
  },
  mainLinkActive: {
    color: theme.colorScheme === 'dark' ? theme.white : theme.black,
    borderLeftColor: 'black',
  },
  dropdown: {
    position: 'absolute',
    top: HEADER_HEIGHT,
    left: 0,
    right: 0,
    zIndex: 1000,
    borderTopRightRadius: 0,
    borderTopLeftRadius: 0,
    borderTopWidth: 0,
    overflow: 'hidden',

    [theme.fn.largerThan('md')]: {
      display: 'none',
    },
  },
}));

// interface LinkProps {
//     label: string;
//     link: string;
// }

const mainLinks = [
  { link: '/', label: 'Single study analysis' },
  { link: '/crossstudy', label: 'Cross-study analysis' },
  { link: '/markergene', label: 'Marker gene search' },
];

function NavBar() {
  const [opened, { toggle }] = useDisclosure(false);
  const { classes, cx } = useStyles();
  const mainItems = mainLinks.map((item) => (
    <NavLink to={item.link} key={item.label} className={({ isActive }) => (isActive ? cx([classes.mainLink, classes.mainLinkActive]) : classes.mainLink)}>
      {item.label}
    </NavLink>
  ));

  return (
    <Header height={HEADER_HEIGHT} zIndex={1000}>
      <Container className={classes.inner} fluid={true}>
        <NavLink to={'/'} style={{ textDecoration: 'none', color: 'black' }}>
          <Group spacing={5}>
            <img src={ProjPlotIcon} alt="proj plot icon" />
            <Title>cellenium</Title>
          </Group>
        </NavLink>
        <div className={classes.links}>
          <Group spacing={0} position="right" className={classes.mainLinks}>
            {mainItems}
          </Group>
        </div>
        <Burger opened={opened} onClick={toggle} className={classes.burger} size="sm" />
        <Transition transition="scale-y" duration={200} mounted={opened}>
          {(styles) => (
            <Paper className={classes.dropdown} withBorder style={styles}>
              <Stack spacing="md" px="md" py="md" align="center">
                {mainItems}
              </Stack>
            </Paper>
          )}
        </Transition>
      </Container>
    </Header>
  );
}

export { NavBar };
