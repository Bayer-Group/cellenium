import { createStyles, Navbar, Stack } from '@mantine/core';
import { ReactNode } from 'react';
import { StudyTitle } from '../StudyTitle/StudyTitle.tsx';

const useStyles = createStyles((theme) => ({
  wrapper: {
    display: 'flex',
    width: '100%',
    height: '100%',
  },
  nav: {
    height: '100%',
    top: 0,
    bottom: 0,
    right: 0,
  },
  main: {
    flex: 1,
    borderLeft: '1px solid #e9efef',
    paddingLeft: theme.spacing.md,
    paddingTop: theme.spacing.md,
    backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[6] : theme.colors.gray[0],
    overflowY: 'auto',
    overflowX: 'hidden',
  },
}));

export function RightSidePanel({ children }: { children?: ReactNode }) {
  const { classes, cx } = useStyles();

  return (
    <Navbar miw={300} zIndex={150} w={300} m={0} className={classes.nav}>
      <Navbar.Section grow className={classes.wrapper}>
        <Stack className={cx(classes.main, ['no-scrollbar'])} spacing="md" p="md">
          <StudyTitle />
          {children}
        </Stack>
      </Navbar.Section>
    </Navbar>
  );
}
