import { createStyles, Navbar, Stack } from '@mantine/core';
import { ReactNode } from 'react';

const useStyles = createStyles((theme) => ({
  wrapper: {
    display: 'flex',
  },

  main: {
    flex: 1,
    borderLeft: '1px solid #e9efef',
    paddingLeft: theme.spacing.md,
    paddingTop: theme.spacing.md,
    backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[6] : theme.colors.gray[0],
  },
}));

export function RightSidePanel({ children }: { children?: ReactNode }) {
  const { classes } = useStyles();

  return (
    <Navbar height="100vh" width={{ sm: 300 }} style={{ overflowY: 'scroll', overflowX: 'hidden', minWidth: 300 }}>
      <Navbar.Section grow className={classes.wrapper}>
        <Stack className={classes.main} spacing="md" p={10}>
          {children}
        </Stack>
      </Navbar.Section>
    </Navbar>
  );
}
