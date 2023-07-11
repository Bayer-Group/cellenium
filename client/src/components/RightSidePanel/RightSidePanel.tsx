import { createStyles, Navbar, Stack } from "@mantine/core";

const useStyles = createStyles((theme) => ({
  wrapper: {
    display: "flex",
  },

  main: {
    flex: 1,
    borderLeft: "1px solid #e9efef",
    paddingLeft: theme.spacing.md,
    paddingTop: theme.spacing.md,
    backgroundColor:
      theme.colorScheme === "dark"
        ? theme.colors.dark[6]
        : theme.colors.gray[0],
  },
}));

type Props = {
  children?: JSX.Element[] | JSX.Element;
};

function RightSidePanel({ children }: Props) {
  const { classes } = useStyles();

  return (
    <Navbar
      height={"100vh"}
      width={{ sm: 300 }}
      style={{ overflowY: "scroll", minWidth: 300 }}
    >
      <Navbar.Section grow className={classes.wrapper}>
        <Stack className={classes.main} spacing={"md"} p={10}>
          {children}
        </Stack>
      </Navbar.Section>
    </Navbar>
  );
}

export { RightSidePanel };
