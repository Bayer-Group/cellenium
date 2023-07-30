import React from "react";
import {
  Anchor,
  createStyles,
  Group,
  Navbar,
  Title,
  Tooltip,
  UnstyledButton,
} from "@mantine/core";
import CellTypeMarkerIcon from "../../icons/study_analysis.svg";
import ExpressionAnalysisIcon from "../../icons/expression_analysis.svg";
import CoExpressionAnalysisIcon from "../../icons/coexpression_analysis.svg";
import CompareAnnotationsIcon from "../../icons/annotation_comparison.svg";
import UserAnnotationIcon from "../../icons/user_annotation.svg";

import CelleniumLogo from "../../images/logo.svg";
import { useRecoilState } from "recoil";
import { pageState } from "../../atoms";

const useStyles = createStyles((theme) => ({
  wrapper: {
    display: "flex",
  },

  aside: {
    flex: "0 0 60px",
    backgroundColor:
      theme.colorScheme === "dark" ? theme.colors.dark[7] : theme.white,
    display: "flex",
    flexDirection: "column",
    alignItems: "center",
    borderRight: `1px solid ${
      theme.colorScheme === "dark" ? theme.colors.dark[7] : theme.colors.gray[3]
    }`,
  },

  main: {
    flex: 1,
    backgroundColor:
      theme.colorScheme === "dark"
        ? theme.colors.dark[6]
        : theme.colors.gray[0],
  },

  navigation: {
    padding: 10,
  },
  mainLink: {
    width: 44,
    height: 44,
    borderRadius: theme.radius.md,
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
    color:
      theme.colorScheme === "dark"
        ? theme.colors.dark[0]
        : theme.colors.gray[7],

    "&:hover": {
      backgroundColor:
        theme.colorScheme === "dark"
          ? theme.colors.dark[5]
          : theme.colors.gray[0],
    },
  },

  mainLinkActive: {
    "&, &:hover": {
      backgroundColor: theme.fn.variant({
        variant: "light",
        color: theme.colors.gray[4],
      }).background,
      color: theme.fn.variant({ variant: "light", color: theme.colors.gray[4] })
        .color,
    },
  },

  title: {
    boxSizing: "border-box",
    fontFamily: `Exo-bold, ${theme.fontFamily}`,
    fontSize: "1.8rem",
    marginBottom: theme.spacing.md,
    backgroundColor:
      theme.colorScheme === "dark" ? theme.colors.dark[7] : theme.white,
    paddingTop: 10,
    paddingLeft: 15,
    height: 60,
    borderBottom: `1px solid ${
      theme.colorScheme === "dark" ? theme.colors.dark[7] : theme.colors.gray[3]
    }`,
  },

  logo: {
    boxSizing: "border-box",
    width: "100%",
    display: "flex",
    justifyContent: "center",
    height: 60,
    paddingTop: theme.spacing.md,
    marginBottom: theme.spacing.md + 4,

    borderBottom: `1px solid ${
      theme.colorScheme === "dark" ? theme.colors.dark[7] : theme.colors.gray[3]
    }`,
  },
}));

const viewLinks = [
  {
    icon: <img src={CellTypeMarkerIcon} alt="cell type marker icon" />,
    label: "Cell type marker analysis",
    link: "CellMarkerAnalysis",
  },
  {
    icon: <img src={ExpressionAnalysisIcon} alt="expression analysis icon" />,
    label: "Expression analysis",
    link: "ExpressionAnalysis",
  },
  {
    icon: <img src={CoExpressionAnalysisIcon} alt="co-expression icon" />,
    label: "Co-Expression analysis",
    link: "CoexpressionAnalysis",
  },
  {
    icon: <img src={UserAnnotationIcon} alt="user annotations icon" />,
    label: "Interactive cell type annotation",
    link: "CelltypeDiscovery",
  },
  {
    icon: <img src={CompareAnnotationsIcon} alt="compare annotations icon" />,
    label: "Compare annotations",
    link: "AnnotationComparison",
  },
];

type Props = {
  children?: React.ReactNode;
};

function LeftSidePanel({ children }: Props) {
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
    <Navbar style={{ height: "100vh", minWidth: 300 }} width={{ sm: 300 }}>
      <Navbar.Section grow className={classes.wrapper}>
        <div className={classes.aside}>
          <div className={classes.logo}>
            <Anchor href={"/"}>
              <img src={CelleniumLogo} alt="cellenium logo" />
            </Anchor>
          </div>
          {mainLinks}
        </div>
        <div className={classes.main}>
          <Title order={4} className={classes.title}>
            cellenium
          </Title>
          <Group className={classes.navigation} spacing={"md"}>
            {children}
          </Group>
        </div>
      </Navbar.Section>
    </Navbar>
  );
}

export { LeftSidePanel };
