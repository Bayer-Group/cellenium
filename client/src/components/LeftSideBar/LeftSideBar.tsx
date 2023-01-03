import {useState} from 'react';
import {createStyles, Navbar, Title, Tooltip, UnstyledButton} from '@mantine/core';
import {ReactComponent as StudyAnalysisIcon} from "../../icons/study_analysis.svg";
import {ReactComponent as ExpressionAnalysisIcon} from "../../icons/expression_analysis.svg";
import {ReactComponent as CoExpressionAnalysisIcon} from "../../icons/coexpression_analysis.svg";
import {ReactComponent as CompareAnnotationsIcon} from "../../icons/annotation_comparison.svg";
import {ReactComponent as UserAnnotationIcon} from "../../icons/user_annotation.svg";

import {ReactComponent as CelleniumLogo} from "../../images/logo.svg";
import {IconEdit} from "@tabler/icons";

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
        borderRight: `1px solid ${
            theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.colors.gray[3]
        }`,
    },

    main: {
        flex: 1,
        backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[6] : theme.colors.gray[0],
    },

    mainLink: {
        width: 44,
        height: 44,
        borderRadius: theme.radius.md,
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        color: theme.colorScheme === 'dark' ? theme.colors.dark[0] : theme.colors.gray[7],

        '&:hover': {
            backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[5] : theme.colors.gray[0],
        },
    },

    mainLinkActive: {
        '&, &:hover': {
            backgroundColor: theme.fn.variant({variant: 'light', color: theme.primaryColor}).background,
            color: theme.fn.variant({variant: 'light', color: theme.primaryColor}).color,
        },
    },

    title: {
        boxSizing: 'border-box',
        fontFamily: `Exo-bold, ${theme.fontFamily}`,
        fontSize: '1.8rem',
        marginBottom: theme.spacing.xl,
        backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.white,
        paddingTop: 10,
        paddingLeft: 15,
        height: 60,
        borderBottom: `1px solid ${
            theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.colors.gray[3]
        }`,
    },

    logo: {
        boxSizing: 'border-box',
        width: '100%',
        display: 'flex',
        justifyContent: 'center',
        height: 60,
        paddingTop: theme.spacing.md,
        borderBottom: `1px solid ${
            theme.colorScheme === 'dark' ? theme.colors.dark[7] : theme.colors.gray[3]
        }`,
        marginBottom: theme.spacing.xl,
    },

    link: {
        boxSizing: 'border-box',
        display: 'block',
        textDecoration: 'none',
        borderTopRightRadius: theme.radius.md,
        borderBottomRightRadius: theme.radius.md,
        color: theme.colorScheme === 'dark' ? theme.colors.dark[0] : theme.colors.gray[7],
        padding: `0 ${theme.spacing.md}px`,
        fontSize: theme.fontSizes.sm,
        marginRight: theme.spacing.md,
        fontWeight: 500,
        height: 44,
        lineHeight: '44px',

        '&:hover': {
            backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[5] : theme.colors.gray[1],
            color: theme.colorScheme === 'dark' ? theme.white : theme.black,
        },
    },

    linkActive: {
        '&, &:hover': {
            borderLeftColor: theme.fn.variant({variant: 'filled', color: theme.primaryColor})
                .background,
            backgroundColor: theme.fn.variant({variant: 'filled', color: theme.primaryColor})
                .background,
            color: theme.white,
        },
    },
}));

const mainLinksMockdata = [
    {icon: StudyAnalysisIcon, label: 'Marker analysis'},
    {icon: ExpressionAnalysisIcon, label: 'Expression analysis'},
    {icon: CoExpressionAnalysisIcon, label: 'Co-Expression analysis'},
    {icon: CompareAnnotationsIcon, label: 'Compare annotations'},
    {icon: UserAnnotationIcon, label: 'User annotations'},


];

const linksMockdata = [
    'Security',
    'Settings',
    'Dashboard',
    'Releases',
    'Account',
    'Orders',
    'Clients',
    'Databases',
    'Pull Requests',
    'Open Issues',
    'Wiki pages',
];

function LeftSideNavBar() {
    const {classes, cx} = useStyles();
    const [active, setActive] = useState('asdf');
    const [activeLink, setActiveLink] = useState('Settings');

    const mainLinks = mainLinksMockdata.map((link) => (
        <Tooltip label={link.label} position="right" withArrow transitionDuration={0} key={link.label}>
            <UnstyledButton
                onClick={() => setActive(link.label)}
                className={cx(classes.mainLink, {[classes.mainLinkActive]: link.label === active})}
            >
                <link.icon/>
            </UnstyledButton>
        </Tooltip>
    ));

    const links = linksMockdata.map((link) => (
        <a
            className={cx(classes.link, {[classes.linkActive]: activeLink === link})}
            href="/"
            onClick={(event) => {
                event.preventDefault();
                setActiveLink(link);
            }}
            key={link}
        >
            {link}
        </a>
    ));

    return (
        <Navbar height={'100vh'} width={{sm: 300}}>
            <Navbar.Section grow className={classes.wrapper}>
                <div className={classes.aside}>
                    <div className={classes.logo}>
                        <CelleniumLogo/>
                    </div>
                    {mainLinks}
                </div>
                <div className={classes.main}>
                    <Title order={4} className={classes.title}>
                        cellenium
                    </Title>

                    {links}
                </div>
            </Navbar.Section>
        </Navbar>
    );
}

export {LeftSideNavBar}