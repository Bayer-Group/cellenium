import {ReactComponent as ProjPlotIcon} from "../../images/logo.svg";
import {useState} from 'react';
import {Anchor, Burger, Container, createStyles, Group, Title, Header} from '@mantine/core';
import {useDisclosure} from '@mantine/hooks';

const HEADER_HEIGHT = 60;

const useStyles = createStyles((theme) => ({
    inner: {
        height: HEADER_HEIGHT,
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'space-between',
    },

    burger: {
        [theme.fn.largerThan('sm')]: {
            display: 'none',
        },
    },

    links: {
        paddingTop: '2.5rem',
        height: HEADER_HEIGHT,
        display: 'flex',
        flexDirection: 'column',
        justifyContent: 'space-between',
        [theme.fn.smallerThan('sm')]: {
            display: 'none',
        },
    },

    mainLinks: {
        marginRight: -theme.spacing.sm,
    },

    mainLink: {
        textTransform: 'uppercase',
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
}));

interface LinkProps {
    label: string;
    link: string;
}

interface DoubleHeaderProps {
    mainLinks: LinkProps[];
}

function NavBar({mainLinks}: DoubleHeaderProps) {
    const [opened, {toggle}] = useDisclosure(false);
    const {classes, cx} = useStyles();
    const [active, setActive] = useState(0);

    const mainItems = mainLinks.map((item, index) => (
        <Anchor<'a'>
            href={item.link}
            key={item.label}
            className={cx(classes.mainLink, {[classes.mainLinkActive]: index === active})}
            onClick={(event) => {
                event.preventDefault();
                setActive(index);
            }}
        >
            {item.label}
        </Anchor>
    ));


    return (
        <Header height={HEADER_HEIGHT}>
            <Container className={classes.inner}>
                <Group spacing={5}>
                    <ProjPlotIcon/>
                    <Title >cellenium</Title>
                </Group>
                <div className={classes.links}>
                    <Group spacing={0} position="right" className={classes.mainLinks}>
                        {mainItems}
                    </Group>
                </div>
                <Burger opened={opened} onClick={toggle} className={classes.burger} size="sm"/>
            </Container>
        </Header>
    );
}

export {NavBar}