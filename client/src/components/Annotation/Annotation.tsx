import React from 'react';
import {ColorSwatch, createStyles, Grid, Text} from "@mantine/core";

type Props = {
    onSelect: Function;
    label: string;
    color: string;
    selected: Boolean;
}
const useStyles = createStyles((theme) => ({
    main: {'&:hover': {
      backgroundColor: theme.colorScheme === 'dark' ? theme.colors.dark[5] : theme.colors.gray[1],
      color: theme.colorScheme === 'dark' ? theme.white : theme.black,
    },
    },
}));


function Annotation({onSelect, label, color, selected}: Props) {
    const {classes, cx} = useStyles();

    return (
        <Grid pl={10} gutter={0} sx={{cursor: 'pointer'}} justify={'space-between'} align={'center'} onClick={() => onSelect(label)}
        className={classes.main}
        onMouseOver={()=>console.log(`hovering over ${label}`)}>
            <Grid.Col span={10}>
                <Text size={'md'} fw={selected ? 800 : "md"} truncate>
                    {label}
                </Text>
            </Grid.Col>
            <Grid.Col span={2} style={{justifyContent: 'center'}}>
                <ColorSwatch key={color} color={color} size={15} style={{border: selected ? '2px solid black' : ''}}/>
            </Grid.Col>
        </Grid>
    );
};

export {Annotation};