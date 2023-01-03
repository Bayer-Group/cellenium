import React from 'react';
import {ColorSwatch, Grid, Text} from "@mantine/core";

type Props = {
    onSelect: Function;
    label: string;
    color: string;
    selected: Boolean;
}

function Annotation({onSelect, label, color, selected}: Props) {
    return (
        <Grid sx={{cursor: 'pointer'}} justify={'space-between'} align={'center'} onClick={()=>onSelect(label)} >
            <Grid.Col span={10}>
                <Text size={'md'} fw={selected ? 800 : "md"} truncate>
                    {label}
                </Text>
            </Grid.Col>
            <Grid.Col span={2} style={{justifyContent: 'center'}}>
                <ColorSwatch key={color} color={color} size={15} style={{border: selected?'2px solid black':''}}/>
            </Grid.Col>
        </Grid>
    );
};

export {Annotation};