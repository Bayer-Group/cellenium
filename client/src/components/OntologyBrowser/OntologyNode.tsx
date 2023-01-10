import React from 'react';
import {OntologyItem} from "../../model";
import {ActionIcon, createStyles, Group, Text, useMantineTheme} from "@mantine/core";
import {IconCaretRight, IconCaretDown} from "@tabler/icons";



const useStyles = createStyles((theme) => ({
    main: {
        paddingLeft: '2px',
        paddingRight: '2px',
        '&:hover':{
            backgroundColor: theme.colors.gray[1],
            borderRadius: '2px'
        }
    },
}));

type Props = {
    item: OntologyItem;
    selected: Boolean;
    hasChildren: Boolean | undefined;
    level: number;
    onToggle: () => void;
    handleAddOntologyItem: Function;
}
const OntologyNode = ({item, selected, hasChildren, level, onToggle, handleAddOntologyItem}: Props) => {
    const theme = useMantineTheme();
    const {cx, classes} = useStyles();
    return (
        <Group pl={`${level * 10}px`} spacing={0}>
            {hasChildren && <ActionIcon size='xs' variant={'subtle'} onClick={onToggle}>
                {selected?<IconCaretDown color={theme.colors.dark[9]}/>:<IconCaretRight color={theme.colors.dark[9]}/>}
            </ActionIcon>}
            <Text onClick={()=>handleAddOntologyItem(item.label)} className={classes.main} size={'xs'} style={{cursor:'pointer'}}>{item.label}</Text>
        </Group>
    );
};

export default OntologyNode;