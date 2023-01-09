import React from 'react';
import {OntologyItem} from "../../model";
import {ActionIcon, Group, Text} from "@mantine/core";
import {IconCaretRight, IconCaretDown} from "@tabler/icons";

type Props = {
    item: OntologyItem;
    selected: Boolean;
    hasChildren: Boolean | undefined;
    level: number;
    onToggle: () => void;
}
const OntologyNode = ({item, selected, hasChildren, level, onToggle}: Props) => {
    return (
        <Group style={{paddingLeft: `${level * 10}px`}} spacing={0}>
            {hasChildren && <ActionIcon size='xs' variant={'transparent'} onClick={onToggle}>
                {selected?<IconCaretDown/>:<IconCaretRight/>}
            </ActionIcon>}
            <Text size={'xs'}>{item.label}</Text>
        </Group>
    );
};

export default OntologyNode;