import React from 'react';
import {Badge, ActionIcon} from "@mantine/core";
import {IconX} from "@tabler/icons";

type Props = {
    color: string;
    onRemove: Function;
    label: string;
}
const SearchBadge = ({color, onRemove,label}: Props) => {
    const removeButton = (
        <ActionIcon onClick={()=>onRemove(label)} size="xs" color="blue" radius="xl" variant="transparent">
            <IconX size={10}/>
        </ActionIcon>
    );
    return (
        <div>
            <Badge radius={4} size={'xl'} variant="outline" sx={{paddingRight: 3}} rightSection={removeButton}>
                {label}
            </Badge>
        </div>
    );
};

export default SearchBadge;