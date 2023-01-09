import React from 'react';
import {Badge, ActionIcon} from "@mantine/core";
import {IconX} from "@tabler/icons";
import {OfferingItem} from "../SearchBar/SearchBar";
type Props = {
    color: string;
    onRemove: Function;
    item: OfferingItem;
}
const SearchBadge = ({color, onRemove, item}: Props) => {
    const removeButton = (
        <ActionIcon onClick={()=>onRemove(item)} size="xs" color="blue" radius="xl" variant="transparent">
            <IconX size={10}/>
        </ActionIcon>
    );
    return (
        <div>
            <Badge radius={4} size={'xl'} variant="outline" sx={{paddingRight: 3}} rightSection={removeButton}>
                {item.label}
            </Badge>
        </div>
    );
};

export default SearchBadge;