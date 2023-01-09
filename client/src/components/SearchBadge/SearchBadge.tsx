import React from 'react';
import {ActionIcon, Badge, useMantineTheme} from "@mantine/core";
import {IconX} from "@tabler/icons";
import {OfferingItem} from "../SearchBar/SearchBar";

type Props = {
    onRemove: Function;
    item: OfferingItem;
}

const SearchBadge = ({onRemove, item}: Props) => {

    const theme = useMantineTheme();

    function ontology2Color(ontology: string) {
        let color = '';
        switch (ontology) {
            case 'NCBI_taxonomy':
                color = 'yellow'//theme.colors.yellow[5];
                break;
            case 'NCIT':
                color = 'red'//theme.colors.red[5];
                break;
            case 'MeSH':
                color = 'violet'//theme.colors.violet[5];
                break;
            default:
                color = 'gray'//theme.colors.gray[5];
                break
        }
        return color;
    }

    const removeButton = (
        <ActionIcon onClick={() => onRemove(item)} size="xs" radius="xl" variant="transparent">
            <IconX color='white' size={15}/>
        </ActionIcon>
    );
    return (
        <div>
            <Badge color={ontology2Color(item.ontology)} radius={4} size={'xl'} variant="filled" sx={{paddingRight: 3, paddingLeft: 8}}
                   rightSection={removeButton}>
                {item.value}
            </Badge>
        </div>
    );
};

export default SearchBadge;