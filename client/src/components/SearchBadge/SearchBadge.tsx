import {ActionIcon, Badge} from "@mantine/core";
import {IconX} from "@tabler/icons-react";
import {OfferingItem} from "../SearchBar/SearchBar";
import {Omics} from "../../model";
import {ontology2Color} from "../../pages/helper";

type Props = {
    onRemove: Function;
    item: OfferingItem | Omics;
}

const SearchBadge = ({onRemove, item}: Props) => {


    const removeButton = (
        <ActionIcon onClick={() => onRemove(item)} size="xs" radius="xl" variant="transparent">
            <IconX color='white' size={15}/>
        </ActionIcon>
    );
    return (
        <div>
            <Badge color={ontology2Color(item.ontology ? item.ontology : '')} radius={4} size={'xl'} variant="filled"
                   sx={{paddingRight: 3, paddingLeft: 8}}
                   rightSection={removeButton}>
                {item.value}
            </Badge>
        </div>
    );
};

export default SearchBadge;