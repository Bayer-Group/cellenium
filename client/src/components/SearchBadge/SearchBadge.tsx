import { ActionIcon, Badge } from '@mantine/core';
import { IconX } from '@tabler/icons-react';
import { Omics } from '../../model';
import { ontology2Color } from '../../utils/helper';
import { OfferingItem } from '../SearchBar/interfaces';

export function SearchBadge({ onRemove, item }: { onRemove: (item: OfferingItem | Omics) => void; item: OfferingItem | Omics }) {
  const removeButton = (
    <ActionIcon onClick={() => onRemove(item)} size="xs" radius="xl" variant="transparent">
      <IconX color="white" size={15} />
    </ActionIcon>
  );
  return (
    <div>
      <Badge
        color={ontology2Color(item.ontology ? item.ontology : '')}
        radius={4}
        size="xl"
        variant="filled"
        sx={{ paddingRight: 3, paddingLeft: 8 }}
        rightSection={removeButton}
      >
        {item.value}
      </Badge>
    </div>
  );
}
