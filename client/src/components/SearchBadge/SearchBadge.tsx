import { ActionIcon, Badge } from '@mantine/core';
import { IconX } from '@tabler/icons-react';
import { useMemo } from 'react';
import { Omics } from '../../model';
import { ontology2Color } from '../../utils/helper';
import { OfferingItem } from '../SearchBar/interfaces';

export function SearchBadge({ onRemove, item }: { onRemove: (item: OfferingItem | Omics) => void; item: OfferingItem | Omics }) {
  const removeButton = useMemo(() => {
    return (
      <ActionIcon onClick={() => onRemove(item)} size="xs" radius="xl" variant="transparent">
        <IconX color="white" size={15} />
      </ActionIcon>
    );
  }, [item, onRemove]);

  const color = useMemo(() => {
    return ontology2Color(item.ontology ? item.ontology : '');
  }, [item.ontology]);

  return (
    <Badge color={color} radius={4} size="xl" variant="filled" pl={8} pr={3} rightSection={removeButton}>
      {item.value}
    </Badge>
  );
}
