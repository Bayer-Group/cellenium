import { Group, Stack, Text } from '@mantine/core';
import { OntologyItem } from '../../model';
import { OntologyBranch } from './OntologyBranch';

export function OntologyTree({ ontology, handleAddOntologyItem }: { ontology: OntologyItem; handleAddOntologyItem: (item: OntologyItem) => void }) {
  return (
    <Group position="center" align="center" grow>
      <Stack spacing={0} justify="center" align="flex-start">
        <Text size="xs" color="dimmed">
          Choose from:
        </Text>
        {[ontology].map((item: OntologyItem) => (
          <OntologyBranch handleAddOntologyItem={handleAddOntologyItem} key={item.id} item={item} level={0} />
        ))}
      </Stack>
    </Group>
  );
}
