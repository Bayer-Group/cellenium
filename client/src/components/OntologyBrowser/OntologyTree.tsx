import { Group, Stack, Text } from '@mantine/core';
import { OntologyItem } from '../../model';
import OntologyBranch from './OntologyBranch';

type Props = {
  ontology: OntologyItem;
  handleAddOntologyItem: Function;
};

const OntologyTree = ({ ontology, handleAddOntologyItem }: Props) => {
  return (
    <Group position={'center'} align={'center'} grow>
      <Stack spacing={0} justify={'center'} align={'flex-start'}>
        <Text size={'xs'} color={'dimmed'}>
          Choose from:
        </Text>
        {[ontology].map((item: OntologyItem) => (
          <OntologyBranch handleAddOntologyItem={handleAddOntologyItem} key={item.id} item={item} level={0}></OntologyBranch>
        ))}
      </Stack>
    </Group>
  );
};

export { OntologyTree };
