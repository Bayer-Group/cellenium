import { useState } from 'react';
import { Stack } from '@mantine/core';
import { useOntologiesQuery } from '../../generated/types';
import { OntologyTree } from './OntologyTree';
import { OntologySelect } from './OntologySelect';
import { OntologyItem } from '../../model';

function OntologyBrowser({
  handleAddOntologyItem,
  ontologyTrees,
}: {
  handleAddOntologyItem: (item: OntologyItem) => void;
  ontologyTrees: Map<string, OntologyItem>;
}) {
  const { data: ontologyData } = useOntologiesQuery();
  const [selectedOntology, setSelectedOntology] = useState<string>('');

  return (
    <Stack>
      {ontologyData && <OntologySelect handleChange={setSelectedOntology} ontologies={ontologyData.ontologiesList} />}
      {ontologyTrees && ontologyTrees.get(selectedOntology) !== undefined && selectedOntology && (
        <OntologyTree ontology={ontologyTrees.get(selectedOntology) as OntologyItem} handleAddOntologyItem={handleAddOntologyItem} />
      )}
    </Stack>
  );
}

export { OntologyBrowser };
