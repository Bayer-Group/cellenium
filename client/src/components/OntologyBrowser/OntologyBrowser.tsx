import { useState } from 'react';
import { Stack } from '@mantine/core';
import { useOntologiesQuery } from '../../generated/types';
import { OntologyTree } from './OntologyTree';
import OntologySelect from './OntologySelect';

const OntologyBrowser = ({ handleAddOntologyItem, ontologyTrees }: { handleAddOntologyItem: Function; ontologyTrees: any }) => {
  const { data: ontologyData } = useOntologiesQuery();
  const [selectedOntology, setSelectedOntology] = useState<string>();

  return (
    <Stack>
      {ontologyData && <OntologySelect handleChange={setSelectedOntology} ontologies={ontologyData.ontologiesList} />}
      {ontologyTrees && ontologyTrees.get(selectedOntology) !== undefined && selectedOntology && (
        <OntologyTree ontology={ontologyTrees.get(selectedOntology) as any} handleAddOntologyItem={handleAddOntologyItem} />
      )}
    </Stack>
  );
};

export { OntologyBrowser };
