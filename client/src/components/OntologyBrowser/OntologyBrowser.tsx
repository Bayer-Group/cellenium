import { useState } from 'react';
import { Stack } from '@mantine/core';
import { StudyOverview, useOntologiesQuery } from '../../generated/types';
import { OntologyTree } from './OntologyTree';
import { OntologySelect } from './OntologySelect';
import { OntologyItem } from '../../model';

function OntologyBrowser({
  handleAddOntologyItem,
  ontologyTrees,
  studies,
}: {
  handleAddOntologyItem: (item: OntologyItem) => void;
  ontologyTrees: Map<string, OntologyItem>;
  studies?: StudyOverview & { allOntCodes: string[] }[];
}) {
  const { data: ontologyData } = useOntologiesQuery();
  const [selectedOntology, setSelectedOntology] = useState<string>('');

  return (
    <Stack>
      {ontologyData && <OntologySelect handleChange={setSelectedOntology} ontologies={ontologyData.ontologiesList} />}
      {ontologyTrees && ontologyTrees.get(selectedOntology) !== undefined && selectedOntology && (
        <OntologyTree studies={studies} ontology={ontologyTrees.get(selectedOntology) as OntologyItem} handleAddOntologyItem={handleAddOntologyItem} />
      )}
    </Stack>
  );
}

export { OntologyBrowser };
