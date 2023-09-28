import { useCallback, useEffect, useState } from 'react';
import { Center, SegmentedControl } from '@mantine/core';
import { SelectBoxItem } from '../../model';
import { OntologyOverviewFragment } from '../../generated/types';

export function OntologySelect({ ontologies, handleChange }: { handleChange: (item: string) => void; ontologies: OntologyOverviewFragment[] }) {
  const [value, setValue] = useState<string>();

  const selectOntologies: SelectBoxItem[] = ontologies.map((ele) => {
    return {
      value: ele.name,
      label: ele.name,
    };
  });

  const updateChange = useCallback(
    (item: string) => {
      setValue(item);
      handleChange(item);
    },
    [handleChange],
  );

  useEffect(() => {
    updateChange(selectOntologies[0].value);
  }, []);

  return (
    <Center>
      <SegmentedControl value={value} onChange={updateChange} data={selectOntologies} />
    </Center>
  );
}
