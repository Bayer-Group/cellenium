import { useEffect, useState } from 'react';
import { SelectBoxItem } from '../../model';
import { Center, SegmentedControl } from '@mantine/core';
import { OntologyOverviewFragment } from '../../generated/types';

type OSProps = {
  handleChange: any;
  ontologies: OntologyOverviewFragment[];
};

const OntologySelect = ({ ontologies, handleChange }: OSProps) => {
  const [value, setValue] = useState<string>();

  let selectOntologies: SelectBoxItem[] = ontologies.map((ele) => {
    return {
      value: ele.name,
      label: ele.name,
    };
  });
  useEffect(() => {
    updateChange(selectOntologies[0].value);
  }, []);

  function updateChange(item: string) {
    setValue(item);
    handleChange(item);
  }

  return (
    <Center>
      <SegmentedControl value={value} onChange={updateChange} data={selectOntologies} />
    </Center>
  );
};
export default OntologySelect;
