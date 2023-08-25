import { Select, SelectItem } from '@mantine/core';

interface Props {
  species: string;
  handleChange: Function;
  data: SelectItem[];
}

const SpeciesSelect = ({ handleChange, species, data }: Props) => {
  return (
    <Select
      style={{ borderColor: '#000' }}
      onChange={(value: string) => handleChange(value)}
      value={species}
      variant={'default'}
      size={'md'}
      data={data}
      label={'Select species'}
      labelProps={{ fw: 700, size: 'xs' }}
    />
  );
};

export { SpeciesSelect };
