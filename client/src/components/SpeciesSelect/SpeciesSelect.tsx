import { Select, SelectItem } from '@mantine/core';

export function SpeciesSelect({ handleChange, species, data }: { species: string; handleChange: (item: string) => void; data: SelectItem[] }) {
  return <Select onChange={handleChange} value={species} variant="default" size="md" data={data} label="Select species" labelProps={{ fw: 700, size: 'xs' }} />;
}
