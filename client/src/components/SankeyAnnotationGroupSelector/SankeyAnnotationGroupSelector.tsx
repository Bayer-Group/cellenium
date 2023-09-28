import { Select, Stack } from '@mantine/core';
import { SelectBoxItem } from '../../model';

export function SankeyAnnotationGroupSelector({
  annotationGroups,
  value1,
  value2,
  handleChange1,
  handleChange2,
}: {
  annotationGroups: SelectBoxItem[];
  handleChange1: (v: string) => void;
  value1: string | undefined;
  handleChange2: (v: string) => void;
  value2: string | undefined;
}) {
  return (
    <Stack spacing={10}>
      <Select
        value={value1}
        onChange={handleChange1}
        label="Select annotation group 1"
        labelProps={{ size: 'xs' }}
        placeholder="Pick one"
        transitionProps={{
          duration: 80,
          timingFunction: 'ease',
        }}
        data={annotationGroups}
      />
      <Select
        value={value2}
        onChange={handleChange2}
        label="Select annotation group 2"
        labelProps={{ size: 'xs' }}
        placeholder="Pick one"
        transitionProps={{
          duration: 80,
          timingFunction: 'ease',
        }}
        data={annotationGroups.filter((ele) => ele.value !== value1)}
      />
    </Stack>
  );
}
