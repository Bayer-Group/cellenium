import { Select } from '@mantine/core';
import { SelectBoxItem } from '../../model';

export function ExpressionAnalysisTypeSelectBox({
  selection,
  options,
  handleSelection,
}: {
  handleSelection: (value: unknown) => void;
  selection: string;
  options: SelectBoxItem[];
}) {
  return <Select w="100%" labelProps={{ size: 'xs' }} label="Select plot type" value={selection} onChange={handleSelection} data={options} />;
}
