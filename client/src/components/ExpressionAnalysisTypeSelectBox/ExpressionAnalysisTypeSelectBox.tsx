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
  return (
    <Select
      style={{ minWidth: 210, maxWidth: 210, width: 210 }}
      labelProps={{ size: 'xs' }}
      label="Select plot type"
      value={selection}
      onChange={handleSelection}
      data={options}
    />
  );
}
