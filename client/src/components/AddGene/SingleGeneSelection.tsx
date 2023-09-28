import { ActionIcon, Autocomplete, AutocompleteItem } from '@mantine/core';
import { useCallback, useEffect, useState } from 'react';
import { IconX } from '@tabler/icons-react';
import { useRecoilValue } from 'recoil';
import * as aq from 'arquero';
import { studyState } from '../../atoms';
import { Omics } from '../../model';

export function SingleGeneSelection({ selection, onSelectionChange }: { selection: Omics | null; onSelectionChange: (omics: Omics | null) => void }) {
  const [offerings, setOfferings] = useState<Omics[]>([]);
  const [value, setValue] = useState(selection ? selection.displaySymbol : '');
  useEffect(() => {
    if (value !== selection?.displaySymbol) {
      if (selection && value === '') {
        setValue(selection.displaySymbol);
      }
    }
  }, [selection]);
  const study = useRecoilValue(studyState);

  const handleChange = useCallback(
    (inputString: string) => {
      let newOfferings: Omics[] = [];
      if (inputString.length > 0) {
        newOfferings = study?.studyOmicsTable
          .filter(aq.escape((t: { displaySymbol: string }) => aq.op.includes(t.displaySymbol.toLowerCase(), inputString.toLowerCase(), 0)))
          .objects() as Omics[];
      }
      setOfferings(newOfferings);
      setValue(inputString);
    },
    [study?.studyOmicsTable],
  );

  const handleItemSubmit = useCallback(
    (item: AutocompleteItem) => {
      setOfferings([]);
      if (item) {
        onSelectionChange(item as Omics);
      }
    },
    [onSelectionChange],
  );

  const onBlur = () => {
    handleItemSubmit(offerings[0] as AutocompleteItem);
  };

  return (
    <Autocomplete
      value={value}
      radius={0}
      size="xs"
      onChange={handleChange}
      onBlur={onBlur}
      data={offerings as AutocompleteItem[]}
      onItemSubmit={handleItemSubmit}
      rightSection={
        <ActionIcon
          onClick={() => {
            setValue('');
            setOfferings([]);
            onSelectionChange(null);
          }}
        >
          <IconX />
        </ActionIcon>
      }
    />
  );
}
