import { ActionIcon, Autocomplete, AutocompleteItem, createStyles, Group, Loader, Stack, Text, useMantineTheme } from '@mantine/core';
import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { IconSearch, IconX } from '@tabler/icons-react';
import _ from 'lodash';
import { useRecoilValue } from 'recoil';
import { Omics } from '../../model';
import { SearchBadge } from '../SearchBadge/SearchBadge';
import { allGenesState } from '../../atoms';
import { SpeciesSelect } from '../SpeciesSelect/SpeciesSelect';
import { OfferingItem } from './interfaces';

const sortAlphaNum = (a: Omics, b: Omics) => a.displaySymbol.localeCompare(b.displaySymbol, 'en', { numeric: true });

const SPECIES = [
  { value: '9606', label: 'Homo sapiens' },
  { value: '10090', label: 'Mus musculus' },
  { value: '10116', label: 'Rattus norvegicus' },
];

const useStyles = createStyles((theme) => ({
  grow: {
    flexGrow: 1,
  },
  searchBarWrapper: {
    border: `1px ${theme.colors.gray[3]} solid`,
    borderRadius: 5,
    paddingLeft: 4,
  },
  autocomplete: {
    height: 40,
    width: '100%',
  },
}));

function GeneSearchBar({ humanOnly, onGeneSelection }: { humanOnly: boolean; onGeneSelection: (omicsIds: number[]) => void }) {
  const theme = useMantineTheme();
  const { classes } = useStyles();
  const [value, setValue] = useState<string>('');
  const [offerings, setOfferings] = useState<Omics[]>([]);
  const [selectedFilters, setSelectedFilters] = useState<Omics[]>([]);
  const allGenesRecoilState = useRecoilValue(allGenesState);
  const allGenes = useMemo(() => {
    return allGenesRecoilState || new Map();
  }, [allGenesRecoilState]);
  const [species, setSpecies] = useState<string>(SPECIES[0].value);
  const inputRef = useRef<HTMLInputElement>(null);
  const speciesList = useMemo(() => (humanOnly ? SPECIES.filter((s) => s.value === '9606') : SPECIES), [humanOnly]);

  const handleSubmit = useCallback(
    (item: AutocompleteItem) => {
      setValue('');
      const newFilters = [...selectedFilters, item as Omics];
      setSelectedFilters(newFilters);
      onGeneSelection(newFilters.map((f) => f.omicsId));
      if (inputRef && inputRef.current !== null) {
        inputRef.current.focus();
      }
    },
    [onGeneSelection, selectedFilters],
  );

  const handleChange = useCallback(
    (input: string) => {
      if (input === '') {
        setOfferings([]);
        setValue('');
        if (inputRef && inputRef.current !== null) {
          inputRef.current.focus();
        }
        return;
      }
      const newOfferings = _.uniqBy(
        Array.from(allGenes.values())
          .filter((gene) => gene.taxId === parseInt(species, 10))
          .filter((gene) => gene.displaySymbol.toLowerCase().startsWith(input.toLowerCase()))
          .filter((gene) => !humanOnly || gene.taxId === 9606)
          .sort(sortAlphaNum),
        'displaySymbol',
      )
        .slice(0, 20)
        .map((gene) => {
          return { ...gene, ontology: 'GENE', value: gene.displaySymbol };
        });
      if (newOfferings !== undefined) setOfferings(newOfferings);
      setValue(input);
    },
    [allGenes, humanOnly, species],
  );

  const handleFilterRemove = useCallback(
    (filter: Omics | OfferingItem) => {
      if ('ontcode' in filter) {
        return;
      }
      const newFilters = selectedFilters.filter((f) => !(f.omicsId === (filter as Omics).omicsId));
      if (newFilters.length > 0) {
        setSelectedFilters(newFilters);
      } else {
        setSelectedFilters([]);
      }
      if (inputRef && inputRef.current !== null) {
        inputRef.current.focus();
      }
      setOfferings([]);
      onGeneSelection(newFilters.map((f) => f.omicsId));
    },
    [onGeneSelection, selectedFilters],
  );

  const clearInput = useCallback(() => {
    setValue('');
    setOfferings([]);
    setSelectedFilters([]);
    onGeneSelection([]);
  }, [onGeneSelection]);

  const autoCompleteFocusChange = useCallback(() => {
    handleChange(value);
  }, [handleChange, value]);

  useEffect(() => {
    setValue('');
    setOfferings([]);
    setSelectedFilters([]);
    onGeneSelection([]);
  }, [species]);

  return (
    <Group position="left" align="flex-end" spacing={4} noWrap>
      <SpeciesSelect data={speciesList} species={species} handleChange={setSpecies} />
      <Stack spacing={0} className={classes.grow}>
        <Text size="xs" weight={800}>
          Enter identifier(s)
        </Text>

        <Group spacing={4} position="left" align="center" className={classes.searchBarWrapper} noWrap>
          {!allGenes ? <Loader variant="dots" size={25} color="blue" /> : <IconSearch size={25} color={theme.colors.gray[3]} />}
          <Group spacing={2}>
            {selectedFilters.map((filter) => {
              return <SearchBadge key={`${filter.omicsId}`} onRemove={handleFilterRemove} item={filter} />;
            })}
          </Group>
          <Autocomplete
            ref={inputRef}
            autoFocus
            className={classes.autocomplete}
            disabled={allGenes.size === 0}
            onFocus={autoCompleteFocusChange}
            onChange={handleChange}
            onItemSubmit={handleSubmit}
            value={value}
            variant="unstyled"
            data={offerings as AutocompleteItem[]}
            size="md"
            placeholder="EGFR, KLK3, CDK2, ..."
            rightSection={
              <ActionIcon onClick={clearInput}>
                <IconX />
              </ActionIcon>
            }
          />
        </Group>
      </Stack>
    </Group>
  );
}

export { GeneSearchBar };
