import { ActionIcon, Autocomplete, AutocompleteItem, createStyles, Group, Loader, Stack, Text, useMantineTheme } from '@mantine/core';
import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { IconSearch, IconX } from '@tabler/icons-react';
import { useRecoilState } from 'recoil';
import { Omics } from '../../model';
import { SearchBadge } from '../SearchBadge/SearchBadge';
import { SpeciesSelect } from '../SpeciesSelect/SpeciesSelect';
import { OfferingItem } from './interfaces';
import { OmicsType, useOmicsAutocompleteLazyQuery } from '../../generated/types';
import { GeneSearchSelection, GeneSearchSpeciesSelection, SPECIES } from '../../atoms';

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
  const [species, setSpecies] = useRecoilState(GeneSearchSpeciesSelection);
  const [selectedFilters, setSelectedFilters] = useRecoilState(GeneSearchSelection);
  const inputRef = useRef<HTMLInputElement>(null);
  const speciesList = useMemo(() => (humanOnly ? SPECIES.filter((s) => s.value === '9606') : SPECIES), [humanOnly]);
  const [getAutocomplete, { data: autocompleteSuggestions, loading }] = useOmicsAutocompleteLazyQuery();

  useEffect(() => {
    const timeout = setTimeout(async () => {
      if (value.trim().length > 0) {
        await getAutocomplete({
          variables: {
            taxId: parseInt(species.value, 10),
            omicsType: OmicsType.Gene,
            searchQuery: value,
          },
        });
      }
    }, 50);
    return () => clearTimeout(timeout);
  }, [value, species, getAutocomplete]);

  const offerings = useMemo(() => {
    if (!autocompleteSuggestions || !autocompleteSuggestions.omicsAutocompleteList || value.trim().length === 0) {
      return [];
    }
    return autocompleteSuggestions.omicsAutocompleteList.map((o) => {
      return { ...o, value: o.displaySymbol, ontology: 'GENE', displayName: o.displaySymbol };
    });
  }, [autocompleteSuggestions, value]);

  const handleSubmit = useCallback(
    (item: AutocompleteItem) => {
      const omicsItem = offerings.find((o) => o.value === item.value);
      if (omicsItem) {
        setValue('');
        if (selectedFilters.find((f) => f.value === omicsItem.value) === undefined) {
          const newFilters = [...selectedFilters, omicsItem];
          setSelectedFilters(newFilters);
          onGeneSelection(newFilters.map((f) => f.omicsId).flat());
        }

        if (inputRef && inputRef.current !== null) {
          inputRef.current.focus();
        }
      }
    },
    [offerings, onGeneSelection, selectedFilters, setSelectedFilters],
  );

  const handleChange = useCallback((input: string) => {
    if (input === '') {
      if (inputRef && inputRef.current !== null) {
        inputRef.current.focus();
      }
    }
    setValue(input);
  }, []);

  const updateSpecies = useCallback(
    (item: string) => {
      const localSpecies = SPECIES.find((s) => s.value === item);
      if (localSpecies) {
        setSpecies(localSpecies);
        setSelectedFilters([]);
        setValue('');
      }
    },
    [setSelectedFilters, setSpecies],
  );

  const handleFilterRemove = useCallback(
    (item: Omics | OfferingItem) => {
      if ('ontcode' in item) {
        return;
      }
      const newFilters = selectedFilters.filter((f) => !(f.value === (item as Omics).value));
      if (newFilters.length > 0) {
        setSelectedFilters(newFilters);
      } else {
        setSelectedFilters([]);
      }
      if (inputRef && inputRef.current !== null) {
        inputRef.current.focus();
      }
      onGeneSelection(newFilters.map((f) => f.omicsId).flat());
    },
    [onGeneSelection, selectedFilters, setSelectedFilters],
  );

  const clearInput = useCallback(() => {
    setValue('');
    setSelectedFilters([]);
    onGeneSelection([]);
  }, [onGeneSelection, setSelectedFilters]);

  const onKeyDown = useCallback(
    (event: React.KeyboardEvent<HTMLInputElement>) => {
      if (event.key === 'Backspace' && value.length === 0 && selectedFilters.length > 0) {
        handleFilterRemove(selectedFilters[selectedFilters.length - 1]);
      }
    },
    [handleFilterRemove, selectedFilters, value],
  );

  return (
    <Group position="left" align="flex-end" spacing={4} noWrap>
      <SpeciesSelect data={speciesList} species={species.value} handleChange={updateSpecies} />
      <Stack spacing={0} className={classes.grow}>
        <Text size="xs" weight={800}>
          Enter identifier(s)
        </Text>

        <Group spacing={4} position="left" align="center" className={classes.searchBarWrapper} noWrap>
          {loading ? <Loader size={25} color="blue" /> : <IconSearch size={25} color={theme.colors.gray[3]} />}
          <Group spacing={2} align="center" noWrap>
            {selectedFilters.map((filter) => {
              return <SearchBadge key={`${filter.omicsId}`} onRemove={handleFilterRemove} item={filter as unknown as Omics} />;
            })}
          </Group>
          <Autocomplete
            ref={inputRef}
            autoFocus
            className={classes.autocomplete}
            onItemSubmit={handleSubmit}
            onChange={handleChange}
            onKeyDown={onKeyDown}
            value={value}
            variant="unstyled"
            data={offerings.map((o) => ({ value: o.value, key: o.value })) as AutocompleteItem[]}
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
