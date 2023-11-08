import { ActionIcon, Autocomplete, Badge, createStyles, Group, Loader, Stack, Text, useMantineTheme } from '@mantine/core';
import React, { forwardRef, useCallback, useEffect, useMemo, useState } from 'react';
import { IconBinaryTree, IconSearch, IconX } from '@tabler/icons-react';
import { closeModal, openModal } from '@mantine/modals';
import { useRecoilState } from 'recoil';
import { StudyOverview, useAutocompleteLazyQuery } from '../../generated/types';
import { OntologyBrowser } from '../OntologyBrowser/OntologyBrowser';
import { Omics, OntologyItem } from '../../model';
import { SearchBadge } from '../SearchBadge/SearchBadge';
import { ontology2Color } from '../../utils/helper';
import { OfferingItem } from './interfaces';
import { StudySearchSelection } from '../../atoms';

interface ItemProps extends React.ComponentPropsWithoutRef<'div'> {
  value: string;
  ontology: string;
  ontcode?: string;
  preferredLabel?: string;
}

const SelectItem = forwardRef<HTMLDivElement, ItemProps>(function SelectItem({ value, ontology, ontcode, preferredLabel, ...others }: ItemProps, ref) {
  return (
    <div ref={ref} {...others}>
      <Group position="apart" align="center" noWrap>
        <Text>{preferredLabel || value}</Text>
        <Badge color={ontology2Color(ontology)}>{ontology}</Badge>
      </Group>
    </div>
  );
});

const useStyles = createStyles(() => ({
  grow: {
    flexGrow: 1,
  },
  filterGroup: {
    border: '1px lightgray solid',
    borderRadius: 5,
  },
}));

export function SearchBar({
  ontologies,
  onSearchFiltersUpdate,
  studies,
  ontologyLoading,
}: {
  ontologyLoading?: boolean;
  ontologies?: Map<string, OntologyItem>;
  onSearchFiltersUpdate: (filters: OfferingItem[]) => void;
  studies?: StudyOverview[];
}) {
  const { classes } = useStyles();
  const theme = useMantineTheme();
  const [value, setValue] = useState<string>('');
  const [selectedFilters, setSelectedFilters] = useRecoilState(StudySearchSelection);
  const [getAutocomplete, { data: autocompleteSuggestions, loading }] = useAutocompleteLazyQuery();

  const offerings = useMemo(() => {
    const newOfferings = [] as OfferingItem[];
    if (value.trim().length > 0) {
      newOfferings.push({
        ontology: 'FREETEXT',
        value,
        ontcode: value,
      });
    }
    if (autocompleteSuggestions) {
      newOfferings.push(
        ...autocompleteSuggestions.autocompleteList.map((e) => {
          return {
            ontcode: e.ontCode,
            ontology: e.ontology,
            value: e.label,
            preferredLabel: e.isSynonymOfPreferredTerm,
          };
        }),
      );
    }
    return newOfferings;
  }, [autocompleteSuggestions, value]);

  useEffect(() => onSearchFiltersUpdate(selectedFilters), [onSearchFiltersUpdate, selectedFilters]);

  const handleSubmit = useCallback(
    (item: OfferingItem) => {
      setValue('');

      const check = selectedFilters.filter((e) => e.ontology === item.ontology && e.ontcode === item.ontcode);
      if (check.length === 0) {
        setSelectedFilters([...selectedFilters, item]);
      }
    },
    [selectedFilters, setSelectedFilters],
  );

  const handleChange = useCallback(
    (input: string) => {
      setValue(input);
      (async () => {
        await getAutocomplete({
          variables: {
            query: input,
          },
        });
      })();
    },
    [getAutocomplete],
  );

  const handleFilterRemove = useCallback(
    (filter: OfferingItem | Omics) => {
      const newFilters = selectedFilters.filter((f) => !(f.ontcode === (filter as OfferingItem).ontcode && f.ontology === filter.ontology));
      setSelectedFilters(newFilters);
    },
    [selectedFilters, setSelectedFilters],
  );

  const showOntologyBrowser = useCallback(() => {
    if (ontologies) {
      openModal({
        modalId: 'ontologyBrowser',
        title: <Text weight={800}>Ontology browser</Text>,
        children: (
          <OntologyBrowser
            ontologyTrees={ontologies}
            studies={studies}
            handleAddOntologyItem={(item: OntologyItem) => {
              setSelectedFilters([
                ...selectedFilters,
                {
                  value: item.label,
                  ontcode: item.id,
                  ontology: item.ontology,
                },
              ]);

              closeModal('ontologyBrowser');
            }}
          />
        ),
      });
    }
  }, [ontologies, selectedFilters, setSelectedFilters, studies]);

  const clearInput = useCallback(() => {
    setValue('');
    setSelectedFilters([]);
  }, [setSelectedFilters]);

  const onKeyDown = useCallback(
    (event: React.KeyboardEvent<HTMLInputElement>) => {
      if (event.key === 'Enter') {
        handleSubmit({
          ontology: 'FREETEXT',
          value,
          ontcode: value,
        });
      } else if (event.key === 'Backspace' && value.length === 0 && selectedFilters.length > 0) {
        handleFilterRemove(selectedFilters[selectedFilters.length - 1]);
      }
    },
    [handleFilterRemove, handleSubmit, selectedFilters, value],
  );

  return (
    <Group position="left" align="flex-end" w="100%" spacing={4} noWrap>
      <ActionIcon onClick={showOntologyBrowser} size="xl" variant="default" disabled={ontologyLoading}>
        {ontologyLoading ? <Loader size="sm" /> : <IconBinaryTree color="lightgray" />}
      </ActionIcon>

      <Stack spacing={0} w="100%">
        <Text size="xs" weight={800}>
          Filter studies by disease (MESH), tissue (NCIT), species, cell type (CO) or title / description
        </Text>

        <Group spacing={4} position="left" align="center" className={classes.filterGroup} pl={4} noWrap>
          {loading ? <Loader size={25} color="blue" /> : <IconSearch size={25} color={theme.colors.gray[3]} />}
          <Group spacing={2} align="center" noWrap>
            {selectedFilters.map((filter) => {
              return <SearchBadge key={`${filter.ontology}_${filter.ontcode}`} onRemove={handleFilterRemove} item={filter} />;
            })}
          </Group>
          <Autocomplete
            onChange={handleChange}
            onKeyDown={onKeyDown}
            onItemSubmit={handleSubmit}
            value={value}
            itemComponent={SelectItem}
            variant="unstyled"
            data={offerings}
            filter={() => true}
            size="md"
            w="100%"
            placeholder="lung, cancer, heart, ..."
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
