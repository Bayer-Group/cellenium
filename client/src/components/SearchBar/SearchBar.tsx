import { ActionIcon, Autocomplete, Badge, Group, Loader, Stack, Text, useMantineTheme } from '@mantine/core';
import React, { forwardRef, useEffect, useState } from 'react';
import { IconBinaryTree, IconSearch, IconX } from '@tabler/icons-react';
import { useAutocompleteLazyQuery } from '../../generated/types';
import { closeModal, openModal } from '@mantine/modals';
import { OntologyBrowser } from '../OntologyBrowser/OntologyBrowser';
import { OntologyItem } from '../../model';
import SearchBadge from '../SearchBadge/SearchBadge';
import { ontology2Color } from '../../pages/helper';

export type OfferingItem = {
  value: string;
  ontcode?: string;
  ontology: string;
  preferredLabel?: string;
};

interface ItemProps extends React.ComponentPropsWithoutRef<'div'> {
  value: string;
  ontology: string;
  ontcode?: string;
  preferredLabel?: string;
}

const SelectItem = forwardRef<HTMLDivElement, ItemProps>(({ value, ontology, ontcode, preferredLabel, ...others }: ItemProps, ref) => (
  <div ref={ref} {...others}>
    <Group position={'apart'} align={'center'} noWrap>
      <Text>{preferredLabel || value}</Text>
      <Badge color={ontology2Color(ontology)}>{ontology}</Badge>
    </Group>
  </div>
));

interface Props {
  ontologies?: Map<string, OntologyItem>;
  onSearchFiltersUpdate: (filters: OfferingItem[]) => void;
}

function SearchBar({ ontologies, onSearchFiltersUpdate }: Props) {
  const theme = useMantineTheme();
  const [value, setValue] = useState<string>('');
  const [offerings, setOfferings] = useState<OfferingItem[]>([]);
  const [selectedFilters, setSelectedFilters] = useState<OfferingItem[]>([]);
  const [getAutocomplete, { data: autocompleteSuggestions, loading }] = useAutocompleteLazyQuery();

  useEffect(() => {
    const newOfferings = [] as OfferingItem[];
    if (value.length > 0) {
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
      setOfferings(newOfferings);
    }
  }, [value, autocompleteSuggestions]);

  useEffect(() => onSearchFiltersUpdate(selectedFilters), [selectedFilters]);

  function handleSubmit(item: OfferingItem) {
    setValue('');

    let check = selectedFilters.filter((e) => e.ontology === item.ontology && e.ontcode === item.ontcode);
    if (check.length === 0) {
      setSelectedFilters([...selectedFilters, item]);
    }
    setOfferings([]);
  }

  function handleChange(input: string) {
    setOfferings([]);
    setValue(input);
    getAutocomplete({
      variables: {
        query: input,
      },
    });
  }

  function handleFilterRemove(filter: OfferingItem) {
    let newFilters = selectedFilters.filter((f) => !(f.ontcode === filter.ontcode && f.ontology === f.ontology));
    setSelectedFilters(newFilters);
    setOfferings([]);
  }

  function showOntologyBrowser() {
    if (ontologies) {
      openModal({
        modalId: 'ontologyBrowser',
        title: <Text weight={800}>Ontology browser</Text>,
        children: (
          <OntologyBrowser
            ontologyTrees={ontologies}
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
  }

  return (
    <Group position={'left'} align={'flex-end'} spacing={4}>
      <ActionIcon onClick={showOntologyBrowser} size={'xl'} variant={'default'}>
        <IconBinaryTree color={'lightgray'} />
      </ActionIcon>

      <Stack spacing={0} style={{ flexGrow: 1 }}>
        <Text size={'xs'} weight={800}>
          Filter studies by disease (MESH), tissue (NCIT), species, cell type (CO) or title / description
        </Text>

        <Group
          spacing={4}
          position={'left'}
          align={'center'}
          style={{
            border: '1px lightgray solid',
            borderRadius: 5,
            paddingLeft: 4,
          }}
        >
          {loading ? <Loader size={25} color={theme.colors.gray[5]} /> : <IconSearch size={25} color={theme.colors.gray[3]} />}
          <Group spacing={2}>
            {selectedFilters.map((filter) => {
              return <SearchBadge key={`${filter.ontology}_${filter.ontcode}`} onRemove={handleFilterRemove} item={filter} />;
            })}
          </Group>
          <div style={{ flexGrow: 1 }}>
            <Autocomplete
              onChange={handleChange}
              onItemSubmit={handleSubmit}
              value={value}
              itemComponent={SelectItem}
              variant="unstyled"
              styles={{
                label: {
                  fontWeight: 500,
                  fontSize: '0.8rem',
                  display: 'inline-block',
                },
              }}
              data={offerings}
              filter={() => true}
              size="md"
              placeholder={'lung, cancer, heart'}
              rightSection={
                <ActionIcon
                  onClick={() => {
                    setValue('');
                    setOfferings([]);
                    setSelectedFilters([]);
                  }}
                >
                  <IconX />
                </ActionIcon>
              }
            />
          </div>
        </Group>
      </Stack>
    </Group>
  );
}

export { SearchBar };
export type { OfferingItem };
