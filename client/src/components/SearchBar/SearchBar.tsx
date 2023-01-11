import {ActionIcon, Autocomplete, Group, Loader, Stack, Text, useMantineTheme} from '@mantine/core';
import React, {useEffect, useState} from "react";
import SearchBadge from "../SearchBadge/SearchBadge";
import {IconBinaryTree, IconSearch, IconX} from "@tabler/icons";
import { useAutocompleteLazyQuery} from "../../generated/types";
import {closeModal, openModal} from "@mantine/modals";
import {OntologyBrowser} from "../OntologyBrowser/OntologyBrowser";
import {OntologyItem} from "../../model";

type OfferingItem = {
    value: string;
    ontcode: string;
    ontology: string;
}

type Props = {
    ontologies: Map<string,OntologyItem>;
}

function SearchBar({ontologies}: Props) {
    const theme = useMantineTheme();
    const [value, setValue] = useState<string>('')
    const [offerings, setOfferings] = useState<OfferingItem[]>([]);
    const [selectedFilters, setSelectedFilters] = useState<OfferingItem[]>([]);
    const [getAutocomplete, {data: autocompleteSuggestions, error, loading}] = useAutocompleteLazyQuery()


    useEffect(() => {
        if (autocompleteSuggestions) {
            let newOfferings: OfferingItem[] = autocompleteSuggestions.autocompleteList.map((e) => {
                    return {
                        ontcode: e.ontCode,
                        ontology: e.ontology,
                        value: e.label
                    }
                }
            )
            setOfferings(newOfferings);
        }
    }, [autocompleteSuggestions])

    function handleSubmit(item: OfferingItem) {
        setValue('');

        let check = selectedFilters.filter((e) => ((e.ontology === item.ontology) && (e.ontcode === item.ontcode)))
        if (check.length === 0) {
            setSelectedFilters([...selectedFilters, item])
        }
        setOfferings([]);
    }

    function handleChange(input: string) {
        setOfferings([])
        setValue(input)
        getAutocomplete({
            variables: {
                query: input
            }
        })
    }

    function handleFilterRemove(filter: OfferingItem) {
        let newFilters = selectedFilters.filter((f) => !((f.ontcode === filter.ontcode) && (f.ontology === f.ontology)))
        setSelectedFilters(newFilters)
        setOfferings([])
    }

    function showOntologyBrowser() {
        openModal({

            modalId: 'ontologyBrowser',
            title: 'Ontology browser',
            children: <OntologyBrowser ontologyTrees={ontologies} handleAddOntologyItem={(item: OntologyItem) => {
                setSelectedFilters([...selectedFilters, {
                    value: item.label,
                    ontcode: item.id,
                    ontology: item.ontology
                }]);

                closeModal('ontologyBrowser');
            }}/>
        })

    }

    return (
        <Group position={'left'} align={'flex-end'} spacing={4}>
            <ActionIcon onClick={showOntologyBrowser} size={'xl'} variant={'light'}>
                <IconBinaryTree/>
            </ActionIcon>

            <Stack spacing={0} style={{flexGrow: 1}}>
                <Text size={'xs'} weight={800}>
                    Filter studies by disease, tissue, species, title/description
                </Text>

                <Group spacing={4} position={'left'} align={'center'}
                       style={{'border': '1px lightgray solid', borderRadius: 5, paddingLeft: 4}}
                >
                    {loading ? <Loader size={25} color={theme.colors.blue[5]}/> :
                        <IconSearch size={25} color={theme.colors.gray[3]}/>}
                    <Group spacing={2}>
                        {selectedFilters.map((filter) => {
                            return (
                                <SearchBadge key={`${filter.ontology}_${filter.ontcode}`} onRemove={handleFilterRemove}
                                             item={filter}/>
                            )

                        })}
                    </Group>
                    <div style={{flexGrow: 1}}>
                        <Autocomplete
                            onChange={handleChange}
                            onItemSubmit={handleSubmit}
                            value={value}
                            variant='unstyled'
                            styles={{
                                label: {fontWeight: 500, fontSize: '0.8rem', display: 'inline-block'},
                            }}
                            data={offerings}
                            size="md"
                            placeholder='lung, "multiple myelome", heart, mouse'
                            rightSection={
                                <ActionIcon onClick={() => {
                                    setValue('');
                                    setOfferings([]);
                                    setSelectedFilters([]);
                                }
                                }>
                                    <IconX/>
                                </ActionIcon>}
                        />
                    </div>
                </Group>
            </Stack>
        </Group>
    );
}

export {SearchBar};
export type {OfferingItem};