import {ActionIcon, Autocomplete, Group, Loader, Stack, Text, TextInputProps, useMantineTheme} from '@mantine/core';
import React, {useEffect, useState} from "react";
import SearchBadge from "../SearchBadge/SearchBadge";
import {IconSearch, IconX} from "@tabler/icons";
import {useAutocompleteLazyQuery} from "../../generated/types";

type OfferingItem = {
    value: string;
    isSynonymOfPreferredTerm: string | null;
    label: string;
    labelHighlight: string;
    ontCode: string;
    ontology: string;
}

function SearchBar(props: TextInputProps) {
    const theme = useMantineTheme();
    const [value, setValue] = useState<string>('')
    const [offerings, setOfferings] = useState<OfferingItem[]>([]);
    const [selectedFilters, setSelectedFilters] = useState<OfferingItem[]>([]);
    const [getAutocomplete, {data: autocompleteSuggestions, error, loading}] = useAutocompleteLazyQuery()


    useEffect(() => {
        if (autocompleteSuggestions) {
            let newOfferings: OfferingItem[] = autocompleteSuggestions.autocompleteList.map((e) => {
                    return {
                        ...e,
                        value: e.label
                    }
                }
            )
            setOfferings(newOfferings);
        }
    }, [autocompleteSuggestions])

    function handleSubmit(item: OfferingItem) {
        let check = selectedFilters.filter((e) => ((e.ontology === item.ontology) && (e.ontCode === item.ontCode)))
        console.log({check})
        console.log({selectedFilters})
        console.log({item})
        if (check.length === 0) {
            setSelectedFilters([...selectedFilters, item])
        }
        setValue('')
        setOfferings([])
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
        let newFilters = selectedFilters.filter((f) => !((f.ontCode === filter.ontCode) && (f.ontology === f.ontology)))
        setSelectedFilters(newFilters)
    }


    return (
        <Stack spacing={0}>
            <Text size={'xs'} weight={800}>
                Filter studies by disease, tissue, species, title/description
            </Text>

            <Group spacing={4} position={'left'} align={'center'}
                   style={{'border': '1px lightgray solid', borderRadius: 5, paddingLeft: 4}}
            >
                <IconSearch color={theme.colors.gray[3]}/>
                <Group spacing={2}>
                    {selectedFilters.map((filter) => {
                        return (
                            <SearchBadge key={`${filter.ontology}_${filter.ontCode}`} onRemove={handleFilterRemove}
                                         item={filter}
                                         color={'blue'}/>
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
                        rightSection={loading ? <Loader size={'xs'} color={theme.colors.blue[3]}/> :
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
    );
}

export {SearchBar};
export type {OfferingItem};