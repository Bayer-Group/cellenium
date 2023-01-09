import {ActionIcon, Autocomplete, Group, Loader, Stack, Text, TextInputProps, useMantineTheme} from '@mantine/core';
import React, {useEffect, useState} from "react";
import SearchBadge from "../SearchBadge/SearchBadge";
import {IconSearch, IconX} from "@tabler/icons";
import {useAutocompleteLazyQuery} from "../../generated/types";
import {DropdownItem} from "../../components";

function SearchBar(props: TextInputProps) {
    const theme = useMantineTheme();
    const [value, setValue] = useState<string>('')
    const [offerings, setOfferings] = useState<string[]>([]);
    const [selectedFilters, setSelectedFilters] = useState<string[]>(['Lung', 'Cough']);
    const [getAutocomplete, {data: autocompleteSuggestions, error, loading}] = useAutocompleteLazyQuery()

    function handleChange(input: string) {
        setOfferings([])
        setValue(input)
        getAutocomplete({
            variables: {
                query: input
            }
        })
    }

    useEffect(() => {
        if (autocompleteSuggestions) {
            let newOfferings:string[] = autocompleteSuggestions.autocompleteList.map((e) => e.label)
            setOfferings(newOfferings);
        }
    }, [autocompleteSuggestions])

    function handleFilterRemove(filter: string) {
        let newFilters = selectedFilters.filter((f) => f !== filter)
        setSelectedFilters(newFilters)
    }
    const dropdownStyled = ()=>{
        return <DropdownItem/>
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
                            <SearchBadge key={filter} onRemove={() => handleFilterRemove(filter)} label={filter}
                                         color={'blue'}/>
                        )

                    })}
                </Group>
                <div style={{flexGrow: 1}}>
                    <Autocomplete
                        value={value}
                        onChange={handleChange}
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