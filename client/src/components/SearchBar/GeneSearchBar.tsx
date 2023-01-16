import {ActionIcon, Autocomplete, AutocompleteItem, Group, Loader, Stack, Text, useMantineTheme} from '@mantine/core';
import React, {useEffect, useState} from "react";
import {IconSearch, IconX} from "@tabler/icons";
import {useHumanGeneAutocompleteQuery, useStudiesWithMarkerGenesLazyQuery} from "../../generated/types";
import {Omics} from "../../model";
import SearchBadge from "../SearchBadge/SearchBadge";
import _ from 'lodash';

const sortAlphaNum = (a: Omics, b: Omics) => a.displaySymbol.localeCompare(b.displaySymbol, 'en', {numeric: true})

interface Props {
    handleNewFilters: Function;
}

function GeneSearchBar({handleNewFilters}: Props) {
    const theme = useMantineTheme();
    const [value, setValue] = useState<string>('')
    const [offerings, setOfferings] = useState<Omics[]>([]);
    const [selectedFilters, setSelectedFilters] = useState<Omics[]>([]);
    const {data, error, loading} = useHumanGeneAutocompleteQuery()
    const [getCellTypes, {
        data: markerData,
        error: markerError,
        loading: markerLoading
    }] = useStudiesWithMarkerGenesLazyQuery();

    function handleSubmit(item: AutocompleteItem) {
        setValue('');
        let newFilters = [...selectedFilters, item as Omics]
        setSelectedFilters(newFilters)
        getCellTypes({
            variables: {
                omicsIds: newFilters.map((ele) => ele.omicsId)
            }
        })
    }

    useEffect(() => {
        if (markerData)
            handleNewFilters(markerData.differentialExpressionsList)
    }, [markerData])

    function handleChange(input: string) {
        if (input === '') {
            setOfferings([])
            setValue('')
            return
        }

        const newOfferings = _.uniqBy(data?.omicsBasesList
            .filter((gene) => gene.displaySymbol.toLowerCase().startsWith(input.toLowerCase()))
            .sort(sortAlphaNum), 'displaySymbol').slice(0, 20);
        if (newOfferings !== undefined)
            setOfferings(newOfferings)
        setValue(input)


    }

    function handleFilterRemove(filter: Omics) {
        let newFilters = selectedFilters.filter((f) => !((f.omicsId === filter.omicsId)));
        if (newFilters.length > 0) {
            setSelectedFilters(newFilters)
            handleNewFilters(newFilters)
        } else {
            setSelectedFilters([])
            handleNewFilters([])
        }

        setOfferings([])
    }


    return (
        <Group position={'left'} align={'flex-end'} spacing={4}>
            <Stack spacing={0} style={{flexGrow: 1}}>
                <Text size={'xs'} weight={800}>
                    Filter by human gene symbol(s)
                </Text>

                < Group spacing={4} position={'left'} align={'center'}
                        style={{'border': '1px lightgray solid', borderRadius: 5, paddingLeft: 4}}
                >
                    {loading ? <Loader size={25} color={theme.colors.blue[5]}/> :
                        <IconSearch size={25} color={theme.colors.gray[3]}/>}
                    <Group spacing={2}>
                        {selectedFilters.map((filter) => {
                            return (
                                <SearchBadge key={`${filter.omicsId}`} onRemove={handleFilterRemove}
                                             item={filter}/>
                            )

                        })}
                    </Group>
                    <div style={{flexGrow: 1}}>
                        <Autocomplete
                            onFocus={() => handleChange(value)}
                            onChange={handleChange}
                            onItemSubmit={handleSubmit}
                            value={value}
                            variant='unstyled'
                            styles={{
                                label: {fontWeight: 500, fontSize: '0.8rem', display: 'inline-block'},
                            }}
                            data={offerings as AutocompleteItem[]}
                            size="md"
                            placeholder={'EGFR, KLK3, CDK2'}
                            rightSection={
                                <ActionIcon onClick={() => {
                                    setValue('');
                                    setOfferings([]);
                                    setSelectedFilters([]);
                                    handleNewFilters([]);
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

export {GeneSearchBar};
