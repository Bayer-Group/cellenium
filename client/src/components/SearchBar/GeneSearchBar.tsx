import {ActionIcon, Autocomplete, AutocompleteItem, Group, Loader, Stack, Text, useMantineTheme} from '@mantine/core';
import React, {useEffect, useMemo, useRef, useState} from "react";
import {IconSearch, IconX} from "@tabler/icons";
import {useStudiesWithMarkerGenesLazyQuery} from "../../generated/types";
import {Omics} from "../../model";
import SearchBadge from "../SearchBadge/SearchBadge";
import _ from 'lodash';
import {useRecoilValue} from "recoil";
import {allGenesState} from "../../atoms";
import {SpeciesSelect} from "../SpeciesSelect/SpeciesSelect";

const sortAlphaNum = (a: Omics, b: Omics) => a.displaySymbol.localeCompare(b.displaySymbol, 'en', {numeric: true})

interface Props {
    humanOnly: boolean,
    handleNewFilters?: Function;
    onGeneSelection?: (omicsIds: number[]) => void;
}

const SPECIES = [
    {value: "9606", label: 'Homo sapiens'},
    {value: "10090", label: 'Mus musculus'},
    {value: "10116", label: 'Rattus norvegicus'},
]

function GeneSearchBar({humanOnly, handleNewFilters, onGeneSelection}: Props) {
    // either just a gene input (onGeneSelection parameter is used) or
    // also providing the studies that include the gene as differentially expressed (handleNewFilters parameter)

    const theme = useMantineTheme();
    const [value, setValue] = useState<string>('')
    const [offerings, setOfferings] = useState<Omics[]>([]);
    const [selectedFilters, setSelectedFilters] = useState<Omics[]>([]);
    const allGenes = useRecoilValue(allGenesState) || new Map();
    const [species, setSpecies] = useState<string>(SPECIES[0].value);
    const inputRef = useRef<HTMLInputElement>(null);
    const speciesList = useMemo(() => humanOnly ? SPECIES.filter(s => s.value === "9606") : SPECIES, []);

    const [getCellTypes, {
        data: markerData,
        error: markerError,
        loading: markerLoading
    }] = useStudiesWithMarkerGenesLazyQuery();

    function handleSubmit(item: AutocompleteItem) {
        setValue('');
        let newFilters = [...selectedFilters, item as Omics];
        setSelectedFilters(newFilters);
        onGeneSelection && onGeneSelection(newFilters.map(f => f.omicsId));
        if (handleNewFilters) {
            getCellTypes({
                variables: {
                    omicsIds: newFilters.map((ele) => ele.omicsId)
                }
            })
        }
        if (inputRef && inputRef.current !== null) {
            inputRef.current.focus()
        }
    }

    useEffect(() => {
        if (markerData) {
            handleNewFilters && handleNewFilters(markerData.differentialExpressionsList);
        }
    }, [markerData])

    function handleChange(input: string) {
        if (input === '') {
            setOfferings([])
            setValue('')
            if (inputRef && inputRef.current !== null) {
                inputRef.current.focus()
            }
            return
        }
        const newOfferings = _.uniqBy(Array.from(allGenes.values())
            .filter((gene) => gene.taxId === parseInt(species))
            .filter((gene) => gene.displaySymbol.toLowerCase().startsWith(input.toLowerCase()))
            .filter(gene => humanOnly === false || gene.taxId === 9606)
            .sort(sortAlphaNum), 'displaySymbol').slice(0, 20).map((gene) => {
            return {...gene, ontology: 'GENE', value: gene.displaySymbol}
        });
        if (newOfferings !== undefined)
            setOfferings(newOfferings)
        setValue(input)

    }

    function handleFilterRemove(filter: Omics) {
        let newFilters = selectedFilters.filter((f) => !((f.omicsId === filter.omicsId)));
        if (newFilters.length > 0) {
            setSelectedFilters(newFilters)
            handleNewFilters && handleNewFilters(newFilters)
        } else {
            setSelectedFilters([])
            handleNewFilters && handleNewFilters([])
        }
        if (inputRef && inputRef.current !== null) {
            inputRef.current.focus()
        }
        setOfferings([])
    }


    return (
        <Group position={'left'} align={'flex-end'} spacing={4} noWrap>
            <SpeciesSelect data={speciesList} species={species} handleChange={setSpecies}/>
            <Stack spacing={0} style={{flexGrow: 1}}>
                <Text size={'xs'} weight={800}>
                    Enter gene symbol(s)
                </Text>

                < Group spacing={4} position={'left'} align={'center'}
                        style={{'border': `1px ${theme.colors.gray[3]} solid`, borderRadius: 5, paddingLeft: 4}}
                        noWrap
                >
                    {!allGenes ? <Loader variant={'dots'} size={25} color={theme.colors.gray[5]}/> :
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
                            ref={inputRef}
                            autoFocus
                            style={{
                                height: 40
                            }}
                            disabled={allGenes.size === 0}
                            onFocus={() => handleChange(value)}
                            onChange={handleChange}
                            onItemSubmit={handleSubmit}
                            value={value}
                            variant='unstyled'
                            data={offerings as AutocompleteItem[]}
                            styles={{
                                label: {fontWeight: 500, fontSize: '0.8rem', display: 'inline-block'},
                            }}
                            size={'md'}
                            placeholder={'EGFR, KLK3, CDK2'}
                            rightSection={
                                <ActionIcon onClick={() => {
                                    setValue('');
                                    setOfferings([]);
                                    setSelectedFilters([]);
                                    handleNewFilters && handleNewFilters([]);
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
