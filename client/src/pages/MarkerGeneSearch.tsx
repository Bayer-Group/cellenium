import React, {useEffect, useState} from 'react';
import {MarkerCard, NavBar} from "../components";
import {Container, Space, Text} from "@mantine/core";
import {GeneSearchBar} from "../components/SearchBar/GeneSearchBar";
import {Omics} from "../model";
import {DifferentialMarkerFragment, useStudiesWithMarkerGenesLazyQuery} from "../generated/types";


const MarkerGeneSearch = () => {
    const [getCellTypes, {data, error, loading}] = useStudiesWithMarkerGenesLazyQuery();
    const [searchResults, setSearchResults] = useState<DifferentialMarkerFragment[]>([]);

    function findMarkerStudies(genes: Omics[]) {
        getCellTypes({
            variables: {
                omicsIds: genes.map((g) => g.omicsId)
            }
        })
    }

    useEffect(() => {
        if (data && data.differentialExpressionsList.length > 0) {
            setSearchResults([...data.differentialExpressionsList])
        }
    }, [data]);

    useEffect(() => {
        console.log({searchResults})
    }, [searchResults])
    return (
        <Container fluid={true}>
            <NavBar/>
            <Space h="xl"/>
            <Container size={'xl'} style={{paddingBottom: '2rem'}}>
                <GeneSearchBar triggerSearch={findMarkerStudies} onDelete={setSearchResults}/>
            </Container>
            <Container size={'xl'}>
                {(searchResults.length > 0) && searchResults.map((sr: any, i) => {
                    return <MarkerCard key={i} data={sr}/>
                })}
                {searchResults.length === 0 && <Text color={'dimmed'}>Please enter your genes of interest.</Text>}
            </Container>
        </Container>
    );
};

export default MarkerGeneSearch;