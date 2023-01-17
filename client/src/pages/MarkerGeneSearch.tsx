import React, {useEffect, useState} from 'react';
import {MarkerCard, NavBar} from "../components";
import {Container, Space, Text} from "@mantine/core";
import {GeneSearchBar} from "../components/SearchBar/GeneSearchBar";
import {DifferentialMarkerFragment} from "../generated/types";


const MarkerGeneSearch = () => {
    const [searchResults, setSearchResults] = useState<DifferentialMarkerFragment[]>([]);


    useEffect(() => {
        console.log({searchResults})
    }, [searchResults])
    return (
        <Container fluid={true}>
            <NavBar/>
            <Space h="xl"/>
            <Container size={'xl'} style={{paddingBottom: '2rem'}}>
                <GeneSearchBar humanOnly={false} handleNewFilters={setSearchResults}/>
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