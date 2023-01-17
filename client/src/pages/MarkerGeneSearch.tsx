import React, {useEffect, useState} from 'react';
import {MarkerCard, NavBar} from "../components";
import {Container, Grid, Space, Text} from "@mantine/core";
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
                <GeneSearchBar handleNewFilters={setSearchResults}/>
            </Container>
            <Container size={'xl'}>

                <Grid>
                    {(searchResults.length > 0) && searchResults.map((sr: DifferentialMarkerFragment) => {
                        return <Grid.Col span={4} key={`${sr.study}_${sr.annotationValue}`}><MarkerCard
                            data={sr}/></Grid.Col>
                    })}
                </Grid>
                <Container size={'xl'}>
                    {searchResults.length === 0 &&
                        <Text color={'dimmed'}>Please enter your genes of interest.</Text>}
                </Container>
            </Container>
        </Container>
    );
};

export default MarkerGeneSearch;