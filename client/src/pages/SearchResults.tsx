import React from 'react';
import {NavBar, SearchBar, StudyCard} from "../components";
import {Container, Grid, Space} from "@mantine/core";
import STUDIES from '../data/studies.json';

const SearchResults = () => {
    return (
        <Container fluid={true}>
            <NavBar mainLinks={[{link: 'single_studies', label: 'Single study analysis'},
                {link: 'cross_study', label: 'Cross study analysis'},
                {link: 'marker_gene', label: 'Marker gene search'}
            ]}/>
            <Space h="xl"/>
            <Container size={'md'} style={{paddingBottom: '2rem'}}>
                <SearchBar/>
            </Container>
            <Container size={'md'}>
                <Grid>
                    {STUDIES.map((study: any) => {
                        return <Grid.Col span={12}><StudyCard content={study}/></Grid.Col>
                    })}
                </Grid>
            </Container>
        </Container>
    );
};

export default SearchResults;