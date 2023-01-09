import React from 'react';
import {NavBar, SearchBar, StudyCard} from "../components";
import {Container, Grid, Space} from "@mantine/core";
import {StudyOverviewFragment, useStudiesQuery} from "../generated/types";

const SearchResults = () => {
    const {data, error, loading} = useStudiesQuery()

    return (
        <Container fluid={true}>
            <NavBar mainLinks={[{link: 'single_studies', label: 'Single study analysis'},
                {link: 'cross_study', label: 'Cross study analysis'},
                {link: 'marker_gene', label: 'Marker gene search'}
            ]}/>
            <Space h="xl"/>
            <Container size={'xl'} style={{paddingBottom: '2rem'}}>
                {data && data.treeDiseasesList && data.treeTissuesList && <SearchBar diseases={data.treeDiseasesList} tissues={data.treeTissuesList}/>}
            </Container>
            <Container size={'xl'}>
                <Grid>
                    {data && data.studiesList && data.studiesList.map((study: StudyOverviewFragment) => {
                        return <Grid.Col span={12} key={study.studyId}><StudyCard study={study}
                                                                                  tissues={data.treeTissuesList}
                                                                                  diseases={data.treeDiseasesList}/></Grid.Col>
                    })}
                </Grid>
            </Container>
        </Container>
    );
};

export default SearchResults;