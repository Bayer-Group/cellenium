import React, {useEffect, useState} from 'react';
import {NavBar, SearchBar, StudyCard} from "../components";
import {Container, Grid, Space} from "@mantine/core";
import {StudyInfoFragment, TreeOntologyOverviewFragment, useStudiesQuery} from "../generated/types";
import {OntologyItem} from "../model";
import {generateOntologyTrees} from "./helper";


const SearchResults = () => {
    const {data, error, loading} = useStudiesQuery()
    const [ontologyTrees, setOntologyTrees] = useState<Map<string, OntologyItem>>();
    useEffect(()=>{
        if (data)
            setOntologyTrees(generateOntologyTrees(data.treeOntologiesList))
    },[data])
    return (
        <Container fluid={true}>
            <NavBar mainLinks={[{link: 'single_studies', label: 'Single study analysis'},
                {link: 'cross_study', label: 'Cross study analysis'},
                {link: 'marker_gene', label: 'Marker gene search'}
            ]}/>
            <Space h="xl"/>
            <Container size={'xl'} style={{paddingBottom: '2rem'}}>
                {ontologyTrees && <SearchBar ontologies={ontologyTrees}/>}
            </Container>
            <Container size={'xl'}>
                <Grid>
                    {data && data.studyOverviewsList && data.studyOverviewsList.map((study: StudyInfoFragment) => {

                        return <Grid.Col span={12} key={study.studyId}><StudyCard study={study}/></Grid.Col>
                    })}
                </Grid>
            </Container>
        </Container>
    );
};

export default SearchResults;