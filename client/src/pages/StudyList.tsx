import React, {useEffect, useState} from 'react';
import {NavBar, SearchBar, StudyCard} from "../components";
import {Container, Grid, Space} from "@mantine/core";
import {StudyInfoFragment, useStudiesQuery} from "../generated/types";
import {OntologyItem} from "../model";
import {generateOntologyTrees} from "./helper";


const StudyList = () => {
    const {data, error, loading} = useStudiesQuery()
    const [ontologyTrees, setOntologyTrees] = useState<Map<string, OntologyItem>>();
    useEffect(() => {
        if (data)
            setOntologyTrees(generateOntologyTrees(data.treeOntologiesList))
    }, [data])
    return (
        <Container fluid={true}>
            <NavBar/>
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

export default StudyList;