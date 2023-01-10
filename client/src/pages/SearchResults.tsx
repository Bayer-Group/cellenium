import React, {useEffect, useState} from 'react';
import {NavBar, SearchBar, StudyCard} from "../components";
import {Container, Grid, Space} from "@mantine/core";
import {StudyOverviewFragment, TreeOntologyOverviewFragment, useStudiesQuery} from "../generated/types";
import {OntologyItem} from "../model";


function generateOntologyTrees(nodeList: TreeOntologyOverviewFragment[]) {
    // @ts-ignore
    const ontologies = [...new Set(nodeList.map(item => item.ontology))];
    const ontologyItemMap = new Map<string, OntologyItem>()
    // setup the hash
    nodeList.map((nd) => {
        ontologyItemMap.set(`${nd.ontology}_${nd.ontCode}`, {
            id: nd.ontCode,
            label: nd.label,
            ontology: nd.ontology,
            children:[]});
    })
    // fill the children
     nodeList.map((nd)=>{

     });
    console.log({ontologyItemMap})
}

const SearchResults = () => {
    const {data, error, loading} = useStudiesQuery()
    const [ontologyTrees, setOntologyTrees] = useState<Map<string, OntologyItem>>();

    useEffect(() => {
        if (data && data.treeOntologiesList) {
            generateOntologyTrees(data.treeOntologiesList)
        }
    }, [data])
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
                    {data && data.studiesList && data.studiesList.map((study: StudyOverviewFragment) => {
                        return <Grid.Col span={12} key={study.studyId}><StudyCard study={study}
                                                                                  tissues={study.tissueNcitIds}
                                                                                  diseases={study.diseaseMeshIds}/></Grid.Col>
                    })}
                </Grid>
            </Container>
        </Container>
    );
};

export default SearchResults;