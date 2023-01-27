import React, {useEffect, useMemo, useState} from 'react';
import {NavBar, SearchBar, StudyCard} from "../components";
import {Container, Grid, Space} from "@mantine/core";
import {StudyInfoFragment, useStudiesQuery} from "../generated/types";
import {OntologyItem} from "../model";
import {generateOntologyTrees} from "./helper";
import {flatten} from "lodash";


const StudyList = () => {
    const {data, error, loading} = useStudiesQuery();
    const allStudies = useMemo(() => data?.studyOverviewsList && data?.studyOverviewsList.map(study => ({
        ...study,
        allOntCodes: [
            study.studyOntologyList.map(ont => ont.ontCodes),
            study.studyOntologyList.map(ont => ont.parentIds)
        ].flat(2).filter(ontCode => !!ontCode)
    })), [data]);

    const ontologyTrees: Map<string, OntologyItem> | undefined = useMemo(
        () => data?.treeOntologiesList && generateOntologyTrees(data.treeOntologiesList),
        [data]);
    const [searchOntCodes, setSearchOntCodes] = useState<string[]>([]);

    const filteredStudies = useMemo(() => {
        if (!allStudies) {
            return undefined;
        }
        if (searchOntCodes.length === 0) {
            return allStudies;
        }
        return allStudies.filter(study => searchOntCodes.find(searchOntCode => study.allOntCodes.indexOf(searchOntCode) > -1));
    }, [data?.studyOverviewsList, searchOntCodes]);

    return (
        <Container fluid={true}>
            <NavBar/>
            <Space h="xl"/>
            <Container size={'xl'} style={{paddingBottom: '2rem'}}>
                {ontologyTrees && <SearchBar ontologies={ontologyTrees} onSearchElementsUpdate={setSearchOntCodes}/>}
            </Container>
            <Container size={'xl'}>
                <Grid>
                    {filteredStudies && filteredStudies.map((study: StudyInfoFragment) => {
                        return <Grid.Col span={12} key={study.studyId}><StudyCard study={study}/></Grid.Col>
                    })}
                </Grid>
            </Container>
        </Container>
    );
};

export default StudyList;