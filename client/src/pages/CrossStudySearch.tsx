import React, {useEffect, useState} from 'react';
import {NavBar, SearchBar, StudyCard} from "../components";
import {Container, Grid, Space} from "@mantine/core";
import {StudyInfoFragment, TreeOntologyOverviewFragment, useStudiesQuery} from "../generated/types";
import {OntologyItem} from "../model";
import {generateOntologyTrees} from "./helper";


const CrossStudySearch = () => {
    const {data, error, loading} = useStudiesQuery()

    return (
        <Container fluid={true}>
            <NavBar/>
            <Space h="xl"/>

        </Container>
    );
};

export default CrossStudySearch;