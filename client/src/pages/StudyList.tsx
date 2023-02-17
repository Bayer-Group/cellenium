import React, {useEffect, useMemo, useState} from 'react';
import {NavBar, StudyCard} from "../components";
import {Container, Grid, Space} from "@mantine/core";
import {StudyInfoFragment} from "../generated/types";
import StudySearchBar from "../components/SearchBar/StudySearchBar";


const StudyList = () => {
    const [studyList, setStudyList] = useState<StudyInfoFragment[]>([]);

    return (
        <Container fluid={true}>
            <NavBar/>
            <Space h="xl"/>
            <Container size={'xl'} style={{paddingBottom: '2rem'}}>
                <StudySearchBar onStudyListUpdate={setStudyList}/>
            </Container>
            <Container size={'xl'}>
                <Grid>
                    {studyList.map(study => {
                        return <Grid.Col span={12} key={study.studyId}><StudyCard study={study}/></Grid.Col>
                    })}
                </Grid>
            </Container>
        </Container>
    );
};

export default StudyList;