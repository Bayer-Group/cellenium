import React, {useState} from 'react';
import {MarkerCard, NavBar} from "../components";
import {Center, Container, Grid, Loader, Space, Text, useMantineTheme} from "@mantine/core";
import {GeneSearchBar} from "../components/SearchBar/GeneSearchBar";
import {
    DifferentialMarkerFragment,
    useStudiesWithMarkerGenesLazyQuery,
    useStudiesWithMarkerGenesQuery
} from "../generated/types";


const MarkerGeneSearch = () => {
    const theme = useMantineTheme();
    const [omicsIds, setOmicsIds] = useState<number[]>([]);
    const {data, loading} = useStudiesWithMarkerGenesQuery({
        variables: {
            omicsIds
        },
        skip: omicsIds.length === 0
    });
    return (
        <Container fluid={true}>
            <NavBar/>
            <Space h="xl"/>
            <Container size={'xl'} style={{paddingBottom: '2rem'}}>
                <GeneSearchBar humanOnly={false} onGeneSelection={ids => setOmicsIds(ids)}/>
            </Container>
            <Container size={'xl'}>
                {loading && <Loader variant={'dots'} color={theme.colors.gray[5]} size={25}/>}
                <Grid>
                    {data?.differentialExpressionsList && data.differentialExpressionsList.map((sr: DifferentialMarkerFragment) => (
                        <Grid.Col span={4}
                                  key={`${sr.study.studyId}_${sr.annotationValueId}`}>
                            <MarkerCard data={sr}/>
                        </Grid.Col>)
                    )}
                </Grid>
                <Space h={'xl'}/>
                <Center style={{'height': '100%', 'width': '100%'}}>
                    {!data?.differentialExpressionsList &&
                        <Text color={'dimmed'}>Please enter your genes of interest. Cellenium will search for studies
                            and cell annotation clusters that show the entered gene as differentially expressed.</Text>}
                </Center>
            </Container>
        </Container>
    );
};

export default MarkerGeneSearch;