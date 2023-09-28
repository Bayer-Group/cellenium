import { useState } from 'react';
import { Center, Container, Grid, Loader, Space, Text, useMantineTheme } from '@mantine/core';
import { DifferentialMarkerFragment, useStudiesWithMarkerGenesQuery } from '../generated/types';
import { NavBar } from '../components/NavBar/NavBar';
import { GeneSearchBar } from '../components/SearchBar/GeneSearchBar';
import { MarkerCard } from '../components/MarkerCard/MarkerCard';

function MarkerGeneSearch() {
  const theme = useMantineTheme();
  const [omicsIds, setOmicsIds] = useState<number[]>([]);
  const { data, loading } = useStudiesWithMarkerGenesQuery({
    variables: {
      omicsIds,
    },
    skip: omicsIds.length === 0,
  });

  return (
    <Container fluid>
      <NavBar />
      <Space h="xl" />
      <Container size="xl" style={{ paddingBottom: '2rem' }}>
        <GeneSearchBar humanOnly={false} onGeneSelection={(ids: number[]) => setOmicsIds(ids)} />
      </Container>
      <Container size="xl">
        <Center style={{ width: '100%' }} mb="sm">
          <Text color="dimmed">
            Please enter your genes of interest. Cellenium will search for studies and cell annotation clusters that show the entered gene as differentially
            expressed.
          </Text>
          {loading && <Loader variant="dots" color={theme.colors.gray[5]} size={25} />}
        </Center>
        <Grid>
          {data?.differentialExpressionsList &&
            data.differentialExpressionsList.map((sr: DifferentialMarkerFragment) => (
              <Grid.Col span={4} key={`${sr.study.studyId}_${sr.annotationValueId}`}>
                <MarkerCard data={sr} />
              </Grid.Col>
            ))}
        </Grid>
        <Space h="xl" />
      </Container>
    </Container>
  );
}

export default MarkerGeneSearch;
