import { useMemo, useState } from 'react';
import { Center, Container, Loader, Space, Text, Box, Accordion } from '@mantine/core';
import { DifferentialMarkerFragment, useStudiesWithMarkerGenesQuery } from '../generated/types';
import { NavBar } from '../components/NavBar/NavBar';
import { GeneSearchBar } from '../components/SearchBar/GeneSearchBar';
import { MarkerCard } from '../components/MarkerCard/MarkerCard';

interface StudySummary {
  studyId: number;
  studyName: string;
  differentialExpressionsList: DifferentialMarkerFragment[];
}

function StudySummaryLabel({ studySummary }: { studySummary: StudySummary }) {
  return (
    <div>
      <Text>{studySummary.studyName}</Text>
      <Text size={'xs'}>max log2FC: {Math.max(...studySummary.differentialExpressionsList.map((de) => de.log2Foldchange)).toFixed(2)}</Text>
    </div>
  );
}

function StudyMarkerGeneDetails({ studySummary }: { studySummary: StudySummary }) {
  return (
    <Box pl="xs" style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit,26em)' }}>
      {studySummary.differentialExpressionsList.map((sr: DifferentialMarkerFragment) => (
        <MarkerCard data={sr} key={`${sr.study.studyId}_${sr.annotationValueId}`} />
      ))}
    </Box>
  );
}

function MarkerGeneSearch() {
  const [omicsIds, setOmicsIds] = useState<number[]>([]);
  const { data, loading } = useStudiesWithMarkerGenesQuery({
    variables: {
      omicsIds,
    },
    skip: omicsIds.length === 0,
  });

  const studySummaries = useMemo(() => {
    if (data?.differentialExpressionsList) {
      const buildStudySummaries: StudySummary[] = new Array();
      data.differentialExpressionsList.forEach((sr) => {
        let studySummary = buildStudySummaries.find((s) => s.studyId === sr.study.studyId);
        if (!studySummary) {
          studySummary = {
            studyId: sr.study.studyId,
            studyName: sr.study.studyName,
            differentialExpressionsList: [],
          } as StudySummary;
          buildStudySummaries.push(studySummary);
        }
        studySummary.differentialExpressionsList.push(sr);
      });
      buildStudySummaries.forEach((s) => {
        s.differentialExpressionsList.sort((a, b) => b.log2Foldchange - a.log2Foldchange);
      });
      buildStudySummaries.sort((a, b) => b.differentialExpressionsList[0].log2Foldchange - a.differentialExpressionsList[0].log2Foldchange);
      return buildStudySummaries;
    }
  }, [data]);

  return (
    <Container fluid>
      <NavBar />
      <Space h="xl" />
      <Container size="xl" style={{ paddingBottom: '2rem' }}>
        <GeneSearchBar humanOnly={false} onGeneSelection={(ids: number[]) => setOmicsIds(ids)} />

        <Center style={{ width: '100%' }} m="sm">
          {loading && <Loader variant="dots" color="blue" size={25} />}
          {omicsIds.length === 0 && (
            <Text color="dimmed">
              Please enter your genes of interest. Cellenium will search for studies and cell annotation clusters that show the entered gene as differentially
              expressed.
            </Text>
          )}
        </Center>
      </Container>

      {studySummaries && studySummaries.length > 0 && (
        <Accordion chevronPosition="left">
          <Center style={{ width: '100%' }}>
            <Text color="dimmed">
              {omicsIds.length === 1
                ? 'The gene is differentially expressed in the studies below.'
                : 'At least one of the genes is differentially expressed in the studies below.'}
            </Text>
          </Center>
          <Center style={{ width: '100%' }} mb="sm">
            <Text color="dimmed">The half volcano plot highlights the gene among other differentially expressed genes in the same study and annotation.</Text>
          </Center>
          {studySummaries.map((s) => (
            <Accordion.Item value={s.studyName} key={s.studyId}>
              <Accordion.Control>
                <StudySummaryLabel studySummary={s} />
              </Accordion.Control>
              <Accordion.Panel>
                <StudyMarkerGeneDetails studySummary={s} />
              </Accordion.Panel>
            </Accordion.Item>
          ))}
        </Accordion>
      )}
    </Container>
  );
}

export default MarkerGeneSearch;
