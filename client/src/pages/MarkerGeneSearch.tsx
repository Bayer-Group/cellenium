import { useMemo, useState } from 'react';
import { Anchor, Badge, Card, Center, createStyles, Grid, Group, Loader, SimpleGrid, Stack, Text } from '@mantine/core';
import { Link } from 'react-router-dom';
import { DifferentialMarkerFragment, useStudiesWithMarkerGenesQuery } from '../generated/types';
import { NavBarProvider } from '../components/NavBar/NavBar';
import { GeneSearchBar } from '../components/SearchBar/GeneSearchBar';
import { MarkerCard } from '../components/MarkerCard/MarkerCard';

interface StudySummary {
  studyId: number;
  studyName: string;
  cellCount: number;
  description: string;
  differentialExpressionsList: DifferentialMarkerFragment[];
}

const useStyles = createStyles(() => ({
  studyName: { textOverflow: 'ellipsis', overflow: 'hidden' },
  diffExpGrid: { maxHeight: '20rem', overflowY: 'scroll', justifyItems: 'center' },
}));

function StudyMarkerGeneCard({ study }: { study: StudySummary }) {
  const { classes, cx } = useStyles();
  const newStudyUrl = `/study/${study.studyId}`;

  return (
    <Card shadow="sm" p="lg" radius="md" withBorder>
      <Card.Section withBorder inheritPadding py="xs">
        <Grid columns={12}>
          <Grid.Col span={8}>
            <Anchor component={Link} to={newStudyUrl} color="dark">
              <Text align="left" lineClamp={1} className={classes.studyName} weight={800}>
                {study.studyName}
              </Text>
            </Anchor>
          </Grid.Col>
          <Grid.Col span={4}>
            <Group position="right">
              <Badge variant="light" color="gray">
                {Math.round(study.cellCount / 1000)}k cells
              </Badge>
              <Badge variant="light" color="gray">
                max log2FC: {Math.max(...study.differentialExpressionsList.map((de) => de.log2Foldchange)).toFixed(2)}
              </Badge>
            </Group>
          </Grid.Col>
        </Grid>
      </Card.Section>
      <Text mt="sm" mb="sm" color="dimmed" size="sm" lineClamp={3} align="left">
        {study.description}
      </Text>
      <Card.Section withBorder inheritPadding py="xs">
        <SimpleGrid
          breakpoints={[
            { maxWidth: 'md', cols: 1 },
            { maxWidth: 'xl', cols: 2 },
            { minWidth: 'xl', cols: 3 },
          ]}
          className={cx(classes.diffExpGrid, 'no-scrollbar')}
          spacing="xs"
        >
          {study.differentialExpressionsList.map((sr: DifferentialMarkerFragment) => (
            <MarkerCard data={sr} key={`${sr.study.studyId}_${sr.annotationValueId}`} />
          ))}
        </SimpleGrid>
      </Card.Section>
    </Card>
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
      const buildStudySummaries: StudySummary[] = [];
      data.differentialExpressionsList.forEach((sr) => {
        let studySummary = buildStudySummaries.find((s) => s.studyId === sr.study.studyId);
        if (!studySummary) {
          studySummary = {
            studyId: sr.study.studyId,
            studyName: sr.study.studyName,
            cellCount: sr.study.cellCount,
            description: sr.study.description,
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
    return undefined;
  }, [data]);

  return (
    <NavBarProvider scrollable>
      <Stack p="md" spacing={0}>
        <GeneSearchBar humanOnly={false} onGeneSelection={(ids: number[]) => setOmicsIds(ids)} />
        <Center w="100%" m="sm">
          {loading && <Loader variant="dots" color="blue" />}
          {omicsIds.length === 0 && (
            <Text color="dimmed">
              Please enter your genes of interest. Cellenium will search for studies and cell annotation clusters that show the entered gene as differentially
              expressed.
            </Text>
          )}
        </Center>

        {studySummaries && studySummaries.length > 0 && (
          <>
            <Center w="100%">
              <Text color="dimmed" align="center">
                {omicsIds.length === 1
                  ? 'The gene is differentially expressed in the studies below.'
                  : 'At least one of the genes is differentially expressed in the studies below.'}
              </Text>
            </Center>
            <Center w="100%" mb="sm">
              <Text align="center" color="dimmed">
                The half volcano plot highlights the gene among other differentially expressed genes in the same study and annotation.
              </Text>
            </Center>
            <Stack>
              {studySummaries.map((s) => (
                <StudyMarkerGeneCard study={s} key={`marke-gene-study-card-${s.studyId}`} />
              ))}
            </Stack>
          </>
        )}
      </Stack>
    </NavBarProvider>
  );
}

export default MarkerGeneSearch;
