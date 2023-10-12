import { useEffect } from 'react';
import { Center, createStyles, Divider, Group, Loader, Stack, Text } from '@mantine/core';
import { useRecoilValue, useSetRecoilState } from 'recoil';
import {
  correlationOmicsIdState,
  selectedAnnotationFilterState,
  selectedGenesState,
  studyIdState,
  studyLayerIdState,
  studyState,
  userGeneStoreOpenState,
} from '../atoms';
import { useExpressionCorrelationTrianglePlotQuery } from '../generated/types';
import { AnnotationFilterDisplay } from '../components/AnnotationFilterDisplay/AnnotationFilterDisplay';
import { RightSidePanel } from '../components/RightSidePanel/RightSidePanel';
import { UserGeneStore } from '../components/UserGeneStore/UserGeneStore';
import { CorrelationTable } from '../components/CorrelationTable/CorrelationTable';
import { LeftSidePanel } from '../components/LeftSidePanel/LeftSidePanel';

const useStyles = createStyles(() => ({
  coexpressionPlotImg: {
    width: '100%',
    height: '75%',
    objectFit: 'fill',
    overflow: 'hidden',
  },
}));

function CoexpressionAnalysisPlot() {
  const { classes } = useStyles();
  const selectedGenes = useRecoilValue(selectedGenesState);
  const studyId = useRecoilValue(studyIdState);
  const studyLayerId = useRecoilValue(studyLayerIdState);
  const annotationFilter = useRecoilValue(selectedAnnotationFilterState);
  const setOpened = useSetRecoilState(userGeneStoreOpenState);
  useEffect(() => {
    setOpened(true);
  }, []);

  const { data, loading } = useExpressionCorrelationTrianglePlotQuery({
    variables: {
      studyId,
      studyLayerId,
      omicsIds: selectedGenes.map((g) => g.omicsId),
      excludeAnnotationValueIds: annotationFilter,
    },
    skip: selectedGenes.length < 2,
  });

  if (loading) {
    return (
      <Center w="100%" h="100%">
        <Loader variant="dots" color="blue" size="xl" />
      </Center>
    );
  }
  if (!data?.correlationTrianglePlot) {
    return (
      <Center w="100%" h="100%">
        <Text color="dimmed" size="md">
          Please select at least 2 genes from the{' '}
          <Text span weight={800}>
            gene store
          </Text>
          .
        </Text>
      </Center>
    );
  }
  return (
    <Center w="100%" h="100%">
      <img alt="correlation triangle plot" src={data.correlationTrianglePlot} className={classes.coexpressionPlotImg} />
    </Center>
  );
}

export default function CoexpressionAnalysis() {
  const study = useRecoilValue(studyState);
  const setOpened = useSetRecoilState(userGeneStoreOpenState);
  const correlationOmicsId = useRecoilValue(correlationOmicsIdState);
  useEffect(() => {
    setOpened(true);
  }, []);

  if (!study) {
    return (
      <Center w="100%" h="100%">
        <Loader variant="dots" color="blue" />{' '}
      </Center>
    );
  }

  return (
    <Group h="100%" w="100%" position="apart" spacing="xs" noWrap>
      <LeftSidePanel>
        <Stack pt={5}>
          <AnnotationFilterDisplay />
          <Divider />
          <Text size="xs" color="gray">
            Pearson correlation is displayed in each pair-wise scatter plot.
          </Text>
        </Stack>
      </LeftSidePanel>
      <Stack h="100%" w="100%">
        <CoexpressionAnalysisPlot />
      </Stack>
      <RightSidePanel>
        <UserGeneStore multiple findCoexpressors />
        <Divider size="xs" label="Correlated genes" />
        {correlationOmicsId === undefined && (
          <Text size="xs" color="gray">
            No correlation exploration triggered yet.
          </Text>
        )}
        {correlationOmicsId !== undefined && <CorrelationTable omicsId={correlationOmicsId} studyId={study.studyId} />}
      </RightSidePanel>
    </Group>
  );
}
