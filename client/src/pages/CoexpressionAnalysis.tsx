import { useEffect } from 'react';
import { Center, Divider, Group, Loader, Space, Stack, Text, useMantineTheme } from '@mantine/core';
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
import { StudyTitle } from '../components/StudyTitle/StudyTitle';

// interface PreparedPlot {
//     message?: string;
//     allSameSampleExprValues: ColumnTable[];
//     combinations: { dimX: number; dimY: number }[];
//     plotlyData: Partial<Plotly.PlotData>[];
//     plotlyLayout: Partial<Plotly.Layout>;
// }

function CoexpressionAnalysisPlot() {
  const theme = useMantineTheme();
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
      <Center style={{ height: '100%', width: '100%' }}>
        <Loader variant="dots" color={theme.colors.gray[5]} size="xl" />
      </Center>
    );
  }
  if (!data?.correlationTrianglePlot) {
    return (
      <Center style={{ height: '100%', width: '100%' }}>
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
    <Center style={{ height: '100%', width: '100%' }}>
      <img
        style={{
          width: '100%',
          height: selectedGenes.length > 3 ? '100%' : '',
          objectFit: 'fill',
          overflow: 'hidden',
        }}
        alt="correlation triangle plot"
        src={data.correlationTrianglePlot}
      />
    </Center>
  );
}

function CoexpressionAnalysis() {
  const study = useRecoilValue(studyState);
  const setOpened = useSetRecoilState(userGeneStoreOpenState);
  const correlationOmicsId = useRecoilValue(correlationOmicsIdState);
  useEffect(() => {
    setOpened(true);
  }, []);
  if (!study) {
    return (
      <Center style={{ height: '100%', width: '100%' }}>
        <Loader variant="dots" color="gray" />{' '}
      </Center>
    );
  }

  return (
    <Group style={{ height: '100vh' }} align="flex-start" position="apart" spacing="xs" noWrap>
      <LeftSidePanel>
        <Stack pt={5}>
          <AnnotationFilterDisplay />
          <Divider />
          <Text size={'xs'} color={'gray'}>
            Pearson correlation is displayed in each pair-wise scatter plot.
          </Text>
        </Stack>
      </LeftSidePanel>
      <main style={{ height: '100vh' }}>
        <CoexpressionAnalysisPlot />
      </main>
      <RightSidePanel>
        <Stack>
          <StudyTitle />
          <Divider size="xs" label="User gene store" />
          <UserGeneStore multiple findCoexpressors />
          <Space />
          <Divider size="xs" label="Correlated genes" />
          {correlationOmicsId === undefined && (
            <Text size="xs" color="gray">
              No correlation exploration triggered yet.
            </Text>
          )}
          {correlationOmicsId !== undefined && <CorrelationTable omicsId={correlationOmicsId} studyId={study.studyId} />}
        </Stack>
      </RightSidePanel>
    </Group>
  );
}

export default CoexpressionAnalysis;
