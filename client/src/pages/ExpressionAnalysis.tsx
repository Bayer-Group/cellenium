import { useEffect, useMemo, useState } from 'react';
import { Center, Divider, Group, Loader, Stack, Text, Title } from '@mantine/core';
import { useRecoilValue, useSetRecoilState } from 'recoil';
import { Params, Struct } from 'arquero/dist/types/table/transformable';
import { ExpressionAnalysisTypeSelectBox } from '../components/ExpressionAnalysisTypeSelectBox/ExpressionAnalysisTypeSelectBox';
import {
  annotationGroupIdState,
  annotationSecondaryGroupIdState,
  selectedAnnotationFilterState,
  selectedGenesState,
  studyIdState,
  studyLayerIdState,
  studyState,
  userGeneStoreOpenState,
} from '../atoms';
import { ProjectionPlot } from '../components/ProjectionPlot/ProjectionPlot';
import { useExpressionValues } from '../hooks';
import { useExpressionByAnnotationQuery, useExpressionViolinPlotQuery } from '../generated/types';
import { ExpressionDotPlot } from '../components/ExpressionDotPlot/ExpressionDotPlot';
import { ProjectionSelectBox } from '../components/ProjectionSelectBox/ProjectionSelectBox';
import { AnnotationGroupSelectBox, AnnotationSecondGroupSelectBox } from '../components/AnnotationGroupSelectBox/AnnotationGroupSelectBox';
import { LeftSidePanel } from '../components/LeftSidePanel/LeftSidePanel';
import { AnnotationGroupDisplay } from '../components/AnnotationGroupDisplay/AnnotationGroupDisplay';
import { AnnotationFilterDisplay } from '../components/AnnotationFilterDisplay/AnnotationFilterDisplay';
import { RightSidePanel } from '../components/RightSidePanel/RightSidePanel';
import { UserGeneStore } from '../components/UserGeneStore/UserGeneStore';
import { StudyTitle } from '../components/StudyTitle/StudyTitle';

const analysisTypes = [
  { value: 'violinplot', label: 'Violin plot' },
  { value: 'projection', label: 'Projection plot' },
  { value: 'dotplot', label: 'Dot plot' },
  /*
        {value: 'boxplot', label: 'Boxplot'},
            {value: 'dot', label: 'Dotplot'},

         */
];

function ViolinPlot({ omicsId }: { omicsId: number }) {
  const studyId = useRecoilValue(studyIdState);
  const studyLayerId = useRecoilValue(studyLayerIdState);
  const annotationGroupId = useRecoilValue(annotationGroupIdState);
  const annotationSecondaryGroupId = useRecoilValue(annotationSecondaryGroupIdState);
  const annotationFilter = useRecoilValue(selectedAnnotationFilterState);
  const { data, loading } = useExpressionViolinPlotQuery({
    variables: {
      studyId,
      studyLayerId,
      omicsId,
      annotationGroupId: annotationGroupId || -1,
      annotationSecondaryGroupId: annotationSecondaryGroupId as number,
      excludeAnnotationValueIds: annotationFilter,
    },
    skip: !annotationGroupId || !studyId,
  });

  if (data?.violinPlot) {
    return (
      <Stack align="flex-start" h="28rem" pos="relative">
        <img style={{ height: '100%', display: 'block' }} alt="violin plot" src={data.violinPlot} />
      </Stack>
    );
  }
  return (
    <Stack w="100%" h="100%" align="center">
      {loading && <Loader variant="dots" color="blue" size="xl" />}
    </Stack>
  );
}

function ViolinPlots() {
  const selectedGenes = useRecoilValue(selectedGenesState);

  return (
    <Stack align="center" pos="relative" style={{ width: '100%', overflowX: 'hidden' }} bg="white">
      {[...selectedGenes].reverse().map((g) => (
        <Stack key={g.omicsId} align="center" w="100%">
          <Title order={3}>{g.displaySymbol}</Title>
          <Stack h="30rem" w="100%" style={{ overflowX: 'scroll', overflowY: 'hidden' }}>
            <ViolinPlot omicsId={g.omicsId} />
          </Stack>
        </Stack>
      ))}
    </Stack>
  );
}

function ProjectionPlots() {
  const selectedGenes = useRecoilValue(selectedGenesState);
  const { table, loading } = useExpressionValues(
    selectedGenes.map((g) => g.omicsId),
    true,
  );
  const tablePerGene = useMemo(() => {
    if (selectedGenes.length === 0 || !table) {
      return undefined;
    }
    return [...selectedGenes].reverse().map((g) => table.params({ omicsId: g.omicsId }).filter((d: Struct, p: Params) => d.omicsId === p.omicsId));
  }, [selectedGenes, table]);

  if (loading) {
    return (
      <Center style={{ height: '100%', width: '100%' }}>
        <Loader variant="dots" color="blue" size={25} />
      </Center>
    );
  }

  return (
    <Stack w="100%" pt="1rem" spacing="md">
      {tablePerGene &&
        [...selectedGenes].reverse().map((g, i) => (
          <Stack key={g.omicsId} align="center" mih="50vh">
            <Title order={3}>{g.displaySymbol}</Title>
            <ProjectionPlot colorBy="expression" expressionTable={tablePerGene[i]} />
          </Stack>
        ))}
    </Stack>
  );
}

function DotPlots() {
  const study = useRecoilValue(studyState);
  const studyLayerId = useRecoilValue(studyLayerIdState);
  const annotationGroupId = useRecoilValue(annotationGroupIdState);
  const selectedGenes = useRecoilValue(selectedGenesState);
  const annotationFilter = useRecoilValue(selectedAnnotationFilterState);
  const { data, loading } = useExpressionByAnnotationQuery({
    variables: {
      omicsIds: selectedGenes.map((g) => g.omicsId),
      studyLayerIds: [studyLayerId],
      annotationGroupId: annotationGroupId || -1,
      excludeAnnotationValueIds: annotationFilter,
    },
    skip: selectedGenes.length === 0 || !annotationGroupId || !study,
  });
  const heatmapDisplayData = useMemo(() => {
    if (!data?.expressionByAnnotationList) {
      return undefined;
    }
    return data.expressionByAnnotationList.map((dotPlotElement) => ({
      ...dotPlotElement,
      displaySymbol: study?.studyOmicsMap?.get(dotPlotElement.omicsId)?.displaySymbol || 'nn',
    }));
  }, [data?.expressionByAnnotationList, study?.studyOmicsMap]);

  if (loading) {
    return (
      <Center style={{ height: '100%', width: '100%' }}>
        <Loader variant="dots" color="blue" size={25} />
      </Center>
    );
  }
  return (
    heatmapDisplayData && (
      <Stack style={{ height: '100%', width: '100%' }} align="center" justify="center">
        <ExpressionDotPlot data={heatmapDisplayData} xAxis="displaySymbol" />
      </Stack>
    )
  );
}

function ExpressionAnalysis() {
  const [analysisType, setAnalysisType] = useState<string>(analysisTypes[0].value);
  const setOpened = useSetRecoilState(userGeneStoreOpenState);
  useEffect(() => {
    setOpened(() => true);
  }, [setOpened]);

  // const [selectedAnnotationGroup, setSelectedAnnotationGroup] = useState<number>();
  const selectedGenes = useRecoilValue(selectedGenesState);

  const study = useRecoilValue(studyState);
  const annotationGroupId = useRecoilValue(annotationGroupIdState);
  const annotationSecondaryGroupId = useRecoilValue(annotationSecondaryGroupIdState);
  if (!study) {
    return null;
  }
  return (
    <Group align="flex-start" position="apart" spacing="xs" noWrap h="100%" w="100%">
      <LeftSidePanel>
        <Stack spacing="md">
          <ExpressionAnalysisTypeSelectBox handleSelection={setAnalysisType as (v: unknown) => void} selection={analysisType} options={analysisTypes} />
          {(analysisType === 'violinplot' || analysisType === 'dotplot') && (
            <>
              <AnnotationGroupSelectBox />
              {analysisType === 'violinplot' && annotationGroupId && annotationSecondaryGroupId && (
                <>
                  <Text size="xs">Violins are colored by annotation</Text>
                  <AnnotationGroupDisplay disableSelection />
                </>
              )}

              {analysisType === 'violinplot' && <AnnotationSecondGroupSelectBox />}

              <Divider my="sm" />
              <AnnotationFilterDisplay />
            </>
          )}

          {analysisType === 'projection' && <ProjectionSelectBox />}
        </Stack>
      </LeftSidePanel>
      <main style={{ height: '100%', width: '100%', overflowY: 'scroll', flexGrow: 1 }} className="plotContainer">
        {analysisType === 'violinplot' && <ViolinPlots />}
        {analysisType === 'projection' && <ProjectionPlots />}
        {analysisType === 'dotplot' && <DotPlots />}
        {selectedGenes.length === 0 && (
          <Center style={{ height: '100%', width: '100%' }}>
            <Text c="dimmed">
              Please select gene(s) from the{' '}
              <Text span weight={800}>
                gene store
              </Text>
              .
            </Text>
          </Center>
        )}
      </main>
      <RightSidePanel>
        <Stack w="100%">
          <StudyTitle />
          <Divider size="xs" label="User gene store" />
          <UserGeneStore multiple />
        </Stack>
      </RightSidePanel>
    </Group>
  );
}

export default ExpressionAnalysis;
