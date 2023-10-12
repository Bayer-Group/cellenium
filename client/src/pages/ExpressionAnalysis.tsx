import { useEffect, useMemo, useState } from 'react';
import { Center, createStyles, Divider, Group, Loader, Stack, Text, Title } from '@mantine/core';
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

const analysisTypes = [
  { value: 'violinplot', label: 'Violin Plot' },
  { value: 'projection', label: 'Projection Plot' },
  { value: 'dotplot', label: 'Dot Plot' },
  /*
        {value: 'boxplot', label: 'Boxplot'},
            {value: 'dot', label: 'Dotplot'},

         */
];

const useStyles = createStyles(() => ({
  violinPlotImg: { height: '100%', display: 'block', marginLeft: 'auto', marginRight: 'auto' },
  violinPlotWrapper: { overflowX: 'auto', overflowY: 'hidden' },
  expressionAnalysisMain: { overflowX: 'hidden', overflowY: 'auto' },
}));

function ViolinPlot({ omicsId }: { omicsId: number }) {
  const { classes } = useStyles();

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
      <Stack align="flex-start" h="28rem" pos="relative" px="md">
        <img alt="violin plot" src={data.violinPlot} className={classes.violinPlotImg} />
      </Stack>
    );
  }
  return (
    <Stack w="100%" h="100%" align="center">
      {loading && <Loader variant="dots" color="blue" />}
    </Stack>
  );
}

function ViolinPlots() {
  const { classes } = useStyles();
  const selectedGenes = useRecoilValue(selectedGenesState);

  return (
    <>
      {[...selectedGenes].reverse().map((g) => (
        <Stack key={g.omicsId} align="center" w="100%">
          <Title order={3}>{g.displaySymbol}</Title>
          <Stack h="30rem" w="100%" className={classes.violinPlotWrapper}>
            <ViolinPlot omicsId={g.omicsId} />
          </Stack>
        </Stack>
      ))}
    </>
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
      <Center h="100%" w="100%">
        <Loader variant="dots" color="blue" />
      </Center>
    );
  }

  return (
    tablePerGene &&
    [...selectedGenes].reverse().map((g, i) => (
      <Stack key={g.omicsId} align="center" mih="50vh">
        <Title order={3}>{g.displaySymbol}</Title>
        <ProjectionPlot colorBy="expression" expressionTable={tablePerGene[i]} />
      </Stack>
    ))
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
      <Center h="100%" w="100%">
        <Loader variant="dots" color="blue" />
      </Center>
    );
  }
  return (
    heatmapDisplayData && (
      <Stack h="100%" w="100%" align="center" justify="center" p="md">
        <ExpressionDotPlot data={heatmapDisplayData} xAxis="displaySymbol" responsiveHeight />
      </Stack>
    )
  );
}

function ExpressionAnalysis() {
  const { classes, cx } = useStyles();
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
    <Group h="100%" w="100%" position="apart" spacing={0} noWrap>
      <LeftSidePanel>
        {analysisType === 'projection' && <ProjectionSelectBox />}
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
      </LeftSidePanel>
      <Stack h="100%" w="100%" pos="relative" py="md" spacing="md" className={cx(classes.expressionAnalysisMain, 'no-scrollbar')}>
        {analysisType === 'violinplot' && <ViolinPlots />}
        {analysisType === 'projection' && <ProjectionPlots />}
        {analysisType === 'dotplot' && <DotPlots />}
        {selectedGenes.length === 0 && (
          <Center h="100%" w="100%">
            <Text c="dimmed">
              Please select gene(s) from the&nbsp;
              <Text span weight={800}>
                gene store
              </Text>
              .
            </Text>
          </Center>
        )}
      </Stack>
      <RightSidePanel>
        <UserGeneStore multiple />
      </RightSidePanel>
    </Group>
  );
}

export default ExpressionAnalysis;
