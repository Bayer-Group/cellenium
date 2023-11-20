import { useCallback, useEffect, useMemo, useState } from 'react';
import { ActionIcon, Button, Center, ColorSwatch, Container, createStyles, Divider, Group, Loader, SimpleGrid, Stack, Table, Text, Title } from '@mantine/core';
import { useRecoilValue, useSetRecoilState } from 'recoil';
import { Params, Struct } from 'arquero/dist/types/table/transformable';
import { IconTable, IconX, IconInfoCircle } from '@tabler/icons-react';
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
import {
  StudyInfoFragment,
  useExpressionByAnnotationQuery,
  useExpressionGroupTableQuery,
  useExpressionTTestLazyQuery,
  useExpressionViolinPlotQuery,
  useGeneSpecificityStudyQuery,
} from '../generated/types';
import { ExpressionDotPlot } from '../components/ExpressionDotPlot/ExpressionDotPlot';
import { ProjectionSelectBox } from '../components/ProjectionSelectBox/ProjectionSelectBox';
import { AnnotationGroupSelectBox, AnnotationSecondGroupSelectBox } from '../components/AnnotationGroupSelectBox/AnnotationGroupSelectBox';
import { LeftSidePanel } from '../components/LeftSidePanel/LeftSidePanel';
import { AnnotationGroupDisplay } from '../components/AnnotationGroupDisplay/AnnotationGroupDisplay';
import { AnnotationFilterDisplay } from '../components/AnnotationFilterDisplay/AnnotationFilterDisplay';
import { RightSidePanel } from '../components/RightSidePanel/RightSidePanel';
import { UserGeneStore } from '../components/UserGeneStore/UserGeneStore';
import { GeneSpecificityPlot } from '../components/GeneSpecificityPlot/GeneSpecificityPlot';
import { CrossStudySelector } from '../components/CrossStudySelector/CrossStudySelector';
import { GlobalLoading } from './GlobalLoading';
import { StudyInfoModal } from '../components/StudyTitle/StudyTitle';

const analysisTypes = [
  { value: 'violinplot', label: 'Violin Plot' },
  { value: 'projection', label: 'Projection Plot' },
  { value: 'dotplot', label: 'Dot Plot' },
  { value: 'specificity', label: 'Gene Specificity Plot' },
  /*
            {value: 'boxplot', label: 'Boxplot'},
                {value: 'dot', label: 'Dotplot'},
  
             */
];

const useStyles = createStyles((theme) => ({
  violinPlotImg: { display: 'block', marginLeft: 'auto', marginRight: 'auto' },
  violinPlotWrapper: { overflowX: 'auto', overflowY: 'hidden' },
  expressionAnalysisMain: { overflowX: 'hidden', overflowY: 'auto' },
  active: {
    backgroundColor: theme.colors.blue[6],
  },
  main: {
    '&:hover': {
      backgroundColor: theme.colors.blue[4],
    },
  },
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
      <Stack align="flex-start" pos="relative" px="md">
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

function ExpressionTable({ omicsId }: { omicsId: number }) {
  const study = useRecoilValue(studyState);
  const studyId = useRecoilValue(studyIdState);
  const studyLayerId = useRecoilValue(studyLayerIdState);
  const annotationGroupId = useRecoilValue(annotationGroupIdState);
  const annotationSecondaryGroupId = useRecoilValue(annotationSecondaryGroupIdState);
  const annotationFilter = useRecoilValue(selectedAnnotationFilterState);
  const { data, loading } = useExpressionGroupTableQuery({
    variables: {
      studyId,
      studyLayerId,
      omicsId,
      annotationGroupId: annotationGroupId || -1,
      annotationSecondaryGroupId: annotationSecondaryGroupId || annotationGroupId || -1,
      excludeAnnotationValueIds: annotationFilter,
    },
    skip: !annotationGroupId || !studyId,
  });
  const sortedData = useMemo(
    () =>
      data?.expressionByTwoAnnotationsList &&
      [...data.expressionByTwoAnnotationsList].sort((a, b) =>
        a.secondAnnotationDisplayValue !== b.secondAnnotationDisplayValue
          ? a.secondAnnotationDisplayValue.localeCompare(b.secondAnnotationDisplayValue)
          : a.annotationDisplayValue.localeCompare(b.annotationDisplayValue),
      ),
    [data],
  );

  const annotations = study?.annotationGroupMap.get(annotationGroupId || -1)?.annotationValuesList;
  const [selectedTTestIndexes, setSelectedTTestIndexes] = useState<number[]>([]);
  const [tTestResult, setTTestResult] = useState<string | undefined>();
  const [runTTest, { loading: ttestLoading }] = useExpressionTTestLazyQuery();
  useEffect(() => {
    setSelectedTTestIndexes([]);
    setTTestResult(undefined);
  }, [sortedData]);

  const onTTestClick = (index: number) => {
    if (selectedTTestIndexes.length === 0) {
      setSelectedTTestIndexes([index]);
    } else {
      if (!sortedData) {
        return;
      }

      runTTest({
        variables: {
          studyId,
          studyLayerId,
          omicsId,
          annotationGroupId: annotationGroupId || -1,
          // @ts-ignore
          secondAnnotationGroupId: annotationSecondaryGroupId,
          excludeAnnotationValueIds: annotationFilter,
          sample1AnnotationValueId: sortedData[selectedTTestIndexes[0]].annotationValueId,
          sample1SecondAnnotationValueId: sortedData[selectedTTestIndexes[0]].secondAnnotationValueId,
          sample2AnnotationValueId: sortedData[index].annotationValueId,
          sample2SecondAnnotationValueId: sortedData[index].secondAnnotationValueId,
        },
      }).then((r) => {
        setSelectedTTestIndexes([selectedTTestIndexes[0], index]);
        setTTestResult(r.data?.expressionTtest);
      });
    }
  };

  if (sortedData) {
    return (
      <Container>
        <div style={{ width: annotationSecondaryGroupId ? '70em' : '45em' }}>
          {tTestResult && (
            <Text italic color="dimmed">
              t test: no correction for multiple testing or technical covariables is implemented.
            </Text>
          )}
          <Table>
            <thead>
              <tr>
                <th>Annotation Group</th>
                {annotationSecondaryGroupId && <th>Second Group</th>}
                <th>Median</th>
                <th title="amount of cells in this annotation group"># in group</th>
                <th title="amount of cells with gene expression measured (i.e. no dropout)"># expressed</th>
                <th title="test probability of both groups having the same expression level">Run t test</th>
              </tr>
            </thead>
            <tbody>
              {sortedData.map((r, index) => (
                <tr key={`${r.annotationValueId}_${r.secondAnnotationValueId}`}>
                  <td>
                    <Group spacing="xs">
                      <ColorSwatch color={annotations?.find((a) => a.annotationValueId === r.annotationValueId)?.color || ''} size={12} />
                      <Text> {r.annotationDisplayValue}</Text>
                    </Group>
                  </td>
                  {annotationSecondaryGroupId && <td>{r.secondAnnotationDisplayValue}</td>}
                  <td>{r.median.toFixed(2)}</td>
                  <td>{r.valueCount}</td>
                  <td>{r.nonZeroValueCount}</td>
                  <td>
                    {tTestResult === undefined && !ttestLoading && (
                      <Button size="xs" h="1.5em" disabled={selectedTTestIndexes.indexOf(index) !== -1} variant="default" onClick={() => onTTestClick(index)}>
                        Group {selectedTTestIndexes.length > 0 && selectedTTestIndexes.indexOf(index) === -1 ? 2 : 1}
                      </Button>
                    )}
                    {selectedTTestIndexes.indexOf(index) !== -1 && (
                      <>
                        {ttestLoading && <Loader variant="dots" color="blue" />}
                        {tTestResult && (
                          <Group spacing="xs">
                            <Text title="probability of both groups having the same expression level">p={tTestResult}</Text>
                            <ActionIcon
                              size="xs"
                              onClick={() => {
                                setSelectedTTestIndexes([]);
                                setTTestResult(undefined);
                              }}
                            >
                              <IconX />
                            </ActionIcon>
                          </Group>
                        )}
                      </>
                    )}
                  </td>
                </tr>
              ))}
            </tbody>
          </Table>
        </div>
      </Container>
    );
  }
  return <Container>{loading && <Loader variant="dots" color="blue" />}</Container>;
}

function ViolinPlotAndTable({ displaySymbol, omicsId }: { displaySymbol: string; omicsId: number }) {
  const { cx, classes } = useStyles();
  const [showTable, setShowTable] = useState(false);
  return (
    <>
      <Group>
        <Title order={3}>{displaySymbol}</Title>
        <ActionIcon
          className={cx(classes.main, {
            [classes.active]: showTable,
          })}
          onClick={() => setShowTable((p) => !p)}
          title="Cell count table"
        >
          <IconTable color={showTable ? 'white' : 'gray'} />
        </ActionIcon>
      </Group>
      <Stack w="100%" className={classes.violinPlotWrapper}>
        <ViolinPlot omicsId={omicsId} />
      </Stack>
      {showTable && <ExpressionTable omicsId={omicsId} />}
    </>
  );
}

function ViolinPlots() {
  const selectedGenes = useRecoilValue(selectedGenesState);

  return (
    <>
      {[...selectedGenes].reverse().map((g) => (
        <Stack key={g.omicsId} align="center" w="100%">
          <ViolinPlotAndTable displaySymbol={g.displaySymbol} omicsId={g.omicsId} />
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

function GeneSpecificity() {
  const study = useRecoilValue(studyState);
  const studyLayerId = useRecoilValue(studyLayerIdState);
  const annotationGroupId = useRecoilValue(annotationGroupIdState);
  const annotationSecondaryGroupId = useRecoilValue(annotationSecondaryGroupIdState);
  const annotationFilter = useRecoilValue(selectedAnnotationFilterState);
  const selectedGenes = useRecoilValue(selectedGenesState);
  const [secondSelectedStudy, setSelectedStudy] = useState<number | undefined>(undefined);
  const [modalOpen, setModalOpen] = useState(false);

  const [minMaxOne, setMinMaxOne] = useState<[number, number, number, number] | undefined>(undefined);
  const [minMaxTwo, setMinMaxTwo] = useState<[number, number, number, number] | undefined>(undefined);

  const minMaxXY = useMemo(() => {
    if (minMaxOne && minMaxTwo) {
      return {
        minX: Math.min(minMaxOne[0], minMaxTwo[0]),
        maxX: Math.max(minMaxOne[1], minMaxTwo[1]),
        minY: Math.min(minMaxOne[2], minMaxTwo[2]),
        maxY: Math.max(minMaxOne[3], minMaxTwo[3]),
      };
    }
    return undefined;
  }, [minMaxOne, minMaxTwo]);

  const modalClick = useCallback(() => {
    setModalOpen(!modalOpen);
  }, [modalOpen]);

  const { data: secondStudyData, loading } = useGeneSpecificityStudyQuery({
    variables: {
      studyId: secondSelectedStudy ?? 0,
    },
    skip: !secondSelectedStudy,
  });

  const setMinMaxCallback = useCallback((minX: number, maxX: number, minY: number, maxY: number) => setMinMaxOne([minX, maxX, minY, maxY]), [setMinMaxOne]);
  const setMinMaxTwoCallback = useCallback((minX: number, maxX: number, minY: number, maxY: number) => setMinMaxTwo([minX, maxX, minY, maxY]), [setMinMaxOne]);

  return (
    selectedGenes.length > 0 &&
    annotationSecondaryGroupId &&
    study && (
      <SimpleGrid cols={1} h="100%" w="100%" breakpoints={[{ minWidth: '150rem', cols: 2 }]}>
        {secondStudyData?.referenceStudyOverviewsList[0] && (
          <StudyInfoModal opened={modalOpen} onClose={modalClick} study={secondStudyData?.referenceStudyOverviewsList[0] as StudyInfoFragment} />
        )}
        <GeneSpecificityPlot
          studyName={study?.studyName || ''}
          study_id={study?.studyId ?? 0}
          studyLayerId={studyLayerId}
          omicsIds={selectedGenes.map((e) => e.omicsId)}
          annotationGroupId={annotationGroupId ?? 0}
          secondAnnotationGroupId={annotationSecondaryGroupId ?? 0}
          excludeAnnotationValueIds={annotationFilter}
          showAnnotationLegend={false}
          reportMinMaxXY={setMinMaxCallback}
          minMaxXY={minMaxXY}
        />
        {secondSelectedStudy === undefined ? (
          <CrossStudySelector selectStudy={setSelectedStudy} taxId={study.organismTaxId} />
        ) : loading || !secondStudyData ? (
          <GlobalLoading />
        ) : (
          <Stack pos="relative" w="100%" h="100%">
            <Group spacing="xs" pos="absolute" top={20} right={20} style={{ zIndex: 190 }}>
              <ActionIcon size="lg" onClick={modalClick}>
                <IconInfoCircle height="100%" width="auto" />
              </ActionIcon>
              <ActionIcon size="lg" onClick={() => setSelectedStudy(undefined)}>
                <IconX height="100%" width="auto" />
              </ActionIcon>
            </Group>
            <GeneSpecificityPlot
              study_id={secondSelectedStudy}
              omicsIds={selectedGenes.map((e) => e.omicsId)}
              studyName={secondStudyData?.referenceStudyOverviewsList[0].studyName || ''}
              studyLayerId={secondStudyData?.referenceStudyOverviewsList[0].defaultStudyLayerId || 0}
              annotationGroupId={secondStudyData?.referenceStudyOverviewsList[0].referenceStudyInfoList[0].tissueAnnotationGroupId}
              secondAnnotationGroupId={secondStudyData?.referenceStudyOverviewsList[0].referenceStudyInfoList[0].celltypeAnnotationGroupId}
              excludeAnnotationValueIds={[]}
              reportMinMaxXY={setMinMaxTwoCallback}
              minMaxXY={minMaxXY}
            />
          </Stack>
        )}
      </SimpleGrid>
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
        <ExpressionAnalysisTypeSelectBox handleSelection={setAnalysisType as (v: unknown) => void} selection={analysisType} options={analysisTypes} />
        {analysisType === 'projection' && <ProjectionSelectBox />}
        {(analysisType === 'violinplot' || analysisType === 'dotplot' || analysisType === 'specificity') && (
          <>
            <Stack spacing={0}>
              {analysisType === 'specificity' && (
                <Text size="xs" weight="bold">
                  Please select a tissue annotation group in this dropdown
                </Text>
              )}
              <AnnotationGroupSelectBox />
              {((analysisType === 'violinplot' && annotationGroupId && annotationSecondaryGroupId) || analysisType === 'specificity') && (
                <Stack spacing="sm" pt="sm">
                  <Text size="xs">{analysisType === 'violinplot' ? 'Violins' : 'Bubbles'} are colored by annotation</Text>
                  <AnnotationGroupDisplay disableSelection />
                </Stack>
              )}
            </Stack>
            <Stack spacing={0}>
              {analysisType === 'specificity' && (
                <Text size="xs" weight="bold">
                  Please select a celltype annotation group in this dropdown
                </Text>
              )}
              {(analysisType === 'violinplot' || analysisType === 'specificity') && <AnnotationSecondGroupSelectBox />}
            </Stack>
            <Divider my="sm" />
            <AnnotationFilterDisplay />
          </>
        )}
      </LeftSidePanel>
      <Stack h="100%" w="100%" pos="relative" py="md" spacing="md" className={cx(classes.expressionAnalysisMain, 'no-scrollbar')}>
        {analysisType === 'violinplot' && <ViolinPlots />}
        {analysisType === 'projection' && <ProjectionPlots />}
        {analysisType === 'dotplot' && <DotPlots />}
        {analysisType === 'specificity' && <GeneSpecificity />}
        {analysisType === 'specificity' && (selectedGenes.length === 0 || annotationSecondaryGroupId === undefined) && (
          <Center h="100%" w="100%">
            <Text c="dimmed">
              Please select a gene from the&nbsp;
              <Text span weight={800}>
                gene store
              </Text>
              &nbsp;and both&nbsp;
              <Text span weight={800}>
                annotation groups
              </Text>
              .
            </Text>
          </Center>
        )}
        {analysisType !== 'specificity' && selectedGenes.length === 0 && (
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
        <UserGeneStore multiple={analysisType !== 'specificity'} />
      </RightSidePanel>
    </Group>
  );
}

export default ExpressionAnalysis;
