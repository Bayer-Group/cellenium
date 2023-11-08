import { useCallback, useMemo, useState } from 'react';
import { Center, createStyles, Group, Loader, Space, Stack, Text } from '@mantine/core';
import { useRecoilState, useRecoilValue } from 'recoil';
import { useNavigate } from 'react-router-dom';
import { ScenegraphEvent } from 'vega';
import { DotPlotElementFragment, StudyInfoFragment, useExpressionByAnnotationQuery } from '../generated/types';
import { cellOAnnotationGroupIdState, GeneSearchSelection } from '../atoms';
import { ExpressionDotPlot } from '../components/ExpressionDotPlot/ExpressionDotPlot';
import { StudySearchBar } from '../components/SearchBar/StudySearchBar';
import { GeneSearchBar } from '../components/SearchBar/GeneSearchBar';

const useStyles = createStyles(() => ({
  plotContainer: {
    width: '100%',
    overflowX: 'auto',
    transform: 'rotateX(180deg)', // this is a hack to display the scrollbar on the top
  },
  plot: {
    transform: 'rotateX(180deg)', // this is a hack to display the scrollbar on the top
  },
}));

function CrossStudySearch() {
  const { classes } = useStyles();
  const [studyList, setStudyList] = useState<StudyInfoFragment[]>([]);
  const [omicsIds, setOmicsIds] = useState<number[]>([]);
  const cellOAnnotationGroupId = useRecoilValue(cellOAnnotationGroupIdState);
  const navigate = useNavigate();
  const [selectedFilters] = useRecoilState(GeneSearchSelection);

  const { data, loading } = useExpressionByAnnotationQuery({
    variables: {
      studyLayerIds: studyList?.map((s) => s.defaultStudyLayerId) || [],
      omicsIds,
      annotationGroupId: cellOAnnotationGroupId || -1,
      excludeAnnotationValueIds: [],
    },
    skip: omicsIds.length === 0 || studyList.length === 0 || !cellOAnnotationGroupId,
  });

  const heatmapDisplayData = useMemo(() => {
    if (studyList.length === 0 || !data?.expressionByAnnotationList) {
      return undefined;
    }
    const studyLayerIdMap = new Map(studyList.map((s) => [s.defaultStudyLayerId, s]));
    const allData = data.expressionByAnnotationList.map((o) => ({
      ...o,
      studyName: studyLayerIdMap.get(o.studyLayerId)?.studyName,
      studyId: studyLayerIdMap.get(o.studyLayerId)?.studyId,
    }));
    return omicsIds.map((id) => ({
      omicsId: id,
      heatmapData: allData.filter((d) => d.omicsId === id),
    }));
  }, [studyList, data, omicsIds]);

  const onHeatmapClick = useCallback(
    (dotPlotElement: DotPlotElementFragment, event: ScenegraphEvent) => {
      const dpe = dotPlotElement as DotPlotElementFragment & { studyId: number };
      const newStudyUrl = `/study/${dpe.studyId}?page=CellMarkerAnalysis&annotationGroupId=${cellOAnnotationGroupId}&annotationValueId=${dotPlotElement.annotationValueId}&omicsId=${dotPlotElement.omicsId}`;
      if (event.shiftKey || event.altKey) {
        const parsedUrl = new URL(window.location.href);
        const url = `${parsedUrl.origin}${newStudyUrl}`;
        window.open(url, '_blank');
      } else {
        navigate(newStudyUrl);
      }
    },
    [cellOAnnotationGroupId, navigate],
  );

  const onGeneSelection = useCallback((ids: number[]) => setOmicsIds(ids), []);

  return (
    <Stack p="md" spacing={0}>
      <GeneSearchBar humanOnly onGeneSelection={onGeneSelection} />
      <Space h="xl" />
      <StudySearchBar onStudyListUpdate={setStudyList} />
      <Center w="100%">
        <Text color="dimmed" align="center" mt="1rem">
          Please enter your genes of interest. Cellenium will show the gene&apos;s expression in human studies with standardized cell annotation (CellO). As the
          study data is processed and normalized independently, this is a qualitative direction for which studies to explore independently. Click in the chart
          to open a study.
        </Text>
      </Center>
      <Group position="center" align="start" w="100%">
        {loading && <Loader mt="1rem" variant="dots" color="blue" />}
        {heatmapDisplayData &&
          heatmapDisplayData.map((heatmap) => (
            <Stack key={`${heatmap.omicsId}-expression-dot-plot`} w="100%" align="center" mt="1rem">
              <Text weight="bold">{selectedFilters.find((i) => i.omicsId.includes(heatmap.omicsId))?.displaySymbol}</Text>

              <Stack w="100%" align="center" className={classes.plotContainer}>
                <Stack w="100%" className={classes.plot}>
                  <ExpressionDotPlot data={heatmap.heatmapData} xAxis="studyName" onClick={onHeatmapClick} responsiveHeight={false} responsiveWidth={false} />
                </Stack>
              </Stack>
            </Stack>
          ))}
      </Group>
    </Stack>
  );
}

export default CrossStudySearch;
