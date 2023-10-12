import { useCallback, useMemo, useState } from 'react';
import { Center, Group, Loader, Space, Stack, Text } from '@mantine/core';
import { useRecoilValue } from 'recoil';
import { useNavigate } from 'react-router-dom';
import { ScenegraphEvent } from 'vega';
import { DotPlotElementFragment, StudyInfoFragment, useExpressionByAnnotationQuery } from '../generated/types';
import { allGenesState, cellOAnnotationGroupIdState } from '../atoms';
import { ExpressionDotPlot } from '../components/ExpressionDotPlot/ExpressionDotPlot';
import { StudySearchBar } from '../components/SearchBar/StudySearchBar';
import { NavBarProvider } from '../components/NavBar/NavBar';
import { GeneSearchBar } from '../components/SearchBar/GeneSearchBar';

function CrossStudySearch() {
  const [studyList, setStudyList] = useState<StudyInfoFragment[]>([]);
  const [omicsIds, setOmicsIds] = useState<number[]>([]);
  const cellOAnnotationGroupId = useRecoilValue(cellOAnnotationGroupIdState);
  const allGenes = useRecoilValue(allGenesState);
  const navigate = useNavigate();

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

  return (
    <NavBarProvider scrollable>
      <Stack p="md" spacing={0}>
        <GeneSearchBar humanOnly onGeneSelection={(ids) => setOmicsIds(ids)} />
        <Space h="xl" />
        <StudySearchBar onStudyListUpdate={setStudyList} />
        <Center w="100%">
          <Text color="dimmed" align="center" mt="1rem">
            Please enter your genes of interest. Cellenium will show the gene&apos;s expression in human studies with standardized cell annotation (CellO). As
            the study data is processed and normalized independently, this is a qualitative direction for which studies to explore independently. Click in the
            chart to open a study.
          </Text>
        </Center>
        <Group position="center" align="start" w="100%">
          {loading && <Loader mt="1rem" variant="dots" color="blue" size={25} />}
          {heatmapDisplayData &&
            heatmapDisplayData.map((heatmap) => (
              <Stack key={`${heatmap.omicsId}-expression-dot-plot`} w="100%" align="center" mt="1rem">
                <Text weight="bold">{allGenes?.get(heatmap.omicsId)?.displaySymbol}</Text>
                <ExpressionDotPlot data={heatmap.heatmapData} xAxis="studyName" onClick={onHeatmapClick} responsiveHeight={false} />
              </Stack>
            ))}
        </Group>
      </Stack>
    </NavBarProvider>
  );
}

export default CrossStudySearch;
