import { useEffect } from 'react';
import { Divider, Group, Loader, Stack, Text } from '@mantine/core';
import { useRecoilState, useRecoilValue } from 'recoil';
import { annotationGroupIdState, selectedAnnotationState, selectedDEGComparisonAnnotationState, selectedGenesState, studyState } from '../atoms';
import { ProjectionPlot } from '../components/ProjectionPlot/ProjectionPlot';
import { useExpressionValues } from '../hooks';
import { ProjectionSelectBox } from '../components/ProjectionSelectBox/ProjectionSelectBox';
import { AnnotationGroupSelectBox } from '../components/AnnotationGroupSelectBox/AnnotationGroupSelectBox';
import { LeftSidePanel } from '../components/LeftSidePanel/LeftSidePanel';
import { AnnotationGroupDisplay } from '../components/AnnotationGroupDisplay/AnnotationGroupDisplay';
import { RightSidePanel } from '../components/RightSidePanel/RightSidePanel';
import { UserGeneStore } from '../components/UserGeneStore/UserGeneStore';
import { DEGTable } from '../components/DEGTable/DEGTable';

function ProjectionPlotWithOptionalExpression() {
  const selectedGenes = useRecoilValue(selectedGenesState);
  const { table, loading } = useExpressionValues(
    selectedGenes.map((g) => g.omicsId),
    true,
  );

  if (loading) {
    return (
      <Stack w="100%" h="100%" align="center" justify="center">
        <Loader variant="dots" color="blue" size="xl" />
      </Stack>
    );
  }

  return (
    <Stack w="100%" h="100%" align="center" justify="center">
      {table && <ProjectionPlot colorBy="annotation" expressionTable={table} />}
      {table === undefined && <ProjectionPlot colorBy="annotation" />}
    </Stack>
  );
}

function DifferentialExpressionAnalysis() {
  const annotationGroupId = useRecoilValue(annotationGroupIdState);
  const selectedAnnotation = useRecoilValue(selectedAnnotationState);
  const [selectedGenes, setSelectedGenes] = useRecoilState(selectedGenesState);
  const selectedDEGComparison = useRecoilValue(selectedDEGComparisonAnnotationState);
  const study = useRecoilValue(studyState);
  const pairwiseDifferentialExpressionCalculated =
    annotationGroupId && study?.annotationGroupMap.get(annotationGroupId)?.pairwiseDifferentialExpressionCalculated;

  useEffect(() => {
    if (selectedGenes.length > 1) setSelectedGenes(selectedGenes.slice(0, 1));
  }, [selectedGenes, setSelectedGenes]);

  if (!study) {
    return null;
  }
  return (
    <Group h="100%" w="100%" position="apart" spacing="xs" noWrap>
      <LeftSidePanel>
        <ProjectionSelectBox />
        {annotationGroupId && <AnnotationGroupSelectBox />}
        {annotationGroupId && study.annotationGroupMap.get(annotationGroupId)?.differentialExpressionCalculated ? null : (
          <Text color="red" size="xs">
            No DEGs calculated for selected group.
          </Text>
        )}
        {pairwiseDifferentialExpressionCalculated && (
          <Text size="xs">Pairwise attributes differentially expressed genes available, hold Shift to select a second annotation.</Text>
        )}
        {annotationGroupId && <AnnotationGroupDisplay enableDEGComparisonSelection={pairwiseDifferentialExpressionCalculated === true} />}
      </LeftSidePanel>
      <Stack w="100%" h="100%">
        <ProjectionPlotWithOptionalExpression />
      </Stack>
      <RightSidePanel>
        <UserGeneStore multiple={false} />
        <Stack w="100%" mt="xs">
          <Divider size="xs" label="Differential expression table" />
          {annotationGroupId && study.annotationGroupMap.get(annotationGroupId)?.differentialExpressionCalculated ? null : (
            <Text color="red" size="xs">
              No DEGs calculated for selected group.
            </Text>
          )}
          {!selectedAnnotation && study.annotationGroupMap.get(annotationGroupId as number)?.differentialExpressionCalculated === true && (
            <Text size="xs" color="dimmed">
              Select cells in the plot or via the selection panel on the left-hand side.
            </Text>
          )}
          {selectedAnnotation !== undefined && <DEGTable annotationId={selectedAnnotation} selectedDEGComparisonAnnotationId={selectedDEGComparison} />}
        </Stack>
      </RightSidePanel>
    </Group>
  );
}

export default DifferentialExpressionAnalysis;
