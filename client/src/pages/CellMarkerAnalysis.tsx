import { useEffect } from 'react';
import { Divider, Group, Loader, Space, Stack, Text, useMantineTheme } from '@mantine/core';
import { useRecoilState, useRecoilValue } from 'recoil';
import { annotationGroupIdState, selectedAnnotationState, selectedGenesState, studyState } from '../atoms';
import { ProjectionPlot } from '../components/ProjectionPlot/ProjectionPlot';
import { useExpressionValues } from '../hooks';
import { ProjectionSelectBox } from '../components/ProjectionSelectBox/ProjectionSelectBox';
import { AnnotationGroupSelectBox } from '../components/AnnotationGroupSelectBox/AnnotationGroupSelectBox';
import { LeftSidePanel } from '../components/LeftSidePanel/LeftSidePanel';
import { AnnotationGroupDisplay } from '../components/AnnotationGroupDisplay/AnnotationGroupDisplay';
import { RightSidePanel } from '../components/RightSidePanel/RightSidePanel';
import { UserGeneStore } from '../components/UserGeneStore/UserGeneStore';
import { DEGTable } from '../components/DEGTable/DEGTable';
import { StudyTitle } from '../components/StudyTitle/StudyTitle';

// const ANNOTATIONS = [
//     {label: "bone cell", color: "#1f77b4"},
//     {label: "chondrocyte", color: "#ff7f0e"},
//     {label: "endothelial cell", color: "#2ca02c"},
//     {label: "endothelial cell of artery", color: "#d62728"},
//     {label: "fibroblast", color: "#9467bd"},
//     {label: "mesenchymal stem cell", color: "#8c564b"},
//     {label: "pericyte cell", color: "#e377c2"}
// ]

// interface PreparedPlot {
//     message?: string;
//     plotlyData: Partial<Plotly.PlotData>[];
//     plotlyLayout: Partial<Plotly.Layout>;
// }

function ProjectionPlotWithOptionalExpression() {
  const theme = useMantineTheme();
  const selectedGenes = useRecoilValue(selectedGenesState);
  const { table, loading } = useExpressionValues(
    selectedGenes.map((g) => g.omicsId),
    true,
  );

  if (loading) {
    return (
      <div>
        <Loader variant="dots" color={theme.colors.gray[5]} size="xl" />
      </div>
    );
  }

  return (
    <div style={{ width: '100%', height: '100%', position: 'relative' }}>
      {table && <ProjectionPlot colorBy="annotation" expressionTable={table} />}
      {table === undefined && <ProjectionPlot colorBy="annotation" />}
    </div>
  );
}

function DifferentialExpressionAnalysis() {
  const annotationGroupId = useRecoilValue(annotationGroupIdState);
  const selectedAnnotation = useRecoilValue(selectedAnnotationState);
  const [selectedGenes, setSelectedGenes] = useRecoilState(selectedGenesState);
  const study = useRecoilValue(studyState);

  useEffect(() => {
    if (selectedGenes.length > 1) setSelectedGenes(selectedGenes.slice(0, 1));
  }, [selectedGenes, setSelectedGenes]);

  if (!study) {
    return null;
  }
  return (
    <Group style={{ height: '100vh' }} align="flex-start" position="apart" spacing="xs" noWrap>
      <LeftSidePanel>
        <Stack>
          <ProjectionSelectBox />
          {annotationGroupId && <AnnotationGroupSelectBox />}
          {annotationGroupId && study.annotationGroupMap.get(annotationGroupId)?.differentialExpressionCalculated ? null : (
            <Text color="red" size="xs">
              No DEGs calculated for selected group.
            </Text>
          )}
          {annotationGroupId && <AnnotationGroupDisplay />}
        </Stack>
      </LeftSidePanel>
      <main style={{ width: '100%', height: '100%' }}>
        <ProjectionPlotWithOptionalExpression />
      </main>
      <RightSidePanel>
        <Stack>
          <StudyTitle />
          <Divider size="xs" label="User gene store" />
          <UserGeneStore multiple={false} />
          <Space h="xs" />
          <Stack>
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
            {selectedAnnotation !== undefined && <DEGTable annotationId={selectedAnnotation} />}
          </Stack>
        </Stack>
      </RightSidePanel>
    </Group>
  );
}

export default DifferentialExpressionAnalysis;
