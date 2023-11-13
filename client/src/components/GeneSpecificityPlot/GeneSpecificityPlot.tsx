import Plot from 'react-plotly.js';
import * as Plotly from 'plotly.js';
import { Stack } from '@mantine/core';
import { useMemo } from 'react';
import { useGeneSpecificityQuery } from '../../generated/types';
import { GlobalLoading } from '../../pages/GlobalLoading';

const LAYOUT = {
  autosize: true,
};

const CONFIG = {
  responsive: true,
};

export function GeneSpecificityPlot({
  study_id,
  studyName,
  studyLayerId,
  omicsIds,
  annotationGroupId,
  secondAnnotationGroupId,
  excludeAnnotationValueIds,
}: {
  studyName: string;
  study_id: number;
  studyLayerId: number;
  omicsIds: number[];
  annotationGroupId: number;
  secondAnnotationGroupId: number;
  excludeAnnotationValueIds: number[] | number;
}) {
  const { data, loading } = useGeneSpecificityQuery({
    variables: {
      studyId: study_id,
      studyLayerId,
      omicsIds,
      annotationGroupId,
      secondAnnotationGroupId,
      excludeAnnotationValueIds,
    },
  });

  const { trace, layout } = useMemo(() => {
    if (!data) return { trace: undefined, layout: undefined };
    const colorSet = [...new Set(data.expressionByTwoAnnotationsList.map((e) => e.annotationDisplayValue))];
    const tColors = data.expressionByTwoAnnotationsList.map((e) => colorSet.indexOf(e.annotationDisplayValue));
    const tTooltips = data.expressionByTwoAnnotationsList.map((e) => `${e.annotationDisplayValue}:${e.secondAnnotationDisplayValue}`);
    const tRadius = data.expressionByTwoAnnotationsList.map((e) => e.valueCount);
    const tX = data.expressionByTwoAnnotationsList.map((e) => e.mean);
    const tY = data.expressionByTwoAnnotationsList.map((e) => e.exprSamplesFraction);
    return {
      trace: {
        x: tX,
        y: tY,
        text: tTooltips,
        mode: 'markers',
        type: 'scattergl',
        marker: { size: tRadius, color: tColors },
      } as Partial<Plotly.Data>,
      layout: {
        ...LAYOUT,
        xaxis: {
          range: [-1, Math.max(...tX) + 1],
        },
        yaxis: {
          range: [-0.2, Math.max(...tY) + 1],
        },
        title: studyName,
      },
    };
  }, [data]);

  if (loading || !data) {
    return (
      <Stack w="100%" h="100%">
        <GlobalLoading />
      </Stack>
    );
  }

  return <Plot data={trace ? [trace] : []} layout={layout || {}} config={CONFIG} style={{ width: '100%' }} />;
}
