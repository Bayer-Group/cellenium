import Plot from 'react-plotly.js';
import * as Plotly from 'plotly.js';
import { Loader, Stack } from '@mantine/core';
import { useEffect, useMemo } from 'react';
import { useGeneSpecificityQuery } from '../../generated/types';

const LAYOUT: Partial<Plotly.Layout> = {
  autosize: true,
  legend: {
    xanchor: 'left',
    yanchor: 'top',
    itemsizing: 'constant',
    x: 1,
    y: 0.8,
  },
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
  reportMinMaxXY,
  showAnnotationLegend = true,
  minMaxXY,
}: {
  studyName: string;
  study_id: number;
  studyLayerId: number;
  omicsIds: number[];
  annotationGroupId: number;
  secondAnnotationGroupId: number;
  excludeAnnotationValueIds: number[] | number;
  showAnnotationLegend?: boolean;
  reportMinMaxXY?: (minX: number, maxX: number, minY: number, maxY: number) => void;
  minMaxXY?: { minX: number; maxX: number; minY: number; maxY: number };
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

  const { traces, layout } = useMemo(() => {
    if (!data) return { trace: undefined, layout: undefined };

    const annotationSet = [...new Set(data.expressionByTwoAnnotationsList.map((e) => e.annotationDisplayValue))];

    const t = annotationSet.map((cs) => {
      const d = data.expressionByTwoAnnotationsList.filter((e) => e.annotationDisplayValue === cs);
      const tTooltips = d.map(
        (e) =>
          `X:${e.annotationDisplayValue}(mean expression)</br></br>Y:${e.secondAnnotationDisplayValue}(${e.exprSamplesFraction}% of samples have expression values)</br></br>Size: ${e.valueCount} samples (log10)`,
      );
      const tRadius = d.map((e) => Math.max(e.valueCount, 5));
      const tX = d.map((e) => e.mean);
      const tY = d.map((e) => e.exprSamplesFraction);
      const colors = d.map((e) => e.color);

      return {
        minX: Math.min(...tX),
        maxX: Math.max(...tX),
        minY: Math.min(...tY),
        maxY: Math.max(...tY),
        trace: {
          x: tX,
          y: tY,
          hoverinfo: 'text',
          name: cs,
          text: tTooltips,
          mode: 'markers',
          type: 'scattergl',
          marker: { size: tRadius, color: colors },
        } as Partial<Plotly.Data>,
      };
    });

    const minX = Math.min(...t.map((e) => e.minX));
    const maxX = Math.max(...t.map((e) => e.maxX));
    const minY = Math.min(...t.map((e) => e.minY));
    const maxY = Math.max(...t.map((e) => e.maxY));

    return {
      traces: [...t.map((e) => e.trace)],
      layout: {
        ...LAYOUT,
        showlegend: showAnnotationLegend,
        xaxis: {
          range: [minX - 1, maxX + 1],
        },
        yaxis: {
          range: [minY - 0.2, maxY + 1],
        },
        title: studyName.length > 70 ? `${studyName.slice(0, 70)}...` : studyName,
      },
    };
  }, [data, showAnnotationLegend, studyName]);

  useEffect(() => {
    if (layout && reportMinMaxXY && !minMaxXY && layout.xaxis?.range && layout.yaxis?.range) {
      reportMinMaxXY(layout.xaxis?.range[0], layout.xaxis?.range[1], layout.yaxis?.range[0], layout.yaxis?.range[1]);
    }
  }, [minMaxXY, reportMinMaxXY, layout]);

  const calcLayout: Partial<Plotly.Layout> = useMemo(() => {
    if (!minMaxXY) {
      return layout || ({} as Partial<Plotly.Layout>);
    }

    return {
      ...layout,
      xaxis: {
        range: [minMaxXY.minX, minMaxXY.maxX],
        title: {
          text: data?.annotationGroupsList.find((e) => e.annotationGroupId === annotationGroupId)?.displayGroup || '',
        },
      },
      yaxis: {
        range: [minMaxXY.minY, minMaxXY.maxY],
        title: {
          text: data?.annotationGroupsList.find((e) => e.annotationGroupId === secondAnnotationGroupId)?.displayGroup || '',
        },
      },
    } as Partial<Plotly.Layout>;
  }, [annotationGroupId, data?.annotationGroupsList, layout, minMaxXY, secondAnnotationGroupId]);

  if (loading || !data) {
    return (
      <Stack w="100%" h="100%" justify="center" align="center">
        <Loader color="blue" />
      </Stack>
    );
  }
  return (
    <Stack w="100%" h="100%" mih="40rem">
      <Plot data={traces || []} layout={calcLayout} config={CONFIG} style={{ width: '100%', height: '100%' }} />
    </Stack>
  );
}
