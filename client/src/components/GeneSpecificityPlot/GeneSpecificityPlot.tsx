import Plot from 'react-plotly.js';
import * as Plotly from 'plotly.js';
import { Loader, Stack, Text } from '@mantine/core';
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
      const tTooltips = d.map((e) => {
        return `${e.annotationDisplayValue}:${e.secondAnnotationDisplayValue}</br></br>X: ${
          e.mean > 0.0099 ? e.mean.toFixed(2) : e.mean
        } mean expression</br>Y: ${((e.nonZeroValueCount / e.valueCount) * 100).toFixed(2)}% of samples have expression values</br>Size: ${
          e.nonZeroValueCount
        } samples`;
      });
      const tRadius = d.map((e) => Math.max(e.nonZeroValueCount, 5));
      const tX = d.map((e) => e.mean);
      const tY = d.map((e) => (e.nonZeroValueCount / e.valueCount) * 100);
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
          range: [minX - maxX * 0.05, maxX + 1],
        },
        yaxis: {
          range: [minY - maxY * 0.05, maxY + 1],
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
          text: 'Normalized Mean Gene Expression',
        },
      },
      yaxis: {
        range: [minMaxXY.minY, minMaxXY.maxY],
        title: {
          text: 'The percentage of cells where the gene is detected (PCT %)',
        },
      },
    } as Partial<Plotly.Layout>;
  }, [layout, minMaxXY]);

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
      {omicsIds.length == 0 && (
        <Stack align="center" justify="center" mt="4rem" h="100%" w="100%" pos="absolute" bg="rgb(1,1,1, 0.25)" style={{ backdropFilter: 'blur(2px)' }}>
          <Text color="white" weight="bold" size="lg">
            No Gene Selected
          </Text>
        </Stack>
      )}
    </Stack>
  );
}
