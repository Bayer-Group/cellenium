import { useMemo } from 'react';
import { VegaLite, View, VisualizationSpec } from 'react-vega';
import { DotPlotElementFragment } from '../../generated/types';
import { ScenegraphEvent } from 'vega';

function createSpec(xAxis: 'studyName' | 'displaySymbol') {
  return {
    $schema: 'https://vega.github.io/schema/vega-lite/v5.json',
    data: { name: 'table' },
    mark: { type: 'point', filled: true },
    encoding: {
      y: {
        field: 'annotationDisplayValue',
        type: 'ordinal',
        title: '',
        axis: { offset: 10 },
      },
      x:
        xAxis === 'studyName'
          ? {
              field: 'studyName',
              type: 'nominal',
              title: '',
              axis: { offset: 10 },
            }
          : {
              field: 'displaySymbol',
              type: 'nominal',
              title: '',
              axis: { offset: 10 },
            },
      size: {
        field: 'exprSamplesFraction',
        type: 'quantitative',
        scale: { domain: [0.0, 1.0] },
        title: 'Expr. fraction',
      },
      color: {
        field: 'median',
        type: 'quantitative',
        scale: { scheme: 'viridis', reverse: true },
      },
    },
  } as VisualizationSpec;
}

export function ExpressionDotPlot({
  data,
  annotationTitle,
  xAxis,
  onClick,
}: {
  data: DotPlotElementFragment[];
  annotationTitle: string;
  xAxis: 'studyName' | 'displaySymbol';
  onClick?: (dotPlotElement: DotPlotElementFragment, event: ScenegraphEvent) => void;
}) {
  const spec = useMemo(() => createSpec(xAxis), [annotationTitle, xAxis]);

  const setUpSelectionListener = (view: View) => {
    view.addEventListener('click', (event, item) => {
      if (item && onClick) {
        const dotPlotElement = item.datum as DotPlotElementFragment;
        onClick(dotPlotElement, event);
      }
    });
  };

  return (
    <VegaLite
      spec={spec}
      onNewView={(view) => setUpSelectionListener(view)}
      config={{
        view: { stroke: 'transparent' },
      }}
      actions={false}
      data={{
        table: data,
      }}
    />
  );
}
