import { lazy, Suspense, useMemo } from 'react';
import { View, VisualizationSpec } from 'react-vega';
import { DotPlotElementFragment } from '../../generated/types';
import { ScenegraphEvent } from 'vega';
import { Center, Loader } from '@mantine/core';

function createSpec(xAxis: 'studyName' | 'displaySymbol') {
  return {
    $schema: 'https://vega.github.io/schema/vega-lite/v5.json',
    data: { name: 'table' },
    transform: [
      {
        impute: 'exprSamplesFraction',
        key: 'annotationDisplayValue',
        value: 0,
        groupby: [xAxis],
      },
    ],
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
    },
    layer: [
      {
        mark: { type: 'point', filled: true },
        encoding: {
          size: {
            field: 'exprSamplesFraction',
            type: 'quantitative',
            scale: { domain: [0.0, 1.0] },
            title: ['Expr. fraction', 'x: not measured'],
          },
          color: {
            field: 'median',
            type: 'quantitative',
            scale: { scheme: 'viridis', reverse: true },
            title: 'Median',
          },
        },
      },
      {
        mark: 'text',
        encoding: {
          text: {
            condition: { test: "datum['exprSamplesFraction'] === 0", value: 'x' },
            value: '',
          },
          color: {
            value: 'lightgray',
          },
        },
      },
    ],
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
  const VegaLite = lazy(() => import('react-vega/lib/VegaLite'));
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
    <Suspense
      fallback={
        <Center style={{ height: '100%', width: '100%' }}>
          <Loader variant={'dots'} color={'gray'} />
        </Center>
      }
    >
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
    </Suspense>
  );
}
