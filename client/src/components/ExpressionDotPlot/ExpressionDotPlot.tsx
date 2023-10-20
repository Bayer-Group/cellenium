import { lazy, Suspense, useMemo } from 'react';
import { View, VisualizationSpec } from 'react-vega';
import { ScenegraphEvent } from 'vega';
import { Center, createStyles, Loader } from '@mantine/core';
import { DotPlotElementFragment } from '../../generated/types';

const useStyles = createStyles(() => ({
  plot: {
    height: '100%',
    width: '100%',
    maxHeight: '75%',
    margin: 'auto',
    '& canvas': {
      margin: 'auto',
      display: 'block!important',
    },
  },
}));

function createSpec(xAxis: 'studyName' | 'displaySymbol', responsiveHeight: boolean, responsiveWidth: boolean) {
  return {
    $schema: 'https://vega.github.io/schema/vega-lite/v5.json',
    data: { name: 'table' },
    width: responsiveWidth ? 'container' : undefined,
    height: responsiveHeight ? 'container' : undefined,
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
              axis: { offset: 10, orient: 'top' },
            }
          : {
              field: 'displaySymbol',
              type: 'nominal',
              title: '',
              axis: { offset: 10, orient: 'top' },
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
  xAxis,
  onClick,
  responsiveHeight = false,
  responsiveWidth = true,
}: {
  data: DotPlotElementFragment[];
  xAxis: 'studyName' | 'displaySymbol';
  onClick?: (dotPlotElement: DotPlotElementFragment, event: ScenegraphEvent) => void;
  responsiveHeight?: boolean;
  responsiveWidth?: boolean;
}) {
  const { classes } = useStyles();
  const VegaLite = lazy(() => import('react-vega/lib/VegaLite'));
  const spec = useMemo(() => createSpec(xAxis, responsiveHeight, responsiveWidth), [responsiveHeight, responsiveWidth, xAxis]);

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
        <Center h="100%" w="100%">
          <Loader variant="dots" color="blue" />
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
        className={classes.plot}
      />
    </Suspense>
  );
}
