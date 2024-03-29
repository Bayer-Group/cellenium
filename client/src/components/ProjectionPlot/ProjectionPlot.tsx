import React, { useCallback } from 'react';
import { useRecoilState, useRecoilValue } from 'recoil';
import Plot from 'react-plotly.js';
import * as aq from 'arquero';
import * as Plotly from 'plotly.js';
import { Params, Struct } from 'arquero/dist/types/table/transformable';
import { Center, createStyles, Text } from '@mantine/core';
import { annotationGroupIdState, highlightAnnotationState, selectedAnnotationState, selectedProjectionState, studyState } from '../../atoms';
import { ExpressionTable } from '../../model';

interface PreparedPlot {
  message?: string;
  plotlyData: Partial<Plotly.PlotData>[];
  plotlyLayout: Partial<Plotly.Layout>;
}

const plotlyConfig: Partial<Plotly.Config> = {
  responsive: true,
  displayModeBar: false,
};

const useStyles = createStyles(() => ({
  full: {
    width: '100%',
    height: '100%',
  },
}));

export function ProjectionPlot({
  colorBy,
  expressionTable,
  showSampleIds,
  disableSelection,
}: {
  colorBy: 'annotation' | 'expression';
  expressionTable?: ExpressionTable;
  showSampleIds?: number[] | null;
  disableSelection?: boolean;
}) {
  const { classes } = useStyles();
  const annotationGroupId = useRecoilValue(annotationGroupIdState);

  const study = useRecoilValue(studyState);
  const projection = useRecoilValue(selectedProjectionState);
  const [highlightAnnotation, setHighlightAnnotation] = useRecoilState(highlightAnnotationState);
  const [selectedAnnotation, setSelectedAnnotation] = useRecoilState(selectedAnnotationState);
  const isSelectable = disableSelection ? false : (study?.annotationGroupMap.get(annotationGroupId as number)?.differentialExpressionCalculated as boolean);

  const annotationProjectionData = React.useMemo(() => {
    if (!study || !study.samplesProjectionTables.get(projection)) {
      return undefined;
    }

    let samplesAnnotationProjectionTable = study.samplesAnnotationTable
      .params({ annotationGroupId })
      .filter((d: Struct, p: Params) => d.annotationGroupId === p.annotationGroupId);
    samplesAnnotationProjectionTable = samplesAnnotationProjectionTable
      .join_right(study.samplesProjectionTables.get(projection), 'studySampleId')
      .impute({ annotationValueId: () => -1 })
      .reify();
    const rangeX = [
      aq.agg(samplesAnnotationProjectionTable, aq.op.min('projectionX')) * 1.05,
      aq.agg(samplesAnnotationProjectionTable, aq.op.max('projectionX')) * 1.05,
    ];
    const rangeY = [
      aq.agg(samplesAnnotationProjectionTable, aq.op.min('projectionY')) * 1.05,
      aq.agg(samplesAnnotationProjectionTable, aq.op.max('projectionY')) * 1.05,
    ];
    const distinctAnnotationValueIds: number[] = samplesAnnotationProjectionTable
      .rollup({
        annotationValueIds: aq.op.array_agg_distinct('annotationValueId'),
      })
      .array('annotationValueIds')[0];

    return {
      samplesAnnotationProjectionTable,
      rangeX,
      rangeY,
      distinctAnnotationValueIds,
    };
  }, [annotationGroupId, study, projection]);

  // one plotly data trace per category, so that we can assign categorical colors
  const annotationTraces = React.useMemo(() => {
    if (!study || !annotationProjectionData || colorBy !== 'annotation') {
      return undefined;
    }

    // the cells in selected annotation color
    return annotationProjectionData.distinctAnnotationValueIds.map((annotationValueId) => {
      const tableForAnnotation = annotationProjectionData.samplesAnnotationProjectionTable
        .params({ annotationValueId })
        .filter((d: Struct, p: Params) => d.annotationValueId === p.annotationValueId);
      return {
        type: 'scattergl',
        x: tableForAnnotation.array('projectionX', Float32Array),
        y: tableForAnnotation.array('projectionY', Float32Array),
        customdata: study.annotationValueMap.get(annotationValueId) && tableForAnnotation.array('annotationValueId', Int32Array),
        text: study.annotationValueMap.get(annotationValueId)?.displayValue,
        mode: 'markers',
        marker: {
          size: 3,
          opacity: expressionTable || (showSampleIds && showSampleIds.length > 0) ? 0.2 : 0.7,
          color: study.annotationValueMap.get(annotationValueId)?.color || '#d7d5d5',
        },
        showlegend: false,
        hoverinfo: 'text',
      } as Partial<Plotly.PlotData>;
    });
  }, [study, annotationProjectionData, expressionTable, colorBy, showSampleIds]);

  // the hovered cells highlighted
  const annotationHighlightTrace = React.useMemo(() => {
    // Don't show the hovered cluster highlight in the expression view - the small expressionTrace points
    // get overlayed with the annotationHighlightTrace points a lot, causing an unHover event, which
    // causes flickering and web browser halt for large datasets.
    if (!study || !annotationProjectionData || highlightAnnotation === 0 || colorBy !== 'annotation') {
      return undefined;
    }
    const tableForAnnotation = annotationProjectionData.samplesAnnotationProjectionTable
      .params({ highlightAnnotation })
      .filter((d: Struct, p: Params) => d.annotationValueId === p.highlightAnnotation);
    return {
      type: 'scattergl',
      x: tableForAnnotation.array('projectionX', Float32Array),
      y: tableForAnnotation.array('projectionY', Float32Array),
      mode: 'markers',
      text: study.annotationValueMap.get(highlightAnnotation)?.displayValue,
      marker: {
        size: 10,
        opacity: 0.1,
        color: 'yellow',
      },
      showlegend: false,
      hoverinfo: 'text',
    } as Partial<Plotly.PlotData>;
  }, [study, annotationProjectionData, highlightAnnotation, colorBy]);

  // the selected cells highlighted
  const selectedAnnotationHighlightTrace = React.useMemo(() => {
    // Don't show the hovered cluster highlight in the expression view - the small expressionTrace points
    // get overlayed with the annotationHighlightTrace points a lot, causing an unHover event, which
    // causes flickering and web browser halt for large datasets.
    if (!study || !annotationProjectionData || selectedAnnotation === 0 || colorBy !== 'annotation') {
      return undefined;
    }
    const tableForAnnotation = annotationProjectionData.samplesAnnotationProjectionTable
      .params({ selectedAnnotation })
      .filter((d: Struct, p: Params) => d.annotationValueId === p.selectedAnnotation);
    return {
      type: 'scattergl',
      x: tableForAnnotation.array('projectionX', Float32Array),
      y: tableForAnnotation.array('projectionY', Float32Array),
      mode: 'markers',
      text: study.annotationValueMap.get(selectedAnnotation as number)?.displayValue,
      marker: {
        size: 10,
        opacity: 0.1,
        color: colorBy === 'annotation' ? study.annotationValueMap.get(selectedAnnotation)?.color : '#dddddd',
      },
      showlegend: false,
      hoverinfo: 'text',
    } as Partial<Plotly.PlotData>;
  }, [study, annotationProjectionData, selectedAnnotation, colorBy]);

  // the cells colored according to gene expression for the expression analysis page, include 0 values
  const expressionTrace = React.useMemo(() => {
    if (!study || !expressionTable || !annotationProjectionData || colorBy !== 'expression') {
      return undefined;
    }
    // right join: we want to plot all cells, even those without expression data for this gene (set to 0 for plotting)
    const joinedTable = expressionTable
      .join_right(annotationProjectionData.samplesAnnotationProjectionTable, 'studySampleId')
      .impute({ value: () => 0 })
      .orderby(['value'])
      .reify();
    return {
      type: 'scattergl',
      x: joinedTable.array('projectionX', Float32Array),
      y: joinedTable.array('projectionY', Float32Array),
      customdata: joinedTable.array('annotationValueId', Int32Array),
      mode: 'markers',
      marker: {
        size: 3,
        opacity: 0.7,
        color: joinedTable.array('value', Float32Array),
        colorscale: 'Viridis',
        reversescale: true,
        colorbar: {
          tickfont: {
            size: '14',
          },
          thickness: 15,
          outlinewidth: 0,
          len: 0.2,
          lenmode: 'fraction',
          yanchor: 'top',
          y: 1,
          yref: 'paper',
        },
      },
      showlegend: false,
      hoverinfo: 'none',
    } as Partial<Plotly.PlotData>;
  }, [study, colorBy, annotationProjectionData, expressionTable]);

  // gene expression of selected gene projected on top of normal annotation plot
  const selectedGeneExpressionTrace = React.useMemo(() => {
    if (!study || !expressionTable || !annotationProjectionData || colorBy !== 'annotation') {
      return undefined;
    }
    const joinedTable = expressionTable.join(annotationProjectionData.samplesAnnotationProjectionTable, 'studySampleId').orderby(['value']).reify();
    return {
      type: 'scattergl',
      x: joinedTable.array('projectionX', Float32Array),
      y: joinedTable.array('projectionY', Float32Array),
      customdata: joinedTable.array('annotationValueId', Int32Array),
      mode: 'markers',
      marker: {
        size: 5,
        opacity: 1,
        color: joinedTable.array('value', Float32Array),
        colorscale: 'Viridis',
        cmin: 0,
        cmax: Math.max(...joinedTable.array('value')),
        reversescale: true,
        colorbar: {
          tickfont: {
            size: '14',
          },
          thickness: 15,
          outlinewidth: 0,
          len: 0.2,
          lenmode: 'fraction',
          yanchor: 'top',
          y: 1,
          yref: 'paper',
        },
      },
      showlegend: false,
      hoverinfo: 'none',
    } as Partial<Plotly.PlotData>;
  }, [study, annotationProjectionData, expressionTable, colorBy]);

  // traces from the interactive user cell selection
  const selectedSamplesTrace = React.useMemo(() => {
    if (!study || !annotationProjectionData || !showSampleIds) {
      return undefined;
    }
    const filteredTable = annotationProjectionData.samplesAnnotationProjectionTable
      .params({ showSampleIds })
      .filter((d: Struct, p: Params) => aq.op.includes(p.showSampleIds, d.studySampleId, 0))
      .reify();
    return {
      type: 'scattergl',
      x: filteredTable.array('projectionX', Float32Array),
      y: filteredTable.array('projectionY', Float32Array),
      customdata: filteredTable.array('annotationValueId', Int32Array),
      mode: 'markers',
      marker: {
        size: 10,
        opacity: 0.8,
        color: 'yellow',
      },
      showlegend: false,
      hoverinfo: 'none',
    } as Partial<Plotly.PlotData>;
  }, [annotationProjectionData, study, showSampleIds]);

  const preparedPlot: PreparedPlot | undefined = React.useMemo(() => {
    if (!annotationProjectionData || (!annotationTraces && !expressionTrace)) {
      return undefined;
    }
    const plotlyData = [
      ...(annotationTraces || []),
      ...(annotationHighlightTrace ? [annotationHighlightTrace] : []),
      ...(selectedAnnotationHighlightTrace ? [selectedAnnotationHighlightTrace] : []),
      ...(expressionTrace ? [expressionTrace] : []),
      ...(selectedGeneExpressionTrace ? [selectedGeneExpressionTrace] : []),
      ...(selectedSamplesTrace ? [selectedSamplesTrace] : []),
    ];
    return {
      plotlyData,
      plotlyLayout: {
        // width: 850,
        // height: window.innerHeight - 50,
        margin: { l: 0, r: 0, t: 0, b: 0 },
        xaxis: {
          visible: false,
          fixedrange: true,
          range: annotationProjectionData.rangeX,
        },
        yaxis: {
          visible: false,
          fixedrange: true,
          range: annotationProjectionData.rangeY,
        },
      },
    };
  }, [
    annotationProjectionData,
    annotationTraces,
    expressionTrace,
    annotationHighlightTrace,
    selectedAnnotationHighlightTrace,
    selectedGeneExpressionTrace,
    selectedSamplesTrace,
  ]);

  const onHover = useCallback(
    (event: Readonly<Plotly.PlotHoverEvent>) => {
      if (event.points.length > 0 && event.points[0].customdata) {
        const annotationValueId = event.points[0].customdata as number;
        setHighlightAnnotation(annotationValueId);
      }
    },
    [setHighlightAnnotation],
  );

  const onClick = useCallback(
    (event: Readonly<Plotly.PlotMouseEvent>) => {
      if (!isSelectable) {
        return;
      }
      if (event.points.length > 0 && event.points[0].customdata) {
        const annotationValueId = event.points[0].customdata as number;
        setSelectedAnnotation(annotationValueId);
      }
    },
    [isSelectable, setSelectedAnnotation],
  );

  const onDoubleClick = useCallback(() => {
    setSelectedAnnotation(0);
  }, [setSelectedAnnotation]);

  const onUnHover = useCallback(() => {
    setHighlightAnnotation(0);
  }, [setHighlightAnnotation]);
  if (preparedPlot) {
    return (
      <Plot
        data={preparedPlot.plotlyData}
        layout={preparedPlot.plotlyLayout}
        config={plotlyConfig}
        onHover={onHover}
        onClick={onClick}
        onDoubleClick={onDoubleClick}
        onUnhover={onUnHover}
        className={classes.full}
      />
    );
  }
  return (
    <Center w="100%" h="100%">
      <Text weight="bold" color="gray">
        No plot available
      </Text>
    </Center>
  );
}
