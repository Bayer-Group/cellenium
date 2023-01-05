import React, {useEffect, useState} from 'react';
import {useRecoilState, useRecoilValue} from "recoil";
import {studyIdState, studyState} from "../../atoms";
import Plot from 'react-plotly.js';
import * as aq from 'arquero';
import * as Plotly from "plotly.js";
import ColumnTable from "arquero/dist/types/table/column-table";

interface PreparedPlot {
    message?: string;
    plotlyData: Partial<Plotly.PlotData>[];
    plotlyLayout: Partial<Plotly.Layout>;
}

type Props = {
    colorBy: 'annotation' | 'expression';
    expressionTable?: ColumnTable;
}

const ProjectionPlot = ({
                            colorBy,
                            expressionTable
                        }: Props) => {
    // TODO from recoil
    const selectedAnnotationGroupId = 1;
    const selectedOmicsIds = [116, 264];

    const study = useRecoilValue(studyState);
    const [highlightAnnotation, setHighlightAnnotation] = useState(0);


    const annotationProjectionData = React.useMemo(() => {
        if (!study) {
            return undefined;
        }

        // @ts-ignore
        let samplesAnnotationProjectionTable = study.samplesAnnotationTable.params({selectedAnnotationGroupId}).filter((d, p) => d.annotationGroupId === p.selectedAnnotationGroupId);
        samplesAnnotationProjectionTable = samplesAnnotationProjectionTable.join(study.samplesProjectionTable, 'studySampleId').reify();
        const rangeX = [aq.agg(samplesAnnotationProjectionTable, aq.op.min('projectionX')) * 1.05, aq.agg(samplesAnnotationProjectionTable, aq.op.max('projectionX')) * 1.05];
        const rangeY = [aq.agg(samplesAnnotationProjectionTable, aq.op.min('projectionY')) * 1.05, aq.agg(samplesAnnotationProjectionTable, aq.op.max('projectionY')) * 1.05];
        const distinctAnnotationValueIds: number[] = samplesAnnotationProjectionTable.rollup({annotationValueIds: aq.op.array_agg_distinct('annotationValueId')}).array('annotationValueIds')[0];

        return {
            samplesAnnotationProjectionTable,
            rangeX,
            rangeY,
            distinctAnnotationValueIds
        }
    }, [selectedAnnotationGroupId, study]);

    // one plotly data trace per category, so that we can assign categorical colors
    const annotationTraces = React.useMemo(() => {
        if (!study || !annotationProjectionData || colorBy !== 'annotation') {
            return undefined;
        }
        return annotationProjectionData.distinctAnnotationValueIds.map(annotationValueId => {
            const tableForAnnotation = annotationProjectionData.samplesAnnotationProjectionTable.params({annotationValueId}).filter((d: any, p: any) => d.annotationValueId === p.annotationValueId);
            return {
                type: 'scattergl',
                x: tableForAnnotation.array('projectionX', Float32Array),
                y: tableForAnnotation.array('projectionY', Float32Array),
                customdata: tableForAnnotation.array('annotationValueId', Int32Array),
                text: study.annotationValueMap.get(annotationValueId)?.displayValue,
                mode: 'markers',
                marker: {
                    size: 3,
                    opacity: 0.7,
                    color: study.annotationValueMap.get(annotationValueId)?.color,
                },
                showlegend: false,
                hoverinfo: "text"
            } as Partial<Plotly.PlotData>;
        });
    }, [study, annotationProjectionData, colorBy]);

    const annotationHighlightTrace = React.useMemo(() => {
        if (!study || !annotationProjectionData || highlightAnnotation === 0) {
            return undefined;
        }
        const tableForAnnotation = annotationProjectionData.samplesAnnotationProjectionTable.params({highlightAnnotation}).filter((d: any, p: any) => d.annotationValueId === p.highlightAnnotation);
        return {
            type: 'scattergl',
            x: tableForAnnotation.array('projectionX', Float32Array),
            y: tableForAnnotation.array('projectionY', Float32Array),
            mode: 'markers',
            text: study.annotationValueMap.get(highlightAnnotation)?.displayValue,
            marker: {
                size: 30,
                opacity: 0.02,
                color: colorBy === 'annotation'
                    ? study.annotationValueMap.get(highlightAnnotation)?.color
                    : '#dddddd',
            },
            showlegend: false,
            hoverinfo: "text"
        } as Partial<Plotly.PlotData>;
    }, [study, annotationProjectionData, highlightAnnotation]);

    const expressionTrace = React.useMemo(() => {
        if (!study || !expressionTable || !annotationProjectionData || colorBy !== 'expression') {
            return undefined;
        }
        const joinedTable = expressionTable.join(annotationProjectionData.samplesAnnotationProjectionTable, 'studySampleId').reify();
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
                colorscale: 'YlGnBu',
                reversescale: true,
                colorbar: {
                    thickness: 5,
                },
            },
            showlegend: false,
            hoverinfo: 'none'
        } as Partial<Plotly.PlotData>;
    }, [study, colorBy, annotationProjectionData, expressionTable]);

    const preparedPlot: PreparedPlot | undefined = React.useMemo(() => {
        if (!annotationProjectionData || (!annotationTraces && !expressionTrace)) {
            return undefined;
        }
        const plotlyData = [
            ...(annotationTraces || []),
            ...(annotationHighlightTrace ? [annotationHighlightTrace] : []),
            ...(expressionTrace ? [expressionTrace] : [])];
        return {
            plotlyData,
            plotlyLayout: {
                width: 850,
                height: 700,
                margin: {l: 0, r: 0, t: 0, b: 0},
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
        }
    }, [annotationProjectionData, annotationTraces, annotationHighlightTrace, expressionTrace]);

    const onHover = (event: Readonly<Plotly.PlotHoverEvent>) => {
        if (event.points.length > 0 && event.points[0].customdata) {
            const annotationValueId = event.points[0].customdata as number;
            setHighlightAnnotation(annotationValueId);
        }
    };

    if (preparedPlot) {
        return (<Plot data={preparedPlot.plotlyData}
                      layout={preparedPlot.plotlyLayout}
                      onHover={onHover}/>);
    } else {
        return (<div>no plot</div>)
    }
};

export default ProjectionPlot;
