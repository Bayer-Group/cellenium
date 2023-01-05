import React, {useEffect, useState} from 'react';
import {LeftSidePanel, RightSidePanel} from "../components";
import {Group, Stack} from "@mantine/core";
import {useRecoilState, useRecoilValue} from "recoil";
import {studyIdState, studyState} from "../atoms";
import {useStudyBasicsQuery} from "../generated/types";
import {useExpressionValues} from "../hooks";
import Plot from 'react-plotly.js';
import * as aq from 'arquero';
import * as Plotly from "plotly.js";

interface PreparedPlot {
    message?: string;
    plotlyData: Partial<Plotly.PlotData>[];
    plotlyLayout: Partial<Plotly.Layout>;
}

const AnnotationsInUmapScatterplotTestPage = () => {
    const [studyId, setStudyId] = useRecoilState(studyIdState);
    useEffect(() => {
        setStudyId(1)
    });
    const study = useRecoilValue(studyState);
    const [highlightAnnotation, setHighlightAnnotation] = useState(0);

    const preparedPlot: PreparedPlot | undefined = React.useMemo(() => {
        if (!study) {
            return undefined;
        }

        // @ts-ignore
        let t = study.samplesAnnotationTable.filter(r => r.annotationGroupId === 1);
        t = t.join(study.samplesProjectionTable, 'studySampleId').reify();
        const rangeX = [aq.agg(t, aq.op.min('projectionX')) * 1.05, aq.agg(t, aq.op.max('projectionX')) * 1.05];
        const rangeY = [aq.agg(t, aq.op.min('projectionY')) * 1.05, aq.agg(t, aq.op.max('projectionY')) * 1.05];
        const distinctAnnotationValueIds: number[] = t.rollup({annotationValueIds: aq.op.array_agg_distinct('annotationValueId')}).array('annotationValueIds')[0];

        // one plotly data track per category, so that we can assign categorical colors
        const plotlyData = distinctAnnotationValueIds.map(annotationValueId => {
            const tableForAnnotation = t.params({annotationValueId}).filter((d: any, p: any) => d.annotationValueId === p.annotationValueId);
            const samplesTrace = {
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
            if (highlightAnnotation === annotationValueId) {
                return [
                    samplesTrace,
                    {
                        type: 'scattergl',
                        x: tableForAnnotation.array('projectionX', Float32Array),
                        y: tableForAnnotation.array('projectionY', Float32Array),
                        mode: 'markers',
                        marker: {
                            size: 30,
                            opacity: 0.02,
                            color: study.annotationValueMap.get(annotationValueId)?.color,
                        },
                        showlegend: false,
                    } as Partial<Plotly.PlotData>
                ];
            } else {
                return [samplesTrace];
            }
        }).flat(2);


        return {
            plotlyData,
            plotlyLayout: {
                width: 850,
                height: 700,
                margin: {l: 0, r: 0, t: 0, b: 0},
                xaxis: {
                    visible: false,
                    fixedrange: true,
                    range: rangeX,
                },
                yaxis: {
                    visible: false,
                    fixedrange: true,
                    range: rangeY,
                },
            },
        }
    }, [study, highlightAnnotation]);

    const onHover = (event: Readonly<Plotly.PlotHoverEvent>) => {
        if (event.points.length > 0 && event.points[0].customdata) {
            const annotationValueId = event.points[0].customdata as number;
            setHighlightAnnotation(annotationValueId);
        }
    };

    console.log(preparedPlot?.plotlyLayout)

    return (
        <Group position={'apart'}>
            <LeftSidePanel/>
            <Stack>
                {preparedPlot && <Plot data={preparedPlot.plotlyData}
                                       layout={preparedPlot.plotlyLayout}
                                       onHover={onHover}
                />}
            </Stack>
            <RightSidePanel/>
        </Group>
    );
};

export default AnnotationsInUmapScatterplotTestPage;
