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
        setStudyId(2)
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
        const distinctAnnotationValueIds: number[] = t.rollup({annotationValueIds: aq.op.array_agg_distinct('annotationValueId')}).array('annotationValueIds')[0];

        // one plotly data track per category, so that we can assign categorical colors
        const plotlyData = distinctAnnotationValueIds.map(annotationValueId => {
            const color = study.annotationValueMap.get(annotationValueId)?.color || '';
            const tableForAnnotation = t.params({annotationValueId}).filter((d: any, p: any) => d.annotationValueId === p.annotationValueId);
            const samplesTrace = {
                type: 'scattergl',
                x: tableForAnnotation.array('projectionX', Float32Array),
                y: tableForAnnotation.array('projectionY', Float32Array),
                customdata: tableForAnnotation.array('annotationValueId', Int32Array),
                mode: 'markers',
                marker: {
                    size: 3,
                    opacity: 0.7,
                    color,
                },
                showlegend: false,
            } as Partial<Plotly.PlotData>;
            if (highlightAnnotation === annotationValueId) {
                console.log('highlight', color, annotationValueId)
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
                            color,
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
                },
                yaxis: {
                    visible: false,
                    fixedrange: true,
                },
            },
        }
    }, [study, highlightAnnotation]);

    const onHover = (event: Readonly<Plotly.PlotHoverEvent>) => {
        if (event.points.length > 0 && event.points[0].customdata) {
            const an = event.points[0].customdata as number;
            setHighlightAnnotation(an);
        }
    };

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
