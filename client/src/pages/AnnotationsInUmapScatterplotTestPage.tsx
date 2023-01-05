import React, {useEffect} from 'react';
import {LeftSidePanel, RightSidePanel} from "../components";
import {Group, Stack} from "@mantine/core";
import {useRecoilState, useRecoilValue} from "recoil";
import {studyIdState, studyState} from "../atoms";
import {useStudyBasicsQuery} from "../generated/types";
import {useExpressionValues} from "../hooks";
import Plot from 'react-plotly.js';
import * as aq from 'arquero';

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
            const color = study.annotationValueMap.get(annotationValueId);
            const tableForAnnotation = t.params({annotationValueId}).filter((d: any, p: any) => d.annotationValueId === p.annotationValueId);
            return {
                type: 'scattergl',
                x: tableForAnnotation.array('projectionX', Float32Array),
                y: tableForAnnotation.array('projectionY', Float32Array),
                customdata: tableForAnnotation.array('studySampleId', Int32Array),
                mode: 'markers',
                marker: {
                    size: 3,
                    opacity: 0.7,
                    color, //: t.array('annotationValueId', Int32Array),
                },
                showlegend: false,
            } as Partial<Plotly.PlotData>;
        });
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
    }, [study]);


    return (
        <Group position={'apart'}>
            <LeftSidePanel/>
            <Stack>
                {preparedPlot && <Plot data={preparedPlot.plotlyData} layout={preparedPlot.plotlyLayout}/>}
            </Stack>
            <RightSidePanel/>
        </Group>
    );
};

export default AnnotationsInUmapScatterplotTestPage;
