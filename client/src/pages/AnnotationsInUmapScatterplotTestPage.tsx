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

        const plotlyData = [
            {
                type: 'scattergl',
                x: t.array('projectionX', Float32Array),
                y: t.array('projectionY', Float32Array),
                customdata: t.array('studySampleId', Int32Array),
                mode: 'markers',
                marker: {
                    size: 3,
                    opacity: 0.7,
                    color: t.array('annotationValueId', Int32Array),
                    cmin: 1.0,
                    cmax: 15.0,

                    // TODO we need to provide the color strings in color:, plotly doesn't support a categorical palette yet.
                    colorscale: [[0.0, 'rgb(0,0,255)'], [1.0, 'rgb(255,0,255)']],
                },
            } as Partial<Plotly.PlotData>
        ];
        return {
            plotlyData,
            plotlyLayout: {
                width: 850,
                height: 700,
                margin: {l: 0, r: 0, t: 0, b: 0},
                legend: {
                    // Make the legend marker size larger
                    // @ts-ignore
                    itemsizing: 'constant',
                },
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
