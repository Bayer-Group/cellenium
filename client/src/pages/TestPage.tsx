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

const TestPage = () => {
    const [studyId, setStudyId] = useRecoilState(studyIdState);
    useEffect(() => {
        setStudyId(1)
    });
    const study = useRecoilValue(studyState);

    const {table, loading} = useExpressionValues();

    const preparedPlot: PreparedPlot | undefined = React.useMemo(() => {
        if (!table || !study) {
            return;
        }

        // const joinedTable = table.join(study.samplesProjectionTable, 'studySampleId').reify();
        const joinedTable = table.join(study.samplesProjectionTable, 'studySampleId').reify();
        const plotlyData = [
            {
                type: 'scattergl',
                x: joinedTable.array('projectionX', Float32Array),
                y: joinedTable.array('projectionY', Float32Array),
                customdata: joinedTable.array('studySampleId', Int32Array),
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
                selected: {
                    marker: {
                        color: '#ff0000',
                    },
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
    }, [table, study]);


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

export default TestPage;
