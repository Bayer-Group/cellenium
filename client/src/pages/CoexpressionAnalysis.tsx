import React, {useEffect, useMemo, useState} from 'react';
import ExpressionAnalysisTypeSelectBox
    from "../components/ExpressionAnalysisTypeSelectBox/ExpressionAnalysisTypeSelectBox";
import {Group, Space, Stack, Text} from "@mantine/core";
import {AddGene, AnnotationGroupSelectBox, LeftSidePanel, RightSidePanel} from "../components";
import UserGene from "../components/UserGene/UserGene";
import {useRecoilState, useRecoilValue} from "recoil";
import {selectedGenesState, studyState} from "../atoms";
import {useExpressionValues} from "../hooks";
import * as aq from 'arquero';
import ColumnTable from 'arquero/dist/types/table/column-table';
import Plot from "react-plotly.js";

interface PreparedPlot {
    message?: string;
    allSameSampleExprValues: ColumnTable[];
    combinations: { dimX: number; dimY: number }[];
    plotlyData: Partial<Plotly.PlotData>[];
    plotlyLayout: Partial<Plotly.Layout>;
}

const CoexpressionAnalysis = () => {
    const selectedGenes = useRecoilValue(selectedGenesState);
    const {table, loading} = useExpressionValues(selectedGenes.map(g => g.omicsId));
    const study = useRecoilValue(studyState);

    const preparedPlot = useMemo(() => {
        const distinctValues: number[] =
            table
                ?.rollup({omicsIds: aq.op.array_agg_distinct('omicsId')})
                .array('omicsIds')[0]
            /*.sort((a, b) => a.toLowerCase().localeCompare(b.toLowerCase()))*/ || [];

        if (!table || distinctValues.length < 2 || distinctValues.length > 8) {
            return {
                message: 'To show a co-expression plot, at least two genes (and not more than eight genes) must be selected',
                allSameSampleExprValues: [],
                combinations: [],
                plotlyData: [],
                plotlyLayout: {},
            };
        }

        const dimsX = distinctValues.slice(0, distinctValues.length - 1);
        const dimsY = distinctValues.reverse().slice(0, distinctValues.length - 1);
        const combinations = dimsX
            .map((x) => dimsY.map((y) => ({dimX: x, dimY: y})))
            .flat(2)
            .filter((d) => distinctValues.indexOf(d.dimX) > distinctValues.indexOf(d.dimY));

        const allSameSampleExprValues = combinations.map((combination, i) => {
            const subTableA = table
                .params({plotFilter: combination.dimX})
                .filter((d: any, p: any) => d.omicsId === p.plotFilter)
                .select({studySampleId: 'studySampleId', value: 'valueA'});
            const subTableB = table
                .params({plotFilter: combination.dimY})
                .filter((d: any, p: any) => d.omicsId === p.plotFilter)
                .select({studySampleId: 'studySampleId', value: 'valueB'});
            const sameSampleExprValues = subTableA.join(subTableB, ['studySampleId', 'studySampleId']).derive({index: () => aq.op.row_number() - 1});
            return sameSampleExprValues;
        });

        const markerSize = Math.max(Math.min(1000 / dimsX.length / table.numRows(), 0.5), 3.0);

        const subplots: Partial<Plotly.PlotData>[] = combinations.map((combination, i) => {
            const sameSampleExprValues = allSameSampleExprValues[i];
            return {
                type: 'scattergl',
                x: sameSampleExprValues.array('valueA', Float32Array),
                y: sameSampleExprValues.array('valueB', Float32Array),
                customdata: sameSampleExprValues.array('studySampleId', Int32Array),
                xaxis: `x${dimsX.indexOf(combination.dimX) === 0 ? '' : dimsX.indexOf(combination.dimX) + 1}`,
                yaxis: `y${dimsY.indexOf(combination.dimY) === 0 ? '' : dimsY.indexOf(combination.dimY) + 1}`,
                domain: {
                    column: dimsX.indexOf(combination.dimX),
                    row: dimsY.indexOf(combination.dimY),
                },
                mode: 'markers',
                marker: {
                    size: markerSize,
                    color: '#415370',
                    opacity: 0.6,
                },
                hoverinfo: 'none',
                showlegend: false,
            } as Partial<Plotly.PlotData>;
        });

        const axisConfig = Object();
        dimsX.forEach((g, i) => {
            axisConfig[`xaxis${i === 0 ? '' : i + 1}`] = {
                // TODO display gene symbol from omics lookup map
                title: `omics ID ${g} ...`,
                side: 'top',
                // put on top of the plot of the 'y' row (i.e. not y2 etc)
                anchor: 'y',
                color: '#737373',
                zerolinecolor: '#d5d5d5',
            };
        });
        dimsY.forEach((g, i) => {
            axisConfig[`yaxis${i === 0 ? '' : i + 1}`] = {
                title: `omics ID ${g} ...`,
                color: '#737373',
                zerolinecolor: '#d5d5d5',
            };
        });

        return {
            allSameSampleExprValues,
            plotlyData: subplots,
            combinations,
            plotlyLayout: {
                autosize: true,
                height: 1000,
                // margin: { l: 0, r: 0, t: 0, b: 0 },
                ...axisConfig,
                grid: {
                    columns: dimsX.length,
                    rows: dimsY.length,
                    // pattern: 'independent',
                    xgap: 0,
                    ygap: 0,
                },
                dragmode: false,
                clickmode: 'none',
            } as Partial<Plotly.Layout>,
        } as PreparedPlot;
    }, [table]);


    if (!study) {
        return <></>;
    }
    return (
        <Group position={'apart'}>
            <LeftSidePanel>
            </LeftSidePanel>
            <main>
                {preparedPlot?.message && <Text>{preparedPlot.message}</Text>}
                {preparedPlot && !preparedPlot.message &&
                    <Plot data={preparedPlot.plotlyData}
                          layout={preparedPlot.plotlyLayout}
                          config={{
                              responsive: true
                          }}
                    />
                }
            </main>
            <RightSidePanel>
                <Stack align={'flex-start'} justify={'flex-start'} spacing={'md'}>
                    <div style={{width: '80%'}}>
                        <AddGene/>
                    </div>
                    <Space h={'md'}/>
                    {/* TODO reusable gene panel, from CellMarkerAnalysis.tsx */}
                </Stack>
            </RightSidePanel>
        </Group>
    );
};

export default CoexpressionAnalysis;