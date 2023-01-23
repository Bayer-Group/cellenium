import React, {useEffect, useMemo, useState} from 'react';
import {Group, Space, Stack, Text} from "@mantine/core";
import {LeftSidePanel, RightSidePanel} from "../components";
import {useRecoilState, useRecoilValue} from "recoil";
import {celltypeDiscoveryGenesState, selectedGenesState, studyState, userGenesState} from "../atoms";
import {useExpressionValues} from "../hooks";
import * as aq from 'arquero';
import Plot from "react-plotly.js";
import {SingleGeneSelection} from "../components/AddGene/SingleGeneSelection";
import {Omics} from "../model";
import ProjectionPlot from "../components/ProjectionPlot/ProjectionPlot";
import * as Plotly from "plotly.js";

interface PreparedPlot {
    message?: string;
    plotlyData: Partial<Plotly.PlotData>[];
    plotlyLayout: Partial<Plotly.Layout>;
}

const plotlyConfig: Partial<Plotly.Config> = {
    responsive: true
};

const CelltypeDiscovery = () => {
    const [omicsAll, setOmicsAll] = useRecoilState(celltypeDiscoveryGenesState);

    // convenient default for gene input
    const userGenes = useRecoilValue(userGenesState);
    useEffect(() => {
        if (!omicsAll[0] && !omicsAll[1]) {
            if (userGenes.length > 1) {
                setOmicsAll([userGenes[0], userGenes[1]]);
            } else if (userGenes.length > 0) {
                setOmicsAll([userGenes[0], null]);
            }
        }
    }, [omicsAll, userGenes]);

    const omicsX = omicsAll[0];
    const omicsY = omicsAll[1];
    const setOmicsX = (o: Omics | null) => setOmicsAll(old => [o, old[1]]);
    const setOmicsY = (o: Omics | null) => setOmicsAll(old => [old[0], o]);
    const {table, loading} = useExpressionValues((omicsX && omicsY) ? [omicsX.omicsId, omicsY.omicsId] : [], false);
    const study = useRecoilValue(studyState);

    const preparedPlot = useMemo(() => {
        if (!table || !omicsX || !omicsY) {
            return undefined;
        }

        const subTableA = table
            .params({plotFilter: omicsX.omicsId})
            .filter((d: any, p: any) => d.omicsId === p.plotFilter)
            .select({studySampleId: 'studySampleId', value: 'valueA'});
        const subTableB = table
            .params({plotFilter: omicsY.omicsId})
            .filter((d: any, p: any) => d.omicsId === p.plotFilter)
            .select({studySampleId: 'studySampleId', value: 'valueB'});
        const sameSampleExprValues = subTableA.join(subTableB, ['studySampleId', 'studySampleId']).derive({index: () => aq.op.row_number() - 1});

        const plot: Partial<Plotly.PlotData> = {
            type: 'scattergl',
            x: sameSampleExprValues.array('valueA', Float32Array),
            y: sameSampleExprValues.array('valueB', Float32Array),
            customdata: sameSampleExprValues.array('studySampleId', Int32Array),
            mode: 'markers',
            marker: {
                //size: markerSize,
                color: '#415370',
                opacity: 0.6,
            },
            hoverinfo: 'none',
            showlegend: false,
        } as Partial<Plotly.PlotData>;

        return {
            //allSameSampleExprValues,
            plotlyData: [plot],
            //combinations,
            plotlyLayout: {
                //autosize: true,
                width: 250,
                height: 250,
                margin: {l: 40, r: 0, t: 0, b: 40},
                dragmode: 'lasso',
                clickmode: 'none',
            } as Partial<Plotly.Layout>,
        } as PreparedPlot;
    }, [table]);

    const [selectedSampleIds, setSelectedSampleIds] = useState<number[] | undefined>(undefined);

    const onCoexpressionSelection = (event: Readonly<Plotly.PlotSelectionEvent>) => {
        const theelectedSampleIds = event.points.map(p => p.customdata) as number[];
        setSelectedSampleIds(theelectedSampleIds);
    };


    if (!study) {
        return <></>;
    }
    return (
        <Group position={'apart'}>
            <LeftSidePanel>
            </LeftSidePanel>
            <main>
                <ProjectionPlot colorBy={'annotation'} showSampleIds={selectedSampleIds}/>
            </main>
            <RightSidePanel>
                <Stack align={'flex-start'} justify={'flex-start'} spacing={'md'}>
                    <div style={{width: '80%'}}>
                        <Group>
                            <SingleGeneSelection selection={omicsX} onSelectionChange={setOmicsX}/>
                            <SingleGeneSelection selection={omicsY} onSelectionChange={setOmicsY}/>
                        </Group>
                    </div>
                    {preparedPlot && !preparedPlot.message &&
                        <Plot data={preparedPlot.plotlyData}
                              layout={preparedPlot.plotlyLayout}
                              config={plotlyConfig}
                              onSelected={onCoexpressionSelection}
                        />
                    }
                </Stack>
            </RightSidePanel>
        </Group>
    );
};

export default CelltypeDiscovery;
