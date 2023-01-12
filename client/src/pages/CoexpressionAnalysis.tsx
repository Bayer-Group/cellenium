import React, {useEffect, useMemo, useState} from 'react';
import {Group, Loader, Space, Stack, Text} from "@mantine/core";
import {
    AddGene,
    LeftSidePanel,
    RightSidePanel
} from "../components";
import UserGene from "../components/UserGene/UserGene";
import {useRecoilState, useRecoilValue} from "recoil";
import {
    selectedGenesState,
    studyLayerIdState,
    studyState, userGenesState
} from "../atoms";
import ColumnTable from 'arquero/dist/types/table/column-table';
import {useExpressionCorrelationTrianglePlotQuery} from "../generated/types";

interface PreparedPlot {
    message?: string;
    allSameSampleExprValues: ColumnTable[];
    combinations: { dimX: number; dimY: number }[];
    plotlyData: Partial<Plotly.PlotData>[];
    plotlyLayout: Partial<Plotly.Layout>;
}

const CoexpressionAnalysisPlot = () => {
    const selectedGenes = useRecoilValue(selectedGenesState);
    const studyLayerId = useRecoilValue(studyLayerIdState);

    const {data, loading} = useExpressionCorrelationTrianglePlotQuery({
        variables: {
            studyLayerId,
            omicsIds: selectedGenes.map(g => g.omicsId)
        },
        skip: selectedGenes.length < 2
    });

    if (loading) {
        return <Loader size={25}/>;
    }
    if (!data?.correlationTrianglePlot) {
        return <></>;
    }
    return <img src={data.correlationTrianglePlot}/>;
};


function CoexpressionAnalysis() {
    const userGenes = useRecoilValue(userGenesState);
    const study = useRecoilValue(studyState);

    if (!study) {
        return <></>;
    }

    return (
        <Group style={{height: '100vh'}} align={'flex-start'} position={'apart'} spacing={'xs'}>
            <LeftSidePanel>
            </LeftSidePanel>
            <main>
                <CoexpressionAnalysisPlot/>
            </main>
            <RightSidePanel>
                <Stack align={'flex-start'} justify={'flex-start'} spacing={'md'}>
                    <div>
                        <AddGene/>
                    </div>
                    <Stack spacing={'xs'}>
                        {userGenes.length > 0 ? userGenes.map((gene) => <UserGene key={`ug_${gene.displaySymbol}`}
                                                                                  gene={gene}/>) :
                            <Text color={'gray'} size={'xs'}>Nothing added yet.</Text>}
                    </Stack>
                </Stack>
            </RightSidePanel>
        </Group>
    );
};

export default CoexpressionAnalysis;