import React, {useEffect, useMemo, useState} from 'react';
import {Divider, Group, Loader, Space, Stack, Text} from "@mantine/core";
import {
    AddGene, AnnotationFilterDisplay,
    LeftSidePanel,
    RightSidePanel, UserGeneStore
} from "../components";
import UserGene from "../components/UserGene/UserGene";
import {useRecoilState, useRecoilValue} from "recoil";
import {
    annotationGroupIdState,
    selectedAnnotationFilterState,
    selectedGenesState, studyIdState,
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
    const studyId = useRecoilValue(studyIdState);
    const studyLayerId = useRecoilValue(studyLayerIdState);
    const annotationFilter = useRecoilValue(selectedAnnotationFilterState);

    const {data, loading} = useExpressionCorrelationTrianglePlotQuery({
        variables: {
            studyId,
            studyLayerId,
            omicsIds: selectedGenes.map(g => g.omicsId),
            excludeAnnotationValueIds: annotationFilter
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
                <AnnotationFilterDisplay/>
            </LeftSidePanel>
            <main>
                <CoexpressionAnalysisPlot/>
            </main>
            <RightSidePanel>
                <Stack>
                    <Divider size={"xs"} label={'User gene store'}/>
                    <UserGeneStore multiple={true} opened={true}/>
                </Stack>
            </RightSidePanel>
        </Group>
    );
};

export default CoexpressionAnalysis;