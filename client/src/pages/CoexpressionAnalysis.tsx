import React, {useEffect} from 'react';
import {Center, Divider, Group, Loader, Stack, Text, useMantineTheme} from "@mantine/core";
import {AnnotationFilterDisplay, LeftSidePanel, RightSidePanel, UserGeneStore} from "../components";
import {useRecoilState, useRecoilValue} from "recoil";
import {
    selectedAnnotationFilterState,
    selectedGenesState,
    studyIdState,
    studyLayerIdState,
    studyState,
    userGenesState,
    userGeneStoreOpenState
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
    const theme = useMantineTheme()
    const selectedGenes = useRecoilValue(selectedGenesState);
    const studyId = useRecoilValue(studyIdState);
    const studyLayerId = useRecoilValue(studyLayerIdState);
    const annotationFilter = useRecoilValue(selectedAnnotationFilterState);
    const [storeOpened, setOpened] = useRecoilState(userGeneStoreOpenState);
    useEffect(() => {
        setOpened(true)
    }, [])

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
        return <Center style={{height: '100%', width: '100%'}}><Loader variant={'dots'} color={theme.colors.gray[5]}
                                                                       size={'xl'}/></Center>;

    }
    if (!data?.correlationTrianglePlot) {
        return <Center style={{height: '100%', width: '100%'}}><Text color={'dimmed'} size={'md'}>Please select at least
            2 genes
            from the gene
            store.</Text></Center>;
    }
    return <Center style={{height: '100%', width: '100%'}}><img
        style={{width: '100%', height: selectedGenes.length > 3 ? '100%' : '', objectFit: 'fill', overflow: 'hidden'}}
        src={data.correlationTrianglePlot}/></Center>;
};


function CoexpressionAnalysis() {
    const userGenes = useRecoilValue(userGenesState);
    const study = useRecoilValue(studyState);
    const [storeOpened, setOpened] = useRecoilState(userGeneStoreOpenState);

    useEffect(() => {
        setOpened(true)
    }, [])
    if (!study) {
        return <Center style={{height: '100%', width: '100%'}}><Loader variant={'dots'} color={'gray'}/> </Center>;
    }

    return (
        <Group style={{height: '100vh'}} align={'flex-start'} position={'apart'} spacing={'xs'} noWrap={true}>
            <LeftSidePanel>
                <AnnotationFilterDisplay/>
            </LeftSidePanel>
            <main style={{height: '100vh'}}>
                <CoexpressionAnalysisPlot/>
            </main>
            <RightSidePanel>
                <Stack>
                    <Divider size={"xs"} label={'Gene store'}/>
                    <UserGeneStore multiple={true}/>
                </Stack>
            </RightSidePanel>
        </Group>
    );
};

export default CoexpressionAnalysis;