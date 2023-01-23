import React, {useEffect} from 'react';
import {
    AnnotationGroupDisplay,
    AnnotationGroupSelectBox,
    DEGTable,
    LeftSidePanel,
    RightSidePanel,
    UserGeneStore
} from "../components";
import {Divider, Group, Loader, Space, Stack, Text, useMantineTheme} from "@mantine/core";
import {useRecoilState, useRecoilValue} from "recoil";
import {
    annotationGroupIdState,
    selectedAnnotationState,
    selectedGenesState,
    studyState,
    userGenesState
} from "../atoms";
import ProjectionPlot from "../components/ProjectionPlot/ProjectionPlot";
import {useExpressionValues} from "../hooks";

const ANNOTATIONS = [
    {label: "bone cell", color: "#1f77b4"},
    {label: "chondrocyte", color: "#ff7f0e"},
    {label: "endothelial cell", color: "#2ca02c"},
    {label: "endothelial cell of artery", color: "#d62728"},
    {label: "fibroblast", color: "#9467bd"},
    {label: "mesenchymal stem cell", color: "#8c564b"},
    {label: "pericyte cell", color: "#e377c2"}
]


interface PreparedPlot {
    message?: string;
    plotlyData: Partial<Plotly.PlotData>[];
    plotlyLayout: Partial<Plotly.Layout>;
}

function ProjectionPlotWithOptionalExpression() {
    const theme = useMantineTheme();
    const selectedGenes = useRecoilValue(selectedGenesState);
    const {table, loading} = useExpressionValues(selectedGenes.map(g => g.omicsId), true);

    if (loading) {
        return <div><Loader variant={'dots'} color={theme.colors.gray[5]} size={'xl'}/></div>;
    }

    return (<div>
        {table && <ProjectionPlot colorBy={'annotation'} expressionTable={table}/>}
        {table === undefined && <ProjectionPlot colorBy={'annotation'}/>}
    </div>)
}

function DifferentialExpressionAnalysis() {
    const [annotationGroupId, setAnnotationGroupId] = useRecoilState(annotationGroupIdState);
    const [selectedAnnotation, setSelectedAnnotation] = useRecoilState(selectedAnnotationState);
    const [selectedGenes, setSelectedGenes] = useRecoilState(selectedGenesState);
    const [userGenes, setUserGenes] = useRecoilState(userGenesState);
    const study = useRecoilValue(studyState);
    useEffect(() => {
        if (selectedGenes.length > 1)
            setSelectedGenes(selectedGenes.slice(0, 1))
    }, [])

    if (!study) {
        return <></>;
    }
    return (
        <Group style={{height: '100vh'}} align={'flex-start'} position={'apart'} spacing={'xs'} noWrap>
            <LeftSidePanel>
                <Stack>
                    {annotationGroupId && <AnnotationGroupSelectBox/>}
                    {annotationGroupId && <AnnotationGroupDisplay/>}
                </Stack>
            </LeftSidePanel>
            <main>
                {<ProjectionPlotWithOptionalExpression/>}
            </main>
            <RightSidePanel>
                <Stack>
                    <Divider size={"xs"} label={'User gene store'}/>
                    <UserGeneStore multiple={false}/>
                    <Space h={'xs'}/>
                    <Stack>
                        <Divider size={"xs"} label={'Differentially expressed genes'}/>
                        {!selectedAnnotation &&
                            study.annotationGroupMap.get(annotationGroupId as number)?.differentialExpressionCalculated === true &&
                            <Text size='xs' color='dimmed'>Select cells in the plot or via the selection panel on the
                                left-hand side.</Text>}
                        {selectedAnnotation !== undefined && <DEGTable annotationId={selectedAnnotation}/>}
                        {study.annotationGroupMap.get(annotationGroupId as number)?.differentialExpressionCalculated === false &&
                            <Text size={'xs'} color={'dimmed'}>
                                No differential gene expression calculated for selected annotation group.
                            </Text>}
                    </Stack>
                </Stack>
            </RightSidePanel>
        </Group>
    );
};

export default DifferentialExpressionAnalysis;