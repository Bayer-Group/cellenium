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
import ProjectionSelectBox from "../components/ProjectionSelectBox/ProjectionSelectBox";

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
    const annotationGroupId = useRecoilValue(annotationGroupIdState);
    const selectedAnnotation = useRecoilValue(selectedAnnotationState);
    const [selectedGenes, setSelectedGenes] = useRecoilState(selectedGenesState);
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
                    <ProjectionSelectBox />
                    {annotationGroupId && <AnnotationGroupSelectBox/>}
                    {annotationGroupId && study.annotationGroupMap.get(annotationGroupId)?.differentialExpressionCalculated?
                    null:<Text color={'red'} size={'xs'}>No DEGs calculated for selected group.</Text>}
                    {annotationGroupId && <AnnotationGroupDisplay/>}
                </Stack>
            </LeftSidePanel>
            <main>
                {<ProjectionPlotWithOptionalExpression/>}
            </main>
            <RightSidePanel>
                <Stack>
                    <Divider size={"xs"} label={'Gene store'}/>
                    <UserGeneStore multiple={false}/>
                    <Space h={'xs'}/>
                    <Stack>
                        <Divider size={"xs"} label={'Differential expression table'}/>
                        {annotationGroupId && study.annotationGroupMap.get(annotationGroupId)?.differentialExpressionCalculated?
                    null:<Text color={'red'} size={'xs'}>No DEGs calculated for selected group.</Text>}
                        {!selectedAnnotation &&
                            study.annotationGroupMap.get(annotationGroupId as number)?.differentialExpressionCalculated === true &&
                            <Text size='xs' color='dimmed'>Select cells in the plot or via the selection panel on the
                                left-hand side.</Text>}
                        {selectedAnnotation !== undefined && <DEGTable annotationId={selectedAnnotation}/>}
                    </Stack>
                </Stack>
            </RightSidePanel>
        </Group>
    );
};

export default DifferentialExpressionAnalysis;