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
import {annotationGroupIdState, selectedAnnotationState, selectedGenesState, studyState} from "../atoms";
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
    console.log({selectedGenes});
    const {table, loading} = useExpressionValues(selectedGenes.map(g => g.omicsId));

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
    const study = useRecoilValue(studyState);


    useEffect(() => {
        setSelectedGenes([])
    }, [])

    if (!study) {
        return <></>;
    }

    return (
        <Group style={{height: '100vh'}} align={'flex-start'} position={'apart'} spacing={'xs'}>
            <LeftSidePanel>
                <Stack>
                    <AnnotationGroupSelectBox changeHandler={(value: number) => {
                        setAnnotationGroupId(value);
                        setSelectedAnnotation(undefined);
                    }}/>
                    {annotationGroupId && study?.annotationGroupMap.get(annotationGroupId) !== undefined &&
                        <AnnotationGroupDisplay
                            annotationGroup={study.annotationGroupMap.get(annotationGroupId)}/>}
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
                        {selectedAnnotation ? <DEGTable annotationId={selectedAnnotation}/> :
                            <Text size='xs' color='dimmed'>Select cells in the plot or via the selection panel on the
                                left-hand side.</Text>}
                    </Stack>
                </Stack>
            </RightSidePanel>
        </Group>
    );
};

export default DifferentialExpressionAnalysis;