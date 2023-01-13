import React from 'react';
import {
    AnnotationGroupDisplay,
    AnnotationGroupSelectBox,
    DEGTable,
    LeftSidePanel,
    RightSidePanel,
    UserGeneStore
} from "../components";
import {Divider, Group, Space, Stack, Text} from "@mantine/core";
import {useRecoilState, useRecoilValue} from "recoil";
import {annotationGroupIdState, selectedAnnotationState, studyState, userGenesState} from "../atoms";
import ProjectionPlot from "../components/ProjectionPlot/ProjectionPlot";

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

function DifferentialExpressionAnalysis() {
    const [annotationGroupId, setAnnotationGroupId] = useRecoilState(annotationGroupIdState);
    const [selectedAnnotation, setSelectedAnnotation] = useRecoilState(selectedAnnotationState);
    const userGenes = useRecoilValue(userGenesState);
    const study = useRecoilValue(studyState);

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
                <ProjectionPlot colorBy={'annotation'}/>
            </main>
            <RightSidePanel>
                <Stack>
                    <Divider size={"xs"} label={'User gene store'}/>
                    <UserGeneStore/>
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