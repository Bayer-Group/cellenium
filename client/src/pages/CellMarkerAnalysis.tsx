import React, {useEffect, useState} from 'react';
import {
    AddGene,
    AnnotationGroupDisplay,
    AnnotationGroupSelectBox,
    DEGTable,
    LeftSidePanel,
    RightSidePanel
} from "../components";
import {Group, Space, Stack, Text} from "@mantine/core";
import UserGene from "../components/UserGene/UserGene";
import {useRecoilState, useRecoilValue} from "recoil";
import {annotationGroupIdState, selectedAnnotationState, studyIdState, studyState, userGenesState} from "../atoms";
import {SelectBoxItem} from '../model';
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

    //const [selectedAnnotationGroup, setSelectedAnnotationGroup] = useState<number>();
    const [selectedAnnotation, setSelectedAnnotation] = useRecoilState(selectedAnnotationState);
    const userGenes = useRecoilValue(userGenesState);
    const [annotationGroups, setAnnotationGroups] = useState<SelectBoxItem[]>([])
    const study = useRecoilValue(studyState);
    useEffect(() => {
        if (study) {
            const anns: SelectBoxItem[] = []
            study.annotationGroupMap.forEach((value, key) => {
                anns.push({
                    value: key.toString(),
                    label: value.displayGroup
                })
            })
            if (anns.length > 0) {
                setAnnotationGroups(anns)
                setAnnotationGroupId(parseInt(anns[0].value))
            }

        }

    }, [study])

    return (
        <Group position={'apart'} spacing={'xs'}>
            <LeftSidePanel>

                {annotationGroups.length > 0 &&
                    <Stack>
                        <AnnotationGroupSelectBox annotations={annotationGroups}
                                                  changeHandler={(value: number) => {
                                                      setAnnotationGroupId(value);
                                                      setSelectedAnnotation(undefined);
                                                  }}/>
                        {annotationGroupId && study?.annotationGroupMap.get(annotationGroupId) !== undefined &&
                            <AnnotationGroupDisplay
                                annotationGroup={study.annotationGroupMap.get(annotationGroupId)}/>}
                    </Stack>
                }

            </LeftSidePanel>
            <main>
                <ProjectionPlot colorBy={'annotation'}/>
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
                    <Space h={'md'}/>
                    <Stack>
                        <Text weight={800} size={'xs'}>Differentially expressed genes</Text>
                        {selectedAnnotation ? <DEGTable annotationId={selectedAnnotation}/> :
                            <Text size='xs' color='gray'>Nothing selected yet.</Text>}

                    </Stack>
                </Stack>
            </RightSidePanel>
        </Group>
    );
};

export default DifferentialExpressionAnalysis;