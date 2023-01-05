import React, {useEffect, useState} from 'react';
import {
    AddGene,
    AnnotationGroupDisplay,
    AnnotationGroupSelectBox,
    DEGTable,
    LeftSidePanel,
    RightSidePanel
} from "../components";
import {Group, Space, Stack} from "@mantine/core";
import UserGene from "../components/UserGene/UserGene";
import {useRecoilState, useRecoilValue} from "recoil";
import {studyIdState, studyState} from "../atoms";
import {SelectBoxItem} from '../model';

const ANNOTATIONS = [
    {label: "bone cell", color: "#1f77b4"},
    {label: "chondrocyte", color: "#ff7f0e"},
    {label: "endothelial cell", color: "#2ca02c"},
    {label: "endothelial cell of artery", color: "#d62728"},
    {label: "fibroblast", color: "#9467bd"},
    {label: "mesenchymal stem cell", color: "#8c564b"},
    {label: "pericyte cell", color: "#e377c2"}
]

const DEG: object[] = [
    {
        symbol: 'BRD4',
        padj: 0.001342,
        log2fc: 2132.23
    },
    {
        symbol: 'PTK2',
        padj: 0.001342,
        log2fc: 2.23
    },
    {
        symbol: 'CDK2',
        padj: 0.001342,
        log2fc: 24
    },
    {
        symbol: 'EGFR',
        padj: 0.001342,
        log2fc: 1
    },
    {
        symbol: 'KRAS',
        padj: 0.001342,
        log2fc: 3.324
    },
    {
        symbol: 'KLK3',
        padj: 0.00000000000000001342,
        log2fc: 432.3231322
    },
    {
        symbol: 'PLK1',
        padj: 0.001342,
        log2fc: 2132.23
    },
]

const genes = [
    {display_symbol: 'CDK2'},
    {display_symbol: 'KRAS'},
    {display_symbol: 'BRD4'},
    {display_symbol: 'KLK3'},
    {display_symbol: 'ATAD2'},

]

interface PreparedPlot {
    message?: string;
    plotlyData: Partial<Plotly.PlotData>[];
    plotlyLayout: Partial<Plotly.Layout>;
}

function DifferentialExpressionAnalysis() {
    const [selectedAnnotationGroup, setSelectedAnnotationGroup] = useState<number>();
    const [selectedAnnotation, setSelectedAnnotation] = useState<string|undefined>();
    const [studyId, setStudyId] = useRecoilState(studyIdState);
    const [annotationGroups, setAnnotationGroups] = useState<SelectBoxItem[]>([])
    useEffect(() => {
        setStudyId(1)
    });
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
                setSelectedAnnotationGroup(parseInt(anns[0].value))
            }

        }

    }, [study])
    useEffect(() => {
        if (selectedAnnotationGroup && study) {
            console.log(selectedAnnotationGroup);
            console.log(study.annotationGroupMap.get(selectedAnnotationGroup))
        }
    }, [selectedAnnotationGroup])
    return (
        <Group position={'apart'}>
            <LeftSidePanel>

                {annotationGroups.length > 0 &&
                    <Stack>
                        <AnnotationGroupSelectBox annotations={annotationGroups}
                                                  changeHandler={(value:number)=>{
                                                      setSelectedAnnotationGroup(value);
                                                      setSelectedAnnotation(undefined);
                                                  }}/>
                        {selectedAnnotationGroup && study?.annotationGroupMap.get(selectedAnnotationGroup) !== undefined &&
                            <AnnotationGroupDisplay
                                handleSelection={setSelectedAnnotation}
                                selectedAnnotation={selectedAnnotation}
                                annotationGroup={study.annotationGroupMap.get(selectedAnnotationGroup)}/>}
                    </Stack>
                }

            </LeftSidePanel>

            <RightSidePanel>
                <Stack align={'flex-start'} justify={'flex-start'} spacing={'md'}>
                    <div style={{width: '80%'}}>
                        <AddGene/>
                    </div>
                    <Stack style={{width: '100%'}} spacing={'xs'}>
                        {genes.map((gene) => <UserGene key={`ug_${gene.display_symbol}`} gene={gene}/>)}
                    </Stack>
                    <Space h={'md'}/>
                    <DEGTable data={DEG}/>
                </Stack>
            </RightSidePanel>
        </Group>
    );
};

export default DifferentialExpressionAnalysis;