import React, {useEffect, useState} from 'react';
import ExpressionAnalysisTypeSelectBox
    from "../components/ExpressionAnalysisTypeSelectBox/ExpressionAnalysisTypeSelectBox";
import {Group, Space, Stack} from "@mantine/core";
import {AddGene, AnnotationGroupSelectBox, LeftSidePanel, RightSidePanel} from "../components";
import UserGene from "../components/UserGene/UserGene";
import {useRecoilState, useRecoilValue} from "recoil";
import {studyIdState, studyState, userGenesState} from "../atoms";
import {SelectBoxItem} from "../model";
import ProjectionPlot from "../components/ProjectionPlot/ProjectionPlot";
import {useExpressionValues} from "../hooks";

const analysisTypes = [
    {value: 'violinplot', label: 'Violinplot'},
    {value: 'boxplot', label: 'Boxplot'},
    {value: 'projection', label: 'Projectionplot'},
    {value: 'dot', label: 'Dotplot'},
]
const genes = [
    {display_symbol: 'CDK2'},
    {display_symbol: 'KRAS'},
    {display_symbol: 'BRD4'},
    {display_symbol: 'KLK3'},
    {display_symbol: 'ATAD2'},

]
const ExpressionAnalysis = () => {
    const [analysisType, setAnalysisType] = useState<string>(analysisTypes[0].value)
    const [selectedAnnotationGroup, setSelectedAnnotationGroup] = useState<number>();
    const [studyId, setStudyId] = useRecoilState(studyIdState);
    const userGenes = useRecoilValue(userGenesState);
    const {table, loading} = useExpressionValues();
    useEffect(() => {
        console.log({table})
    }, [table])

    const [annotationGroups, setAnnotationGroups] = useState<SelectBoxItem[]>([])
    useEffect(() => {
        setStudyId(1)
    });
    const study = useRecoilValue(studyState);
    console.log({study})
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
    return (
        <Group position={'apart'}>
            <LeftSidePanel>
                <Stack>


                    <ExpressionAnalysisTypeSelectBox handleSelection={setAnalysisType} selection={analysisType}
                                                     options={analysisTypes}/>

                    {annotationGroups.length > 0 &&
                        <Stack>
                            <AnnotationGroupSelectBox annotations={annotationGroups}
                                                      changeHandler={(value: number) => {
                                                          setSelectedAnnotationGroup(value);
                                                      }}/>

                        </Stack>
                    }
                </Stack>

            </LeftSidePanel>
            <main>
                <ProjectionPlot colorBy={'expression'} expressionTable={table}/>
            </main>
            <RightSidePanel>
                <Stack align={'flex-start'} justify={'flex-start'} spacing={'md'}>
                    <AddGene/>
                    <Stack style={{width: '100%'}} spacing={'xs'}>
                        {userGenes.map((gene) => <UserGene key={`ug_${gene.displaySymbol}`} gene={gene}/>)}
                    </Stack>
                </Stack>
            </RightSidePanel>
        </Group>

    );
};

export default ExpressionAnalysis;