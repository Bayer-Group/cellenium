import React, {useEffect, useState} from 'react';
import ExpressionAnalysisTypeSelectBox
    from "../components/ExpressionAnalysisTypeSelectBox/ExpressionAnalysisTypeSelectBox";
import {Group, Space, Stack} from "@mantine/core";
import {AddGene, AnnotationGroupSelectBox, LeftSidePanel, RightSidePanel} from "../components";
import UserGene from "../components/UserGene/UserGene";
import {useRecoilState, useRecoilValue} from "recoil";
import {studyIdState, studyState} from "../atoms";
import {SelectBoxItem} from "../model";

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
const CoexpressionAnalysis = () => {
    const [analysisType, setAnalysisType] = useState<string>(analysisTypes[0].value)
    const [selectedAnnotationGroup, setSelectedAnnotationGroup] = useState<number>();
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

            <RightSidePanel>
                <Stack align={'flex-start'} justify={'flex-start'} spacing={'md'}>
                    <div style={{width: '80%'}}>
                        <AddGene/>
                    </div>

                    <Space h={'md'}/>

                </Stack>
            </RightSidePanel>
        </Group>

    );
};

export default CoexpressionAnalysis;