import React, {useEffect, useState} from 'react';
import ExpressionAnalysisTypeSelectBox
    from "../components/ExpressionAnalysisTypeSelectBox/ExpressionAnalysisTypeSelectBox";
import {Group, Space, Stack} from "@mantine/core";
import {AddGene, AnnotationGroupSelectBox, LeftSidePanel, RightSidePanel} from "../components";
import UserGene from "../components/UserGene/UserGene";
import {useRecoilState, useRecoilValue} from "recoil";
import {studyIdState, studyState} from "../atoms";
import {SelectBoxItem} from "../model";


const CoexpressionAnalysis = () => {

    return (
        <Group position={'apart'}>
            <LeftSidePanel>
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