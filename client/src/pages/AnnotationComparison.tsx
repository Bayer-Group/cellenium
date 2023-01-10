import React, {useEffect, useState} from 'react';
import ExpressionAnalysisTypeSelectBox
    from "../components/ExpressionAnalysisTypeSelectBox/ExpressionAnalysisTypeSelectBox";
import {Group, Space, Stack} from "@mantine/core";
import {AddGene, AnnotationGroupSelectBox, LeftSidePanel, RightSidePanel} from "../components";
import UserGene from "../components/UserGene/UserGene";
import {useRecoilState, useRecoilValue} from "recoil";
import {studyIdState, studyState} from "../atoms";
import {SelectBoxItem} from "../model";

const AnnotationComparison = () => {

    return (
        <Group position={'apart'}>
            <LeftSidePanel>

            </LeftSidePanel>

            <RightSidePanel>

            </RightSidePanel>
        </Group>

    );
};

export default AnnotationComparison;