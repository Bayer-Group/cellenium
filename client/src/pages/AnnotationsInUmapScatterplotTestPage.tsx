import React, {useEffect, useState} from 'react';
import {LeftSidePanel, RightSidePanel} from "../components";
import {Group, Stack} from "@mantine/core";
import {useRecoilState, useRecoilValue} from "recoil";
import {studyIdState, studyState} from "../atoms";

import * as Plotly from "plotly.js";
import ProjectionPlot from "../components/ProjectionPlot/ProjectionPlot";

const AnnotationsInUmapScatterplotTestPage = () => {
    const [studyId, setStudyId] = useRecoilState(studyIdState);
    useEffect(() => {
        setStudyId(1)
    });

    return (
        <Group position={'apart'}>
            <LeftSidePanel/>
            <Stack>
                <ProjectionPlot colorBy={'annotation'}/>
            </Stack>
            <RightSidePanel/>
        </Group>
    );
};

export default AnnotationsInUmapScatterplotTestPage;
