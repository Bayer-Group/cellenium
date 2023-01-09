import React, {useEffect} from 'react';
import {LeftSidePanel, RightSidePanel} from "../components";
import {Group, Stack} from "@mantine/core";
import {useRecoilState, useRecoilValue} from "recoil";
import {studyIdState, studyState} from "../atoms";
import {useExpressionValues} from "../hooks";
import ProjectionPlot from "../components/ProjectionPlot/ProjectionPlot";

const SingleGeneExpressionInUmapScatterplotTestPage = () => {
    const [studyId, setStudyId] = useRecoilState(studyIdState);
    console.log('studyId', studyId);
    const {table, loading} = useExpressionValues();

    return (
        <Group position={'apart'}>
            <LeftSidePanel/>
            <Stack>
                <ProjectionPlot colorBy={'expression'} expressionTable={table}/>
            </Stack>
            <RightSidePanel/>
        </Group>
    );
};

export default SingleGeneExpressionInUmapScatterplotTestPage;
