import React, {useEffect} from 'react';
import {LeftSidePanel, RightSidePanel} from "../components";
import {Group, Stack} from "@mantine/core";
import {useRecoilState, useRecoilValue} from "recoil";
import {studyIdState, studyState} from "../atoms";
import {useStudyBasicsQuery} from "../generated/types";
import {useExpressionValues} from "../hooks";

const TestPage = () => {
    const [studyId, setStudyId] = useRecoilState(studyIdState);
    useEffect(() => {
        setStudyId(1)
    });
    const study = useRecoilValue(studyState);
    console.log('study', study);
    study.samplesProjectionTable.print();

    const {table, loading} = useExpressionValues();
    if(table) {
        table.print();
    }


    return (
        <Group position={'apart'}>
            <LeftSidePanel/>
            <Stack>

            </Stack>

            <RightSidePanel/>
        </Group>
    );
};

export default TestPage;
