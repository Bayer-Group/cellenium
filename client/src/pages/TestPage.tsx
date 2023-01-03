import React, {useEffect} from 'react';
import {LeftSidePanel, RightSidePanel} from "../components";
import {Group, Stack} from "@mantine/core";
import {useRecoilState, useRecoilValue} from "recoil";
import {studyIdState, studyState} from "../atoms";
import {useStudyBasicsQuery} from "../generated/types";

const TestPage = () => {
    const [studyId, setStudyId] = useRecoilState(studyIdState);
    useEffect(() => {
        setStudyId(1)
    });
    const study = useRecoilValue(studyState);
    console.log('study', study);

    const {data} = useStudyBasicsQuery({variables: {studyId: 1}});
    console.log('data', data);


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
