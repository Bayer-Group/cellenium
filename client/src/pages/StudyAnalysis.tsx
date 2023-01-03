import React from 'react';
import {LeftSidePanel, RightSidePanel} from "../components";
import {Group} from "@mantine/core";

const StudyAnalysis = () => {
    return (
        <Group position={'apart'}>
            <LeftSidePanel />

            <RightSidePanel />
        </Group>
    );
};

export default StudyAnalysis;