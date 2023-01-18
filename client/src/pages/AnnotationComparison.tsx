import React from 'react';
import {Group} from "@mantine/core";
import {LeftSidePanel, RightSidePanel} from "../components";
import {SankeyPlot} from "../components/";

const AnnotationComparison = () => {

    return (
        <Group position={'apart'}>
            <LeftSidePanel>

            </LeftSidePanel>
            <main style={{height: '100vh', overflowY: 'scroll', flexGrow: 1, paddingTop: 60}}>
                <SankeyPlot/>
            </main>

            <RightSidePanel>

            </RightSidePanel>
        </Group>

    );
};

export default AnnotationComparison;