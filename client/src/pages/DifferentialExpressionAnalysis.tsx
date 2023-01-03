import React from 'react';
import {AnnotationGroup, AnnotationGroupSelectBox, LeftSidePanel, RightSidePanel} from "../components";
import {Group} from "@mantine/core";

const ANNOTATIONS = [
    {label: "bone cell", color:"#1f77b4"},
    {label: "chondrocyte", color:"#ff7f0e"},
    {label: "endothelial cell", color:"#2ca02c"},
    {label: "endothelial cell of artery", color:"#d62728"},
    {label: "fibroblast", color:"#9467bd"},
    {label: "mesenchymal stem cell", color:"#8c564b"},
    {label: "pericyte cell", color:"#e377c2"}
]
function DifferentialExpressionAnalysis() {
    return (
        <Group position={'apart'}>
            <LeftSidePanel>
                <AnnotationGroupSelectBox/>
                <AnnotationGroup annotations={ANNOTATIONS} />
            </LeftSidePanel>

            <RightSidePanel/>
        </Group>
    );
};

export default DifferentialExpressionAnalysis;