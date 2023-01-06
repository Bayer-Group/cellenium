import React, {useState} from 'react';
import {Annotation} from "../Annotation/Annotation";
import {Stack} from "@mantine/core";
import {AnnotationGroup} from "../../generated/types";
import {useRecoilState} from "recoil";
import {highlightAnnotationState} from "../../atoms";

type Props = {
    annotationGroup: AnnotationGroup|undefined;
}

function AnnotationGroupDisplay({annotationGroup}: Props) {
    const [highlightAnnotation, setHighlightAnnotation] = useRecoilState(highlightAnnotationState);
    return (
        <Stack spacing={2}>
            {annotationGroup && annotationGroup.annotationValuesList.map((annot) => {
                return <Annotation label={annot.displayValue} color={annot.color} annotationId={annot.annotationValueId}/>
            })}
        </Stack>
    );
};

export {AnnotationGroupDisplay};