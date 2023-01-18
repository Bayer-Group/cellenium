import React from 'react';
import {Annotation} from "../Annotation/Annotation";
import {Stack} from "@mantine/core";
import {useRecoilState, useRecoilValue} from "recoil";
import {annotationGroupIdState, highlightAnnotationState, studyState} from "../../atoms";


function AnnotationGroupDisplay() {
    const [highlightAnnotation, setHighlightAnnotation] = useRecoilState(highlightAnnotationState);
    const annotationGroupId = useRecoilValue(annotationGroupIdState);
    const study = useRecoilValue(studyState);

    if (!study || !annotationGroupId) {
        return <></>;
    }
    let annotations = study.annotationGroupMap.get(annotationGroupId)?.annotationValuesList;
    console.log("HH", annotationGroupId)
    return (
        <Stack spacing={2} onMouseLeave={() => setHighlightAnnotation(0)}>
            {annotations !== undefined && annotations.map((annot) => {
                return <Annotation key={annot.annotationValueId} label={annot.displayValue}
                                   sampleCount={annot.sampleCount} color={annot.color}
                                   annotationId={annot.annotationValueId}/>
            })}
        </Stack>
    );
};

export {AnnotationGroupDisplay};