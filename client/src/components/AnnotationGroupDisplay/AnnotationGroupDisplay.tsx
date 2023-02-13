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
    const annotations = study.annotationGroupMap.get(annotationGroupId)?.annotationValuesList;
    const isSelectable = study.annotationGroupMap.get(annotationGroupId)?.differentialExpressionCalculated;
    return (
        <Stack spacing={2} onMouseLeave={() => setHighlightAnnotation(0)}
            style={{maxWidth: 205}}
        >
            {annotations !== undefined && annotations.map((annot) => {
                return <Annotation key={annot.annotationValueId} label={annot.displayValue}
                                   sampleCount={annot.sampleCount} color={annot.color}
                                   annotationId={annot.annotationValueId}
                isSelectable={isSelectable as boolean}/>
            })}
        </Stack>
    );
};

export {AnnotationGroupDisplay};