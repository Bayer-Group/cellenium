import React, {useState} from 'react';
import {Annotation} from "../Annotation/Annotation";
import {Stack} from "@mantine/core";
import {AnnotationGroup} from "../../generated/types";
import {useRecoilState} from "recoil";
import {highlightAnnotationState} from "../../atoms";

type Props = {
    annotationGroup: AnnotationGroup|undefined;
    handleSelection: Function;
    selectedAnnotation: string|undefined;
}

function AnnotationGroupDisplay({annotationGroup,handleSelection, selectedAnnotation}: Props) {
    const [highlightAnnotation, setHighlightAnnotation] = useRecoilState(highlightAnnotationState);
    return (
        <Stack spacing={'xs'}>
            {annotationGroup && annotationGroup.annotationValuesList.map((annot) => {
                return <Annotation onSelect={()=>{
                    annot.annotationValueId===highlightAnnotation?setHighlightAnnotation(0):setHighlightAnnotation(annot.annotationValueId);
                }

                } label={annot.displayValue} color={annot.color}
                                   selected={annot.annotationValueId === highlightAnnotation}/>
            })}
        </Stack>
    );
};

export {AnnotationGroupDisplay};