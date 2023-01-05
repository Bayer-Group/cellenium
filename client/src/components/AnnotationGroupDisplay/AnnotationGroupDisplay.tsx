import React, {useState} from 'react';
import {Annotation} from "../Annotation/Annotation";
import {Stack} from "@mantine/core";
import {AnnotationGroup} from "../../generated/types";

type Props = {
    annotationGroup: AnnotationGroup|undefined;
    handleSelection: Function;
    selectedAnnotation: string|undefined;
}

function AnnotationGroupDisplay({annotationGroup,handleSelection, selectedAnnotation}: Props) {

    return (
        <Stack spacing={'xs'}>
            {annotationGroup && annotationGroup.annotationValuesList.map((annot) => {
                return <Annotation onSelect={()=>{
                    annot.displayValue===selectedAnnotation?handleSelection():handleSelection(annot.displayValue);
                }

                } label={annot.displayValue} color={annot.color}
                                   selected={annot.displayValue === selectedAnnotation}/>
            })}
        </Stack>
    );
};

export {AnnotationGroupDisplay};