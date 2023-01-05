import React, {useState} from 'react';
import {Annotation} from "../Annotation/Annotation";
import {Stack} from "@mantine/core";
import {AnnotationGroup} from "../../generated/types";

type Props = {
    annotationGroup: AnnotationGroup|undefined;
}

function AnnotationGroupDisplay({annotationGroup}: Props) {
    const [selected, setSelected] = useState<string>()

    return (
        <Stack spacing={'xs'}>
            {annotationGroup && annotationGroup.annotationValuesList.map((annot) => {
                return <Annotation onSelect={setSelected} label={annot.displayValue} color={annot.color}
                                   selected={annot.displayValue === selected}/>
            })}
        </Stack>
    );
};

export {AnnotationGroupDisplay};