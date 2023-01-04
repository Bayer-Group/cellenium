import React, {useState} from 'react';
import {Annotation} from "../Annotation/Annotation";
import {Stack} from "@mantine/core";

type Props = {
    annotations: any[];
}
function AnnotationGroup({annotations}:Props) {
    const [selected, setSelected] = useState<string>(annotations[0].label)
    return (
        <Stack>
            {annotations.map((annot)=>{
                return <Annotation onSelect={setSelected} label={annot.label} color={annot.color} selected={annot.label===selected}/>
            })}
        </Stack>
    );
};

export {AnnotationGroup};