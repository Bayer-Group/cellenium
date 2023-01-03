import React, {useState} from 'react';
import {Annotation} from "../Annotation/Annotation";

type Props = {
    annotations: any[];
}
function AnnotationGroup({annotations}:Props) {
    const [selected, setSelected] = useState<string>(annotations[0].label)
    return (
        <div>
            {annotations.map((annot)=>{
                return <Annotation onSelect={setSelected} label={annot.label} color={annot.color} selected={annot.label===selected}/>
            })}
        </div>
    );
};

export {AnnotationGroup};