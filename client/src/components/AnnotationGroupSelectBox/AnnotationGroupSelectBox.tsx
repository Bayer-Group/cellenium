import React, {useMemo, useState} from 'react';
import {Select} from '@mantine/core';
import {SelectBoxItem} from "../../model";
import {useRecoilState, useRecoilValue} from "recoil";
import {annotationGroupIdState, selectedAnnotationState, studyState} from "../../atoms";


function AnnotationGroupSelectBox() {
    const study = useRecoilValue(studyState);
    const [annotationGroupId, setAnnotationGroupId] = useRecoilState(annotationGroupIdState);
    const [selectedAnnotation, setSelectedAnnotation] = useRecoilState(selectedAnnotationState);
    const [value, setValue] = useState<string | undefined>();
    const annotations: SelectBoxItem[] = useMemo(() => {
        const anns: SelectBoxItem[] = [];
        if (study) {
            study.annotationGroupMap.forEach((value, key) => {
                anns.push({
                    value: key.toString(),
                    label: value.displayGroup
                })
            });
        }
        setValue(annotationGroupId?.toString())
        return anns;
    }, [study, annotationGroupId]);

    function update(value: string | null) {
        if (value) {
            setValue(value)
            setAnnotationGroupId(parseInt(value))
            setSelectedAnnotation(undefined)
        }
    }

    return (
        <Select
            value={value} onChange={(value) => update(value)}
            label="Select annotation group"
            labelProps={{size: 'xs'}}
            placeholder="Pick one"
            transitionDuration={80}
            transitionTimingFunction="ease"
            data={annotations as any}
        />
    );
}

export {AnnotationGroupSelectBox};