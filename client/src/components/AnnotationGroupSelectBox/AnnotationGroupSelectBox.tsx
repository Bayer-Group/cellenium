import React from 'react';
import {Select} from '@mantine/core';
import {useMemo, useState} from "react";
import {SelectBoxItem} from "../../model";
import {useRecoilValue} from "recoil";
import {studyState} from "../../atoms";

type Props = {
    changeHandler: Function;
}

function AnnotationGroupSelectBox({changeHandler}: Props) {
    const study = useRecoilValue(studyState);
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
        return anns;
    }, [study]);

    const [value, setValue] = useState<string | null>(annotations[0].value);

    function update(value: string | null) {
        if (value) {
            setValue(value)
            changeHandler(parseInt(value))
        }
    }

    return (
        <Select
            value={value} onChange={(value) => update(value)}
            label="Select annotation group"
            labelProps={{size:'xs'}}
            placeholder="Pick one"
            transitionDuration={80}
            transitionTimingFunction="ease"
            data={annotations as any}
        />
    );
}

export {AnnotationGroupSelectBox};