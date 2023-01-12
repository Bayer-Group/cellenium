import React, {useMemo, useState} from 'react';
import {MultiSelect, Select, SelectItem} from '@mantine/core';
import {useRecoilState, useRecoilValue} from "recoil";
import {selectedAnnotationFilterState, studyState} from "../../atoms";
import _ from 'lodash';

const AnnotationFilterSelectBox = () => {
    const study = useRecoilValue(studyState);
    const [value, setValue] = useState<string[]|undefined>();

    const [selectedAnnotationFilter, setSelectedAnnotationFilter] = useRecoilState(selectedAnnotationFilterState);
    const annotations: SelectItem[] = useMemo(() => {
        const anns: SelectItem[] = [];
        if (study) {
            study.annotationGroupMap.forEach((value, key) => {
                const groupName = value.displayGroup;
                value.annotationValuesList.map((av) => anns.push({
                    value: av.annotationValueId.toString(),
                    label: av.displayValue,
                    group: groupName
                }))
            });
        }
        return anns;
    }, [study]);
    return (
        <MultiSelect
            value={value}
            onChange={(value) => {
                if (value) {
                    setSelectedAnnotationFilter(_.union([...selectedAnnotationFilter], value.map((v)=>parseInt(v))))
                }
            }}
            searchable
            placeholder="Select"
            data={annotations}
            onSubmit={(event) => console.log(event)}
        />
    );
};

export default AnnotationFilterSelectBox;