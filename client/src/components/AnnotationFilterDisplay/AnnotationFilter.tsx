import React from 'react';
import {useRecoilState, useRecoilValue} from "recoil";
import {selectedAnnotationFilterState, selectedAnnotationState, studyState} from "../../atoms";
import {Stack, Text} from "@mantine/core";

const AnnotationFilter = () => {
    const [selectedAnnotationFilter, setSelectedAnnotationFilter] = useRecoilState(selectedAnnotationFilterState);
    const study = useRecoilValue(studyState);
    console.log(study?.annotationValueMap)
    return (
        <Stack spacing={1}>
            {study && study.annotationValueMap && selectedAnnotationFilter.map((av:number)=>{
                const annotationValue = study.annotationValueMap.get(av)?.displayValue;
                return (<Text size={'xs'}>{annotationValue}</Text>)
            })}
        </Stack>
    );
};

export default AnnotationFilter;