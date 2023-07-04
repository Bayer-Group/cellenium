import {useRecoilState, useRecoilValue} from "recoil";
import {selectedAnnotationFilterState, studyState} from "../../atoms";
import {Stack, Text} from "@mantine/core";

const AnnotationFilter = () => {
    const [selectedAnnotationFilter, ] = useRecoilState(selectedAnnotationFilterState);
    const study = useRecoilValue(studyState);
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