import React, {useMemo} from 'react';
import {SelectItem, Stack, Text} from "@mantine/core";
import AnnotationFilterSelectBox from "./AnnotationFilterSelectBox";
import {SelectBoxItem} from "../../model";
import {useRecoilValue} from "recoil";
import {studyState} from "../../atoms";
import AnnotationFilter from "./AnnotationFilter";

const AnnotationFilterDisplay = () => {

    return (
        <Stack spacing={0}>
            <Text size={'xs'}>Filter out cells with attribute(s)</Text>
            <AnnotationFilterSelectBox/>
        </Stack>
    );
};

export  {AnnotationFilterDisplay};
