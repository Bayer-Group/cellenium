import {Select, Stack} from '@mantine/core';
import {SelectBoxItem} from "../../model";

interface Props {
    annotationGroups: SelectBoxItem[];
    handleChange1: Function;
    value1: string|undefined;
    handleChange2: Function;
    value2: string|undefined;
}

function SankeyAnnotationGroupSelector({annotationGroups, value1, value2, handleChange1, handleChange2}: Props) {

    return (
        <Stack spacing={10}>
            <Select
                value={value1} onChange={(value) => handleChange1(value)}
                label="Select annotation group 1"
                labelProps={{size: 'xs'}}
                placeholder="Pick one"
                transitionDuration={80}
                transitionTimingFunction="ease"
                data={annotationGroups}
            />
            <Select
                value={value2} onChange={(value) => handleChange2(value)}
                label="Select annotation group 2"
                labelProps={{size: 'xs'}}
                placeholder="Pick one"
                transitionDuration={80}
                transitionTimingFunction="ease"
                data={annotationGroups.filter((ele)=>ele.value!==value1)}
            />
        </Stack>
    );
}

export {SankeyAnnotationGroupSelector};