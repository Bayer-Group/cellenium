import {Select} from '@mantine/core';
import {useState} from "react";
import {SelectBoxItem} from "../../model";

type Props = {
    annotations: SelectBoxItem[];
    changeHandler: Function;
}

function AnnotationGroupSelectBox({annotations, changeHandler}: Props) {
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
            placeholder="Pick one"
            transitionDuration={80}
            transitionTimingFunction="ease"
            data={annotations as any}
        />
    );
}

export {AnnotationGroupSelectBox};