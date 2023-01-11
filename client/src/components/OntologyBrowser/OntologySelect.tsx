import React, {useState} from "react";
import {SelectBoxItem} from "../../model";
import {Select} from "@mantine/core";
import {OntologyOverviewFragment} from "../../generated/types";

type OSProps = {
    handleChange: any;
    ontologies: OntologyOverviewFragment[];
}

const OntologySelect = ({ontologies, handleChange}: OSProps) => {
    const [value, setValue] = useState<string>();

    let selectOntologies: SelectBoxItem[] = ontologies.map((ele) => {
        return {
            value: ele.name,
            label: ele.name
        }
    })

    function updateChange(item: string) {
        setValue(item)
        handleChange(item)
    }

    return (
        <Select
            value={value}
            onChange={updateChange}
            label={'Select ontology'}
            placeholder={'Pick one'}
            data={selectOntologies}
        />


    )
}
export default OntologySelect;
