import React, {useState} from 'react';
import {Select, Stack, useMantineTheme, Text} from "@mantine/core";
import {OntologyOverviewFragment, useOntologiesQuery,} from "../../generated/types";
import {OntologyTree} from "./OntologyTree";
import {SelectBoxItem} from "../../model";

type OSProps = {
    handleChange: any;
    ontologies: OntologyOverviewFragment[];
}
const OntologySelect = ({ontologies, handleChange}: OSProps) => {
    const [value, setValue] = useState<string>();

    let selectOntologies: SelectBoxItem[] = ontologies.map((ele) => {
        return {
            value: ele.ontid.toString(),
            label: ele.name
        }
    })
    function updateChange(item: string){
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
const OntologyBrowser = () => {
    const {data, error, loading} = useOntologiesQuery()
    const [ontology, setOntology] = useState<string | undefined>();
    const theme = useMantineTheme()
    return (
        <Stack>
            {data && data.ontologiesList &&
                <OntologySelect handleChange={setOntology} ontologies={data.ontologiesList}/>}
            {ontology ? <OntologyTree ontid={parseInt(ontology)}/> :
                <Text color={theme.colors.gray[3]}>Please select an ontology for browsing!</Text>}
        </Stack>
    );
};

export {OntologyBrowser};