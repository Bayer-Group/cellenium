import React, {useState} from 'react';
import {Select, Stack, useMantineTheme, Text} from "@mantine/core";
import {useOntologiesQuery, useStudiesQuery,} from "../../generated/types";
import {OntologyTree} from "./OntologyTree";
import OntologySelect from "./OntologySelect";
import {OntologyItem} from "../../model";


const OntologyBrowser = ({handleAddOntologyItem, ontologyTrees}:{handleAddOntologyItem: Function, ontologyTrees: any}) => {
    const {data:ontologyData, error:ontologyError, loading:ontologyLoading} = useOntologiesQuery()
    const [selectedOntology, setSelectedOntology] = useState<string>();
    const theme = useMantineTheme()
    return (
        <Stack>
            {ontologyData && <OntologySelect handleChange={setSelectedOntology} ontologies={ontologyData.ontologiesList}/>}
            {ontologyTrees && ontologyTrees.get(selectedOntology)!==undefined && selectedOntology &&
                <OntologyTree ontology={ontologyTrees.get(selectedOntology) as any} handleAddOntologyItem={handleAddOntologyItem}/>}
        </Stack>
    );
};

export {OntologyBrowser};