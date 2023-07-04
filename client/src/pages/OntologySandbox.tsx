import {useEffect, useState} from 'react';
import {OntologyTree} from "../components";
import {useOntologiesQuery, useStudiesQuery} from "../generated/types";
import {OntologyItem} from "../model";
import {Stack} from "@mantine/core";
import OntologySelect from "../components/OntologyBrowser/OntologySelect";
import {generateOntologyTrees} from "./helper";

const OntologySandbox = () => {
    const {data: ontologyData} = useOntologiesQuery()
    const {data} = useStudiesQuery()
    const [ontologyTrees, setOntologyTrees] = useState<Map<string, OntologyItem>>();
    const [selectedOntology, setSelectedOntology] = useState<string>('NCIT');
    useEffect(() => {
        if (data && data.treeOntologiesList) {
            setOntologyTrees(generateOntologyTrees(data.treeOntologiesList));
        }
    }, [data])

    return (
        <Stack>
            {ontologyData &&
                <OntologySelect handleChange={setSelectedOntology} ontologies={ontologyData.ontologiesList}/>}
            {ontologyTrees && ontologyTrees.get(selectedOntology) !== undefined && selectedOntology &&
                <OntologyTree ontology={ontologyTrees.get(selectedOntology) as any} handleAddOntologyItem={() => {
                }}/>}
        </Stack>
    );
};

export default OntologySandbox;