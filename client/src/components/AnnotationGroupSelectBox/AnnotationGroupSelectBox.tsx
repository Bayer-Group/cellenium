import {Select} from '@mantine/core';
import {useState} from "react";

const DATA = [
                {value: 'cellOntology', label: 'Cell ontology names', group: 'Annotation groups'},
                {value: 'clusterName', label: 'Cluster names', group: 'Annotation groups'},
                {value: 'CellO Bayer celltype', label: 'CellO cell types', group: 'Advanced annotation groups'},
                {value: 'cellOntologyID', label: 'Cell ontology IDs', group: 'Advanced annotation groups'},
            ];
function AnnotationGroupSelectBox() {
    const [value, setValue] = useState<string|null>(DATA[0].value);

    return (
        <Select
            value={value} onChange={setValue}
            label="Select annotation group"
            placeholder="Pick one"
            transitionDuration={80}
            transitionTimingFunction="ease"
            data={DATA}
        />
    );
}

export {AnnotationGroupSelectBox};