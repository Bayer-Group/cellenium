import React, {useState} from 'react';
import {Select} from "@mantine/core";

const SPECIES = [
    {value: "9606", label: 'Homo sapiens'},
    {value: "10090", label: 'Mus musculus'},
    {value: "10116", label: 'Rattus norvegicus'},
]
const SpeciesSelect = () => {
    const [value, setValue] = useState(SPECIES[0].value)
    return (
        <Select style={{borderColor: '#000'}} onChange={
            (value: string) => setValue(value)
        } value={value} variant={'default'} size={'md'} data={SPECIES as any} label={'Select species'}
                labelProps={{fw: 700, size: 'xs'}}/>
    );
};

export {SpeciesSelect};