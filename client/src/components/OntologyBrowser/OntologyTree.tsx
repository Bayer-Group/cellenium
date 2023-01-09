import React from 'react';
import DATA from "../../data/ontology.json";
import {Container, Group, Stack, Text} from "@mantine/core";
import {OntologyItem} from "../../model";
import OntologyBranch from "./OntologyBranch";

type Props = {
    ontid: number;
}


const OntologyTree = ({ontid}: Props) => {
    return (
        <Group position={'center'} align={'center'} grow>
            <Stack spacing={0} justify={'center'} align={'flex-start'}>
                <Text>Ontology Browser</Text>
                {DATA.map((item: OntologyItem) => <OntologyBranch key={item.id} item={item} level={0}></OntologyBranch>)}
            </Stack>
        </Group>
    );
};

export {OntologyTree};