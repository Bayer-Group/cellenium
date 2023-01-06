import React from 'react';
import {ActionIcon, CloseButton, Group, Text} from "@mantine/core";
import {IconDroplet} from "@tabler/icons";
import {Gene} from "../../model";
import {useRecoilState} from "recoil";
import {selectedGenesState, userGenesState} from "../../atoms";
import {useNavigate} from "react-router-dom";

type Props = {
    gene: Gene
}
const UserGene = ({gene}: Props) => {
    const [geneStore, setGeneStore] = useRecoilState(userGenesState);
    const [selectedGenes, setSelectedGenesStore] = useRecoilState(selectedGenesState);
    const navigate = useNavigate()
    function handleRemove(gene: Gene) {
        let removed = geneStore.filter((g)=>g.omicsId!==gene.omicsId)
        setGeneStore(removed)
    }
    function handleColorClick(gene: Gene) {
        setSelectedGenesStore([...selectedGenes,gene])
        navigate('/expressionanalysis?plotType=projectionplot')
    }
    return (
        <Group position={'apart'}>
            <Group>
                <CloseButton onClick={()=>handleRemove(gene)} size={'xs'} iconSize={15}/>
                <Text size={'xs'}>{gene.displaySymbol}</Text>
            </Group>
            {/* eslint-disable-next-line react/jsx-no-undef */}
            <ActionIcon onClick={()=>handleColorClick(gene)} variant="default" size={'xs'} mr={5}>
                <IconDroplet/>
            </ActionIcon>
        </Group>
    );
};

export default UserGene;
