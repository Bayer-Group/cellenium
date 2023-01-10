import React from 'react';
import {ActionIcon, CloseButton, createStyles, Group, Text} from "@mantine/core";
import {IconDroplet} from "@tabler/icons";
import {Gene} from "../../model";
import {useRecoilState, useSetRecoilState} from "recoil";
import {pageState, selectedGenesState, userGenesState} from "../../atoms";

type Props = {
    gene: Gene
}
const useStyles = createStyles((theme) => ({
    active: {
        backgroundColor: theme.colors.blue[3]

    },
    main: {
        "&:hover": {
            backgroundColor: theme.colors.blue[1]
        }
    }
}));

const UserGene = ({gene}: Props) => {
    const {cx, classes} = useStyles();
    const [geneStore, setGeneStore] = useRecoilState(userGenesState);
    const [selectedGenes, setSelectedGenesStore] = useRecoilState(selectedGenesState);
    const setPage = useSetRecoilState(pageState);

    function handleRemove(gene: Gene) {

        // remove from geneStore
        let removed = geneStore.filter((g) => g.omicsId !== gene.omicsId)
        setGeneStore(removed)
        // remove from Selection
        removed = selectedGenes.filter((g) => g.omicsId !== gene.omicsId)
        setSelectedGenesStore(removed)
    }

    function handleColorClick(gene: Gene) {
        if (selectedGenes.filter((g) => g.omicsId === gene.omicsId).length > 0) {
            // remove
            let removed = selectedGenes.filter((g) => g.omicsId !== gene.omicsId)
            setSelectedGenesStore(removed)
        } else {
            // add
            setSelectedGenesStore([...selectedGenes, gene]);
            setPage('ExpressionAnalysis');
        }

    }

    return (
        <Group position={'apart'}>
            <Group>
                <CloseButton onClick={() => handleRemove(gene)} size={'xs'} iconSize={15}/>
                <Text size={'xs'}>{gene.displaySymbol}</Text>
            </Group>
            {/* eslint-disable-next-line react/jsx-no-undef */}
            <ActionIcon
                className={cx(classes.main, {[classes.active]: selectedGenes.filter((g) => g.omicsId == gene.omicsId).length !== 0})}
                onClick={() => handleColorClick(gene)} variant="default" size={'xs'} mr={5}>
                <IconDroplet/>
            </ActionIcon>
        </Group>
    );
};

export default UserGene;
