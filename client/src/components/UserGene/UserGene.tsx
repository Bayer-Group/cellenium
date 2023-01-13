import React from 'react';
import {ActionIcon, CloseButton, createStyles, Grid, Group, Text} from "@mantine/core";
import {IconEye} from "@tabler/icons";
import {Omics} from "../../model";
import {useRecoilState, useSetRecoilState} from "recoil";
import {pageState, selectedGenesState, userGenesState} from "../../atoms";


const useStyles = createStyles((theme) => ({
    active: {
        backgroundColor: theme.colors.blue[3]

    }, main: {
        "&:hover": {
            backgroundColor: theme.colors.blue[1]
        }
    }

}));

interface Props {
    gene: Omics;
    buttons?: string[];
}

const UserGene = ({gene, buttons = []}: Props) => {
    const {cx, classes} = useStyles();
    const [geneStore, setGeneStore] = useRecoilState(userGenesState);
    const [selectedGenes, setSelectedGenesStore] = useRecoilState(selectedGenesState);
    const setPage = useSetRecoilState(pageState);

    function handleRemove(gene: Omics) {
        // remove from geneStore
        let removed = geneStore.filter((g) => g.omicsId !== gene.omicsId)
        setGeneStore(removed)
        // remove from Selection
        removed = selectedGenes.filter((g) => g.omicsId !== gene.omicsId)
        setSelectedGenesStore(removed)
    }

    function handleColorClick(gene: Omics) {
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
        <Group align={'center'} position={'left'} spacing={2} style={{width: '100%'}}>
            <Grid columns={12} style={{width: '100%'}} gutter={'md'}>
                <Grid.Col span={1}>
                    <ActionIcon variant={'subtle'} size={'xs'}>
                        <CloseButton onClick={() => handleRemove(gene)} size={'xs'} iconSize={15}/>
                    </ActionIcon>
                </Grid.Col>
                <Grid.Col span={1}>
                    <ActionIcon
                        className={cx(classes.main, {[classes.active]: selectedGenes.filter((g) => g.omicsId === gene.omicsId).length !== 0})}
                        onClick={() => handleColorClick(gene)} variant="subtle" size={'xs'} mr={5}>
                        <IconEye/>
                    </ActionIcon>
                </Grid.Col>
                <Grid.Col span={3}>
                    <Text size={'xs'}>{gene.displaySymbol}</Text>
                </Grid.Col>
                <Grid.Col span={7}>
                    <Text size={'xs'} lineClamp={1}>{gene.displayName}</Text>
                </Grid.Col>

            </Grid>
        </Group>
    );
};

export default UserGene;
