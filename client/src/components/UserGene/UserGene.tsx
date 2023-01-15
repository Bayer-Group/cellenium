import React, {useState} from 'react';
import {ActionIcon, CloseButton, createStyles, Grid, Group, Text, Tooltip} from "@mantine/core";
import {IconEye, IconInfoCircle} from "@tabler/icons";
import {Omics} from "../../model";
import {useRecoilState, useSetRecoilState} from "recoil";
import {pageState, selectedGenesState, userGenesState} from "../../atoms";


const useStyles = createStyles((theme) => ({
    active: {
        backgroundColor: theme.colors.blue[6]

    }, main: {
        "&:hover": {
            backgroundColor: theme.colors.blue[4]
        }
    }

}));

interface Props {
    gene: Omics;
    multiple?: boolean;
}

const UserGene = ({gene, multiple = false}: Props) => {
    const {cx, classes} = useStyles();
    const [geneStore, setGeneStore] = useRecoilState(userGenesState);
    const [selectedGenes, setSelectedGenesStore] = useRecoilState(selectedGenesState);
    const setPage = useSetRecoilState(pageState);
    const [showInfo, setShowInfo] = useState(false)
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
            if (multiple)
                setSelectedGenesStore([...selectedGenes, gene]);
            else {
                setSelectedGenesStore([gene]);
            }
            //setPage('ExpressionAnalysis');
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
                        <IconEye style={{color: 'green'}}
                                 color={(selectedGenes.filter((g) => g.omicsId === gene.omicsId).length !== 0) ? 'white' : 'gray'}/>
                    </ActionIcon>
                </Grid.Col>
                <Grid.Col span={1}>
                    <Tooltip label={`${gene.displayName}`} opened={showInfo}>
                        <ActionIcon variant={'subtle'} size={'xs'} onClick={()=>setShowInfo(true)} onMouseLeave={()=>setShowInfo(false)}>
                            <IconInfoCircle/>
                        </ActionIcon>
                    </Tooltip>
                </Grid.Col>
                <Grid.Col span={5}>
                    <Text size={'xs'}>{gene.displaySymbol}</Text>
                </Grid.Col>


            </Grid>
        </Group>
    );
};

export default UserGene;
