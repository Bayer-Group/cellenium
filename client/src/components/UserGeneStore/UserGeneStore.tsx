import React from 'react';
import {ActionIcon, Collapse, Group, Indicator, Stack, Text, useMantineTheme} from "@mantine/core";
import {AddGene} from "../AddGene/AddGene";
import {IconChevronDown, IconChevronRight, IconTrashX} from "@tabler/icons";
import {useRecoilState, useRecoilValue} from "recoil";
import {selectedGenesState, userGenesState, userGeneStoreCounterColor, userGeneStoreOpenState} from "../../atoms";
import UserGene from "../UserGene/UserGene";

interface Props {
    multiple?: boolean;
    findCoexpressors?: boolean;
}

const UserGeneStore = ({multiple = false, findCoexpressors=false}: Props) => {
    const [storeOpened, setOpened] = useRecoilState(userGeneStoreOpenState);
    const indicatorColor = useRecoilValue(userGeneStoreCounterColor);
    const theme = useMantineTheme()
    const [userGeneStore, setUserGeneStore] = useRecoilState(userGenesState);
    const [selectedGenes, setSelectedGenes] = useRecoilState(selectedGenesState);
    return (
        <Stack>
            <AddGene multipleSelected={multiple}/>
            <Group onClick={() => {
                setOpened(!storeOpened)
            }} style={{cursor: 'pointer'}}>
                <Group spacing={0}>
                    <ActionIcon size='xs' variant={'subtle'}>
                        {storeOpened ? <IconChevronDown color={theme.colors.dark[9]}/> :
                            <IconChevronRight color={theme.colors.dark[9]}/>}
                    </ActionIcon>
                    <Indicator color={indicatorColor} position={"middle-end"} inline offset={-20}
                               label={`${userGeneStore.length}`} size={20}>
                        <Text size={'xs'}>Stored gene(s)</Text>
                    </Indicator>
                </Group>
            </Group>
            <Collapse in={storeOpened} transitionDuration={0} transitionTimingFunction="linear">
                {userGeneStore.length === 0 ? <Text size={'xs'} color={'dimmed'}>No genes added yet.</Text> :
                    <Stack> <Group position={'left'} spacing={'xs'}>
                        <Text size={'xs'} color={'dimmed'}>Remove all genes from store</Text>
                        <ActionIcon size={'xs'} onClick={() => {
                            setSelectedGenes([]);
                            setUserGeneStore([]);
                        }}>
                            <IconTrashX/>
                        </ActionIcon>


                    </Group>
                        {[...userGeneStore].reverse().map((omics) => {
                            return (<Group key={omics.omicsId} align={'center'} position={'left'}>
                                <UserGene multiple={multiple} gene={omics} findCoexpressors={findCoexpressors}></UserGene>
                            </Group>)
                        })}</Stack>}

            </Collapse>


        </Stack>
    );
};

export {UserGeneStore};