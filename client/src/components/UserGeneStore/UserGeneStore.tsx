import React, {useState} from 'react';
import {ActionIcon, Collapse, Group, Indicator, Stack, Text, useMantineTheme} from "@mantine/core";
import {AddGene} from "../AddGene/AddGene";
import {IconChevronDown, IconChevronRight, IconTrashX} from "@tabler/icons";
import {useRecoilState, useRecoilValue} from "recoil";
import {selectedGenesState, useGeneStoreCounterColor, userGenesState} from "../../atoms";
import UserGene from "../UserGene/UserGene";

interface Props {
    multiple?: boolean;
    opened?: boolean;
}

const UserGeneStore = ({multiple = false, opened = false}: Props) => {
    const [storeOpened, setOpened] = useState(opened);
    const indicatorColor = useRecoilValue(useGeneStoreCounterColor);
    const theme = useMantineTheme()
    const [userGeneStore, setUserGeneStore] = useRecoilState(userGenesState);
    const [selectedGenes, setSelectedGenes] = useRecoilState(selectedGenesState);
    return (
        <Stack>
            <AddGene multipleSelected={multiple}/>
            <Group onClick={() => setOpened(!storeOpened)} style={{cursor: 'pointer'}}>
                <Group spacing={0}>
                    <ActionIcon size='xs' variant={'subtle'}>
                        {storeOpened ? <IconChevronDown color={theme.colors.dark[9]}/> :
                            <IconChevronRight color={theme.colors.dark[9]}/>}
                    </ActionIcon>
                    <Indicator color={indicatorColor} position={"middle-end"} inline offset={-20}
                               label={`${userGeneStore.length}`} size={20}>
                        <Text size={'xs'}>User gene(s)</Text>
                    </Indicator>
                </Group>
            </Group>
            <Collapse in={storeOpened} transitionDuration={0} transitionTimingFunction="linear">
                {userGeneStore.length === 0 ? <Text size={'xs'} color={'dimmed'}>No genes added yet.</Text> :
                    <Stack> <Group>
                        <ActionIcon size={'xs'} onClick={() => {
                            setSelectedGenes([]);
                            setUserGeneStore([]);

                        }}>
                            <IconTrashX/>
                        </ActionIcon>
                        <Text size={'xs'} color={'dimmed'}>Remove all genes from store</Text>

                    </Group>
                        {userGeneStore.map((omics) => {
                            return (<Group align={'center'} position={'left'}>
                                <UserGene multiple={multiple} gene={omics}></UserGene>
                            </Group>)
                        })}</Stack>}
            </Collapse>


        </Stack>
    );
};

export {UserGeneStore};