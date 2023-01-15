import React, {useState} from 'react';
import {ActionIcon, Collapse, Group, Indicator, Stack, Text, useMantineTheme} from "@mantine/core";
import {AddGene} from "../AddGene/AddGene";
import {IconChevronDown, IconChevronRight} from "@tabler/icons";
import {useRecoilValue} from "recoil";
import {userGenesState} from "../../atoms";
import UserGene from "../UserGene/UserGene";

interface Props {
    multiple?: boolean;
    opened?: boolean;
}

const UserGeneStore = ({multiple=false, opened=false}:Props) => {
    const [storeOpened, setOpened] = useState(opened);
    const theme = useMantineTheme()
    const userGeneStore = useRecoilValue(userGenesState);
    return (
        <Stack>
            <AddGene/>
            <Group onClick={() => setOpened(!storeOpened)} style={{cursor: 'pointer'}}>
                <Group spacing={0}><ActionIcon size='xs' variant={'subtle'}>
                    {storeOpened ? <IconChevronDown color={theme.colors.dark[9]}/> :
                        <IconChevronRight color={theme.colors.dark[9]}/>}
                </ActionIcon>
                    <Indicator position={"middle-end"} inline offset={-20} label={`${userGeneStore.length}`} size={20}>
                        <Text size={'xs'}>User gene(s)</Text>
                    </Indicator>
                </Group>
            </Group>
            <Collapse in={storeOpened} transitionDuration={0} transitionTimingFunction="linear">
                {userGeneStore.length===0?<Text size={'xs'} color={'dimmed'}>No genes added yet.</Text>:
                    userGeneStore.map((omics)=>{
                        return (<Group align={'center'} position={'left'}>
                            <UserGene multiple={multiple} gene={omics}></UserGene>
                        </Group>)
                    })}
            </Collapse>


        </Stack>
    );
};

export {UserGeneStore};