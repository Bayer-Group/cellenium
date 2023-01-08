import {ActionIcon, Autocomplete, Group, Stack, Text, TextInputProps, useMantineTheme} from '@mantine/core';
import React, {FormEvent, useState} from "react";
import _ from 'lodash';
import {Gene} from "../../model";
import {IconArrowRight, IconX} from "@tabler/icons";
import {useRecoilState} from "recoil";
import {userGenesState} from "../../atoms";
import {showNotification} from '@mantine/notifications';
import {useForm} from '@mantine/form';

let GENES: Gene[] = [
    {'displayName': 'CDK2', 'displaySymbol': 'CDK2', omicsId: 1},
    {'displayName': 'ATAD2', 'displaySymbol': 'ATAD2', omicsId: 2},
    {'displayName': 'CDK3', 'displaySymbol': 'CDK3', omicsId: 3},
    {'displayName': 'BRAF', 'displaySymbol': 'BRAF', omicsId: 4},
    {'displayName': 'KRAS', 'displaySymbol': 'KRAS', omicsId: 5},
    {'displayName': 'PTK2', 'displaySymbol': 'PTK2', omicsId: 6},
    {'displayName': 'ADAR1', 'displaySymbol': 'ADAR1', omicsId: 7},
    {'displayName': 'CDK6', 'displaySymbol': 'CDK6', omicsId: 8},
    {'displayName': 'CDK5', 'displaySymbol': 'CDK5', omicsId: 9},
    {'displayName': 'CDK4', 'displaySymbol': 'CDK4', omicsId: 10},
];

GENES = _.orderBy(GENES, ['displaySymbol'], ['asc']);

type OfferingGene = {
    displaySymbol: string,
    displayName: string,
    omicsId: number,
    value: string
}

function AddGene(props: TextInputProps) {
    const [offerings, setOfferings] = useState<OfferingGene[]>([])
    const [value, setValue] = useState('');
    const theme = useMantineTheme();
    const [userGenes, setUserGenes] = useRecoilState(userGenesState);

    const form = useForm();

    function handleChange(input: string) {

        let newOfferings: any = [];
        if (input.length > 0) {
            for (let g of GENES) {
                if (g.displaySymbol.toLowerCase().startsWith(input.toLowerCase())) {
                    newOfferings.push({
                        ...g,
                        value: g.displaySymbol
                    })
                }
                if (newOfferings.length === 5) {
                    break
                }
            }
        }
        setOfferings(newOfferings)
        setValue(input)
    }

    function handleItemSubmit(item: OfferingGene) {
        setOfferings([])
        setValue('')
        if (userGenes.filter((g) => g.omicsId === item.omicsId).length === 1) {
            showNotification({
                title: 'Your input is already in the store!',
                message: "It's not a problem, really!",
                color: 'red',
                autoClose: 5000
            })
        } else {
            setUserGenes([...userGenes, {
                displaySymbol: item.displaySymbol,
                displayName: item.displayName,
                omicsId: item.omicsId
            }])

        }

    }

    function handleSubmit(event: React.MouseEvent | FormEvent) {
        event.preventDefault();
        const addGene: OfferingGene[] = offerings.filter((g) => g.displaySymbol.toLowerCase() === value.toLowerCase())

        if (value === '')
            return

        setValue('');
        setOfferings([]);
        if (addGene.length === 0) {
            showNotification({
                title: 'Provide a valid selection!',
                message: 'Please choose from the autocompletion list!',
                color: 'red',
                autoClose: 5000
            })
        } else if (userGenes.filter((g) => g.omicsId === addGene[0].omicsId).length === 1) {
            showNotification({
                title: 'Your input is already in the store!',
                message: "It's not a problem, really!",
                color: 'red',
                autoClose: 5000
            })
        } else if (addGene.length === 1) {

            setUserGenes([...userGenes, {
                displaySymbol: addGene[0].displaySymbol,
                displayName: addGene[0].displayName,
                omicsId: addGene[0].omicsId
            }])
        }
    }

    return (

        <Stack spacing={0}>
            <Text weight={800} size={'xs'}>
                Enter gene(s)
            </Text>
            <Group align={'center'} spacing={3}>
                <form onSubmit={(event) => handleSubmit(event)}>
                    <Autocomplete
                        value={value}
                        onChange={handleChange}
                        data={offerings}
                        radius={'md'}
                        onItemSubmit={(item: OfferingGene) => {
                            handleItemSubmit(item)
                        }}
                        rightSection={
                            <ActionIcon onClick={() => {
                                setValue('')
                                setOfferings([])
                            }
                            }>
                                <IconX/>
                            </ActionIcon>
                        }
                    />
                </form>

                <ActionIcon onClick={(event) => {
                    handleSubmit(event);
                }} size={32} radius="md" color={theme.primaryColor}
                            variant="filled">
                    <IconArrowRight/>
                </ActionIcon>
            </Group>
        </Stack>
    );
}

export {AddGene}

/*
                        <ActionIcon onClick={(event) => {
                            handleClickSubmit(event);
                            setValue('');
                            setOfferings([]);
                        }} size={32} radius="md" color={theme.primaryColor}
                                    variant="filled">
                            <IconPlus/>
                        </ActionIcon>
 */