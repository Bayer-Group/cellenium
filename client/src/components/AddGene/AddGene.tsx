import {
    ActionIcon,
    Autocomplete,
    AutocompleteItem,
    Group,
    Stack,
    Text,
    TextInputProps,
    useMantineTheme
} from '@mantine/core';
import React, {FormEvent, useState} from "react";
import {IconArrowRight, IconX} from "@tabler/icons";
import {useRecoilState, useRecoilValue} from "recoil";
import {studyState, userGenesState} from "../../atoms";
import {showNotification} from '@mantine/notifications';
import {useForm} from '@mantine/form';
import * as aq from 'arquero';
import {Omics} from "../../model";




function AddGene(props: TextInputProps) {
    const [offerings, setOfferings] = useState<Omics[]>([])
    const [value, setValue] = useState('');
    const theme = useMantineTheme();
    const [userGenes, setUserGenes] = useRecoilState(userGenesState);
    const study = useRecoilValue(studyState);
    const form = useForm();

    function handleChange(inputString: string) {
        console.log({inputString})
        let newOfferings: Omics[] = [];
        if (inputString.length > 0) {
            // @ts-ignore
            newOfferings = study?.studyOmicsTable.filter(aq.escape(t => aq.op.includes(t.displaySymbol, inputString, 0))).objects();
        }
        console.log({newOfferings});
        setOfferings(newOfferings)
        setValue(inputString)
    }

    function handleItemSubmit(item: Omics) {
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
            setUserGenes([...userGenes, item])
        }

    }

    function handleSubmit(event: React.MouseEvent | FormEvent) {
        event.preventDefault();
        const addGene: Omics[] = offerings.filter((g) => g.displaySymbol.toLowerCase() === value.toLowerCase())

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

            setUserGenes([...userGenes, addGene[0]])
        }
    }

    return (
        <Stack spacing={0}>
            <Text size={'xs'}>
                Enter gene(s)
            </Text>
            <Group align={'center'} spacing={3}>
                <form onSubmit={(event) => handleSubmit(event)}>
                    <Autocomplete
                        value={value}
                        onChange={handleChange}
                        data={offerings as AutocompleteItem[]}
                        radius={'md'}
                        onItemSubmit={(item: any) => {
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
