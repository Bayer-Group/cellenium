import { ActionIcon, Autocomplete, AutocompleteItem, Group, Stack, Text, useMantineTheme } from '@mantine/core';
import React, { FormEvent, useState } from 'react';
import { IconArrowRight, IconX } from '@tabler/icons-react';
import { useRecoilState, useRecoilValue, useSetRecoilState } from 'recoil';
import { selectedGenesState, studyState, userGenesState, userGeneStoreCounterColor, userGeneStoreOpenState } from '../../atoms';
import { showNotification } from '@mantine/notifications';
import * as aq from 'arquero';
import { Omics } from '../../model';

interface Props {
  multipleSelected?: boolean;
}

function AddGene({ multipleSelected = false }: Props) {
  const [offerings, setOfferings] = useState<Omics[]>([]);
  const [value, setValue] = useState('');
  const theme = useMantineTheme();
  const [userGenes, setUserGenes] = useRecoilState(userGenesState);
  const [selectedGenes, setSelectedGenes] = useRecoilState(selectedGenesState);
  const setIndicatorColor = useSetRecoilState(userGeneStoreCounterColor);
  const setOpened = useSetRecoilState(userGeneStoreOpenState);

  const study = useRecoilValue(studyState);
  // const form = useForm();

  function handleChange(inputString: string) {
    let newOfferings: Omics[] = [];
    if (inputString.length > 0) {
      // @ts-ignore
      newOfferings = study?.studyOmicsTable
        .filter(aq.escape((t: any) => aq.op.startswith(t.displaySymbol.toLowerCase(), inputString.toLowerCase(), 0)))
        .objects();
      newOfferings = newOfferings.sort((a, b) => a.displaySymbol.localeCompare(b.displaySymbol));
    }
    setOfferings(newOfferings);
    setValue(inputString);
  }

  function handleItemSubmit(item: Omics) {
    setOfferings([]);
    setValue('');
    if (userGenes.filter((g) => g.omicsId === item.omicsId).length === 1) {
      showNotification({
        title: 'Your input is already in the store!',
        message: "It's not a problem, really!",
        color: 'red',
        autoClose: 1000,
      });
    } else {
      if (multipleSelected) setSelectedGenes([...selectedGenes, item]);
      else setSelectedGenes([item]);
      setUserGenes([...userGenes, item]);
      setIndicatorColor('pink');
      setTimeout(() => {
        setIndicatorColor('blue');
      }, 100);
      setOpened(true);
    }
  }

  function handleSubmit(event: React.MouseEvent | FormEvent) {
    event.preventDefault();
    const addGene: Omics[] = offerings.filter((g) => g.displaySymbol.toLowerCase() === value.toLowerCase());

    if (value === '') return;

    setValue('');
    setOfferings([]);
    if (addGene.length === 0) {
      showNotification({
        title: 'Provide a valid selection!',
        message: 'Please choose from the autocompletion list!',
        color: 'red',
        autoClose: 5000,
      });
    } else if (userGenes.filter((g) => g.omicsId === addGene[0].omicsId).length === 1) {
      showNotification({
        title: 'Your input is already in the store!',
        message: "It's not a problem, really!",
        color: 'red',
        autoClose: 1000,
      });
    } else if (addGene.length === 1) {
      if (multipleSelected) setSelectedGenes([...selectedGenes, ...addGene]);
      else setSelectedGenes(addGene);
      setUserGenes([...userGenes, ...addGene]);
      setOpened(true);
    }
  }

  return (
    <Stack spacing={0}>
      <Text size={'xs'}>Enter identifiers(s)</Text>
      <Group align={'center'} spacing={0}>
        <form onSubmit={(event) => handleSubmit(event)}>
          <Autocomplete
            value={value}
            radius={0}
            onChange={handleChange}
            data={offerings as AutocompleteItem[]}
            limit={15}
            onItemSubmit={(item: any) => {
              handleItemSubmit(item);
            }}
            rightSection={
              <ActionIcon
                onClick={() => {
                  setValue('');
                  setOfferings([]);
                }}
              >
                <IconX />
              </ActionIcon>
            }
          />
        </form>

        <ActionIcon
          onClick={(event) => {
            handleSubmit(event);
          }}
          size={36}
          style={{ borderBottomLeftRadius: 0, borderTopLeftRadius: 0 }}
          color={theme.primaryColor}
          variant="filled"
        >
          <IconArrowRight />
        </ActionIcon>
      </Group>
    </Stack>
  );
}

export { AddGene };
