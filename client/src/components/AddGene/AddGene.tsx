import { ActionIcon, Autocomplete, AutocompleteItem, createStyles, Group, Stack, Text, useMantineTheme } from '@mantine/core';
import React, { FormEvent, useCallback, useState } from 'react';
import { IconArrowRight, IconX } from '@tabler/icons-react';
import { useRecoilState, useRecoilValue, useSetRecoilState } from 'recoil';
import { showNotification } from '@mantine/notifications';
import * as aq from 'arquero';
import { selectedGenesState, studyState, userGenesState, userGeneStoreCounterColor, userGeneStoreOpenState } from '../../atoms';
import { Omics } from '../../model';

const useStyles = createStyles(() => ({
  form: {
    width: '100%',
  },
  icon: {
    borderBottomLeftRadius: 0,
    borderTopLeftRadius: 0,
  },
}));

export function AddGene({ multipleSelected = false }: { multipleSelected?: boolean }) {
  const { classes } = useStyles();
  const [offerings, setOfferings] = useState<Omics[]>([]);
  const [value, setValue] = useState('');
  const theme = useMantineTheme();
  const [userGenes, setUserGenes] = useRecoilState(userGenesState);
  const [selectedGenes, setSelectedGenes] = useRecoilState(selectedGenesState);
  const setIndicatorColor = useSetRecoilState(userGeneStoreCounterColor);
  const setOpened = useSetRecoilState(userGeneStoreOpenState);
  const study = useRecoilValue(studyState);

  const handleChange = useCallback(
    (inputString: string) => {
      let newOfferings: Omics[] = [];
      if (inputString.length > 0) {
        newOfferings = study?.studyOmicsTable
          .filter(aq.escape((t: { displaySymbol: string }) => aq.op.startswith(t.displaySymbol.toLowerCase(), inputString.toLowerCase(), 0)))
          .objects() as Omics[];
        newOfferings = newOfferings.sort((a, b) => a.displaySymbol.localeCompare(b.displaySymbol));
      }
      setOfferings(newOfferings);
      setValue(inputString);
    },
    [study?.studyOmicsTable],
  );

  const handleItemSubmit = useCallback(
    (item: AutocompleteItem) => {
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
        if (multipleSelected) setSelectedGenes([...selectedGenes, item as Omics]);
        else setSelectedGenes([item as Omics]);
        setUserGenes([...userGenes, item as Omics]);
        setIndicatorColor('pink');
        setTimeout(() => {
          setIndicatorColor('blue');
        }, 100);
        setOpened(true);
      }
    },
    [multipleSelected, selectedGenes, setIndicatorColor, setOpened, setSelectedGenes, setUserGenes, userGenes],
  );

  const handleSubmit = useCallback(
    (event: React.MouseEvent | FormEvent) => {
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
    },
    [multipleSelected, offerings, selectedGenes, setOpened, setSelectedGenes, setUserGenes, userGenes, value],
  );

  return (
    <Stack spacing={0}>
      <Text size="xs">Enter identifiers(s)</Text>
      <Group spacing={0} noWrap>
        <form className={classes.form} onSubmit={(event) => handleSubmit(event)}>
          <Autocomplete
            value={value}
            radius={0}
            onChange={handleChange}
            data={offerings as AutocompleteItem[]}
            limit={15}
            onItemSubmit={handleItemSubmit}
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
          className={classes.icon}
          color={theme.primaryColor}
          variant="filled"
        >
          <IconArrowRight />
        </ActionIcon>
      </Group>
    </Stack>
  );
}
