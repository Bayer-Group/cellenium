import { ActionIcon, Collapse, createStyles, Divider, Group, Indicator, Stack, Text, useMantineTheme } from '@mantine/core';
import { IconChevronDown, IconChevronRight, IconTrashX } from '@tabler/icons-react';
import { useRecoilState, useRecoilValue } from 'recoil';
import { useCallback, useMemo } from 'react';
import { AddGene } from '../AddGene/AddGene';
import { selectedGenesState, userGenesState, userGeneStoreCounterColor, userGeneStoreOpenState } from '../../atoms';
import { UserGene } from '../UserGene/UserGene';

const useStyles = createStyles(() => ({
  cursor: {
    cursor: 'pointer',
  },
}));

export function UserGeneStore({ multiple = false, findCoexpressors = false }: { multiple?: boolean; findCoexpressors?: boolean }) {
  const { classes } = useStyles();
  const [storeOpened, setOpened] = useRecoilState(userGeneStoreOpenState);
  const indicatorColor = useRecoilValue(userGeneStoreCounterColor);
  const theme = useMantineTheme();
  const [userGeneStore, setUserGeneStore] = useRecoilState(userGenesState);
  const [, setSelectedGenes] = useRecoilState(selectedGenesState);

  const handleOpened = useCallback(() => {
    setOpened(!storeOpened);
  }, [setOpened, storeOpened]);

  const resetStore = useCallback(() => {
    setSelectedGenes([]);
    setUserGeneStore([]);
  }, [setSelectedGenes, setUserGeneStore]);

  const userGenes = useMemo(() => {
    return [...userGeneStore].reverse().map((omics) => {
      return (
        <Group key={omics.omicsId} align="center" position="left">
          <UserGene multiple={multiple} gene={omics} findCoexpressors={findCoexpressors} />
        </Group>
      );
    });
  }, [findCoexpressors, multiple, userGeneStore]);

  return (
    <Stack spacing="sm">
      <Divider size="xs" label="User gene store" />
      <AddGene multipleSelected={multiple} />
      <Group onClick={handleOpened} className={classes.cursor}>
        <Group spacing={0}>
          <ActionIcon size="xs" variant="subtle">
            {storeOpened ? <IconChevronDown color={theme.colors.dark[9]} /> : <IconChevronRight color={theme.colors.dark[9]} />}
          </ActionIcon>
          <Indicator color={indicatorColor} position="middle-end" inline offset={-20} label={`${userGeneStore.length}`} size={20}>
            <Text size="xs">Stored genes/proteins/regions</Text>
          </Indicator>
        </Group>
      </Group>
      <Collapse in={storeOpened} transitionDuration={0} transitionTimingFunction="linear">
        {userGeneStore.length === 0 ? (
          <Text size="xs" color="dimmed" pl="md">
            No genes added yet.
          </Text>
        ) : (
          <Stack pl="md">
            <Group position="left" spacing="xs">
              <Text size="xs" color="dimmed">
                Remove everything from store
              </Text>
              <ActionIcon size="xs" onClick={resetStore}>
                <IconTrashX />
              </ActionIcon>
            </Group>
            {userGenes}
          </Stack>
        )}
      </Collapse>
    </Stack>
  );
}
