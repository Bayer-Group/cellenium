import { ActionIcon, Collapse, Group, Indicator, Stack, Text, useMantineTheme } from '@mantine/core';
import { IconChevronDown, IconChevronRight, IconTrashX } from '@tabler/icons-react';
import { useRecoilState, useRecoilValue } from 'recoil';
import { AddGene } from '../AddGene/AddGene';
import { selectedGenesState, userGenesState, userGeneStoreCounterColor, userGeneStoreOpenState } from '../../atoms';
import { UserGene } from '../UserGene/UserGene';

export function UserGeneStore({ multiple = false, findCoexpressors = false }: { multiple?: boolean; findCoexpressors?: boolean }) {
  const [storeOpened, setOpened] = useRecoilState(userGeneStoreOpenState);
  const indicatorColor = useRecoilValue(userGeneStoreCounterColor);
  const theme = useMantineTheme();
  const [userGeneStore, setUserGeneStore] = useRecoilState(userGenesState);
  const [, setSelectedGenes] = useRecoilState(selectedGenesState);
  return (
    <Stack w="100%">
      <AddGene multipleSelected={multiple} />
      <Group
        onClick={() => {
          setOpened(!storeOpened);
        }}
        style={{ cursor: 'pointer' }}
      >
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
          <Text size="xs" color="dimmed">
            No genes added yet.
          </Text>
        ) : (
          <Stack>
            {' '}
            <Group position="left" spacing="xs">
              <Text size="xs" color="dimmed">
                Remove everything from store
              </Text>
              <ActionIcon
                size="xs"
                onClick={() => {
                  setSelectedGenes([]);
                  setUserGeneStore([]);
                }}
              >
                <IconTrashX />
              </ActionIcon>
            </Group>
            {[...userGeneStore].reverse().map((omics) => {
              return (
                <Group key={omics.omicsId} align="center" position="left">
                  <UserGene multiple={multiple} gene={omics} findCoexpressors={findCoexpressors} />
                </Group>
              );
            })}
          </Stack>
        )}
      </Collapse>
    </Stack>
  );
}
