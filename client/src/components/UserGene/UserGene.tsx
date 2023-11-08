import { useCallback, useState } from 'react';
import { ActionIcon, createStyles, Grid, Group, Text, Tooltip } from '@mantine/core';
import { IconCalculator, IconEye, IconInfoCircle, IconX } from '@tabler/icons-react';
import { useRecoilState } from 'recoil';
import { Omics } from '../../model';
import { correlationOmicsIdState, selectedGenesState, userGenesState } from '../../atoms';

const useStyles = createStyles((theme) => ({
  active: {
    backgroundColor: theme.colors.blue[6],
  },
  main: {
    '&:hover': {
      backgroundColor: theme.colors.blue[4],
    },
  },
}));

export function UserGene({ gene, multiple = false, findCoexpressors = false }: { gene: Omics; multiple?: boolean; findCoexpressors?: boolean }) {
  const { cx, classes } = useStyles();
  const [geneStore, setGeneStore] = useRecoilState(userGenesState);
  const [selectedGenes, setSelectedGenesStore] = useRecoilState(selectedGenesState);
  const [showInfo, setShowInfo] = useState(false);
  const [correlationOmicsId, setCorrelationOmicsId] = useRecoilState(correlationOmicsIdState);

  // const study = useRecoilValue(studyState);

  const handleRemove = useCallback(() => {
    // remove from geneStore
    let removed = geneStore.filter((ge) => ge.omicsId !== gene.omicsId);
    setGeneStore(removed);
    // remove from Selection
    removed = selectedGenes.filter((ge) => ge.omicsId !== gene.omicsId);
    setSelectedGenesStore(removed);
  }, [gene.omicsId, geneStore, selectedGenes, setGeneStore, setSelectedGenesStore]);

  const searchCoexpressors = useCallback(() => {
    setCorrelationOmicsId(gene.omicsId);
    setSelectedGenesStore([...selectedGenes, gene]);
  }, [gene, selectedGenes, setCorrelationOmicsId, setSelectedGenesStore]);

  const handleColorClick = useCallback(() => {
    if (selectedGenes.filter((g) => g.omicsId === gene.omicsId).length > 0) {
      // remove
      const removed = selectedGenes.filter((g) => g.omicsId !== gene.omicsId);
      setSelectedGenesStore(removed);
    } else if (multiple) {
      // add
      setSelectedGenesStore([...selectedGenes, gene]);
    } else {
      setSelectedGenesStore([gene]);
    }
  }, [gene, multiple, selectedGenes, setSelectedGenesStore]);

  const onInfoClick = useCallback(() => {
    setShowInfo((o) => !o);
  }, []);

  const onInfoLeave = useCallback(() => {
    setShowInfo(false);
  }, []);

  return (
    <Group align="center" position="left" spacing={2} w="100%">
      <Grid columns={12} w="100%" gutter="md">
        <Grid.Col span={1}>
          <ActionIcon variant="subtle" size="xs" title="remove from user gene list">
            <IconX onClick={handleRemove} />
          </ActionIcon>
        </Grid.Col>
        <Grid.Col span={1}>
          <ActionIcon
            className={cx(classes.main, {
              [classes.active]: selectedGenes.filter((g) => g.omicsId === gene.omicsId).length !== 0,
            })}
            onClick={handleColorClick}
            variant="subtle"
            size="xs"
            title="show/hide gene in plot"
            mr={5}
          >
            <IconEye color={selectedGenes.filter((g) => g.omicsId === gene.omicsId).length !== 0 ? 'white' : 'gray'} />
          </ActionIcon>
        </Grid.Col>
        {findCoexpressors && (
          <Grid.Col span={1}>
            <ActionIcon
              className={cx(classes.main, {
                [classes.active]: gene.omicsId === correlationOmicsId,
              })}
              variant="subtle"
              size="xs"
              title="find other genes with correlating expression"
              onClick={searchCoexpressors}
            >
              <IconCalculator color={correlationOmicsId === gene.omicsId ? 'white' : 'gray'} />
            </ActionIcon>
          </Grid.Col>
        )}
        <Grid.Col span={1}>
          <Tooltip label={`${gene.displayName}`} opened={showInfo}>
            <ActionIcon variant="subtle" size="xs" onClick={onInfoClick} onMouseLeave={onInfoLeave}>
              <IconInfoCircle />
            </ActionIcon>
          </Tooltip>
        </Grid.Col>

        <Grid.Col span={4}>
          <Group noWrap spacing={1}>
            <Text size="xs">{gene.displaySymbol}</Text>
            <Text size={8}>({gene.omicsType})</Text>
          </Group>
        </Grid.Col>
      </Grid>
    </Group>
  );
}
