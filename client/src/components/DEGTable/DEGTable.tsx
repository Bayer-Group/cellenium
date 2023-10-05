import DataTable from 'react-data-table-component';
import { ActionIcon, Center, Group, Loader, Stack, Text } from '@mantine/core';
import { IconChevronDown, IconChevronRight, IconEye, IconPlus } from '@tabler/icons-react';
import memoize from 'memoize-one';
import { useRecoilState, useRecoilValue } from 'recoil';
import _ from 'lodash';
import { showNotification } from '@mantine/notifications';
import { ReactNode, useCallback } from 'react';
import { selectedGenesState, studyState, userGenesState, userGeneStoreCounterColor, userGeneStoreOpenState } from '../../atoms';
import { Omics } from '../../model';
import { DifferentialExpressionV, useDegQuery } from '../../generated/types';

const customStyles = {
  table: {
    style: {
      backgroundColor: 'transparent',
      marginRight: '10',
      overflow: 'hidden',
    },
  },
  header: {
    style: {
      paddingLeft: '2px',
      backgroundColor: 'transparent',
    },
  },
  head: {
    style: {
      paddingLeft: '2px',
      backgroundColor: 'transparent',
    },
  },
  headRow: {
    style: {
      paddingLeft: '2px',
      backgroundColor: 'transparent',
    },
  },
  rows: {
    style: {
      minHeight: '72px', // override the row height
      backgroundColor: 'transparent',
    },
  },
  headCells: {
    style: {
      paddingLeft: '2px', // override the cell padding for head cells
      paddingRight: '2px',
      backgroundColor: 'transparent',
    },
  },
  cells: {
    style: {
      paddingLeft: '2px', // override the cell padding for data cells
      paddingRight: '2px',
      backgroundColor: 'transparent',
    },
  },
};
const columns = memoize((clickHandler, handleColorClick) => [
  {
    name: 'identifier',
    cell: (row: DifferentialExpressionV) => <Text title={row.displaySymbol}>{row.displaySymbol}</Text>,
    sortable: true,
    width: '80px',
  },
  {
    name: 'padj',
    selector: (row: DifferentialExpressionV) => row.pvalueAdj.toFixed(2),
    sortable: true,
    width: '50px',
  },
  {
    name: 'lgFC',
    selector: (row: DifferentialExpressionV) => {
      let ret;
      if (row.log2Foldchange !== -1) ret = +row.log2Foldchange.toFixed(2);
      else ret = '';
      return ret;
    },
    sortable: true,
    width: '50px',
  },
  {
    name: '',
    cell: (row: DifferentialExpressionV) => {
      const gene = {
        omicsId: row.omicsId,
        displayName: row.displayName,
        displaySymbol: row.displaySymbol,
        omicsType: row.omicsType,
        value: row.displaySymbol,
        linkedGenes: row.linkedGenes,
      };

      return (
        <Group position="center" align="center" noWrap spacing={0}>
          <ActionIcon title="superpose expression" onClick={() => handleColorClick(gene)} variant="default" size="xs" mr={5}>
            <IconEye />
          </ActionIcon>
          <ActionIcon title="add to gene store" color="blue.3" onClick={() => clickHandler(gene)} size="xs" variant="default">
            <IconPlus size={12} color="black" />
          </ActionIcon>
        </Group>
      );
    },
    width: '20px',
  },
]);

function ExpandedComponent({ data }: { data: DifferentialExpressionV }) {
  const study = useRecoilValue(studyState);
  const [userGenes, setUserGenes] = useRecoilState(userGenesState);
  const [, setIndicatorColor] = useRecoilState(userGeneStoreCounterColor);
  const [, setStoreOpen] = useRecoilState(userGeneStoreOpenState);
  const [selectedGenes, setSelectedGenesStore] = useRecoilState(selectedGenesState);

  function showExpression(gene: Omics) {
    if (selectedGenes.filter((g) => g.omicsId === gene.omicsId).length > 0) {
      // remove
      const removed = selectedGenes.filter((g) => g.omicsId !== gene.omicsId);
      setSelectedGenesStore(removed);
    } else {
      setSelectedGenesStore([gene]);
    }
  }

  function addToStore(gene: Omics) {
    const check = userGenes.filter((g) => g.omicsId === gene.omicsId);
    if (check.length === 0) {
      setIndicatorColor('pink');
      setUserGenes(_.union(userGenes, [gene]));
      setStoreOpen(false);
      setTimeout(() => {
        setIndicatorColor('blue');
      }, 200);
    } else {
      showNotification({
        title: 'Your selection is already in the store',
        message: '',
        color: 'red',
        autoClose: 1000,
      });
    }
  }

  const linkedGenes: ReactNode[] = data.linkedGenes.map((id: number) => {
    const gene = study?.studyOmicsMap.get(id);
    if (gene !== undefined)
      return (
        <Group key={id} spacing="xs">
          <Text size="xs">{gene.displaySymbol}</Text>
          <ActionIcon title="add to gene store" color="blue.3" onClick={() => showExpression(gene as Omics)} size="xs" variant="default">
            <IconEye size={12} color="black" />
          </ActionIcon>
          <ActionIcon title="add to gene store" color="blue.3" onClick={() => addToStore(gene as Omics)} size="xs" variant="default">
            <IconPlus size={12} color="black" />
          </ActionIcon>
        </Group>
      );
    return null;
  });
  return (
    <pre>
      <Center>
        <Stack>
          {linkedGenes && linkedGenes.length > 0 && (
            <Text weight={800} size="xs">
              Corresponding gene(s)
            </Text>
          )}
          {linkedGenes && linkedGenes.length > 0 && linkedGenes}
        </Stack>
      </Center>
    </pre>
  );
}

export function DEGTable({ annotationId }: { annotationId: number }) {
  const [userGenes, setUserGenes] = useRecoilState(userGenesState);
  const [, setIndicatorColor] = useRecoilState(userGeneStoreCounterColor);
  const [selectedGenes, setSelectedGenesStore] = useRecoilState(selectedGenesState);
  // const annotationGroup = useRecoilValue(annotationGroupIdState);
  const study = useRecoilValue(studyState);
  const [, setStoreOpen] = useRecoilState(userGeneStoreOpenState);

  const { data, loading } = useDegQuery({
    variables: {
      annotationValueId: annotationId,
      studyId: study?.studyId || 0,
    },
  });

  const handleColorClick = useCallback(
    (gene: Omics) => {
      if (selectedGenes.filter((g) => g.omicsId === gene.omicsId).length > 0) {
        // remove
        const removed = selectedGenes.filter((g) => g.omicsId !== gene.omicsId);
        setSelectedGenesStore(removed);
      } else {
        setSelectedGenesStore([gene]);
      }
    },
    [selectedGenes, setSelectedGenesStore],
  );

  const handleClick = useCallback(
    (gene: Omics) => {
      const check = userGenes.filter((g) => g.omicsId === gene.omicsId);
      if (check.length === 0) {
        setIndicatorColor('pink');
        setUserGenes(_.union(userGenes, [gene]));
        setStoreOpen(false);
        setTimeout(() => {
          setIndicatorColor('blue');
        }, 200);
      } else {
        showNotification({
          title: 'Your selection is already in the store',
          message: '',
          color: 'red',
          autoClose: 1000,
        });
      }
    },
    [setIndicatorColor, setStoreOpen, setUserGenes, userGenes],
  );

  return (
    <Stack justify="flex-start" align="center" w="100%">
      {/* TODO ExpandedComponent can also link from gene to protein, so for a multi-omics study all omics row types can be expanded */}
      {data && study && data.differentialExpressionVsList.length > 0 && study.omicsTypes.length > 1 && (
        <DataTable
          dense
          columns={columns(handleClick, handleColorClick)}
          data={data.differentialExpressionVsList}
          defaultSortFieldId={3}
          defaultSortAsc={false}
          customStyles={customStyles}
          fixedHeader
          fixedHeaderScrollHeight="100%"
          noDataComponent={<Text>No data.</Text>}
          expandableIcon={{
            collapsed: <IconChevronRight size={10} />,
            expanded: <IconChevronDown size={10} />,
          }}
          expandableRows
          expandableRowsComponent={ExpandedComponent}
        />
      )}
      {data && study && data.differentialExpressionVsList.length > 0 && study.omicsTypes.length === 1 && (
        <DataTable
          dense
          columns={columns(handleClick, handleColorClick)}
          data={data.differentialExpressionVsList}
          defaultSortFieldId={3}
          defaultSortAsc={false}
          customStyles={customStyles}
          fixedHeader
          fixedHeaderScrollHeight="100%"
          noDataComponent={<Text>No data.</Text>}
        />
      )}
      {loading && <Loader variant="dots" color="blue" />}
    </Stack>
  );
}
