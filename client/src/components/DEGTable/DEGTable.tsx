import DataTable from 'react-data-table-component';
import { ActionIcon, Center, Group, Loader, Stack, Text } from '@mantine/core';
import { IconChevronDown, IconChevronRight, IconEye, IconPlus } from '@tabler/icons-react';
import memoize from 'memoize-one';
import { useRecoilState, useRecoilValue } from 'recoil';
import _ from 'lodash';
import { showNotification } from '@mantine/notifications';
import { useCallback } from 'react';
import { selectedGenesState, studyState, userGenesState, userGeneStoreCounterColor, userGeneStoreOpenState } from '../../atoms';
import { Omics } from '../../model';
import { DifferentialExpressionV, DifferentialExpressionVFilter, useDegQuery } from '../../generated/types';

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

function LinkedGene({ gene, showExpression, addToStore }: { gene: Omics; showExpression: (gene: Omics) => void; addToStore: (gene: Omics) => void }) {
  return (
    <Group spacing="xs">
      <Text size="xs">{gene.displaySymbol}</Text>
      <ActionIcon title="add to gene store" color="blue.3" onClick={() => showExpression(gene as Omics)} size="xs" variant="default">
        <IconEye size={12} color="black" />
      </ActionIcon>
      <ActionIcon title="add to gene store" color="blue.3" onClick={() => addToStore(gene as Omics)} size="xs" variant="default">
        <IconPlus size={12} color="black" />
      </ActionIcon>
    </Group>
  );
}

function ExpandedComponent({ data }: { data: DifferentialExpressionV }) {
  const study = useRecoilValue(studyState);
  const [userGenes, setUserGenes] = useRecoilState(userGenesState);
  const [, setIndicatorColor] = useRecoilState(userGeneStoreCounterColor);
  const [, setStoreOpen] = useRecoilState(userGeneStoreOpenState);
  const [selectedGenes, setSelectedGenesStore] = useRecoilState(selectedGenesState);

  const showExpression = useCallback(
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

  const addToStore = useCallback(
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
    [setIndicatorColor, setUserGenes, userGenes, setStoreOpen],
  );

  const linkedGenes: Omics[] = data.linkedGenes
    .map((id: number) => {
      return study?.studyOmicsMap.get(id) as Omics;
    })
    .filter((g) => g !== undefined);

  return (
    <pre>
      <Center>
        <Stack>
          {linkedGenes && linkedGenes.length > 0 && (
            <Text weight={800} size="xs">
              Corresponding gene(s)
            </Text>
          )}
          {linkedGenes &&
            linkedGenes.length > 0 &&
            linkedGenes.map((gene) => {
              return <LinkedGene gene={gene} key={gene.omicsId} showExpression={showExpression} addToStore={addToStore} />;
            })}
        </Stack>
      </Center>
    </pre>
  );
}

export function DEGTable({ annotationId, selectedDEGComparisonAnnotationId }: { annotationId: number; selectedDEGComparisonAnnotationId: number }) {
  const [userGenes, setUserGenes] = useRecoilState(userGenesState);
  const [, setIndicatorColor] = useRecoilState(userGeneStoreCounterColor);
  const [selectedGenes, setSelectedGenesStore] = useRecoilState(selectedGenesState);
  // const annotationGroup = useRecoilValue(annotationGroupIdState);
  const study = useRecoilValue(studyState);
  // const [, setStoreOpen] = useRecoilState(userGeneStoreOpenState);

  const { data, loading } = useDegQuery({
    variables: {
      filter: {
        annotationValueId: { equalTo: annotationId },
        studyId: { equalTo: study?.studyId || 0 },
        otherAnnotationValueId: selectedDEGComparisonAnnotationId === 0 ? { isNull: true } : { equalTo: selectedDEGComparisonAnnotationId },
      } as unknown as DifferentialExpressionVFilter,
    },
    skip: annotationId === 0 || !study,
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
        // setStoreOpen(false);
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
    [setIndicatorColor, setUserGenes, userGenes],
  );

  return (
    <Stack justify="flex-start" align="left" w="100%">
      {data && study && data.differentialExpressionVsList.length > 0 && selectedDEGComparisonAnnotationId ? (
        <Text size="xs" color="dimmed">
          Showing differentially expressed genes of two selected annotations.
        </Text>
      ) : null}
      {/* TODO ExpandedComponent can also link from gene to protein, so for a multi-omics study all omics row types can be expanded */}
      {data && study && data.differentialExpressionVsList.length === 0 && (
        <Text size="xs">No differentially expressed genes found for the current group selection.</Text>
      )}
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
      {loading && (
        <Center>
          <Loader variant="dots" color="blue" />
        </Center>
      )}
    </Stack>
  );
}
