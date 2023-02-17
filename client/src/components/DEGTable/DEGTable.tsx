import React from 'react';
import DataTable from "react-data-table-component";
import {ActionIcon, Group, Stack, Text} from "@mantine/core";
import {IconEye, IconPlus} from "@tabler/icons";
import {useDegQuery} from "../../generated/types";
import memoize from 'memoize-one';
import {useRecoilState, useRecoilValue} from "recoil";
import {
    selectedGenesState,
    studyState,
    userGenesState,
    userGeneStoreCounterColor,
    userGeneStoreOpenState
} from "../../atoms";
import {Omics} from "../../model";
import _ from 'lodash';
import {showNotification} from "@mantine/notifications";

const customStyles = {
    table: {
        style: {
            backgroundColor: 'transparent',
            marginRight: '10',
            overflow: 'hidden'
        }
    },
    header: {
        style: {
            paddingLeft: '2px',
            backgroundColor: 'transparent',

        }
    },
    head: {
        style: {
            paddingLeft: '2px',
            backgroundColor: 'transparent',

        }
    },
    headRow: {
        style: {
            paddingLeft: '2px',
            backgroundColor: 'transparent',
        }
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
        name: 'gene',
        selector: (row: any) => row.displaySymbol,
        sortable: true,
        width: '80px'
    },
    {
        name: 'padj',
        selector: (row: any) => row.pvalueAdj.toFixed(2),
        sortable: true,
        width: '70px'
    },
    {
        name: 'log2FC',
        selector: (row: any) => +row.log2Foldchange.toFixed(2),
        sortable: true,
        width: '70px'
    },
    {
        name: '',
        cell: (row: any) => {
            let gene = {
                omicsId: row.omicsId,
                displayName: row.displayName,
                displaySymbol: row.displaySymbol,
                omicsType: row.omicsType,
                value: row.displaySymbol
            }
            return (
                <Group position={'center'} align={'center'} noWrap={true} spacing={0}>
                    <ActionIcon title={'superpose expression'}
                        onClick={() => handleColorClick(gene)} variant="default" size={'xs'} mr={5}>
                        <IconEye/>
                    </ActionIcon>
                    <ActionIcon title={'add to gene store'} color={'blue.3'} onClick={() => clickHandler(gene)} size='xs'
                                variant={"default"}><IconPlus
                        size={12} color={'black'}/></ActionIcon>
                </Group>
            )
        },
        width: '20px',
    }
]);

type Props = {
    annotationId: number;
}

const DEGTable = ({annotationId}: Props) => {
    const [userGenes, setUserGenes] = useRecoilState(userGenesState);
    const [indicatorColor, setIndicatorColor] = useRecoilState(userGeneStoreCounterColor);
    const [selectedGenes, setSelectedGenesStore] = useRecoilState(selectedGenesState);

    const study = useRecoilValue(studyState);
    const [storeOpen, setStoreOpen] = useRecoilState(userGeneStoreOpenState)
    const {data, error, loading} = useDegQuery({
        variables: {
            annotationValueId: annotationId,
            studyId: study?.studyId || 0
        }
    })

    function handleColorClick(gene: Omics) {
        if (selectedGenes.filter((g) => g.omicsId === gene.omicsId).length > 0) {
            // remove
            let removed = selectedGenes.filter((g) => g.omicsId !== gene.omicsId)
            setSelectedGenesStore(removed)
        } else {
            setSelectedGenesStore([gene]);
        }

    }

    function handleClick(gene: Omics) {
        let check = userGenes.filter((g) => g.omicsId === gene.omicsId)
        if (check.length === 0) {
            setIndicatorColor('pink')
            setUserGenes(_.union(userGenes, [gene]))
            setStoreOpen(false)
            setTimeout(() => {
                setIndicatorColor('blue')
            }, 200)

        } else {
            showNotification({
                title: 'Your selection is already in the store',
                message: "It's not a problem, really!",
                color: 'red',
                autoClose: 1000
            })
        }
    }

    return (
        <Stack justify={'flex-start'} align={'center'} w={'100%'}>
            {data && data.differentialExpressionVsList.length > 0 &&
                <DataTable dense columns={columns(handleClick, handleColorClick)}
                           data={data.differentialExpressionVsList}
                           defaultSortFieldId={3}
                           defaultSortAsc={false}
                           customStyles={customStyles} fixedHeader
                           fixedHeaderScrollHeight="100%"
                           noDataComponent={<Text>No data.</Text>}/>
            }

        </Stack>
    );
};

export {DEGTable};