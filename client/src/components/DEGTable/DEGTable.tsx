import React from 'react';
import DataTable from "react-data-table-component";
import {ActionIcon, Stack, Text} from "@mantine/core";
import {IconPlus} from "@tabler/icons";
import {useDegQuery} from "../../generated/types";
import memoize from 'memoize-one';
import {useRecoilState, useRecoilValue} from "recoil";
import {studyState, userGenesState, userGeneStoreCounterColor, userGeneStoreOpenState} from "../../atoms";
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
const columns = memoize((clickHandler) => [
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
                <ActionIcon color={'blue.3'} onClick={() => clickHandler(gene)} size='xs' variant={"default"}><IconPlus
                    size={12} color={'black'}/></ActionIcon>)
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
    const study = useRecoilValue(studyState);
    const [storeOpen, setStoreOpen] = useRecoilState(userGeneStoreOpenState)
    const {data, error, loading} = useDegQuery({
        variables: {
            annotationValueId: annotationId,
            studyId: study?.studyId || 0
        }
    })

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
                <DataTable dense columns={columns(handleClick)} data={data.differentialExpressionVsList}
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