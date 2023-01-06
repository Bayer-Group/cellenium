import React from 'react';
import DataTable from "react-data-table-component";
import {Space, Stack, Text} from "@mantine/core";
import {IconPlus} from "@tabler/icons";
import { ActionIcon } from '@mantine/core';


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
const columns = [
    {
        name: 'gene',
        selector: (row: any) => row.symbol,
        sortable: true,
        width: '80px'
    },
    {
        name: 'padj',
        selector: (row: any) => row.padj.toExponential(2),
        sortable: true,
        width: '80px'
    },
    {
        name: 'log2FC',
        selector: (row: any) => +row.log2fc.toFixed(2),
        sortable: true,
        width: '80px'
    },
    {
        name: '',
        cell: () => <ActionIcon size='xs' variant={"default"}><IconPlus size={12}/></ActionIcon>,
        width: '20px',
    }
];

type Props = {
    data: object[];
}

const DEGTable = ({data}: Props) => {
    return (
        <Stack justify={'flex-start'} align={'center'} w={'100%'}>
            <div>
                <Text weight={800} size={'xs'}>Differentially expressed genes</Text>
                <Space h={'xs'}/>
                <DataTable dense columns={columns} data={data} defaultSortFieldId={1}
                           customStyles={customStyles}/>
            </div>
        </Stack>
    );
};

export {DEGTable};