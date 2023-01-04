import React from 'react';
import DataTable from "react-data-table-component";
import {Stack} from "@mantine/core";

const columns = [
    {
        name: 'gene',
        selector: (row: any) => row.symbol,
        sortable: true,
        width: '4rem'
    },
    {
        name: 'padj',
        selector: (row: any) => row.padj.toExponential(2),
        sortable: true,
        width: '4rem'
    },
    {
        name: 'log2FC',
        selector: (row: any) =>+row.log2fc.toFixed(2),
        sortable: true,
        width: '6rem'
    },
];

const customStyles = {
    table: {
        style: {
            backgroundColor: 'transparent',
            marginRight: '10',
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

type Props = {
    data: object[];
}

const DEGTable = ({data}: Props) => {
    return (
        <Stack justify={'flex-start'} align={'center'}>
            <div>
                <DataTable dense title="DEG" columns={columns} data={data} defaultSortFieldId={1}
                           customStyles={customStyles}/>
            </div>
        </Stack>
    );
};

export {DEGTable};