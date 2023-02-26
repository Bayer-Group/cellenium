import React from 'react';
import DataTable from "react-data-table-component";
import {ActionIcon, Center, Group, Loader, Stack, Text} from "@mantine/core";
import {IconChevronDown, IconChevronRight, IconEye, IconPlus} from "@tabler/icons";
import {useDegQuery} from "../../generated/types";
import memoize from 'memoize-one';
import {useRecoilState, useRecoilValue} from "recoil";
import {
    annotationGroupIdState,
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
        name: 'identifier',
        selector: (row: any) => <Text title={row.displaySymbol}>{row.displaySymbol}</Text>,
        sortable: true,
        width: '80px'
    },
    {
        name: 'padj',
        selector: (row: any) => row.pvalueAdj.toFixed(2),
        sortable: true,
        width: '50px'
    },
    {
        name: 'lgFC',
        selector: (row: any) => {
            let ret;
            if (row.log2Foldchange !== -1)
                ret = +row.log2Foldchange.toFixed(2);
            else
                ret = '';
            return ret
        },
        sortable: true,
        width: '50px'
    },
    {
        name: '',
        cell: (row: any) => {
            let gene = {
                omicsId: row.omicsId,
                displayName: row.displayName,
                displaySymbol: row.displaySymbol,
                omicsType: row.omicsType,
                value: row.displaySymbol,
                linkedGenes: row.linkedGenes
            }

            return (
                <Group position={'center'} align={'center'} noWrap={true} spacing={0}>
                    <ActionIcon title={'superpose expression'}
                                onClick={() => handleColorClick(gene)} variant="default" size={'xs'} mr={5}>
                        <IconEye/>
                    </ActionIcon>
                    <ActionIcon title={'add to gene store'} color={'blue.3'} onClick={() => clickHandler(gene)}
                                size='xs'
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
const ExpandedComponent = ({data}: any) => {
    const study = useRecoilValue(studyState);
    const [userGenes, setUserGenes] = useRecoilState(userGenesState);
    const [indicatorColor, setIndicatorColor] = useRecoilState(userGeneStoreCounterColor);
    const [storeOpen, setStoreOpen] = useRecoilState(userGeneStoreOpenState)
    const [selectedGenes, setSelectedGenesStore] = useRecoilState(selectedGenesState);

    function showExpression(gene: Omics) {
        if (selectedGenes.filter((g) => g.omicsId === gene.omicsId).length > 0) {
            // remove
            let removed = selectedGenes.filter((g) => g.omicsId !== gene.omicsId)
            setSelectedGenesStore(removed)
        } else {
            setSelectedGenesStore([gene]);
        }

    }

    function addToStore(gene: Omics) {
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
                message: "",
                color: 'red',
                autoClose: 1000
            })
        }
    }

    const linkedGenes: any[] = data.linkedGenes.map((id: number) => {
        let gene = study?.studyOmicsMap.get(id);
        if (gene !== undefined)
            return (
                <Group spacing={'xs'}>
                    <Text size={'xs'}>{gene['displaySymbol']}</Text>
                    <ActionIcon title={'add to gene store'} color={'blue.3'}
                                onClick={() => showExpression(gene as Omics)}
                                size='xs'
                                variant={"default"}><IconEye
                        size={12} color={'black'}/>
                    </ActionIcon>
                    <ActionIcon title={'add to gene store'} color={'blue.3'} onClick={() => addToStore(gene as Omics)}
                                size='xs'
                                variant={"default"}><IconPlus
                        size={12} color={'black'}/></ActionIcon>
                </Group>)
        else
            return undefined
    });
    return (
        <pre>
            <Center>
                <Stack>
            <Text weight={800} size={'xs'}>Corresponding gene(s)</Text>
                {linkedGenes && linkedGenes.length > 0 &&
                    linkedGenes
                }
                </Stack>
        </Center>
        </pre>
    )
}


const DEGTable = ({annotationId}: Props) => {
    const [userGenes, setUserGenes] = useRecoilState(userGenesState);
    const [indicatorColor, setIndicatorColor] = useRecoilState(userGeneStoreCounterColor);
    const [selectedGenes, setSelectedGenesStore] = useRecoilState(selectedGenesState);
    const annotationGroup = useRecoilValue(annotationGroupIdState);
    const study = useRecoilValue(studyState);
    const [storeOpen, setStoreOpen] = useRecoilState(userGeneStoreOpenState)

    let selectedAGM = study?.annotationGroupMap;
    let modality: string | undefined = 'rna';
    if (annotationGroup && selectedAGM && selectedAGM.get(annotationGroup) !== undefined) {
        modality = selectedAGM.get(annotationGroup)?.modality;
    }
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
                message: "",
                color: 'red',
                autoClose: 1000
            })
        }
    }

    return (
        <Stack justify={'flex-start'} align={'center'} w={'100%'}>
            {data && data.differentialExpressionVsList.length > 0 && (modality === 'atac' || modality === "prot") &&
                <DataTable dense columns={columns(handleClick, handleColorClick)}
                           data={data.differentialExpressionVsList}
                           defaultSortFieldId={3}
                           defaultSortAsc={false}
                           customStyles={customStyles} fixedHeader
                           fixedHeaderScrollHeight="100%"
                           noDataComponent={<Text>No data.</Text>}
                           expandableIcon={{
                               collapsed: <IconChevronRight size={10}/>,
                               expanded: <IconChevronDown size={10}/>
                           }}
                           expandableRows
                           expandableRowsComponent={ExpandedComponent}
                />
            }
            {data && data.differentialExpressionVsList.length > 0 && (modality === 'rna') &&
                <DataTable dense columns={columns(handleClick, handleColorClick)}
                           data={data.differentialExpressionVsList}
                           defaultSortFieldId={3}
                           defaultSortAsc={false}
                           customStyles={customStyles} fixedHeader
                           fixedHeaderScrollHeight="100%"
                           noDataComponent={<Text>No data.</Text>}
                />
            }
            {loading && <Loader variant={'dots'} color={'gray'}/>}
        </Stack>
    );
};

export {DEGTable};