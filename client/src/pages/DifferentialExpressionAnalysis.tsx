import React from 'react';
import {
    AddGene,
    AnnotationGroup,
    AnnotationGroupSelectBox,
    DEGTable,
    LeftSidePanel,
    RightSidePanel
} from "../components";
import {Group, Space, Stack} from "@mantine/core";
import UserGene from "../components/UserGene/UserGene";

const ANNOTATIONS = [
    {label: "bone cell", color: "#1f77b4"},
    {label: "chondrocyte", color: "#ff7f0e"},
    {label: "endothelial cell", color: "#2ca02c"},
    {label: "endothelial cell of artery", color: "#d62728"},
    {label: "fibroblast", color: "#9467bd"},
    {label: "mesenchymal stem cell", color: "#8c564b"},
    {label: "pericyte cell", color: "#e377c2"}
]

const DEG: object[] = [
    {
        symbol: 'BRD4',
        padj: 0.001342,
        log2fc: 2132.23
    },
    {
        symbol: 'PTK2',
        padj: 0.001342,
        log2fc: 2.23
    },
    {
        symbol: 'CDK2',
        padj: 0.001342,
        log2fc: 24
    },
    {
        symbol: 'EGFR',
        padj: 0.001342,
        log2fc: 1
    },
    {
        symbol: 'KRAS',
        padj: 0.001342,
        log2fc: 3.324
    },
    {
        symbol: 'KLK3',
        padj: 0.00000000000000001342,
        log2fc: 432.3231322
    },
    {
        symbol: 'PLK1',
        padj: 0.001342,
        log2fc: 2132.23
    },
]

const genes = [
    {display_symbol: 'CDK2'},
    {display_symbol: 'KRAS'},
    {display_symbol: 'BRD4'},
    {display_symbol: 'KLK3'},
    {display_symbol: 'ATAD2'},

]

function DifferentialExpressionAnalysis() {
    return (
        <Group position={'apart'}>
            <LeftSidePanel>
                <Stack>
                    <AnnotationGroupSelectBox/>
                    <AnnotationGroup annotations={ANNOTATIONS}/>
                </Stack>
            </LeftSidePanel>

            <RightSidePanel>
                <Stack align={'flex-start'} justify={'flex-start'} spacing={'md'}>
                    <div style={{width: '80%'}}>
                        <AddGene/>
                    </div>
                    <Stack style={{width: '100%'}} spacing={'xs'}>
                        {genes.map((gene) => <UserGene gene={gene}/>)}
                    </Stack>
                    <Space h={'md'}/>
                    <DEGTable data={DEG}/>
                </Stack>
            </RightSidePanel>
        </Group>
    );
};

export default DifferentialExpressionAnalysis;