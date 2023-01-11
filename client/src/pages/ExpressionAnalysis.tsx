import React, {useMemo, useState} from 'react';
import ExpressionAnalysisTypeSelectBox
    from "../components/ExpressionAnalysisTypeSelectBox/ExpressionAnalysisTypeSelectBox";
import {Group, Loader, Stack, Text, Title} from "@mantine/core";
import {AddGene, AnnotationGroupSelectBox, LeftSidePanel, RightSidePanel} from "../components";
import UserGene from "../components/UserGene/UserGene";
import {useRecoilState, useRecoilValue} from "recoil";
import {
    annotationGroupIdState,
    selectedGenesState,
    studyIdState,
    studyLayerIdState,
    studyState,
    userGenesState
} from "../atoms";
import ProjectionPlot from "../components/ProjectionPlot/ProjectionPlot";
import {useExpressionValues} from "../hooks";
import {useExpressionViolinPlotQuery} from "../generated/types";

const analysisTypes = [
    {value: 'violinplot', label: 'Violinplot'},
    {value: 'boxplot', label: 'Boxplot'},
    {value: 'projection', label: 'Projectionplot'},
    {value: 'dot', label: 'Dotplot'},
]
const genes = [
    {display_symbol: 'CDK2'},
    {display_symbol: 'KRAS'},
    {display_symbol: 'BRD4'},
    {display_symbol: 'KLK3'},
    {display_symbol: 'ATAD2'},

]


function ViolinPlot({omicsId}: { omicsId: number }) {
    const studyId = useRecoilValue(studyIdState);
    const studyLayerId = useRecoilValue(studyLayerIdState);
    const annotationGroupId = useRecoilValue(annotationGroupIdState);
    const {data, loading} = useExpressionViolinPlotQuery({
        variables: {
            studyId,
            studyLayerId,
            omicsId,
            annotationGroupId: annotationGroupId || -1
        },
        skip: !annotationGroupId || !studyId
    })

    if (data?.violinPlot) {
        return <img src={data.violinPlot}/>;
    }
    return <div>{loading && <Loader size={25}/>}</div>
}

function ViolinPlots() {
    const selectedGenes = useRecoilValue(selectedGenesState);

    return <Stack align={'center'}>
        {selectedGenes.map((g, i) => <div key={g.omicsId}>
            <Title order={3}>{g.displaySymbol}</Title>
            <ViolinPlot omicsId={g.omicsId}/>
        </div>)}
    </Stack>;
}


function ProjectionPlots() {
    const selectedGenes = useRecoilValue(selectedGenesState);
    const {table, loading} = useExpressionValues(selectedGenes.map(g => g.omicsId));
    const tablePerGene = useMemo(() => {
        if (selectedGenes.length === 0 || !table) {
            return undefined;
        }
        return selectedGenes.map(g =>
            table.params({omicsId: g.omicsId}).filter((d: any, p: any) => d.omicsId === p.omicsId));
    }, [selectedGenes, table]);

    if (loading) {
        return <div><Loader size={25}/></div>;
    }

    return <Stack align={'center'}>
        {tablePerGene && selectedGenes.map((g, i) => <div key={g.omicsId}>
            <Title order={3}>{g.displaySymbol}</Title>
            <ProjectionPlot colorBy={'expression'} expressionTable={tablePerGene[i]}/>
        </div>)}
    </Stack>;
}


const ExpressionAnalysis = () => {
    const [analysisType, setAnalysisType] = useState<string>(analysisTypes[0].value);
    const [annotationGroupId, setAnnotationGroupId] = useRecoilState(annotationGroupIdState);

    // const [selectedAnnotationGroup, setSelectedAnnotationGroup] = useState<number>();
    const userGenes = useRecoilValue(userGenesState);
    const selectedGenes = useRecoilValue(selectedGenesState);

    const study = useRecoilValue(studyState);
    if (!study) {
        return <></>;
    }

    return (
        <Group align={'flex-start'} position={'apart'} spacing={'xs'}>
            <LeftSidePanel>
                <Stack>
                    <ExpressionAnalysisTypeSelectBox handleSelection={setAnalysisType} selection={analysisType}
                                                     options={analysisTypes}/>

                    <Stack>
                        <AnnotationGroupSelectBox changeHandler={(value: number) => {
                            setAnnotationGroupId(value);
                        }}/>
                    </Stack>
                </Stack>

            </LeftSidePanel>
            <main style={{height: '100vh', overflowY: 'scroll', flexGrow: 1}}
                  className={'plotContainer'}>
                {analysisType === 'violinplot' && <ViolinPlots/>}
                {analysisType === 'projection' && <ProjectionPlots/>}
            </main>
            <RightSidePanel>
                <Stack align={'flex-start'} justify={'flex-start'} spacing={'md'}>
                    <AddGene/>
                    <Stack spacing={'xs'}>
                        {userGenes.length > 0 ? userGenes.map((gene) => <UserGene key={`ug_${gene.displaySymbol}`}
                                                                                  gene={gene}/>) :
                            <Text color={'gray'} size={'xs'}>Nothing added yet.</Text>}
                    </Stack>
                </Stack>
            </RightSidePanel>
        </Group>

    );
};

export default ExpressionAnalysis;