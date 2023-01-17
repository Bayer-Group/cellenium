import React, {useMemo, useState} from 'react';
import ExpressionAnalysisTypeSelectBox
    from "../components/ExpressionAnalysisTypeSelectBox/ExpressionAnalysisTypeSelectBox";
import {Center, Divider, Group, Loader, Stack, Text, Title, useMantineTheme} from "@mantine/core";
import {
    AnnotationFilterDisplay,
    AnnotationGroupSelectBox,
    LeftSidePanel,
    RightSidePanel,
    UserGeneStore
} from "../components";
import {useRecoilState, useRecoilValue} from "recoil";
import {
    annotationGroupIdState,
    selectedAnnotationFilterState,
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
    {value: 'projection', label: 'Projectionplot'}
    /*
    {value: 'boxplot', label: 'Boxplot'},
        {value: 'dot', label: 'Dotplot'},

     */
]

function ViolinPlot({omicsId}: { omicsId: number }) {
    const theme = useMantineTheme();
    const studyId = useRecoilValue(studyIdState);
    const studyLayerId = useRecoilValue(studyLayerIdState);
    const annotationGroupId = useRecoilValue(annotationGroupIdState);
    const annotationFilter = useRecoilValue(selectedAnnotationFilterState);
    const {data, loading} = useExpressionViolinPlotQuery({
        variables: {
            studyId,
            studyLayerId,
            omicsId,
            annotationGroupId: annotationGroupId || -1,
            excludeAnnotationValueIds: annotationFilter
        },
        skip: !annotationGroupId || !studyId
    })

    if (data?.violinPlot) {
        return <img src={data.violinPlot}/>;
    }
    return <div>{loading && <Loader variant={'dots'} color={theme.colors.gray[5]} size={'xl'}/>}</div>
}

function ViolinPlots() {
    const selectedGenes = useRecoilValue(selectedGenesState);

    return <Stack align={'center'}>
        {[...selectedGenes].reverse().map((g, i) => <Stack key={g.omicsId} align={'center'}>
            <Title order={3}>{g.displaySymbol}</Title>
            <ViolinPlot omicsId={g.omicsId}/>
        </Stack>)}
    </Stack>;
}


function ProjectionPlots() {
    const theme = useMantineTheme();
    const selectedGenes = useRecoilValue(selectedGenesState);
    const {table, loading} = useExpressionValues(selectedGenes.map(g => g.omicsId));
    const tablePerGene = useMemo(() => {
        if (selectedGenes.length === 0 || !table) {
            return undefined;
        }
        return [...selectedGenes].reverse().map(g =>
            table.params({omicsId: g.omicsId}).filter((d: any, p: any) => d.omicsId === p.omicsId));
    }, [selectedGenes, table]);

    if (loading) {
        return <Center style={{height: '100%', width: '100%'}}><Loader variant={'dots'} color={theme.colors.gray[5]}
                                                                       size={25}/></Center>;
    }

    return <Stack>
        {tablePerGene && [...selectedGenes].reverse().map((g, i) => <Stack key={g.omicsId} align={'center'}>
            <Title>{g.displaySymbol}</Title>
            <ProjectionPlot colorBy={'expression'} expressionTable={tablePerGene[i]}/>
        </Stack>)}
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
    const showAnnotationSelectors = (
        <div>
            <AnnotationGroupSelectBox changeHandler={(value: number) => {
                setAnnotationGroupId(value);
            }}/>
            <Divider my="sm"/>
            <AnnotationFilterDisplay/>
        </div>
    );
    return (
        <Group align={'flex-start'} position={'apart'} spacing={'xs'}>
            <LeftSidePanel>
                <Stack>
                    <ExpressionAnalysisTypeSelectBox handleSelection={setAnalysisType} selection={analysisType}
                                                     options={analysisTypes}/>
                    {analysisType !== 'projection' && showAnnotationSelectors}
                </Stack>

            </LeftSidePanel>
            <main style={{height: '100vh', overflowY: 'scroll', flexGrow: 1, paddingTop: 60}}
                  className={'plotContainer'}>
                <Stack align={'center'} justify={'center'}>
                    {analysisType === 'violinplot' && <ViolinPlots/>}
                    {analysisType === 'projection' && <ProjectionPlots/>}
                    {selectedGenes.length === 0 &&
                        <Text c={'dimmed'}>Please select gene(s) from the <Text span weight={800}>User gene
                            store</Text></Text>}
                </Stack>
            </main>
            <RightSidePanel>
                <Stack>
                    <Divider size={"xs"} label={'User gene store'}/>
                    <UserGeneStore multiple={true} opened={true}/>
                </Stack>
            </RightSidePanel>
        </Group>

    );
};

export default ExpressionAnalysis;