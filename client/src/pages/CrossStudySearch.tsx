import React, {useEffect, useMemo, useState} from 'react';
import {GeneSearchBar, MarkerCard, NavBar, SearchBar, StudyCard} from "../components";
import {Center, Container, Grid, Group, Loader, Space, Stack, Text, useMantineTheme} from "@mantine/core";
import {
    DotPlotElementFragment, ExpressionByAnnotationFilter,
    useExpressionByAnnotationQuery,
    useStudiesQuery
} from "../generated/types";
import {useRecoilState, useRecoilValue} from "recoil";
import {allGenesState, cellOAnnotationGroupIdState} from "../atoms";
import {useNavigate} from "react-router-dom";
import {ExpressionDotPlot} from "../components/ExpressionDotPlot/ExpressionDotPlot";
import {ScenegraphEvent} from "vega";

const CrossStudySearch = () => {
    const theme = useMantineTheme();
    const {data: studyData, error, loading: studyDataLoading} = useStudiesQuery();
    const [omicsIds, setOmicsIds] = useState<number[]>([]);
    const cellOAnnotationGroupId = useRecoilValue(cellOAnnotationGroupIdState);
    const allGenes = useRecoilValue(allGenesState);
    const navigate = useNavigate();
    const {data, loading} = useExpressionByAnnotationQuery({
        variables: {
            studyLayerIds: studyData?.studyOverviewsList.map(s => s.defaultStudyLayerId) || [],
            omicsIds,
            annotationGroupId: cellOAnnotationGroupId || -1,
            excludeAnnotationValueIds: []
        },
        skip: omicsIds.length === 0 || !studyData?.studyOverviewsList || !cellOAnnotationGroupId
    })

    const heatmapDisplayData = useMemo(() => {
        if (!studyData?.studyOverviewsList || !data?.expressionByAnnotationList) {
            return undefined;
        }
        const studyLayerIdMap = new Map(studyData.studyOverviewsList.map(s => [s.defaultStudyLayerId, s]));
        const allData = data.expressionByAnnotationList.map(o => ({
            ...o,
            studyName: studyLayerIdMap.get(o.studyLayerId)?.studyName,
            studyId: studyLayerIdMap.get(o.studyLayerId)?.studyId
        }));
        return omicsIds.map(id => ({
            omicsId: id,
            heatmapData: allData.filter(d => d.omicsId === id)
        }));
    }, [studyData, data, omicsIds]);

    const onHeatmapClick = (dotPlotElement: DotPlotElementFragment, event: ScenegraphEvent) => {
        const dpe = dotPlotElement as DotPlotElementFragment & { studyId: number };
        const newStudyUrl = `/study/${dpe.studyId}?page=CellMarkerAnalysis&annotationGroupId=${cellOAnnotationGroupId}&annotationValueId=${dotPlotElement.annotationValueId}&omicsId=${dotPlotElement.omicsId}`;
        if (event.shiftKey || event.altKey) {
            const parsedUrl = new URL(window.location.href);
            const url = `${parsedUrl.origin}${newStudyUrl}`
            window.open(url, '_blank')
        } else {
            navigate(newStudyUrl);
        }
    };

    return (
        <Container fluid={true}>
            <NavBar/>
            <Space h="xl"/>
            <Container size={'xl'} style={{paddingBottom: '2rem'}}>
                <GeneSearchBar humanOnly={true} onGeneSelection={ids => setOmicsIds(ids)}/>
            </Container>
            {omicsIds.length === 0 && <Center>
                <Text style={{'width': '50em'}} color={'dimmed'}>Please enter your genes of interest. Cellenium will
                    show the gene's
                    expression in human studies with standardized cell annotation (CellO). As the study data is
                    processed and normalized independently, this is a qualitative direction for which studies
                    to explore independently. Click in the chart to open a study.</Text>
            </Center>}
            <Group position={"center"}>
                {(loading || studyDataLoading) && <Loader variant={'dots'} color={theme.colors.gray[5]} size={25}/>}
                {heatmapDisplayData && heatmapDisplayData.map(heatmap => <Stack>
                    <Text>{allGenes?.get(heatmap.omicsId)?.displaySymbol}</Text>
                    <ExpressionDotPlot data={heatmap.heatmapData}
                                       annotationTitle={"Cell Type"}
                                       xAxis={"studyName"}
                                       onClick={onHeatmapClick}
                    />
                </Stack>)}
            </Group>
        </Container>
    );
};

export default CrossStudySearch;