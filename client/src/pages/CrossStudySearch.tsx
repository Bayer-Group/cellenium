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
            filter: {
                omicsId: {in: omicsIds},
                annotationGroupId: {equalTo: cellOAnnotationGroupId || -1}
            } as ExpressionByAnnotationFilter
        },
        skip: omicsIds.length === 0 || !studyData?.studyOverviewsList || !cellOAnnotationGroupId
    })

    const heatmapDisplayData = useMemo(() => {
        if (!studyData?.studyOverviewsList || !data?.expressionByAnnotationsList) {
            return undefined;
        }
        const studyIdMap = new Map(studyData.studyOverviewsList.map(s => [s.studyId, s.studyName]));
        const allData = data.expressionByAnnotationsList.map(o => ({
            ...o,
            studyName: studyIdMap.get(o.studyId)
        }));
        return omicsIds.map(id => ({
            omicsId: id,
            heatmapData: allData.filter(d => d.omicsId === id)
        }));
    }, [studyData, data, omicsIds]);

    const onHeatmapClick = (dotPlotElement: DotPlotElementFragment, event: ScenegraphEvent) => {
        const newStudyUrl = `/study/${dotPlotElement.studyId}?page=CellMarkerAnalysis&annotationGroupId=${dotPlotElement.annotationGroupId}&annotationValueId=${dotPlotElement.annotationValueId}&omicsId=${dotPlotElement.omicsId}`;
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