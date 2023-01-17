import React, {useEffect, useMemo, useState} from 'react';
import {GeneSearchBar, MarkerCard, NavBar, SearchBar, StudyCard} from "../components";
import {Container, Grid, Space, Text} from "@mantine/core";
import {
    useExpressionByCelltypeQuery,
    useStudiesQuery
} from "../generated/types";
import {VegaLite, VisualizationSpec} from "react-vega";
import {useRecoilValue} from "recoil";
import {allGenesState} from "../atoms";


const createSpec = (multipleGenes: boolean): VisualizationSpec => {
    if (multipleGenes) {
        return {
            "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
            "data": {"name": "table"},
            "facet": {
                "column": {
                    // TODO studyName is not readable due to overlaps, needs to be rotated 45 degrees or so
                    "field": "studyName", "type": "ordinal", "header": {
                        "title": "Study"
                    }
                }
            },
            "spec": {
                "mark": {"type": "point", "filled": true},
                "encoding": {
                    "y": {
                        "field": "celltype",
                        "type": "ordinal",
                        //"sort": heatMapData.getMultiStudyExpressionHeatMap.sorting,
                        "title": "Cell Type"
                    },
                    "x": {"field": "geneSymbol", "type": "ordinal", "title": "Gene"},
                    "size": {
                        "field": "exprCellsFraction",
                        "type": "quantitative",
                        "title": "Expr. fraction"
                    },
                    "color": {
                        "field": "q3",
                        "type": "quantitative",
                        "scale": {"scheme": "viridis", reverse: true}
                    }
                }
            }
        };
    } else {
        return {
            "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
            "data": {"name": "table"},
            "mark": {"type": "point", "filled": true},
            "encoding": {
                "y": {
                    "field": "celltype",
                    "type": "ordinal",
                    //"sort": heatMapData.getMultiStudyExpressionHeatMap.sorting,
                    "title": "Cell Type"
                },
                "x": {
                    "field": "studyName", "type": "nominal", "title": "Study",
                },
                "size": {
                    "field": "exprCellsFraction",
                    "type": "quantitative",
                    "title": "Expr. fraction"
                },
                "color": {
                    "field": "q3",
                    "type": "quantitative",
                    "scale": {"scheme": "viridis", reverse: true}
                }

            }
        };
    }
};

const CrossStudySearch = () => {
    const {data: studyData, error, loading: studyDataLoading} = useStudiesQuery();
    const [omicsIds, setOmicsIds] = useState<number[]>([]);
    const allGenes = useRecoilValue(allGenesState);
    const {data, loading} = useExpressionByCelltypeQuery({
        variables: {
            omicsIds: omicsIds
        },
        skip: omicsIds.length === 0
    });

    const heatmapDisplayData = useMemo(() => {
        if (!studyData?.studyOverviewsList || !data?.expressionByCelltypesList) {
            return undefined;
        }
        const studyIdMap = new Map(studyData.studyOverviewsList.map(s => [s.studyId, s.studyName]));
        return data.expressionByCelltypesList.map(o => ({
            ...o,
            studyName: studyIdMap.get(o.studyId),
            geneSymbol: allGenes?.get(o.omicsId) || String(o.omicsId)
        }));

    }, [studyData, data]);

    return (
        <Container fluid={true}>
            <NavBar/>
            <Space h="xl"/>
            <Container size={'xl'} style={{paddingBottom: '2rem'}}>
                <GeneSearchBar humanOnly={true} onGeneSelection={ids => setOmicsIds(ids)}/>
            </Container>
            <Container size={'xl'}>
                {heatmapDisplayData && <VegaLite
                    spec={createSpec(omicsIds.length > 1)} //heatmapSelectedGenes.length > 1, heatmapChosenStudies.length)}
                    //onNewView={(view) => setUpSelectionListener(view)}
                    data={{
                        "table": heatmapDisplayData
                    }}/>}
            </Container>

        </Container>
    );
};

export default CrossStudySearch;