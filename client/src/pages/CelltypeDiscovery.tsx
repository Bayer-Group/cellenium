import React, {useEffect, useMemo, useState} from 'react';
import {Box, Button, Group, Loader, Overlay, Stack, Text, TextInput} from "@mantine/core";
import {AnnotationGroupDisplay, AnnotationGroupSelectBox, LeftSidePanel, RightSidePanel} from "../components";
import {useRecoilState, useRecoilValue, useSetRecoilState} from "recoil";
import {
    annotationGroupIdState,
    celltypeDiscoveryCoexpressionSamplesState,
    celltypeDiscoveryGenesState,
    selectedAnnotationState,
    studyReloadHelperState,
    studyState,
    userGenesState
} from "../atoms";
import {useExpressionValues} from "../hooks";
import * as aq from 'arquero';
import Plot from "react-plotly.js";
import {SingleGeneSelection} from "../components/AddGene/SingleGeneSelection";
import {Omics} from "../model";
import ProjectionPlot from "../components/ProjectionPlot/ProjectionPlot";
import * as Plotly from "plotly.js";
import {InputMaybe, useSaveUserAnnotationMutation} from "../generated/types";
import {showNotification} from "@mantine/notifications";

interface PreparedPlot {
    message?: string;
    plotlyData: Partial<Plotly.PlotData>[];
    plotlyLayout: Partial<Plotly.Layout>;
}

const plotlyConfig: Partial<Plotly.Config> = {
    responsive: true
};

const UNEXPRESSED_SAMPLE_ID = -100;

function CoexpressionPlot({
                              stateOffset,
                              onSelection,
                              onDoubleClick
                          }: {
    stateOffset: number;
    onSelection: (event: Readonly<Plotly.PlotSelectionEvent>) => void;
    onDoubleClick: () => void;
}) {
    const [omicsAll, setOmicsAll] = useRecoilState(celltypeDiscoveryGenesState);
    const omicsX = omicsAll[stateOffset * 2];
    const omicsY = omicsAll[stateOffset * 2 + 1];
    const setOmicsX = (omics: Omics | null) => setOmicsAll(prev => prev.map((o, i) => i === stateOffset * 2 ? omics : o));
    const setOmicsY = (omics: Omics | null) => setOmicsAll(prev => prev.map((o, i) => i === stateOffset * 2 + 1 ? omics : o));
    const {table, loading} = useExpressionValues((omicsX && omicsY) ? [omicsX.omicsId, omicsY.omicsId] : [], false);

    const celltypeDiscoveryCoexpressionSamples = useRecoilValue(celltypeDiscoveryCoexpressionSamplesState);
    const filterSampleIds = celltypeDiscoveryCoexpressionSamples[stateOffset];

    const preparedPlot = useMemo(() => {
        if (!table || !omicsX || !omicsY) {
            return undefined;
        }

        const subTableA = table
            .params({plotFilter: omicsX.omicsId})
            .filter((d: any, p: any) => d.omicsId === p.plotFilter)
            .select({studySampleId: 'studySampleId', value: 'valueA'});
        const subTableB = table
            .params({plotFilter: omicsY.omicsId})
            .filter((d: any, p: any) => d.omicsId === p.plotFilter)
            .select({studySampleId: 'studySampleId', value: 'valueB'});
        let sameSampleExprValues = subTableA.join_full(subTableB, ['studySampleId', 'studySampleId']).derive({index: () => aq.op.row_number() - 1})
            .impute({valueA: () => 0, valueB: () => 0});
        if (filterSampleIds) {
            // @ts-ignore
            sameSampleExprValues = sameSampleExprValues.params({filterSampleIds}).filter((d, p) => aq.op.includes(p.filterSampleIds, d.studySampleId, 0));
        }

        // to show a dummy sample with 0 expression in both genes (not retrievable from the sparse DB data)
        sameSampleExprValues = sameSampleExprValues.concat(aq.from([{
            valueA: 0,
            valueB: 0,
            studySampleId: UNEXPRESSED_SAMPLE_ID
        }]));


        const plot: Partial<Plotly.PlotData> = {
            type: 'scattergl',
            x: sameSampleExprValues.array('valueA', Float32Array),
            y: sameSampleExprValues.array('valueB', Float32Array),
            customdata: sameSampleExprValues.array('studySampleId', Int32Array),
            mode: 'markers',
            marker: {
                //size: markerSize,
                color: '#415370',
                opacity: 0.6,
            },
            hoverinfo: 'none',
            showlegend: false,
        } as Partial<Plotly.PlotData>;

        return {
            //allSameSampleExprValues,
            plotlyData: [plot],
            //combinations,
            plotlyLayout: {
                //autosize: true,
                width: 250,
                height: 250,
                margin: {l: 40, r: 0, t: 0, b: 40},
                dragmode: 'lasso',
                clickmode: 'none',
            } as Partial<Plotly.Layout>,
        } as PreparedPlot;
    }, [table]);

    return (<div style={{width: '100%'}}>
        <Group>
            <Group><Text>X</Text><SingleGeneSelection selection={omicsX} onSelectionChange={setOmicsX}/></Group>
            <Group><Text>Y</Text><SingleGeneSelection selection={omicsY} onSelectionChange={setOmicsY}/></Group>
        </Group>
        {preparedPlot && !preparedPlot.message && <Plot data={preparedPlot.plotlyData}
                                                        layout={preparedPlot.plotlyLayout}
                                                        config={plotlyConfig}
                                                        onSelected={onSelection}
                                                        onDeselect={onDoubleClick}
        />}
        {!preparedPlot && <div style={{height: 250, width: 250}}>
            {loading && <Loader variant={'dots'}/>}
        </div>}
    </div>);
}

function CelltypeDiscovery() {
    const study = useRecoilValue(studyState);
    const annotationGroupId = useRecoilValue(annotationGroupIdState);
    const [omicsAll, setOmicsAll] = useRecoilState(celltypeDiscoveryGenesState);
    const [celltypeDiscoveryCoexpressionSamples, setCelltypeDiscoveryCoexpressionSamples] = useRecoilState(celltypeDiscoveryCoexpressionSamplesState);
    const [selected, setSelected] = useRecoilState(selectedAnnotationState);

    useEffect(() => {
        setSelected(0)
    }, [])

    // convenient default for gene input
    const userGenes = useRecoilValue(userGenesState);
    useEffect(() => {
        if (!omicsAll[0] && !omicsAll[1]) {
            if (userGenes.length > 1) {
                setOmicsAll([userGenes[0], userGenes[1]]);
            } else if (userGenes.length > 0) {
                setOmicsAll([userGenes[0], null]);
            }
        }
    }, [omicsAll, userGenes]);

    const [selectedSampleIds, setSelectedSampleIds] = useState<number[] | null>(null);

    const onCoexpressionSelection = (event: Readonly<Plotly.PlotSelectionEvent>) => {
        const theelectedSampleIds = event.points.map(p => p.customdata) as number[];
        setSelectedSampleIds(theelectedSampleIds);
    };
    const onCoexpressionDoubleClick = () => {
        setSelectedSampleIds([]);
    };

    const newPlotBasedOnSelectedSamples = () => {
        setOmicsAll(prev => [...prev, prev[prev.length - 2], prev[prev.length - 1]]);
        setCelltypeDiscoveryCoexpressionSamples(prev => [...prev, selectedSampleIds]);
        setSelectedSampleIds(null);
    };
    const removeLastPlot = () => {
        setOmicsAll(prev => prev.slice(0, prev.length - 2));
        setCelltypeDiscoveryCoexpressionSamples(prev => prev.slice(0, prev.length - 1));
    };
    const [annotationName, setAnnotationName] = useState("");
    const [saveUserAnnotationMutation, {
        data: saveUserAnnotationResult,
        loading: saveUserAnnotationLoading
    }] = useSaveUserAnnotationMutation();
    const setStudyReloadHelper = useSetRecoilState(studyReloadHelperState);
    const saveUserAnnotation = () => {
        if (study && selectedSampleIds && selectedSampleIds.length > 0) {

            const omicsX = omicsAll[(celltypeDiscoveryCoexpressionSamples.length - 1) * 2];
            const omicsY = omicsAll[(celltypeDiscoveryCoexpressionSamples.length - 1) * 2 + 1];
            saveUserAnnotationMutation({
                variables: {
                    studyId: study.studyId,
                    annotationGroupName: annotationName,
                    selectedSampleIds: selectedSampleIds.filter(id => id !== UNEXPRESSED_SAMPLE_ID).join(','),
                    unexpressedSamplesOmicsIds: (selectedSampleIds.indexOf(UNEXPRESSED_SAMPLE_ID) !== -1 ?
                            [omicsX?.omicsId, omicsY?.omicsId] : null
                    ) as InputMaybe<number[]>
                }
            }).then(value => {
                showNotification({
                    title: 'Successfully saved user annotation',
                    message: 'Differentially expressed genes are calculated in the background and will appear later.',
                    color: 'green',
                    autoClose: 2000
                });
                setAnnotationName("");
                setStudyReloadHelper(prev => prev + 1);
            }).catch(reason => {
                showNotification({
                    title: 'Could not save user annotation',
                    message: reason.message,
                    color: 'red'
                });
            });
        }
    };

    if (!study) {
        return <></>;
    }
    return (
        <Group position={'apart'}>
            <LeftSidePanel>
                <Stack>
                    {annotationGroupId && <AnnotationGroupSelectBox/>}
                    {annotationGroupId && <AnnotationGroupDisplay disableSelection={true}/>}
                </Stack>
            </LeftSidePanel>
            <main>
                <ProjectionPlot colorBy={'annotation'} showSampleIds={selectedSampleIds} disableSelection={true}/>
            </main>
            <RightSidePanel>
                <Stack align={'flex-start'} justify={'flex-start'} spacing={'md'}>
                    {[...Array(celltypeDiscoveryCoexpressionSamples.length)].map((z, i) => (
                        <Box sx={{height: 350, position: 'relative'}} key={i}>
                            {/*  Disable edit if there is a plot following (which depends on the selection in the current plot)  */}
                            {i < celltypeDiscoveryCoexpressionSamples.length - 1 &&
                                <Overlay opacity={0.6} color="#000" zIndex={5}/>}
                            <CoexpressionPlot
                                stateOffset={i}
                                onSelection={onCoexpressionSelection}
                                onDoubleClick={onCoexpressionDoubleClick}
                            />
                        </Box>
                    ))}
                    <Button onClick={newPlotBasedOnSelectedSamples}
                            disabled={selectedSampleIds === null || selectedSampleIds.length === 0}>Add plot, based on
                        selected samples</Button>
                    <Button onClick={removeLastPlot}
                            disabled={celltypeDiscoveryCoexpressionSamples.length < 2}>Remove last plot</Button>
                    <Group>
                        <TextInput value={annotationName}
                                   placeholder={"My Custom Cell Annotation"}
                                   onChange={(event) => setAnnotationName(event.currentTarget.value)}/>
                        <Button
                            disabled={selectedSampleIds === null || selectedSampleIds.length === 0 || annotationName.length === 0}
                            onClick={saveUserAnnotation}
                            loading={saveUserAnnotationLoading}>Save</Button>
                    </Group>
                </Stack>
            </RightSidePanel>
        </Group>
    );
};

export default CelltypeDiscovery;
