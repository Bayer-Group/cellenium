import React from 'react';
import Plot from "react-plotly.js";
import {useHalfAVolcanoQuery} from "../../generated/types";
import {Loader} from "@mantine/core";


interface Props {
    studyId: number;
    annotationValueId: number;
    pval: number;
    log2fc: number;
};

const InlineFoldChangePlot = ({studyId, annotationValueId, pval, log2fc}: Props) => {
    const {data, error, loading} = useHalfAVolcanoQuery({
        variables: {
            studyId: studyId,
            annotationValueId: annotationValueId
        }
    })
    if (loading)
        return <Loader variant={'dots'}/>
    return (<>
            {
                <Plot
                    config={{
                        modeBarButtons: false,
                        displaylogo: false,
                        staticPlot: true
                    }}
                    data={[
                        {
                            x: data?.differentialExpressionsList.map((e) => e.log2Foldchange),
                            y: data?.differentialExpressionsList.map((e) => {
                                let ppval = -1 * Math.log10(e.pvalueAdj)
                                return ppval?ppval:0
                            }),
                            type: 'scatter',
                            mode: 'markers',
                            marker: {color: 'lightblue', size: 3},
                            showlegend: false,
                        }, {
                            x: [log2fc],
                            y: [-1 * Math.log10(pval)],
                            marker: {
                                color: 'red',
                                size: 5
                            },
                            showlegend: false
                        }]}

                    layout={{
                        width: 160, height: 100, margin: {l: 26, r: 5, t: 5, b: 20},

                        xaxis: {
                            range: [0, 10],
                            visible: true,
                            fixedrange: true,
                            showgrid: false,
                            title: 'log2FC',
                            titlefont: {
                                size: 10,
                                color: 'grey'
                            },
                            tickfont: {
                                size: 6
                            }
                        },
                        yaxis: {
                            visible: true,
                            fixedrange: true,
                            showgrid: false,
                            title: '-log(padj)',
                            titlefont: {
                                family: 'Arial, sans-serif',
                                size: 10,
                                color: 'grey'
                            },
                            tickfont: {
                                size: 6
                            }
                        },
                    }}

                />
            }
        </>
    );
};

export {InlineFoldChangePlot};