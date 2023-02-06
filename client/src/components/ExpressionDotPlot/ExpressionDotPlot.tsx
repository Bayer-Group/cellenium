import React, {useEffect, useMemo, useState} from 'react';
import {VegaLite, View, VisualizationSpec} from "react-vega";
import {DotPlotElementFragment} from "../../generated/types";
import {ScenegraphEvent} from "vega";


function createSpec(annotationTitle: string, xAxis: "studyName" | "displaySymbol") {
    return {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "data": {"name": "table"},
        "mark": {"type": "point", "filled": true},
        "encoding": {
            "y": {
                "field": "annotationDisplayValue",
                "type": "ordinal",
                "title": annotationTitle
            },
            "x": xAxis === "studyName" ? {
                "field": "studyName", "type": "nominal", "title": "Study",
            } : {
                "field": "displaySymbol", "type": "nominal", "title": "Gene",
            },
            "size": {
                "field": "exprSamplesFraction",
                "type": "quantitative",
                "title": "Expr. fraction"
            },
            "color": {
                "field": "q3",
                "type": "quantitative",
                "scale": {"scheme": "viridis", reverse: true}
            }

        }
    } as VisualizationSpec;
}

export function ExpressionDotPlot({
                                      data,
                                      annotationTitle,
                                      xAxis,
                                      onClick
                                  }:
                                      {
                                          data: DotPlotElementFragment[],
                                          annotationTitle: string,
                                          xAxis: "studyName" | "displaySymbol",
                                          onClick?: (dotPlotElement: DotPlotElementFragment, event: ScenegraphEvent) => void
                                      }) {
    const spec = useMemo(() => createSpec(annotationTitle, xAxis), [annotationTitle, xAxis]);

    const setUpSelectionListener = (view: View) => {
        view.addEventListener("click", (event, item) => {
            if (item && onClick) {
                const dotPlotElement = item.datum as DotPlotElementFragment;
                onClick(dotPlotElement, event);
            }
        });
    };

    return <VegaLite
        spec={spec}
        onNewView={(view) => setUpSelectionListener(view)}
        data={{
            "table": data
        }}/>
}

