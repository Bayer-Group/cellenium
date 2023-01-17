import React from 'react';
import Plot from "react-plotly.js";

const DATA = [
    {
        "log2Foldchange": 1.612797498703003,
        "pvalueAdj": 0.0033313693114854475
    },
    {
        "log2Foldchange": 1.5811443328857422,
        "pvalueAdj": 0.006197080265921125
    },
    {
        "log2Foldchange": 0.7544171214103699,
        "pvalueAdj": 0.3117244734657245
    },
    {
        "log2Foldchange": 1.50431227684021,
        "pvalueAdj": 0.052188886673494556
    },
    {
        "log2Foldchange": 1.3831335306167603,
        "pvalueAdj": 0.30767269921571144
    },
    {
        "log2Foldchange": 1.2032649517059326,
        "pvalueAdj": 0.8812961615097907
    },
    {
        "log2Foldchange": 1.8256185054779053,
        "pvalueAdj": 0.24273365095579313
    },
    {
        "log2Foldchange": 8.180359840393066,
        "pvalueAdj": 5.366866685848026e-7
    },
    {
        "log2Foldchange": 4.331547737121582,
        "pvalueAdj": 0.12409382977793353
    },
    {
        "log2Foldchange": 0.813087522983551,
        "pvalueAdj": 0.9938842973611571
    },
    {
        "log2Foldchange": 2.204742908477783,
        "pvalueAdj": 7.984925053884862e-64
    },
    {
        "log2Foldchange": 4.0682244300842285,
        "pvalueAdj": 3.200788274760086e-106
    },
    {
        "log2Foldchange": 2.106973886489868,
        "pvalueAdj": 6.260167157358415e-56
    },
    {
        "log2Foldchange": 3.4032015800476074,
        "pvalueAdj": 3.1389571380136343e-119
    },
    {
        "log2Foldchange": 3.6343913078308105,
        "pvalueAdj": 1.526460567321708e-137
    },
    {
        "log2Foldchange": 3.4973249435424805,
        "pvalueAdj": 2.7533295250144865e-161
    },
    {
        "log2Foldchange": 2.152918815612793,
        "pvalueAdj": 3.777820711758066e-46
    },
    {
        "log2Foldchange": 1.8759461641311646,
        "pvalueAdj": 6.7387520058319996e-46
    },
    {
        "log2Foldchange": 2.062788724899292,
        "pvalueAdj": 1.415123717468114e-46
    },
    {
        "log2Foldchange": 2.086505174636841,
        "pvalueAdj": 7.498988631792026e-66
    },
    {
        "log2Foldchange": 2.8979320526123047,
        "pvalueAdj": 1.0152666382922017e-76
    },
    {
        "log2Foldchange": 2.257690191268921,
        "pvalueAdj": 2.8394856755129314e-92
    },
    {
        "log2Foldchange": 2.02933669090271,
        "pvalueAdj": 1.4381342938891404e-47
    },
    {
        "log2Foldchange": 3.354607582092285,
        "pvalueAdj": 4.8248113746532835e-173
    },
    {
        "log2Foldchange": 3.9728646278381348,
        "pvalueAdj": 9.175482615019964e-184
    },
    {
        "log2Foldchange": 2.531480073928833,
        "pvalueAdj": 1.8897259785929053e-79
    },
    {
        "log2Foldchange": 2.263294219970703,
        "pvalueAdj": 3.661427953996245e-47
    },
    {
        "log2Foldchange": 2.2775838375091553,
        "pvalueAdj": 1.9156301768348698e-54
    },
    {
        "log2Foldchange": 2.079608201980591,
        "pvalueAdj": 7.456436689934556e-58
    },
    {
        "log2Foldchange": 2.6080663204193115,
        "pvalueAdj": 5.2770397277189734e-64
    },
    {
        "log2Foldchange": 2.034295082092285,
        "pvalueAdj": 1.3158091464578723e-65
    },
    {
        "log2Foldchange": 2.0929014682769775,
        "pvalueAdj": 1.7434279730613616e-54
    },
    {
        "log2Foldchange": 2.343421697616577,
        "pvalueAdj": 1.5463259623005745e-110
    },
    {
        "log2Foldchange": 2.2093381881713867,
        "pvalueAdj": 2.9859583730818105e-56
    },
    {
        "log2Foldchange": 1.7461391687393188,
        "pvalueAdj": 3.4895163929658505e-54
    },
    {
        "log2Foldchange": 3.9775924682617188,
        "pvalueAdj": 2.5979576823945325e-118
    },
    {
        "log2Foldchange": 4.020167827606201,
        "pvalueAdj": 3.940262980849446e-132
    },
    {
        "log2Foldchange": 3.503763198852539,
        "pvalueAdj": 1.6528632489049293e-133
    },
    {
        "log2Foldchange": 2.625159740447998,
        "pvalueAdj": 3.615204146836571e-84
    },
    {
        "log2Foldchange": 3.25948166847229,
        "pvalueAdj": 1.9177061314157802e-83
    },
    {
        "log2Foldchange": 3.098356008529663,
        "pvalueAdj": 1.6753480634662904e-78
    },
    {
        "log2Foldchange": 3.0033230781555176,
        "pvalueAdj": 2.884262415729038e-80
    },
    {
        "log2Foldchange": 4.844970226287842,
        "pvalueAdj": 1.5362514388073766e-55
    },
    {
        "log2Foldchange": 4.481184959411621,
        "pvalueAdj": 3.769464504706972e-191
    },
    {
        "log2Foldchange": 2.473616600036621,
        "pvalueAdj": 2.0257640810121384e-107
    },
    {
        "log2Foldchange": 2.4688780307769775,
        "pvalueAdj": 1.0271719626290837e-75
    },
    {
        "log2Foldchange": 13.342408180236816,
        "pvalueAdj": 0
    },
    {
        "log2Foldchange": 3.7583723068237305,
        "pvalueAdj": 2.267565215613455e-197
    },
    {
        "log2Foldchange": 3.447967052459717,
        "pvalueAdj": 7.974114591586608e-150
    },
    {
        "log2Foldchange": 5.5788068771362305,
        "pvalueAdj": 8.240193180552716e-148
    },
    {
        "log2Foldchange": 5.368076801300049,
        "pvalueAdj": 4.5809775780429434e-126
    },
    {
        "log2Foldchange": 6.797235488891602,
        "pvalueAdj": 2.5422683237015346e-125
    },
    {
        "log2Foldchange": 3.6316161155700684,
        "pvalueAdj": 1.2930136532483277e-113
    },
    {
        "log2Foldchange": 2.336184501647949,
        "pvalueAdj": 9.519249831804025e-103
    },
    {
        "log2Foldchange": 3.9643449783325195,
        "pvalueAdj": 1.0710678632047282e-99
    },
    {
        "log2Foldchange": 4.264341354370117,
        "pvalueAdj": 2.277309876471936e-99
    },
    {
        "log2Foldchange": 4.2674479484558105,
        "pvalueAdj": 4.3458251583928827e-97
    },
    {
        "log2Foldchange": 3.2196500301361084,
        "pvalueAdj": 3.800413367519762e-86
    },
    {
        "log2Foldchange": 2.971247673034668,
        "pvalueAdj": 2.798364201395149e-85
    },
    {
        "log2Foldchange": 4.398364543914795,
        "pvalueAdj": 2.711432422595307e-83
    },
    {
        "log2Foldchange": 2.7018558979034424,
        "pvalueAdj": 1.7981927659980787e-82
    },
    {
        "log2Foldchange": 3.3664605617523193,
        "pvalueAdj": 8.900395177046458e-82
    },
    {
        "log2Foldchange": 2.565990447998047,
        "pvalueAdj": 1.3302694297605615e-81
    },
    {
        "log2Foldchange": 3.2112245559692383,
        "pvalueAdj": 4.1161553252823467e-79
    },
    {
        "log2Foldchange": 2.7399978637695312,
        "pvalueAdj": 5.281063931066421e-75
    },
    {
        "log2Foldchange": 2.5879554748535156,
        "pvalueAdj": 1.1369446111073664e-71
    },
    {
        "log2Foldchange": 2.651322364807129,
        "pvalueAdj": 2.5931089929191884e-70
    },
    {
        "log2Foldchange": 5.481463432312012,
        "pvalueAdj": 2.6037333602462146e-65
    },
    {
        "log2Foldchange": 2.4212400913238525,
        "pvalueAdj": 7.626214531571608e-64
    },
    {
        "log2Foldchange": 2.1872122287750244,
        "pvalueAdj": 1.2302997546666002e-63
    },
    {
        "log2Foldchange": 2.477537155151367,
        "pvalueAdj": 4.468848698638319e-61
    },
    {
        "log2Foldchange": 4.037140369415283,
        "pvalueAdj": 5.926375749714031e-61
    },
    {
        "log2Foldchange": 2.8292157649993896,
        "pvalueAdj": 1.2196292673536187e-60
    },
    {
        "log2Foldchange": 2.2851476669311523,
        "pvalueAdj": 2.3976576335702135e-59
    },
    {
        "log2Foldchange": 2.282695770263672,
        "pvalueAdj": 5.400139413442975e-58
    },
    {
        "log2Foldchange": 2.070852279663086,
        "pvalueAdj": 9.658540932977372e-58
    },
    {
        "log2Foldchange": 1.8804264068603516,
        "pvalueAdj": 3.259191858764868e-57
    },
    {
        "log2Foldchange": 1.95876145362854,
        "pvalueAdj": 1.8871511843296134e-56
    },
    {
        "log2Foldchange": 3.6474297046661377,
        "pvalueAdj": 3.1235879284658743e-56
    },
    {
        "log2Foldchange": 2.106504201889038,
        "pvalueAdj": 3.248252405288402e-56
    },
    {
        "log2Foldchange": 2.222846508026123,
        "pvalueAdj": 2.0684628798489026e-55
    },
    {
        "log2Foldchange": 1.6512264013290405,
        "pvalueAdj": 6.674112622471447e-55
    },
    {
        "log2Foldchange": 1.9693562984466553,
        "pvalueAdj": 8.943875286333584e-55
    },
    {
        "log2Foldchange": 2.000840425491333,
        "pvalueAdj": 1.3697031423254848e-53
    },
    {
        "log2Foldchange": 2.0654680728912354,
        "pvalueAdj": 9.895300693969308e-53
    },
    {
        "log2Foldchange": 2.0881412029266357,
        "pvalueAdj": 1.2911469724588967e-52
    },
    {
        "log2Foldchange": 2.113406181335449,
        "pvalueAdj": 3.5853699807997176e-52
    },
    {
        "log2Foldchange": 2.1782634258270264,
        "pvalueAdj": 5.5563456414821955e-52
    },
    {
        "log2Foldchange": 2.4987406730651855,
        "pvalueAdj": 1.6140744589032138e-51
    },
    {
        "log2Foldchange": 2.200211763381958,
        "pvalueAdj": 2.683175256579481e-50
    },
    {
        "log2Foldchange": 1.5154441595077515,
        "pvalueAdj": 9.830856977156838e-50
    },
    {
        "log2Foldchange": 2.644683361053467,
        "pvalueAdj": 1.5012963881409237e-49
    },
    {
        "log2Foldchange": 2.0861330032348633,
        "pvalueAdj": 3.842900426731567e-49
    },
    {
        "log2Foldchange": 2.1288156509399414,
        "pvalueAdj": 2.3427344330700177e-48
    },
    {
        "log2Foldchange": 1.965607762336731,
        "pvalueAdj": 3.7892015854411527e-48
    },
    {
        "log2Foldchange": 1.9021252393722534,
        "pvalueAdj": 4.2353550910546195e-48
    },
    {
        "log2Foldchange": 2.3865034580230713,
        "pvalueAdj": 4.568725115657556e-48
    },
    {
        "log2Foldchange": 3.4802119731903076,
        "pvalueAdj": 8.367951230664752e-48
    },
    {
        "log2Foldchange": 2.117946147918701,
        "pvalueAdj": 3.3670944198543407e-47
    },
    {
        "log2Foldchange": 2.0320141315460205,
        "pvalueAdj": 5.942042480715578e-47
    },
    {
        "log2Foldchange": 2.528416395187378,
        "pvalueAdj": 1.1285202109361255e-46
    },
    {
        "log2Foldchange": 2.558823585510254,
        "pvalueAdj": 2.0587856136790993e-46
    },
    {
        "log2Foldchange": 2.235530376434326,
        "pvalueAdj": 2.3628749548889976e-46
    },
    {
        "log2Foldchange": 2.020297050476074,
        "pvalueAdj": 4.050844796540206e-46
    },
    {
        "log2Foldchange": 3.3248326778411865,
        "pvalueAdj": 8.548433908315337e-46
    },
    {
        "log2Foldchange": 2.941375732421875,
        "pvalueAdj": 1.1949767702364276e-45
    },
    {
        "log2Foldchange": 2.205721378326416,
        "pvalueAdj": 1.4220134955829276e-45
    },
    {
        "log2Foldchange": 1.8769129514694214,
        "pvalueAdj": 1.9807836227041696e-45
    },
    {
        "log2Foldchange": 2.1289784908294678,
        "pvalueAdj": 2.3848309026719735e-45
    }
];
const InlineFoldChangePlot = () => {
    return (
        <Plot
            config={{
                modeBarButtons: false,
                displaylogo: false,
                staticPlot: true
            }}
            data={[
                {
                    x: DATA.map((e) => e.log2Foldchange),
                    y: DATA.map((e) => -1 * Math.log10(e.pvalueAdj)),
                    type: 'scatter',
                    mode: 'markers',
                    marker: {color: 'lightblue'},
                    showlegend: false,
                }, {
                    x: [2],
                    y: [50],
                    marker: {
                        color: 'red',
                        size: 5
                    },
                    showlegend: false
                }]}

            layout={{
                width: 160, height: 100, margin: {l: 26, r: 0, t: 0, b: 20},

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
    );
};

export {InlineFoldChangePlot};