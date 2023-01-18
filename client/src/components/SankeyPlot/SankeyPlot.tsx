import {ResponsiveSankey} from '@nivo/sankey';
import React from 'react';

const data = {
    "nodes": [
        {
            "id": "John",
            "nodeColor": "hsl(71, 70%, 50%)"
        },
        {
            "id": "Raoul",
            "nodeColor": "hsl(59, 70%, 50%)"
        },
        {
            "id": "Jane",
            "nodeColor": "hsl(54, 70%, 50%)"
        },
        {
            "id": "Marcel",
            "nodeColor": "hsl(170, 70%, 50%)"
        },
        {
            "id": "Ibrahim",
            "nodeColor": "hsl(36, 70%, 50%)"
        },
        {
            "id": "Junko",
            "nodeColor": "hsl(13, 70%, 50%)"
        }
    ],
    "links": [
        {
            "source": "Jane",
            "target": "Raoul",
            "value": 130
        },
        {
            "source": "Junko",
            "target": "John",
            "value": 179
        },
        {
            "source": "Junko",
            "target": "Jane",
            "value": 70
        },
        {
            "source": "Junko",
            "target": "Ibrahim",
            "value": 122
        },
        {
            "source": "Ibrahim",
            "target": "Jane",
            "value": 27
        },
        {
            "source": "Marcel",
            "target": "Ibrahim",
            "value": 127
        },
        {
            "source": "Marcel",
            "target": "Junko",
            "value": 155
        },
        {
            "source": "John",
            "target": "Raoul",
            "value": 103
        }
    ]
};
const SankeyPlot = () => {
    return (
        <div style={{width: '100%', height: '100%'}}>
            <ResponsiveSankey
                data={data}
                margin={{top: 40, right: 160, bottom: 40, left: 50}}
                align="justify"
                colors={{scheme: 'category10'}}
                nodeOpacity={1}
                nodeHoverOthersOpacity={0.35}
                nodeThickness={18}
                nodeSpacing={24}
                nodeBorderWidth={0}
                nodeBorderColor={{
                    from: 'color',
                    modifiers: [
                        [
                            'darker',
                            0.8
                        ]
                    ]
                }}
                nodeBorderRadius={3}
                linkOpacity={0.5}
                linkHoverOthersOpacity={0.1}
                linkContract={3}
                enableLinkGradient={true}
                labelPosition="outside"
                labelOrientation="vertical"
                labelPadding={16}
                labelTextColor={{
                    from: 'color',
                    modifiers: [
                        [
                            'darker',
                            1
                        ]
                    ]
                }}
                
            />
        </div>
    );
};

export {SankeyPlot};