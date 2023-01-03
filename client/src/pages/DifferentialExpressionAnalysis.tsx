import React from 'react';
import {AnnotationGroup, AnnotationGroupSelectBox, DEGTable, LeftSidePanel, RightSidePanel} from "../components";
import {Group} from "@mantine/core";

const ANNOTATIONS = [
    {label: "bone cell", color: "#1f77b4"},
    {label: "chondrocyte", color: "#ff7f0e"},
    {label: "endothelial cell", color: "#2ca02c"},
    {label: "endothelial cell of artery", color: "#d62728"},
    {label: "fibroblast", color: "#9467bd"},
    {label: "mesenchymal stem cell", color: "#8c564b"},
    {label: "pericyte cell", color: "#e377c2"}
]

const DEG = [
    {
        "name": "Athena Weissnat",
        "company": "Little - Rippin",
        "email": "Elouise.Prohaska@yahoo.com"
    },
    {
        "name": "Deangelo Runolfsson",
        "company": "Greenfelder - Krajcik",
        "email": "Kadin_Trantow87@yahoo.com"
    },
    {
        "name": "Danny Carter",
        "company": "Kohler and Sons",
        "email": "Marina3@hotmail.com"
    },
    {
        "name": "Trace Tremblay PhD",
        "company": "Crona, Aufderhar and Senger",
        "email": "Antonina.Pouros@yahoo.com"
    },
    {
        "name": "Derek Dibbert",
        "company": "Gottlieb LLC",
        "email": "Abagail29@hotmail.com"
    },
    {
        "name": "Viola Bernhard",
        "company": "Funk, Rohan and Kreiger",
        "email": "Jamie23@hotmail.com"
    },
    {
        "name": "Austin Jacobi",
        "company": "Botsford - Corwin",
        "email": "Genesis42@yahoo.com"
    },
    {
        "name": "Hershel Mosciski",
        "company": "Okuneva, Farrell and Kilback",
        "email": "Idella.Stehr28@yahoo.com"
    },
    {
        "name": "Mylene Ebert",
        "company": "Kirlin and Sons",
        "email": "Hildegard17@hotmail.com"
    },
    {
        "name": "Lou Trantow",
        "company": "Parisian - Lemke",
        "email": "Hillard.Barrows1@hotmail.com"
    },
    {
        "name": "Dariana Weimann",
        "company": "Schowalter - Donnelly",
        "email": "Colleen80@gmail.com"
    },
    {
        "name": "Dr. Christy Herman",
        "company": "VonRueden - Labadie",
        "email": "Lilyan98@gmail.com"
    },
    {
        "name": "Katelin Schuster",
        "company": "Jacobson - Smitham",
        "email": "Erich_Brekke76@gmail.com"
    },
    {
        "name": "Melyna Macejkovic",
        "company": "Schuster LLC",
        "email": "Kylee4@yahoo.com"
    },
    {
        "name": "Pinkie Rice",
        "company": "Wolf, Trantow and Zulauf",
        "email": "Fiona.Kutch@hotmail.com"
    },
    {
        "name": "Brain Kreiger",
        "company": "Lueilwitz Group",
        "email": "Rico98@hotmail.com"
    }
];

function DifferentialExpressionAnalysis() {
    return (
        <Group position={'apart'}>
            <LeftSidePanel>
                <AnnotationGroupSelectBox/>
                <AnnotationGroup annotations={ANNOTATIONS}/>
            </LeftSidePanel>

            <RightSidePanel>
                <DEGTable data={DEG}/>
            </RightSidePanel>
        </Group>
    );
};

export default DifferentialExpressionAnalysis;