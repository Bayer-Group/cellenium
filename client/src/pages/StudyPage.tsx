import React, {useEffect, useState} from 'react';
import {useRecoilState, useRecoilValue} from "recoil";
import {pageState, selectedGenesState} from "../atoms";
import {useSetStudyFromUrl} from "../hooks";
import CellMarkerAnalysis from "./CellMarkerAnalysis";
import ExpressionAnalysis from "./ExpressionAnalysis";
import CoexpressionAnalysis from "./CoexpressionAnalysis";
import AnnotationComparison from "./AnnotationComparison";
import UserAnnotation from "./UserAnnotation";

export function StudyPage() {
    useSetStudyFromUrl();
    const page = useRecoilValue(pageState);

    switch (page) {
        case 'CellMarkerAnalysis':
            return <CellMarkerAnalysis/>;
        case 'ExpressionAnalysis':
            return <ExpressionAnalysis/>;
        case 'CoexpressionAnalysis':
            return <CoexpressionAnalysis/>;
        case 'AnnotationComparison':
            return <AnnotationComparison/>;
        case 'UserAnnotation':
            return <UserAnnotation/>;
    }
    return <div>page {page} not defined</div>;
}
