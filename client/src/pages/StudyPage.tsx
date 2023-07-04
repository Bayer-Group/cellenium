import {useRecoilValue} from "recoil";
import {pageState} from "../atoms";
import {useSetStudyFromUrl} from "../hooks";
import CellMarkerAnalysis from "./CellMarkerAnalysis";
import ExpressionAnalysis from "./ExpressionAnalysis";
import CoexpressionAnalysis from "./CoexpressionAnalysis";
import AnnotationComparison from "./AnnotationComparison";
import UserAnnotation from "./UserAnnotation";
import CelltypeDiscovery from "./CelltypeDiscovery";

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
        case 'CelltypeDiscovery':
            return <CelltypeDiscovery/>;
        case 'AnnotationComparison':
            return <AnnotationComparison/>;
        case 'UserAnnotation':
            return <UserAnnotation/>;
    }
    return <div>page {page} not defined</div>;
}
