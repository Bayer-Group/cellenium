import {AnnotationGroup, AnnotationValue, StudyBasicsFragment} from "./generated/types";
import ColumnTable from 'arquero/dist/types/table/column-table';

export type Study = StudyBasicsFragment & {
    samplesProjectionTable: ColumnTable;
    samplesAnnotationTable: ColumnTable;
    annotationGroupMap: Map<number, AnnotationGroup>;
    annotationValueMap: Map<number, AnnotationValue>;
};

export type SelectBoxItem = {
    value:string,
    label:string
}