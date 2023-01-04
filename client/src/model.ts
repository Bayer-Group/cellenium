import {StudyBasicsFragment} from "./generated/types";
import ColumnTable from 'arquero/dist/types/table/column-table';

export type Study = StudyBasicsFragment & {
    samplesProjectionTable: ColumnTable;
    samplesAnnotationTable: ColumnTable;
};
