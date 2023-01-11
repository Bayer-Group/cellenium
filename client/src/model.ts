import {AnnotationGroup, AnnotationValue, StudyBasicsFragment} from "./generated/types";
import {internal} from 'arquero';


export type OntologyItem = {
    id: string;
    ontology: string;
    unique_id: string;
    label: string;
    parent_unique_id: string|undefined;
    children?: OntologyItem[];
}


export type Gene = {
    omicsId: number;
    displayName: string;
    displaySymbol: string;
}

export type Study = StudyBasicsFragment & {
    samplesProjectionTable: SamplesProjectionTable;
    samplesAnnotationTable: SamplesAnnotationTable;
    annotationGroupMap: Map<number, AnnotationGroup>;
    annotationValueMap: Map<number, AnnotationValue>;

};

export type SelectBoxItem = {
    value: string,
    label: string
}

class DefinedTable extends internal.ColumnTable {
    static checkTable(t: internal.ColumnTable, expectedColumns: string[]) {
        expectedColumns.forEach((c) => {
            if (!t.column(c)) {
                throw Error(`table with columns ${t.columnNames()} does not match expected columns ${expectedColumns}`);
            }
        })
    }
}

export class SamplesProjectionTable extends DefinedTable {
    samplesProjectionTable = true;

    static definedTable(t: internal.ColumnTable): SamplesProjectionTable {
        DefinedTable.checkTable(t, ['studySampleId', 'projectionX', 'projectionY']);
        return t as SamplesProjectionTable;
    }
};

export class SamplesAnnotationTable extends DefinedTable {
    samplesAnnotationTable = true;

    static definedTable(t: internal.ColumnTable): SamplesAnnotationTable {
        DefinedTable.checkTable(t, ['studySampleId', 'annotationValueId', 'annotationGroupId']);
        return t as SamplesAnnotationTable;
    }
};

export class ExpressionTable extends DefinedTable {
    expressionTable = true;

    static definedTable(t: internal.ColumnTable): ExpressionTable {
        DefinedTable.checkTable(t, ['value', 'omicsId', 'studySampleId']);
        return t as ExpressionTable;
    }
};
