import { internal } from 'arquero';
import { AnnotationGrpFragment, OmicsType, StudyBasicsFragment } from './generated/types';

export type OntologyItem = {
  id: string;
  ontology: string;
  unique_id: string;
  label: string;
  parent_unique_id: string | undefined;
  children?: OntologyItem[];
};

export type Omics = {
  omicsId: number;
  omicsType: string;
  displayName: string;
  displaySymbol: string;
  value?: string;
  ontology?: string;
  taxId?: number;
  linkedGenes?: number[];
};

export type Study = StudyBasicsFragment & {
  samplesProjectionTables: Map<string, SamplesProjectionTable>;
  samplesAnnotationTable: SamplesAnnotationTable;
  studyOmicsTable: StudyOmicsTable;
  studyOmicsMap: Map<number, Omics>;
  annotationGroupMap: Map<number, AnnotationGrpFragment>;
  annotationValueMap: Map<number, AnnotationGrpFragment['annotationValuesList'][0]>;
  omicsTypes: OmicsType[];
};

export type SelectBoxItem = {
  value: string;
  label: string;
};

class DefinedTable extends internal.ColumnTable {
  static checkTable(t: internal.ColumnTable, expectedColumns: string[]) {
    expectedColumns.forEach((c) => {
      if (!t.column(c)) {
        throw Error(`table with columns ${t.columnNames()} does not match expected columns ${expectedColumns}`);
      }
    });
  }
}

export class SamplesProjectionTable extends DefinedTable {
  // samplesProjectionTable = true;

  static definedTable(t: internal.ColumnTable): SamplesProjectionTable {
    DefinedTable.checkTable(t, ['studySampleId', 'projectionX', 'projectionY']);
    return t as SamplesProjectionTable;
  }
}

export class SamplesAnnotationTable extends DefinedTable {
  // samplesAnnotationTable = true;

  static definedTable(t: internal.ColumnTable): SamplesAnnotationTable {
    DefinedTable.checkTable(t, ['studySampleId', 'annotationValueId', 'annotationGroupId']);
    return t as SamplesAnnotationTable;
  }
}

export class StudyOmicsTable extends DefinedTable {
  // StudyOmicsTable = true;

  static definedTable(t: internal.ColumnTable): StudyOmicsTable {
    DefinedTable.checkTable(t, ['omicsId', 'omicsType', 'value', 'displaySymbol', 'displayName']);
    return t as StudyOmicsTable;
  }
}

export class ExpressionTable extends DefinedTable {
  // expressionTable = true;

  static definedTable(t: internal.ColumnTable): ExpressionTable {
    DefinedTable.checkTable(t, ['value', 'omicsId', 'studySampleId']);
    return t as ExpressionTable;
  }
}
