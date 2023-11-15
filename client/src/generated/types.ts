/* eslint-disable */
import { gql } from '@apollo/client';
import * as Apollo from '@apollo/client';
export type Maybe<T> = T;
export type InputMaybe<T> = T;
export type Exact<T extends { [key: string]: unknown }> = { [K in keyof T]: T[K] };
export type MakeOptional<T, K extends keyof T> = Omit<T, K> & { [SubKey in K]?: Maybe<T[SubKey]> };
export type MakeMaybe<T, K extends keyof T> = Omit<T, K> & { [SubKey in K]: Maybe<T[SubKey]> };
const defaultOptions = {} as const;
/** All built-in and custom scalars, mapped to their actual values */
export type Scalars = {
  ID: string;
  String: string;
  Boolean: boolean;
  Int: number;
  Float: number;
  BigInt: any;
  Cursor: any;
  Datetime: any;
  JSON: any;
};

export type AnnotationGroup = Node & {
  __typename?: 'AnnotationGroup';
  annotationGroupId: Scalars['Int'];
  /** Reads and enables pagination through a set of `AnnotationValue`. */
  annotationValues: AnnotationValuesConnection;
  /** Reads and enables pagination through a set of `AnnotationValue`. */
  annotationValuesList: Array<AnnotationValue>;
  displayGroup: Scalars['String'];
  h5AdColumn: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads and enables pagination through a set of `ReferenceStudy`. */
  referenceStudiesByCelltypeAnnotationGroupId: ReferenceStudiesConnection;
  /** Reads and enables pagination through a set of `ReferenceStudy`. */
  referenceStudiesByCelltypeAnnotationGroupIdList: Array<ReferenceStudy>;
  /** Reads and enables pagination through a set of `ReferenceStudy`. */
  referenceStudiesByTissueAnnotationGroupId: ReferenceStudiesConnection;
  /** Reads and enables pagination through a set of `ReferenceStudy`. */
  referenceStudiesByTissueAnnotationGroupIdList: Array<ReferenceStudy>;
  /** Reads and enables pagination through a set of `StudyAnnotationGroupUi`. */
  studyAnnotationGroupUis: StudyAnnotationGroupUisConnection;
  /** Reads and enables pagination through a set of `StudyAnnotationGroupUi`. */
  studyAnnotationGroupUisList: Array<StudyAnnotationGroupUi>;
  /** Reads and enables pagination through a set of `UserAnnotationGroup`. */
  userAnnotationGroupsBySavedAsAnnotationGroupId: UserAnnotationGroupsConnection;
  /** Reads and enables pagination through a set of `UserAnnotationGroup`. */
  userAnnotationGroupsBySavedAsAnnotationGroupIdList: Array<UserAnnotationGroup>;
};


export type AnnotationGroupAnnotationValuesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<AnnotationValueCondition>;
  filter: InputMaybe<AnnotationValueFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<AnnotationValuesOrderBy>>;
};


export type AnnotationGroupAnnotationValuesListArgs = {
  condition: InputMaybe<AnnotationValueCondition>;
  filter: InputMaybe<AnnotationValueFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<AnnotationValuesOrderBy>>;
};


export type AnnotationGroupReferenceStudiesByCelltypeAnnotationGroupIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ReferenceStudyCondition>;
  filter: InputMaybe<ReferenceStudyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ReferenceStudiesOrderBy>>;
};


export type AnnotationGroupReferenceStudiesByCelltypeAnnotationGroupIdListArgs = {
  condition: InputMaybe<ReferenceStudyCondition>;
  filter: InputMaybe<ReferenceStudyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ReferenceStudiesOrderBy>>;
};


export type AnnotationGroupReferenceStudiesByTissueAnnotationGroupIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ReferenceStudyCondition>;
  filter: InputMaybe<ReferenceStudyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ReferenceStudiesOrderBy>>;
};


export type AnnotationGroupReferenceStudiesByTissueAnnotationGroupIdListArgs = {
  condition: InputMaybe<ReferenceStudyCondition>;
  filter: InputMaybe<ReferenceStudyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ReferenceStudiesOrderBy>>;
};


export type AnnotationGroupStudyAnnotationGroupUisArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyAnnotationGroupUiCondition>;
  filter: InputMaybe<StudyAnnotationGroupUiFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyAnnotationGroupUisOrderBy>>;
};


export type AnnotationGroupStudyAnnotationGroupUisListArgs = {
  condition: InputMaybe<StudyAnnotationGroupUiCondition>;
  filter: InputMaybe<StudyAnnotationGroupUiFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyAnnotationGroupUisOrderBy>>;
};


export type AnnotationGroupUserAnnotationGroupsBySavedAsAnnotationGroupIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<UserAnnotationGroupCondition>;
  filter: InputMaybe<UserAnnotationGroupFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<UserAnnotationGroupsOrderBy>>;
};


export type AnnotationGroupUserAnnotationGroupsBySavedAsAnnotationGroupIdListArgs = {
  condition: InputMaybe<UserAnnotationGroupCondition>;
  filter: InputMaybe<UserAnnotationGroupFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<UserAnnotationGroupsOrderBy>>;
};

/**
 * A condition to be used against `AnnotationGroup` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type AnnotationGroupCondition = {
  /** Checks for equality with the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `displayGroup` field. */
  displayGroup: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `h5AdColumn` field. */
  h5AdColumn: InputMaybe<Scalars['String']>;
};

/** A filter to be used against `AnnotationGroup` object types. All fields are combined with a logical ‘and.’ */
export type AnnotationGroupFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<AnnotationGroupFilter>>;
  /** Filter by the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<IntFilter>;
  /** Filter by the object’s `displayGroup` field. */
  displayGroup: InputMaybe<StringFilter>;
  /** Filter by the object’s `h5AdColumn` field. */
  h5AdColumn: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<AnnotationGroupFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<AnnotationGroupFilter>>;
};

/** An input for mutations affecting `AnnotationGroup` */
export type AnnotationGroupInput = {
  annotationGroupId: InputMaybe<Scalars['Int']>;
  displayGroup: Scalars['String'];
  h5AdColumn: Scalars['String'];
};

/** Represents an update to a `AnnotationGroup`. Fields that are set will be updated. */
export type AnnotationGroupPatch = {
  annotationGroupId: InputMaybe<Scalars['Int']>;
  displayGroup: InputMaybe<Scalars['String']>;
  h5AdColumn: InputMaybe<Scalars['String']>;
};

/** A connection to a list of `AnnotationGroup` values. */
export type AnnotationGroupsConnection = {
  __typename?: 'AnnotationGroupsConnection';
  /** A list of edges which contains the `AnnotationGroup` and cursor to aid in pagination. */
  edges: Array<AnnotationGroupsEdge>;
  /** A list of `AnnotationGroup` objects. */
  nodes: Array<Maybe<AnnotationGroup>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `AnnotationGroup` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `AnnotationGroup` edge in the connection. */
export type AnnotationGroupsEdge = {
  __typename?: 'AnnotationGroupsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `AnnotationGroup` at the end of the edge. */
  node: Maybe<AnnotationGroup>;
};

/** Methods to use when ordering `AnnotationGroup`. */
export enum AnnotationGroupsOrderBy {
  AnnotationGroupIdAsc = 'ANNOTATION_GROUP_ID_ASC',
  AnnotationGroupIdDesc = 'ANNOTATION_GROUP_ID_DESC',
  DisplayGroupAsc = 'DISPLAY_GROUP_ASC',
  DisplayGroupDesc = 'DISPLAY_GROUP_DESC',
  H5AdColumnAsc = 'H5AD_COLUMN_ASC',
  H5AdColumnDesc = 'H5AD_COLUMN_DESC',
  Natural = 'NATURAL',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC'
}

export type AnnotationValue = Node & {
  __typename?: 'AnnotationValue';
  /** Reads a single `AnnotationGroup` that is related to this `AnnotationValue`. */
  annotationGroup: Maybe<AnnotationGroup>;
  annotationGroupId: Scalars['Int'];
  annotationValueId: Scalars['Int'];
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressions: DifferentialExpressionsConnection;
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressionsList: Array<DifferentialExpression>;
  displayValue: Scalars['String'];
  h5AdValue: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads and enables pagination through a set of `StudySampleAnnotation`. */
  studySampleAnnotations: StudySampleAnnotationsConnection;
  /** Reads and enables pagination through a set of `StudySampleAnnotation`. */
  studySampleAnnotationsList: Array<StudySampleAnnotation>;
};


export type AnnotationValueDifferentialExpressionsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<DifferentialExpressionCondition>;
  filter: InputMaybe<DifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


export type AnnotationValueDifferentialExpressionsListArgs = {
  condition: InputMaybe<DifferentialExpressionCondition>;
  filter: InputMaybe<DifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


export type AnnotationValueStudySampleAnnotationsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudySampleAnnotationCondition>;
  filter: InputMaybe<StudySampleAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySampleAnnotationsOrderBy>>;
};


export type AnnotationValueStudySampleAnnotationsListArgs = {
  condition: InputMaybe<StudySampleAnnotationCondition>;
  filter: InputMaybe<StudySampleAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySampleAnnotationsOrderBy>>;
};

/**
 * A condition to be used against `AnnotationValue` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type AnnotationValueCondition = {
  /** Checks for equality with the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `displayValue` field. */
  displayValue: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `h5AdValue` field. */
  h5AdValue: InputMaybe<Scalars['String']>;
};

/** A connection to a list of `AnnotationValueCoocurrenceRecord` values. */
export type AnnotationValueCoocurrenceConnection = {
  __typename?: 'AnnotationValueCoocurrenceConnection';
  /** A list of edges which contains the `AnnotationValueCoocurrenceRecord` and cursor to aid in pagination. */
  edges: Array<AnnotationValueCoocurrenceEdge>;
  /** A list of `AnnotationValueCoocurrenceRecord` objects. */
  nodes: Array<Maybe<AnnotationValueCoocurrenceRecord>>;
  /** The count of *all* `AnnotationValueCoocurrenceRecord` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `AnnotationValueCoocurrenceRecord` edge in the connection. */
export type AnnotationValueCoocurrenceEdge = {
  __typename?: 'AnnotationValueCoocurrenceEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `AnnotationValueCoocurrenceRecord` at the end of the edge. */
  node: Maybe<AnnotationValueCoocurrenceRecord>;
};

/** The return type of our `annotationValueCoocurrence` query. */
export type AnnotationValueCoocurrenceRecord = {
  __typename?: 'AnnotationValueCoocurrenceRecord';
  annotationValueId1: Maybe<Scalars['Int']>;
  annotationValueId2: Maybe<Scalars['Int']>;
  occurrence: Maybe<Scalars['Int']>;
};

/** A filter to be used against `AnnotationValueCoocurrenceRecord` object types. All fields are combined with a logical ‘and.’ */
export type AnnotationValueCoocurrenceRecordFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<AnnotationValueCoocurrenceRecordFilter>>;
  /** Filter by the object’s `annotationValueId1` field. */
  annotationValueId1: InputMaybe<IntFilter>;
  /** Filter by the object’s `annotationValueId2` field. */
  annotationValueId2: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not: InputMaybe<AnnotationValueCoocurrenceRecordFilter>;
  /** Filter by the object’s `occurrence` field. */
  occurrence: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<AnnotationValueCoocurrenceRecordFilter>>;
};

/** A filter to be used against `AnnotationValue` object types. All fields are combined with a logical ‘and.’ */
export type AnnotationValueFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<AnnotationValueFilter>>;
  /** Filter by the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<IntFilter>;
  /** Filter by the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<IntFilter>;
  /** Filter by the object’s `displayValue` field. */
  displayValue: InputMaybe<StringFilter>;
  /** Filter by the object’s `h5AdValue` field. */
  h5AdValue: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<AnnotationValueFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<AnnotationValueFilter>>;
};

/** An input for mutations affecting `AnnotationValue` */
export type AnnotationValueInput = {
  annotationGroupId: Scalars['Int'];
  annotationValueId: InputMaybe<Scalars['Int']>;
  displayValue: Scalars['String'];
  h5AdValue: Scalars['String'];
};

/** Represents an update to a `AnnotationValue`. Fields that are set will be updated. */
export type AnnotationValuePatch = {
  annotationGroupId: InputMaybe<Scalars['Int']>;
  annotationValueId: InputMaybe<Scalars['Int']>;
  displayValue: InputMaybe<Scalars['String']>;
  h5AdValue: InputMaybe<Scalars['String']>;
};

/** A connection to a list of `AnnotationValue` values. */
export type AnnotationValuesConnection = {
  __typename?: 'AnnotationValuesConnection';
  /** A list of edges which contains the `AnnotationValue` and cursor to aid in pagination. */
  edges: Array<AnnotationValuesEdge>;
  /** A list of `AnnotationValue` objects. */
  nodes: Array<Maybe<AnnotationValue>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `AnnotationValue` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `AnnotationValue` edge in the connection. */
export type AnnotationValuesEdge = {
  __typename?: 'AnnotationValuesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `AnnotationValue` at the end of the edge. */
  node: Maybe<AnnotationValue>;
};

/** Methods to use when ordering `AnnotationValue`. */
export enum AnnotationValuesOrderBy {
  AnnotationGroupIdAsc = 'ANNOTATION_GROUP_ID_ASC',
  AnnotationGroupIdDesc = 'ANNOTATION_GROUP_ID_DESC',
  AnnotationValueIdAsc = 'ANNOTATION_VALUE_ID_ASC',
  AnnotationValueIdDesc = 'ANNOTATION_VALUE_ID_DESC',
  DisplayValueAsc = 'DISPLAY_VALUE_ASC',
  DisplayValueDesc = 'DISPLAY_VALUE_DESC',
  H5AdValueAsc = 'H5AD_VALUE_ASC',
  H5AdValueDesc = 'H5AD_VALUE_DESC',
  Natural = 'NATURAL',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC'
}

export type ApiDifferentialExpression = {
  __typename?: 'ApiDifferentialExpression';
  annotationDisplayValue: Maybe<Scalars['String']>;
  annotationGroupColumn: Maybe<Scalars['String']>;
  annotationGroupDisplay: Maybe<Scalars['String']>;
  annotationValue: Maybe<Scalars['String']>;
  celleniumInternalOmicsId: Maybe<Scalars['Int']>;
  ensemblGeneId: Maybe<Scalars['String']>;
  entrezGeneIds: Maybe<Array<Maybe<Scalars['String']>>>;
  hgncSymbols: Maybe<Array<Maybe<Scalars['String']>>>;
  log2Foldchange: Maybe<Scalars['Float']>;
  omicsDetails: Maybe<ApiOmic>;
  pvalue: Maybe<Scalars['Float']>;
  pvalueAdj: Maybe<Scalars['Float']>;
  score: Maybe<Scalars['Float']>;
  studyId: Maybe<Scalars['Int']>;
  studyName: Maybe<Scalars['String']>;
  studyTaxId: Maybe<Scalars['String']>;
};

/**
 * A condition to be used against `ApiDifferentialExpression` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type ApiDifferentialExpressionCondition = {
  /** Checks for equality with the object’s `annotationDisplayValue` field. */
  annotationDisplayValue: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `annotationGroupColumn` field. */
  annotationGroupColumn: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `annotationGroupDisplay` field. */
  annotationGroupDisplay: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `annotationValue` field. */
  annotationValue: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `celleniumInternalOmicsId` field. */
  celleniumInternalOmicsId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `ensemblGeneId` field. */
  ensemblGeneId: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `entrezGeneIds` field. */
  entrezGeneIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `hgncSymbols` field. */
  hgncSymbols: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `log2Foldchange` field. */
  log2Foldchange: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `omicsDetails` field. */
  omicsDetails: InputMaybe<ApiOmicInput>;
  /** Checks for equality with the object’s `pvalue` field. */
  pvalue: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `pvalueAdj` field. */
  pvalueAdj: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `score` field. */
  score: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyName` field. */
  studyName: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyTaxId` field. */
  studyTaxId: InputMaybe<Scalars['String']>;
};

/** A filter to be used against `ApiDifferentialExpression` object types. All fields are combined with a logical ‘and.’ */
export type ApiDifferentialExpressionFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ApiDifferentialExpressionFilter>>;
  /** Filter by the object’s `annotationDisplayValue` field. */
  annotationDisplayValue: InputMaybe<StringFilter>;
  /** Filter by the object’s `annotationGroupColumn` field. */
  annotationGroupColumn: InputMaybe<StringFilter>;
  /** Filter by the object’s `annotationGroupDisplay` field. */
  annotationGroupDisplay: InputMaybe<StringFilter>;
  /** Filter by the object’s `annotationValue` field. */
  annotationValue: InputMaybe<StringFilter>;
  /** Filter by the object’s `celleniumInternalOmicsId` field. */
  celleniumInternalOmicsId: InputMaybe<IntFilter>;
  /** Filter by the object’s `ensemblGeneId` field. */
  ensemblGeneId: InputMaybe<StringFilter>;
  /** Filter by the object’s `entrezGeneIds` field. */
  entrezGeneIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `hgncSymbols` field. */
  hgncSymbols: InputMaybe<StringListFilter>;
  /** Filter by the object’s `log2Foldchange` field. */
  log2Foldchange: InputMaybe<FloatFilter>;
  /** Negates the expression. */
  not: InputMaybe<ApiDifferentialExpressionFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ApiDifferentialExpressionFilter>>;
  /** Filter by the object’s `pvalue` field. */
  pvalue: InputMaybe<FloatFilter>;
  /** Filter by the object’s `pvalueAdj` field. */
  pvalueAdj: InputMaybe<FloatFilter>;
  /** Filter by the object’s `score` field. */
  score: InputMaybe<FloatFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyName` field. */
  studyName: InputMaybe<StringFilter>;
  /** Filter by the object’s `studyTaxId` field. */
  studyTaxId: InputMaybe<StringFilter>;
};

/** A connection to a list of `ApiDifferentialExpression` values. */
export type ApiDifferentialExpressionsConnection = {
  __typename?: 'ApiDifferentialExpressionsConnection';
  /** A list of edges which contains the `ApiDifferentialExpression` and cursor to aid in pagination. */
  edges: Array<ApiDifferentialExpressionsEdge>;
  /** A list of `ApiDifferentialExpression` objects. */
  nodes: Array<Maybe<ApiDifferentialExpression>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ApiDifferentialExpression` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ApiDifferentialExpression` edge in the connection. */
export type ApiDifferentialExpressionsEdge = {
  __typename?: 'ApiDifferentialExpressionsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ApiDifferentialExpression` at the end of the edge. */
  node: Maybe<ApiDifferentialExpression>;
};

/** Methods to use when ordering `ApiDifferentialExpression`. */
export enum ApiDifferentialExpressionsOrderBy {
  AnnotationDisplayValueAsc = 'ANNOTATION_DISPLAY_VALUE_ASC',
  AnnotationDisplayValueDesc = 'ANNOTATION_DISPLAY_VALUE_DESC',
  AnnotationGroupColumnAsc = 'ANNOTATION_GROUP_COLUMN_ASC',
  AnnotationGroupColumnDesc = 'ANNOTATION_GROUP_COLUMN_DESC',
  AnnotationGroupDisplayAsc = 'ANNOTATION_GROUP_DISPLAY_ASC',
  AnnotationGroupDisplayDesc = 'ANNOTATION_GROUP_DISPLAY_DESC',
  AnnotationValueAsc = 'ANNOTATION_VALUE_ASC',
  AnnotationValueDesc = 'ANNOTATION_VALUE_DESC',
  CelleniumInternalOmicsIdAsc = 'CELLENIUM_INTERNAL_OMICS_ID_ASC',
  CelleniumInternalOmicsIdDesc = 'CELLENIUM_INTERNAL_OMICS_ID_DESC',
  EnsemblGeneIdAsc = 'ENSEMBL_GENE_ID_ASC',
  EnsemblGeneIdDesc = 'ENSEMBL_GENE_ID_DESC',
  EntrezGeneIdsAsc = 'ENTREZ_GENE_IDS_ASC',
  EntrezGeneIdsDesc = 'ENTREZ_GENE_IDS_DESC',
  HgncSymbolsAsc = 'HGNC_SYMBOLS_ASC',
  HgncSymbolsDesc = 'HGNC_SYMBOLS_DESC',
  Log2FoldchangeAsc = 'LOG2_FOLDCHANGE_ASC',
  Log2FoldchangeDesc = 'LOG2_FOLDCHANGE_DESC',
  Natural = 'NATURAL',
  OmicsDetailsAsc = 'OMICS_DETAILS_ASC',
  OmicsDetailsDesc = 'OMICS_DETAILS_DESC',
  PvalueAdjAsc = 'PVALUE_ADJ_ASC',
  PvalueAdjDesc = 'PVALUE_ADJ_DESC',
  PvalueAsc = 'PVALUE_ASC',
  PvalueDesc = 'PVALUE_DESC',
  ScoreAsc = 'SCORE_ASC',
  ScoreDesc = 'SCORE_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudyNameAsc = 'STUDY_NAME_ASC',
  StudyNameDesc = 'STUDY_NAME_DESC',
  StudyTaxIdAsc = 'STUDY_TAX_ID_ASC',
  StudyTaxIdDesc = 'STUDY_TAX_ID_DESC'
}

export type ApiExpressionByAnnotation = {
  __typename?: 'ApiExpressionByAnnotation';
  annotationDisplayValue: Maybe<Scalars['String']>;
  annotationValueId: Maybe<Scalars['Int']>;
  boxplotParams: Maybe<BoxplotValue>;
  ensemblGeneId: Maybe<Scalars['String']>;
  entrezGeneIds: Maybe<Array<Maybe<Scalars['String']>>>;
  exprSamplesFraction: Maybe<Scalars['Float']>;
  hgncSymbols: Maybe<Array<Maybe<Scalars['String']>>>;
  layer: Maybe<Scalars['String']>;
  mean: Maybe<Scalars['Float']>;
  median: Maybe<Scalars['Float']>;
  omicsId: Maybe<Scalars['Int']>;
  q3: Maybe<Scalars['Float']>;
  studyId: Maybe<Scalars['Int']>;
  studyLayerId: Maybe<Scalars['Int']>;
  valueCount: Maybe<Scalars['Int']>;
  values: Maybe<Array<Maybe<Scalars['Float']>>>;
};

/** A filter to be used against `ApiExpressionByAnnotation` object types. All fields are combined with a logical ‘and.’ */
export type ApiExpressionByAnnotationFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ApiExpressionByAnnotationFilter>>;
  /** Filter by the object’s `annotationDisplayValue` field. */
  annotationDisplayValue: InputMaybe<StringFilter>;
  /** Filter by the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<IntFilter>;
  /** Filter by the object’s `boxplotParams` field. */
  boxplotParams: InputMaybe<BoxplotValueFilter>;
  /** Filter by the object’s `ensemblGeneId` field. */
  ensemblGeneId: InputMaybe<StringFilter>;
  /** Filter by the object’s `entrezGeneIds` field. */
  entrezGeneIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `exprSamplesFraction` field. */
  exprSamplesFraction: InputMaybe<FloatFilter>;
  /** Filter by the object’s `hgncSymbols` field. */
  hgncSymbols: InputMaybe<StringListFilter>;
  /** Filter by the object’s `layer` field. */
  layer: InputMaybe<StringFilter>;
  /** Filter by the object’s `mean` field. */
  mean: InputMaybe<FloatFilter>;
  /** Filter by the object’s `median` field. */
  median: InputMaybe<FloatFilter>;
  /** Negates the expression. */
  not: InputMaybe<ApiExpressionByAnnotationFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ApiExpressionByAnnotationFilter>>;
  /** Filter by the object’s `q3` field. */
  q3: InputMaybe<FloatFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId: InputMaybe<IntFilter>;
  /** Filter by the object’s `valueCount` field. */
  valueCount: InputMaybe<IntFilter>;
  /** Filter by the object’s `values` field. */
  values: InputMaybe<FloatListFilter>;
};

/** A connection to a list of `ApiExpressionByAnnotation` values. */
export type ApiExpressionByAnnotationsConnection = {
  __typename?: 'ApiExpressionByAnnotationsConnection';
  /** A list of edges which contains the `ApiExpressionByAnnotation` and cursor to aid in pagination. */
  edges: Array<ApiExpressionByAnnotationsEdge>;
  /** A list of `ApiExpressionByAnnotation` objects. */
  nodes: Array<Maybe<ApiExpressionByAnnotation>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ApiExpressionByAnnotation` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ApiExpressionByAnnotation` edge in the connection. */
export type ApiExpressionByAnnotationsEdge = {
  __typename?: 'ApiExpressionByAnnotationsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ApiExpressionByAnnotation` at the end of the edge. */
  node: Maybe<ApiExpressionByAnnotation>;
};

export type ApiOmic = {
  __typename?: 'ApiOmic';
  displayName: Maybe<Scalars['String']>;
  displaySymbol: Maybe<Scalars['String']>;
  ensemblGeneId: Maybe<Scalars['String']>;
  entrezGeneIds: Maybe<Array<Maybe<Scalars['String']>>>;
  hgncSymbols: Maybe<Array<Maybe<Scalars['String']>>>;
  linkedGenes: Maybe<Array<Maybe<Scalars['Int']>>>;
  omicsId: Maybe<Scalars['Int']>;
  omicsType: Maybe<OmicsType>;
  region: Maybe<Scalars['String']>;
  taxId: Maybe<Scalars['Int']>;
};

/** A condition to be used against `ApiOmic` object types. All fields are tested for equality and combined with a logical ‘and.’ */
export type ApiOmicCondition = {
  /** Checks for equality with the object’s `displayName` field. */
  displayName: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `ensemblGeneId` field. */
  ensemblGeneId: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `entrezGeneIds` field. */
  entrezGeneIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `hgncSymbols` field. */
  hgncSymbols: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `linkedGenes` field. */
  linkedGenes: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `omicsType` field. */
  omicsType: InputMaybe<OmicsType>;
  /** Checks for equality with the object’s `region` field. */
  region: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `taxId` field. */
  taxId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `ApiOmic` object types. All fields are combined with a logical ‘and.’ */
export type ApiOmicFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ApiOmicFilter>>;
  /** Filter by the object’s `displayName` field. */
  displayName: InputMaybe<StringFilter>;
  /** Filter by the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<StringFilter>;
  /** Filter by the object’s `ensemblGeneId` field. */
  ensemblGeneId: InputMaybe<StringFilter>;
  /** Filter by the object’s `entrezGeneIds` field. */
  entrezGeneIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `hgncSymbols` field. */
  hgncSymbols: InputMaybe<StringListFilter>;
  /** Filter by the object’s `linkedGenes` field. */
  linkedGenes: InputMaybe<IntListFilter>;
  /** Negates the expression. */
  not: InputMaybe<ApiOmicFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Filter by the object’s `omicsType` field. */
  omicsType: InputMaybe<OmicsTypeFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ApiOmicFilter>>;
  /** Filter by the object’s `region` field. */
  region: InputMaybe<StringFilter>;
  /** Filter by the object’s `taxId` field. */
  taxId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `ApiOmic` */
export type ApiOmicInput = {
  displayName: InputMaybe<Scalars['String']>;
  displaySymbol: InputMaybe<Scalars['String']>;
  ensemblGeneId: InputMaybe<Scalars['String']>;
  entrezGeneIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  hgncSymbols: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  linkedGenes: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  omicsId: InputMaybe<Scalars['Int']>;
  omicsType: InputMaybe<OmicsType>;
  region: InputMaybe<Scalars['String']>;
  taxId: InputMaybe<Scalars['Int']>;
};

/** A connection to a list of `ApiOmic` values. */
export type ApiOmicsConnection = {
  __typename?: 'ApiOmicsConnection';
  /** A list of edges which contains the `ApiOmic` and cursor to aid in pagination. */
  edges: Array<ApiOmicsEdge>;
  /** A list of `ApiOmic` objects. */
  nodes: Array<Maybe<ApiOmic>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ApiOmic` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ApiOmic` edge in the connection. */
export type ApiOmicsEdge = {
  __typename?: 'ApiOmicsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ApiOmic` at the end of the edge. */
  node: Maybe<ApiOmic>;
};

/** Methods to use when ordering `ApiOmic`. */
export enum ApiOmicsOrderBy {
  DisplayNameAsc = 'DISPLAY_NAME_ASC',
  DisplayNameDesc = 'DISPLAY_NAME_DESC',
  DisplaySymbolAsc = 'DISPLAY_SYMBOL_ASC',
  DisplaySymbolDesc = 'DISPLAY_SYMBOL_DESC',
  EnsemblGeneIdAsc = 'ENSEMBL_GENE_ID_ASC',
  EnsemblGeneIdDesc = 'ENSEMBL_GENE_ID_DESC',
  EntrezGeneIdsAsc = 'ENTREZ_GENE_IDS_ASC',
  EntrezGeneIdsDesc = 'ENTREZ_GENE_IDS_DESC',
  HgncSymbolsAsc = 'HGNC_SYMBOLS_ASC',
  HgncSymbolsDesc = 'HGNC_SYMBOLS_DESC',
  LinkedGenesAsc = 'LINKED_GENES_ASC',
  LinkedGenesDesc = 'LINKED_GENES_DESC',
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  OmicsTypeAsc = 'OMICS_TYPE_ASC',
  OmicsTypeDesc = 'OMICS_TYPE_DESC',
  RegionAsc = 'REGION_ASC',
  RegionDesc = 'REGION_DESC',
  TaxIdAsc = 'TAX_ID_ASC',
  TaxIdDesc = 'TAX_ID_DESC'
}

export type ApiStudiesBase = {
  __typename?: 'ApiStudiesBase';
  annotations: Maybe<Array<Maybe<ApiStudyAnnotationOverview>>>;
  cellCount: Maybe<Scalars['Int']>;
  cellOntologyIds: Maybe<Array<Maybe<Scalars['String']>>>;
  cellOntologyLabels: Maybe<Array<Maybe<Scalars['String']>>>;
  description: Maybe<Scalars['String']>;
  diseaseMeshIds: Maybe<Array<Maybe<Scalars['String']>>>;
  diseaseMeshLabels: Maybe<Array<Maybe<Scalars['String']>>>;
  externalWebsite: Maybe<Scalars['String']>;
  layers: Maybe<Array<Maybe<Scalars['String']>>>;
  organismLabel: Maybe<Scalars['String']>;
  organismTaxId: Maybe<Scalars['String']>;
  studyId: Maybe<Scalars['Int']>;
  studyName: Maybe<Scalars['String']>;
  studyType: Maybe<Scalars['String']>;
  tissueNcitIds: Maybe<Array<Maybe<Scalars['String']>>>;
  tissueNcitLabels: Maybe<Array<Maybe<Scalars['String']>>>;
};

/**
 * A condition to be used against `ApiStudiesBase` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type ApiStudiesBaseCondition = {
  /** Checks for equality with the object’s `annotations` field. */
  annotations: InputMaybe<Array<InputMaybe<ApiStudyAnnotationOverviewInput>>>;
  /** Checks for equality with the object’s `cellCount` field. */
  cellCount: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `cellOntologyIds` field. */
  cellOntologyIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `cellOntologyLabels` field. */
  cellOntologyLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `description` field. */
  description: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `diseaseMeshIds` field. */
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `diseaseMeshLabels` field. */
  diseaseMeshLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `layers` field. */
  layers: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `organismLabel` field. */
  organismLabel: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyName` field. */
  studyName: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyType` field. */
  studyType: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `tissueNcitIds` field. */
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `tissueNcitLabels` field. */
  tissueNcitLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** A filter to be used against `ApiStudiesBase` object types. All fields are combined with a logical ‘and.’ */
export type ApiStudiesBaseFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ApiStudiesBaseFilter>>;
  /** Filter by the object’s `cellCount` field. */
  cellCount: InputMaybe<IntFilter>;
  /** Filter by the object’s `cellOntologyIds` field. */
  cellOntologyIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `cellOntologyLabels` field. */
  cellOntologyLabels: InputMaybe<StringListFilter>;
  /** Filter by the object’s `description` field. */
  description: InputMaybe<StringFilter>;
  /** Filter by the object’s `diseaseMeshIds` field. */
  diseaseMeshIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `diseaseMeshLabels` field. */
  diseaseMeshLabels: InputMaybe<StringListFilter>;
  /** Filter by the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<StringFilter>;
  /** Filter by the object’s `layers` field. */
  layers: InputMaybe<StringListFilter>;
  /** Negates the expression. */
  not: InputMaybe<ApiStudiesBaseFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ApiStudiesBaseFilter>>;
  /** Filter by the object’s `organismLabel` field. */
  organismLabel: InputMaybe<StringFilter>;
  /** Filter by the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<StringFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyName` field. */
  studyName: InputMaybe<StringFilter>;
  /** Filter by the object’s `studyType` field. */
  studyType: InputMaybe<StringFilter>;
  /** Filter by the object’s `tissueNcitIds` field. */
  tissueNcitIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `tissueNcitLabels` field. */
  tissueNcitLabels: InputMaybe<StringListFilter>;
};

/** An input for mutations affecting `ApiStudiesBase` */
export type ApiStudiesBaseInput = {
  annotations: InputMaybe<Array<InputMaybe<ApiStudyAnnotationOverviewInput>>>;
  cellCount: InputMaybe<Scalars['Int']>;
  cellOntologyIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  cellOntologyLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  description: InputMaybe<Scalars['String']>;
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  diseaseMeshLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  externalWebsite: InputMaybe<Scalars['String']>;
  layers: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  organismLabel: InputMaybe<Scalars['String']>;
  organismTaxId: InputMaybe<Scalars['String']>;
  studyId: InputMaybe<Scalars['Int']>;
  studyName: InputMaybe<Scalars['String']>;
  studyType: InputMaybe<Scalars['String']>;
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  tissueNcitLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** A connection to a list of `ApiStudiesBase` values. */
export type ApiStudiesBasesConnection = {
  __typename?: 'ApiStudiesBasesConnection';
  /** A list of edges which contains the `ApiStudiesBase` and cursor to aid in pagination. */
  edges: Array<ApiStudiesBasesEdge>;
  /** A list of `ApiStudiesBase` objects. */
  nodes: Array<Maybe<ApiStudiesBase>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ApiStudiesBase` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ApiStudiesBase` edge in the connection. */
export type ApiStudiesBasesEdge = {
  __typename?: 'ApiStudiesBasesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ApiStudiesBase` at the end of the edge. */
  node: Maybe<ApiStudiesBase>;
};

/** Methods to use when ordering `ApiStudiesBase`. */
export enum ApiStudiesBasesOrderBy {
  AnnotationsAsc = 'ANNOTATIONS_ASC',
  AnnotationsDesc = 'ANNOTATIONS_DESC',
  CellCountAsc = 'CELL_COUNT_ASC',
  CellCountDesc = 'CELL_COUNT_DESC',
  CellOntologyIdsAsc = 'CELL_ONTOLOGY_IDS_ASC',
  CellOntologyIdsDesc = 'CELL_ONTOLOGY_IDS_DESC',
  CellOntologyLabelsAsc = 'CELL_ONTOLOGY_LABELS_ASC',
  CellOntologyLabelsDesc = 'CELL_ONTOLOGY_LABELS_DESC',
  DescriptionAsc = 'DESCRIPTION_ASC',
  DescriptionDesc = 'DESCRIPTION_DESC',
  DiseaseMeshIdsAsc = 'DISEASE_MESH_IDS_ASC',
  DiseaseMeshIdsDesc = 'DISEASE_MESH_IDS_DESC',
  DiseaseMeshLabelsAsc = 'DISEASE_MESH_LABELS_ASC',
  DiseaseMeshLabelsDesc = 'DISEASE_MESH_LABELS_DESC',
  ExternalWebsiteAsc = 'EXTERNAL_WEBSITE_ASC',
  ExternalWebsiteDesc = 'EXTERNAL_WEBSITE_DESC',
  LayersAsc = 'LAYERS_ASC',
  LayersDesc = 'LAYERS_DESC',
  Natural = 'NATURAL',
  OrganismLabelAsc = 'ORGANISM_LABEL_ASC',
  OrganismLabelDesc = 'ORGANISM_LABEL_DESC',
  OrganismTaxIdAsc = 'ORGANISM_TAX_ID_ASC',
  OrganismTaxIdDesc = 'ORGANISM_TAX_ID_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudyNameAsc = 'STUDY_NAME_ASC',
  StudyNameDesc = 'STUDY_NAME_DESC',
  StudyTypeAsc = 'STUDY_TYPE_ASC',
  StudyTypeDesc = 'STUDY_TYPE_DESC',
  TissueNcitIdsAsc = 'TISSUE_NCIT_IDS_ASC',
  TissueNcitIdsDesc = 'TISSUE_NCIT_IDS_DESC',
  TissueNcitLabelsAsc = 'TISSUE_NCIT_LABELS_ASC',
  TissueNcitLabelsDesc = 'TISSUE_NCIT_LABELS_DESC'
}

export type ApiStudiesBulkRna = {
  __typename?: 'ApiStudiesBulkRna';
  annotations: Maybe<Array<Maybe<ApiStudyAnnotationOverview>>>;
  cellCount: Maybe<Scalars['Int']>;
  cellOntologyIds: Maybe<Array<Maybe<Scalars['String']>>>;
  cellOntologyLabels: Maybe<Array<Maybe<Scalars['String']>>>;
  description: Maybe<Scalars['String']>;
  diseaseMeshIds: Maybe<Array<Maybe<Scalars['String']>>>;
  diseaseMeshLabels: Maybe<Array<Maybe<Scalars['String']>>>;
  externalWebsite: Maybe<Scalars['String']>;
  layers: Maybe<Array<Maybe<Scalars['String']>>>;
  organismLabel: Maybe<Scalars['String']>;
  organismTaxId: Maybe<Scalars['String']>;
  studyId: Maybe<Scalars['Int']>;
  studyName: Maybe<Scalars['String']>;
  studyType: Maybe<Scalars['String']>;
  tissueNcitIds: Maybe<Array<Maybe<Scalars['String']>>>;
  tissueNcitLabels: Maybe<Array<Maybe<Scalars['String']>>>;
};

/**
 * A condition to be used against `ApiStudiesBulkRna` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type ApiStudiesBulkRnaCondition = {
  /** Checks for equality with the object’s `annotations` field. */
  annotations: InputMaybe<Array<InputMaybe<ApiStudyAnnotationOverviewInput>>>;
  /** Checks for equality with the object’s `cellCount` field. */
  cellCount: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `cellOntologyIds` field. */
  cellOntologyIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `cellOntologyLabels` field. */
  cellOntologyLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `description` field. */
  description: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `diseaseMeshIds` field. */
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `diseaseMeshLabels` field. */
  diseaseMeshLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `layers` field. */
  layers: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `organismLabel` field. */
  organismLabel: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyName` field. */
  studyName: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyType` field. */
  studyType: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `tissueNcitIds` field. */
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `tissueNcitLabels` field. */
  tissueNcitLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** A filter to be used against `ApiStudiesBulkRna` object types. All fields are combined with a logical ‘and.’ */
export type ApiStudiesBulkRnaFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ApiStudiesBulkRnaFilter>>;
  /** Filter by the object’s `cellCount` field. */
  cellCount: InputMaybe<IntFilter>;
  /** Filter by the object’s `cellOntologyIds` field. */
  cellOntologyIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `cellOntologyLabels` field. */
  cellOntologyLabels: InputMaybe<StringListFilter>;
  /** Filter by the object’s `description` field. */
  description: InputMaybe<StringFilter>;
  /** Filter by the object’s `diseaseMeshIds` field. */
  diseaseMeshIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `diseaseMeshLabels` field. */
  diseaseMeshLabels: InputMaybe<StringListFilter>;
  /** Filter by the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<StringFilter>;
  /** Filter by the object’s `layers` field. */
  layers: InputMaybe<StringListFilter>;
  /** Negates the expression. */
  not: InputMaybe<ApiStudiesBulkRnaFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ApiStudiesBulkRnaFilter>>;
  /** Filter by the object’s `organismLabel` field. */
  organismLabel: InputMaybe<StringFilter>;
  /** Filter by the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<StringFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyName` field. */
  studyName: InputMaybe<StringFilter>;
  /** Filter by the object’s `studyType` field. */
  studyType: InputMaybe<StringFilter>;
  /** Filter by the object’s `tissueNcitIds` field. */
  tissueNcitIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `tissueNcitLabels` field. */
  tissueNcitLabels: InputMaybe<StringListFilter>;
};

/** An input for mutations affecting `ApiStudiesBulkRna` */
export type ApiStudiesBulkRnaInput = {
  annotations: InputMaybe<Array<InputMaybe<ApiStudyAnnotationOverviewInput>>>;
  cellCount: InputMaybe<Scalars['Int']>;
  cellOntologyIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  cellOntologyLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  description: InputMaybe<Scalars['String']>;
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  diseaseMeshLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  externalWebsite: InputMaybe<Scalars['String']>;
  layers: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  organismLabel: InputMaybe<Scalars['String']>;
  organismTaxId: InputMaybe<Scalars['String']>;
  studyId: InputMaybe<Scalars['Int']>;
  studyName: InputMaybe<Scalars['String']>;
  studyType: InputMaybe<Scalars['String']>;
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  tissueNcitLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** A connection to a list of `ApiStudiesBulkRna` values. */
export type ApiStudiesBulkRnasConnection = {
  __typename?: 'ApiStudiesBulkRnasConnection';
  /** A list of edges which contains the `ApiStudiesBulkRna` and cursor to aid in pagination. */
  edges: Array<ApiStudiesBulkRnasEdge>;
  /** A list of `ApiStudiesBulkRna` objects. */
  nodes: Array<Maybe<ApiStudiesBulkRna>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ApiStudiesBulkRna` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ApiStudiesBulkRna` edge in the connection. */
export type ApiStudiesBulkRnasEdge = {
  __typename?: 'ApiStudiesBulkRnasEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ApiStudiesBulkRna` at the end of the edge. */
  node: Maybe<ApiStudiesBulkRna>;
};

/** Methods to use when ordering `ApiStudiesBulkRna`. */
export enum ApiStudiesBulkRnasOrderBy {
  AnnotationsAsc = 'ANNOTATIONS_ASC',
  AnnotationsDesc = 'ANNOTATIONS_DESC',
  CellCountAsc = 'CELL_COUNT_ASC',
  CellCountDesc = 'CELL_COUNT_DESC',
  CellOntologyIdsAsc = 'CELL_ONTOLOGY_IDS_ASC',
  CellOntologyIdsDesc = 'CELL_ONTOLOGY_IDS_DESC',
  CellOntologyLabelsAsc = 'CELL_ONTOLOGY_LABELS_ASC',
  CellOntologyLabelsDesc = 'CELL_ONTOLOGY_LABELS_DESC',
  DescriptionAsc = 'DESCRIPTION_ASC',
  DescriptionDesc = 'DESCRIPTION_DESC',
  DiseaseMeshIdsAsc = 'DISEASE_MESH_IDS_ASC',
  DiseaseMeshIdsDesc = 'DISEASE_MESH_IDS_DESC',
  DiseaseMeshLabelsAsc = 'DISEASE_MESH_LABELS_ASC',
  DiseaseMeshLabelsDesc = 'DISEASE_MESH_LABELS_DESC',
  ExternalWebsiteAsc = 'EXTERNAL_WEBSITE_ASC',
  ExternalWebsiteDesc = 'EXTERNAL_WEBSITE_DESC',
  LayersAsc = 'LAYERS_ASC',
  LayersDesc = 'LAYERS_DESC',
  Natural = 'NATURAL',
  OrganismLabelAsc = 'ORGANISM_LABEL_ASC',
  OrganismLabelDesc = 'ORGANISM_LABEL_DESC',
  OrganismTaxIdAsc = 'ORGANISM_TAX_ID_ASC',
  OrganismTaxIdDesc = 'ORGANISM_TAX_ID_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudyNameAsc = 'STUDY_NAME_ASC',
  StudyNameDesc = 'STUDY_NAME_DESC',
  StudyTypeAsc = 'STUDY_TYPE_ASC',
  StudyTypeDesc = 'STUDY_TYPE_DESC',
  TissueNcitIdsAsc = 'TISSUE_NCIT_IDS_ASC',
  TissueNcitIdsDesc = 'TISSUE_NCIT_IDS_DESC',
  TissueNcitLabelsAsc = 'TISSUE_NCIT_LABELS_ASC',
  TissueNcitLabelsDesc = 'TISSUE_NCIT_LABELS_DESC'
}

export type ApiStudiesSingleCell = {
  __typename?: 'ApiStudiesSingleCell';
  annotations: Maybe<Array<Maybe<ApiStudyAnnotationOverview>>>;
  cellCount: Maybe<Scalars['Int']>;
  cellOntologyIds: Maybe<Array<Maybe<Scalars['String']>>>;
  cellOntologyLabels: Maybe<Array<Maybe<Scalars['String']>>>;
  description: Maybe<Scalars['String']>;
  diseaseMeshIds: Maybe<Array<Maybe<Scalars['String']>>>;
  diseaseMeshLabels: Maybe<Array<Maybe<Scalars['String']>>>;
  externalWebsite: Maybe<Scalars['String']>;
  layers: Maybe<Array<Maybe<Scalars['String']>>>;
  organismLabel: Maybe<Scalars['String']>;
  organismTaxId: Maybe<Scalars['String']>;
  studyId: Maybe<Scalars['Int']>;
  studyName: Maybe<Scalars['String']>;
  studyType: Maybe<Scalars['String']>;
  tissueNcitIds: Maybe<Array<Maybe<Scalars['String']>>>;
  tissueNcitLabels: Maybe<Array<Maybe<Scalars['String']>>>;
};

/**
 * A condition to be used against `ApiStudiesSingleCell` object types. All fields
 * are tested for equality and combined with a logical ‘and.’
 */
export type ApiStudiesSingleCellCondition = {
  /** Checks for equality with the object’s `annotations` field. */
  annotations: InputMaybe<Array<InputMaybe<ApiStudyAnnotationOverviewInput>>>;
  /** Checks for equality with the object’s `cellCount` field. */
  cellCount: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `cellOntologyIds` field. */
  cellOntologyIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `cellOntologyLabels` field. */
  cellOntologyLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `description` field. */
  description: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `diseaseMeshIds` field. */
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `diseaseMeshLabels` field. */
  diseaseMeshLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `layers` field. */
  layers: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `organismLabel` field. */
  organismLabel: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyName` field. */
  studyName: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyType` field. */
  studyType: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `tissueNcitIds` field. */
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `tissueNcitLabels` field. */
  tissueNcitLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** A filter to be used against `ApiStudiesSingleCell` object types. All fields are combined with a logical ‘and.’ */
export type ApiStudiesSingleCellFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ApiStudiesSingleCellFilter>>;
  /** Filter by the object’s `cellCount` field. */
  cellCount: InputMaybe<IntFilter>;
  /** Filter by the object’s `cellOntologyIds` field. */
  cellOntologyIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `cellOntologyLabels` field. */
  cellOntologyLabels: InputMaybe<StringListFilter>;
  /** Filter by the object’s `description` field. */
  description: InputMaybe<StringFilter>;
  /** Filter by the object’s `diseaseMeshIds` field. */
  diseaseMeshIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `diseaseMeshLabels` field. */
  diseaseMeshLabels: InputMaybe<StringListFilter>;
  /** Filter by the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<StringFilter>;
  /** Filter by the object’s `layers` field. */
  layers: InputMaybe<StringListFilter>;
  /** Negates the expression. */
  not: InputMaybe<ApiStudiesSingleCellFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ApiStudiesSingleCellFilter>>;
  /** Filter by the object’s `organismLabel` field. */
  organismLabel: InputMaybe<StringFilter>;
  /** Filter by the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<StringFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyName` field. */
  studyName: InputMaybe<StringFilter>;
  /** Filter by the object’s `studyType` field. */
  studyType: InputMaybe<StringFilter>;
  /** Filter by the object’s `tissueNcitIds` field. */
  tissueNcitIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `tissueNcitLabels` field. */
  tissueNcitLabels: InputMaybe<StringListFilter>;
};

/** An input for mutations affecting `ApiStudiesSingleCell` */
export type ApiStudiesSingleCellInput = {
  annotations: InputMaybe<Array<InputMaybe<ApiStudyAnnotationOverviewInput>>>;
  cellCount: InputMaybe<Scalars['Int']>;
  cellOntologyIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  cellOntologyLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  description: InputMaybe<Scalars['String']>;
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  diseaseMeshLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  externalWebsite: InputMaybe<Scalars['String']>;
  layers: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  organismLabel: InputMaybe<Scalars['String']>;
  organismTaxId: InputMaybe<Scalars['String']>;
  studyId: InputMaybe<Scalars['Int']>;
  studyName: InputMaybe<Scalars['String']>;
  studyType: InputMaybe<Scalars['String']>;
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  tissueNcitLabels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** A connection to a list of `ApiStudiesSingleCell` values. */
export type ApiStudiesSingleCellsConnection = {
  __typename?: 'ApiStudiesSingleCellsConnection';
  /** A list of edges which contains the `ApiStudiesSingleCell` and cursor to aid in pagination. */
  edges: Array<ApiStudiesSingleCellsEdge>;
  /** A list of `ApiStudiesSingleCell` objects. */
  nodes: Array<Maybe<ApiStudiesSingleCell>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ApiStudiesSingleCell` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ApiStudiesSingleCell` edge in the connection. */
export type ApiStudiesSingleCellsEdge = {
  __typename?: 'ApiStudiesSingleCellsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ApiStudiesSingleCell` at the end of the edge. */
  node: Maybe<ApiStudiesSingleCell>;
};

/** Methods to use when ordering `ApiStudiesSingleCell`. */
export enum ApiStudiesSingleCellsOrderBy {
  AnnotationsAsc = 'ANNOTATIONS_ASC',
  AnnotationsDesc = 'ANNOTATIONS_DESC',
  CellCountAsc = 'CELL_COUNT_ASC',
  CellCountDesc = 'CELL_COUNT_DESC',
  CellOntologyIdsAsc = 'CELL_ONTOLOGY_IDS_ASC',
  CellOntologyIdsDesc = 'CELL_ONTOLOGY_IDS_DESC',
  CellOntologyLabelsAsc = 'CELL_ONTOLOGY_LABELS_ASC',
  CellOntologyLabelsDesc = 'CELL_ONTOLOGY_LABELS_DESC',
  DescriptionAsc = 'DESCRIPTION_ASC',
  DescriptionDesc = 'DESCRIPTION_DESC',
  DiseaseMeshIdsAsc = 'DISEASE_MESH_IDS_ASC',
  DiseaseMeshIdsDesc = 'DISEASE_MESH_IDS_DESC',
  DiseaseMeshLabelsAsc = 'DISEASE_MESH_LABELS_ASC',
  DiseaseMeshLabelsDesc = 'DISEASE_MESH_LABELS_DESC',
  ExternalWebsiteAsc = 'EXTERNAL_WEBSITE_ASC',
  ExternalWebsiteDesc = 'EXTERNAL_WEBSITE_DESC',
  LayersAsc = 'LAYERS_ASC',
  LayersDesc = 'LAYERS_DESC',
  Natural = 'NATURAL',
  OrganismLabelAsc = 'ORGANISM_LABEL_ASC',
  OrganismLabelDesc = 'ORGANISM_LABEL_DESC',
  OrganismTaxIdAsc = 'ORGANISM_TAX_ID_ASC',
  OrganismTaxIdDesc = 'ORGANISM_TAX_ID_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudyNameAsc = 'STUDY_NAME_ASC',
  StudyNameDesc = 'STUDY_NAME_DESC',
  StudyTypeAsc = 'STUDY_TYPE_ASC',
  StudyTypeDesc = 'STUDY_TYPE_DESC',
  TissueNcitIdsAsc = 'TISSUE_NCIT_IDS_ASC',
  TissueNcitIdsDesc = 'TISSUE_NCIT_IDS_DESC',
  TissueNcitLabelsAsc = 'TISSUE_NCIT_LABELS_ASC',
  TissueNcitLabelsDesc = 'TISSUE_NCIT_LABELS_DESC'
}

export type ApiStudyAnnotationOverview = {
  __typename?: 'ApiStudyAnnotationOverview';
  annotationDisplayValues: Maybe<Array<Maybe<Scalars['String']>>>;
  annotationGroupColumn: Maybe<Scalars['String']>;
  annotationGroupDisplay: Maybe<Scalars['String']>;
  annotationValues: Maybe<Array<Maybe<Scalars['String']>>>;
  studyId: Maybe<Scalars['Int']>;
};

/**
 * A condition to be used against `ApiStudyAnnotationOverview` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type ApiStudyAnnotationOverviewCondition = {
  /** Checks for equality with the object’s `annotationDisplayValues` field. */
  annotationDisplayValues: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `annotationGroupColumn` field. */
  annotationGroupColumn: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `annotationGroupDisplay` field. */
  annotationGroupDisplay: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `annotationValues` field. */
  annotationValues: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `ApiStudyAnnotationOverview` object types. All fields are combined with a logical ‘and.’ */
export type ApiStudyAnnotationOverviewFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ApiStudyAnnotationOverviewFilter>>;
  /** Filter by the object’s `annotationDisplayValues` field. */
  annotationDisplayValues: InputMaybe<StringListFilter>;
  /** Filter by the object’s `annotationGroupColumn` field. */
  annotationGroupColumn: InputMaybe<StringFilter>;
  /** Filter by the object’s `annotationGroupDisplay` field. */
  annotationGroupDisplay: InputMaybe<StringFilter>;
  /** Filter by the object’s `annotationValues` field. */
  annotationValues: InputMaybe<StringListFilter>;
  /** Negates the expression. */
  not: InputMaybe<ApiStudyAnnotationOverviewFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ApiStudyAnnotationOverviewFilter>>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `ApiStudyAnnotationOverview` */
export type ApiStudyAnnotationOverviewInput = {
  annotationDisplayValues: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  annotationGroupColumn: InputMaybe<Scalars['String']>;
  annotationGroupDisplay: InputMaybe<Scalars['String']>;
  annotationValues: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  studyId: InputMaybe<Scalars['Int']>;
};

/** A connection to a list of `ApiStudyAnnotationOverview` values. */
export type ApiStudyAnnotationOverviewsConnection = {
  __typename?: 'ApiStudyAnnotationOverviewsConnection';
  /** A list of edges which contains the `ApiStudyAnnotationOverview` and cursor to aid in pagination. */
  edges: Array<ApiStudyAnnotationOverviewsEdge>;
  /** A list of `ApiStudyAnnotationOverview` objects. */
  nodes: Array<Maybe<ApiStudyAnnotationOverview>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ApiStudyAnnotationOverview` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ApiStudyAnnotationOverview` edge in the connection. */
export type ApiStudyAnnotationOverviewsEdge = {
  __typename?: 'ApiStudyAnnotationOverviewsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ApiStudyAnnotationOverview` at the end of the edge. */
  node: Maybe<ApiStudyAnnotationOverview>;
};

/** Methods to use when ordering `ApiStudyAnnotationOverview`. */
export enum ApiStudyAnnotationOverviewsOrderBy {
  AnnotationDisplayValuesAsc = 'ANNOTATION_DISPLAY_VALUES_ASC',
  AnnotationDisplayValuesDesc = 'ANNOTATION_DISPLAY_VALUES_DESC',
  AnnotationGroupColumnAsc = 'ANNOTATION_GROUP_COLUMN_ASC',
  AnnotationGroupColumnDesc = 'ANNOTATION_GROUP_COLUMN_DESC',
  AnnotationGroupDisplayAsc = 'ANNOTATION_GROUP_DISPLAY_ASC',
  AnnotationGroupDisplayDesc = 'ANNOTATION_GROUP_DISPLAY_DESC',
  AnnotationValuesAsc = 'ANNOTATION_VALUES_ASC',
  AnnotationValuesDesc = 'ANNOTATION_VALUES_DESC',
  Natural = 'NATURAL',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

export type ApiStudyH5Download = {
  __typename?: 'ApiStudyH5Download';
  presignedUrl: Maybe<Scalars['String']>;
  studyId: Maybe<Scalars['Int']>;
};

/**
 * A condition to be used against `ApiStudyH5Download` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type ApiStudyH5DownloadCondition = {
  /** Checks for equality with the object’s `presignedUrl` field. */
  presignedUrl: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `ApiStudyH5Download` object types. All fields are combined with a logical ‘and.’ */
export type ApiStudyH5DownloadFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ApiStudyH5DownloadFilter>>;
  /** Negates the expression. */
  not: InputMaybe<ApiStudyH5DownloadFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ApiStudyH5DownloadFilter>>;
  /** Filter by the object’s `presignedUrl` field. */
  presignedUrl: InputMaybe<StringFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `ApiStudyH5Download` */
export type ApiStudyH5DownloadInput = {
  presignedUrl: InputMaybe<Scalars['String']>;
  studyId: InputMaybe<Scalars['Int']>;
};

/** A connection to a list of `ApiStudyH5Download` values. */
export type ApiStudyH5DownloadsConnection = {
  __typename?: 'ApiStudyH5DownloadsConnection';
  /** A list of edges which contains the `ApiStudyH5Download` and cursor to aid in pagination. */
  edges: Array<ApiStudyH5DownloadsEdge>;
  /** A list of `ApiStudyH5Download` objects. */
  nodes: Array<Maybe<ApiStudyH5Download>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ApiStudyH5Download` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ApiStudyH5Download` edge in the connection. */
export type ApiStudyH5DownloadsEdge = {
  __typename?: 'ApiStudyH5DownloadsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ApiStudyH5Download` at the end of the edge. */
  node: Maybe<ApiStudyH5Download>;
};

/** Methods to use when ordering `ApiStudyH5Download`. */
export enum ApiStudyH5DownloadsOrderBy {
  Natural = 'NATURAL',
  PresignedUrlAsc = 'PRESIGNED_URL_ASC',
  PresignedUrlDesc = 'PRESIGNED_URL_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

export type AutocompleteResult = {
  __typename?: 'AutocompleteResult';
  isSynonymOfPreferredTerm: Maybe<Scalars['String']>;
  label: Maybe<Scalars['String']>;
  labelHighlight: Maybe<Scalars['String']>;
  ontCode: Maybe<Scalars['String']>;
  ontology: Maybe<Scalars['String']>;
};

/** A filter to be used against `AutocompleteResult` object types. All fields are combined with a logical ‘and.’ */
export type AutocompleteResultFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<AutocompleteResultFilter>>;
  /** Filter by the object’s `isSynonymOfPreferredTerm` field. */
  isSynonymOfPreferredTerm: InputMaybe<StringFilter>;
  /** Filter by the object’s `label` field. */
  label: InputMaybe<StringFilter>;
  /** Filter by the object’s `labelHighlight` field. */
  labelHighlight: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<AutocompleteResultFilter>;
  /** Filter by the object’s `ontCode` field. */
  ontCode: InputMaybe<StringFilter>;
  /** Filter by the object’s `ontology` field. */
  ontology: InputMaybe<StringFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<AutocompleteResultFilter>>;
};

/** A connection to a list of `AutocompleteResult` values. */
export type AutocompleteResultsConnection = {
  __typename?: 'AutocompleteResultsConnection';
  /** A list of edges which contains the `AutocompleteResult` and cursor to aid in pagination. */
  edges: Array<AutocompleteResultsEdge>;
  /** A list of `AutocompleteResult` objects. */
  nodes: Array<Maybe<AutocompleteResult>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `AutocompleteResult` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `AutocompleteResult` edge in the connection. */
export type AutocompleteResultsEdge = {
  __typename?: 'AutocompleteResultsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `AutocompleteResult` at the end of the edge. */
  node: Maybe<AutocompleteResult>;
};

/** A filter to be used against BigInt fields. All fields are combined with a logical ‘and.’ */
export type BigIntFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<Scalars['BigInt']>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<Scalars['BigInt']>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<Scalars['BigInt']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<Scalars['BigInt']>;
  /** Included in the specified list. */
  in: InputMaybe<Array<Scalars['BigInt']>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<Scalars['BigInt']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<Scalars['BigInt']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<Scalars['BigInt']>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<Scalars['BigInt']>;
  /** Not included in the specified list. */
  notIn: InputMaybe<Array<Scalars['BigInt']>>;
};

/** A filter to be used against Boolean fields. All fields are combined with a logical ‘and.’ */
export type BooleanFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<Scalars['Boolean']>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<Scalars['Boolean']>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<Scalars['Boolean']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<Scalars['Boolean']>;
  /** Included in the specified list. */
  in: InputMaybe<Array<Scalars['Boolean']>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<Scalars['Boolean']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<Scalars['Boolean']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<Scalars['Boolean']>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<Scalars['Boolean']>;
  /** Not included in the specified list. */
  notIn: InputMaybe<Array<Scalars['Boolean']>>;
};

export type BoxplotValue = {
  __typename?: 'BoxplotValue';
  median: Maybe<Scalars['Float']>;
  n: Maybe<Scalars['Int']>;
  outliers: Maybe<Array<Maybe<Scalars['Float']>>>;
  q1: Maybe<Scalars['Float']>;
  q1Whisker: Maybe<Scalars['Float']>;
  q3: Maybe<Scalars['Float']>;
  q3Whisker: Maybe<Scalars['Float']>;
};

/** A filter to be used against `BoxplotValue` object types. All fields are combined with a logical ‘and.’ */
export type BoxplotValueFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<BoxplotValueFilter>>;
  /** Filter by the object’s `median` field. */
  median: InputMaybe<FloatFilter>;
  /** Filter by the object’s `n` field. */
  n: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not: InputMaybe<BoxplotValueFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<BoxplotValueFilter>>;
  /** Filter by the object’s `outliers` field. */
  outliers: InputMaybe<FloatListFilter>;
  /** Filter by the object’s `q1` field. */
  q1: InputMaybe<FloatFilter>;
  /** Filter by the object’s `q1Whisker` field. */
  q1Whisker: InputMaybe<FloatFilter>;
  /** Filter by the object’s `q3` field. */
  q3: InputMaybe<FloatFilter>;
  /** Filter by the object’s `q3Whisker` field. */
  q3Whisker: InputMaybe<FloatFilter>;
};

export type Concept = Node & {
  __typename?: 'Concept';
  /** Reads and enables pagination through a set of `Concept`. */
  allChildren: ConceptsConnection;
  /** Reads and enables pagination through a set of `Concept`. */
  allChildrenList: Maybe<Array<Maybe<Concept>>>;
  /** Reads and enables pagination through a set of `Concept`. */
  allParents: ConceptsConnection;
  /** Reads and enables pagination through a set of `Concept`. */
  allParentsList: Maybe<Array<Maybe<Concept>>>;
  allParentsPaths: ConceptAllParentsPathsConnection;
  allParentsPathsList: Maybe<Array<Maybe<ConceptAllParentsPathsRecord>>>;
  /** Reads and enables pagination through a set of `Concept`. */
  childrenDepthLimit: ConceptsConnection;
  /** Reads and enables pagination through a set of `Concept`. */
  childrenDepthLimitList: Maybe<Array<Maybe<Concept>>>;
  childrenPathsDepthLimit: ConceptChildrenPathsDepthLimitConnection;
  childrenPathsDepthLimitList: Maybe<Array<Maybe<ConceptChildrenPathsDepthLimitRecord>>>;
  cid: Scalars['Int'];
  /** Reads and enables pagination through a set of `ConceptHierarchy`. */
  conceptHierarchiesByCid: ConceptHierarchiesConnection;
  /** Reads and enables pagination through a set of `ConceptHierarchy`. */
  conceptHierarchiesByCidList: Array<ConceptHierarchy>;
  /** Reads and enables pagination through a set of `ConceptHierarchy`. */
  conceptHierarchiesByParentCid: ConceptHierarchiesConnection;
  /** Reads and enables pagination through a set of `ConceptHierarchy`. */
  conceptHierarchiesByParentCidList: Array<ConceptHierarchy>;
  /** Reads and enables pagination through a set of `ConceptSynonym`. */
  conceptSynonymsByCid: ConceptSynonymsConnection;
  /** Reads and enables pagination through a set of `ConceptSynonym`. */
  conceptSynonymsByCidList: Array<ConceptSynonym>;
  /** Reads and enables pagination through a set of `Concept`. */
  directChildren: ConceptsConnection;
  /** Reads and enables pagination through a set of `Concept`. */
  directChildrenList: Maybe<Array<Maybe<Concept>>>;
  /** Reads and enables pagination through a set of `Concept`. */
  directParents: ConceptsConnection;
  /** Reads and enables pagination through a set of `Concept`. */
  directParentsList: Maybe<Array<Maybe<Concept>>>;
  label: Maybe<Scalars['String']>;
  labelTsvector: Maybe<Scalars['String']>;
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  ontCode: Maybe<Scalars['String']>;
  ontid: Maybe<Scalars['Int']>;
  /** Reads a single `Ontology` that is related to this `Concept`. */
  ontologyByOntid: Maybe<Ontology>;
  /** Reads and enables pagination through a set of `Concept`. */
  testPotentialChildren: ConceptsConnection;
  /** Reads and enables pagination through a set of `Concept`. */
  testPotentialChildrenList: Maybe<Array<Maybe<Concept>>>;
};


export type ConceptAllChildrenArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptAllChildrenListArgs = {
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptAllParentsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptAllParentsListArgs = {
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptAllParentsPathsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<ConceptAllParentsPathsRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptAllParentsPathsListArgs = {
  filter: InputMaybe<ConceptAllParentsPathsRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptChildrenDepthLimitArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  depthLimit: InputMaybe<Scalars['Int']>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptChildrenDepthLimitListArgs = {
  depthLimit: InputMaybe<Scalars['Int']>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptChildrenPathsDepthLimitArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  depthLimit: InputMaybe<Scalars['Int']>;
  filter: InputMaybe<ConceptChildrenPathsDepthLimitRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptChildrenPathsDepthLimitListArgs = {
  depthLimit: InputMaybe<Scalars['Int']>;
  filter: InputMaybe<ConceptChildrenPathsDepthLimitRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptConceptHierarchiesByCidArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ConceptHierarchyCondition>;
  filter: InputMaybe<ConceptHierarchyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
};


export type ConceptConceptHierarchiesByCidListArgs = {
  condition: InputMaybe<ConceptHierarchyCondition>;
  filter: InputMaybe<ConceptHierarchyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
};


export type ConceptConceptHierarchiesByParentCidArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ConceptHierarchyCondition>;
  filter: InputMaybe<ConceptHierarchyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
};


export type ConceptConceptHierarchiesByParentCidListArgs = {
  condition: InputMaybe<ConceptHierarchyCondition>;
  filter: InputMaybe<ConceptHierarchyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
};


export type ConceptConceptSynonymsByCidArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ConceptSynonymCondition>;
  filter: InputMaybe<ConceptSynonymFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptSynonymsOrderBy>>;
};


export type ConceptConceptSynonymsByCidListArgs = {
  condition: InputMaybe<ConceptSynonymCondition>;
  filter: InputMaybe<ConceptSynonymFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ConceptSynonymsOrderBy>>;
};


export type ConceptDirectChildrenArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptDirectChildrenListArgs = {
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptDirectParentsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptDirectParentsListArgs = {
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptTestPotentialChildrenArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  potentialChildrenCodes: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};


export type ConceptTestPotentialChildrenListArgs = {
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  potentialChildrenCodes: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** All input for the `conceptAllChildrenPaths` mutation. */
export type ConceptAllChildrenPathsInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  conceptInp: InputMaybe<ConceptInput>;
};

/** The output of our `conceptAllChildrenPaths` mutation. */
export type ConceptAllChildrenPathsPayload = {
  __typename?: 'ConceptAllChildrenPathsPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  results: Maybe<Array<Maybe<ConceptAllChildrenPathsRecord>>>;
};

/** The return type of our `conceptAllChildrenPaths` mutation. */
export type ConceptAllChildrenPathsRecord = {
  __typename?: 'ConceptAllChildrenPathsRecord';
  cidPath: Maybe<Array<Maybe<Scalars['Int']>>>;
  ontCodePath: Maybe<Array<Maybe<Scalars['String']>>>;
};

/** A `ConceptAllParentsPathsRecord` edge in the connection. */
export type ConceptAllParentsPathEdge = {
  __typename?: 'ConceptAllParentsPathEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ConceptAllParentsPathsRecord` at the end of the edge. */
  node: Maybe<ConceptAllParentsPathsRecord>;
};

/** A connection to a list of `ConceptAllParentsPathsRecord` values. */
export type ConceptAllParentsPathsConnection = {
  __typename?: 'ConceptAllParentsPathsConnection';
  /** A list of edges which contains the `ConceptAllParentsPathsRecord` and cursor to aid in pagination. */
  edges: Array<ConceptAllParentsPathEdge>;
  /** A list of `ConceptAllParentsPathsRecord` objects. */
  nodes: Array<Maybe<ConceptAllParentsPathsRecord>>;
  /** The count of *all* `ConceptAllParentsPathsRecord` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** The return type of our `allParentsPaths` query. */
export type ConceptAllParentsPathsRecord = {
  __typename?: 'ConceptAllParentsPathsRecord';
  cidPath: Maybe<Array<Maybe<Scalars['Int']>>>;
  ontCodePath: Maybe<Array<Maybe<Scalars['String']>>>;
};

/** A filter to be used against `ConceptAllParentsPathsRecord` object types. All fields are combined with a logical ‘and.’ */
export type ConceptAllParentsPathsRecordFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ConceptAllParentsPathsRecordFilter>>;
  /** Filter by the object’s `cidPath` field. */
  cidPath: InputMaybe<IntListFilter>;
  /** Negates the expression. */
  not: InputMaybe<ConceptAllParentsPathsRecordFilter>;
  /** Filter by the object’s `ontCodePath` field. */
  ontCodePath: InputMaybe<StringListFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ConceptAllParentsPathsRecordFilter>>;
};

/** A connection to a list of `ConceptChildrenPathsDepthLimitRecord` values. */
export type ConceptChildrenPathsDepthLimitConnection = {
  __typename?: 'ConceptChildrenPathsDepthLimitConnection';
  /** A list of edges which contains the `ConceptChildrenPathsDepthLimitRecord` and cursor to aid in pagination. */
  edges: Array<ConceptChildrenPathsDepthLimitEdge>;
  /** A list of `ConceptChildrenPathsDepthLimitRecord` objects. */
  nodes: Array<Maybe<ConceptChildrenPathsDepthLimitRecord>>;
  /** The count of *all* `ConceptChildrenPathsDepthLimitRecord` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ConceptChildrenPathsDepthLimitRecord` edge in the connection. */
export type ConceptChildrenPathsDepthLimitEdge = {
  __typename?: 'ConceptChildrenPathsDepthLimitEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ConceptChildrenPathsDepthLimitRecord` at the end of the edge. */
  node: Maybe<ConceptChildrenPathsDepthLimitRecord>;
};

/** The return type of our `childrenPathsDepthLimit` query. */
export type ConceptChildrenPathsDepthLimitRecord = {
  __typename?: 'ConceptChildrenPathsDepthLimitRecord';
  cidPath: Maybe<Array<Maybe<Scalars['Int']>>>;
  ontCodePath: Maybe<Array<Maybe<Scalars['String']>>>;
};

/** A filter to be used against `ConceptChildrenPathsDepthLimitRecord` object types. All fields are combined with a logical ‘and.’ */
export type ConceptChildrenPathsDepthLimitRecordFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ConceptChildrenPathsDepthLimitRecordFilter>>;
  /** Filter by the object’s `cidPath` field. */
  cidPath: InputMaybe<IntListFilter>;
  /** Negates the expression. */
  not: InputMaybe<ConceptChildrenPathsDepthLimitRecordFilter>;
  /** Filter by the object’s `ontCodePath` field. */
  ontCodePath: InputMaybe<StringListFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ConceptChildrenPathsDepthLimitRecordFilter>>;
};

/** A condition to be used against `Concept` object types. All fields are tested for equality and combined with a logical ‘and.’ */
export type ConceptCondition = {
  /** Checks for equality with the object’s `cid` field. */
  cid: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `label` field. */
  label: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `labelTsvector` field. */
  labelTsvector: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `ontCode` field. */
  ontCode: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `ontid` field. */
  ontid: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `Concept` object types. All fields are combined with a logical ‘and.’ */
export type ConceptFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ConceptFilter>>;
  /** Filter by the object’s `cid` field. */
  cid: InputMaybe<IntFilter>;
  /** Filter by the object’s `label` field. */
  label: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<ConceptFilter>;
  /** Filter by the object’s `ontCode` field. */
  ontCode: InputMaybe<StringFilter>;
  /** Filter by the object’s `ontid` field. */
  ontid: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ConceptFilter>>;
};

/** A connection to a list of `ConceptHierarchy` values. */
export type ConceptHierarchiesConnection = {
  __typename?: 'ConceptHierarchiesConnection';
  /** A list of edges which contains the `ConceptHierarchy` and cursor to aid in pagination. */
  edges: Array<ConceptHierarchiesEdge>;
  /** A list of `ConceptHierarchy` objects. */
  nodes: Array<Maybe<ConceptHierarchy>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ConceptHierarchy` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ConceptHierarchy` edge in the connection. */
export type ConceptHierarchiesEdge = {
  __typename?: 'ConceptHierarchiesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ConceptHierarchy` at the end of the edge. */
  node: Maybe<ConceptHierarchy>;
};

/** Methods to use when ordering `ConceptHierarchy`. */
export enum ConceptHierarchiesOrderBy {
  CidAsc = 'CID_ASC',
  CidDesc = 'CID_DESC',
  Natural = 'NATURAL',
  ParentCidAsc = 'PARENT_CID_ASC',
  ParentCidDesc = 'PARENT_CID_DESC'
}

export type ConceptHierarchy = {
  __typename?: 'ConceptHierarchy';
  cid: Maybe<Scalars['Int']>;
  /** Reads a single `Concept` that is related to this `ConceptHierarchy`. */
  conceptByCid: Maybe<Concept>;
  /** Reads a single `Concept` that is related to this `ConceptHierarchy`. */
  conceptByParentCid: Maybe<Concept>;
  parentCid: Maybe<Scalars['Int']>;
};

/**
 * A condition to be used against `ConceptHierarchy` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type ConceptHierarchyCondition = {
  /** Checks for equality with the object’s `cid` field. */
  cid: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `parentCid` field. */
  parentCid: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `ConceptHierarchy` object types. All fields are combined with a logical ‘and.’ */
export type ConceptHierarchyFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ConceptHierarchyFilter>>;
  /** Filter by the object’s `cid` field. */
  cid: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not: InputMaybe<ConceptHierarchyFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ConceptHierarchyFilter>>;
  /** Filter by the object’s `parentCid` field. */
  parentCid: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `ConceptHierarchy` */
export type ConceptHierarchyInput = {
  cid: InputMaybe<Scalars['Int']>;
  parentCid: InputMaybe<Scalars['Int']>;
};

/** A `ConceptHierarchyMinimumTreesParentsListsRecord` edge in the connection. */
export type ConceptHierarchyMinimumTreesParentsListEdge = {
  __typename?: 'ConceptHierarchyMinimumTreesParentsListEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ConceptHierarchyMinimumTreesParentsListsRecord` at the end of the edge. */
  node: Maybe<ConceptHierarchyMinimumTreesParentsListsRecord>;
};

/** A connection to a list of `ConceptHierarchyMinimumTreesParentsListsRecord` values. */
export type ConceptHierarchyMinimumTreesParentsListsConnection = {
  __typename?: 'ConceptHierarchyMinimumTreesParentsListsConnection';
  /** A list of edges which contains the `ConceptHierarchyMinimumTreesParentsListsRecord` and cursor to aid in pagination. */
  edges: Array<ConceptHierarchyMinimumTreesParentsListEdge>;
  /** A list of `ConceptHierarchyMinimumTreesParentsListsRecord` objects. */
  nodes: Array<Maybe<ConceptHierarchyMinimumTreesParentsListsRecord>>;
  /** The count of *all* `ConceptHierarchyMinimumTreesParentsListsRecord` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** The return type of our `conceptHierarchyMinimumTreesParentsLists` query. */
export type ConceptHierarchyMinimumTreesParentsListsRecord = {
  __typename?: 'ConceptHierarchyMinimumTreesParentsListsRecord';
  cid: Maybe<Scalars['Int']>;
  label: Maybe<Scalars['String']>;
  ontCode: Maybe<Scalars['String']>;
  parentCids: Maybe<Array<Maybe<Scalars['Int']>>>;
  parentOntCodePath: Maybe<Array<Maybe<Scalars['String']>>>;
};

/** A filter to be used against `ConceptHierarchyMinimumTreesParentsListsRecord` object types. All fields are combined with a logical ‘and.’ */
export type ConceptHierarchyMinimumTreesParentsListsRecordFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ConceptHierarchyMinimumTreesParentsListsRecordFilter>>;
  /** Filter by the object’s `cid` field. */
  cid: InputMaybe<IntFilter>;
  /** Filter by the object’s `label` field. */
  label: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<ConceptHierarchyMinimumTreesParentsListsRecordFilter>;
  /** Filter by the object’s `ontCode` field. */
  ontCode: InputMaybe<StringFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ConceptHierarchyMinimumTreesParentsListsRecordFilter>>;
  /** Filter by the object’s `parentCids` field. */
  parentCids: InputMaybe<IntListFilter>;
  /** Filter by the object’s `parentOntCodePath` field. */
  parentOntCodePath: InputMaybe<StringListFilter>;
};

/** An input for mutations affecting `Concept` */
export type ConceptInput = {
  cid: InputMaybe<Scalars['Int']>;
  label: InputMaybe<Scalars['String']>;
  labelTsvector: InputMaybe<Scalars['String']>;
  ontCode: InputMaybe<Scalars['String']>;
  ontid: InputMaybe<Scalars['Int']>;
};

/** Represents an update to a `Concept`. Fields that are set will be updated. */
export type ConceptPatch = {
  cid: InputMaybe<Scalars['Int']>;
  label: InputMaybe<Scalars['String']>;
  labelTsvector: InputMaybe<Scalars['String']>;
  ontCode: InputMaybe<Scalars['String']>;
  ontid: InputMaybe<Scalars['Int']>;
};

export type ConceptPath = {
  __typename?: 'ConceptPath';
  cid: Maybe<Scalars['Int']>;
  parentCids: Maybe<Array<Maybe<Scalars['Int']>>>;
};

export type ConceptSynonym = {
  __typename?: 'ConceptSynonym';
  cid: Scalars['Int'];
  /** Reads a single `Concept` that is related to this `ConceptSynonym`. */
  conceptByCid: Maybe<Concept>;
  synonym: Maybe<Scalars['String']>;
  synonymTsvector: Maybe<Scalars['String']>;
};

/**
 * A condition to be used against `ConceptSynonym` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type ConceptSynonymCondition = {
  /** Checks for equality with the object’s `cid` field. */
  cid: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `synonym` field. */
  synonym: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `synonymTsvector` field. */
  synonymTsvector: InputMaybe<Scalars['String']>;
};

/** A filter to be used against `ConceptSynonym` object types. All fields are combined with a logical ‘and.’ */
export type ConceptSynonymFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ConceptSynonymFilter>>;
  /** Filter by the object’s `cid` field. */
  cid: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not: InputMaybe<ConceptSynonymFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ConceptSynonymFilter>>;
  /** Filter by the object’s `synonym` field. */
  synonym: InputMaybe<StringFilter>;
};

/** An input for mutations affecting `ConceptSynonym` */
export type ConceptSynonymInput = {
  cid: Scalars['Int'];
  synonym: InputMaybe<Scalars['String']>;
  synonymTsvector: InputMaybe<Scalars['String']>;
};

/** A connection to a list of `ConceptSynonym` values. */
export type ConceptSynonymsConnection = {
  __typename?: 'ConceptSynonymsConnection';
  /** A list of edges which contains the `ConceptSynonym` and cursor to aid in pagination. */
  edges: Array<ConceptSynonymsEdge>;
  /** A list of `ConceptSynonym` objects. */
  nodes: Array<Maybe<ConceptSynonym>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ConceptSynonym` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ConceptSynonym` edge in the connection. */
export type ConceptSynonymsEdge = {
  __typename?: 'ConceptSynonymsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ConceptSynonym` at the end of the edge. */
  node: Maybe<ConceptSynonym>;
};

/** Methods to use when ordering `ConceptSynonym`. */
export enum ConceptSynonymsOrderBy {
  CidAsc = 'CID_ASC',
  CidDesc = 'CID_DESC',
  Natural = 'NATURAL',
  SynonymAsc = 'SYNONYM_ASC',
  SynonymDesc = 'SYNONYM_DESC',
  SynonymTsvectorAsc = 'SYNONYM_TSVECTOR_ASC',
  SynonymTsvectorDesc = 'SYNONYM_TSVECTOR_DESC'
}

export type ConceptTreeElement = {
  __typename?: 'ConceptTreeElement';
  cid: Maybe<Scalars['Int']>;
  level: Maybe<Scalars['Int']>;
};

/** An input for mutations affecting `ConceptWeightedParent` */
export type ConceptWeightedParentInput = {
  cid: InputMaybe<Scalars['Int']>;
  parentCid: InputMaybe<Scalars['Int']>;
  semanticRelationshipWeight: InputMaybe<Scalars['Float']>;
};

/** A connection to a list of `Concept` values. */
export type ConceptsConnection = {
  __typename?: 'ConceptsConnection';
  /** A list of edges which contains the `Concept` and cursor to aid in pagination. */
  edges: Array<ConceptsEdge>;
  /** A list of `Concept` objects. */
  nodes: Array<Maybe<Concept>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `Concept` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `Concept` edge in the connection. */
export type ConceptsEdge = {
  __typename?: 'ConceptsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `Concept` at the end of the edge. */
  node: Maybe<Concept>;
};

/** Methods to use when ordering `Concept`. */
export enum ConceptsOrderBy {
  CidAsc = 'CID_ASC',
  CidDesc = 'CID_DESC',
  LabelAsc = 'LABEL_ASC',
  LabelDesc = 'LABEL_DESC',
  LabelTsvectorAsc = 'LABEL_TSVECTOR_ASC',
  LabelTsvectorDesc = 'LABEL_TSVECTOR_DESC',
  Natural = 'NATURAL',
  OntidAsc = 'ONTID_ASC',
  OntidDesc = 'ONTID_DESC',
  OntCodeAsc = 'ONT_CODE_ASC',
  OntCodeDesc = 'ONT_CODE_DESC',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC'
}

/** All input for the create `AnnotationGroup` mutation. */
export type CreateAnnotationGroupInput = {
  /** The `AnnotationGroup` to be created by this mutation. */
  annotationGroup: AnnotationGroupInput;
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
};

/** The output of our create `AnnotationGroup` mutation. */
export type CreateAnnotationGroupPayload = {
  __typename?: 'CreateAnnotationGroupPayload';
  /** The `AnnotationGroup` that was created by this mutation. */
  annotationGroup: Maybe<AnnotationGroup>;
  /** An edge for our `AnnotationGroup`. May be used by Relay 1. */
  annotationGroupEdge: Maybe<AnnotationGroupsEdge>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `AnnotationGroup` mutation. */
export type CreateAnnotationGroupPayloadAnnotationGroupEdgeArgs = {
  orderBy?: InputMaybe<Array<AnnotationGroupsOrderBy>>;
};

/** All input for the create `AnnotationValue` mutation. */
export type CreateAnnotationValueInput = {
  /** The `AnnotationValue` to be created by this mutation. */
  annotationValue: AnnotationValueInput;
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
};

/** The output of our create `AnnotationValue` mutation. */
export type CreateAnnotationValuePayload = {
  __typename?: 'CreateAnnotationValuePayload';
  /** Reads a single `AnnotationGroup` that is related to this `AnnotationValue`. */
  annotationGroup: Maybe<AnnotationGroup>;
  /** The `AnnotationValue` that was created by this mutation. */
  annotationValue: Maybe<AnnotationValue>;
  /** An edge for our `AnnotationValue`. May be used by Relay 1. */
  annotationValueEdge: Maybe<AnnotationValuesEdge>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `AnnotationValue` mutation. */
export type CreateAnnotationValuePayloadAnnotationValueEdgeArgs = {
  orderBy?: InputMaybe<Array<AnnotationValuesOrderBy>>;
};

/** All input for the create `ApiStudiesBase` mutation. */
export type CreateApiStudiesBaseInput = {
  /** The `ApiStudiesBase` to be created by this mutation. */
  apiStudiesBase: ApiStudiesBaseInput;
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
};

/** The output of our create `ApiStudiesBase` mutation. */
export type CreateApiStudiesBasePayload = {
  __typename?: 'CreateApiStudiesBasePayload';
  /** The `ApiStudiesBase` that was created by this mutation. */
  apiStudiesBase: Maybe<ApiStudiesBase>;
  /** An edge for our `ApiStudiesBase`. May be used by Relay 1. */
  apiStudiesBaseEdge: Maybe<ApiStudiesBasesEdge>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `ApiStudiesBase` mutation. */
export type CreateApiStudiesBasePayloadApiStudiesBaseEdgeArgs = {
  orderBy?: InputMaybe<Array<ApiStudiesBasesOrderBy>>;
};

/** All input for the create `ApiStudiesBulkRna` mutation. */
export type CreateApiStudiesBulkRnaInput = {
  /** The `ApiStudiesBulkRna` to be created by this mutation. */
  apiStudiesBulkRna: ApiStudiesBulkRnaInput;
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
};

/** The output of our create `ApiStudiesBulkRna` mutation. */
export type CreateApiStudiesBulkRnaPayload = {
  __typename?: 'CreateApiStudiesBulkRnaPayload';
  /** The `ApiStudiesBulkRna` that was created by this mutation. */
  apiStudiesBulkRna: Maybe<ApiStudiesBulkRna>;
  /** An edge for our `ApiStudiesBulkRna`. May be used by Relay 1. */
  apiStudiesBulkRnaEdge: Maybe<ApiStudiesBulkRnasEdge>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `ApiStudiesBulkRna` mutation. */
export type CreateApiStudiesBulkRnaPayloadApiStudiesBulkRnaEdgeArgs = {
  orderBy?: InputMaybe<Array<ApiStudiesBulkRnasOrderBy>>;
};

/** All input for the create `ApiStudiesSingleCell` mutation. */
export type CreateApiStudiesSingleCellInput = {
  /** The `ApiStudiesSingleCell` to be created by this mutation. */
  apiStudiesSingleCell: ApiStudiesSingleCellInput;
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
};

/** The output of our create `ApiStudiesSingleCell` mutation. */
export type CreateApiStudiesSingleCellPayload = {
  __typename?: 'CreateApiStudiesSingleCellPayload';
  /** The `ApiStudiesSingleCell` that was created by this mutation. */
  apiStudiesSingleCell: Maybe<ApiStudiesSingleCell>;
  /** An edge for our `ApiStudiesSingleCell`. May be used by Relay 1. */
  apiStudiesSingleCellEdge: Maybe<ApiStudiesSingleCellsEdge>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `ApiStudiesSingleCell` mutation. */
export type CreateApiStudiesSingleCellPayloadApiStudiesSingleCellEdgeArgs = {
  orderBy?: InputMaybe<Array<ApiStudiesSingleCellsOrderBy>>;
};

/** All input for the create `ApiStudyH5Download` mutation. */
export type CreateApiStudyH5DownloadInput = {
  /** The `ApiStudyH5Download` to be created by this mutation. */
  apiStudyH5Download: ApiStudyH5DownloadInput;
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
};

/** The output of our create `ApiStudyH5Download` mutation. */
export type CreateApiStudyH5DownloadPayload = {
  __typename?: 'CreateApiStudyH5DownloadPayload';
  /** The `ApiStudyH5Download` that was created by this mutation. */
  apiStudyH5Download: Maybe<ApiStudyH5Download>;
  /** An edge for our `ApiStudyH5Download`. May be used by Relay 1. */
  apiStudyH5DownloadEdge: Maybe<ApiStudyH5DownloadsEdge>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `ApiStudyH5Download` mutation. */
export type CreateApiStudyH5DownloadPayloadApiStudyH5DownloadEdgeArgs = {
  orderBy?: InputMaybe<Array<ApiStudyH5DownloadsOrderBy>>;
};

/** All input for the create `ConceptHierarchy` mutation. */
export type CreateConceptHierarchyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `ConceptHierarchy` to be created by this mutation. */
  conceptHierarchy: ConceptHierarchyInput;
};

/** The output of our create `ConceptHierarchy` mutation. */
export type CreateConceptHierarchyPayload = {
  __typename?: 'CreateConceptHierarchyPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Reads a single `Concept` that is related to this `ConceptHierarchy`. */
  conceptByCid: Maybe<Concept>;
  /** Reads a single `Concept` that is related to this `ConceptHierarchy`. */
  conceptByParentCid: Maybe<Concept>;
  /** The `ConceptHierarchy` that was created by this mutation. */
  conceptHierarchy: Maybe<ConceptHierarchy>;
  /** An edge for our `ConceptHierarchy`. May be used by Relay 1. */
  conceptHierarchyEdge: Maybe<ConceptHierarchiesEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `ConceptHierarchy` mutation. */
export type CreateConceptHierarchyPayloadConceptHierarchyEdgeArgs = {
  orderBy?: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
};

/** All input for the create `Concept` mutation. */
export type CreateConceptInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `Concept` to be created by this mutation. */
  concept: ConceptInput;
};

/** The output of our create `Concept` mutation. */
export type CreateConceptPayload = {
  __typename?: 'CreateConceptPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** The `Concept` that was created by this mutation. */
  concept: Maybe<Concept>;
  /** An edge for our `Concept`. May be used by Relay 1. */
  conceptEdge: Maybe<ConceptsEdge>;
  /** Reads a single `Ontology` that is related to this `Concept`. */
  ontologyByOntid: Maybe<Ontology>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `Concept` mutation. */
export type CreateConceptPayloadConceptEdgeArgs = {
  orderBy?: InputMaybe<Array<ConceptsOrderBy>>;
};

/** All input for the create `ConceptSynonym` mutation. */
export type CreateConceptSynonymInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `ConceptSynonym` to be created by this mutation. */
  conceptSynonym: ConceptSynonymInput;
};

/** The output of our create `ConceptSynonym` mutation. */
export type CreateConceptSynonymPayload = {
  __typename?: 'CreateConceptSynonymPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Reads a single `Concept` that is related to this `ConceptSynonym`. */
  conceptByCid: Maybe<Concept>;
  /** The `ConceptSynonym` that was created by this mutation. */
  conceptSynonym: Maybe<ConceptSynonym>;
  /** An edge for our `ConceptSynonym`. May be used by Relay 1. */
  conceptSynonymEdge: Maybe<ConceptSynonymsEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `ConceptSynonym` mutation. */
export type CreateConceptSynonymPayloadConceptSynonymEdgeArgs = {
  orderBy?: InputMaybe<Array<ConceptSynonymsOrderBy>>;
};

/** All input for the create `DifferentialExpression` mutation. */
export type CreateDifferentialExpressionInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `DifferentialExpression` to be created by this mutation. */
  differentialExpression: DifferentialExpressionInput;
};

/** The output of our create `DifferentialExpression` mutation. */
export type CreateDifferentialExpressionPayload = {
  __typename?: 'CreateDifferentialExpressionPayload';
  /** Reads a single `AnnotationValue` that is related to this `DifferentialExpression`. */
  annotationValue: Maybe<AnnotationValue>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** The `DifferentialExpression` that was created by this mutation. */
  differentialExpression: Maybe<DifferentialExpression>;
  /** An edge for our `DifferentialExpression`. May be used by Relay 1. */
  differentialExpressionEdge: Maybe<DifferentialExpressionsEdge>;
  /** Reads a single `OmicsBase` that is related to this `DifferentialExpression`. */
  omics: Maybe<OmicsBase>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `Study` that is related to this `DifferentialExpression`. */
  study: Maybe<Study>;
};


/** The output of our create `DifferentialExpression` mutation. */
export type CreateDifferentialExpressionPayloadDifferentialExpressionEdgeArgs = {
  orderBy?: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};

/** All input for the create `OmicsBase` mutation. */
export type CreateOmicsBaseInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `OmicsBase` to be created by this mutation. */
  omicsBase: OmicsBaseInput;
};

/** The output of our create `OmicsBase` mutation. */
export type CreateOmicsBasePayload = {
  __typename?: 'CreateOmicsBasePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** The `OmicsBase` that was created by this mutation. */
  omicsBase: Maybe<OmicsBase>;
  /** An edge for our `OmicsBase`. May be used by Relay 1. */
  omicsBaseEdge: Maybe<OmicsBasesEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `OmicsBase` mutation. */
export type CreateOmicsBasePayloadOmicsBaseEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsBasesOrderBy>>;
};

/** All input for the create `OmicsGene` mutation. */
export type CreateOmicsGeneInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `OmicsGene` to be created by this mutation. */
  omicsGene: OmicsGeneInput;
};

/** The output of our create `OmicsGene` mutation. */
export type CreateOmicsGenePayload = {
  __typename?: 'CreateOmicsGenePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Reads a single `OmicsBase` that is related to this `OmicsGene`. */
  gene: Maybe<OmicsBase>;
  /** The `OmicsGene` that was created by this mutation. */
  omicsGene: Maybe<OmicsGene>;
  /** An edge for our `OmicsGene`. May be used by Relay 1. */
  omicsGeneEdge: Maybe<OmicsGenesEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `OmicsGene` mutation. */
export type CreateOmicsGenePayloadOmicsGeneEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsGenesOrderBy>>;
};

/** All input for the create `OmicsProteinAntibodyTagGene` mutation. */
export type CreateOmicsProteinAntibodyTagGeneInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `OmicsProteinAntibodyTagGene` to be created by this mutation. */
  omicsProteinAntibodyTagGene: OmicsProteinAntibodyTagGeneInput;
};

/** The output of our create `OmicsProteinAntibodyTagGene` mutation. */
export type CreateOmicsProteinAntibodyTagGenePayload = {
  __typename?: 'CreateOmicsProteinAntibodyTagGenePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Reads a single `OmicsGene` that is related to this `OmicsProteinAntibodyTagGene`. */
  gene: Maybe<OmicsGene>;
  /** The `OmicsProteinAntibodyTagGene` that was created by this mutation. */
  omicsProteinAntibodyTagGene: Maybe<OmicsProteinAntibodyTagGene>;
  /** An edge for our `OmicsProteinAntibodyTagGene`. May be used by Relay 1. */
  omicsProteinAntibodyTagGeneEdge: Maybe<OmicsProteinAntibodyTagGenesEdge>;
  /** Reads a single `OmicsProteinAntibodyTag` that is related to this `OmicsProteinAntibodyTagGene`. */
  proteinAntibodyTag: Maybe<OmicsProteinAntibodyTag>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `OmicsProteinAntibodyTagGene` mutation. */
export type CreateOmicsProteinAntibodyTagGenePayloadOmicsProteinAntibodyTagGeneEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsProteinAntibodyTagGenesOrderBy>>;
};

/** All input for the create `OmicsProteinAntibodyTag` mutation. */
export type CreateOmicsProteinAntibodyTagInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `OmicsProteinAntibodyTag` to be created by this mutation. */
  omicsProteinAntibodyTag: OmicsProteinAntibodyTagInput;
};

/** The output of our create `OmicsProteinAntibodyTag` mutation. */
export type CreateOmicsProteinAntibodyTagPayload = {
  __typename?: 'CreateOmicsProteinAntibodyTagPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** The `OmicsProteinAntibodyTag` that was created by this mutation. */
  omicsProteinAntibodyTag: Maybe<OmicsProteinAntibodyTag>;
  /** An edge for our `OmicsProteinAntibodyTag`. May be used by Relay 1. */
  omicsProteinAntibodyTagEdge: Maybe<OmicsProteinAntibodyTagsEdge>;
  /** Reads a single `OmicsBase` that is related to this `OmicsProteinAntibodyTag`. */
  proteinAntibodyTag: Maybe<OmicsBase>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `OmicsProteinAntibodyTag` mutation. */
export type CreateOmicsProteinAntibodyTagPayloadOmicsProteinAntibodyTagEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsProteinAntibodyTagsOrderBy>>;
};

/** All input for the create `OmicsRegionGene` mutation. */
export type CreateOmicsRegionGeneInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `OmicsRegionGene` to be created by this mutation. */
  omicsRegionGene: OmicsRegionGeneInput;
};

/** The output of our create `OmicsRegionGene` mutation. */
export type CreateOmicsRegionGenePayload = {
  __typename?: 'CreateOmicsRegionGenePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Reads a single `OmicsGene` that is related to this `OmicsRegionGene`. */
  gene: Maybe<OmicsGene>;
  /** The `OmicsRegionGene` that was created by this mutation. */
  omicsRegionGene: Maybe<OmicsRegionGene>;
  /** An edge for our `OmicsRegionGene`. May be used by Relay 1. */
  omicsRegionGeneEdge: Maybe<OmicsRegionGenesEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `OmicsRegion` that is related to this `OmicsRegionGene`. */
  region: Maybe<OmicsRegion>;
};


/** The output of our create `OmicsRegionGene` mutation. */
export type CreateOmicsRegionGenePayloadOmicsRegionGeneEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsRegionGenesOrderBy>>;
};

/** All input for the create `OmicsRegion` mutation. */
export type CreateOmicsRegionInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `OmicsRegion` to be created by this mutation. */
  omicsRegion: OmicsRegionInput;
};

/** The output of our create `OmicsRegion` mutation. */
export type CreateOmicsRegionPayload = {
  __typename?: 'CreateOmicsRegionPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** The `OmicsRegion` that was created by this mutation. */
  omicsRegion: Maybe<OmicsRegion>;
  /** An edge for our `OmicsRegion`. May be used by Relay 1. */
  omicsRegionEdge: Maybe<OmicsRegionsEdge>;
  /** Reads a single `OmicsBase` that is related to this `OmicsRegion`. */
  omics_region_newNameHere: Maybe<OmicsBase>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `OmicsRegion` mutation. */
export type CreateOmicsRegionPayloadOmicsRegionEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsRegionsOrderBy>>;
};

/** All input for the create `OmicsTranscriptionFactorGene` mutation. */
export type CreateOmicsTranscriptionFactorGeneInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `OmicsTranscriptionFactorGene` to be created by this mutation. */
  omicsTranscriptionFactorGene: OmicsTranscriptionFactorGeneInput;
};

/** The output of our create `OmicsTranscriptionFactorGene` mutation. */
export type CreateOmicsTranscriptionFactorGenePayload = {
  __typename?: 'CreateOmicsTranscriptionFactorGenePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Reads a single `OmicsGene` that is related to this `OmicsTranscriptionFactorGene`. */
  gene: Maybe<OmicsGene>;
  /** The `OmicsTranscriptionFactorGene` that was created by this mutation. */
  omicsTranscriptionFactorGene: Maybe<OmicsTranscriptionFactorGene>;
  /** An edge for our `OmicsTranscriptionFactorGene`. May be used by Relay 1. */
  omicsTranscriptionFactorGeneEdge: Maybe<OmicsTranscriptionFactorGenesEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `OmicsTranscriptionFactor` that is related to this `OmicsTranscriptionFactorGene`. */
  transcriptionFactor: Maybe<OmicsTranscriptionFactor>;
};


/** The output of our create `OmicsTranscriptionFactorGene` mutation. */
export type CreateOmicsTranscriptionFactorGenePayloadOmicsTranscriptionFactorGeneEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsTranscriptionFactorGenesOrderBy>>;
};

/** All input for the create `OmicsTranscriptionFactor` mutation. */
export type CreateOmicsTranscriptionFactorInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `OmicsTranscriptionFactor` to be created by this mutation. */
  omicsTranscriptionFactor: OmicsTranscriptionFactorInput;
};

/** The output of our create `OmicsTranscriptionFactor` mutation. */
export type CreateOmicsTranscriptionFactorPayload = {
  __typename?: 'CreateOmicsTranscriptionFactorPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Reads a single `OmicsBase` that is related to this `OmicsTranscriptionFactor`. */
  omics: Maybe<OmicsBase>;
  /** The `OmicsTranscriptionFactor` that was created by this mutation. */
  omicsTranscriptionFactor: Maybe<OmicsTranscriptionFactor>;
  /** An edge for our `OmicsTranscriptionFactor`. May be used by Relay 1. */
  omicsTranscriptionFactorEdge: Maybe<OmicsTranscriptionFactorsEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `OmicsTranscriptionFactor` mutation. */
export type CreateOmicsTranscriptionFactorPayloadOmicsTranscriptionFactorEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsTranscriptionFactorsOrderBy>>;
};

/** All input for the create `Ontology` mutation. */
export type CreateOntologyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `Ontology` to be created by this mutation. */
  ontology: OntologyInput;
};

/** The output of our create `Ontology` mutation. */
export type CreateOntologyPayload = {
  __typename?: 'CreateOntologyPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** The `Ontology` that was created by this mutation. */
  ontology: Maybe<Ontology>;
  /** An edge for our `Ontology`. May be used by Relay 1. */
  ontologyEdge: Maybe<OntologiesEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our create `Ontology` mutation. */
export type CreateOntologyPayloadOntologyEdgeArgs = {
  orderBy?: InputMaybe<Array<OntologiesOrderBy>>;
};

/** All input for the create `ReferenceStudy` mutation. */
export type CreateReferenceStudyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `ReferenceStudy` to be created by this mutation. */
  referenceStudy: ReferenceStudyInput;
};

/** The output of our create `ReferenceStudy` mutation. */
export type CreateReferenceStudyPayload = {
  __typename?: 'CreateReferenceStudyPayload';
  /** Reads a single `AnnotationGroup` that is related to this `ReferenceStudy`. */
  celltypeAnnotationGroup: Maybe<AnnotationGroup>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `ReferenceStudyOverview` that is related to this `ReferenceStudy`. */
  refStudy: Maybe<ReferenceStudyOverview>;
  /** The `ReferenceStudy` that was created by this mutation. */
  referenceStudy: Maybe<ReferenceStudy>;
  /** An edge for our `ReferenceStudy`. May be used by Relay 1. */
  referenceStudyEdge: Maybe<ReferenceStudiesEdge>;
  /** Reads a single `Study` that is related to this `ReferenceStudy`. */
  study: Maybe<Study>;
  /** Reads a single `AnnotationGroup` that is related to this `ReferenceStudy`. */
  tissueAnnotationGroup: Maybe<AnnotationGroup>;
};


/** The output of our create `ReferenceStudy` mutation. */
export type CreateReferenceStudyPayloadReferenceStudyEdgeArgs = {
  orderBy?: InputMaybe<Array<ReferenceStudiesOrderBy>>;
};

/** All input for the `createS3TempCredentials` mutation. */
export type CreateS3TempCredentialsInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
};

/** The output of our `createS3TempCredentials` mutation. */
export type CreateS3TempCredentialsPayload = {
  __typename?: 'CreateS3TempCredentialsPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  strings: Maybe<Array<Maybe<Scalars['String']>>>;
};

/** All input for the create `StudyAdministrableCurrentuser` mutation. */
export type CreateStudyAdministrableCurrentuserInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `StudyAdministrableCurrentuser` to be created by this mutation. */
  studyAdministrableCurrentuser: StudyAdministrableCurrentuserInput;
};

/** The output of our create `StudyAdministrableCurrentuser` mutation. */
export type CreateStudyAdministrableCurrentuserPayload = {
  __typename?: 'CreateStudyAdministrableCurrentuserPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** The `StudyAdministrableCurrentuser` that was created by this mutation. */
  studyAdministrableCurrentuser: Maybe<StudyAdministrableCurrentuser>;
  /** An edge for our `StudyAdministrableCurrentuser`. May be used by Relay 1. */
  studyAdministrableCurrentuserEdge: Maybe<StudyAdministrableCurrentusersEdge>;
};


/** The output of our create `StudyAdministrableCurrentuser` mutation. */
export type CreateStudyAdministrableCurrentuserPayloadStudyAdministrableCurrentuserEdgeArgs = {
  orderBy?: InputMaybe<Array<StudyAdministrableCurrentusersOrderBy>>;
};

/** All input for the create `StudyAnnotationGroupUi` mutation. */
export type CreateStudyAnnotationGroupUiInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `StudyAnnotationGroupUi` to be created by this mutation. */
  studyAnnotationGroupUi: StudyAnnotationGroupUiInput;
};

/** The output of our create `StudyAnnotationGroupUi` mutation. */
export type CreateStudyAnnotationGroupUiPayload = {
  __typename?: 'CreateStudyAnnotationGroupUiPayload';
  /** Reads a single `AnnotationGroup` that is related to this `StudyAnnotationGroupUi`. */
  annotationGroup: Maybe<AnnotationGroup>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudyAnnotationGroupUi`. */
  study: Maybe<Study>;
  /** The `StudyAnnotationGroupUi` that was created by this mutation. */
  studyAnnotationGroupUi: Maybe<StudyAnnotationGroupUi>;
  /** An edge for our `StudyAnnotationGroupUi`. May be used by Relay 1. */
  studyAnnotationGroupUiEdge: Maybe<StudyAnnotationGroupUisEdge>;
};


/** The output of our create `StudyAnnotationGroupUi` mutation. */
export type CreateStudyAnnotationGroupUiPayloadStudyAnnotationGroupUiEdgeArgs = {
  orderBy?: InputMaybe<Array<StudyAnnotationGroupUisOrderBy>>;
};

/** All input for the `createStudyForCurrentUser` mutation. */
export type CreateStudyForCurrentUserInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  studyName: InputMaybe<Scalars['String']>;
};

/** The output of our `createStudyForCurrentUser` mutation. */
export type CreateStudyForCurrentUserPayload = {
  __typename?: 'CreateStudyForCurrentUserPayload';
  boolean: Maybe<Scalars['Boolean']>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};

/** All input for the `createStudyH5AdPresignedUrl` mutation. */
export type CreateStudyH5AdPresignedUrlInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  fileUrl: InputMaybe<Scalars['String']>;
};

/** The output of our `createStudyH5AdPresignedUrl` mutation. */
export type CreateStudyH5AdPresignedUrlPayload = {
  __typename?: 'CreateStudyH5AdPresignedUrlPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  string: Maybe<Scalars['String']>;
};

/** All input for the create `Study` mutation. */
export type CreateStudyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `Study` to be created by this mutation. */
  study: StudyInput;
};

/** All input for the create `StudyLayer` mutation. */
export type CreateStudyLayerInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `StudyLayer` to be created by this mutation. */
  studyLayer: StudyLayerInput;
};

/** The output of our create `StudyLayer` mutation. */
export type CreateStudyLayerPayload = {
  __typename?: 'CreateStudyLayerPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudyLayer`. */
  study: Maybe<Study>;
  /** The `StudyLayer` that was created by this mutation. */
  studyLayer: Maybe<StudyLayer>;
  /** An edge for our `StudyLayer`. May be used by Relay 1. */
  studyLayerEdge: Maybe<StudyLayersEdge>;
};


/** The output of our create `StudyLayer` mutation. */
export type CreateStudyLayerPayloadStudyLayerEdgeArgs = {
  orderBy?: InputMaybe<Array<StudyLayersOrderBy>>;
};

/** All input for the create `StudyOmic` mutation. */
export type CreateStudyOmicInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `StudyOmic` to be created by this mutation. */
  studyOmic: StudyOmicInput;
};

/** The output of our create `StudyOmic` mutation. */
export type CreateStudyOmicPayload = {
  __typename?: 'CreateStudyOmicPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Reads a single `OmicsBase` that is related to this `StudyOmic`. */
  omics: Maybe<OmicsBase>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudyOmic`. */
  study: Maybe<Study>;
  /** The `StudyOmic` that was created by this mutation. */
  studyOmic: Maybe<StudyOmic>;
  /** An edge for our `StudyOmic`. May be used by Relay 1. */
  studyOmicEdge: Maybe<StudyOmicsEdge>;
};


/** The output of our create `StudyOmic` mutation. */
export type CreateStudyOmicPayloadStudyOmicEdgeArgs = {
  orderBy?: InputMaybe<Array<StudyOmicsOrderBy>>;
};

/** All input for the create `StudyOverview` mutation. */
export type CreateStudyOverviewInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `StudyOverview` to be created by this mutation. */
  studyOverview: StudyOverviewInput;
};

/** The output of our create `StudyOverview` mutation. */
export type CreateStudyOverviewPayload = {
  __typename?: 'CreateStudyOverviewPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** The `StudyOverview` that was created by this mutation. */
  studyOverview: Maybe<StudyOverview>;
  /** An edge for our `StudyOverview`. May be used by Relay 1. */
  studyOverviewEdge: Maybe<StudyOverviewsEdge>;
};


/** The output of our create `StudyOverview` mutation. */
export type CreateStudyOverviewPayloadStudyOverviewEdgeArgs = {
  orderBy?: InputMaybe<Array<StudyOverviewsOrderBy>>;
};

/** The output of our create `Study` mutation. */
export type CreateStudyPayload = {
  __typename?: 'CreateStudyPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** The `Study` that was created by this mutation. */
  study: Maybe<Study>;
  /** An edge for our `Study`. May be used by Relay 1. */
  studyEdge: Maybe<StudiesEdge>;
};


/** The output of our create `Study` mutation. */
export type CreateStudyPayloadStudyEdgeArgs = {
  orderBy?: InputMaybe<Array<StudiesOrderBy>>;
};

/** All input for the create `StudySampleAnnotation` mutation. */
export type CreateStudySampleAnnotationInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `StudySampleAnnotation` to be created by this mutation. */
  studySampleAnnotation: StudySampleAnnotationInput;
};

/** The output of our create `StudySampleAnnotation` mutation. */
export type CreateStudySampleAnnotationPayload = {
  __typename?: 'CreateStudySampleAnnotationPayload';
  /** Reads a single `AnnotationValue` that is related to this `StudySampleAnnotation`. */
  annotationValue: Maybe<AnnotationValue>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudySampleAnnotation`. */
  study: Maybe<Study>;
  /** The `StudySampleAnnotation` that was created by this mutation. */
  studySampleAnnotation: Maybe<StudySampleAnnotation>;
  /** An edge for our `StudySampleAnnotation`. May be used by Relay 1. */
  studySampleAnnotationEdge: Maybe<StudySampleAnnotationsEdge>;
};


/** The output of our create `StudySampleAnnotation` mutation. */
export type CreateStudySampleAnnotationPayloadStudySampleAnnotationEdgeArgs = {
  orderBy?: InputMaybe<Array<StudySampleAnnotationsOrderBy>>;
};

/** All input for the create `StudySample` mutation. */
export type CreateStudySampleInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `StudySample` to be created by this mutation. */
  studySample: StudySampleInput;
};

/** The output of our create `StudySample` mutation. */
export type CreateStudySamplePayload = {
  __typename?: 'CreateStudySamplePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudySample`. */
  study: Maybe<Study>;
  /** The `StudySample` that was created by this mutation. */
  studySample: Maybe<StudySample>;
  /** An edge for our `StudySample`. May be used by Relay 1. */
  studySampleEdge: Maybe<StudySamplesEdge>;
};


/** The output of our create `StudySample` mutation. */
export type CreateStudySamplePayloadStudySampleEdgeArgs = {
  orderBy?: InputMaybe<Array<StudySamplesOrderBy>>;
};

/** All input for the create `StudySampleProjection` mutation. */
export type CreateStudySampleProjectionInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `StudySampleProjection` to be created by this mutation. */
  studySampleProjection: StudySampleProjectionInput;
};

/** The output of our create `StudySampleProjection` mutation. */
export type CreateStudySampleProjectionPayload = {
  __typename?: 'CreateStudySampleProjectionPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** The `StudySampleProjection` that was created by this mutation. */
  studySampleProjection: Maybe<StudySampleProjection>;
  /** An edge for our `StudySampleProjection`. May be used by Relay 1. */
  studySampleProjectionEdge: Maybe<StudySampleProjectionsEdge>;
  /** Reads a single `StudySample` that is related to this `StudySampleProjection`. */
  studyStudySample: Maybe<StudySample>;
};


/** The output of our create `StudySampleProjection` mutation. */
export type CreateStudySampleProjectionPayloadStudySampleProjectionEdgeArgs = {
  orderBy?: InputMaybe<Array<StudySampleProjectionsOrderBy>>;
};

/** All input for the `createStudyUpload` mutation. */
export type CreateStudyUploadInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  filename: InputMaybe<Scalars['String']>;
};

/** The output of our `createStudyUpload` mutation. */
export type CreateStudyUploadPayload = {
  __typename?: 'CreateStudyUploadPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  json: Maybe<Scalars['JSON']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};

/** All input for the create `StudyVisibleCurrentuser` mutation. */
export type CreateStudyVisibleCurrentuserInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `StudyVisibleCurrentuser` to be created by this mutation. */
  studyVisibleCurrentuser: StudyVisibleCurrentuserInput;
};

/** The output of our create `StudyVisibleCurrentuser` mutation. */
export type CreateStudyVisibleCurrentuserPayload = {
  __typename?: 'CreateStudyVisibleCurrentuserPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** The `StudyVisibleCurrentuser` that was created by this mutation. */
  studyVisibleCurrentuser: Maybe<StudyVisibleCurrentuser>;
  /** An edge for our `StudyVisibleCurrentuser`. May be used by Relay 1. */
  studyVisibleCurrentuserEdge: Maybe<StudyVisibleCurrentusersEdge>;
};


/** The output of our create `StudyVisibleCurrentuser` mutation. */
export type CreateStudyVisibleCurrentuserPayloadStudyVisibleCurrentuserEdgeArgs = {
  orderBy?: InputMaybe<Array<StudyVisibleCurrentusersOrderBy>>;
};

/** All input for the create `UserAnnotationGroup` mutation. */
export type CreateUserAnnotationGroupInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The `UserAnnotationGroup` to be created by this mutation. */
  userAnnotationGroup: UserAnnotationGroupInput;
};

/** The output of our create `UserAnnotationGroup` mutation. */
export type CreateUserAnnotationGroupPayload = {
  __typename?: 'CreateUserAnnotationGroupPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `AnnotationGroup` that is related to this `UserAnnotationGroup`. */
  savedAsAnnotationGroup: Maybe<AnnotationGroup>;
  /** Reads a single `Study` that is related to this `UserAnnotationGroup`. */
  study: Maybe<Study>;
  /** The `UserAnnotationGroup` that was created by this mutation. */
  userAnnotationGroup: Maybe<UserAnnotationGroup>;
  /** An edge for our `UserAnnotationGroup`. May be used by Relay 1. */
  userAnnotationGroupEdge: Maybe<UserAnnotationGroupsEdge>;
};


/** The output of our create `UserAnnotationGroup` mutation. */
export type CreateUserAnnotationGroupPayloadUserAnnotationGroupEdgeArgs = {
  orderBy?: InputMaybe<Array<UserAnnotationGroupsOrderBy>>;
};

/** A filter to be used against Datetime fields. All fields are combined with a logical ‘and.’ */
export type DatetimeFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<Scalars['Datetime']>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<Scalars['Datetime']>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<Scalars['Datetime']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<Scalars['Datetime']>;
  /** Included in the specified list. */
  in: InputMaybe<Array<Scalars['Datetime']>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<Scalars['Datetime']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<Scalars['Datetime']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<Scalars['Datetime']>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<Scalars['Datetime']>;
  /** Not included in the specified list. */
  notIn: InputMaybe<Array<Scalars['Datetime']>>;
};

/** All input for the `deleteAllStudyData` mutation. */
export type DeleteAllStudyDataInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  pStudyId: InputMaybe<Scalars['Int']>;
};

/** The output of our `deleteAllStudyData` mutation. */
export type DeleteAllStudyDataPayload = {
  __typename?: 'DeleteAllStudyDataPayload';
  boolean: Maybe<Scalars['Boolean']>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};

/** All input for the `deleteAnnotationGroupByNodeId` mutation. */
export type DeleteAnnotationGroupByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `AnnotationGroup` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteAnnotationGroup` mutation. */
export type DeleteAnnotationGroupInput = {
  annotationGroupId: Scalars['Int'];
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
};

/** The output of our delete `AnnotationGroup` mutation. */
export type DeleteAnnotationGroupPayload = {
  __typename?: 'DeleteAnnotationGroupPayload';
  /** The `AnnotationGroup` that was deleted by this mutation. */
  annotationGroup: Maybe<AnnotationGroup>;
  /** An edge for our `AnnotationGroup`. May be used by Relay 1. */
  annotationGroupEdge: Maybe<AnnotationGroupsEdge>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedAnnotationGroupNodeId: Maybe<Scalars['ID']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our delete `AnnotationGroup` mutation. */
export type DeleteAnnotationGroupPayloadAnnotationGroupEdgeArgs = {
  orderBy?: InputMaybe<Array<AnnotationGroupsOrderBy>>;
};

/** All input for the `deleteAnnotationValueByNodeId` mutation. */
export type DeleteAnnotationValueByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `AnnotationValue` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteAnnotationValue` mutation. */
export type DeleteAnnotationValueInput = {
  annotationValueId: Scalars['Int'];
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
};

/** The output of our delete `AnnotationValue` mutation. */
export type DeleteAnnotationValuePayload = {
  __typename?: 'DeleteAnnotationValuePayload';
  /** Reads a single `AnnotationGroup` that is related to this `AnnotationValue`. */
  annotationGroup: Maybe<AnnotationGroup>;
  /** The `AnnotationValue` that was deleted by this mutation. */
  annotationValue: Maybe<AnnotationValue>;
  /** An edge for our `AnnotationValue`. May be used by Relay 1. */
  annotationValueEdge: Maybe<AnnotationValuesEdge>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedAnnotationValueNodeId: Maybe<Scalars['ID']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our delete `AnnotationValue` mutation. */
export type DeleteAnnotationValuePayloadAnnotationValueEdgeArgs = {
  orderBy?: InputMaybe<Array<AnnotationValuesOrderBy>>;
};

/** All input for the `deleteConceptByNodeId` mutation. */
export type DeleteConceptByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `Concept` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteConcept` mutation. */
export type DeleteConceptInput = {
  cid: Scalars['Int'];
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
};

/** The output of our delete `Concept` mutation. */
export type DeleteConceptPayload = {
  __typename?: 'DeleteConceptPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** The `Concept` that was deleted by this mutation. */
  concept: Maybe<Concept>;
  /** An edge for our `Concept`. May be used by Relay 1. */
  conceptEdge: Maybe<ConceptsEdge>;
  deletedConceptNodeId: Maybe<Scalars['ID']>;
  /** Reads a single `Ontology` that is related to this `Concept`. */
  ontologyByOntid: Maybe<Ontology>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our delete `Concept` mutation. */
export type DeleteConceptPayloadConceptEdgeArgs = {
  orderBy?: InputMaybe<Array<ConceptsOrderBy>>;
};

/** All input for the `deleteOmicsBaseByNodeId` mutation. */
export type DeleteOmicsBaseByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `OmicsBase` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteOmicsBase` mutation. */
export type DeleteOmicsBaseInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  omicsId: Scalars['Int'];
};

/** The output of our delete `OmicsBase` mutation. */
export type DeleteOmicsBasePayload = {
  __typename?: 'DeleteOmicsBasePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedOmicsBaseNodeId: Maybe<Scalars['ID']>;
  /** The `OmicsBase` that was deleted by this mutation. */
  omicsBase: Maybe<OmicsBase>;
  /** An edge for our `OmicsBase`. May be used by Relay 1. */
  omicsBaseEdge: Maybe<OmicsBasesEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our delete `OmicsBase` mutation. */
export type DeleteOmicsBasePayloadOmicsBaseEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsBasesOrderBy>>;
};

/** All input for the `deleteOmicsGeneByNodeId` mutation. */
export type DeleteOmicsGeneByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `OmicsGene` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteOmicsGene` mutation. */
export type DeleteOmicsGeneInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  geneId: Scalars['Int'];
};

/** The output of our delete `OmicsGene` mutation. */
export type DeleteOmicsGenePayload = {
  __typename?: 'DeleteOmicsGenePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedOmicsGeneNodeId: Maybe<Scalars['ID']>;
  /** Reads a single `OmicsBase` that is related to this `OmicsGene`. */
  gene: Maybe<OmicsBase>;
  /** The `OmicsGene` that was deleted by this mutation. */
  omicsGene: Maybe<OmicsGene>;
  /** An edge for our `OmicsGene`. May be used by Relay 1. */
  omicsGeneEdge: Maybe<OmicsGenesEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our delete `OmicsGene` mutation. */
export type DeleteOmicsGenePayloadOmicsGeneEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsGenesOrderBy>>;
};

/** All input for the `deleteOmicsProteinAntibodyTagByNodeId` mutation. */
export type DeleteOmicsProteinAntibodyTagByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `OmicsProteinAntibodyTag` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteOmicsProteinAntibodyTag` mutation. */
export type DeleteOmicsProteinAntibodyTagInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  proteinAntibodyTagId: Scalars['Int'];
};

/** The output of our delete `OmicsProteinAntibodyTag` mutation. */
export type DeleteOmicsProteinAntibodyTagPayload = {
  __typename?: 'DeleteOmicsProteinAntibodyTagPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedOmicsProteinAntibodyTagNodeId: Maybe<Scalars['ID']>;
  /** The `OmicsProteinAntibodyTag` that was deleted by this mutation. */
  omicsProteinAntibodyTag: Maybe<OmicsProteinAntibodyTag>;
  /** An edge for our `OmicsProteinAntibodyTag`. May be used by Relay 1. */
  omicsProteinAntibodyTagEdge: Maybe<OmicsProteinAntibodyTagsEdge>;
  /** Reads a single `OmicsBase` that is related to this `OmicsProteinAntibodyTag`. */
  proteinAntibodyTag: Maybe<OmicsBase>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our delete `OmicsProteinAntibodyTag` mutation. */
export type DeleteOmicsProteinAntibodyTagPayloadOmicsProteinAntibodyTagEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsProteinAntibodyTagsOrderBy>>;
};

/** All input for the `deleteOmicsRegionByNodeId` mutation. */
export type DeleteOmicsRegionByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `OmicsRegion` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteOmicsRegion` mutation. */
export type DeleteOmicsRegionInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  regionId: Scalars['Int'];
};

/** The output of our delete `OmicsRegion` mutation. */
export type DeleteOmicsRegionPayload = {
  __typename?: 'DeleteOmicsRegionPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedOmicsRegionNodeId: Maybe<Scalars['ID']>;
  /** The `OmicsRegion` that was deleted by this mutation. */
  omicsRegion: Maybe<OmicsRegion>;
  /** An edge for our `OmicsRegion`. May be used by Relay 1. */
  omicsRegionEdge: Maybe<OmicsRegionsEdge>;
  /** Reads a single `OmicsBase` that is related to this `OmicsRegion`. */
  omics_region_newNameHere: Maybe<OmicsBase>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our delete `OmicsRegion` mutation. */
export type DeleteOmicsRegionPayloadOmicsRegionEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsRegionsOrderBy>>;
};

/** All input for the `deleteOmicsTranscriptionFactorByNodeId` mutation. */
export type DeleteOmicsTranscriptionFactorByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `OmicsTranscriptionFactor` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteOmicsTranscriptionFactor` mutation. */
export type DeleteOmicsTranscriptionFactorInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  omicsId: Scalars['Int'];
};

/** The output of our delete `OmicsTranscriptionFactor` mutation. */
export type DeleteOmicsTranscriptionFactorPayload = {
  __typename?: 'DeleteOmicsTranscriptionFactorPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedOmicsTranscriptionFactorNodeId: Maybe<Scalars['ID']>;
  /** Reads a single `OmicsBase` that is related to this `OmicsTranscriptionFactor`. */
  omics: Maybe<OmicsBase>;
  /** The `OmicsTranscriptionFactor` that was deleted by this mutation. */
  omicsTranscriptionFactor: Maybe<OmicsTranscriptionFactor>;
  /** An edge for our `OmicsTranscriptionFactor`. May be used by Relay 1. */
  omicsTranscriptionFactorEdge: Maybe<OmicsTranscriptionFactorsEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our delete `OmicsTranscriptionFactor` mutation. */
export type DeleteOmicsTranscriptionFactorPayloadOmicsTranscriptionFactorEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsTranscriptionFactorsOrderBy>>;
};

/** All input for the `deleteOntologyByNodeId` mutation. */
export type DeleteOntologyByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `Ontology` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteOntology` mutation. */
export type DeleteOntologyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  ontid: Scalars['Int'];
};

/** The output of our delete `Ontology` mutation. */
export type DeleteOntologyPayload = {
  __typename?: 'DeleteOntologyPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedOntologyNodeId: Maybe<Scalars['ID']>;
  /** The `Ontology` that was deleted by this mutation. */
  ontology: Maybe<Ontology>;
  /** An edge for our `Ontology`. May be used by Relay 1. */
  ontologyEdge: Maybe<OntologiesEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our delete `Ontology` mutation. */
export type DeleteOntologyPayloadOntologyEdgeArgs = {
  orderBy?: InputMaybe<Array<OntologiesOrderBy>>;
};

/** All input for the `deleteStudyByNodeId` mutation. */
export type DeleteStudyByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `Study` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteStudy` mutation. */
export type DeleteStudyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  studyId: Scalars['Int'];
};

/** All input for the `deleteStudyLayerByNodeId` mutation. */
export type DeleteStudyLayerByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `StudyLayer` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteStudyLayer` mutation. */
export type DeleteStudyLayerInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  studyLayerId: Scalars['Int'];
};

/** The output of our delete `StudyLayer` mutation. */
export type DeleteStudyLayerPayload = {
  __typename?: 'DeleteStudyLayerPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedStudyLayerNodeId: Maybe<Scalars['ID']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudyLayer`. */
  study: Maybe<Study>;
  /** The `StudyLayer` that was deleted by this mutation. */
  studyLayer: Maybe<StudyLayer>;
  /** An edge for our `StudyLayer`. May be used by Relay 1. */
  studyLayerEdge: Maybe<StudyLayersEdge>;
};


/** The output of our delete `StudyLayer` mutation. */
export type DeleteStudyLayerPayloadStudyLayerEdgeArgs = {
  orderBy?: InputMaybe<Array<StudyLayersOrderBy>>;
};

/** The output of our delete `Study` mutation. */
export type DeleteStudyPayload = {
  __typename?: 'DeleteStudyPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedStudyNodeId: Maybe<Scalars['ID']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** The `Study` that was deleted by this mutation. */
  study: Maybe<Study>;
  /** An edge for our `Study`. May be used by Relay 1. */
  studyEdge: Maybe<StudiesEdge>;
};


/** The output of our delete `Study` mutation. */
export type DeleteStudyPayloadStudyEdgeArgs = {
  orderBy?: InputMaybe<Array<StudiesOrderBy>>;
};

/** All input for the `deleteStudySampleByNodeId` mutation. */
export type DeleteStudySampleByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `StudySample` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteStudySample` mutation. */
export type DeleteStudySampleInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  studyId: Scalars['Int'];
  studySampleId: Scalars['Int'];
};

/** The output of our delete `StudySample` mutation. */
export type DeleteStudySamplePayload = {
  __typename?: 'DeleteStudySamplePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedStudySampleNodeId: Maybe<Scalars['ID']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudySample`. */
  study: Maybe<Study>;
  /** The `StudySample` that was deleted by this mutation. */
  studySample: Maybe<StudySample>;
  /** An edge for our `StudySample`. May be used by Relay 1. */
  studySampleEdge: Maybe<StudySamplesEdge>;
};


/** The output of our delete `StudySample` mutation. */
export type DeleteStudySamplePayloadStudySampleEdgeArgs = {
  orderBy?: InputMaybe<Array<StudySamplesOrderBy>>;
};

export type DifferentialExpression = {
  __typename?: 'DifferentialExpression';
  /** Reads a single `AnnotationValue` that is related to this `DifferentialExpression`. */
  annotationValue: Maybe<AnnotationValue>;
  annotationValueId: Scalars['Int'];
  log2Foldchange: Maybe<Scalars['Float']>;
  /** Reads a single `OmicsBase` that is related to this `DifferentialExpression`. */
  omics: Maybe<OmicsBase>;
  omicsId: Scalars['Int'];
  pvalue: Maybe<Scalars['Float']>;
  pvalueAdj: Maybe<Scalars['Float']>;
  score: Maybe<Scalars['Float']>;
  /** Reads a single `Study` that is related to this `DifferentialExpression`. */
  study: Maybe<Study>;
  studyId: Scalars['Int'];
};

/**
 * A condition to be used against `DifferentialExpression` object types. All fields
 * are tested for equality and combined with a logical ‘and.’
 */
export type DifferentialExpressionCondition = {
  /** Checks for equality with the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `log2Foldchange` field. */
  log2Foldchange: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `pvalue` field. */
  pvalue: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `pvalueAdj` field. */
  pvalueAdj: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `score` field. */
  score: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `DifferentialExpression` object types. All fields are combined with a logical ‘and.’ */
export type DifferentialExpressionFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<DifferentialExpressionFilter>>;
  /** Filter by the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<IntFilter>;
  /** Filter by the object’s `log2Foldchange` field. */
  log2Foldchange: InputMaybe<FloatFilter>;
  /** Negates the expression. */
  not: InputMaybe<DifferentialExpressionFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<DifferentialExpressionFilter>>;
  /** Filter by the object’s `pvalue` field. */
  pvalue: InputMaybe<FloatFilter>;
  /** Filter by the object’s `pvalueAdj` field. */
  pvalueAdj: InputMaybe<FloatFilter>;
  /** Filter by the object’s `score` field. */
  score: InputMaybe<FloatFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `DifferentialExpression` */
export type DifferentialExpressionInput = {
  annotationValueId: Scalars['Int'];
  log2Foldchange: InputMaybe<Scalars['Float']>;
  omicsId: Scalars['Int'];
  pvalue: InputMaybe<Scalars['Float']>;
  pvalueAdj: InputMaybe<Scalars['Float']>;
  score: InputMaybe<Scalars['Float']>;
  studyId: Scalars['Int'];
};

export type DifferentialExpressionV = {
  __typename?: 'DifferentialExpressionV';
  annotationValueId: Maybe<Scalars['Int']>;
  displayName: Maybe<Scalars['String']>;
  displaySymbol: Maybe<Scalars['String']>;
  linkedGenes: Maybe<Array<Maybe<Scalars['Int']>>>;
  log2Foldchange: Maybe<Scalars['Float']>;
  omicsId: Maybe<Scalars['Int']>;
  omicsType: Maybe<OmicsType>;
  pvalue: Maybe<Scalars['Float']>;
  pvalueAdj: Maybe<Scalars['Float']>;
  score: Maybe<Scalars['Float']>;
  studyId: Maybe<Scalars['Int']>;
};

/**
 * A condition to be used against `DifferentialExpressionV` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type DifferentialExpressionVCondition = {
  /** Checks for equality with the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `displayName` field. */
  displayName: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `linkedGenes` field. */
  linkedGenes: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `log2Foldchange` field. */
  log2Foldchange: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `omicsType` field. */
  omicsType: InputMaybe<OmicsType>;
  /** Checks for equality with the object’s `pvalue` field. */
  pvalue: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `pvalueAdj` field. */
  pvalueAdj: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `score` field. */
  score: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `DifferentialExpressionV` object types. All fields are combined with a logical ‘and.’ */
export type DifferentialExpressionVFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<DifferentialExpressionVFilter>>;
  /** Filter by the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<IntFilter>;
  /** Filter by the object’s `displayName` field. */
  displayName: InputMaybe<StringFilter>;
  /** Filter by the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<StringFilter>;
  /** Filter by the object’s `linkedGenes` field. */
  linkedGenes: InputMaybe<IntListFilter>;
  /** Filter by the object’s `log2Foldchange` field. */
  log2Foldchange: InputMaybe<FloatFilter>;
  /** Negates the expression. */
  not: InputMaybe<DifferentialExpressionVFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Filter by the object’s `omicsType` field. */
  omicsType: InputMaybe<OmicsTypeFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<DifferentialExpressionVFilter>>;
  /** Filter by the object’s `pvalue` field. */
  pvalue: InputMaybe<FloatFilter>;
  /** Filter by the object’s `pvalueAdj` field. */
  pvalueAdj: InputMaybe<FloatFilter>;
  /** Filter by the object’s `score` field. */
  score: InputMaybe<FloatFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** A connection to a list of `DifferentialExpressionV` values. */
export type DifferentialExpressionVsConnection = {
  __typename?: 'DifferentialExpressionVsConnection';
  /** A list of edges which contains the `DifferentialExpressionV` and cursor to aid in pagination. */
  edges: Array<DifferentialExpressionVsEdge>;
  /** A list of `DifferentialExpressionV` objects. */
  nodes: Array<Maybe<DifferentialExpressionV>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `DifferentialExpressionV` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `DifferentialExpressionV` edge in the connection. */
export type DifferentialExpressionVsEdge = {
  __typename?: 'DifferentialExpressionVsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `DifferentialExpressionV` at the end of the edge. */
  node: Maybe<DifferentialExpressionV>;
};

/** Methods to use when ordering `DifferentialExpressionV`. */
export enum DifferentialExpressionVsOrderBy {
  AnnotationValueIdAsc = 'ANNOTATION_VALUE_ID_ASC',
  AnnotationValueIdDesc = 'ANNOTATION_VALUE_ID_DESC',
  DisplayNameAsc = 'DISPLAY_NAME_ASC',
  DisplayNameDesc = 'DISPLAY_NAME_DESC',
  DisplaySymbolAsc = 'DISPLAY_SYMBOL_ASC',
  DisplaySymbolDesc = 'DISPLAY_SYMBOL_DESC',
  LinkedGenesAsc = 'LINKED_GENES_ASC',
  LinkedGenesDesc = 'LINKED_GENES_DESC',
  Log2FoldchangeAsc = 'LOG2_FOLDCHANGE_ASC',
  Log2FoldchangeDesc = 'LOG2_FOLDCHANGE_DESC',
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  OmicsTypeAsc = 'OMICS_TYPE_ASC',
  OmicsTypeDesc = 'OMICS_TYPE_DESC',
  PvalueAdjAsc = 'PVALUE_ADJ_ASC',
  PvalueAdjDesc = 'PVALUE_ADJ_DESC',
  PvalueAsc = 'PVALUE_ASC',
  PvalueDesc = 'PVALUE_DESC',
  ScoreAsc = 'SCORE_ASC',
  ScoreDesc = 'SCORE_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

/** A connection to a list of `DifferentialExpression` values. */
export type DifferentialExpressionsConnection = {
  __typename?: 'DifferentialExpressionsConnection';
  /** A list of edges which contains the `DifferentialExpression` and cursor to aid in pagination. */
  edges: Array<DifferentialExpressionsEdge>;
  /** A list of `DifferentialExpression` objects. */
  nodes: Array<Maybe<DifferentialExpression>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `DifferentialExpression` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `DifferentialExpression` edge in the connection. */
export type DifferentialExpressionsEdge = {
  __typename?: 'DifferentialExpressionsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `DifferentialExpression` at the end of the edge. */
  node: Maybe<DifferentialExpression>;
};

/** Methods to use when ordering `DifferentialExpression`. */
export enum DifferentialExpressionsOrderBy {
  AnnotationValueIdAsc = 'ANNOTATION_VALUE_ID_ASC',
  AnnotationValueIdDesc = 'ANNOTATION_VALUE_ID_DESC',
  Log2FoldchangeAsc = 'LOG2_FOLDCHANGE_ASC',
  Log2FoldchangeDesc = 'LOG2_FOLDCHANGE_DESC',
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  PvalueAdjAsc = 'PVALUE_ADJ_ASC',
  PvalueAdjDesc = 'PVALUE_ADJ_DESC',
  PvalueAsc = 'PVALUE_ASC',
  PvalueDesc = 'PVALUE_DESC',
  ScoreAsc = 'SCORE_ASC',
  ScoreDesc = 'SCORE_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

export type ExpressionByAnnotation = {
  __typename?: 'ExpressionByAnnotation';
  annotationDisplayValue: Maybe<Scalars['String']>;
  annotationValueId: Maybe<Scalars['Int']>;
  boxplotParams: Maybe<BoxplotValue>;
  exprSamplesFraction: Maybe<Scalars['Float']>;
  mean: Maybe<Scalars['Float']>;
  median: Maybe<Scalars['Float']>;
  omicsId: Maybe<Scalars['Int']>;
  q3: Maybe<Scalars['Float']>;
  studyLayerId: Maybe<Scalars['Int']>;
  valueCount: Maybe<Scalars['Int']>;
  values: Maybe<Array<Maybe<Scalars['Float']>>>;
};

/** A filter to be used against `ExpressionByAnnotation` object types. All fields are combined with a logical ‘and.’ */
export type ExpressionByAnnotationFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ExpressionByAnnotationFilter>>;
  /** Filter by the object’s `annotationDisplayValue` field. */
  annotationDisplayValue: InputMaybe<StringFilter>;
  /** Filter by the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<IntFilter>;
  /** Filter by the object’s `boxplotParams` field. */
  boxplotParams: InputMaybe<BoxplotValueFilter>;
  /** Filter by the object’s `exprSamplesFraction` field. */
  exprSamplesFraction: InputMaybe<FloatFilter>;
  /** Filter by the object’s `mean` field. */
  mean: InputMaybe<FloatFilter>;
  /** Filter by the object’s `median` field. */
  median: InputMaybe<FloatFilter>;
  /** Negates the expression. */
  not: InputMaybe<ExpressionByAnnotationFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ExpressionByAnnotationFilter>>;
  /** Filter by the object’s `q3` field. */
  q3: InputMaybe<FloatFilter>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId: InputMaybe<IntFilter>;
  /** Filter by the object’s `valueCount` field. */
  valueCount: InputMaybe<IntFilter>;
  /** Filter by the object’s `values` field. */
  values: InputMaybe<FloatListFilter>;
};

/** A connection to a list of `ExpressionByAnnotation` values. */
export type ExpressionByAnnotationsConnection = {
  __typename?: 'ExpressionByAnnotationsConnection';
  /** A list of edges which contains the `ExpressionByAnnotation` and cursor to aid in pagination. */
  edges: Array<ExpressionByAnnotationsEdge>;
  /** A list of `ExpressionByAnnotation` objects. */
  nodes: Array<Maybe<ExpressionByAnnotation>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ExpressionByAnnotation` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ExpressionByAnnotation` edge in the connection. */
export type ExpressionByAnnotationsEdge = {
  __typename?: 'ExpressionByAnnotationsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ExpressionByAnnotation` at the end of the edge. */
  node: Maybe<ExpressionByAnnotation>;
};

export type ExpressionByOmic = {
  __typename?: 'ExpressionByOmic';
  omicsId: Maybe<Scalars['Int']>;
  studySampleIds: Maybe<Array<Maybe<Scalars['Int']>>>;
  values: Maybe<Array<Maybe<Scalars['Float']>>>;
};

/** A filter to be used against `ExpressionByOmic` object types. All fields are combined with a logical ‘and.’ */
export type ExpressionByOmicFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ExpressionByOmicFilter>>;
  /** Negates the expression. */
  not: InputMaybe<ExpressionByOmicFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ExpressionByOmicFilter>>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds: InputMaybe<IntListFilter>;
  /** Filter by the object’s `values` field. */
  values: InputMaybe<FloatListFilter>;
};

/** A connection to a list of `ExpressionByOmic` values. */
export type ExpressionByOmicsConnection = {
  __typename?: 'ExpressionByOmicsConnection';
  /** A list of edges which contains the `ExpressionByOmic` and cursor to aid in pagination. */
  edges: Array<ExpressionByOmicsEdge>;
  /** A list of `ExpressionByOmic` objects. */
  nodes: Array<Maybe<ExpressionByOmic>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ExpressionByOmic` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ExpressionByOmic` edge in the connection. */
export type ExpressionByOmicsEdge = {
  __typename?: 'ExpressionByOmicsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ExpressionByOmic` at the end of the edge. */
  node: Maybe<ExpressionByOmic>;
};

export type ExpressionByTwoAnnotation = {
  __typename?: 'ExpressionByTwoAnnotation';
  annotationDisplayValue: Maybe<Scalars['String']>;
  annotationValueId: Maybe<Scalars['Int']>;
  boxplotParams: Maybe<BoxplotValue>;
  color: Maybe<Scalars['String']>;
  exprSamplesFraction: Maybe<Scalars['Float']>;
  mean: Maybe<Scalars['Float']>;
  median: Maybe<Scalars['Float']>;
  nonZeroValueCount: Maybe<Scalars['Int']>;
  omicsId: Maybe<Scalars['Int']>;
  q3: Maybe<Scalars['Float']>;
  secondAnnotationDisplayValue: Maybe<Scalars['String']>;
  secondAnnotationValueId: Maybe<Scalars['Int']>;
  valueCount: Maybe<Scalars['Int']>;
  values: Maybe<Array<Maybe<Scalars['Float']>>>;
};

/** A filter to be used against `ExpressionByTwoAnnotation` object types. All fields are combined with a logical ‘and.’ */
export type ExpressionByTwoAnnotationFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ExpressionByTwoAnnotationFilter>>;
  /** Filter by the object’s `annotationDisplayValue` field. */
  annotationDisplayValue: InputMaybe<StringFilter>;
  /** Filter by the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<IntFilter>;
  /** Filter by the object’s `boxplotParams` field. */
  boxplotParams: InputMaybe<BoxplotValueFilter>;
  /** Filter by the object’s `color` field. */
  color: InputMaybe<StringFilter>;
  /** Filter by the object’s `exprSamplesFraction` field. */
  exprSamplesFraction: InputMaybe<FloatFilter>;
  /** Filter by the object’s `mean` field. */
  mean: InputMaybe<FloatFilter>;
  /** Filter by the object’s `median` field. */
  median: InputMaybe<FloatFilter>;
  /** Filter by the object’s `nonZeroValueCount` field. */
  nonZeroValueCount: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not: InputMaybe<ExpressionByTwoAnnotationFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ExpressionByTwoAnnotationFilter>>;
  /** Filter by the object’s `q3` field. */
  q3: InputMaybe<FloatFilter>;
  /** Filter by the object’s `secondAnnotationDisplayValue` field. */
  secondAnnotationDisplayValue: InputMaybe<StringFilter>;
  /** Filter by the object’s `secondAnnotationValueId` field. */
  secondAnnotationValueId: InputMaybe<IntFilter>;
  /** Filter by the object’s `valueCount` field. */
  valueCount: InputMaybe<IntFilter>;
  /** Filter by the object’s `values` field. */
  values: InputMaybe<FloatListFilter>;
};

/** A connection to a list of `ExpressionByTwoAnnotation` values. */
export type ExpressionByTwoAnnotationsConnection = {
  __typename?: 'ExpressionByTwoAnnotationsConnection';
  /** A list of edges which contains the `ExpressionByTwoAnnotation` and cursor to aid in pagination. */
  edges: Array<ExpressionByTwoAnnotationsEdge>;
  /** A list of `ExpressionByTwoAnnotation` objects. */
  nodes: Array<Maybe<ExpressionByTwoAnnotation>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ExpressionByTwoAnnotation` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ExpressionByTwoAnnotation` edge in the connection. */
export type ExpressionByTwoAnnotationsEdge = {
  __typename?: 'ExpressionByTwoAnnotationsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ExpressionByTwoAnnotation` at the end of the edge. */
  node: Maybe<ExpressionByTwoAnnotation>;
};

/** A filter to be used against Float fields. All fields are combined with a logical ‘and.’ */
export type FloatFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<Scalars['Float']>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<Scalars['Float']>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<Scalars['Float']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<Scalars['Float']>;
  /** Included in the specified list. */
  in: InputMaybe<Array<Scalars['Float']>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<Scalars['Float']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<Scalars['Float']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<Scalars['Float']>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<Scalars['Float']>;
  /** Not included in the specified list. */
  notIn: InputMaybe<Array<Scalars['Float']>>;
};

/** A filter to be used against Float List fields. All fields are combined with a logical ‘and.’ */
export type FloatListFilter = {
  /** Any array item is equal to the specified value. */
  anyEqualTo: InputMaybe<Scalars['Float']>;
  /** Any array item is greater than the specified value. */
  anyGreaterThan: InputMaybe<Scalars['Float']>;
  /** Any array item is greater than or equal to the specified value. */
  anyGreaterThanOrEqualTo: InputMaybe<Scalars['Float']>;
  /** Any array item is less than the specified value. */
  anyLessThan: InputMaybe<Scalars['Float']>;
  /** Any array item is less than or equal to the specified value. */
  anyLessThanOrEqualTo: InputMaybe<Scalars['Float']>;
  /** Any array item is not equal to the specified value. */
  anyNotEqualTo: InputMaybe<Scalars['Float']>;
  /** Contained by the specified list of values. */
  containedBy: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Contains the specified list of values. */
  contains: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Overlaps the specified list of values. */
  overlaps: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A `GetCorrelatedGenesRecord` edge in the connection. */
export type GetCorrelatedGeneEdge = {
  __typename?: 'GetCorrelatedGeneEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `GetCorrelatedGenesRecord` at the end of the edge. */
  node: Maybe<GetCorrelatedGenesRecord>;
};

/** A connection to a list of `GetCorrelatedGenesRecord` values. */
export type GetCorrelatedGenesConnection = {
  __typename?: 'GetCorrelatedGenesConnection';
  /** A list of edges which contains the `GetCorrelatedGenesRecord` and cursor to aid in pagination. */
  edges: Array<GetCorrelatedGeneEdge>;
  /** A list of `GetCorrelatedGenesRecord` objects. */
  nodes: Array<Maybe<GetCorrelatedGenesRecord>>;
  /** The count of *all* `GetCorrelatedGenesRecord` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** The return type of our `getCorrelatedGenes` query. */
export type GetCorrelatedGenesRecord = {
  __typename?: 'GetCorrelatedGenesRecord';
  displayName: Maybe<Scalars['String']>;
  displaySymbol: Maybe<Scalars['String']>;
  omicsId: Maybe<Scalars['Int']>;
  r: Maybe<Scalars['Float']>;
};

/** A filter to be used against `GetCorrelatedGenesRecord` object types. All fields are combined with a logical ‘and.’ */
export type GetCorrelatedGenesRecordFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<GetCorrelatedGenesRecordFilter>>;
  /** Filter by the object’s `displayName` field. */
  displayName: InputMaybe<StringFilter>;
  /** Filter by the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<GetCorrelatedGenesRecordFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<GetCorrelatedGenesRecordFilter>>;
  /** Filter by the object’s `r` field. */
  r: InputMaybe<FloatFilter>;
};

/** A filter to be used against Int fields. All fields are combined with a logical ‘and.’ */
export type IntFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<Scalars['Int']>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<Scalars['Int']>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<Scalars['Int']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<Scalars['Int']>;
  /** Included in the specified list. */
  in: InputMaybe<Array<Scalars['Int']>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<Scalars['Int']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<Scalars['Int']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<Scalars['Int']>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<Scalars['Int']>;
  /** Not included in the specified list. */
  notIn: InputMaybe<Array<Scalars['Int']>>;
};

/** A filter to be used against Int List fields. All fields are combined with a logical ‘and.’ */
export type IntListFilter = {
  /** Any array item is equal to the specified value. */
  anyEqualTo: InputMaybe<Scalars['Int']>;
  /** Any array item is greater than the specified value. */
  anyGreaterThan: InputMaybe<Scalars['Int']>;
  /** Any array item is greater than or equal to the specified value. */
  anyGreaterThanOrEqualTo: InputMaybe<Scalars['Int']>;
  /** Any array item is less than the specified value. */
  anyLessThan: InputMaybe<Scalars['Int']>;
  /** Any array item is less than or equal to the specified value. */
  anyLessThanOrEqualTo: InputMaybe<Scalars['Int']>;
  /** Any array item is not equal to the specified value. */
  anyNotEqualTo: InputMaybe<Scalars['Int']>;
  /** Contained by the specified list of values. */
  containedBy: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Contains the specified list of values. */
  contains: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Overlaps the specified list of values. */
  overlaps: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
};

/** A filter to be used against JSON fields. All fields are combined with a logical ‘and.’ */
export type JsonFilter = {
  /** Contained by the specified JSON. */
  containedBy: InputMaybe<Scalars['JSON']>;
  /** Contains the specified JSON. */
  contains: InputMaybe<Scalars['JSON']>;
  /** Contains all of the specified keys. */
  containsAllKeys: InputMaybe<Array<Scalars['String']>>;
  /** Contains any of the specified keys. */
  containsAnyKeys: InputMaybe<Array<Scalars['String']>>;
  /** Contains the specified key. */
  containsKey: InputMaybe<Scalars['String']>;
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<Scalars['JSON']>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<Scalars['JSON']>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<Scalars['JSON']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<Scalars['JSON']>;
  /** Included in the specified list. */
  in: InputMaybe<Array<Scalars['JSON']>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<Scalars['JSON']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<Scalars['JSON']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<Scalars['JSON']>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<Scalars['JSON']>;
  /** Not included in the specified list. */
  notIn: InputMaybe<Array<Scalars['JSON']>>;
};

export type MinimumTreesResult = {
  __typename?: 'MinimumTreesResult';
  conceptPaths: Maybe<Array<Maybe<ConceptPath>>>;
  conceptTreeElements: Maybe<Array<Maybe<ConceptTreeElement>>>;
};

/** The root mutation type which contains root level fields which mutate data. */
export type Mutation = {
  __typename?: 'Mutation';
  conceptAllChildrenPaths: Maybe<ConceptAllChildrenPathsPayload>;
  /** Creates a single `AnnotationGroup`. */
  createAnnotationGroup: Maybe<CreateAnnotationGroupPayload>;
  /** Creates a single `AnnotationValue`. */
  createAnnotationValue: Maybe<CreateAnnotationValuePayload>;
  /** Creates a single `ApiStudiesBase`. */
  createApiStudiesBase: Maybe<CreateApiStudiesBasePayload>;
  /** Creates a single `ApiStudiesBulkRna`. */
  createApiStudiesBulkRna: Maybe<CreateApiStudiesBulkRnaPayload>;
  /** Creates a single `ApiStudiesSingleCell`. */
  createApiStudiesSingleCell: Maybe<CreateApiStudiesSingleCellPayload>;
  /** Creates a single `ApiStudyH5Download`. */
  createApiStudyH5Download: Maybe<CreateApiStudyH5DownloadPayload>;
  /** Creates a single `Concept`. */
  createConcept: Maybe<CreateConceptPayload>;
  /** Creates a single `ConceptHierarchy`. */
  createConceptHierarchy: Maybe<CreateConceptHierarchyPayload>;
  /** Creates a single `ConceptSynonym`. */
  createConceptSynonym: Maybe<CreateConceptSynonymPayload>;
  /** Creates a single `DifferentialExpression`. */
  createDifferentialExpression: Maybe<CreateDifferentialExpressionPayload>;
  /** Creates a single `OmicsBase`. */
  createOmicsBase: Maybe<CreateOmicsBasePayload>;
  /** Creates a single `OmicsGene`. */
  createOmicsGene: Maybe<CreateOmicsGenePayload>;
  /** Creates a single `OmicsProteinAntibodyTag`. */
  createOmicsProteinAntibodyTag: Maybe<CreateOmicsProteinAntibodyTagPayload>;
  /** Creates a single `OmicsProteinAntibodyTagGene`. */
  createOmicsProteinAntibodyTagGene: Maybe<CreateOmicsProteinAntibodyTagGenePayload>;
  /** Creates a single `OmicsRegion`. */
  createOmicsRegion: Maybe<CreateOmicsRegionPayload>;
  /** Creates a single `OmicsRegionGene`. */
  createOmicsRegionGene: Maybe<CreateOmicsRegionGenePayload>;
  /** Creates a single `OmicsTranscriptionFactor`. */
  createOmicsTranscriptionFactor: Maybe<CreateOmicsTranscriptionFactorPayload>;
  /** Creates a single `OmicsTranscriptionFactorGene`. */
  createOmicsTranscriptionFactorGene: Maybe<CreateOmicsTranscriptionFactorGenePayload>;
  /** Creates a single `Ontology`. */
  createOntology: Maybe<CreateOntologyPayload>;
  /** Creates a single `ReferenceStudy`. */
  createReferenceStudy: Maybe<CreateReferenceStudyPayload>;
  createS3TempCredentials: Maybe<CreateS3TempCredentialsPayload>;
  /** Creates a single `Study`. */
  createStudy: Maybe<CreateStudyPayload>;
  /** Creates a single `StudyAdministrableCurrentuser`. */
  createStudyAdministrableCurrentuser: Maybe<CreateStudyAdministrableCurrentuserPayload>;
  /** Creates a single `StudyAnnotationGroupUi`. */
  createStudyAnnotationGroupUi: Maybe<CreateStudyAnnotationGroupUiPayload>;
  createStudyForCurrentUser: Maybe<CreateStudyForCurrentUserPayload>;
  createStudyH5AdPresignedUrl: Maybe<CreateStudyH5AdPresignedUrlPayload>;
  /** Creates a single `StudyLayer`. */
  createStudyLayer: Maybe<CreateStudyLayerPayload>;
  /** Creates a single `StudyOmic`. */
  createStudyOmic: Maybe<CreateStudyOmicPayload>;
  /** Creates a single `StudyOverview`. */
  createStudyOverview: Maybe<CreateStudyOverviewPayload>;
  /** Creates a single `StudySample`. */
  createStudySample: Maybe<CreateStudySamplePayload>;
  /** Creates a single `StudySampleAnnotation`. */
  createStudySampleAnnotation: Maybe<CreateStudySampleAnnotationPayload>;
  /** Creates a single `StudySampleProjection`. */
  createStudySampleProjection: Maybe<CreateStudySampleProjectionPayload>;
  createStudyUpload: Maybe<CreateStudyUploadPayload>;
  /** Creates a single `StudyVisibleCurrentuser`. */
  createStudyVisibleCurrentuser: Maybe<CreateStudyVisibleCurrentuserPayload>;
  /** Creates a single `UserAnnotationGroup`. */
  createUserAnnotationGroup: Maybe<CreateUserAnnotationGroupPayload>;
  deleteAllStudyData: Maybe<DeleteAllStudyDataPayload>;
  /** Deletes a single `AnnotationGroup` using a unique key. */
  deleteAnnotationGroup: Maybe<DeleteAnnotationGroupPayload>;
  /** Deletes a single `AnnotationGroup` using its globally unique id. */
  deleteAnnotationGroupByNodeId: Maybe<DeleteAnnotationGroupPayload>;
  /** Deletes a single `AnnotationValue` using a unique key. */
  deleteAnnotationValue: Maybe<DeleteAnnotationValuePayload>;
  /** Deletes a single `AnnotationValue` using its globally unique id. */
  deleteAnnotationValueByNodeId: Maybe<DeleteAnnotationValuePayload>;
  /** Deletes a single `Concept` using a unique key. */
  deleteConcept: Maybe<DeleteConceptPayload>;
  /** Deletes a single `Concept` using its globally unique id. */
  deleteConceptByNodeId: Maybe<DeleteConceptPayload>;
  /** Deletes a single `OmicsBase` using a unique key. */
  deleteOmicsBase: Maybe<DeleteOmicsBasePayload>;
  /** Deletes a single `OmicsBase` using its globally unique id. */
  deleteOmicsBaseByNodeId: Maybe<DeleteOmicsBasePayload>;
  /** Deletes a single `OmicsGene` using a unique key. */
  deleteOmicsGene: Maybe<DeleteOmicsGenePayload>;
  /** Deletes a single `OmicsGene` using its globally unique id. */
  deleteOmicsGeneByNodeId: Maybe<DeleteOmicsGenePayload>;
  /** Deletes a single `OmicsProteinAntibodyTag` using a unique key. */
  deleteOmicsProteinAntibodyTag: Maybe<DeleteOmicsProteinAntibodyTagPayload>;
  /** Deletes a single `OmicsProteinAntibodyTag` using its globally unique id. */
  deleteOmicsProteinAntibodyTagByNodeId: Maybe<DeleteOmicsProteinAntibodyTagPayload>;
  /** Deletes a single `OmicsRegion` using a unique key. */
  deleteOmicsRegion: Maybe<DeleteOmicsRegionPayload>;
  /** Deletes a single `OmicsRegion` using its globally unique id. */
  deleteOmicsRegionByNodeId: Maybe<DeleteOmicsRegionPayload>;
  /** Deletes a single `OmicsTranscriptionFactor` using a unique key. */
  deleteOmicsTranscriptionFactor: Maybe<DeleteOmicsTranscriptionFactorPayload>;
  /** Deletes a single `OmicsTranscriptionFactor` using its globally unique id. */
  deleteOmicsTranscriptionFactorByNodeId: Maybe<DeleteOmicsTranscriptionFactorPayload>;
  /** Deletes a single `Ontology` using a unique key. */
  deleteOntology: Maybe<DeleteOntologyPayload>;
  /** Deletes a single `Ontology` using its globally unique id. */
  deleteOntologyByNodeId: Maybe<DeleteOntologyPayload>;
  /** Deletes a single `Study` using a unique key. */
  deleteStudy: Maybe<DeleteStudyPayload>;
  /** Deletes a single `Study` using its globally unique id. */
  deleteStudyByNodeId: Maybe<DeleteStudyPayload>;
  /** Deletes a single `StudyLayer` using a unique key. */
  deleteStudyLayer: Maybe<DeleteStudyLayerPayload>;
  /** Deletes a single `StudyLayer` using its globally unique id. */
  deleteStudyLayerByNodeId: Maybe<DeleteStudyLayerPayload>;
  /** Deletes a single `StudySample` using a unique key. */
  deleteStudySample: Maybe<DeleteStudySamplePayload>;
  /** Deletes a single `StudySample` using its globally unique id. */
  deleteStudySampleByNodeId: Maybe<DeleteStudySamplePayload>;
  studyDefinitionUpdate: Maybe<StudyDefinitionUpdatePayload>;
  /** Updates a single `AnnotationGroup` using a unique key and a patch. */
  updateAnnotationGroup: Maybe<UpdateAnnotationGroupPayload>;
  /** Updates a single `AnnotationGroup` using its globally unique id and a patch. */
  updateAnnotationGroupByNodeId: Maybe<UpdateAnnotationGroupPayload>;
  /** Updates a single `AnnotationValue` using a unique key and a patch. */
  updateAnnotationValue: Maybe<UpdateAnnotationValuePayload>;
  /** Updates a single `AnnotationValue` using its globally unique id and a patch. */
  updateAnnotationValueByNodeId: Maybe<UpdateAnnotationValuePayload>;
  /** Updates a single `Concept` using a unique key and a patch. */
  updateConcept: Maybe<UpdateConceptPayload>;
  /** Updates a single `Concept` using its globally unique id and a patch. */
  updateConceptByNodeId: Maybe<UpdateConceptPayload>;
  /** Updates a single `OmicsBase` using a unique key and a patch. */
  updateOmicsBase: Maybe<UpdateOmicsBasePayload>;
  /** Updates a single `OmicsBase` using its globally unique id and a patch. */
  updateOmicsBaseByNodeId: Maybe<UpdateOmicsBasePayload>;
  /** Updates a single `OmicsGene` using a unique key and a patch. */
  updateOmicsGene: Maybe<UpdateOmicsGenePayload>;
  /** Updates a single `OmicsGene` using its globally unique id and a patch. */
  updateOmicsGeneByNodeId: Maybe<UpdateOmicsGenePayload>;
  /** Updates a single `OmicsProteinAntibodyTag` using a unique key and a patch. */
  updateOmicsProteinAntibodyTag: Maybe<UpdateOmicsProteinAntibodyTagPayload>;
  /** Updates a single `OmicsProteinAntibodyTag` using its globally unique id and a patch. */
  updateOmicsProteinAntibodyTagByNodeId: Maybe<UpdateOmicsProteinAntibodyTagPayload>;
  /** Updates a single `OmicsRegion` using a unique key and a patch. */
  updateOmicsRegion: Maybe<UpdateOmicsRegionPayload>;
  /** Updates a single `OmicsRegion` using its globally unique id and a patch. */
  updateOmicsRegionByNodeId: Maybe<UpdateOmicsRegionPayload>;
  /** Updates a single `OmicsTranscriptionFactor` using a unique key and a patch. */
  updateOmicsTranscriptionFactor: Maybe<UpdateOmicsTranscriptionFactorPayload>;
  /** Updates a single `OmicsTranscriptionFactor` using its globally unique id and a patch. */
  updateOmicsTranscriptionFactorByNodeId: Maybe<UpdateOmicsTranscriptionFactorPayload>;
  /** Updates a single `Ontology` using a unique key and a patch. */
  updateOntology: Maybe<UpdateOntologyPayload>;
  /** Updates a single `Ontology` using its globally unique id and a patch. */
  updateOntologyByNodeId: Maybe<UpdateOntologyPayload>;
  /** Updates a single `Study` using a unique key and a patch. */
  updateStudy: Maybe<UpdateStudyPayload>;
  /** Updates a single `Study` using its globally unique id and a patch. */
  updateStudyByNodeId: Maybe<UpdateStudyPayload>;
  /** Updates a single `StudyLayer` using a unique key and a patch. */
  updateStudyLayer: Maybe<UpdateStudyLayerPayload>;
  /** Updates a single `StudyLayer` using its globally unique id and a patch. */
  updateStudyLayerByNodeId: Maybe<UpdateStudyLayerPayload>;
  /** Updates a single `StudySample` using a unique key and a patch. */
  updateStudySample: Maybe<UpdateStudySamplePayload>;
  /** Updates a single `StudySample` using its globally unique id and a patch. */
  updateStudySampleByNodeId: Maybe<UpdateStudySamplePayload>;
  userAnnotationDefine: Maybe<UserAnnotationDefinePayload>;
  userAnnotationDelete: Maybe<UserAnnotationDeletePayload>;
  userAnnotationEdit: Maybe<UserAnnotationEditPayload>;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationConceptAllChildrenPathsArgs = {
  input: ConceptAllChildrenPathsInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateAnnotationGroupArgs = {
  input: CreateAnnotationGroupInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateAnnotationValueArgs = {
  input: CreateAnnotationValueInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateApiStudiesBaseArgs = {
  input: CreateApiStudiesBaseInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateApiStudiesBulkRnaArgs = {
  input: CreateApiStudiesBulkRnaInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateApiStudiesSingleCellArgs = {
  input: CreateApiStudiesSingleCellInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateApiStudyH5DownloadArgs = {
  input: CreateApiStudyH5DownloadInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateConceptArgs = {
  input: CreateConceptInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateConceptHierarchyArgs = {
  input: CreateConceptHierarchyInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateConceptSynonymArgs = {
  input: CreateConceptSynonymInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateDifferentialExpressionArgs = {
  input: CreateDifferentialExpressionInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateOmicsBaseArgs = {
  input: CreateOmicsBaseInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateOmicsGeneArgs = {
  input: CreateOmicsGeneInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateOmicsProteinAntibodyTagArgs = {
  input: CreateOmicsProteinAntibodyTagInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateOmicsProteinAntibodyTagGeneArgs = {
  input: CreateOmicsProteinAntibodyTagGeneInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateOmicsRegionArgs = {
  input: CreateOmicsRegionInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateOmicsRegionGeneArgs = {
  input: CreateOmicsRegionGeneInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateOmicsTranscriptionFactorArgs = {
  input: CreateOmicsTranscriptionFactorInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateOmicsTranscriptionFactorGeneArgs = {
  input: CreateOmicsTranscriptionFactorGeneInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateOntologyArgs = {
  input: CreateOntologyInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateReferenceStudyArgs = {
  input: CreateReferenceStudyInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateS3TempCredentialsArgs = {
  input: CreateS3TempCredentialsInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudyArgs = {
  input: CreateStudyInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudyAdministrableCurrentuserArgs = {
  input: CreateStudyAdministrableCurrentuserInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudyAnnotationGroupUiArgs = {
  input: CreateStudyAnnotationGroupUiInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudyForCurrentUserArgs = {
  input: CreateStudyForCurrentUserInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudyH5AdPresignedUrlArgs = {
  input: CreateStudyH5AdPresignedUrlInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudyLayerArgs = {
  input: CreateStudyLayerInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudyOmicArgs = {
  input: CreateStudyOmicInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudyOverviewArgs = {
  input: CreateStudyOverviewInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudySampleArgs = {
  input: CreateStudySampleInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudySampleAnnotationArgs = {
  input: CreateStudySampleAnnotationInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudySampleProjectionArgs = {
  input: CreateStudySampleProjectionInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudyUploadArgs = {
  input: CreateStudyUploadInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudyVisibleCurrentuserArgs = {
  input: CreateStudyVisibleCurrentuserInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateUserAnnotationGroupArgs = {
  input: CreateUserAnnotationGroupInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteAllStudyDataArgs = {
  input: DeleteAllStudyDataInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteAnnotationGroupArgs = {
  input: DeleteAnnotationGroupInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteAnnotationGroupByNodeIdArgs = {
  input: DeleteAnnotationGroupByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteAnnotationValueArgs = {
  input: DeleteAnnotationValueInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteAnnotationValueByNodeIdArgs = {
  input: DeleteAnnotationValueByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteConceptArgs = {
  input: DeleteConceptInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteConceptByNodeIdArgs = {
  input: DeleteConceptByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOmicsBaseArgs = {
  input: DeleteOmicsBaseInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOmicsBaseByNodeIdArgs = {
  input: DeleteOmicsBaseByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOmicsGeneArgs = {
  input: DeleteOmicsGeneInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOmicsGeneByNodeIdArgs = {
  input: DeleteOmicsGeneByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOmicsProteinAntibodyTagArgs = {
  input: DeleteOmicsProteinAntibodyTagInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOmicsProteinAntibodyTagByNodeIdArgs = {
  input: DeleteOmicsProteinAntibodyTagByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOmicsRegionArgs = {
  input: DeleteOmicsRegionInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOmicsRegionByNodeIdArgs = {
  input: DeleteOmicsRegionByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOmicsTranscriptionFactorArgs = {
  input: DeleteOmicsTranscriptionFactorInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOmicsTranscriptionFactorByNodeIdArgs = {
  input: DeleteOmicsTranscriptionFactorByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOntologyArgs = {
  input: DeleteOntologyInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOntologyByNodeIdArgs = {
  input: DeleteOntologyByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteStudyArgs = {
  input: DeleteStudyInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteStudyByNodeIdArgs = {
  input: DeleteStudyByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteStudyLayerArgs = {
  input: DeleteStudyLayerInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteStudyLayerByNodeIdArgs = {
  input: DeleteStudyLayerByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteStudySampleArgs = {
  input: DeleteStudySampleInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteStudySampleByNodeIdArgs = {
  input: DeleteStudySampleByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationStudyDefinitionUpdateArgs = {
  input: StudyDefinitionUpdateInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateAnnotationGroupArgs = {
  input: UpdateAnnotationGroupInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateAnnotationGroupByNodeIdArgs = {
  input: UpdateAnnotationGroupByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateAnnotationValueArgs = {
  input: UpdateAnnotationValueInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateAnnotationValueByNodeIdArgs = {
  input: UpdateAnnotationValueByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateConceptArgs = {
  input: UpdateConceptInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateConceptByNodeIdArgs = {
  input: UpdateConceptByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOmicsBaseArgs = {
  input: UpdateOmicsBaseInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOmicsBaseByNodeIdArgs = {
  input: UpdateOmicsBaseByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOmicsGeneArgs = {
  input: UpdateOmicsGeneInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOmicsGeneByNodeIdArgs = {
  input: UpdateOmicsGeneByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOmicsProteinAntibodyTagArgs = {
  input: UpdateOmicsProteinAntibodyTagInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOmicsProteinAntibodyTagByNodeIdArgs = {
  input: UpdateOmicsProteinAntibodyTagByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOmicsRegionArgs = {
  input: UpdateOmicsRegionInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOmicsRegionByNodeIdArgs = {
  input: UpdateOmicsRegionByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOmicsTranscriptionFactorArgs = {
  input: UpdateOmicsTranscriptionFactorInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOmicsTranscriptionFactorByNodeIdArgs = {
  input: UpdateOmicsTranscriptionFactorByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOntologyArgs = {
  input: UpdateOntologyInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOntologyByNodeIdArgs = {
  input: UpdateOntologyByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateStudyArgs = {
  input: UpdateStudyInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateStudyByNodeIdArgs = {
  input: UpdateStudyByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateStudyLayerArgs = {
  input: UpdateStudyLayerInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateStudyLayerByNodeIdArgs = {
  input: UpdateStudyLayerByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateStudySampleArgs = {
  input: UpdateStudySampleInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateStudySampleByNodeIdArgs = {
  input: UpdateStudySampleByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUserAnnotationDefineArgs = {
  input: UserAnnotationDefineInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUserAnnotationDeleteArgs = {
  input: UserAnnotationDeleteInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUserAnnotationEditArgs = {
  input: UserAnnotationEditInput;
};

/** An object with a globally unique `ID`. */
export type Node = {
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
};

export type OmicsAll = {
  __typename?: 'OmicsAll';
  displayName: Maybe<Scalars['String']>;
  displaySymbol: Maybe<Scalars['String']>;
  ensemblGeneId: Maybe<Scalars['String']>;
  entrezGeneIds: Maybe<Array<Maybe<Scalars['String']>>>;
  hgncSymbols: Maybe<Array<Maybe<Scalars['String']>>>;
  linkedGenes: Maybe<Array<Maybe<Scalars['Int']>>>;
  omicsId: Maybe<Scalars['Int']>;
  omicsType: Maybe<OmicsType>;
  region: Maybe<Scalars['String']>;
  taxId: Maybe<Scalars['Int']>;
};

/**
 * A condition to be used against `OmicsAll` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type OmicsAllCondition = {
  /** Checks for equality with the object’s `displayName` field. */
  displayName: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `ensemblGeneId` field. */
  ensemblGeneId: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `entrezGeneIds` field. */
  entrezGeneIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `hgncSymbols` field. */
  hgncSymbols: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `linkedGenes` field. */
  linkedGenes: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `omicsType` field. */
  omicsType: InputMaybe<OmicsType>;
  /** Checks for equality with the object’s `region` field. */
  region: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `taxId` field. */
  taxId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `OmicsAll` object types. All fields are combined with a logical ‘and.’ */
export type OmicsAllFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<OmicsAllFilter>>;
  /** Filter by the object’s `displayName` field. */
  displayName: InputMaybe<StringFilter>;
  /** Filter by the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<StringFilter>;
  /** Filter by the object’s `ensemblGeneId` field. */
  ensemblGeneId: InputMaybe<StringFilter>;
  /** Filter by the object’s `entrezGeneIds` field. */
  entrezGeneIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `hgncSymbols` field. */
  hgncSymbols: InputMaybe<StringListFilter>;
  /** Filter by the object’s `linkedGenes` field. */
  linkedGenes: InputMaybe<IntListFilter>;
  /** Negates the expression. */
  not: InputMaybe<OmicsAllFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Filter by the object’s `omicsType` field. */
  omicsType: InputMaybe<OmicsTypeFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<OmicsAllFilter>>;
  /** Filter by the object’s `region` field. */
  region: InputMaybe<StringFilter>;
  /** Filter by the object’s `taxId` field. */
  taxId: InputMaybe<IntFilter>;
};

/** A connection to a list of `OmicsAll` values. */
export type OmicsAllsConnection = {
  __typename?: 'OmicsAllsConnection';
  /** A list of edges which contains the `OmicsAll` and cursor to aid in pagination. */
  edges: Array<OmicsAllsEdge>;
  /** A list of `OmicsAll` objects. */
  nodes: Array<Maybe<OmicsAll>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `OmicsAll` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `OmicsAll` edge in the connection. */
export type OmicsAllsEdge = {
  __typename?: 'OmicsAllsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `OmicsAll` at the end of the edge. */
  node: Maybe<OmicsAll>;
};

/** Methods to use when ordering `OmicsAll`. */
export enum OmicsAllsOrderBy {
  DisplayNameAsc = 'DISPLAY_NAME_ASC',
  DisplayNameDesc = 'DISPLAY_NAME_DESC',
  DisplaySymbolAsc = 'DISPLAY_SYMBOL_ASC',
  DisplaySymbolDesc = 'DISPLAY_SYMBOL_DESC',
  EnsemblGeneIdAsc = 'ENSEMBL_GENE_ID_ASC',
  EnsemblGeneIdDesc = 'ENSEMBL_GENE_ID_DESC',
  EntrezGeneIdsAsc = 'ENTREZ_GENE_IDS_ASC',
  EntrezGeneIdsDesc = 'ENTREZ_GENE_IDS_DESC',
  HgncSymbolsAsc = 'HGNC_SYMBOLS_ASC',
  HgncSymbolsDesc = 'HGNC_SYMBOLS_DESC',
  LinkedGenesAsc = 'LINKED_GENES_ASC',
  LinkedGenesDesc = 'LINKED_GENES_DESC',
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  OmicsTypeAsc = 'OMICS_TYPE_ASC',
  OmicsTypeDesc = 'OMICS_TYPE_DESC',
  RegionAsc = 'REGION_ASC',
  RegionDesc = 'REGION_DESC',
  TaxIdAsc = 'TAX_ID_ASC',
  TaxIdDesc = 'TAX_ID_DESC'
}

export type OmicsAutocompleteResult = {
  __typename?: 'OmicsAutocompleteResult';
  displaySymbol: Maybe<Scalars['String']>;
  labelHighlight: Maybe<Scalars['String']>;
  omicsId: Maybe<Array<Maybe<Scalars['Int']>>>;
  omicsType: Maybe<Array<Maybe<OmicsType>>>;
};

/** A filter to be used against `OmicsAutocompleteResult` object types. All fields are combined with a logical ‘and.’ */
export type OmicsAutocompleteResultFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<OmicsAutocompleteResultFilter>>;
  /** Filter by the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<StringFilter>;
  /** Filter by the object’s `labelHighlight` field. */
  labelHighlight: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<OmicsAutocompleteResultFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntListFilter>;
  /** Filter by the object’s `omicsType` field. */
  omicsType: InputMaybe<OmicsTypeListFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<OmicsAutocompleteResultFilter>>;
};

/** A connection to a list of `OmicsAutocompleteResult` values. */
export type OmicsAutocompleteResultsConnection = {
  __typename?: 'OmicsAutocompleteResultsConnection';
  /** A list of edges which contains the `OmicsAutocompleteResult` and cursor to aid in pagination. */
  edges: Array<OmicsAutocompleteResultsEdge>;
  /** A list of `OmicsAutocompleteResult` objects. */
  nodes: Array<Maybe<OmicsAutocompleteResult>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `OmicsAutocompleteResult` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `OmicsAutocompleteResult` edge in the connection. */
export type OmicsAutocompleteResultsEdge = {
  __typename?: 'OmicsAutocompleteResultsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `OmicsAutocompleteResult` at the end of the edge. */
  node: Maybe<OmicsAutocompleteResult>;
};

export type OmicsBase = Node & {
  __typename?: 'OmicsBase';
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressionsByOmicsId: DifferentialExpressionsConnection;
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressionsByOmicsIdList: Array<DifferentialExpression>;
  displayName: Maybe<Scalars['String']>;
  displayNameTsvector: Maybe<Scalars['String']>;
  displaySymbol: Scalars['String'];
  displaySymbolTsvector: Maybe<Scalars['String']>;
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads a single `OmicsGene` that is related to this `OmicsBase`. */
  omicsGeneByGeneId: Maybe<OmicsGene>;
  /**
   * Reads and enables pagination through a set of `OmicsGene`.
   * @deprecated Please use omicsGeneByGeneId instead
   */
  omicsGenesByGeneId: OmicsGenesConnection;
  omicsId: Scalars['Int'];
  /** Reads a single `OmicsProteinAntibodyTag` that is related to this `OmicsBase`. */
  omicsProteinAntibodyTagByProteinAntibodyTagId: Maybe<OmicsProteinAntibodyTag>;
  /**
   * Reads and enables pagination through a set of `OmicsProteinAntibodyTag`.
   * @deprecated Please use omicsProteinAntibodyTagByProteinAntibodyTagId instead
   */
  omicsProteinAntibodyTagsByProteinAntibodyTagId: OmicsProteinAntibodyTagsConnection;
  /**
   * Reads and enables pagination through a set of `OmicsRegion`.
   * @deprecated Please use omics_region_newNameHere instead
   */
  omicsRegionsByRegionId: OmicsRegionsConnection;
  /** Reads a single `OmicsTranscriptionFactor` that is related to this `OmicsBase`. */
  omicsTranscriptionFactorByOmicsId: Maybe<OmicsTranscriptionFactor>;
  /**
   * Reads and enables pagination through a set of `OmicsTranscriptionFactor`.
   * @deprecated Please use omicsTranscriptionFactorByOmicsId instead
   */
  omicsTranscriptionFactorsByOmicsId: OmicsTranscriptionFactorsConnection;
  omicsType: OmicsType;
  /** Reads a single `OmicsRegion` that is related to this `OmicsBase`. */
  omics_region_newNameHere: Maybe<OmicsRegion>;
  /** Reads and enables pagination through a set of `StudyOmic`. */
  studyOmicsByOmicsId: StudyOmicsConnection;
  /** Reads and enables pagination through a set of `StudyOmic`. */
  studyOmicsByOmicsIdList: Array<StudyOmic>;
  taxId: Scalars['Int'];
};


export type OmicsBaseDifferentialExpressionsByOmicsIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<DifferentialExpressionCondition>;
  filter: InputMaybe<DifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


export type OmicsBaseDifferentialExpressionsByOmicsIdListArgs = {
  condition: InputMaybe<DifferentialExpressionCondition>;
  filter: InputMaybe<DifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


export type OmicsBaseOmicsGenesByGeneIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsGeneCondition>;
  filter: InputMaybe<OmicsGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsGenesOrderBy>>;
};


export type OmicsBaseOmicsProteinAntibodyTagsByProteinAntibodyTagIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsProteinAntibodyTagCondition>;
  filter: InputMaybe<OmicsProteinAntibodyTagFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsProteinAntibodyTagsOrderBy>>;
};


export type OmicsBaseOmicsRegionsByRegionIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsRegionCondition>;
  filter: InputMaybe<OmicsRegionFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsRegionsOrderBy>>;
};


export type OmicsBaseOmicsTranscriptionFactorsByOmicsIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsTranscriptionFactorCondition>;
  filter: InputMaybe<OmicsTranscriptionFactorFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsTranscriptionFactorsOrderBy>>;
};


export type OmicsBaseStudyOmicsByOmicsIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyOmicCondition>;
  filter: InputMaybe<StudyOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyOmicsOrderBy>>;
};


export type OmicsBaseStudyOmicsByOmicsIdListArgs = {
  condition: InputMaybe<StudyOmicCondition>;
  filter: InputMaybe<StudyOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyOmicsOrderBy>>;
};

/**
 * A condition to be used against `OmicsBase` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type OmicsBaseCondition = {
  /** Checks for equality with the object’s `displayName` field. */
  displayName: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `displayNameTsvector` field. */
  displayNameTsvector: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `displaySymbolTsvector` field. */
  displaySymbolTsvector: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `omicsType` field. */
  omicsType: InputMaybe<OmicsType>;
  /** Checks for equality with the object’s `taxId` field. */
  taxId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `OmicsBase` object types. All fields are combined with a logical ‘and.’ */
export type OmicsBaseFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<OmicsBaseFilter>>;
  /** Filter by the object’s `displayName` field. */
  displayName: InputMaybe<StringFilter>;
  /** Filter by the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<OmicsBaseFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Filter by the object’s `omicsType` field. */
  omicsType: InputMaybe<OmicsTypeFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<OmicsBaseFilter>>;
  /** Filter by the object’s `taxId` field. */
  taxId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `OmicsBase` */
export type OmicsBaseInput = {
  displayName: InputMaybe<Scalars['String']>;
  displayNameTsvector: InputMaybe<Scalars['String']>;
  displaySymbol: Scalars['String'];
  displaySymbolTsvector: InputMaybe<Scalars['String']>;
  omicsId: InputMaybe<Scalars['Int']>;
  omicsType: OmicsType;
  taxId: Scalars['Int'];
};

/** Represents an update to a `OmicsBase`. Fields that are set will be updated. */
export type OmicsBasePatch = {
  displayName: InputMaybe<Scalars['String']>;
  displayNameTsvector: InputMaybe<Scalars['String']>;
  displaySymbol: InputMaybe<Scalars['String']>;
  displaySymbolTsvector: InputMaybe<Scalars['String']>;
  omicsId: InputMaybe<Scalars['Int']>;
  omicsType: InputMaybe<OmicsType>;
  taxId: InputMaybe<Scalars['Int']>;
};

/** A connection to a list of `OmicsBase` values. */
export type OmicsBasesConnection = {
  __typename?: 'OmicsBasesConnection';
  /** A list of edges which contains the `OmicsBase` and cursor to aid in pagination. */
  edges: Array<OmicsBasesEdge>;
  /** A list of `OmicsBase` objects. */
  nodes: Array<Maybe<OmicsBase>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `OmicsBase` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `OmicsBase` edge in the connection. */
export type OmicsBasesEdge = {
  __typename?: 'OmicsBasesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `OmicsBase` at the end of the edge. */
  node: Maybe<OmicsBase>;
};

/** Methods to use when ordering `OmicsBase`. */
export enum OmicsBasesOrderBy {
  DisplayNameAsc = 'DISPLAY_NAME_ASC',
  DisplayNameDesc = 'DISPLAY_NAME_DESC',
  DisplayNameTsvectorAsc = 'DISPLAY_NAME_TSVECTOR_ASC',
  DisplayNameTsvectorDesc = 'DISPLAY_NAME_TSVECTOR_DESC',
  DisplaySymbolAsc = 'DISPLAY_SYMBOL_ASC',
  DisplaySymbolDesc = 'DISPLAY_SYMBOL_DESC',
  DisplaySymbolTsvectorAsc = 'DISPLAY_SYMBOL_TSVECTOR_ASC',
  DisplaySymbolTsvectorDesc = 'DISPLAY_SYMBOL_TSVECTOR_DESC',
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  OmicsTypeAsc = 'OMICS_TYPE_ASC',
  OmicsTypeDesc = 'OMICS_TYPE_DESC',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  TaxIdAsc = 'TAX_ID_ASC',
  TaxIdDesc = 'TAX_ID_DESC'
}

export type OmicsGene = Node & {
  __typename?: 'OmicsGene';
  ensemblGeneId: Scalars['String'];
  entrezGeneIds: Maybe<Array<Maybe<Scalars['String']>>>;
  /** Reads a single `OmicsBase` that is related to this `OmicsGene`. */
  gene: Maybe<OmicsBase>;
  geneId: Scalars['Int'];
  hgncSymbols: Maybe<Array<Maybe<Scalars['String']>>>;
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads and enables pagination through a set of `OmicsProteinAntibodyTagGene`. */
  omicsProteinAntibodyTagGenesByGeneId: OmicsProteinAntibodyTagGenesConnection;
  /** Reads and enables pagination through a set of `OmicsProteinAntibodyTagGene`. */
  omicsProteinAntibodyTagGenesByGeneIdList: Array<OmicsProteinAntibodyTagGene>;
  /** Reads and enables pagination through a set of `OmicsRegionGene`. */
  omicsRegionGenesByGeneId: OmicsRegionGenesConnection;
  /** Reads and enables pagination through a set of `OmicsRegionGene`. */
  omicsRegionGenesByGeneIdList: Array<OmicsRegionGene>;
  /** Reads and enables pagination through a set of `OmicsTranscriptionFactorGene`. */
  omicsTranscriptionFactorGenesByGeneId: OmicsTranscriptionFactorGenesConnection;
  /** Reads and enables pagination through a set of `OmicsTranscriptionFactorGene`. */
  omicsTranscriptionFactorGenesByGeneIdList: Array<OmicsTranscriptionFactorGene>;
};


export type OmicsGeneOmicsProteinAntibodyTagGenesByGeneIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsProteinAntibodyTagGeneCondition>;
  filter: InputMaybe<OmicsProteinAntibodyTagGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsProteinAntibodyTagGenesOrderBy>>;
};


export type OmicsGeneOmicsProteinAntibodyTagGenesByGeneIdListArgs = {
  condition: InputMaybe<OmicsProteinAntibodyTagGeneCondition>;
  filter: InputMaybe<OmicsProteinAntibodyTagGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsProteinAntibodyTagGenesOrderBy>>;
};


export type OmicsGeneOmicsRegionGenesByGeneIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsRegionGeneCondition>;
  filter: InputMaybe<OmicsRegionGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsRegionGenesOrderBy>>;
};


export type OmicsGeneOmicsRegionGenesByGeneIdListArgs = {
  condition: InputMaybe<OmicsRegionGeneCondition>;
  filter: InputMaybe<OmicsRegionGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsRegionGenesOrderBy>>;
};


export type OmicsGeneOmicsTranscriptionFactorGenesByGeneIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsTranscriptionFactorGeneCondition>;
  filter: InputMaybe<OmicsTranscriptionFactorGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsTranscriptionFactorGenesOrderBy>>;
};


export type OmicsGeneOmicsTranscriptionFactorGenesByGeneIdListArgs = {
  condition: InputMaybe<OmicsTranscriptionFactorGeneCondition>;
  filter: InputMaybe<OmicsTranscriptionFactorGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsTranscriptionFactorGenesOrderBy>>;
};

/**
 * A condition to be used against `OmicsGene` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type OmicsGeneCondition = {
  /** Checks for equality with the object’s `ensemblGeneId` field. */
  ensemblGeneId: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `entrezGeneIds` field. */
  entrezGeneIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `geneId` field. */
  geneId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `hgncSymbols` field. */
  hgncSymbols: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** A filter to be used against `OmicsGene` object types. All fields are combined with a logical ‘and.’ */
export type OmicsGeneFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<OmicsGeneFilter>>;
  /** Filter by the object’s `ensemblGeneId` field. */
  ensemblGeneId: InputMaybe<StringFilter>;
  /** Filter by the object’s `entrezGeneIds` field. */
  entrezGeneIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `geneId` field. */
  geneId: InputMaybe<IntFilter>;
  /** Filter by the object’s `hgncSymbols` field. */
  hgncSymbols: InputMaybe<StringListFilter>;
  /** Negates the expression. */
  not: InputMaybe<OmicsGeneFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<OmicsGeneFilter>>;
};

/** An input for mutations affecting `OmicsGene` */
export type OmicsGeneInput = {
  ensemblGeneId: Scalars['String'];
  entrezGeneIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  geneId: Scalars['Int'];
  hgncSymbols: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** Represents an update to a `OmicsGene`. Fields that are set will be updated. */
export type OmicsGenePatch = {
  ensemblGeneId: InputMaybe<Scalars['String']>;
  entrezGeneIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  geneId: InputMaybe<Scalars['Int']>;
  hgncSymbols: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** A connection to a list of `OmicsGene` values. */
export type OmicsGenesConnection = {
  __typename?: 'OmicsGenesConnection';
  /** A list of edges which contains the `OmicsGene` and cursor to aid in pagination. */
  edges: Array<OmicsGenesEdge>;
  /** A list of `OmicsGene` objects. */
  nodes: Array<Maybe<OmicsGene>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `OmicsGene` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `OmicsGene` edge in the connection. */
export type OmicsGenesEdge = {
  __typename?: 'OmicsGenesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `OmicsGene` at the end of the edge. */
  node: Maybe<OmicsGene>;
};

/** Methods to use when ordering `OmicsGene`. */
export enum OmicsGenesOrderBy {
  EnsemblGeneIdAsc = 'ENSEMBL_GENE_ID_ASC',
  EnsemblGeneIdDesc = 'ENSEMBL_GENE_ID_DESC',
  EntrezGeneIdsAsc = 'ENTREZ_GENE_IDS_ASC',
  EntrezGeneIdsDesc = 'ENTREZ_GENE_IDS_DESC',
  GeneIdAsc = 'GENE_ID_ASC',
  GeneIdDesc = 'GENE_ID_DESC',
  HgncSymbolsAsc = 'HGNC_SYMBOLS_ASC',
  HgncSymbolsDesc = 'HGNC_SYMBOLS_DESC',
  Natural = 'NATURAL',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC'
}

export type OmicsProteinAntibodyTag = Node & {
  __typename?: 'OmicsProteinAntibodyTag';
  antibodySymbol: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads and enables pagination through a set of `OmicsProteinAntibodyTagGene`. */
  omicsProteinAntibodyTagGenesByProteinAntibodyTagId: OmicsProteinAntibodyTagGenesConnection;
  /** Reads and enables pagination through a set of `OmicsProteinAntibodyTagGene`. */
  omicsProteinAntibodyTagGenesByProteinAntibodyTagIdList: Array<OmicsProteinAntibodyTagGene>;
  /** Reads a single `OmicsBase` that is related to this `OmicsProteinAntibodyTag`. */
  proteinAntibodyTag: Maybe<OmicsBase>;
  proteinAntibodyTagId: Scalars['Int'];
  taxId: Scalars['Int'];
};


export type OmicsProteinAntibodyTagOmicsProteinAntibodyTagGenesByProteinAntibodyTagIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsProteinAntibodyTagGeneCondition>;
  filter: InputMaybe<OmicsProteinAntibodyTagGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsProteinAntibodyTagGenesOrderBy>>;
};


export type OmicsProteinAntibodyTagOmicsProteinAntibodyTagGenesByProteinAntibodyTagIdListArgs = {
  condition: InputMaybe<OmicsProteinAntibodyTagGeneCondition>;
  filter: InputMaybe<OmicsProteinAntibodyTagGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsProteinAntibodyTagGenesOrderBy>>;
};

/**
 * A condition to be used against `OmicsProteinAntibodyTag` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type OmicsProteinAntibodyTagCondition = {
  /** Checks for equality with the object’s `antibodySymbol` field. */
  antibodySymbol: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `proteinAntibodyTagId` field. */
  proteinAntibodyTagId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `taxId` field. */
  taxId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `OmicsProteinAntibodyTag` object types. All fields are combined with a logical ‘and.’ */
export type OmicsProteinAntibodyTagFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<OmicsProteinAntibodyTagFilter>>;
  /** Filter by the object’s `antibodySymbol` field. */
  antibodySymbol: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<OmicsProteinAntibodyTagFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<OmicsProteinAntibodyTagFilter>>;
  /** Filter by the object’s `proteinAntibodyTagId` field. */
  proteinAntibodyTagId: InputMaybe<IntFilter>;
  /** Filter by the object’s `taxId` field. */
  taxId: InputMaybe<IntFilter>;
};

export type OmicsProteinAntibodyTagGene = {
  __typename?: 'OmicsProteinAntibodyTagGene';
  /** Reads a single `OmicsGene` that is related to this `OmicsProteinAntibodyTagGene`. */
  gene: Maybe<OmicsGene>;
  geneId: Scalars['Int'];
  /** Reads a single `OmicsProteinAntibodyTag` that is related to this `OmicsProteinAntibodyTagGene`. */
  proteinAntibodyTag: Maybe<OmicsProteinAntibodyTag>;
  proteinAntibodyTagId: Scalars['Int'];
};

/**
 * A condition to be used against `OmicsProteinAntibodyTagGene` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type OmicsProteinAntibodyTagGeneCondition = {
  /** Checks for equality with the object’s `geneId` field. */
  geneId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `proteinAntibodyTagId` field. */
  proteinAntibodyTagId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `OmicsProteinAntibodyTagGene` object types. All fields are combined with a logical ‘and.’ */
export type OmicsProteinAntibodyTagGeneFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<OmicsProteinAntibodyTagGeneFilter>>;
  /** Filter by the object’s `geneId` field. */
  geneId: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not: InputMaybe<OmicsProteinAntibodyTagGeneFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<OmicsProteinAntibodyTagGeneFilter>>;
  /** Filter by the object’s `proteinAntibodyTagId` field. */
  proteinAntibodyTagId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `OmicsProteinAntibodyTagGene` */
export type OmicsProteinAntibodyTagGeneInput = {
  geneId: Scalars['Int'];
  proteinAntibodyTagId: Scalars['Int'];
};

/** A connection to a list of `OmicsProteinAntibodyTagGene` values. */
export type OmicsProteinAntibodyTagGenesConnection = {
  __typename?: 'OmicsProteinAntibodyTagGenesConnection';
  /** A list of edges which contains the `OmicsProteinAntibodyTagGene` and cursor to aid in pagination. */
  edges: Array<OmicsProteinAntibodyTagGenesEdge>;
  /** A list of `OmicsProteinAntibodyTagGene` objects. */
  nodes: Array<Maybe<OmicsProteinAntibodyTagGene>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `OmicsProteinAntibodyTagGene` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `OmicsProteinAntibodyTagGene` edge in the connection. */
export type OmicsProteinAntibodyTagGenesEdge = {
  __typename?: 'OmicsProteinAntibodyTagGenesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `OmicsProteinAntibodyTagGene` at the end of the edge. */
  node: Maybe<OmicsProteinAntibodyTagGene>;
};

/** Methods to use when ordering `OmicsProteinAntibodyTagGene`. */
export enum OmicsProteinAntibodyTagGenesOrderBy {
  GeneIdAsc = 'GENE_ID_ASC',
  GeneIdDesc = 'GENE_ID_DESC',
  Natural = 'NATURAL',
  ProteinAntibodyTagIdAsc = 'PROTEIN_ANTIBODY_TAG_ID_ASC',
  ProteinAntibodyTagIdDesc = 'PROTEIN_ANTIBODY_TAG_ID_DESC'
}

/** An input for mutations affecting `OmicsProteinAntibodyTag` */
export type OmicsProteinAntibodyTagInput = {
  antibodySymbol: Scalars['String'];
  proteinAntibodyTagId: Scalars['Int'];
  taxId: Scalars['Int'];
};

/** Represents an update to a `OmicsProteinAntibodyTag`. Fields that are set will be updated. */
export type OmicsProteinAntibodyTagPatch = {
  antibodySymbol: InputMaybe<Scalars['String']>;
  proteinAntibodyTagId: InputMaybe<Scalars['Int']>;
  taxId: InputMaybe<Scalars['Int']>;
};

/** A connection to a list of `OmicsProteinAntibodyTag` values. */
export type OmicsProteinAntibodyTagsConnection = {
  __typename?: 'OmicsProteinAntibodyTagsConnection';
  /** A list of edges which contains the `OmicsProteinAntibodyTag` and cursor to aid in pagination. */
  edges: Array<OmicsProteinAntibodyTagsEdge>;
  /** A list of `OmicsProteinAntibodyTag` objects. */
  nodes: Array<Maybe<OmicsProteinAntibodyTag>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `OmicsProteinAntibodyTag` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `OmicsProteinAntibodyTag` edge in the connection. */
export type OmicsProteinAntibodyTagsEdge = {
  __typename?: 'OmicsProteinAntibodyTagsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `OmicsProteinAntibodyTag` at the end of the edge. */
  node: Maybe<OmicsProteinAntibodyTag>;
};

/** Methods to use when ordering `OmicsProteinAntibodyTag`. */
export enum OmicsProteinAntibodyTagsOrderBy {
  AntibodySymbolAsc = 'ANTIBODY_SYMBOL_ASC',
  AntibodySymbolDesc = 'ANTIBODY_SYMBOL_DESC',
  Natural = 'NATURAL',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  ProteinAntibodyTagIdAsc = 'PROTEIN_ANTIBODY_TAG_ID_ASC',
  ProteinAntibodyTagIdDesc = 'PROTEIN_ANTIBODY_TAG_ID_DESC',
  TaxIdAsc = 'TAX_ID_ASC',
  TaxIdDesc = 'TAX_ID_DESC'
}

export type OmicsRegion = Node & {
  __typename?: 'OmicsRegion';
  chromosome: Scalars['String'];
  endPosition: Scalars['Int'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads and enables pagination through a set of `OmicsRegionGene`. */
  omicsRegionGenesByRegionId: OmicsRegionGenesConnection;
  /** Reads and enables pagination through a set of `OmicsRegionGene`. */
  omicsRegionGenesByRegionIdList: Array<OmicsRegionGene>;
  /** Reads a single `OmicsBase` that is related to this `OmicsRegion`. */
  omics_region_newNameHere: Maybe<OmicsBase>;
  region: Scalars['String'];
  regionId: Scalars['Int'];
  startPosition: Scalars['Int'];
};


export type OmicsRegionOmicsRegionGenesByRegionIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsRegionGeneCondition>;
  filter: InputMaybe<OmicsRegionGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsRegionGenesOrderBy>>;
};


export type OmicsRegionOmicsRegionGenesByRegionIdListArgs = {
  condition: InputMaybe<OmicsRegionGeneCondition>;
  filter: InputMaybe<OmicsRegionGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsRegionGenesOrderBy>>;
};

/**
 * A condition to be used against `OmicsRegion` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type OmicsRegionCondition = {
  /** Checks for equality with the object’s `chromosome` field. */
  chromosome: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `endPosition` field. */
  endPosition: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `region` field. */
  region: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `regionId` field. */
  regionId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `startPosition` field. */
  startPosition: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `OmicsRegion` object types. All fields are combined with a logical ‘and.’ */
export type OmicsRegionFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<OmicsRegionFilter>>;
  /** Filter by the object’s `chromosome` field. */
  chromosome: InputMaybe<StringFilter>;
  /** Filter by the object’s `endPosition` field. */
  endPosition: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not: InputMaybe<OmicsRegionFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<OmicsRegionFilter>>;
  /** Filter by the object’s `region` field. */
  region: InputMaybe<StringFilter>;
  /** Filter by the object’s `regionId` field. */
  regionId: InputMaybe<IntFilter>;
  /** Filter by the object’s `startPosition` field. */
  startPosition: InputMaybe<IntFilter>;
};

export type OmicsRegionGene = {
  __typename?: 'OmicsRegionGene';
  /** Reads a single `OmicsGene` that is related to this `OmicsRegionGene`. */
  gene: Maybe<OmicsGene>;
  geneId: Scalars['Int'];
  /** Reads a single `OmicsRegion` that is related to this `OmicsRegionGene`. */
  region: Maybe<OmicsRegion>;
  regionId: Scalars['Int'];
};

/**
 * A condition to be used against `OmicsRegionGene` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type OmicsRegionGeneCondition = {
  /** Checks for equality with the object’s `geneId` field. */
  geneId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `regionId` field. */
  regionId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `OmicsRegionGene` object types. All fields are combined with a logical ‘and.’ */
export type OmicsRegionGeneFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<OmicsRegionGeneFilter>>;
  /** Filter by the object’s `geneId` field. */
  geneId: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not: InputMaybe<OmicsRegionGeneFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<OmicsRegionGeneFilter>>;
  /** Filter by the object’s `regionId` field. */
  regionId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `OmicsRegionGene` */
export type OmicsRegionGeneInput = {
  geneId: Scalars['Int'];
  regionId: Scalars['Int'];
};

/** A connection to a list of `OmicsRegionGene` values. */
export type OmicsRegionGenesConnection = {
  __typename?: 'OmicsRegionGenesConnection';
  /** A list of edges which contains the `OmicsRegionGene` and cursor to aid in pagination. */
  edges: Array<OmicsRegionGenesEdge>;
  /** A list of `OmicsRegionGene` objects. */
  nodes: Array<Maybe<OmicsRegionGene>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `OmicsRegionGene` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `OmicsRegionGene` edge in the connection. */
export type OmicsRegionGenesEdge = {
  __typename?: 'OmicsRegionGenesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `OmicsRegionGene` at the end of the edge. */
  node: Maybe<OmicsRegionGene>;
};

/** Methods to use when ordering `OmicsRegionGene`. */
export enum OmicsRegionGenesOrderBy {
  GeneIdAsc = 'GENE_ID_ASC',
  GeneIdDesc = 'GENE_ID_DESC',
  Natural = 'NATURAL',
  RegionIdAsc = 'REGION_ID_ASC',
  RegionIdDesc = 'REGION_ID_DESC'
}

/** An input for mutations affecting `OmicsRegion` */
export type OmicsRegionInput = {
  chromosome: Scalars['String'];
  endPosition: Scalars['Int'];
  region: Scalars['String'];
  regionId: Scalars['Int'];
  startPosition: Scalars['Int'];
};

/** Represents an update to a `OmicsRegion`. Fields that are set will be updated. */
export type OmicsRegionPatch = {
  chromosome: InputMaybe<Scalars['String']>;
  endPosition: InputMaybe<Scalars['Int']>;
  region: InputMaybe<Scalars['String']>;
  regionId: InputMaybe<Scalars['Int']>;
  startPosition: InputMaybe<Scalars['Int']>;
};

/** A connection to a list of `OmicsRegion` values. */
export type OmicsRegionsConnection = {
  __typename?: 'OmicsRegionsConnection';
  /** A list of edges which contains the `OmicsRegion` and cursor to aid in pagination. */
  edges: Array<OmicsRegionsEdge>;
  /** A list of `OmicsRegion` objects. */
  nodes: Array<Maybe<OmicsRegion>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `OmicsRegion` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `OmicsRegion` edge in the connection. */
export type OmicsRegionsEdge = {
  __typename?: 'OmicsRegionsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `OmicsRegion` at the end of the edge. */
  node: Maybe<OmicsRegion>;
};

/** Methods to use when ordering `OmicsRegion`. */
export enum OmicsRegionsOrderBy {
  ChromosomeAsc = 'CHROMOSOME_ASC',
  ChromosomeDesc = 'CHROMOSOME_DESC',
  EndPositionAsc = 'END_POSITION_ASC',
  EndPositionDesc = 'END_POSITION_DESC',
  Natural = 'NATURAL',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  RegionAsc = 'REGION_ASC',
  RegionDesc = 'REGION_DESC',
  RegionIdAsc = 'REGION_ID_ASC',
  RegionIdDesc = 'REGION_ID_DESC',
  StartPositionAsc = 'START_POSITION_ASC',
  StartPositionDesc = 'START_POSITION_DESC'
}

export type OmicsTranscriptionFactor = Node & {
  __typename?: 'OmicsTranscriptionFactor';
  jasparMatrixId: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads a single `OmicsBase` that is related to this `OmicsTranscriptionFactor`. */
  omics: Maybe<OmicsBase>;
  omicsId: Scalars['Int'];
  /** Reads and enables pagination through a set of `OmicsTranscriptionFactorGene`. */
  omicsTranscriptionFactorGenesByTranscriptionFactorId: OmicsTranscriptionFactorGenesConnection;
  /** Reads and enables pagination through a set of `OmicsTranscriptionFactorGene`. */
  omicsTranscriptionFactorGenesByTranscriptionFactorIdList: Array<OmicsTranscriptionFactorGene>;
};


export type OmicsTranscriptionFactorOmicsTranscriptionFactorGenesByTranscriptionFactorIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsTranscriptionFactorGeneCondition>;
  filter: InputMaybe<OmicsTranscriptionFactorGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsTranscriptionFactorGenesOrderBy>>;
};


export type OmicsTranscriptionFactorOmicsTranscriptionFactorGenesByTranscriptionFactorIdListArgs = {
  condition: InputMaybe<OmicsTranscriptionFactorGeneCondition>;
  filter: InputMaybe<OmicsTranscriptionFactorGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsTranscriptionFactorGenesOrderBy>>;
};

/**
 * A condition to be used against `OmicsTranscriptionFactor` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type OmicsTranscriptionFactorCondition = {
  /** Checks for equality with the object’s `jasparMatrixId` field. */
  jasparMatrixId: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `OmicsTranscriptionFactor` object types. All fields are combined with a logical ‘and.’ */
export type OmicsTranscriptionFactorFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<OmicsTranscriptionFactorFilter>>;
  /** Filter by the object’s `jasparMatrixId` field. */
  jasparMatrixId: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<OmicsTranscriptionFactorFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<OmicsTranscriptionFactorFilter>>;
};

export type OmicsTranscriptionFactorGene = {
  __typename?: 'OmicsTranscriptionFactorGene';
  /** Reads a single `OmicsGene` that is related to this `OmicsTranscriptionFactorGene`. */
  gene: Maybe<OmicsGene>;
  geneId: Scalars['Int'];
  /** Reads a single `OmicsTranscriptionFactor` that is related to this `OmicsTranscriptionFactorGene`. */
  transcriptionFactor: Maybe<OmicsTranscriptionFactor>;
  transcriptionFactorId: Scalars['Int'];
};

/**
 * A condition to be used against `OmicsTranscriptionFactorGene` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type OmicsTranscriptionFactorGeneCondition = {
  /** Checks for equality with the object’s `geneId` field. */
  geneId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `transcriptionFactorId` field. */
  transcriptionFactorId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `OmicsTranscriptionFactorGene` object types. All fields are combined with a logical ‘and.’ */
export type OmicsTranscriptionFactorGeneFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<OmicsTranscriptionFactorGeneFilter>>;
  /** Filter by the object’s `geneId` field. */
  geneId: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not: InputMaybe<OmicsTranscriptionFactorGeneFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<OmicsTranscriptionFactorGeneFilter>>;
  /** Filter by the object’s `transcriptionFactorId` field. */
  transcriptionFactorId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `OmicsTranscriptionFactorGene` */
export type OmicsTranscriptionFactorGeneInput = {
  geneId: Scalars['Int'];
  transcriptionFactorId: Scalars['Int'];
};

/** A connection to a list of `OmicsTranscriptionFactorGene` values. */
export type OmicsTranscriptionFactorGenesConnection = {
  __typename?: 'OmicsTranscriptionFactorGenesConnection';
  /** A list of edges which contains the `OmicsTranscriptionFactorGene` and cursor to aid in pagination. */
  edges: Array<OmicsTranscriptionFactorGenesEdge>;
  /** A list of `OmicsTranscriptionFactorGene` objects. */
  nodes: Array<Maybe<OmicsTranscriptionFactorGene>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `OmicsTranscriptionFactorGene` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `OmicsTranscriptionFactorGene` edge in the connection. */
export type OmicsTranscriptionFactorGenesEdge = {
  __typename?: 'OmicsTranscriptionFactorGenesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `OmicsTranscriptionFactorGene` at the end of the edge. */
  node: Maybe<OmicsTranscriptionFactorGene>;
};

/** Methods to use when ordering `OmicsTranscriptionFactorGene`. */
export enum OmicsTranscriptionFactorGenesOrderBy {
  GeneIdAsc = 'GENE_ID_ASC',
  GeneIdDesc = 'GENE_ID_DESC',
  Natural = 'NATURAL',
  TranscriptionFactorIdAsc = 'TRANSCRIPTION_FACTOR_ID_ASC',
  TranscriptionFactorIdDesc = 'TRANSCRIPTION_FACTOR_ID_DESC'
}

/** An input for mutations affecting `OmicsTranscriptionFactor` */
export type OmicsTranscriptionFactorInput = {
  jasparMatrixId: Scalars['String'];
  omicsId: Scalars['Int'];
};

/** Represents an update to a `OmicsTranscriptionFactor`. Fields that are set will be updated. */
export type OmicsTranscriptionFactorPatch = {
  jasparMatrixId: InputMaybe<Scalars['String']>;
  omicsId: InputMaybe<Scalars['Int']>;
};

/** A connection to a list of `OmicsTranscriptionFactor` values. */
export type OmicsTranscriptionFactorsConnection = {
  __typename?: 'OmicsTranscriptionFactorsConnection';
  /** A list of edges which contains the `OmicsTranscriptionFactor` and cursor to aid in pagination. */
  edges: Array<OmicsTranscriptionFactorsEdge>;
  /** A list of `OmicsTranscriptionFactor` objects. */
  nodes: Array<Maybe<OmicsTranscriptionFactor>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `OmicsTranscriptionFactor` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `OmicsTranscriptionFactor` edge in the connection. */
export type OmicsTranscriptionFactorsEdge = {
  __typename?: 'OmicsTranscriptionFactorsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `OmicsTranscriptionFactor` at the end of the edge. */
  node: Maybe<OmicsTranscriptionFactor>;
};

/** Methods to use when ordering `OmicsTranscriptionFactor`. */
export enum OmicsTranscriptionFactorsOrderBy {
  JasparMatrixIdAsc = 'JASPAR_MATRIX_ID_ASC',
  JasparMatrixIdDesc = 'JASPAR_MATRIX_ID_DESC',
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC'
}

export enum OmicsType {
  Gene = 'GENE',
  ProteinAntibodyTag = 'PROTEIN_ANTIBODY_TAG',
  Region = 'REGION',
  TranscriptionFactor = 'TRANSCRIPTION_FACTOR'
}

/** A filter to be used against OmicsType fields. All fields are combined with a logical ‘and.’ */
export type OmicsTypeFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<OmicsType>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<OmicsType>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<OmicsType>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<OmicsType>;
  /** Included in the specified list. */
  in: InputMaybe<Array<OmicsType>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<OmicsType>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<OmicsType>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<OmicsType>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<OmicsType>;
  /** Not included in the specified list. */
  notIn: InputMaybe<Array<OmicsType>>;
};

/** A filter to be used against OmicsType List fields. All fields are combined with a logical ‘and.’ */
export type OmicsTypeListFilter = {
  /** Any array item is equal to the specified value. */
  anyEqualTo: InputMaybe<OmicsType>;
  /** Any array item is greater than the specified value. */
  anyGreaterThan: InputMaybe<OmicsType>;
  /** Any array item is greater than or equal to the specified value. */
  anyGreaterThanOrEqualTo: InputMaybe<OmicsType>;
  /** Any array item is less than the specified value. */
  anyLessThan: InputMaybe<OmicsType>;
  /** Any array item is less than or equal to the specified value. */
  anyLessThanOrEqualTo: InputMaybe<OmicsType>;
  /** Any array item is not equal to the specified value. */
  anyNotEqualTo: InputMaybe<OmicsType>;
  /** Contained by the specified list of values. */
  containedBy: InputMaybe<Array<InputMaybe<OmicsType>>>;
  /** Contains the specified list of values. */
  contains: InputMaybe<Array<InputMaybe<OmicsType>>>;
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<Array<InputMaybe<OmicsType>>>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<Array<InputMaybe<OmicsType>>>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<Array<InputMaybe<OmicsType>>>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<Array<InputMaybe<OmicsType>>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<Array<InputMaybe<OmicsType>>>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<Array<InputMaybe<OmicsType>>>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<Array<InputMaybe<OmicsType>>>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<Array<InputMaybe<OmicsType>>>;
  /** Overlaps the specified list of values. */
  overlaps: InputMaybe<Array<InputMaybe<OmicsType>>>;
};

/** A connection to a list of `OntCodesInfoRecord` values. */
export type OntCodesInfoConnection = {
  __typename?: 'OntCodesInfoConnection';
  /** A list of edges which contains the `OntCodesInfoRecord` and cursor to aid in pagination. */
  edges: Array<OntCodesInfoEdge>;
  /** A list of `OntCodesInfoRecord` objects. */
  nodes: Array<Maybe<OntCodesInfoRecord>>;
  /** The count of *all* `OntCodesInfoRecord` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `OntCodesInfoRecord` edge in the connection. */
export type OntCodesInfoEdge = {
  __typename?: 'OntCodesInfoEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `OntCodesInfoRecord` at the end of the edge. */
  node: Maybe<OntCodesInfoRecord>;
};

/** The return type of our `ontCodesInfo` query. */
export type OntCodesInfoRecord = {
  __typename?: 'OntCodesInfoRecord';
  labels: Maybe<Array<Maybe<Scalars['String']>>>;
  parentIds: Maybe<Array<Maybe<Scalars['String']>>>;
};

/** A filter to be used against `OntCodesInfoRecord` object types. All fields are combined with a logical ‘and.’ */
export type OntCodesInfoRecordFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<OntCodesInfoRecordFilter>>;
  /** Filter by the object’s `labels` field. */
  labels: InputMaybe<StringListFilter>;
  /** Negates the expression. */
  not: InputMaybe<OntCodesInfoRecordFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<OntCodesInfoRecordFilter>>;
  /** Filter by the object’s `parentIds` field. */
  parentIds: InputMaybe<StringListFilter>;
};

/** A connection to a list of `Ontology` values. */
export type OntologiesConnection = {
  __typename?: 'OntologiesConnection';
  /** A list of edges which contains the `Ontology` and cursor to aid in pagination. */
  edges: Array<OntologiesEdge>;
  /** A list of `Ontology` objects. */
  nodes: Array<Maybe<Ontology>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `Ontology` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `Ontology` edge in the connection. */
export type OntologiesEdge = {
  __typename?: 'OntologiesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `Ontology` at the end of the edge. */
  node: Maybe<Ontology>;
};

/** Methods to use when ordering `Ontology`. */
export enum OntologiesOrderBy {
  NameAsc = 'NAME_ASC',
  NameDesc = 'NAME_DESC',
  Natural = 'NATURAL',
  OntidAsc = 'ONTID_ASC',
  OntidDesc = 'ONTID_DESC',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC'
}

export type Ontology = Node & {
  __typename?: 'Ontology';
  /** Reads and enables pagination through a set of `Concept`. */
  conceptsByOntid: ConceptsConnection;
  /** Reads and enables pagination through a set of `Concept`. */
  conceptsByOntidList: Array<Concept>;
  name: Maybe<Scalars['String']>;
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  ontid: Scalars['Int'];
};


export type OntologyConceptsByOntidArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ConceptCondition>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptsOrderBy>>;
};


export type OntologyConceptsByOntidListArgs = {
  condition: InputMaybe<ConceptCondition>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ConceptsOrderBy>>;
};

/**
 * A condition to be used against `Ontology` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type OntologyCondition = {
  /** Checks for equality with the object’s `name` field. */
  name: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `ontid` field. */
  ontid: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `Ontology` object types. All fields are combined with a logical ‘and.’ */
export type OntologyFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<OntologyFilter>>;
  /** Filter by the object’s `name` field. */
  name: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<OntologyFilter>;
  /** Filter by the object’s `ontid` field. */
  ontid: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<OntologyFilter>>;
};

/** An input for mutations affecting `Ontology` */
export type OntologyInput = {
  name: InputMaybe<Scalars['String']>;
  ontid: Scalars['Int'];
};

/** Represents an update to a `Ontology`. Fields that are set will be updated. */
export type OntologyPatch = {
  name: InputMaybe<Scalars['String']>;
  ontid: InputMaybe<Scalars['Int']>;
};

/** Information about pagination in a connection. */
export type PageInfo = {
  __typename?: 'PageInfo';
  /** When paginating forwards, the cursor to continue. */
  endCursor: Maybe<Scalars['Cursor']>;
  /** When paginating forwards, are there more items? */
  hasNextPage: Scalars['Boolean'];
  /** When paginating backwards, are there more items? */
  hasPreviousPage: Scalars['Boolean'];
  /** When paginating backwards, the cursor to continue. */
  startCursor: Maybe<Scalars['Cursor']>;
};

/** The root query type which gives access points into the data universe. */
export type Query = Node & {
  __typename?: 'Query';
  /** Reads and enables pagination through a set of `_AllUsedOntologyId`. */
  _allUsedOntologyIds: Maybe<_AllUsedOntologyIdsConnection>;
  /** Reads a set of `_AllUsedOntologyId`. */
  _allUsedOntologyIdsList: Maybe<Array<_AllUsedOntologyId>>;
  _conceptHierarchyMinimumTreesImpl: Maybe<MinimumTreesResult>;
  _finalBoxplot: Maybe<BoxplotValue>;
  _semanticOrderImpl: Maybe<Array<Maybe<Scalars['Int']>>>;
  annotationGroup: Maybe<AnnotationGroup>;
  /** Reads a single `AnnotationGroup` using its globally unique `ID`. */
  annotationGroupByNodeId: Maybe<AnnotationGroup>;
  /** Reads and enables pagination through a set of `AnnotationGroup`. */
  annotationGroups: Maybe<AnnotationGroupsConnection>;
  /** Reads a set of `AnnotationGroup`. */
  annotationGroupsList: Maybe<Array<AnnotationGroup>>;
  annotationValue: Maybe<AnnotationValue>;
  /** Reads a single `AnnotationValue` using its globally unique `ID`. */
  annotationValueByNodeId: Maybe<AnnotationValue>;
  annotationValueCoocurrence: Maybe<AnnotationValueCoocurrenceConnection>;
  annotationValueCoocurrenceList: Maybe<Array<Maybe<AnnotationValueCoocurrenceRecord>>>;
  /** Reads and enables pagination through a set of `AnnotationValue`. */
  annotationValues: Maybe<AnnotationValuesConnection>;
  /** Reads a set of `AnnotationValue`. */
  annotationValuesList: Maybe<Array<AnnotationValue>>;
  /** Reads and enables pagination through a set of `ApiDifferentialExpression`. */
  apiDifferentialExpressions: Maybe<ApiDifferentialExpressionsConnection>;
  /** Reads a set of `ApiDifferentialExpression`. */
  apiDifferentialExpressionsList: Maybe<Array<ApiDifferentialExpression>>;
  /** Reads and enables pagination through a set of `ApiExpressionByAnnotation`. */
  apiExpressionByAnnotation: Maybe<ApiExpressionByAnnotationsConnection>;
  /** Reads and enables pagination through a set of `ApiExpressionByAnnotation`. */
  apiExpressionByAnnotationList: Maybe<Array<Maybe<ApiExpressionByAnnotation>>>;
  /** Reads and enables pagination through a set of `ApiOmic`. */
  apiOmics: Maybe<ApiOmicsConnection>;
  /** Reads a set of `ApiOmic`. */
  apiOmicsList: Maybe<Array<ApiOmic>>;
  /** Reads and enables pagination through a set of `ApiStudiesBase`. */
  apiStudiesBases: Maybe<ApiStudiesBasesConnection>;
  /** Reads a set of `ApiStudiesBase`. */
  apiStudiesBasesList: Maybe<Array<ApiStudiesBase>>;
  /** Reads and enables pagination through a set of `ApiStudiesBulkRna`. */
  apiStudiesBulkRnas: Maybe<ApiStudiesBulkRnasConnection>;
  /** Reads a set of `ApiStudiesBulkRna`. */
  apiStudiesBulkRnasList: Maybe<Array<ApiStudiesBulkRna>>;
  /** Reads and enables pagination through a set of `ApiStudiesSingleCell`. */
  apiStudiesSingleCells: Maybe<ApiStudiesSingleCellsConnection>;
  /** Reads a set of `ApiStudiesSingleCell`. */
  apiStudiesSingleCellsList: Maybe<Array<ApiStudiesSingleCell>>;
  /** Reads and enables pagination through a set of `ApiStudyAnnotationOverview`. */
  apiStudyAnnotationOverviews: Maybe<ApiStudyAnnotationOverviewsConnection>;
  /** Reads a set of `ApiStudyAnnotationOverview`. */
  apiStudyAnnotationOverviewsList: Maybe<Array<ApiStudyAnnotationOverview>>;
  /** Reads and enables pagination through a set of `ApiStudyH5Download`. */
  apiStudyH5Downloads: Maybe<ApiStudyH5DownloadsConnection>;
  /** Reads a set of `ApiStudyH5Download`. */
  apiStudyH5DownloadsList: Maybe<Array<ApiStudyH5Download>>;
  /** Reads and enables pagination through a set of `AutocompleteResult`. */
  autocomplete: Maybe<AutocompleteResultsConnection>;
  /** Reads and enables pagination through a set of `AutocompleteResult`. */
  autocompleteList: Maybe<Array<Maybe<AutocompleteResult>>>;
  concept: Maybe<Concept>;
  /** Reads a single `Concept` using its globally unique `ID`. */
  conceptByNodeId: Maybe<Concept>;
  conceptCidArrayToCodes: Maybe<Array<Maybe<Scalars['String']>>>;
  /** Reads and enables pagination through a set of `ConceptHierarchy`. */
  conceptHierarchies: Maybe<ConceptHierarchiesConnection>;
  /** Reads a set of `ConceptHierarchy`. */
  conceptHierarchiesList: Maybe<Array<ConceptHierarchy>>;
  conceptHierarchyMinimumTreesParentsLists: Maybe<ConceptHierarchyMinimumTreesParentsListsConnection>;
  conceptHierarchyMinimumTreesParentsListsList: Maybe<Array<Maybe<ConceptHierarchyMinimumTreesParentsListsRecord>>>;
  /** Reads and enables pagination through a set of `ConceptSynonym`. */
  conceptSynonyms: Maybe<ConceptSynonymsConnection>;
  /** Reads a set of `ConceptSynonym`. */
  conceptSynonymsList: Maybe<Array<ConceptSynonym>>;
  /** Reads and enables pagination through a set of `Concept`. */
  concepts: Maybe<ConceptsConnection>;
  /** Reads and enables pagination through a set of `Concept`. */
  conceptsInSemanticOrder: Maybe<ConceptsConnection>;
  /** Reads and enables pagination through a set of `Concept`. */
  conceptsInSemanticOrderList: Maybe<Array<Maybe<Concept>>>;
  /** Reads a set of `Concept`. */
  conceptsList: Maybe<Array<Concept>>;
  correlationTrianglePlot: Maybe<Scalars['String']>;
  currentUserEmail: Maybe<Scalars['String']>;
  currentUserGroups: Maybe<Array<Maybe<Scalars['String']>>>;
  /** Reads and enables pagination through a set of `DifferentialExpressionV`. */
  differentialExpressionVs: Maybe<DifferentialExpressionVsConnection>;
  /** Reads a set of `DifferentialExpressionV`. */
  differentialExpressionVsList: Maybe<Array<DifferentialExpressionV>>;
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressions: Maybe<DifferentialExpressionsConnection>;
  /** Reads a set of `DifferentialExpression`. */
  differentialExpressionsList: Maybe<Array<DifferentialExpression>>;
  ensureTextArray: Maybe<Array<Maybe<Scalars['String']>>>;
  /** Reads and enables pagination through a set of `ExpressionByAnnotation`. */
  expressionByAnnotation: Maybe<ExpressionByAnnotationsConnection>;
  /** Reads and enables pagination through a set of `ExpressionByAnnotation`. */
  expressionByAnnotationList: Maybe<Array<Maybe<ExpressionByAnnotation>>>;
  /** Reads and enables pagination through a set of `ExpressionByOmic`. */
  expressionByOmicsIds: Maybe<ExpressionByOmicsConnection>;
  /** Reads and enables pagination through a set of `ExpressionByOmic`. */
  expressionByOmicsIdsList: Maybe<Array<Maybe<ExpressionByOmic>>>;
  /** Reads and enables pagination through a set of `ExpressionByTwoAnnotation`. */
  expressionByTwoAnnotations: Maybe<ExpressionByTwoAnnotationsConnection>;
  /** Reads and enables pagination through a set of `ExpressionByTwoAnnotation`. */
  expressionByTwoAnnotationsList: Maybe<Array<Maybe<ExpressionByTwoAnnotation>>>;
  expressionTtest: Maybe<Scalars['String']>;
  getCorrelatedGenes: Maybe<GetCorrelatedGenesConnection>;
  getCorrelatedGenesList: Maybe<Array<Maybe<GetCorrelatedGenesRecord>>>;
  /** Fetches an object given its globally unique `ID`. */
  node: Maybe<Node>;
  /** The root query type must be a `Node` to work well with Relay 1 mutations. This just resolves to `query`. */
  nodeId: Scalars['ID'];
  /** Reads and enables pagination through a set of `OmicsAll`. */
  omicsAlls: Maybe<OmicsAllsConnection>;
  /** Reads a set of `OmicsAll`. */
  omicsAllsList: Maybe<Array<OmicsAll>>;
  /** Reads and enables pagination through a set of `OmicsAutocompleteResult`. */
  omicsAutocomplete: Maybe<OmicsAutocompleteResultsConnection>;
  /** Reads and enables pagination through a set of `OmicsAutocompleteResult`. */
  omicsAutocompleteList: Maybe<Array<Maybe<OmicsAutocompleteResult>>>;
  omicsBase: Maybe<OmicsBase>;
  /** Reads a single `OmicsBase` using its globally unique `ID`. */
  omicsBaseByNodeId: Maybe<OmicsBase>;
  /** Reads and enables pagination through a set of `OmicsBase`. */
  omicsBases: Maybe<OmicsBasesConnection>;
  /** Reads a set of `OmicsBase`. */
  omicsBasesList: Maybe<Array<OmicsBase>>;
  omicsGene: Maybe<OmicsGene>;
  /** Reads a single `OmicsGene` using its globally unique `ID`. */
  omicsGeneByNodeId: Maybe<OmicsGene>;
  /** Reads and enables pagination through a set of `OmicsGene`. */
  omicsGenes: Maybe<OmicsGenesConnection>;
  /** Reads a set of `OmicsGene`. */
  omicsGenesList: Maybe<Array<OmicsGene>>;
  omicsProteinAntibodyTag: Maybe<OmicsProteinAntibodyTag>;
  /** Reads a single `OmicsProteinAntibodyTag` using its globally unique `ID`. */
  omicsProteinAntibodyTagByNodeId: Maybe<OmicsProteinAntibodyTag>;
  /** Reads and enables pagination through a set of `OmicsProteinAntibodyTagGene`. */
  omicsProteinAntibodyTagGenes: Maybe<OmicsProteinAntibodyTagGenesConnection>;
  /** Reads a set of `OmicsProteinAntibodyTagGene`. */
  omicsProteinAntibodyTagGenesList: Maybe<Array<OmicsProteinAntibodyTagGene>>;
  /** Reads and enables pagination through a set of `OmicsProteinAntibodyTag`. */
  omicsProteinAntibodyTags: Maybe<OmicsProteinAntibodyTagsConnection>;
  /** Reads a set of `OmicsProteinAntibodyTag`. */
  omicsProteinAntibodyTagsList: Maybe<Array<OmicsProteinAntibodyTag>>;
  omicsRegion: Maybe<OmicsRegion>;
  /** Reads a single `OmicsRegion` using its globally unique `ID`. */
  omicsRegionByNodeId: Maybe<OmicsRegion>;
  /** Reads and enables pagination through a set of `OmicsRegionGene`. */
  omicsRegionGenes: Maybe<OmicsRegionGenesConnection>;
  /** Reads a set of `OmicsRegionGene`. */
  omicsRegionGenesList: Maybe<Array<OmicsRegionGene>>;
  /** Reads and enables pagination through a set of `OmicsRegion`. */
  omicsRegions: Maybe<OmicsRegionsConnection>;
  /** Reads a set of `OmicsRegion`. */
  omicsRegionsList: Maybe<Array<OmicsRegion>>;
  omicsTranscriptionFactor: Maybe<OmicsTranscriptionFactor>;
  /** Reads a single `OmicsTranscriptionFactor` using its globally unique `ID`. */
  omicsTranscriptionFactorByNodeId: Maybe<OmicsTranscriptionFactor>;
  /** Reads and enables pagination through a set of `OmicsTranscriptionFactorGene`. */
  omicsTranscriptionFactorGenes: Maybe<OmicsTranscriptionFactorGenesConnection>;
  /** Reads a set of `OmicsTranscriptionFactorGene`. */
  omicsTranscriptionFactorGenesList: Maybe<Array<OmicsTranscriptionFactorGene>>;
  /** Reads and enables pagination through a set of `OmicsTranscriptionFactor`. */
  omicsTranscriptionFactors: Maybe<OmicsTranscriptionFactorsConnection>;
  /** Reads a set of `OmicsTranscriptionFactor`. */
  omicsTranscriptionFactorsList: Maybe<Array<OmicsTranscriptionFactor>>;
  ontCodesInfo: Maybe<OntCodesInfoConnection>;
  ontCodesInfoList: Maybe<Array<Maybe<OntCodesInfoRecord>>>;
  /** Reads and enables pagination through a set of `Ontology`. */
  ontologies: Maybe<OntologiesConnection>;
  /** Reads a set of `Ontology`. */
  ontologiesList: Maybe<Array<Ontology>>;
  ontology: Maybe<Ontology>;
  /** Reads a single `Ontology` using its globally unique `ID`. */
  ontologyByNodeId: Maybe<Ontology>;
  /**
   * Exposes the root query type nested one level down. This is helpful for Relay 1
   * which can only query top level fields if they are in a particular form.
   */
  query: Query;
  /** Reads and enables pagination through a set of `ReferenceStudy`. */
  referenceStudies: Maybe<ReferenceStudiesConnection>;
  /** Reads a set of `ReferenceStudy`. */
  referenceStudiesList: Maybe<Array<ReferenceStudy>>;
  /** Reads and enables pagination through a set of `ReferenceStudyOverview`. */
  referenceStudyOverviews: Maybe<ReferenceStudyOverviewsConnection>;
  /** Reads a set of `ReferenceStudyOverview`. */
  referenceStudyOverviewsList: Maybe<Array<ReferenceStudyOverview>>;
  /** Reads and enables pagination through a set of `Study`. */
  studies: Maybe<StudiesConnection>;
  /** Reads a set of `Study`. */
  studiesList: Maybe<Array<Study>>;
  study: Maybe<Study>;
  /** Reads and enables pagination through a set of `StudyAdminDetail`. */
  studyAdminDetails: Maybe<StudyAdminDetailsConnection>;
  /** Reads a set of `StudyAdminDetail`. */
  studyAdminDetailsList: Maybe<Array<StudyAdminDetail>>;
  /** Reads and enables pagination through a set of `StudyAdministrableCurrentuser`. */
  studyAdministrableCurrentusers: Maybe<StudyAdministrableCurrentusersConnection>;
  /** Reads a set of `StudyAdministrableCurrentuser`. */
  studyAdministrableCurrentusersList: Maybe<Array<StudyAdministrableCurrentuser>>;
  /** Reads and enables pagination through a set of `StudyAnnotationFrontendGroup`. */
  studyAnnotationFrontendGroups: Maybe<StudyAnnotationFrontendGroupsConnection>;
  /** Reads a set of `StudyAnnotationFrontendGroup`. */
  studyAnnotationFrontendGroupsList: Maybe<Array<StudyAnnotationFrontendGroup>>;
  /** Reads and enables pagination through a set of `StudyAnnotationFrontendValue`. */
  studyAnnotationFrontendValues: Maybe<StudyAnnotationFrontendValuesConnection>;
  /** Reads a set of `StudyAnnotationFrontendValue`. */
  studyAnnotationFrontendValuesList: Maybe<Array<StudyAnnotationFrontendValue>>;
  /** Reads and enables pagination through a set of `StudyAnnotationGroupUi`. */
  studyAnnotationGroupUis: Maybe<StudyAnnotationGroupUisConnection>;
  /** Reads a set of `StudyAnnotationGroupUi`. */
  studyAnnotationGroupUisList: Maybe<Array<StudyAnnotationGroupUi>>;
  /** Reads a single `Study` using its globally unique `ID`. */
  studyByNodeId: Maybe<Study>;
  /** Reads and enables pagination through a set of `StudyImportLog`. */
  studyImportLogs: Maybe<StudyImportLogsConnection>;
  /** Reads a set of `StudyImportLog`. */
  studyImportLogsList: Maybe<Array<StudyImportLog>>;
  studyLayer: Maybe<StudyLayer>;
  /** Reads a single `StudyLayer` using its globally unique `ID`. */
  studyLayerByNodeId: Maybe<StudyLayer>;
  /** Reads and enables pagination through a set of `StudyLayer`. */
  studyLayers: Maybe<StudyLayersConnection>;
  /** Reads a set of `StudyLayer`. */
  studyLayersList: Maybe<Array<StudyLayer>>;
  /** Reads and enables pagination through a set of `StudyOmic`. */
  studyOmics: Maybe<StudyOmicsConnection>;
  /** Reads a set of `StudyOmic`. */
  studyOmicsList: Maybe<Array<StudyOmic>>;
  /** Reads and enables pagination through a set of `StudyOmicsTransposed`. */
  studyOmicsTransposeds: Maybe<StudyOmicsTransposedsConnection>;
  /** Reads a set of `StudyOmicsTransposed`. */
  studyOmicsTransposedsList: Maybe<Array<StudyOmicsTransposed>>;
  /** Reads and enables pagination through a set of `StudyOverviewOntology`. */
  studyOverviewOntologies: Maybe<StudyOverviewOntologiesConnection>;
  /** Reads a set of `StudyOverviewOntology`. */
  studyOverviewOntologiesList: Maybe<Array<StudyOverviewOntology>>;
  /** Reads and enables pagination through a set of `StudyOverview`. */
  studyOverviews: Maybe<StudyOverviewsConnection>;
  /** Reads a set of `StudyOverview`. */
  studyOverviewsList: Maybe<Array<StudyOverview>>;
  studySample: Maybe<StudySample>;
  /** Reads and enables pagination through a set of `StudySampleAnnotationSubsampling`. */
  studySampleAnnotationSubsamplings: Maybe<StudySampleAnnotationSubsamplingsConnection>;
  /** Reads a set of `StudySampleAnnotationSubsampling`. */
  studySampleAnnotationSubsamplingsList: Maybe<Array<StudySampleAnnotationSubsampling>>;
  /** Reads and enables pagination through a set of `StudySampleAnnotation`. */
  studySampleAnnotations: Maybe<StudySampleAnnotationsConnection>;
  /** Reads a set of `StudySampleAnnotation`. */
  studySampleAnnotationsList: Maybe<Array<StudySampleAnnotation>>;
  /** Reads a single `StudySample` using its globally unique `ID`. */
  studySampleByNodeId: Maybe<StudySample>;
  /** Reads and enables pagination through a set of `StudySampleProjectionSubsamplingTransposed`. */
  studySampleProjectionSubsamplingTransposeds: Maybe<StudySampleProjectionSubsamplingTransposedsConnection>;
  /** Reads a set of `StudySampleProjectionSubsamplingTransposed`. */
  studySampleProjectionSubsamplingTransposedsList: Maybe<Array<StudySampleProjectionSubsamplingTransposed>>;
  /** Reads and enables pagination through a set of `StudySampleProjection`. */
  studySampleProjections: Maybe<StudySampleProjectionsConnection>;
  /** Reads a set of `StudySampleProjection`. */
  studySampleProjectionsList: Maybe<Array<StudySampleProjection>>;
  /** Reads and enables pagination through a set of `StudySample`. */
  studySamples: Maybe<StudySamplesConnection>;
  /** Reads a set of `StudySample`. */
  studySamplesList: Maybe<Array<StudySample>>;
  /** Reads and enables pagination through a set of `StudyVisibleCurrentuser`. */
  studyVisibleCurrentusers: Maybe<StudyVisibleCurrentusersConnection>;
  /** Reads a set of `StudyVisibleCurrentuser`. */
  studyVisibleCurrentusersList: Maybe<Array<StudyVisibleCurrentuser>>;
  /** Reads and enables pagination through a set of `TreeOntology`. */
  treeOntologies: Maybe<TreeOntologiesConnection>;
  /** Reads a set of `TreeOntology`. */
  treeOntologiesList: Maybe<Array<TreeOntology>>;
  /** Reads and enables pagination through a set of `UserAnnotationGroup`. */
  userAnnotationGroups: Maybe<UserAnnotationGroupsConnection>;
  /** Reads a set of `UserAnnotationGroup`. */
  userAnnotationGroupsList: Maybe<Array<UserAnnotationGroup>>;
  userStudyUploadConfigured: Maybe<Scalars['Boolean']>;
  violinPlot: Maybe<Scalars['String']>;
};


/** The root query type which gives access points into the data universe. */
export type Query_AllUsedOntologyIdsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<_AllUsedOntologyIdCondition>;
  filter: InputMaybe<_AllUsedOntologyIdFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<_AllUsedOntologyIdsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type Query_AllUsedOntologyIdsListArgs = {
  condition: InputMaybe<_AllUsedOntologyIdCondition>;
  filter: InputMaybe<_AllUsedOntologyIdFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<_AllUsedOntologyIdsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type Query_ConceptHierarchyMinimumTreesImplArgs = {
  cidsLeaves: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  fullgraph: InputMaybe<Array<InputMaybe<ConceptWeightedParentInput>>>;
};


/** The root query type which gives access points into the data universe. */
export type Query_FinalBoxplotArgs = {
  a: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};


/** The root query type which gives access points into the data universe. */
export type Query_SemanticOrderImplArgs = {
  cidsToOrder: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  fullgraph: InputMaybe<Array<InputMaybe<ConceptWeightedParentInput>>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryAnnotationGroupArgs = {
  annotationGroupId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryAnnotationGroupByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryAnnotationGroupsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<AnnotationGroupCondition>;
  filter: InputMaybe<AnnotationGroupFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<AnnotationGroupsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryAnnotationGroupsListArgs = {
  condition: InputMaybe<AnnotationGroupCondition>;
  filter: InputMaybe<AnnotationGroupFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<AnnotationGroupsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryAnnotationValueArgs = {
  annotationValueId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryAnnotationValueByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryAnnotationValueCoocurrenceArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  annotationGroupId1: InputMaybe<Scalars['Int']>;
  annotationGroupId2: InputMaybe<Scalars['Int']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<AnnotationValueCoocurrenceRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  studyId: InputMaybe<Scalars['Int']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryAnnotationValueCoocurrenceListArgs = {
  annotationGroupId1: InputMaybe<Scalars['Int']>;
  annotationGroupId2: InputMaybe<Scalars['Int']>;
  filter: InputMaybe<AnnotationValueCoocurrenceRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  studyId: InputMaybe<Scalars['Int']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryAnnotationValuesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<AnnotationValueCondition>;
  filter: InputMaybe<AnnotationValueFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<AnnotationValuesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryAnnotationValuesListArgs = {
  condition: InputMaybe<AnnotationValueCondition>;
  filter: InputMaybe<AnnotationValueFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<AnnotationValuesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiDifferentialExpressionsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ApiDifferentialExpressionCondition>;
  filter: InputMaybe<ApiDifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ApiDifferentialExpressionsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiDifferentialExpressionsListArgs = {
  condition: InputMaybe<ApiDifferentialExpressionCondition>;
  filter: InputMaybe<ApiDifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ApiDifferentialExpressionsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiExpressionByAnnotationArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<ApiExpressionByAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  pAnnotationGroupColumn: InputMaybe<Scalars['String']>;
  pEnsemblGeneIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  pEntrezGeneIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  pHgncSymbols: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  pLayerName: InputMaybe<Scalars['String']>;
  pOmicsIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pStudyIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiExpressionByAnnotationListArgs = {
  filter: InputMaybe<ApiExpressionByAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  pAnnotationGroupColumn: InputMaybe<Scalars['String']>;
  pEnsemblGeneIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  pEntrezGeneIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  pHgncSymbols: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  pLayerName: InputMaybe<Scalars['String']>;
  pOmicsIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pStudyIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiOmicsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ApiOmicCondition>;
  filter: InputMaybe<ApiOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ApiOmicsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiOmicsListArgs = {
  condition: InputMaybe<ApiOmicCondition>;
  filter: InputMaybe<ApiOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ApiOmicsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiStudiesBasesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ApiStudiesBaseCondition>;
  filter: InputMaybe<ApiStudiesBaseFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ApiStudiesBasesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiStudiesBasesListArgs = {
  condition: InputMaybe<ApiStudiesBaseCondition>;
  filter: InputMaybe<ApiStudiesBaseFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ApiStudiesBasesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiStudiesBulkRnasArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ApiStudiesBulkRnaCondition>;
  filter: InputMaybe<ApiStudiesBulkRnaFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ApiStudiesBulkRnasOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiStudiesBulkRnasListArgs = {
  condition: InputMaybe<ApiStudiesBulkRnaCondition>;
  filter: InputMaybe<ApiStudiesBulkRnaFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ApiStudiesBulkRnasOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiStudiesSingleCellsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ApiStudiesSingleCellCondition>;
  filter: InputMaybe<ApiStudiesSingleCellFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ApiStudiesSingleCellsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiStudiesSingleCellsListArgs = {
  condition: InputMaybe<ApiStudiesSingleCellCondition>;
  filter: InputMaybe<ApiStudiesSingleCellFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ApiStudiesSingleCellsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiStudyAnnotationOverviewsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ApiStudyAnnotationOverviewCondition>;
  filter: InputMaybe<ApiStudyAnnotationOverviewFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ApiStudyAnnotationOverviewsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiStudyAnnotationOverviewsListArgs = {
  condition: InputMaybe<ApiStudyAnnotationOverviewCondition>;
  filter: InputMaybe<ApiStudyAnnotationOverviewFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ApiStudyAnnotationOverviewsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiStudyH5DownloadsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ApiStudyH5DownloadCondition>;
  filter: InputMaybe<ApiStudyH5DownloadFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ApiStudyH5DownloadsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryApiStudyH5DownloadsListArgs = {
  condition: InputMaybe<ApiStudyH5DownloadCondition>;
  filter: InputMaybe<ApiStudyH5DownloadFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ApiStudyH5DownloadsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryAutocompleteArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<AutocompleteResultFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  searchQuery: InputMaybe<Scalars['String']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryAutocompleteListArgs = {
  filter: InputMaybe<AutocompleteResultFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  searchQuery: InputMaybe<Scalars['String']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptArgs = {
  cid: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptCidArrayToCodesArgs = {
  arg0: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptHierarchiesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ConceptHierarchyCondition>;
  filter: InputMaybe<ConceptHierarchyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptHierarchiesListArgs = {
  condition: InputMaybe<ConceptHierarchyCondition>;
  filter: InputMaybe<ConceptHierarchyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptHierarchyMinimumTreesParentsListsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<ConceptHierarchyMinimumTreesParentsListsRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  ontologyCodes: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  queryOntology: InputMaybe<Scalars['String']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptHierarchyMinimumTreesParentsListsListArgs = {
  filter: InputMaybe<ConceptHierarchyMinimumTreesParentsListsRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  ontologyCodes: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  queryOntology: InputMaybe<Scalars['String']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptSynonymsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ConceptSynonymCondition>;
  filter: InputMaybe<ConceptSynonymFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptSynonymsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptSynonymsListArgs = {
  condition: InputMaybe<ConceptSynonymCondition>;
  filter: InputMaybe<ConceptSynonymFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ConceptSynonymsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ConceptCondition>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptsInSemanticOrderArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  ontologyCodes: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  queryOntology: InputMaybe<Scalars['String']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptsInSemanticOrderListArgs = {
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  ontologyCodes: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  queryOntology: InputMaybe<Scalars['String']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptsListArgs = {
  condition: InputMaybe<ConceptCondition>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ConceptsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryCorrelationTrianglePlotArgs = {
  pExcludeAnnotationValueIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pOmicsIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pStudyId: InputMaybe<Scalars['Int']>;
  pStudyLayerId: InputMaybe<Scalars['Int']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryDifferentialExpressionVsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<DifferentialExpressionVCondition>;
  filter: InputMaybe<DifferentialExpressionVFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<DifferentialExpressionVsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryDifferentialExpressionVsListArgs = {
  condition: InputMaybe<DifferentialExpressionVCondition>;
  filter: InputMaybe<DifferentialExpressionVFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<DifferentialExpressionVsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryDifferentialExpressionsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<DifferentialExpressionCondition>;
  filter: InputMaybe<DifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryDifferentialExpressionsListArgs = {
  condition: InputMaybe<DifferentialExpressionCondition>;
  filter: InputMaybe<DifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryEnsureTextArrayArgs = {
  a: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpressionByAnnotationArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<ExpressionByAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  pAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pExcludeAnnotationValueIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pOmicsIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pStudyLayerIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpressionByAnnotationListArgs = {
  filter: InputMaybe<ExpressionByAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  pAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pExcludeAnnotationValueIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pOmicsIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pStudyLayerIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpressionByOmicsIdsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<ExpressionByOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  pOmicsIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pStudyLayerId: InputMaybe<Scalars['Int']>;
  pSubsamplingProjection: InputMaybe<Scalars['String']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpressionByOmicsIdsListArgs = {
  filter: InputMaybe<ExpressionByOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  pOmicsIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pStudyLayerId: InputMaybe<Scalars['Int']>;
  pSubsamplingProjection: InputMaybe<Scalars['String']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpressionByTwoAnnotationsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<ExpressionByTwoAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  pAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pDropoutsAsZero: InputMaybe<Scalars['Boolean']>;
  pExcludeAnnotationValueIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pOmicsIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pSecondAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pStudyId: InputMaybe<Scalars['Int']>;
  pStudyLayerId: InputMaybe<Scalars['Int']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpressionByTwoAnnotationsListArgs = {
  filter: InputMaybe<ExpressionByTwoAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  pAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pDropoutsAsZero: InputMaybe<Scalars['Boolean']>;
  pExcludeAnnotationValueIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pOmicsIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pSecondAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pStudyId: InputMaybe<Scalars['Int']>;
  pStudyLayerId: InputMaybe<Scalars['Int']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpressionTtestArgs = {
  pAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pExcludeAnnotationValueIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pOmicsId: InputMaybe<Scalars['Int']>;
  pSample1AnnotationValueId: InputMaybe<Scalars['Int']>;
  pSample1SecondAnnotationValueId: InputMaybe<Scalars['Int']>;
  pSample2AnnotationValueId: InputMaybe<Scalars['Int']>;
  pSample2SecondAnnotationValueId: InputMaybe<Scalars['Int']>;
  pSecondaryAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pStudyId: InputMaybe<Scalars['Int']>;
  pStudyLayerId: InputMaybe<Scalars['Int']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryGetCorrelatedGenesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<GetCorrelatedGenesRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  omicsId: InputMaybe<Scalars['Int']>;
  studyId: InputMaybe<Scalars['Int']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryGetCorrelatedGenesListArgs = {
  filter: InputMaybe<GetCorrelatedGenesRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  omicsId: InputMaybe<Scalars['Int']>;
  studyId: InputMaybe<Scalars['Int']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryNodeArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsAllsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsAllCondition>;
  filter: InputMaybe<OmicsAllFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsAllsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsAllsListArgs = {
  condition: InputMaybe<OmicsAllCondition>;
  filter: InputMaybe<OmicsAllFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsAllsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsAutocompleteArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<OmicsAutocompleteResultFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  omicsTypeFilter: InputMaybe<OmicsType>;
  searchQuery: InputMaybe<Scalars['String']>;
  taxIdFilter: InputMaybe<Scalars['Int']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsAutocompleteListArgs = {
  filter: InputMaybe<OmicsAutocompleteResultFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  omicsTypeFilter: InputMaybe<OmicsType>;
  searchQuery: InputMaybe<Scalars['String']>;
  taxIdFilter: InputMaybe<Scalars['Int']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsBaseArgs = {
  omicsId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsBaseByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsBasesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsBaseCondition>;
  filter: InputMaybe<OmicsBaseFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsBasesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsBasesListArgs = {
  condition: InputMaybe<OmicsBaseCondition>;
  filter: InputMaybe<OmicsBaseFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsBasesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsGeneArgs = {
  geneId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsGeneByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsGenesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsGeneCondition>;
  filter: InputMaybe<OmicsGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsGenesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsGenesListArgs = {
  condition: InputMaybe<OmicsGeneCondition>;
  filter: InputMaybe<OmicsGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsGenesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsProteinAntibodyTagArgs = {
  proteinAntibodyTagId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsProteinAntibodyTagByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsProteinAntibodyTagGenesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsProteinAntibodyTagGeneCondition>;
  filter: InputMaybe<OmicsProteinAntibodyTagGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsProteinAntibodyTagGenesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsProteinAntibodyTagGenesListArgs = {
  condition: InputMaybe<OmicsProteinAntibodyTagGeneCondition>;
  filter: InputMaybe<OmicsProteinAntibodyTagGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsProteinAntibodyTagGenesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsProteinAntibodyTagsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsProteinAntibodyTagCondition>;
  filter: InputMaybe<OmicsProteinAntibodyTagFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsProteinAntibodyTagsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsProteinAntibodyTagsListArgs = {
  condition: InputMaybe<OmicsProteinAntibodyTagCondition>;
  filter: InputMaybe<OmicsProteinAntibodyTagFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsProteinAntibodyTagsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsRegionArgs = {
  regionId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsRegionByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsRegionGenesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsRegionGeneCondition>;
  filter: InputMaybe<OmicsRegionGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsRegionGenesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsRegionGenesListArgs = {
  condition: InputMaybe<OmicsRegionGeneCondition>;
  filter: InputMaybe<OmicsRegionGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsRegionGenesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsRegionsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsRegionCondition>;
  filter: InputMaybe<OmicsRegionFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsRegionsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsRegionsListArgs = {
  condition: InputMaybe<OmicsRegionCondition>;
  filter: InputMaybe<OmicsRegionFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsRegionsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsTranscriptionFactorArgs = {
  omicsId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsTranscriptionFactorByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsTranscriptionFactorGenesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsTranscriptionFactorGeneCondition>;
  filter: InputMaybe<OmicsTranscriptionFactorGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsTranscriptionFactorGenesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsTranscriptionFactorGenesListArgs = {
  condition: InputMaybe<OmicsTranscriptionFactorGeneCondition>;
  filter: InputMaybe<OmicsTranscriptionFactorGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsTranscriptionFactorGenesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsTranscriptionFactorsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OmicsTranscriptionFactorCondition>;
  filter: InputMaybe<OmicsTranscriptionFactorFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsTranscriptionFactorsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsTranscriptionFactorsListArgs = {
  condition: InputMaybe<OmicsTranscriptionFactorCondition>;
  filter: InputMaybe<OmicsTranscriptionFactorFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsTranscriptionFactorsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOntCodesInfoArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  filter: InputMaybe<OntCodesInfoRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  pOntCodes: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  pOntology: InputMaybe<Scalars['String']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOntCodesInfoListArgs = {
  filter: InputMaybe<OntCodesInfoRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  pOntCodes: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  pOntology: InputMaybe<Scalars['String']>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOntologiesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<OntologyCondition>;
  filter: InputMaybe<OntologyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OntologiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOntologiesListArgs = {
  condition: InputMaybe<OntologyCondition>;
  filter: InputMaybe<OntologyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OntologiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOntologyArgs = {
  ontid: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOntologyByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryReferenceStudiesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ReferenceStudyCondition>;
  filter: InputMaybe<ReferenceStudyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ReferenceStudiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryReferenceStudiesListArgs = {
  condition: InputMaybe<ReferenceStudyCondition>;
  filter: InputMaybe<ReferenceStudyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ReferenceStudiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryReferenceStudyOverviewsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ReferenceStudyOverviewCondition>;
  filter: InputMaybe<ReferenceStudyOverviewFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ReferenceStudyOverviewsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryReferenceStudyOverviewsListArgs = {
  condition: InputMaybe<ReferenceStudyOverviewCondition>;
  filter: InputMaybe<ReferenceStudyOverviewFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ReferenceStudyOverviewsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudiesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyCondition>;
  filter: InputMaybe<StudyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudiesListArgs = {
  condition: InputMaybe<StudyCondition>;
  filter: InputMaybe<StudyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyArgs = {
  studyId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyAdminDetailsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyAdminDetailCondition>;
  filter: InputMaybe<StudyAdminDetailFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyAdminDetailsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyAdminDetailsListArgs = {
  condition: InputMaybe<StudyAdminDetailCondition>;
  filter: InputMaybe<StudyAdminDetailFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyAdminDetailsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyAdministrableCurrentusersArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyAdministrableCurrentuserCondition>;
  filter: InputMaybe<StudyAdministrableCurrentuserFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyAdministrableCurrentusersOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyAdministrableCurrentusersListArgs = {
  condition: InputMaybe<StudyAdministrableCurrentuserCondition>;
  filter: InputMaybe<StudyAdministrableCurrentuserFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyAdministrableCurrentusersOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyAnnotationFrontendGroupsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyAnnotationFrontendGroupCondition>;
  filter: InputMaybe<StudyAnnotationFrontendGroupFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyAnnotationFrontendGroupsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyAnnotationFrontendGroupsListArgs = {
  condition: InputMaybe<StudyAnnotationFrontendGroupCondition>;
  filter: InputMaybe<StudyAnnotationFrontendGroupFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyAnnotationFrontendGroupsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyAnnotationFrontendValuesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyAnnotationFrontendValueCondition>;
  filter: InputMaybe<StudyAnnotationFrontendValueFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyAnnotationFrontendValuesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyAnnotationFrontendValuesListArgs = {
  condition: InputMaybe<StudyAnnotationFrontendValueCondition>;
  filter: InputMaybe<StudyAnnotationFrontendValueFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyAnnotationFrontendValuesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyAnnotationGroupUisArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyAnnotationGroupUiCondition>;
  filter: InputMaybe<StudyAnnotationGroupUiFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyAnnotationGroupUisOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyAnnotationGroupUisListArgs = {
  condition: InputMaybe<StudyAnnotationGroupUiCondition>;
  filter: InputMaybe<StudyAnnotationGroupUiFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyAnnotationGroupUisOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyImportLogsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyImportLogCondition>;
  filter: InputMaybe<StudyImportLogFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyImportLogsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyImportLogsListArgs = {
  condition: InputMaybe<StudyImportLogCondition>;
  filter: InputMaybe<StudyImportLogFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyImportLogsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyLayerArgs = {
  studyLayerId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyLayerByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyLayersArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyLayerCondition>;
  filter: InputMaybe<StudyLayerFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyLayersOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyLayersListArgs = {
  condition: InputMaybe<StudyLayerCondition>;
  filter: InputMaybe<StudyLayerFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyLayersOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyOmicsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyOmicCondition>;
  filter: InputMaybe<StudyOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyOmicsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyOmicsListArgs = {
  condition: InputMaybe<StudyOmicCondition>;
  filter: InputMaybe<StudyOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyOmicsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyOmicsTransposedsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyOmicsTransposedCondition>;
  filter: InputMaybe<StudyOmicsTransposedFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyOmicsTransposedsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyOmicsTransposedsListArgs = {
  condition: InputMaybe<StudyOmicsTransposedCondition>;
  filter: InputMaybe<StudyOmicsTransposedFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyOmicsTransposedsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyOverviewOntologiesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyOverviewOntologyCondition>;
  filter: InputMaybe<StudyOverviewOntologyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyOverviewOntologiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyOverviewOntologiesListArgs = {
  condition: InputMaybe<StudyOverviewOntologyCondition>;
  filter: InputMaybe<StudyOverviewOntologyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyOverviewOntologiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyOverviewsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyOverviewCondition>;
  filter: InputMaybe<StudyOverviewFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyOverviewsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyOverviewsListArgs = {
  condition: InputMaybe<StudyOverviewCondition>;
  filter: InputMaybe<StudyOverviewFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyOverviewsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySampleArgs = {
  studyId: Scalars['Int'];
  studySampleId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySampleAnnotationSubsamplingsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudySampleAnnotationSubsamplingCondition>;
  filter: InputMaybe<StudySampleAnnotationSubsamplingFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySampleAnnotationSubsamplingsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySampleAnnotationSubsamplingsListArgs = {
  condition: InputMaybe<StudySampleAnnotationSubsamplingCondition>;
  filter: InputMaybe<StudySampleAnnotationSubsamplingFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySampleAnnotationSubsamplingsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySampleAnnotationsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudySampleAnnotationCondition>;
  filter: InputMaybe<StudySampleAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySampleAnnotationsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySampleAnnotationsListArgs = {
  condition: InputMaybe<StudySampleAnnotationCondition>;
  filter: InputMaybe<StudySampleAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySampleAnnotationsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySampleByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySampleProjectionSubsamplingTransposedsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudySampleProjectionSubsamplingTransposedCondition>;
  filter: InputMaybe<StudySampleProjectionSubsamplingTransposedFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySampleProjectionSubsamplingTransposedsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySampleProjectionSubsamplingTransposedsListArgs = {
  condition: InputMaybe<StudySampleProjectionSubsamplingTransposedCondition>;
  filter: InputMaybe<StudySampleProjectionSubsamplingTransposedFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySampleProjectionSubsamplingTransposedsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySampleProjectionsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudySampleProjectionCondition>;
  filter: InputMaybe<StudySampleProjectionFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySampleProjectionsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySampleProjectionsListArgs = {
  condition: InputMaybe<StudySampleProjectionCondition>;
  filter: InputMaybe<StudySampleProjectionFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySampleProjectionsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySamplesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudySampleCondition>;
  filter: InputMaybe<StudySampleFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySamplesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySamplesListArgs = {
  condition: InputMaybe<StudySampleCondition>;
  filter: InputMaybe<StudySampleFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySamplesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyVisibleCurrentusersArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyVisibleCurrentuserCondition>;
  filter: InputMaybe<StudyVisibleCurrentuserFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyVisibleCurrentusersOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyVisibleCurrentusersListArgs = {
  condition: InputMaybe<StudyVisibleCurrentuserCondition>;
  filter: InputMaybe<StudyVisibleCurrentuserFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyVisibleCurrentusersOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryTreeOntologiesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<TreeOntologyCondition>;
  filter: InputMaybe<TreeOntologyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<TreeOntologiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryTreeOntologiesListArgs = {
  condition: InputMaybe<TreeOntologyCondition>;
  filter: InputMaybe<TreeOntologyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<TreeOntologiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryUserAnnotationGroupsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<UserAnnotationGroupCondition>;
  filter: InputMaybe<UserAnnotationGroupFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<UserAnnotationGroupsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryUserAnnotationGroupsListArgs = {
  condition: InputMaybe<UserAnnotationGroupCondition>;
  filter: InputMaybe<UserAnnotationGroupFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<UserAnnotationGroupsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryViolinPlotArgs = {
  pAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pExcludeAnnotationValueIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pOmicsId: InputMaybe<Scalars['Int']>;
  pSecondaryAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pStudyId: InputMaybe<Scalars['Int']>;
  pStudyLayerId: InputMaybe<Scalars['Int']>;
};

/** A connection to a list of `ReferenceStudy` values. */
export type ReferenceStudiesConnection = {
  __typename?: 'ReferenceStudiesConnection';
  /** A list of edges which contains the `ReferenceStudy` and cursor to aid in pagination. */
  edges: Array<ReferenceStudiesEdge>;
  /** A list of `ReferenceStudy` objects. */
  nodes: Array<Maybe<ReferenceStudy>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ReferenceStudy` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ReferenceStudy` edge in the connection. */
export type ReferenceStudiesEdge = {
  __typename?: 'ReferenceStudiesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ReferenceStudy` at the end of the edge. */
  node: Maybe<ReferenceStudy>;
};

/** Methods to use when ordering `ReferenceStudy`. */
export enum ReferenceStudiesOrderBy {
  CelltypeAnnotationGroupIdAsc = 'CELLTYPE_ANNOTATION_GROUP_ID_ASC',
  CelltypeAnnotationGroupIdDesc = 'CELLTYPE_ANNOTATION_GROUP_ID_DESC',
  Natural = 'NATURAL',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  TissueAnnotationGroupIdAsc = 'TISSUE_ANNOTATION_GROUP_ID_ASC',
  TissueAnnotationGroupIdDesc = 'TISSUE_ANNOTATION_GROUP_ID_DESC'
}

export type ReferenceStudy = {
  __typename?: 'ReferenceStudy';
  /** Reads a single `AnnotationGroup` that is related to this `ReferenceStudy`. */
  celltypeAnnotationGroup: Maybe<AnnotationGroup>;
  celltypeAnnotationGroupId: Scalars['Int'];
  /** Reads a single `ReferenceStudyOverview` that is related to this `ReferenceStudy`. */
  refStudy: Maybe<ReferenceStudyOverview>;
  /** Reads a single `Study` that is related to this `ReferenceStudy`. */
  study: Maybe<Study>;
  studyId: Scalars['Int'];
  /** Reads a single `AnnotationGroup` that is related to this `ReferenceStudy`. */
  tissueAnnotationGroup: Maybe<AnnotationGroup>;
  tissueAnnotationGroupId: Scalars['Int'];
};

/**
 * A condition to be used against `ReferenceStudy` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type ReferenceStudyCondition = {
  /** Checks for equality with the object’s `celltypeAnnotationGroupId` field. */
  celltypeAnnotationGroupId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `tissueAnnotationGroupId` field. */
  tissueAnnotationGroupId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `ReferenceStudy` object types. All fields are combined with a logical ‘and.’ */
export type ReferenceStudyFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ReferenceStudyFilter>>;
  /** Filter by the object’s `celltypeAnnotationGroupId` field. */
  celltypeAnnotationGroupId: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not: InputMaybe<ReferenceStudyFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ReferenceStudyFilter>>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `tissueAnnotationGroupId` field. */
  tissueAnnotationGroupId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `ReferenceStudy` */
export type ReferenceStudyInput = {
  celltypeAnnotationGroupId: Scalars['Int'];
  studyId: Scalars['Int'];
  tissueAnnotationGroupId: Scalars['Int'];
};

export type ReferenceStudyOverview = {
  __typename?: 'ReferenceStudyOverview';
  cellCount: Maybe<Scalars['Int']>;
  defaultStudyLayerId: Maybe<Scalars['Int']>;
  description: Maybe<Scalars['String']>;
  externalWebsite: Maybe<Scalars['String']>;
  metadata: Maybe<Scalars['JSON']>;
  organismTaxId: Maybe<Scalars['String']>;
  readerPermissions: Maybe<Array<Maybe<Scalars['String']>>>;
  /** Reads and enables pagination through a set of `ReferenceStudy`. */
  referenceStudyInfo: ReferenceStudiesConnection;
  /** Reads and enables pagination through a set of `ReferenceStudy`. */
  referenceStudyInfoList: Array<ReferenceStudy>;
  studyId: Maybe<Scalars['Int']>;
  studyName: Maybe<Scalars['String']>;
  /** Reads and enables pagination through a set of `StudyOverviewOntology`. */
  studyOntology: StudyOverviewOntologiesConnection;
  /** Reads and enables pagination through a set of `StudyOverviewOntology`. */
  studyOntologyList: Array<StudyOverviewOntology>;
};


export type ReferenceStudyOverviewReferenceStudyInfoArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ReferenceStudyCondition>;
  filter: InputMaybe<ReferenceStudyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ReferenceStudiesOrderBy>>;
};


export type ReferenceStudyOverviewReferenceStudyInfoListArgs = {
  condition: InputMaybe<ReferenceStudyCondition>;
  filter: InputMaybe<ReferenceStudyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ReferenceStudiesOrderBy>>;
};


export type ReferenceStudyOverviewStudyOntologyArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyOverviewOntologyCondition>;
  filter: InputMaybe<StudyOverviewOntologyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyOverviewOntologiesOrderBy>>;
};


export type ReferenceStudyOverviewStudyOntologyListArgs = {
  condition: InputMaybe<StudyOverviewOntologyCondition>;
  filter: InputMaybe<StudyOverviewOntologyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyOverviewOntologiesOrderBy>>;
};

/**
 * A condition to be used against `ReferenceStudyOverview` object types. All fields
 * are tested for equality and combined with a logical ‘and.’
 */
export type ReferenceStudyOverviewCondition = {
  /** Checks for equality with the object’s `cellCount` field. */
  cellCount: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `defaultStudyLayerId` field. */
  defaultStudyLayerId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `description` field. */
  description: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `metadata` field. */
  metadata: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `readerPermissions` field. */
  readerPermissions: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyName` field. */
  studyName: InputMaybe<Scalars['String']>;
};

/** A filter to be used against `ReferenceStudyOverview` object types. All fields are combined with a logical ‘and.’ */
export type ReferenceStudyOverviewFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ReferenceStudyOverviewFilter>>;
  /** Filter by the object’s `cellCount` field. */
  cellCount: InputMaybe<IntFilter>;
  /** Filter by the object’s `defaultStudyLayerId` field. */
  defaultStudyLayerId: InputMaybe<IntFilter>;
  /** Filter by the object’s `description` field. */
  description: InputMaybe<StringFilter>;
  /** Filter by the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<StringFilter>;
  /** Filter by the object’s `metadata` field. */
  metadata: InputMaybe<JsonFilter>;
  /** Negates the expression. */
  not: InputMaybe<ReferenceStudyOverviewFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ReferenceStudyOverviewFilter>>;
  /** Filter by the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<StringFilter>;
  /** Filter by the object’s `readerPermissions` field. */
  readerPermissions: InputMaybe<StringListFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyName` field. */
  studyName: InputMaybe<StringFilter>;
};

/** A connection to a list of `ReferenceStudyOverview` values. */
export type ReferenceStudyOverviewsConnection = {
  __typename?: 'ReferenceStudyOverviewsConnection';
  /** A list of edges which contains the `ReferenceStudyOverview` and cursor to aid in pagination. */
  edges: Array<ReferenceStudyOverviewsEdge>;
  /** A list of `ReferenceStudyOverview` objects. */
  nodes: Array<Maybe<ReferenceStudyOverview>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `ReferenceStudyOverview` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `ReferenceStudyOverview` edge in the connection. */
export type ReferenceStudyOverviewsEdge = {
  __typename?: 'ReferenceStudyOverviewsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `ReferenceStudyOverview` at the end of the edge. */
  node: Maybe<ReferenceStudyOverview>;
};

/** Methods to use when ordering `ReferenceStudyOverview`. */
export enum ReferenceStudyOverviewsOrderBy {
  CellCountAsc = 'CELL_COUNT_ASC',
  CellCountDesc = 'CELL_COUNT_DESC',
  DefaultStudyLayerIdAsc = 'DEFAULT_STUDY_LAYER_ID_ASC',
  DefaultStudyLayerIdDesc = 'DEFAULT_STUDY_LAYER_ID_DESC',
  DescriptionAsc = 'DESCRIPTION_ASC',
  DescriptionDesc = 'DESCRIPTION_DESC',
  ExternalWebsiteAsc = 'EXTERNAL_WEBSITE_ASC',
  ExternalWebsiteDesc = 'EXTERNAL_WEBSITE_DESC',
  MetadataAsc = 'METADATA_ASC',
  MetadataDesc = 'METADATA_DESC',
  Natural = 'NATURAL',
  OrganismTaxIdAsc = 'ORGANISM_TAX_ID_ASC',
  OrganismTaxIdDesc = 'ORGANISM_TAX_ID_DESC',
  ReaderPermissionsAsc = 'READER_PERMISSIONS_ASC',
  ReaderPermissionsDesc = 'READER_PERMISSIONS_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudyNameAsc = 'STUDY_NAME_ASC',
  StudyNameDesc = 'STUDY_NAME_DESC'
}

/** A filter to be used against String fields. All fields are combined with a logical ‘and.’ */
export type StringFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<Scalars['String']>;
  /** Not equal to the specified value, treating null like an ordinary value (case-insensitive). */
  distinctFromInsensitive: InputMaybe<Scalars['String']>;
  /** Ends with the specified string (case-sensitive). */
  endsWith: InputMaybe<Scalars['String']>;
  /** Ends with the specified string (case-insensitive). */
  endsWithInsensitive: InputMaybe<Scalars['String']>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<Scalars['String']>;
  /** Equal to the specified value (case-insensitive). */
  equalToInsensitive: InputMaybe<Scalars['String']>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<Scalars['String']>;
  /** Greater than the specified value (case-insensitive). */
  greaterThanInsensitive: InputMaybe<Scalars['String']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<Scalars['String']>;
  /** Greater than or equal to the specified value (case-insensitive). */
  greaterThanOrEqualToInsensitive: InputMaybe<Scalars['String']>;
  /** Included in the specified list. */
  in: InputMaybe<Array<Scalars['String']>>;
  /** Included in the specified list (case-insensitive). */
  inInsensitive: InputMaybe<Array<Scalars['String']>>;
  /** Contains the specified string (case-sensitive). */
  includes: InputMaybe<Scalars['String']>;
  /** Contains the specified string (case-insensitive). */
  includesInsensitive: InputMaybe<Scalars['String']>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<Scalars['String']>;
  /** Less than the specified value (case-insensitive). */
  lessThanInsensitive: InputMaybe<Scalars['String']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<Scalars['String']>;
  /** Less than or equal to the specified value (case-insensitive). */
  lessThanOrEqualToInsensitive: InputMaybe<Scalars['String']>;
  /** Matches the specified pattern (case-sensitive). An underscore (_) matches any single character; a percent sign (%) matches any sequence of zero or more characters. */
  like: InputMaybe<Scalars['String']>;
  /** Matches the specified pattern (case-insensitive). An underscore (_) matches any single character; a percent sign (%) matches any sequence of zero or more characters. */
  likeInsensitive: InputMaybe<Scalars['String']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<Scalars['String']>;
  /** Equal to the specified value, treating null like an ordinary value (case-insensitive). */
  notDistinctFromInsensitive: InputMaybe<Scalars['String']>;
  /** Does not end with the specified string (case-sensitive). */
  notEndsWith: InputMaybe<Scalars['String']>;
  /** Does not end with the specified string (case-insensitive). */
  notEndsWithInsensitive: InputMaybe<Scalars['String']>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<Scalars['String']>;
  /** Not equal to the specified value (case-insensitive). */
  notEqualToInsensitive: InputMaybe<Scalars['String']>;
  /** Not included in the specified list. */
  notIn: InputMaybe<Array<Scalars['String']>>;
  /** Not included in the specified list (case-insensitive). */
  notInInsensitive: InputMaybe<Array<Scalars['String']>>;
  /** Does not contain the specified string (case-sensitive). */
  notIncludes: InputMaybe<Scalars['String']>;
  /** Does not contain the specified string (case-insensitive). */
  notIncludesInsensitive: InputMaybe<Scalars['String']>;
  /** Does not match the specified pattern (case-sensitive). An underscore (_) matches any single character; a percent sign (%) matches any sequence of zero or more characters. */
  notLike: InputMaybe<Scalars['String']>;
  /** Does not match the specified pattern (case-insensitive). An underscore (_) matches any single character; a percent sign (%) matches any sequence of zero or more characters. */
  notLikeInsensitive: InputMaybe<Scalars['String']>;
  /** Does not start with the specified string (case-sensitive). */
  notStartsWith: InputMaybe<Scalars['String']>;
  /** Does not start with the specified string (case-insensitive). */
  notStartsWithInsensitive: InputMaybe<Scalars['String']>;
  /** Starts with the specified string (case-sensitive). */
  startsWith: InputMaybe<Scalars['String']>;
  /** Starts with the specified string (case-insensitive). */
  startsWithInsensitive: InputMaybe<Scalars['String']>;
};

/** A filter to be used against String List fields. All fields are combined with a logical ‘and.’ */
export type StringListFilter = {
  /** Any array item is equal to the specified value. */
  anyEqualTo: InputMaybe<Scalars['String']>;
  /** Any array item is greater than the specified value. */
  anyGreaterThan: InputMaybe<Scalars['String']>;
  /** Any array item is greater than or equal to the specified value. */
  anyGreaterThanOrEqualTo: InputMaybe<Scalars['String']>;
  /** Any array item is less than the specified value. */
  anyLessThan: InputMaybe<Scalars['String']>;
  /** Any array item is less than or equal to the specified value. */
  anyLessThanOrEqualTo: InputMaybe<Scalars['String']>;
  /** Any array item is not equal to the specified value. */
  anyNotEqualTo: InputMaybe<Scalars['String']>;
  /** Contained by the specified list of values. */
  containedBy: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Contains the specified list of values. */
  contains: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Overlaps the specified list of values. */
  overlaps: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** A connection to a list of `Study` values. */
export type StudiesConnection = {
  __typename?: 'StudiesConnection';
  /** A list of edges which contains the `Study` and cursor to aid in pagination. */
  edges: Array<StudiesEdge>;
  /** A list of `Study` objects. */
  nodes: Array<Maybe<Study>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `Study` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `Study` edge in the connection. */
export type StudiesEdge = {
  __typename?: 'StudiesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `Study` at the end of the edge. */
  node: Maybe<Study>;
};

/** Methods to use when ordering `Study`. */
export enum StudiesOrderBy {
  AdminPermissionsAsc = 'ADMIN_PERMISSIONS_ASC',
  AdminPermissionsDesc = 'ADMIN_PERMISSIONS_DESC',
  CellCountAsc = 'CELL_COUNT_ASC',
  CellCountDesc = 'CELL_COUNT_DESC',
  CellOntologyIdsAsc = 'CELL_ONTOLOGY_IDS_ASC',
  CellOntologyIdsDesc = 'CELL_ONTOLOGY_IDS_DESC',
  DescriptionAsc = 'DESCRIPTION_ASC',
  DescriptionDesc = 'DESCRIPTION_DESC',
  DiseaseMeshIdsAsc = 'DISEASE_MESH_IDS_ASC',
  DiseaseMeshIdsDesc = 'DISEASE_MESH_IDS_DESC',
  ExternalWebsiteAsc = 'EXTERNAL_WEBSITE_ASC',
  ExternalWebsiteDesc = 'EXTERNAL_WEBSITE_DESC',
  FilenameAsc = 'FILENAME_ASC',
  FilenameDesc = 'FILENAME_DESC',
  ImportFailedAsc = 'IMPORT_FAILED_ASC',
  ImportFailedDesc = 'IMPORT_FAILED_DESC',
  ImportFileAsc = 'IMPORT_FILE_ASC',
  ImportFileDesc = 'IMPORT_FILE_DESC',
  ImportFinishedAsc = 'IMPORT_FINISHED_ASC',
  ImportFinishedDesc = 'IMPORT_FINISHED_DESC',
  ImportLogAsc = 'IMPORT_LOG_ASC',
  ImportLogDesc = 'IMPORT_LOG_DESC',
  ImportStartedAsc = 'IMPORT_STARTED_ASC',
  ImportStartedDesc = 'IMPORT_STARTED_DESC',
  LegacyConfigAsc = 'LEGACY_CONFIG_ASC',
  LegacyConfigDesc = 'LEGACY_CONFIG_DESC',
  MetadataAsc = 'METADATA_ASC',
  MetadataDesc = 'METADATA_DESC',
  Natural = 'NATURAL',
  OrganismTaxIdAsc = 'ORGANISM_TAX_ID_ASC',
  OrganismTaxIdDesc = 'ORGANISM_TAX_ID_DESC',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  ProjectionsAsc = 'PROJECTIONS_ASC',
  ProjectionsDesc = 'PROJECTIONS_DESC',
  ReaderPermissionsAsc = 'READER_PERMISSIONS_ASC',
  ReaderPermissionsDesc = 'READER_PERMISSIONS_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudyNameAsc = 'STUDY_NAME_ASC',
  StudyNameDesc = 'STUDY_NAME_DESC',
  TissueNcitIdsAsc = 'TISSUE_NCIT_IDS_ASC',
  TissueNcitIdsDesc = 'TISSUE_NCIT_IDS_DESC',
  VisibleAsc = 'VISIBLE_ASC',
  VisibleDesc = 'VISIBLE_DESC'
}

export type Study = Node & {
  __typename?: 'Study';
  adminPermissions: Maybe<Array<Maybe<Scalars['String']>>>;
  cellCount: Maybe<Scalars['Int']>;
  cellOntologyIds: Maybe<Array<Maybe<Scalars['String']>>>;
  description: Maybe<Scalars['String']>;
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressions: DifferentialExpressionsConnection;
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressionsList: Array<DifferentialExpression>;
  diseaseMeshIds: Maybe<Array<Maybe<Scalars['String']>>>;
  externalWebsite: Maybe<Scalars['String']>;
  filename: Maybe<Scalars['String']>;
  importFailed: Maybe<Scalars['Boolean']>;
  importFile: Maybe<Scalars['String']>;
  importFinished: Maybe<Scalars['Boolean']>;
  importLog: Maybe<Scalars['String']>;
  importStarted: Maybe<Scalars['Boolean']>;
  legacyConfig: Maybe<Scalars['JSON']>;
  metadata: Maybe<Scalars['JSON']>;
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  organismTaxId: Maybe<Scalars['String']>;
  projections: Maybe<Array<Maybe<Scalars['String']>>>;
  readerPermissions: Maybe<Array<Maybe<Scalars['String']>>>;
  /** Reads and enables pagination through a set of `ReferenceStudy`. */
  referenceStudies: ReferenceStudiesConnection;
  /** Reads and enables pagination through a set of `ReferenceStudy`. */
  referenceStudiesList: Array<ReferenceStudy>;
  /** Reads and enables pagination through a set of `StudyAnnotationGroupUi`. */
  studyAnnotationGroupUis: StudyAnnotationGroupUisConnection;
  /** Reads and enables pagination through a set of `StudyAnnotationGroupUi`. */
  studyAnnotationGroupUisList: Array<StudyAnnotationGroupUi>;
  studyId: Scalars['Int'];
  /** Reads and enables pagination through a set of `StudyLayer`. */
  studyLayers: StudyLayersConnection;
  /** Reads and enables pagination through a set of `StudyLayer`. */
  studyLayersList: Array<StudyLayer>;
  studyName: Scalars['String'];
  /** Reads and enables pagination through a set of `StudyOmic`. */
  studyOmics: StudyOmicsConnection;
  /** Reads and enables pagination through a set of `StudyOmic`. */
  studyOmicsList: Array<StudyOmic>;
  /** Reads and enables pagination through a set of `StudySampleAnnotation`. */
  studySampleAnnotations: StudySampleAnnotationsConnection;
  /** Reads and enables pagination through a set of `StudySampleAnnotation`. */
  studySampleAnnotationsList: Array<StudySampleAnnotation>;
  /** Reads and enables pagination through a set of `StudySample`. */
  studySamples: StudySamplesConnection;
  /** Reads and enables pagination through a set of `StudySample`. */
  studySamplesList: Array<StudySample>;
  tissueNcitIds: Maybe<Array<Maybe<Scalars['String']>>>;
  /** Reads and enables pagination through a set of `UserAnnotationGroup`. */
  userAnnotationGroups: UserAnnotationGroupsConnection;
  /** Reads and enables pagination through a set of `UserAnnotationGroup`. */
  userAnnotationGroupsList: Array<UserAnnotationGroup>;
  visible: Maybe<Scalars['Boolean']>;
};


export type StudyDifferentialExpressionsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<DifferentialExpressionCondition>;
  filter: InputMaybe<DifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


export type StudyDifferentialExpressionsListArgs = {
  condition: InputMaybe<DifferentialExpressionCondition>;
  filter: InputMaybe<DifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


export type StudyReferenceStudiesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<ReferenceStudyCondition>;
  filter: InputMaybe<ReferenceStudyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ReferenceStudiesOrderBy>>;
};


export type StudyReferenceStudiesListArgs = {
  condition: InputMaybe<ReferenceStudyCondition>;
  filter: InputMaybe<ReferenceStudyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ReferenceStudiesOrderBy>>;
};


export type StudyStudyAnnotationGroupUisArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyAnnotationGroupUiCondition>;
  filter: InputMaybe<StudyAnnotationGroupUiFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyAnnotationGroupUisOrderBy>>;
};


export type StudyStudyAnnotationGroupUisListArgs = {
  condition: InputMaybe<StudyAnnotationGroupUiCondition>;
  filter: InputMaybe<StudyAnnotationGroupUiFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyAnnotationGroupUisOrderBy>>;
};


export type StudyStudyLayersArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyLayerCondition>;
  filter: InputMaybe<StudyLayerFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyLayersOrderBy>>;
};


export type StudyStudyLayersListArgs = {
  condition: InputMaybe<StudyLayerCondition>;
  filter: InputMaybe<StudyLayerFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyLayersOrderBy>>;
};


export type StudyStudyOmicsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyOmicCondition>;
  filter: InputMaybe<StudyOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyOmicsOrderBy>>;
};


export type StudyStudyOmicsListArgs = {
  condition: InputMaybe<StudyOmicCondition>;
  filter: InputMaybe<StudyOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyOmicsOrderBy>>;
};


export type StudyStudySampleAnnotationsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudySampleAnnotationCondition>;
  filter: InputMaybe<StudySampleAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySampleAnnotationsOrderBy>>;
};


export type StudyStudySampleAnnotationsListArgs = {
  condition: InputMaybe<StudySampleAnnotationCondition>;
  filter: InputMaybe<StudySampleAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySampleAnnotationsOrderBy>>;
};


export type StudyStudySamplesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudySampleCondition>;
  filter: InputMaybe<StudySampleFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySamplesOrderBy>>;
};


export type StudyStudySamplesListArgs = {
  condition: InputMaybe<StudySampleCondition>;
  filter: InputMaybe<StudySampleFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySamplesOrderBy>>;
};


export type StudyUserAnnotationGroupsArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<UserAnnotationGroupCondition>;
  filter: InputMaybe<UserAnnotationGroupFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<UserAnnotationGroupsOrderBy>>;
};


export type StudyUserAnnotationGroupsListArgs = {
  condition: InputMaybe<UserAnnotationGroupCondition>;
  filter: InputMaybe<UserAnnotationGroupFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<UserAnnotationGroupsOrderBy>>;
};

export type StudyAdminDetail = {
  __typename?: 'StudyAdminDetail';
  adminPermissionGranted: Maybe<Scalars['Boolean']>;
  adminPermissions: Maybe<Array<Maybe<Scalars['String']>>>;
  cellCount: Maybe<Scalars['Int']>;
  description: Maybe<Scalars['String']>;
  diseaseMeshIds: Maybe<Array<Maybe<Scalars['String']>>>;
  externalWebsite: Maybe<Scalars['String']>;
  filename: Maybe<Scalars['String']>;
  hasImportLog: Maybe<Scalars['Boolean']>;
  importFailed: Maybe<Scalars['Boolean']>;
  importFinished: Maybe<Scalars['Boolean']>;
  importStarted: Maybe<Scalars['Boolean']>;
  readerPermissionGranted: Maybe<Scalars['Boolean']>;
  readerPermissions: Maybe<Array<Maybe<Scalars['String']>>>;
  studyId: Maybe<Scalars['Int']>;
  studyName: Maybe<Scalars['String']>;
  tissueNcitIds: Maybe<Array<Maybe<Scalars['String']>>>;
  visible: Maybe<Scalars['Boolean']>;
};

/**
 * A condition to be used against `StudyAdminDetail` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type StudyAdminDetailCondition = {
  /** Checks for equality with the object’s `adminPermissionGranted` field. */
  adminPermissionGranted: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `adminPermissions` field. */
  adminPermissions: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `cellCount` field. */
  cellCount: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `description` field. */
  description: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `diseaseMeshIds` field. */
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `filename` field. */
  filename: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `hasImportLog` field. */
  hasImportLog: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `importFailed` field. */
  importFailed: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `importFinished` field. */
  importFinished: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `importStarted` field. */
  importStarted: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `readerPermissionGranted` field. */
  readerPermissionGranted: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `readerPermissions` field. */
  readerPermissions: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyName` field. */
  studyName: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `tissueNcitIds` field. */
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `visible` field. */
  visible: InputMaybe<Scalars['Boolean']>;
};

/** A filter to be used against `StudyAdminDetail` object types. All fields are combined with a logical ‘and.’ */
export type StudyAdminDetailFilter = {
  /** Filter by the object’s `adminPermissionGranted` field. */
  adminPermissionGranted: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `adminPermissions` field. */
  adminPermissions: InputMaybe<StringListFilter>;
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyAdminDetailFilter>>;
  /** Filter by the object’s `cellCount` field. */
  cellCount: InputMaybe<IntFilter>;
  /** Filter by the object’s `description` field. */
  description: InputMaybe<StringFilter>;
  /** Filter by the object’s `diseaseMeshIds` field. */
  diseaseMeshIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<StringFilter>;
  /** Filter by the object’s `filename` field. */
  filename: InputMaybe<StringFilter>;
  /** Filter by the object’s `hasImportLog` field. */
  hasImportLog: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `importFailed` field. */
  importFailed: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `importFinished` field. */
  importFinished: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `importStarted` field. */
  importStarted: InputMaybe<BooleanFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyAdminDetailFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyAdminDetailFilter>>;
  /** Filter by the object’s `readerPermissionGranted` field. */
  readerPermissionGranted: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `readerPermissions` field. */
  readerPermissions: InputMaybe<StringListFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyName` field. */
  studyName: InputMaybe<StringFilter>;
  /** Filter by the object’s `tissueNcitIds` field. */
  tissueNcitIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `visible` field. */
  visible: InputMaybe<BooleanFilter>;
};

/** A connection to a list of `StudyAdminDetail` values. */
export type StudyAdminDetailsConnection = {
  __typename?: 'StudyAdminDetailsConnection';
  /** A list of edges which contains the `StudyAdminDetail` and cursor to aid in pagination. */
  edges: Array<StudyAdminDetailsEdge>;
  /** A list of `StudyAdminDetail` objects. */
  nodes: Array<Maybe<StudyAdminDetail>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudyAdminDetail` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudyAdminDetail` edge in the connection. */
export type StudyAdminDetailsEdge = {
  __typename?: 'StudyAdminDetailsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudyAdminDetail` at the end of the edge. */
  node: Maybe<StudyAdminDetail>;
};

/** Methods to use when ordering `StudyAdminDetail`. */
export enum StudyAdminDetailsOrderBy {
  AdminPermissionsAsc = 'ADMIN_PERMISSIONS_ASC',
  AdminPermissionsDesc = 'ADMIN_PERMISSIONS_DESC',
  AdminPermissionGrantedAsc = 'ADMIN_PERMISSION_GRANTED_ASC',
  AdminPermissionGrantedDesc = 'ADMIN_PERMISSION_GRANTED_DESC',
  CellCountAsc = 'CELL_COUNT_ASC',
  CellCountDesc = 'CELL_COUNT_DESC',
  DescriptionAsc = 'DESCRIPTION_ASC',
  DescriptionDesc = 'DESCRIPTION_DESC',
  DiseaseMeshIdsAsc = 'DISEASE_MESH_IDS_ASC',
  DiseaseMeshIdsDesc = 'DISEASE_MESH_IDS_DESC',
  ExternalWebsiteAsc = 'EXTERNAL_WEBSITE_ASC',
  ExternalWebsiteDesc = 'EXTERNAL_WEBSITE_DESC',
  FilenameAsc = 'FILENAME_ASC',
  FilenameDesc = 'FILENAME_DESC',
  HasImportLogAsc = 'HAS_IMPORT_LOG_ASC',
  HasImportLogDesc = 'HAS_IMPORT_LOG_DESC',
  ImportFailedAsc = 'IMPORT_FAILED_ASC',
  ImportFailedDesc = 'IMPORT_FAILED_DESC',
  ImportFinishedAsc = 'IMPORT_FINISHED_ASC',
  ImportFinishedDesc = 'IMPORT_FINISHED_DESC',
  ImportStartedAsc = 'IMPORT_STARTED_ASC',
  ImportStartedDesc = 'IMPORT_STARTED_DESC',
  Natural = 'NATURAL',
  ReaderPermissionsAsc = 'READER_PERMISSIONS_ASC',
  ReaderPermissionsDesc = 'READER_PERMISSIONS_DESC',
  ReaderPermissionGrantedAsc = 'READER_PERMISSION_GRANTED_ASC',
  ReaderPermissionGrantedDesc = 'READER_PERMISSION_GRANTED_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudyNameAsc = 'STUDY_NAME_ASC',
  StudyNameDesc = 'STUDY_NAME_DESC',
  TissueNcitIdsAsc = 'TISSUE_NCIT_IDS_ASC',
  TissueNcitIdsDesc = 'TISSUE_NCIT_IDS_DESC',
  VisibleAsc = 'VISIBLE_ASC',
  VisibleDesc = 'VISIBLE_DESC'
}

export type StudyAdministrableCurrentuser = {
  __typename?: 'StudyAdministrableCurrentuser';
  studyId: Maybe<Scalars['Int']>;
};

/**
 * A condition to be used against `StudyAdministrableCurrentuser` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type StudyAdministrableCurrentuserCondition = {
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudyAdministrableCurrentuser` object types. All fields are combined with a logical ‘and.’ */
export type StudyAdministrableCurrentuserFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyAdministrableCurrentuserFilter>>;
  /** Negates the expression. */
  not: InputMaybe<StudyAdministrableCurrentuserFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyAdministrableCurrentuserFilter>>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `StudyAdministrableCurrentuser` */
export type StudyAdministrableCurrentuserInput = {
  studyId: InputMaybe<Scalars['Int']>;
};

/** A connection to a list of `StudyAdministrableCurrentuser` values. */
export type StudyAdministrableCurrentusersConnection = {
  __typename?: 'StudyAdministrableCurrentusersConnection';
  /** A list of edges which contains the `StudyAdministrableCurrentuser` and cursor to aid in pagination. */
  edges: Array<StudyAdministrableCurrentusersEdge>;
  /** A list of `StudyAdministrableCurrentuser` objects. */
  nodes: Array<Maybe<StudyAdministrableCurrentuser>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudyAdministrableCurrentuser` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudyAdministrableCurrentuser` edge in the connection. */
export type StudyAdministrableCurrentusersEdge = {
  __typename?: 'StudyAdministrableCurrentusersEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudyAdministrableCurrentuser` at the end of the edge. */
  node: Maybe<StudyAdministrableCurrentuser>;
};

/** Methods to use when ordering `StudyAdministrableCurrentuser`. */
export enum StudyAdministrableCurrentusersOrderBy {
  Natural = 'NATURAL',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

export type StudyAnnotationFrontendGroup = {
  __typename?: 'StudyAnnotationFrontendGroup';
  annotationGroupId: Maybe<Scalars['Int']>;
  /** Reads and enables pagination through a set of `StudyAnnotationFrontendValue`. */
  annotationValues: StudyAnnotationFrontendValuesConnection;
  /** Reads and enables pagination through a set of `StudyAnnotationFrontendValue`. */
  annotationValuesList: Array<StudyAnnotationFrontendValue>;
  createdByUser: Maybe<Scalars['String']>;
  currentUserIsOwner: Maybe<Scalars['Boolean']>;
  differentialExpressionCalculated: Maybe<Scalars['Boolean']>;
  displayGroup: Maybe<Scalars['String']>;
  isPrimary: Maybe<Scalars['Boolean']>;
  ordering: Maybe<Scalars['Int']>;
  privateToUser: Maybe<Scalars['Boolean']>;
  studyId: Maybe<Scalars['Int']>;
};


export type StudyAnnotationFrontendGroupAnnotationValuesArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyAnnotationFrontendValueCondition>;
  filter: InputMaybe<StudyAnnotationFrontendValueFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyAnnotationFrontendValuesOrderBy>>;
};


export type StudyAnnotationFrontendGroupAnnotationValuesListArgs = {
  condition: InputMaybe<StudyAnnotationFrontendValueCondition>;
  filter: InputMaybe<StudyAnnotationFrontendValueFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyAnnotationFrontendValuesOrderBy>>;
};

/**
 * A condition to be used against `StudyAnnotationFrontendGroup` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type StudyAnnotationFrontendGroupCondition = {
  /** Checks for equality with the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `createdByUser` field. */
  createdByUser: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `currentUserIsOwner` field. */
  currentUserIsOwner: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `differentialExpressionCalculated` field. */
  differentialExpressionCalculated: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `displayGroup` field. */
  displayGroup: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `isPrimary` field. */
  isPrimary: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `ordering` field. */
  ordering: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `privateToUser` field. */
  privateToUser: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudyAnnotationFrontendGroup` object types. All fields are combined with a logical ‘and.’ */
export type StudyAnnotationFrontendGroupFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyAnnotationFrontendGroupFilter>>;
  /** Filter by the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<IntFilter>;
  /** Filter by the object’s `createdByUser` field. */
  createdByUser: InputMaybe<StringFilter>;
  /** Filter by the object’s `currentUserIsOwner` field. */
  currentUserIsOwner: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `differentialExpressionCalculated` field. */
  differentialExpressionCalculated: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `displayGroup` field. */
  displayGroup: InputMaybe<StringFilter>;
  /** Filter by the object’s `isPrimary` field. */
  isPrimary: InputMaybe<BooleanFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyAnnotationFrontendGroupFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyAnnotationFrontendGroupFilter>>;
  /** Filter by the object’s `ordering` field. */
  ordering: InputMaybe<IntFilter>;
  /** Filter by the object’s `privateToUser` field. */
  privateToUser: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** A connection to a list of `StudyAnnotationFrontendGroup` values. */
export type StudyAnnotationFrontendGroupsConnection = {
  __typename?: 'StudyAnnotationFrontendGroupsConnection';
  /** A list of edges which contains the `StudyAnnotationFrontendGroup` and cursor to aid in pagination. */
  edges: Array<StudyAnnotationFrontendGroupsEdge>;
  /** A list of `StudyAnnotationFrontendGroup` objects. */
  nodes: Array<Maybe<StudyAnnotationFrontendGroup>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudyAnnotationFrontendGroup` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudyAnnotationFrontendGroup` edge in the connection. */
export type StudyAnnotationFrontendGroupsEdge = {
  __typename?: 'StudyAnnotationFrontendGroupsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudyAnnotationFrontendGroup` at the end of the edge. */
  node: Maybe<StudyAnnotationFrontendGroup>;
};

/** Methods to use when ordering `StudyAnnotationFrontendGroup`. */
export enum StudyAnnotationFrontendGroupsOrderBy {
  AnnotationGroupIdAsc = 'ANNOTATION_GROUP_ID_ASC',
  AnnotationGroupIdDesc = 'ANNOTATION_GROUP_ID_DESC',
  CreatedByUserAsc = 'CREATED_BY_USER_ASC',
  CreatedByUserDesc = 'CREATED_BY_USER_DESC',
  CurrentUserIsOwnerAsc = 'CURRENT_USER_IS_OWNER_ASC',
  CurrentUserIsOwnerDesc = 'CURRENT_USER_IS_OWNER_DESC',
  DifferentialExpressionCalculatedAsc = 'DIFFERENTIAL_EXPRESSION_CALCULATED_ASC',
  DifferentialExpressionCalculatedDesc = 'DIFFERENTIAL_EXPRESSION_CALCULATED_DESC',
  DisplayGroupAsc = 'DISPLAY_GROUP_ASC',
  DisplayGroupDesc = 'DISPLAY_GROUP_DESC',
  IsPrimaryAsc = 'IS_PRIMARY_ASC',
  IsPrimaryDesc = 'IS_PRIMARY_DESC',
  Natural = 'NATURAL',
  OrderingAsc = 'ORDERING_ASC',
  OrderingDesc = 'ORDERING_DESC',
  PrivateToUserAsc = 'PRIVATE_TO_USER_ASC',
  PrivateToUserDesc = 'PRIVATE_TO_USER_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

export type StudyAnnotationFrontendValue = {
  __typename?: 'StudyAnnotationFrontendValue';
  annotationGroupId: Maybe<Scalars['Int']>;
  annotationValueId: Maybe<Scalars['Int']>;
  color: Maybe<Scalars['String']>;
  displayValue: Maybe<Scalars['String']>;
  /** Reads a single `StudyAnnotationFrontendGroup` that is related to this `StudyAnnotationFrontendValue`. */
  group: Maybe<StudyAnnotationFrontendGroup>;
  sampleCount: Maybe<Scalars['Int']>;
  studyId: Maybe<Scalars['Int']>;
};

/**
 * A condition to be used against `StudyAnnotationFrontendValue` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type StudyAnnotationFrontendValueCondition = {
  /** Checks for equality with the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `color` field. */
  color: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `displayValue` field. */
  displayValue: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `sampleCount` field. */
  sampleCount: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudyAnnotationFrontendValue` object types. All fields are combined with a logical ‘and.’ */
export type StudyAnnotationFrontendValueFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyAnnotationFrontendValueFilter>>;
  /** Filter by the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<IntFilter>;
  /** Filter by the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<IntFilter>;
  /** Filter by the object’s `color` field. */
  color: InputMaybe<StringFilter>;
  /** Filter by the object’s `displayValue` field. */
  displayValue: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyAnnotationFrontendValueFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyAnnotationFrontendValueFilter>>;
  /** Filter by the object’s `sampleCount` field. */
  sampleCount: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** A connection to a list of `StudyAnnotationFrontendValue` values. */
export type StudyAnnotationFrontendValuesConnection = {
  __typename?: 'StudyAnnotationFrontendValuesConnection';
  /** A list of edges which contains the `StudyAnnotationFrontendValue` and cursor to aid in pagination. */
  edges: Array<StudyAnnotationFrontendValuesEdge>;
  /** A list of `StudyAnnotationFrontendValue` objects. */
  nodes: Array<Maybe<StudyAnnotationFrontendValue>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudyAnnotationFrontendValue` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudyAnnotationFrontendValue` edge in the connection. */
export type StudyAnnotationFrontendValuesEdge = {
  __typename?: 'StudyAnnotationFrontendValuesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudyAnnotationFrontendValue` at the end of the edge. */
  node: Maybe<StudyAnnotationFrontendValue>;
};

/** Methods to use when ordering `StudyAnnotationFrontendValue`. */
export enum StudyAnnotationFrontendValuesOrderBy {
  AnnotationGroupIdAsc = 'ANNOTATION_GROUP_ID_ASC',
  AnnotationGroupIdDesc = 'ANNOTATION_GROUP_ID_DESC',
  AnnotationValueIdAsc = 'ANNOTATION_VALUE_ID_ASC',
  AnnotationValueIdDesc = 'ANNOTATION_VALUE_ID_DESC',
  ColorAsc = 'COLOR_ASC',
  ColorDesc = 'COLOR_DESC',
  DisplayValueAsc = 'DISPLAY_VALUE_ASC',
  DisplayValueDesc = 'DISPLAY_VALUE_DESC',
  Natural = 'NATURAL',
  SampleCountAsc = 'SAMPLE_COUNT_ASC',
  SampleCountDesc = 'SAMPLE_COUNT_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

export type StudyAnnotationGroupUi = {
  __typename?: 'StudyAnnotationGroupUi';
  /** Reads a single `AnnotationGroup` that is related to this `StudyAnnotationGroupUi`. */
  annotationGroup: Maybe<AnnotationGroup>;
  annotationGroupId: Scalars['Int'];
  differentialExpressionCalculated: Scalars['Boolean'];
  isPrimary: Scalars['Boolean'];
  ordering: Scalars['Int'];
  /** Reads a single `Study` that is related to this `StudyAnnotationGroupUi`. */
  study: Maybe<Study>;
  studyId: Scalars['Int'];
};

/**
 * A condition to be used against `StudyAnnotationGroupUi` object types. All fields
 * are tested for equality and combined with a logical ‘and.’
 */
export type StudyAnnotationGroupUiCondition = {
  /** Checks for equality with the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `differentialExpressionCalculated` field. */
  differentialExpressionCalculated: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `isPrimary` field. */
  isPrimary: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `ordering` field. */
  ordering: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudyAnnotationGroupUi` object types. All fields are combined with a logical ‘and.’ */
export type StudyAnnotationGroupUiFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyAnnotationGroupUiFilter>>;
  /** Filter by the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<IntFilter>;
  /** Filter by the object’s `differentialExpressionCalculated` field. */
  differentialExpressionCalculated: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `isPrimary` field. */
  isPrimary: InputMaybe<BooleanFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyAnnotationGroupUiFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyAnnotationGroupUiFilter>>;
  /** Filter by the object’s `ordering` field. */
  ordering: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `StudyAnnotationGroupUi` */
export type StudyAnnotationGroupUiInput = {
  annotationGroupId: Scalars['Int'];
  differentialExpressionCalculated: Scalars['Boolean'];
  isPrimary: Scalars['Boolean'];
  ordering: Scalars['Int'];
  studyId: Scalars['Int'];
};

/** A connection to a list of `StudyAnnotationGroupUi` values. */
export type StudyAnnotationGroupUisConnection = {
  __typename?: 'StudyAnnotationGroupUisConnection';
  /** A list of edges which contains the `StudyAnnotationGroupUi` and cursor to aid in pagination. */
  edges: Array<StudyAnnotationGroupUisEdge>;
  /** A list of `StudyAnnotationGroupUi` objects. */
  nodes: Array<Maybe<StudyAnnotationGroupUi>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudyAnnotationGroupUi` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudyAnnotationGroupUi` edge in the connection. */
export type StudyAnnotationGroupUisEdge = {
  __typename?: 'StudyAnnotationGroupUisEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudyAnnotationGroupUi` at the end of the edge. */
  node: Maybe<StudyAnnotationGroupUi>;
};

/** Methods to use when ordering `StudyAnnotationGroupUi`. */
export enum StudyAnnotationGroupUisOrderBy {
  AnnotationGroupIdAsc = 'ANNOTATION_GROUP_ID_ASC',
  AnnotationGroupIdDesc = 'ANNOTATION_GROUP_ID_DESC',
  DifferentialExpressionCalculatedAsc = 'DIFFERENTIAL_EXPRESSION_CALCULATED_ASC',
  DifferentialExpressionCalculatedDesc = 'DIFFERENTIAL_EXPRESSION_CALCULATED_DESC',
  IsPrimaryAsc = 'IS_PRIMARY_ASC',
  IsPrimaryDesc = 'IS_PRIMARY_DESC',
  Natural = 'NATURAL',
  OrderingAsc = 'ORDERING_ASC',
  OrderingDesc = 'ORDERING_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

/** A condition to be used against `Study` object types. All fields are tested for equality and combined with a logical ‘and.’ */
export type StudyCondition = {
  /** Checks for equality with the object’s `adminPermissions` field. */
  adminPermissions: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `cellCount` field. */
  cellCount: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `cellOntologyIds` field. */
  cellOntologyIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `description` field. */
  description: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `diseaseMeshIds` field. */
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `filename` field. */
  filename: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `importFailed` field. */
  importFailed: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `importFile` field. */
  importFile: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `importFinished` field. */
  importFinished: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `importLog` field. */
  importLog: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `importStarted` field. */
  importStarted: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `legacyConfig` field. */
  legacyConfig: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `metadata` field. */
  metadata: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `projections` field. */
  projections: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `readerPermissions` field. */
  readerPermissions: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyName` field. */
  studyName: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `tissueNcitIds` field. */
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `visible` field. */
  visible: InputMaybe<Scalars['Boolean']>;
};

/** All input for the `studyDefinitionUpdate` mutation. */
export type StudyDefinitionUpdateInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
};

/** The output of our `studyDefinitionUpdate` mutation. */
export type StudyDefinitionUpdatePayload = {
  __typename?: 'StudyDefinitionUpdatePayload';
  boolean: Maybe<Scalars['Boolean']>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};

/** A filter to be used against `Study` object types. All fields are combined with a logical ‘and.’ */
export type StudyFilter = {
  /** Filter by the object’s `adminPermissions` field. */
  adminPermissions: InputMaybe<StringListFilter>;
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyFilter>>;
  /** Filter by the object’s `cellCount` field. */
  cellCount: InputMaybe<IntFilter>;
  /** Filter by the object’s `cellOntologyIds` field. */
  cellOntologyIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `description` field. */
  description: InputMaybe<StringFilter>;
  /** Filter by the object’s `diseaseMeshIds` field. */
  diseaseMeshIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<StringFilter>;
  /** Filter by the object’s `filename` field. */
  filename: InputMaybe<StringFilter>;
  /** Filter by the object’s `importFailed` field. */
  importFailed: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `importFile` field. */
  importFile: InputMaybe<StringFilter>;
  /** Filter by the object’s `importFinished` field. */
  importFinished: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `importLog` field. */
  importLog: InputMaybe<StringFilter>;
  /** Filter by the object’s `importStarted` field. */
  importStarted: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `legacyConfig` field. */
  legacyConfig: InputMaybe<JsonFilter>;
  /** Filter by the object’s `metadata` field. */
  metadata: InputMaybe<JsonFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyFilter>>;
  /** Filter by the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<StringFilter>;
  /** Filter by the object’s `projections` field. */
  projections: InputMaybe<StringListFilter>;
  /** Filter by the object’s `readerPermissions` field. */
  readerPermissions: InputMaybe<StringListFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyName` field. */
  studyName: InputMaybe<StringFilter>;
  /** Filter by the object’s `tissueNcitIds` field. */
  tissueNcitIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `visible` field. */
  visible: InputMaybe<BooleanFilter>;
};

export type StudyImportLog = {
  __typename?: 'StudyImportLog';
  importFile: Maybe<Scalars['String']>;
  importLog: Maybe<Scalars['String']>;
  studyId: Maybe<Scalars['Int']>;
};

/**
 * A condition to be used against `StudyImportLog` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type StudyImportLogCondition = {
  /** Checks for equality with the object’s `importFile` field. */
  importFile: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `importLog` field. */
  importLog: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudyImportLog` object types. All fields are combined with a logical ‘and.’ */
export type StudyImportLogFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyImportLogFilter>>;
  /** Filter by the object’s `importFile` field. */
  importFile: InputMaybe<StringFilter>;
  /** Filter by the object’s `importLog` field. */
  importLog: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyImportLogFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyImportLogFilter>>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** A connection to a list of `StudyImportLog` values. */
export type StudyImportLogsConnection = {
  __typename?: 'StudyImportLogsConnection';
  /** A list of edges which contains the `StudyImportLog` and cursor to aid in pagination. */
  edges: Array<StudyImportLogsEdge>;
  /** A list of `StudyImportLog` objects. */
  nodes: Array<Maybe<StudyImportLog>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudyImportLog` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudyImportLog` edge in the connection. */
export type StudyImportLogsEdge = {
  __typename?: 'StudyImportLogsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudyImportLog` at the end of the edge. */
  node: Maybe<StudyImportLog>;
};

/** Methods to use when ordering `StudyImportLog`. */
export enum StudyImportLogsOrderBy {
  ImportFileAsc = 'IMPORT_FILE_ASC',
  ImportFileDesc = 'IMPORT_FILE_DESC',
  ImportLogAsc = 'IMPORT_LOG_ASC',
  ImportLogDesc = 'IMPORT_LOG_DESC',
  Natural = 'NATURAL',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

/** An input for mutations affecting `Study` */
export type StudyInput = {
  adminPermissions: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  cellCount: InputMaybe<Scalars['Int']>;
  cellOntologyIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  description: InputMaybe<Scalars['String']>;
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  externalWebsite: InputMaybe<Scalars['String']>;
  filename: InputMaybe<Scalars['String']>;
  importFailed: InputMaybe<Scalars['Boolean']>;
  importFile: InputMaybe<Scalars['String']>;
  importFinished: InputMaybe<Scalars['Boolean']>;
  importLog: InputMaybe<Scalars['String']>;
  importStarted: InputMaybe<Scalars['Boolean']>;
  legacyConfig: InputMaybe<Scalars['JSON']>;
  metadata: InputMaybe<Scalars['JSON']>;
  organismTaxId: InputMaybe<Scalars['String']>;
  projections: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  readerPermissions: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  studyId: InputMaybe<Scalars['Int']>;
  studyName: Scalars['String'];
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  visible: InputMaybe<Scalars['Boolean']>;
};

export type StudyLayer = Node & {
  __typename?: 'StudyLayer';
  layer: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  omicsType: Maybe<OmicsType>;
  /** Reads a single `Study` that is related to this `StudyLayer`. */
  study: Maybe<Study>;
  studyId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
};

/**
 * A condition to be used against `StudyLayer` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type StudyLayerCondition = {
  /** Checks for equality with the object’s `layer` field. */
  layer: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `omicsType` field. */
  omicsType: InputMaybe<OmicsType>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudyLayer` object types. All fields are combined with a logical ‘and.’ */
export type StudyLayerFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyLayerFilter>>;
  /** Filter by the object’s `layer` field. */
  layer: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyLayerFilter>;
  /** Filter by the object’s `omicsType` field. */
  omicsType: InputMaybe<OmicsTypeFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyLayerFilter>>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `StudyLayer` */
export type StudyLayerInput = {
  layer: Scalars['String'];
  omicsType: InputMaybe<OmicsType>;
  studyId: Scalars['Int'];
  studyLayerId: InputMaybe<Scalars['Int']>;
};

/** Represents an update to a `StudyLayer`. Fields that are set will be updated. */
export type StudyLayerPatch = {
  layer: InputMaybe<Scalars['String']>;
  omicsType: InputMaybe<OmicsType>;
  studyId: InputMaybe<Scalars['Int']>;
  studyLayerId: InputMaybe<Scalars['Int']>;
};

/** A connection to a list of `StudyLayer` values. */
export type StudyLayersConnection = {
  __typename?: 'StudyLayersConnection';
  /** A list of edges which contains the `StudyLayer` and cursor to aid in pagination. */
  edges: Array<StudyLayersEdge>;
  /** A list of `StudyLayer` objects. */
  nodes: Array<Maybe<StudyLayer>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudyLayer` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudyLayer` edge in the connection. */
export type StudyLayersEdge = {
  __typename?: 'StudyLayersEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudyLayer` at the end of the edge. */
  node: Maybe<StudyLayer>;
};

/** Methods to use when ordering `StudyLayer`. */
export enum StudyLayersOrderBy {
  LayerAsc = 'LAYER_ASC',
  LayerDesc = 'LAYER_DESC',
  Natural = 'NATURAL',
  OmicsTypeAsc = 'OMICS_TYPE_ASC',
  OmicsTypeDesc = 'OMICS_TYPE_DESC',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudyLayerIdAsc = 'STUDY_LAYER_ID_ASC',
  StudyLayerIdDesc = 'STUDY_LAYER_ID_DESC'
}

export type StudyOmic = {
  __typename?: 'StudyOmic';
  h5AdVarIndex: Scalars['Int'];
  /** Reads a single `OmicsBase` that is related to this `StudyOmic`. */
  omics: Maybe<OmicsBase>;
  omicsId: Scalars['Int'];
  regionEnd: Maybe<Scalars['Int']>;
  regionStart: Maybe<Scalars['Int']>;
  /** Reads a single `Study` that is related to this `StudyOmic`. */
  study: Maybe<Study>;
  studyId: Scalars['Int'];
};

/**
 * A condition to be used against `StudyOmic` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type StudyOmicCondition = {
  /** Checks for equality with the object’s `h5AdVarIndex` field. */
  h5AdVarIndex: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `regionEnd` field. */
  regionEnd: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `regionStart` field. */
  regionStart: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudyOmic` object types. All fields are combined with a logical ‘and.’ */
export type StudyOmicFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyOmicFilter>>;
  /** Filter by the object’s `h5AdVarIndex` field. */
  h5AdVarIndex: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyOmicFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyOmicFilter>>;
  /** Filter by the object’s `regionEnd` field. */
  regionEnd: InputMaybe<IntFilter>;
  /** Filter by the object’s `regionStart` field. */
  regionStart: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `StudyOmic` */
export type StudyOmicInput = {
  h5AdVarIndex: Scalars['Int'];
  omicsId: Scalars['Int'];
  regionEnd: InputMaybe<Scalars['Int']>;
  regionStart: InputMaybe<Scalars['Int']>;
  studyId: Scalars['Int'];
};

/** A connection to a list of `StudyOmic` values. */
export type StudyOmicsConnection = {
  __typename?: 'StudyOmicsConnection';
  /** A list of edges which contains the `StudyOmic` and cursor to aid in pagination. */
  edges: Array<StudyOmicsEdge>;
  /** A list of `StudyOmic` objects. */
  nodes: Array<Maybe<StudyOmic>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudyOmic` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudyOmic` edge in the connection. */
export type StudyOmicsEdge = {
  __typename?: 'StudyOmicsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudyOmic` at the end of the edge. */
  node: Maybe<StudyOmic>;
};

/** Methods to use when ordering `StudyOmic`. */
export enum StudyOmicsOrderBy {
  H5AdVarIndexAsc = 'H5AD_VAR_INDEX_ASC',
  H5AdVarIndexDesc = 'H5AD_VAR_INDEX_DESC',
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  RegionEndAsc = 'REGION_END_ASC',
  RegionEndDesc = 'REGION_END_DESC',
  RegionStartAsc = 'REGION_START_ASC',
  RegionStartDesc = 'REGION_START_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

export type StudyOmicsTransposed = {
  __typename?: 'StudyOmicsTransposed';
  displayName: Maybe<Array<Maybe<Scalars['String']>>>;
  displaySymbol: Maybe<Array<Maybe<Scalars['String']>>>;
  omicsId: Maybe<Array<Maybe<Scalars['Int']>>>;
  omicsType: Maybe<Array<Maybe<OmicsType>>>;
  rowNumber: Maybe<Scalars['BigInt']>;
  studyId: Maybe<Scalars['Int']>;
};

/**
 * A condition to be used against `StudyOmicsTransposed` object types. All fields
 * are tested for equality and combined with a logical ‘and.’
 */
export type StudyOmicsTransposedCondition = {
  /** Checks for equality with the object’s `displayName` field. */
  displayName: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `omicsType` field. */
  omicsType: InputMaybe<Array<InputMaybe<OmicsType>>>;
  /** Checks for equality with the object’s `rowNumber` field. */
  rowNumber: InputMaybe<Scalars['BigInt']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudyOmicsTransposed` object types. All fields are combined with a logical ‘and.’ */
export type StudyOmicsTransposedFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyOmicsTransposedFilter>>;
  /** Filter by the object’s `displayName` field. */
  displayName: InputMaybe<StringListFilter>;
  /** Filter by the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<StringListFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyOmicsTransposedFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntListFilter>;
  /** Filter by the object’s `omicsType` field. */
  omicsType: InputMaybe<OmicsTypeListFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyOmicsTransposedFilter>>;
  /** Filter by the object’s `rowNumber` field. */
  rowNumber: InputMaybe<BigIntFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** A connection to a list of `StudyOmicsTransposed` values. */
export type StudyOmicsTransposedsConnection = {
  __typename?: 'StudyOmicsTransposedsConnection';
  /** A list of edges which contains the `StudyOmicsTransposed` and cursor to aid in pagination. */
  edges: Array<StudyOmicsTransposedsEdge>;
  /** A list of `StudyOmicsTransposed` objects. */
  nodes: Array<Maybe<StudyOmicsTransposed>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudyOmicsTransposed` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudyOmicsTransposed` edge in the connection. */
export type StudyOmicsTransposedsEdge = {
  __typename?: 'StudyOmicsTransposedsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudyOmicsTransposed` at the end of the edge. */
  node: Maybe<StudyOmicsTransposed>;
};

/** Methods to use when ordering `StudyOmicsTransposed`. */
export enum StudyOmicsTransposedsOrderBy {
  DisplayNameAsc = 'DISPLAY_NAME_ASC',
  DisplayNameDesc = 'DISPLAY_NAME_DESC',
  DisplaySymbolAsc = 'DISPLAY_SYMBOL_ASC',
  DisplaySymbolDesc = 'DISPLAY_SYMBOL_DESC',
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  OmicsTypeAsc = 'OMICS_TYPE_ASC',
  OmicsTypeDesc = 'OMICS_TYPE_DESC',
  RowNumberAsc = 'ROW_NUMBER_ASC',
  RowNumberDesc = 'ROW_NUMBER_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

export type StudyOverview = {
  __typename?: 'StudyOverview';
  cellCount: Maybe<Scalars['Int']>;
  defaultStudyLayerId: Maybe<Scalars['Int']>;
  description: Maybe<Scalars['String']>;
  externalWebsite: Maybe<Scalars['String']>;
  metadata: Maybe<Scalars['JSON']>;
  organismTaxId: Maybe<Scalars['String']>;
  studyId: Maybe<Scalars['Int']>;
  studyName: Maybe<Scalars['String']>;
  /** Reads and enables pagination through a set of `StudyOverviewOntology`. */
  studyOntology: StudyOverviewOntologiesConnection;
  /** Reads and enables pagination through a set of `StudyOverviewOntology`. */
  studyOntologyList: Array<StudyOverviewOntology>;
};


export type StudyOverviewStudyOntologyArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudyOverviewOntologyCondition>;
  filter: InputMaybe<StudyOverviewOntologyFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyOverviewOntologiesOrderBy>>;
};


export type StudyOverviewStudyOntologyListArgs = {
  condition: InputMaybe<StudyOverviewOntologyCondition>;
  filter: InputMaybe<StudyOverviewOntologyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyOverviewOntologiesOrderBy>>;
};

/**
 * A condition to be used against `StudyOverview` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type StudyOverviewCondition = {
  /** Checks for equality with the object’s `cellCount` field. */
  cellCount: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `defaultStudyLayerId` field. */
  defaultStudyLayerId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `description` field. */
  description: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `metadata` field. */
  metadata: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyName` field. */
  studyName: InputMaybe<Scalars['String']>;
};

/** A filter to be used against `StudyOverview` object types. All fields are combined with a logical ‘and.’ */
export type StudyOverviewFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyOverviewFilter>>;
  /** Filter by the object’s `cellCount` field. */
  cellCount: InputMaybe<IntFilter>;
  /** Filter by the object’s `defaultStudyLayerId` field. */
  defaultStudyLayerId: InputMaybe<IntFilter>;
  /** Filter by the object’s `description` field. */
  description: InputMaybe<StringFilter>;
  /** Filter by the object’s `externalWebsite` field. */
  externalWebsite: InputMaybe<StringFilter>;
  /** Filter by the object’s `metadata` field. */
  metadata: InputMaybe<JsonFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyOverviewFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyOverviewFilter>>;
  /** Filter by the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<StringFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyName` field. */
  studyName: InputMaybe<StringFilter>;
};

/** An input for mutations affecting `StudyOverview` */
export type StudyOverviewInput = {
  cellCount: InputMaybe<Scalars['Int']>;
  defaultStudyLayerId: InputMaybe<Scalars['Int']>;
  description: InputMaybe<Scalars['String']>;
  externalWebsite: InputMaybe<Scalars['String']>;
  metadata: InputMaybe<Scalars['JSON']>;
  organismTaxId: InputMaybe<Scalars['String']>;
  studyId: InputMaybe<Scalars['Int']>;
  studyName: InputMaybe<Scalars['String']>;
};

/** A connection to a list of `StudyOverviewOntology` values. */
export type StudyOverviewOntologiesConnection = {
  __typename?: 'StudyOverviewOntologiesConnection';
  /** A list of edges which contains the `StudyOverviewOntology` and cursor to aid in pagination. */
  edges: Array<StudyOverviewOntologiesEdge>;
  /** A list of `StudyOverviewOntology` objects. */
  nodes: Array<Maybe<StudyOverviewOntology>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudyOverviewOntology` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudyOverviewOntology` edge in the connection. */
export type StudyOverviewOntologiesEdge = {
  __typename?: 'StudyOverviewOntologiesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudyOverviewOntology` at the end of the edge. */
  node: Maybe<StudyOverviewOntology>;
};

/** Methods to use when ordering `StudyOverviewOntology`. */
export enum StudyOverviewOntologiesOrderBy {
  LabelsAsc = 'LABELS_ASC',
  LabelsDesc = 'LABELS_DESC',
  Natural = 'NATURAL',
  OntologyAsc = 'ONTOLOGY_ASC',
  OntologyDesc = 'ONTOLOGY_DESC',
  OntCodesAsc = 'ONT_CODES_ASC',
  OntCodesDesc = 'ONT_CODES_DESC',
  ParentIdsAsc = 'PARENT_IDS_ASC',
  ParentIdsDesc = 'PARENT_IDS_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

export type StudyOverviewOntology = {
  __typename?: 'StudyOverviewOntology';
  labels: Maybe<Array<Maybe<Scalars['String']>>>;
  ontCodes: Maybe<Array<Maybe<Scalars['String']>>>;
  ontology: Maybe<Scalars['String']>;
  parentIds: Maybe<Array<Maybe<Scalars['String']>>>;
  /** Reads a single `StudyOverview` that is related to this `StudyOverviewOntology`. */
  study: Maybe<StudyOverview>;
  studyId: Maybe<Scalars['Int']>;
  /** Reads a single `ReferenceStudyOverview` that is related to this `StudyOverviewOntology`. */
  studyOntology: Maybe<ReferenceStudyOverview>;
};

/**
 * A condition to be used against `StudyOverviewOntology` object types. All fields
 * are tested for equality and combined with a logical ‘and.’
 */
export type StudyOverviewOntologyCondition = {
  /** Checks for equality with the object’s `labels` field. */
  labels: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `ontCodes` field. */
  ontCodes: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `ontology` field. */
  ontology: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `parentIds` field. */
  parentIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudyOverviewOntology` object types. All fields are combined with a logical ‘and.’ */
export type StudyOverviewOntologyFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyOverviewOntologyFilter>>;
  /** Filter by the object’s `labels` field. */
  labels: InputMaybe<StringListFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyOverviewOntologyFilter>;
  /** Filter by the object’s `ontCodes` field. */
  ontCodes: InputMaybe<StringListFilter>;
  /** Filter by the object’s `ontology` field. */
  ontology: InputMaybe<StringFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyOverviewOntologyFilter>>;
  /** Filter by the object’s `parentIds` field. */
  parentIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** A connection to a list of `StudyOverview` values. */
export type StudyOverviewsConnection = {
  __typename?: 'StudyOverviewsConnection';
  /** A list of edges which contains the `StudyOverview` and cursor to aid in pagination. */
  edges: Array<StudyOverviewsEdge>;
  /** A list of `StudyOverview` objects. */
  nodes: Array<Maybe<StudyOverview>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudyOverview` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudyOverview` edge in the connection. */
export type StudyOverviewsEdge = {
  __typename?: 'StudyOverviewsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudyOverview` at the end of the edge. */
  node: Maybe<StudyOverview>;
};

/** Methods to use when ordering `StudyOverview`. */
export enum StudyOverviewsOrderBy {
  CellCountAsc = 'CELL_COUNT_ASC',
  CellCountDesc = 'CELL_COUNT_DESC',
  DefaultStudyLayerIdAsc = 'DEFAULT_STUDY_LAYER_ID_ASC',
  DefaultStudyLayerIdDesc = 'DEFAULT_STUDY_LAYER_ID_DESC',
  DescriptionAsc = 'DESCRIPTION_ASC',
  DescriptionDesc = 'DESCRIPTION_DESC',
  ExternalWebsiteAsc = 'EXTERNAL_WEBSITE_ASC',
  ExternalWebsiteDesc = 'EXTERNAL_WEBSITE_DESC',
  MetadataAsc = 'METADATA_ASC',
  MetadataDesc = 'METADATA_DESC',
  Natural = 'NATURAL',
  OrganismTaxIdAsc = 'ORGANISM_TAX_ID_ASC',
  OrganismTaxIdDesc = 'ORGANISM_TAX_ID_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudyNameAsc = 'STUDY_NAME_ASC',
  StudyNameDesc = 'STUDY_NAME_DESC'
}

/** Represents an update to a `Study`. Fields that are set will be updated. */
export type StudyPatch = {
  adminPermissions: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  cellCount: InputMaybe<Scalars['Int']>;
  cellOntologyIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  description: InputMaybe<Scalars['String']>;
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  externalWebsite: InputMaybe<Scalars['String']>;
  filename: InputMaybe<Scalars['String']>;
  importFailed: InputMaybe<Scalars['Boolean']>;
  importFile: InputMaybe<Scalars['String']>;
  importFinished: InputMaybe<Scalars['Boolean']>;
  importLog: InputMaybe<Scalars['String']>;
  importStarted: InputMaybe<Scalars['Boolean']>;
  legacyConfig: InputMaybe<Scalars['JSON']>;
  metadata: InputMaybe<Scalars['JSON']>;
  organismTaxId: InputMaybe<Scalars['String']>;
  projections: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  readerPermissions: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  studyId: InputMaybe<Scalars['Int']>;
  studyName: InputMaybe<Scalars['String']>;
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  visible: InputMaybe<Scalars['Boolean']>;
};

export type StudySample = Node & {
  __typename?: 'StudySample';
  h5AdObsIndex: Scalars['Int'];
  h5AdObsKey: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads a single `Study` that is related to this `StudySample`. */
  study: Maybe<Study>;
  studyId: Scalars['Int'];
  studySampleId: Scalars['Int'];
  /** Reads and enables pagination through a set of `StudySampleProjection`. */
  studySampleProjectionsByStudyIdAndStudySampleId: StudySampleProjectionsConnection;
  /** Reads and enables pagination through a set of `StudySampleProjection`. */
  studySampleProjectionsByStudyIdAndStudySampleIdList: Array<StudySampleProjection>;
};


export type StudySampleStudySampleProjectionsByStudyIdAndStudySampleIdArgs = {
  after: InputMaybe<Scalars['Cursor']>;
  before: InputMaybe<Scalars['Cursor']>;
  condition: InputMaybe<StudySampleProjectionCondition>;
  filter: InputMaybe<StudySampleProjectionFilter>;
  first: InputMaybe<Scalars['Int']>;
  last: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySampleProjectionsOrderBy>>;
};


export type StudySampleStudySampleProjectionsByStudyIdAndStudySampleIdListArgs = {
  condition: InputMaybe<StudySampleProjectionCondition>;
  filter: InputMaybe<StudySampleProjectionFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySampleProjectionsOrderBy>>;
};

export type StudySampleAnnotation = {
  __typename?: 'StudySampleAnnotation';
  /** Reads a single `AnnotationValue` that is related to this `StudySampleAnnotation`. */
  annotationValue: Maybe<AnnotationValue>;
  annotationValueId: Scalars['Int'];
  color: Maybe<Scalars['String']>;
  /** Reads a single `Study` that is related to this `StudySampleAnnotation`. */
  study: Maybe<Study>;
  studyId: Scalars['Int'];
  studySampleIds: Array<Maybe<Scalars['Int']>>;
};

/**
 * A condition to be used against `StudySampleAnnotation` object types. All fields
 * are tested for equality and combined with a logical ‘and.’
 */
export type StudySampleAnnotationCondition = {
  /** Checks for equality with the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `color` field. */
  color: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleIds` field. */
  studySampleIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
};

/** A filter to be used against `StudySampleAnnotation` object types. All fields are combined with a logical ‘and.’ */
export type StudySampleAnnotationFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudySampleAnnotationFilter>>;
  /** Filter by the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<IntFilter>;
  /** Filter by the object’s `color` field. */
  color: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudySampleAnnotationFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudySampleAnnotationFilter>>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds: InputMaybe<IntListFilter>;
};

/** An input for mutations affecting `StudySampleAnnotation` */
export type StudySampleAnnotationInput = {
  annotationValueId: Scalars['Int'];
  color: InputMaybe<Scalars['String']>;
  studyId: Scalars['Int'];
  studySampleIds: Array<InputMaybe<Scalars['Int']>>;
};

export type StudySampleAnnotationSubsampling = {
  __typename?: 'StudySampleAnnotationSubsampling';
  annotationValueId: Maybe<Scalars['Int']>;
  rowNumber: Maybe<Scalars['BigInt']>;
  studyId: Maybe<Scalars['Int']>;
  studySampleIds: Maybe<Array<Maybe<Scalars['Int']>>>;
};

/**
 * A condition to be used against `StudySampleAnnotationSubsampling` object types.
 * All fields are tested for equality and combined with a logical ‘and.’
 */
export type StudySampleAnnotationSubsamplingCondition = {
  /** Checks for equality with the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `rowNumber` field. */
  rowNumber: InputMaybe<Scalars['BigInt']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleIds` field. */
  studySampleIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
};

/** A filter to be used against `StudySampleAnnotationSubsampling` object types. All fields are combined with a logical ‘and.’ */
export type StudySampleAnnotationSubsamplingFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudySampleAnnotationSubsamplingFilter>>;
  /** Filter by the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudySampleAnnotationSubsamplingFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudySampleAnnotationSubsamplingFilter>>;
  /** Filter by the object’s `rowNumber` field. */
  rowNumber: InputMaybe<BigIntFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds: InputMaybe<IntListFilter>;
};

/** A connection to a list of `StudySampleAnnotationSubsampling` values. */
export type StudySampleAnnotationSubsamplingsConnection = {
  __typename?: 'StudySampleAnnotationSubsamplingsConnection';
  /** A list of edges which contains the `StudySampleAnnotationSubsampling` and cursor to aid in pagination. */
  edges: Array<StudySampleAnnotationSubsamplingsEdge>;
  /** A list of `StudySampleAnnotationSubsampling` objects. */
  nodes: Array<Maybe<StudySampleAnnotationSubsampling>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudySampleAnnotationSubsampling` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudySampleAnnotationSubsampling` edge in the connection. */
export type StudySampleAnnotationSubsamplingsEdge = {
  __typename?: 'StudySampleAnnotationSubsamplingsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudySampleAnnotationSubsampling` at the end of the edge. */
  node: Maybe<StudySampleAnnotationSubsampling>;
};

/** Methods to use when ordering `StudySampleAnnotationSubsampling`. */
export enum StudySampleAnnotationSubsamplingsOrderBy {
  AnnotationValueIdAsc = 'ANNOTATION_VALUE_ID_ASC',
  AnnotationValueIdDesc = 'ANNOTATION_VALUE_ID_DESC',
  Natural = 'NATURAL',
  RowNumberAsc = 'ROW_NUMBER_ASC',
  RowNumberDesc = 'ROW_NUMBER_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC'
}

/** A connection to a list of `StudySampleAnnotation` values. */
export type StudySampleAnnotationsConnection = {
  __typename?: 'StudySampleAnnotationsConnection';
  /** A list of edges which contains the `StudySampleAnnotation` and cursor to aid in pagination. */
  edges: Array<StudySampleAnnotationsEdge>;
  /** A list of `StudySampleAnnotation` objects. */
  nodes: Array<Maybe<StudySampleAnnotation>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudySampleAnnotation` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudySampleAnnotation` edge in the connection. */
export type StudySampleAnnotationsEdge = {
  __typename?: 'StudySampleAnnotationsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudySampleAnnotation` at the end of the edge. */
  node: Maybe<StudySampleAnnotation>;
};

/** Methods to use when ordering `StudySampleAnnotation`. */
export enum StudySampleAnnotationsOrderBy {
  AnnotationValueIdAsc = 'ANNOTATION_VALUE_ID_ASC',
  AnnotationValueIdDesc = 'ANNOTATION_VALUE_ID_DESC',
  ColorAsc = 'COLOR_ASC',
  ColorDesc = 'COLOR_DESC',
  Natural = 'NATURAL',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC'
}

/**
 * A condition to be used against `StudySample` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type StudySampleCondition = {
  /** Checks for equality with the object’s `h5AdObsIndex` field. */
  h5AdObsIndex: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `h5AdObsKey` field. */
  h5AdObsKey: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleId` field. */
  studySampleId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudySample` object types. All fields are combined with a logical ‘and.’ */
export type StudySampleFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudySampleFilter>>;
  /** Filter by the object’s `h5AdObsIndex` field. */
  h5AdObsIndex: InputMaybe<IntFilter>;
  /** Filter by the object’s `h5AdObsKey` field. */
  h5AdObsKey: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudySampleFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudySampleFilter>>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleId` field. */
  studySampleId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `StudySample` */
export type StudySampleInput = {
  h5AdObsIndex: Scalars['Int'];
  h5AdObsKey: Scalars['String'];
  studyId: Scalars['Int'];
  studySampleId: Scalars['Int'];
};

/** Represents an update to a `StudySample`. Fields that are set will be updated. */
export type StudySamplePatch = {
  h5AdObsIndex: InputMaybe<Scalars['Int']>;
  h5AdObsKey: InputMaybe<Scalars['String']>;
  studyId: InputMaybe<Scalars['Int']>;
  studySampleId: InputMaybe<Scalars['Int']>;
};

export type StudySampleProjection = {
  __typename?: 'StudySampleProjection';
  displaySubsampling: Scalars['Boolean'];
  modality: Maybe<Scalars['String']>;
  projection: Array<Maybe<Scalars['Float']>>;
  projectionType: Scalars['String'];
  studyId: Scalars['Int'];
  studySampleId: Scalars['Int'];
  /** Reads a single `StudySample` that is related to this `StudySampleProjection`. */
  studyStudySample: Maybe<StudySample>;
};

/**
 * A condition to be used against `StudySampleProjection` object types. All fields
 * are tested for equality and combined with a logical ‘and.’
 */
export type StudySampleProjectionCondition = {
  /** Checks for equality with the object’s `displaySubsampling` field. */
  displaySubsampling: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `modality` field. */
  modality: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `projection` field. */
  projection: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Checks for equality with the object’s `projectionType` field. */
  projectionType: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleId` field. */
  studySampleId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudySampleProjection` object types. All fields are combined with a logical ‘and.’ */
export type StudySampleProjectionFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudySampleProjectionFilter>>;
  /** Filter by the object’s `displaySubsampling` field. */
  displaySubsampling: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `modality` field. */
  modality: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudySampleProjectionFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudySampleProjectionFilter>>;
  /** Filter by the object’s `projection` field. */
  projection: InputMaybe<FloatListFilter>;
  /** Filter by the object’s `projectionType` field. */
  projectionType: InputMaybe<StringFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleId` field. */
  studySampleId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `StudySampleProjection` */
export type StudySampleProjectionInput = {
  displaySubsampling: Scalars['Boolean'];
  modality: InputMaybe<Scalars['String']>;
  projection: Array<InputMaybe<Scalars['Float']>>;
  projectionType: Scalars['String'];
  studyId: Scalars['Int'];
  studySampleId: Scalars['Int'];
};

export type StudySampleProjectionSubsamplingTransposed = {
  __typename?: 'StudySampleProjectionSubsamplingTransposed';
  modality: Maybe<Scalars['String']>;
  projection: Maybe<Array<Maybe<Scalars['Float']>>>;
  projectionType: Maybe<Scalars['String']>;
  rowNumber: Maybe<Scalars['BigInt']>;
  studyId: Maybe<Scalars['Int']>;
  studySampleId: Maybe<Array<Maybe<Scalars['Int']>>>;
};

/**
 * A condition to be used against `StudySampleProjectionSubsamplingTransposed`
 * object types. All fields are tested for equality and combined with a logical ‘and.’
 */
export type StudySampleProjectionSubsamplingTransposedCondition = {
  /** Checks for equality with the object’s `modality` field. */
  modality: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `projection` field. */
  projection: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Checks for equality with the object’s `projectionType` field. */
  projectionType: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `rowNumber` field. */
  rowNumber: InputMaybe<Scalars['BigInt']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleId` field. */
  studySampleId: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
};

/** A filter to be used against `StudySampleProjectionSubsamplingTransposed` object types. All fields are combined with a logical ‘and.’ */
export type StudySampleProjectionSubsamplingTransposedFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudySampleProjectionSubsamplingTransposedFilter>>;
  /** Filter by the object’s `modality` field. */
  modality: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudySampleProjectionSubsamplingTransposedFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudySampleProjectionSubsamplingTransposedFilter>>;
  /** Filter by the object’s `projection` field. */
  projection: InputMaybe<FloatListFilter>;
  /** Filter by the object’s `projectionType` field. */
  projectionType: InputMaybe<StringFilter>;
  /** Filter by the object’s `rowNumber` field. */
  rowNumber: InputMaybe<BigIntFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleId` field. */
  studySampleId: InputMaybe<IntListFilter>;
};

/** A connection to a list of `StudySampleProjectionSubsamplingTransposed` values. */
export type StudySampleProjectionSubsamplingTransposedsConnection = {
  __typename?: 'StudySampleProjectionSubsamplingTransposedsConnection';
  /** A list of edges which contains the `StudySampleProjectionSubsamplingTransposed` and cursor to aid in pagination. */
  edges: Array<StudySampleProjectionSubsamplingTransposedsEdge>;
  /** A list of `StudySampleProjectionSubsamplingTransposed` objects. */
  nodes: Array<Maybe<StudySampleProjectionSubsamplingTransposed>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudySampleProjectionSubsamplingTransposed` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudySampleProjectionSubsamplingTransposed` edge in the connection. */
export type StudySampleProjectionSubsamplingTransposedsEdge = {
  __typename?: 'StudySampleProjectionSubsamplingTransposedsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudySampleProjectionSubsamplingTransposed` at the end of the edge. */
  node: Maybe<StudySampleProjectionSubsamplingTransposed>;
};

/** Methods to use when ordering `StudySampleProjectionSubsamplingTransposed`. */
export enum StudySampleProjectionSubsamplingTransposedsOrderBy {
  ModalityAsc = 'MODALITY_ASC',
  ModalityDesc = 'MODALITY_DESC',
  Natural = 'NATURAL',
  ProjectionAsc = 'PROJECTION_ASC',
  ProjectionDesc = 'PROJECTION_DESC',
  ProjectionTypeAsc = 'PROJECTION_TYPE_ASC',
  ProjectionTypeDesc = 'PROJECTION_TYPE_DESC',
  RowNumberAsc = 'ROW_NUMBER_ASC',
  RowNumberDesc = 'ROW_NUMBER_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudySampleIdAsc = 'STUDY_SAMPLE_ID_ASC',
  StudySampleIdDesc = 'STUDY_SAMPLE_ID_DESC'
}

/** A connection to a list of `StudySampleProjection` values. */
export type StudySampleProjectionsConnection = {
  __typename?: 'StudySampleProjectionsConnection';
  /** A list of edges which contains the `StudySampleProjection` and cursor to aid in pagination. */
  edges: Array<StudySampleProjectionsEdge>;
  /** A list of `StudySampleProjection` objects. */
  nodes: Array<Maybe<StudySampleProjection>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudySampleProjection` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudySampleProjection` edge in the connection. */
export type StudySampleProjectionsEdge = {
  __typename?: 'StudySampleProjectionsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudySampleProjection` at the end of the edge. */
  node: Maybe<StudySampleProjection>;
};

/** Methods to use when ordering `StudySampleProjection`. */
export enum StudySampleProjectionsOrderBy {
  DisplaySubsamplingAsc = 'DISPLAY_SUBSAMPLING_ASC',
  DisplaySubsamplingDesc = 'DISPLAY_SUBSAMPLING_DESC',
  ModalityAsc = 'MODALITY_ASC',
  ModalityDesc = 'MODALITY_DESC',
  Natural = 'NATURAL',
  ProjectionAsc = 'PROJECTION_ASC',
  ProjectionDesc = 'PROJECTION_DESC',
  ProjectionTypeAsc = 'PROJECTION_TYPE_ASC',
  ProjectionTypeDesc = 'PROJECTION_TYPE_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudySampleIdAsc = 'STUDY_SAMPLE_ID_ASC',
  StudySampleIdDesc = 'STUDY_SAMPLE_ID_DESC'
}

/** A connection to a list of `StudySample` values. */
export type StudySamplesConnection = {
  __typename?: 'StudySamplesConnection';
  /** A list of edges which contains the `StudySample` and cursor to aid in pagination. */
  edges: Array<StudySamplesEdge>;
  /** A list of `StudySample` objects. */
  nodes: Array<Maybe<StudySample>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudySample` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudySample` edge in the connection. */
export type StudySamplesEdge = {
  __typename?: 'StudySamplesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudySample` at the end of the edge. */
  node: Maybe<StudySample>;
};

/** Methods to use when ordering `StudySample`. */
export enum StudySamplesOrderBy {
  H5AdObsIndexAsc = 'H5AD_OBS_INDEX_ASC',
  H5AdObsIndexDesc = 'H5AD_OBS_INDEX_DESC',
  H5AdObsKeyAsc = 'H5AD_OBS_KEY_ASC',
  H5AdObsKeyDesc = 'H5AD_OBS_KEY_DESC',
  Natural = 'NATURAL',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudySampleIdAsc = 'STUDY_SAMPLE_ID_ASC',
  StudySampleIdDesc = 'STUDY_SAMPLE_ID_DESC'
}

export type StudyVisibleCurrentuser = {
  __typename?: 'StudyVisibleCurrentuser';
  studyId: Maybe<Scalars['Int']>;
};

/**
 * A condition to be used against `StudyVisibleCurrentuser` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type StudyVisibleCurrentuserCondition = {
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudyVisibleCurrentuser` object types. All fields are combined with a logical ‘and.’ */
export type StudyVisibleCurrentuserFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyVisibleCurrentuserFilter>>;
  /** Negates the expression. */
  not: InputMaybe<StudyVisibleCurrentuserFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyVisibleCurrentuserFilter>>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `StudyVisibleCurrentuser` */
export type StudyVisibleCurrentuserInput = {
  studyId: InputMaybe<Scalars['Int']>;
};

/** A connection to a list of `StudyVisibleCurrentuser` values. */
export type StudyVisibleCurrentusersConnection = {
  __typename?: 'StudyVisibleCurrentusersConnection';
  /** A list of edges which contains the `StudyVisibleCurrentuser` and cursor to aid in pagination. */
  edges: Array<StudyVisibleCurrentusersEdge>;
  /** A list of `StudyVisibleCurrentuser` objects. */
  nodes: Array<Maybe<StudyVisibleCurrentuser>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `StudyVisibleCurrentuser` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `StudyVisibleCurrentuser` edge in the connection. */
export type StudyVisibleCurrentusersEdge = {
  __typename?: 'StudyVisibleCurrentusersEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `StudyVisibleCurrentuser` at the end of the edge. */
  node: Maybe<StudyVisibleCurrentuser>;
};

/** Methods to use when ordering `StudyVisibleCurrentuser`. */
export enum StudyVisibleCurrentusersOrderBy {
  Natural = 'NATURAL',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

/** A connection to a list of `TreeOntology` values. */
export type TreeOntologiesConnection = {
  __typename?: 'TreeOntologiesConnection';
  /** A list of edges which contains the `TreeOntology` and cursor to aid in pagination. */
  edges: Array<TreeOntologiesEdge>;
  /** A list of `TreeOntology` objects. */
  nodes: Array<Maybe<TreeOntology>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `TreeOntology` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `TreeOntology` edge in the connection. */
export type TreeOntologiesEdge = {
  __typename?: 'TreeOntologiesEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `TreeOntology` at the end of the edge. */
  node: Maybe<TreeOntology>;
};

/** Methods to use when ordering `TreeOntology`. */
export enum TreeOntologiesOrderBy {
  CidAsc = 'CID_ASC',
  CidDesc = 'CID_DESC',
  LabelAsc = 'LABEL_ASC',
  LabelDesc = 'LABEL_DESC',
  Natural = 'NATURAL',
  OntologyAsc = 'ONTOLOGY_ASC',
  OntologyDesc = 'ONTOLOGY_DESC',
  OntCodeAsc = 'ONT_CODE_ASC',
  OntCodeDesc = 'ONT_CODE_DESC',
  ParentCidsAsc = 'PARENT_CIDS_ASC',
  ParentCidsDesc = 'PARENT_CIDS_DESC',
  ParentOntCodePathAsc = 'PARENT_ONT_CODE_PATH_ASC',
  ParentOntCodePathDesc = 'PARENT_ONT_CODE_PATH_DESC'
}

export type TreeOntology = {
  __typename?: 'TreeOntology';
  cid: Maybe<Scalars['Int']>;
  label: Maybe<Scalars['String']>;
  ontCode: Maybe<Scalars['String']>;
  ontology: Maybe<Scalars['String']>;
  parentCids: Maybe<Array<Maybe<Scalars['Int']>>>;
  parentOntCodePath: Maybe<Array<Maybe<Scalars['String']>>>;
};

/**
 * A condition to be used against `TreeOntology` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type TreeOntologyCondition = {
  /** Checks for equality with the object’s `cid` field. */
  cid: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `label` field. */
  label: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `ontCode` field. */
  ontCode: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `ontology` field. */
  ontology: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `parentCids` field. */
  parentCids: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `parentOntCodePath` field. */
  parentOntCodePath: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** A filter to be used against `TreeOntology` object types. All fields are combined with a logical ‘and.’ */
export type TreeOntologyFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<TreeOntologyFilter>>;
  /** Filter by the object’s `cid` field. */
  cid: InputMaybe<IntFilter>;
  /** Filter by the object’s `label` field. */
  label: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<TreeOntologyFilter>;
  /** Filter by the object’s `ontCode` field. */
  ontCode: InputMaybe<StringFilter>;
  /** Filter by the object’s `ontology` field. */
  ontology: InputMaybe<StringFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<TreeOntologyFilter>>;
  /** Filter by the object’s `parentCids` field. */
  parentCids: InputMaybe<IntListFilter>;
  /** Filter by the object’s `parentOntCodePath` field. */
  parentOntCodePath: InputMaybe<StringListFilter>;
};

/** All input for the `updateAnnotationGroupByNodeId` mutation. */
export type UpdateAnnotationGroupByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `AnnotationGroup` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `AnnotationGroup` being updated. */
  patch: AnnotationGroupPatch;
};

/** All input for the `updateAnnotationGroup` mutation. */
export type UpdateAnnotationGroupInput = {
  annotationGroupId: Scalars['Int'];
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** An object where the defined keys will be set on the `AnnotationGroup` being updated. */
  patch: AnnotationGroupPatch;
};

/** The output of our update `AnnotationGroup` mutation. */
export type UpdateAnnotationGroupPayload = {
  __typename?: 'UpdateAnnotationGroupPayload';
  /** The `AnnotationGroup` that was updated by this mutation. */
  annotationGroup: Maybe<AnnotationGroup>;
  /** An edge for our `AnnotationGroup`. May be used by Relay 1. */
  annotationGroupEdge: Maybe<AnnotationGroupsEdge>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our update `AnnotationGroup` mutation. */
export type UpdateAnnotationGroupPayloadAnnotationGroupEdgeArgs = {
  orderBy?: InputMaybe<Array<AnnotationGroupsOrderBy>>;
};

/** All input for the `updateAnnotationValueByNodeId` mutation. */
export type UpdateAnnotationValueByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `AnnotationValue` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `AnnotationValue` being updated. */
  patch: AnnotationValuePatch;
};

/** All input for the `updateAnnotationValue` mutation. */
export type UpdateAnnotationValueInput = {
  annotationValueId: Scalars['Int'];
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** An object where the defined keys will be set on the `AnnotationValue` being updated. */
  patch: AnnotationValuePatch;
};

/** The output of our update `AnnotationValue` mutation. */
export type UpdateAnnotationValuePayload = {
  __typename?: 'UpdateAnnotationValuePayload';
  /** Reads a single `AnnotationGroup` that is related to this `AnnotationValue`. */
  annotationGroup: Maybe<AnnotationGroup>;
  /** The `AnnotationValue` that was updated by this mutation. */
  annotationValue: Maybe<AnnotationValue>;
  /** An edge for our `AnnotationValue`. May be used by Relay 1. */
  annotationValueEdge: Maybe<AnnotationValuesEdge>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our update `AnnotationValue` mutation. */
export type UpdateAnnotationValuePayloadAnnotationValueEdgeArgs = {
  orderBy?: InputMaybe<Array<AnnotationValuesOrderBy>>;
};

/** All input for the `updateConceptByNodeId` mutation. */
export type UpdateConceptByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `Concept` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `Concept` being updated. */
  patch: ConceptPatch;
};

/** All input for the `updateConcept` mutation. */
export type UpdateConceptInput = {
  cid: Scalars['Int'];
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** An object where the defined keys will be set on the `Concept` being updated. */
  patch: ConceptPatch;
};

/** The output of our update `Concept` mutation. */
export type UpdateConceptPayload = {
  __typename?: 'UpdateConceptPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** The `Concept` that was updated by this mutation. */
  concept: Maybe<Concept>;
  /** An edge for our `Concept`. May be used by Relay 1. */
  conceptEdge: Maybe<ConceptsEdge>;
  /** Reads a single `Ontology` that is related to this `Concept`. */
  ontologyByOntid: Maybe<Ontology>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our update `Concept` mutation. */
export type UpdateConceptPayloadConceptEdgeArgs = {
  orderBy?: InputMaybe<Array<ConceptsOrderBy>>;
};

/** All input for the `updateOmicsBaseByNodeId` mutation. */
export type UpdateOmicsBaseByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `OmicsBase` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `OmicsBase` being updated. */
  patch: OmicsBasePatch;
};

/** All input for the `updateOmicsBase` mutation. */
export type UpdateOmicsBaseInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  omicsId: Scalars['Int'];
  /** An object where the defined keys will be set on the `OmicsBase` being updated. */
  patch: OmicsBasePatch;
};

/** The output of our update `OmicsBase` mutation. */
export type UpdateOmicsBasePayload = {
  __typename?: 'UpdateOmicsBasePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** The `OmicsBase` that was updated by this mutation. */
  omicsBase: Maybe<OmicsBase>;
  /** An edge for our `OmicsBase`. May be used by Relay 1. */
  omicsBaseEdge: Maybe<OmicsBasesEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our update `OmicsBase` mutation. */
export type UpdateOmicsBasePayloadOmicsBaseEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsBasesOrderBy>>;
};

/** All input for the `updateOmicsGeneByNodeId` mutation. */
export type UpdateOmicsGeneByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `OmicsGene` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `OmicsGene` being updated. */
  patch: OmicsGenePatch;
};

/** All input for the `updateOmicsGene` mutation. */
export type UpdateOmicsGeneInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  geneId: Scalars['Int'];
  /** An object where the defined keys will be set on the `OmicsGene` being updated. */
  patch: OmicsGenePatch;
};

/** The output of our update `OmicsGene` mutation. */
export type UpdateOmicsGenePayload = {
  __typename?: 'UpdateOmicsGenePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Reads a single `OmicsBase` that is related to this `OmicsGene`. */
  gene: Maybe<OmicsBase>;
  /** The `OmicsGene` that was updated by this mutation. */
  omicsGene: Maybe<OmicsGene>;
  /** An edge for our `OmicsGene`. May be used by Relay 1. */
  omicsGeneEdge: Maybe<OmicsGenesEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our update `OmicsGene` mutation. */
export type UpdateOmicsGenePayloadOmicsGeneEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsGenesOrderBy>>;
};

/** All input for the `updateOmicsProteinAntibodyTagByNodeId` mutation. */
export type UpdateOmicsProteinAntibodyTagByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `OmicsProteinAntibodyTag` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `OmicsProteinAntibodyTag` being updated. */
  patch: OmicsProteinAntibodyTagPatch;
};

/** All input for the `updateOmicsProteinAntibodyTag` mutation. */
export type UpdateOmicsProteinAntibodyTagInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** An object where the defined keys will be set on the `OmicsProteinAntibodyTag` being updated. */
  patch: OmicsProteinAntibodyTagPatch;
  proteinAntibodyTagId: Scalars['Int'];
};

/** The output of our update `OmicsProteinAntibodyTag` mutation. */
export type UpdateOmicsProteinAntibodyTagPayload = {
  __typename?: 'UpdateOmicsProteinAntibodyTagPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** The `OmicsProteinAntibodyTag` that was updated by this mutation. */
  omicsProteinAntibodyTag: Maybe<OmicsProteinAntibodyTag>;
  /** An edge for our `OmicsProteinAntibodyTag`. May be used by Relay 1. */
  omicsProteinAntibodyTagEdge: Maybe<OmicsProteinAntibodyTagsEdge>;
  /** Reads a single `OmicsBase` that is related to this `OmicsProteinAntibodyTag`. */
  proteinAntibodyTag: Maybe<OmicsBase>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our update `OmicsProteinAntibodyTag` mutation. */
export type UpdateOmicsProteinAntibodyTagPayloadOmicsProteinAntibodyTagEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsProteinAntibodyTagsOrderBy>>;
};

/** All input for the `updateOmicsRegionByNodeId` mutation. */
export type UpdateOmicsRegionByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `OmicsRegion` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `OmicsRegion` being updated. */
  patch: OmicsRegionPatch;
};

/** All input for the `updateOmicsRegion` mutation. */
export type UpdateOmicsRegionInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** An object where the defined keys will be set on the `OmicsRegion` being updated. */
  patch: OmicsRegionPatch;
  regionId: Scalars['Int'];
};

/** The output of our update `OmicsRegion` mutation. */
export type UpdateOmicsRegionPayload = {
  __typename?: 'UpdateOmicsRegionPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** The `OmicsRegion` that was updated by this mutation. */
  omicsRegion: Maybe<OmicsRegion>;
  /** An edge for our `OmicsRegion`. May be used by Relay 1. */
  omicsRegionEdge: Maybe<OmicsRegionsEdge>;
  /** Reads a single `OmicsBase` that is related to this `OmicsRegion`. */
  omics_region_newNameHere: Maybe<OmicsBase>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our update `OmicsRegion` mutation. */
export type UpdateOmicsRegionPayloadOmicsRegionEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsRegionsOrderBy>>;
};

/** All input for the `updateOmicsTranscriptionFactorByNodeId` mutation. */
export type UpdateOmicsTranscriptionFactorByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `OmicsTranscriptionFactor` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `OmicsTranscriptionFactor` being updated. */
  patch: OmicsTranscriptionFactorPatch;
};

/** All input for the `updateOmicsTranscriptionFactor` mutation. */
export type UpdateOmicsTranscriptionFactorInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  omicsId: Scalars['Int'];
  /** An object where the defined keys will be set on the `OmicsTranscriptionFactor` being updated. */
  patch: OmicsTranscriptionFactorPatch;
};

/** The output of our update `OmicsTranscriptionFactor` mutation. */
export type UpdateOmicsTranscriptionFactorPayload = {
  __typename?: 'UpdateOmicsTranscriptionFactorPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Reads a single `OmicsBase` that is related to this `OmicsTranscriptionFactor`. */
  omics: Maybe<OmicsBase>;
  /** The `OmicsTranscriptionFactor` that was updated by this mutation. */
  omicsTranscriptionFactor: Maybe<OmicsTranscriptionFactor>;
  /** An edge for our `OmicsTranscriptionFactor`. May be used by Relay 1. */
  omicsTranscriptionFactorEdge: Maybe<OmicsTranscriptionFactorsEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our update `OmicsTranscriptionFactor` mutation. */
export type UpdateOmicsTranscriptionFactorPayloadOmicsTranscriptionFactorEdgeArgs = {
  orderBy?: InputMaybe<Array<OmicsTranscriptionFactorsOrderBy>>;
};

/** All input for the `updateOntologyByNodeId` mutation. */
export type UpdateOntologyByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `Ontology` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `Ontology` being updated. */
  patch: OntologyPatch;
};

/** All input for the `updateOntology` mutation. */
export type UpdateOntologyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  ontid: Scalars['Int'];
  /** An object where the defined keys will be set on the `Ontology` being updated. */
  patch: OntologyPatch;
};

/** The output of our update `Ontology` mutation. */
export type UpdateOntologyPayload = {
  __typename?: 'UpdateOntologyPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** The `Ontology` that was updated by this mutation. */
  ontology: Maybe<Ontology>;
  /** An edge for our `Ontology`. May be used by Relay 1. */
  ontologyEdge: Maybe<OntologiesEdge>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};


/** The output of our update `Ontology` mutation. */
export type UpdateOntologyPayloadOntologyEdgeArgs = {
  orderBy?: InputMaybe<Array<OntologiesOrderBy>>;
};

/** All input for the `updateStudyByNodeId` mutation. */
export type UpdateStudyByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `Study` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `Study` being updated. */
  patch: StudyPatch;
};

/** All input for the `updateStudy` mutation. */
export type UpdateStudyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** An object where the defined keys will be set on the `Study` being updated. */
  patch: StudyPatch;
  studyId: Scalars['Int'];
};

/** All input for the `updateStudyLayerByNodeId` mutation. */
export type UpdateStudyLayerByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `StudyLayer` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `StudyLayer` being updated. */
  patch: StudyLayerPatch;
};

/** All input for the `updateStudyLayer` mutation. */
export type UpdateStudyLayerInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** An object where the defined keys will be set on the `StudyLayer` being updated. */
  patch: StudyLayerPatch;
  studyLayerId: Scalars['Int'];
};

/** The output of our update `StudyLayer` mutation. */
export type UpdateStudyLayerPayload = {
  __typename?: 'UpdateStudyLayerPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudyLayer`. */
  study: Maybe<Study>;
  /** The `StudyLayer` that was updated by this mutation. */
  studyLayer: Maybe<StudyLayer>;
  /** An edge for our `StudyLayer`. May be used by Relay 1. */
  studyLayerEdge: Maybe<StudyLayersEdge>;
};


/** The output of our update `StudyLayer` mutation. */
export type UpdateStudyLayerPayloadStudyLayerEdgeArgs = {
  orderBy?: InputMaybe<Array<StudyLayersOrderBy>>;
};

/** The output of our update `Study` mutation. */
export type UpdateStudyPayload = {
  __typename?: 'UpdateStudyPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** The `Study` that was updated by this mutation. */
  study: Maybe<Study>;
  /** An edge for our `Study`. May be used by Relay 1. */
  studyEdge: Maybe<StudiesEdge>;
};


/** The output of our update `Study` mutation. */
export type UpdateStudyPayloadStudyEdgeArgs = {
  orderBy?: InputMaybe<Array<StudiesOrderBy>>;
};

/** All input for the `updateStudySampleByNodeId` mutation. */
export type UpdateStudySampleByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `StudySample` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `StudySample` being updated. */
  patch: StudySamplePatch;
};

/** All input for the `updateStudySample` mutation. */
export type UpdateStudySampleInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  /** An object where the defined keys will be set on the `StudySample` being updated. */
  patch: StudySamplePatch;
  studyId: Scalars['Int'];
  studySampleId: Scalars['Int'];
};

/** The output of our update `StudySample` mutation. */
export type UpdateStudySamplePayload = {
  __typename?: 'UpdateStudySamplePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudySample`. */
  study: Maybe<Study>;
  /** The `StudySample` that was updated by this mutation. */
  studySample: Maybe<StudySample>;
  /** An edge for our `StudySample`. May be used by Relay 1. */
  studySampleEdge: Maybe<StudySamplesEdge>;
};


/** The output of our update `StudySample` mutation. */
export type UpdateStudySamplePayloadStudySampleEdgeArgs = {
  orderBy?: InputMaybe<Array<StudySamplesOrderBy>>;
};

/** All input for the `userAnnotationDefine` mutation. */
export type UserAnnotationDefineInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  pAnnotationGroupName: InputMaybe<Scalars['String']>;
  pSelectedSampleIds: InputMaybe<Scalars['String']>;
  pStudyId: InputMaybe<Scalars['Int']>;
  pUnexpressedSamplesOmicsIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
};

/** The output of our `userAnnotationDefine` mutation. */
export type UserAnnotationDefinePayload = {
  __typename?: 'UserAnnotationDefinePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  integer: Maybe<Scalars['Int']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};

/** All input for the `userAnnotationDelete` mutation. */
export type UserAnnotationDeleteInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  pAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pStudyId: InputMaybe<Scalars['Int']>;
};

/** The output of our `userAnnotationDelete` mutation. */
export type UserAnnotationDeletePayload = {
  __typename?: 'UserAnnotationDeletePayload';
  boolean: Maybe<Scalars['Boolean']>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};

/** All input for the `userAnnotationEdit` mutation. */
export type UserAnnotationEditInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId: InputMaybe<Scalars['String']>;
  pAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pPrivateToUser: InputMaybe<Scalars['Boolean']>;
  pStudyId: InputMaybe<Scalars['Int']>;
};

/** The output of our `userAnnotationEdit` mutation. */
export type UserAnnotationEditPayload = {
  __typename?: 'UserAnnotationEditPayload';
  boolean: Maybe<Scalars['Boolean']>;
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
};

export type UserAnnotationGroup = {
  __typename?: 'UserAnnotationGroup';
  calculateDifferentialExpression: Scalars['Boolean'];
  createdByUser: Maybe<Scalars['String']>;
  creationTimestamp: Scalars['Datetime'];
  privateToUser: Scalars['Boolean'];
  /** Reads a single `AnnotationGroup` that is related to this `UserAnnotationGroup`. */
  savedAsAnnotationGroup: Maybe<AnnotationGroup>;
  savedAsAnnotationGroupId: Scalars['Int'];
  /** Reads a single `Study` that is related to this `UserAnnotationGroup`. */
  study: Maybe<Study>;
  studyId: Scalars['Int'];
};

/**
 * A condition to be used against `UserAnnotationGroup` object types. All fields
 * are tested for equality and combined with a logical ‘and.’
 */
export type UserAnnotationGroupCondition = {
  /** Checks for equality with the object’s `calculateDifferentialExpression` field. */
  calculateDifferentialExpression: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `createdByUser` field. */
  createdByUser: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `creationTimestamp` field. */
  creationTimestamp: InputMaybe<Scalars['Datetime']>;
  /** Checks for equality with the object’s `privateToUser` field. */
  privateToUser: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `savedAsAnnotationGroupId` field. */
  savedAsAnnotationGroupId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `UserAnnotationGroup` object types. All fields are combined with a logical ‘and.’ */
export type UserAnnotationGroupFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<UserAnnotationGroupFilter>>;
  /** Filter by the object’s `calculateDifferentialExpression` field. */
  calculateDifferentialExpression: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `createdByUser` field. */
  createdByUser: InputMaybe<StringFilter>;
  /** Filter by the object’s `creationTimestamp` field. */
  creationTimestamp: InputMaybe<DatetimeFilter>;
  /** Negates the expression. */
  not: InputMaybe<UserAnnotationGroupFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<UserAnnotationGroupFilter>>;
  /** Filter by the object’s `privateToUser` field. */
  privateToUser: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `savedAsAnnotationGroupId` field. */
  savedAsAnnotationGroupId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `UserAnnotationGroup` */
export type UserAnnotationGroupInput = {
  calculateDifferentialExpression: Scalars['Boolean'];
  createdByUser: InputMaybe<Scalars['String']>;
  creationTimestamp: InputMaybe<Scalars['Datetime']>;
  privateToUser: InputMaybe<Scalars['Boolean']>;
  savedAsAnnotationGroupId: Scalars['Int'];
  studyId: Scalars['Int'];
};

/** A connection to a list of `UserAnnotationGroup` values. */
export type UserAnnotationGroupsConnection = {
  __typename?: 'UserAnnotationGroupsConnection';
  /** A list of edges which contains the `UserAnnotationGroup` and cursor to aid in pagination. */
  edges: Array<UserAnnotationGroupsEdge>;
  /** A list of `UserAnnotationGroup` objects. */
  nodes: Array<Maybe<UserAnnotationGroup>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `UserAnnotationGroup` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `UserAnnotationGroup` edge in the connection. */
export type UserAnnotationGroupsEdge = {
  __typename?: 'UserAnnotationGroupsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `UserAnnotationGroup` at the end of the edge. */
  node: Maybe<UserAnnotationGroup>;
};

/** Methods to use when ordering `UserAnnotationGroup`. */
export enum UserAnnotationGroupsOrderBy {
  CalculateDifferentialExpressionAsc = 'CALCULATE_DIFFERENTIAL_EXPRESSION_ASC',
  CalculateDifferentialExpressionDesc = 'CALCULATE_DIFFERENTIAL_EXPRESSION_DESC',
  CreatedByUserAsc = 'CREATED_BY_USER_ASC',
  CreatedByUserDesc = 'CREATED_BY_USER_DESC',
  CreationTimestampAsc = 'CREATION_TIMESTAMP_ASC',
  CreationTimestampDesc = 'CREATION_TIMESTAMP_DESC',
  Natural = 'NATURAL',
  PrivateToUserAsc = 'PRIVATE_TO_USER_ASC',
  PrivateToUserDesc = 'PRIVATE_TO_USER_DESC',
  SavedAsAnnotationGroupIdAsc = 'SAVED_AS_ANNOTATION_GROUP_ID_ASC',
  SavedAsAnnotationGroupIdDesc = 'SAVED_AS_ANNOTATION_GROUP_ID_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

export type _AllUsedOntologyId = {
  __typename?: '_AllUsedOntologyId';
  ontCode: Maybe<Scalars['String']>;
  ontology: Maybe<Scalars['String']>;
};

/**
 * A condition to be used against `_AllUsedOntologyId` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type _AllUsedOntologyIdCondition = {
  /** Checks for equality with the object’s `ontCode` field. */
  ontCode: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `ontology` field. */
  ontology: InputMaybe<Scalars['String']>;
};

/** A filter to be used against `_AllUsedOntologyId` object types. All fields are combined with a logical ‘and.’ */
export type _AllUsedOntologyIdFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<_AllUsedOntologyIdFilter>>;
  /** Negates the expression. */
  not: InputMaybe<_AllUsedOntologyIdFilter>;
  /** Filter by the object’s `ontCode` field. */
  ontCode: InputMaybe<StringFilter>;
  /** Filter by the object’s `ontology` field. */
  ontology: InputMaybe<StringFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<_AllUsedOntologyIdFilter>>;
};

/** A connection to a list of `_AllUsedOntologyId` values. */
export type _AllUsedOntologyIdsConnection = {
  __typename?: '_AllUsedOntologyIdsConnection';
  /** A list of edges which contains the `_AllUsedOntologyId` and cursor to aid in pagination. */
  edges: Array<_AllUsedOntologyIdsEdge>;
  /** A list of `_AllUsedOntologyId` objects. */
  nodes: Array<Maybe<_AllUsedOntologyId>>;
  /** Information to aid in pagination. */
  pageInfo: PageInfo;
  /** The count of *all* `_AllUsedOntologyId` you could get from the connection. */
  totalCount: Scalars['Int'];
};

/** A `_AllUsedOntologyId` edge in the connection. */
export type _AllUsedOntologyIdsEdge = {
  __typename?: '_AllUsedOntologyIdsEdge';
  /** A cursor for use in pagination. */
  cursor: Maybe<Scalars['Cursor']>;
  /** The `_AllUsedOntologyId` at the end of the edge. */
  node: Maybe<_AllUsedOntologyId>;
};

/** Methods to use when ordering `_AllUsedOntologyId`. */
export enum _AllUsedOntologyIdsOrderBy {
  Natural = 'NATURAL',
  OntologyAsc = 'ONTOLOGY_ASC',
  OntologyDesc = 'ONTOLOGY_DESC',
  OntCodeAsc = 'ONT_CODE_ASC',
  OntCodeDesc = 'ONT_CODE_DESC'
}

export type CorrelatedgenesQueryVariables = Exact<{
  studyId: Scalars['Int'];
  omicsId: Scalars['Int'];
}>;


export type CorrelatedgenesQuery = { __typename?: 'Query', getCorrelatedGenesList: Array<{ __typename?: 'GetCorrelatedGenesRecord', displayName: string, displaySymbol: string, omicsId: number, r: number }> };

export type StudyInfoFragment = { __typename?: 'StudyOverview', studyId: number, studyName: string, description: string, cellCount: number, externalWebsite: string, metadata: any, defaultStudyLayerId: number, studyOntologyList: Array<{ __typename?: 'StudyOverviewOntology', ontCodes: Array<string>, labels: Array<string>, ontology: string, parentIds: Array<string> }> };

export type ReferenceStudyInfoFragment = { __typename?: 'ReferenceStudyOverview', organismTaxId: string, studyName: string, studyId: number, cellCount: number, defaultStudyLayerId: number, externalWebsite: string, description: string, metadata: any, studyOntologyList: Array<{ __typename?: 'StudyOverviewOntology', ontCodes: Array<string>, ontology: string, labels: Array<string>, parentIds: Array<string> }>, referenceStudyInfoList: Array<{ __typename?: 'ReferenceStudy', celltypeAnnotationGroupId: number, tissueAnnotationGroupId: number }> };

export type TreeOntologyOverviewFragment = { __typename?: 'TreeOntology', label: string, ontCode: string, ontology: string, parentOntCodePath: Array<string> };

export type DegQueryVariables = Exact<{
  studyId: Scalars['Int'];
  annotationValueId: Scalars['Int'];
}>;


export type DegQuery = { __typename?: 'Query', differentialExpressionVsList: Array<{ __typename?: 'DifferentialExpressionV', omicsId: number, studyId: number, annotationValueId: number, omicsType: OmicsType, displayName: string, displaySymbol: string, pvalueAdj: number, log2Foldchange: number, linkedGenes: Array<number> }> };

export type StudiesQueryVariables = Exact<{ [key: string]: never; }>;


export type StudiesQuery = { __typename?: 'Query', studyOverviewsList: Array<{ __typename?: 'StudyOverview', studyId: number, studyName: string, description: string, cellCount: number, externalWebsite: string, metadata: any, defaultStudyLayerId: number, studyOntologyList: Array<{ __typename?: 'StudyOverviewOntology', ontCodes: Array<string>, labels: Array<string>, ontology: string, parentIds: Array<string> }> }> };

export type StudiesTreeOntologiesQueryVariables = Exact<{ [key: string]: never; }>;


export type StudiesTreeOntologiesQuery = { __typename?: 'Query', treeOntologiesList: Array<{ __typename?: 'TreeOntology', label: string, ontCode: string, ontology: string, parentOntCodePath: Array<string> }> };

export type SingleStudyInfoQueryVariables = Exact<{
  studyId: Scalars['Int'];
}>;


export type SingleStudyInfoQuery = { __typename?: 'Query', studyOverviewsList: Array<{ __typename?: 'StudyOverview', studyId: number, studyName: string, description: string, cellCount: number, externalWebsite: string, metadata: any, defaultStudyLayerId: number, studyOntologyList: Array<{ __typename?: 'StudyOverviewOntology', ontCodes: Array<string>, labels: Array<string>, ontology: string, parentIds: Array<string> }> }> };

export type ReferenceStudiesOverviewQueryVariables = Exact<{
  organismTaxId: Scalars['String'];
}>;


export type ReferenceStudiesOverviewQuery = { __typename?: 'Query', referenceStudyOverviewsList: Array<{ __typename?: 'ReferenceStudyOverview', organismTaxId: string, studyName: string, studyId: number, cellCount: number, defaultStudyLayerId: number, externalWebsite: string, description: string, metadata: any, studyOntologyList: Array<{ __typename?: 'StudyOverviewOntology', ontCodes: Array<string>, ontology: string, labels: Array<string>, parentIds: Array<string> }>, referenceStudyInfoList: Array<{ __typename?: 'ReferenceStudy', celltypeAnnotationGroupId: number, tissueAnnotationGroupId: number }> }> };

export type GeneSpecificityStudyQueryVariables = Exact<{
  studyId: Scalars['Int'];
}>;


export type GeneSpecificityStudyQuery = { __typename?: 'Query', referenceStudyOverviewsList: Array<{ __typename?: 'ReferenceStudyOverview', organismTaxId: string, studyName: string, studyId: number, cellCount: number, defaultStudyLayerId: number, externalWebsite: string, description: string, metadata: any, studyOntologyList: Array<{ __typename?: 'StudyOverviewOntology', ontCodes: Array<string>, ontology: string, labels: Array<string>, parentIds: Array<string> }>, referenceStudyInfoList: Array<{ __typename?: 'ReferenceStudy', celltypeAnnotationGroupId: number, tissueAnnotationGroupId: number }> }>, studyAnnotationFrontendGroupsList: Array<{ __typename?: 'StudyAnnotationFrontendGroup', annotationGroupId: number, annotationValuesList: Array<{ __typename?: 'StudyAnnotationFrontendValue', annotationValueId: number, displayValue: string, color: string, sampleCount: number }> }> };

export type AnnotationGrpFragment = { __typename?: 'StudyAnnotationFrontendGroup', annotationGroupId: number, isPrimary: boolean, ordering: number, displayGroup: string, differentialExpressionCalculated: boolean, createdByUser: string, currentUserIsOwner: boolean, privateToUser: boolean, annotationValuesList: Array<{ __typename?: 'StudyAnnotationFrontendValue', annotationValueId: number, displayValue: string, color: string, sampleCount: number }> };

export type DifferentialMarkerFragment = { __typename?: 'DifferentialExpression', annotationValueId: number, log2Foldchange: number, pvalueAdj: number, score: number, study: { __typename?: 'Study', studyName: string, studyId: number, cellCount: number, description: string }, annotationValue: { __typename?: 'AnnotationValue', displayValue: string, annotationGroup: { __typename?: 'AnnotationGroup', displayGroup: string, annotationGroupId: number } }, omics: { __typename?: 'OmicsBase', displaySymbol: string, taxId: number, omicsId: number, omicsType: OmicsType, displayName: string } };

export type StudiesWithMarkerGenesQueryVariables = Exact<{
  omicsIds: Array<Scalars['Int']> | Scalars['Int'];
}>;


export type StudiesWithMarkerGenesQuery = { __typename?: 'Query', differentialExpressionsList: Array<{ __typename?: 'DifferentialExpression', annotationValueId: number, log2Foldchange: number, pvalueAdj: number, score: number, study: { __typename?: 'Study', studyName: string, studyId: number, cellCount: number, description: string }, annotationValue: { __typename?: 'AnnotationValue', displayValue: string, annotationGroup: { __typename?: 'AnnotationGroup', displayGroup: string, annotationGroupId: number } }, omics: { __typename?: 'OmicsBase', displaySymbol: string, taxId: number, omicsId: number, omicsType: OmicsType, displayName: string } }> };

export type OmicsGeneFragment = { __typename?: 'OmicsBase', displayName: string, displaySymbol: string, omicsId: number, taxId: number };

export type AllGenesQueryVariables = Exact<{ [key: string]: never; }>;


export type AllGenesQuery = { __typename?: 'Query', omicsBasesList: Array<{ __typename?: 'OmicsBase', displayName: string, displaySymbol: string, omicsId: number, taxId: number, value: string, ontology: OmicsType }> };

export type StudyOmicsQueryVariables = Exact<{
  studyId: Scalars['Int'];
}>;


export type StudyOmicsQuery = { __typename?: 'Query', studyOmicsList: Array<{ __typename?: 'StudyOmic', omics: { __typename?: 'OmicsBase', omicsId: number, displayName: string, displaySymbol: string } }> };

export type StudyBasicsFragment = { __typename?: 'Study', studyId: number, studyName: string, organismTaxId: string, cellCount: number, projections: Array<string>, studyLayersList: Array<{ __typename?: 'StudyLayer', layer: string, studyLayerId: number }> };

export type StudyBasicsQueryVariables = Exact<{
  studyId: Scalars['Int'];
}>;


export type StudyBasicsQuery = { __typename?: 'Query', study: { __typename?: 'Study', studyId: number, studyName: string, cellCount: number, projections: Array<string>, organismTaxId: string }, studyLayersList: Array<{ __typename?: 'StudyLayer', layer: string, studyLayerId: number }> };

export type StudyBasics2QueryVariables = Exact<{
  studyId: Scalars['Int'];
}>;


export type StudyBasics2Query = { __typename?: 'Query', studyAnnotationFrontendGroupsList: Array<{ __typename?: 'StudyAnnotationFrontendGroup', annotationGroupId: number, isPrimary: boolean, ordering: number, displayGroup: string, differentialExpressionCalculated: boolean, createdByUser: string, currentUserIsOwner: boolean, privateToUser: boolean, annotationValuesList: Array<{ __typename?: 'StudyAnnotationFrontendValue', annotationValueId: number, displayValue: string, color: string, sampleCount: number }> }>, studySampleProjectionSubsamplingTransposedsList: Array<{ __typename?: 'StudySampleProjectionSubsamplingTransposed', projectionType: string, studySampleId: Array<number>, projection: Array<number>, modality: string }> };

export type StudyBasics3QueryVariables = Exact<{
  studyId: Scalars['Int'];
}>;


export type StudyBasics3Query = { __typename?: 'Query', studySampleAnnotationSubsamplingsList: Array<{ __typename?: 'StudySampleAnnotationSubsampling', annotationValueId: number, studySampleIds: Array<number> }>, studyOmicsTransposedsList: Array<{ __typename?: 'StudyOmicsTransposed', displayName: Array<string>, displaySymbol: Array<string>, omicsId: Array<number>, omicsType: Array<OmicsType> }> };

export type ExpressionByOmicsIdsQueryVariables = Exact<{
  studyLayerId: Scalars['Int'];
  omicsIds: Array<Scalars['Int']> | Scalars['Int'];
  subsamplingProjection: InputMaybe<Scalars['String']>;
}>;


export type ExpressionByOmicsIdsQuery = { __typename?: 'Query', expressionByOmicsIdsList: Array<{ __typename?: 'ExpressionByOmic', omicsId: number, studySampleIds: Array<number>, values: Array<number> }> };

export type ExpressionViolinPlotQueryVariables = Exact<{
  studyId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  omicsId: Scalars['Int'];
  annotationGroupId: Scalars['Int'];
  excludeAnnotationValueIds: Array<Scalars['Int']> | Scalars['Int'];
  annotationSecondaryGroupId: InputMaybe<Scalars['Int']>;
}>;


export type ExpressionViolinPlotQuery = { __typename?: 'Query', violinPlot: string };

export type ExpressionGroupTableQueryVariables = Exact<{
  studyId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  omicsId: Scalars['Int'];
  annotationGroupId: Scalars['Int'];
  excludeAnnotationValueIds: Array<Scalars['Int']> | Scalars['Int'];
  annotationSecondaryGroupId: Scalars['Int'];
}>;


export type ExpressionGroupTableQuery = { __typename?: 'Query', expressionByTwoAnnotationsList: Array<{ __typename?: 'ExpressionByTwoAnnotation', annotationValueId: number, annotationDisplayValue: string, secondAnnotationValueId: number, secondAnnotationDisplayValue: string, valueCount: number, median: number, nonZeroValueCount: number }> };

export type ExpressionTTestQueryVariables = Exact<{
  studyId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  omicsId: Scalars['Int'];
  annotationGroupId: Scalars['Int'];
  excludeAnnotationValueIds: Array<Scalars['Int']> | Scalars['Int'];
  secondAnnotationGroupId: InputMaybe<Scalars['Int']>;
  sample1AnnotationValueId: Scalars['Int'];
  sample1SecondAnnotationValueId: InputMaybe<Scalars['Int']>;
  sample2AnnotationValueId: Scalars['Int'];
  sample2SecondAnnotationValueId: InputMaybe<Scalars['Int']>;
}>;


export type ExpressionTTestQuery = { __typename?: 'Query', expressionTtest: string };

export type ExpressionCorrelationTrianglePlotQueryVariables = Exact<{
  studyId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  omicsIds: Array<Scalars['Int']> | Scalars['Int'];
  excludeAnnotationValueIds: Array<Scalars['Int']> | Scalars['Int'];
}>;


export type ExpressionCorrelationTrianglePlotQuery = { __typename?: 'Query', correlationTrianglePlot: string };

export type AutocompleteQueryVariables = Exact<{
  query: Scalars['String'];
}>;


export type AutocompleteQuery = { __typename?: 'Query', autocompleteList: Array<{ __typename?: 'AutocompleteResult', isSynonymOfPreferredTerm: string, label: string, labelHighlight: string, ontCode: string, ontology: string }> };

export type OmicsAutocompleteQueryVariables = Exact<{
  searchQuery: Scalars['String'];
  omicsType: OmicsType;
  taxId: Scalars['Int'];
}>;


export type OmicsAutocompleteQuery = { __typename?: 'Query', omicsAutocompleteList: Array<{ __typename?: 'OmicsAutocompleteResult', displaySymbol: string, labelHighlight: string, omicsId: Array<number>, omicsType: Array<OmicsType> }> };

export type OntologyOverviewFragment = { __typename?: 'Ontology', name: string, ontid: number, nodeId: string };

export type OntologiesQueryVariables = Exact<{ [key: string]: never; }>;


export type OntologiesQuery = { __typename?: 'Query', ontologiesList: Array<{ __typename?: 'Ontology', name: string, ontid: number, nodeId: string }> };

export type DotPlotElementFragment = { __typename?: 'ExpressionByAnnotation', studyLayerId: number, omicsId: number, annotationValueId: number, annotationDisplayValue: string, q3: number, median: number, exprSamplesFraction: number };

export type ExpressionByAnnotationQueryVariables = Exact<{
  studyLayerIds: Array<Scalars['Int']> | Scalars['Int'];
  omicsIds: Array<Scalars['Int']> | Scalars['Int'];
  annotationGroupId: Scalars['Int'];
  excludeAnnotationValueIds: Array<Scalars['Int']> | Scalars['Int'];
}>;


export type ExpressionByAnnotationQuery = { __typename?: 'Query', expressionByAnnotationList: Array<{ __typename?: 'ExpressionByAnnotation', studyLayerId: number, omicsId: number, annotationValueId: number, annotationDisplayValue: string, q3: number, median: number, exprSamplesFraction: number }> };

export type CellOAnnotationGroupIdQueryVariables = Exact<{ [key: string]: never; }>;


export type CellOAnnotationGroupIdQuery = { __typename?: 'Query', annotationGroupsList: Array<{ __typename?: 'AnnotationGroup', annotationGroupId: number }> };

export type HalfAVolcanoQueryVariables = Exact<{
  annotationValueId: Scalars['Int'];
  studyId: Scalars['Int'];
}>;


export type HalfAVolcanoQuery = { __typename?: 'Query', differentialExpressionsList: Array<{ __typename?: 'DifferentialExpression', log2Foldchange: number, pvalueAdj: number }> };

export type AnnotationValueCoocurrenceQueryVariables = Exact<{
  studyId: Scalars['Int'];
  annotationGroupId1: Scalars['Int'];
  annotationGroupId2: Scalars['Int'];
}>;


export type AnnotationValueCoocurrenceQuery = { __typename?: 'Query', annotationValueCoocurrenceList: Array<{ __typename?: 'AnnotationValueCoocurrenceRecord', annotationValueId1: number, annotationValueId2: number, occurrence: number }> };

export type SaveUserAnnotationMutationVariables = Exact<{
  studyId: Scalars['Int'];
  annotationGroupName: Scalars['String'];
  selectedSampleIds: Scalars['String'];
  unexpressedSamplesOmicsIds: InputMaybe<Array<Scalars['Int']> | Scalars['Int']>;
}>;


export type SaveUserAnnotationMutation = { __typename?: 'Mutation', userAnnotationDefine: { __typename?: 'UserAnnotationDefinePayload', clientMutationId: string, integer: number } };

export type EditUserAnnotationMutationVariables = Exact<{
  studyId: Scalars['Int'];
  annotationGroupId: Scalars['Int'];
  privateToUser: Scalars['Boolean'];
}>;


export type EditUserAnnotationMutation = { __typename?: 'Mutation', userAnnotationEdit: { __typename?: 'UserAnnotationEditPayload', boolean: boolean } };

export type DeleteUserAnnotationMutationVariables = Exact<{
  studyId: Scalars['Int'];
  annotationGroupId: Scalars['Int'];
}>;


export type DeleteUserAnnotationMutation = { __typename?: 'Mutation', userAnnotationDelete: { __typename?: 'UserAnnotationDeletePayload', boolean: boolean } };

export type StudyAdminDetailsFragment = { __typename?: 'StudyAdminDetail', studyId: number, studyName: string, description: string, filename: string, cellCount: number, tissueNcitIds: Array<string>, diseaseMeshIds: Array<string>, visible: boolean, externalWebsite: string, readerPermissions: Array<string>, readerPermissionGranted: boolean, adminPermissions: Array<string>, adminPermissionGranted: boolean, importStarted: boolean, importFailed: boolean, importFinished: boolean, hasImportLog: boolean };

export type StudyAdminListQueryVariables = Exact<{ [key: string]: never; }>;


export type StudyAdminListQuery = { __typename?: 'Query', userStudyUploadConfigured: boolean, studyAdminDetailsList: Array<{ __typename?: 'StudyAdminDetail', studyId: number, studyName: string, description: string, filename: string, cellCount: number, tissueNcitIds: Array<string>, diseaseMeshIds: Array<string>, visible: boolean, externalWebsite: string, readerPermissions: Array<string>, readerPermissionGranted: boolean, adminPermissions: Array<string>, adminPermissionGranted: boolean, importStarted: boolean, importFailed: boolean, importFinished: boolean, hasImportLog: boolean }> };

export type StudyLogsQueryVariables = Exact<{
  studyId: Scalars['Int'];
}>;


export type StudyLogsQuery = { __typename?: 'Query', studyImportLogsList: Array<{ __typename?: 'StudyImportLog', importFile: string, importLog: string }> };

export type GeneSpecificityQueryVariables = Exact<{
  studyId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  omicsIds: Array<Scalars['Int']> | Scalars['Int'];
  annotationGroupId: Scalars['Int'];
  secondAnnotationGroupId: Scalars['Int'];
  excludeAnnotationValueIds: Array<Scalars['Int']> | Scalars['Int'];
}>;


export type GeneSpecificityQuery = { __typename?: 'Query', expressionByTwoAnnotationsList: Array<{ __typename?: 'ExpressionByTwoAnnotation', omicsId: number, annotationValueId: number, annotationDisplayValue: string, secondAnnotationDisplayValue: string, secondAnnotationValueId: number, valueCount: number, exprSamplesFraction: number, mean: number, color: string }>, annotationGroupsList: Array<{ __typename?: 'AnnotationGroup', displayGroup: string, annotationGroupId: number }> };

export type StudyUpdateMutationVariables = Exact<{
  studyId: Scalars['Int'];
  studyName: Scalars['String'];
  description: InputMaybe<Scalars['String']>;
  readerPermissions: InputMaybe<Array<Scalars['String']> | Scalars['String']>;
  adminPermissions: InputMaybe<Array<Scalars['String']> | Scalars['String']>;
  tissueNcitIds: InputMaybe<Array<Scalars['String']> | Scalars['String']>;
  diseaseMeshIds: InputMaybe<Array<Scalars['String']> | Scalars['String']>;
  visible: Scalars['Boolean'];
  externalWebsite: InputMaybe<Scalars['String']>;
}>;


export type StudyUpdateMutation = { __typename?: 'Mutation', updateStudy: { __typename?: 'UpdateStudyPayload', clientMutationId: string } };

export type StudyDeleteMutationVariables = Exact<{
  studyId: Scalars['Int'];
}>;


export type StudyDeleteMutation = { __typename?: 'Mutation', deleteAllStudyData: { __typename?: 'DeleteAllStudyDataPayload', boolean: boolean } };

export type CreateStudyUploadMutationVariables = Exact<{
  filename: Scalars['String'];
}>;


export type CreateStudyUploadMutation = { __typename?: 'Mutation', createStudyUpload: { __typename?: 'CreateStudyUploadPayload', json: any } };

export type StudyDefinitionUpdateMutationVariables = Exact<{ [key: string]: never; }>;


export type StudyDefinitionUpdateMutation = { __typename?: 'Mutation', studyDefinitionUpdate: { __typename?: 'StudyDefinitionUpdatePayload', boolean: boolean } };

export const StudyInfoFragmentDoc = gql`
    fragment StudyInfo on StudyOverview {
  studyId
  studyName
  description
  cellCount
  externalWebsite
  metadata
  defaultStudyLayerId
  studyOntologyList {
    ontCodes
    labels
    ontology
    parentIds
  }
}
    `;
export const ReferenceStudyInfoFragmentDoc = gql`
    fragment ReferenceStudyInfo on ReferenceStudyOverview {
  organismTaxId
  studyName
  studyId
  cellCount
  defaultStudyLayerId
  externalWebsite
  description
  metadata
  studyOntologyList {
    ontCodes
    ontology
    labels
    parentIds
  }
  referenceStudyInfoList {
    celltypeAnnotationGroupId
    tissueAnnotationGroupId
  }
}
    `;
export const TreeOntologyOverviewFragmentDoc = gql`
    fragment TreeOntologyOverview on TreeOntology {
  label
  ontCode
  ontology
  parentOntCodePath
}
    `;
export const AnnotationGrpFragmentDoc = gql`
    fragment AnnotationGrp on StudyAnnotationFrontendGroup {
  annotationGroupId
  isPrimary
  ordering
  displayGroup
  differentialExpressionCalculated
  createdByUser
  currentUserIsOwner
  privateToUser
  annotationValuesList {
    annotationValueId
    displayValue
    color
    sampleCount
  }
}
    `;
export const DifferentialMarkerFragmentDoc = gql`
    fragment DifferentialMarker on DifferentialExpression {
  annotationValueId
  log2Foldchange
  pvalueAdj
  score
  study {
    studyName
    studyId
    cellCount
    description
  }
  annotationValue {
    annotationGroup {
      displayGroup
      annotationGroupId
    }
    displayValue
  }
  omics {
    displaySymbol
    taxId
    omicsId
    omicsType
    displayName
  }
}
    `;
export const OmicsGeneFragmentDoc = gql`
    fragment OmicsGene on OmicsBase {
  displayName
  displaySymbol
  omicsId
  taxId
}
    `;
export const StudyBasicsFragmentDoc = gql`
    fragment StudyBasics on Study {
  studyId
  studyName
  organismTaxId
  studyLayersList {
    layer
    studyLayerId
  }
  cellCount
  projections
}
    `;
export const OntologyOverviewFragmentDoc = gql`
    fragment ontologyOverview on Ontology {
  name
  ontid
  nodeId
}
    `;
export const DotPlotElementFragmentDoc = gql`
    fragment DotPlotElement on ExpressionByAnnotation {
  studyLayerId
  omicsId
  annotationValueId
  annotationDisplayValue
  q3
  median
  exprSamplesFraction
}
    `;
export const StudyAdminDetailsFragmentDoc = gql`
    fragment StudyAdminDetails on StudyAdminDetail {
  studyId
  studyName
  description
  filename
  cellCount
  tissueNcitIds
  diseaseMeshIds
  visible
  externalWebsite
  readerPermissions
  readerPermissionGranted
  adminPermissions
  adminPermissionGranted
  importStarted
  importFailed
  importFinished
  hasImportLog
}
    `;
export const CorrelatedgenesDocument = gql`
    query correlatedgenes($studyId: Int!, $omicsId: Int!) {
  getCorrelatedGenesList(studyId: $studyId, omicsId: $omicsId) {
    displayName
    displaySymbol
    omicsId
    r
  }
}
    `;

/**
 * __useCorrelatedgenesQuery__
 *
 * To run a query within a React component, call `useCorrelatedgenesQuery` and pass it any options that fit your needs.
 * When your component renders, `useCorrelatedgenesQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useCorrelatedgenesQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *      omicsId: // value for 'omicsId'
 *   },
 * });
 */
export function useCorrelatedgenesQuery(baseOptions: Apollo.QueryHookOptions<CorrelatedgenesQuery, CorrelatedgenesQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<CorrelatedgenesQuery, CorrelatedgenesQueryVariables>(CorrelatedgenesDocument, options);
      }
export function useCorrelatedgenesLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<CorrelatedgenesQuery, CorrelatedgenesQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<CorrelatedgenesQuery, CorrelatedgenesQueryVariables>(CorrelatedgenesDocument, options);
        }
export type CorrelatedgenesQueryHookResult = ReturnType<typeof useCorrelatedgenesQuery>;
export type CorrelatedgenesLazyQueryHookResult = ReturnType<typeof useCorrelatedgenesLazyQuery>;
export type CorrelatedgenesQueryResult = Apollo.QueryResult<CorrelatedgenesQuery, CorrelatedgenesQueryVariables>;
export const DegDocument = gql`
    query deg($studyId: Int!, $annotationValueId: Int!) {
  differentialExpressionVsList(
    filter: {annotationValueId: {equalTo: $annotationValueId}, studyId: {equalTo: $studyId}}
  ) {
    omicsId
    studyId
    annotationValueId
    omicsType
    displayName
    displaySymbol
    pvalueAdj
    log2Foldchange
    linkedGenes
  }
}
    `;

/**
 * __useDegQuery__
 *
 * To run a query within a React component, call `useDegQuery` and pass it any options that fit your needs.
 * When your component renders, `useDegQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useDegQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *      annotationValueId: // value for 'annotationValueId'
 *   },
 * });
 */
export function useDegQuery(baseOptions: Apollo.QueryHookOptions<DegQuery, DegQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<DegQuery, DegQueryVariables>(DegDocument, options);
      }
export function useDegLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<DegQuery, DegQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<DegQuery, DegQueryVariables>(DegDocument, options);
        }
export type DegQueryHookResult = ReturnType<typeof useDegQuery>;
export type DegLazyQueryHookResult = ReturnType<typeof useDegLazyQuery>;
export type DegQueryResult = Apollo.QueryResult<DegQuery, DegQueryVariables>;
export const StudiesDocument = gql`
    query studies {
  studyOverviewsList {
    ...StudyInfo
  }
}
    ${StudyInfoFragmentDoc}`;

/**
 * __useStudiesQuery__
 *
 * To run a query within a React component, call `useStudiesQuery` and pass it any options that fit your needs.
 * When your component renders, `useStudiesQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useStudiesQuery({
 *   variables: {
 *   },
 * });
 */
export function useStudiesQuery(baseOptions?: Apollo.QueryHookOptions<StudiesQuery, StudiesQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<StudiesQuery, StudiesQueryVariables>(StudiesDocument, options);
      }
export function useStudiesLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<StudiesQuery, StudiesQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<StudiesQuery, StudiesQueryVariables>(StudiesDocument, options);
        }
export type StudiesQueryHookResult = ReturnType<typeof useStudiesQuery>;
export type StudiesLazyQueryHookResult = ReturnType<typeof useStudiesLazyQuery>;
export type StudiesQueryResult = Apollo.QueryResult<StudiesQuery, StudiesQueryVariables>;
export const StudiesTreeOntologiesDocument = gql`
    query studiesTreeOntologies {
  treeOntologiesList {
    ...TreeOntologyOverview
  }
}
    ${TreeOntologyOverviewFragmentDoc}`;

/**
 * __useStudiesTreeOntologiesQuery__
 *
 * To run a query within a React component, call `useStudiesTreeOntologiesQuery` and pass it any options that fit your needs.
 * When your component renders, `useStudiesTreeOntologiesQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useStudiesTreeOntologiesQuery({
 *   variables: {
 *   },
 * });
 */
export function useStudiesTreeOntologiesQuery(baseOptions?: Apollo.QueryHookOptions<StudiesTreeOntologiesQuery, StudiesTreeOntologiesQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<StudiesTreeOntologiesQuery, StudiesTreeOntologiesQueryVariables>(StudiesTreeOntologiesDocument, options);
      }
export function useStudiesTreeOntologiesLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<StudiesTreeOntologiesQuery, StudiesTreeOntologiesQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<StudiesTreeOntologiesQuery, StudiesTreeOntologiesQueryVariables>(StudiesTreeOntologiesDocument, options);
        }
export type StudiesTreeOntologiesQueryHookResult = ReturnType<typeof useStudiesTreeOntologiesQuery>;
export type StudiesTreeOntologiesLazyQueryHookResult = ReturnType<typeof useStudiesTreeOntologiesLazyQuery>;
export type StudiesTreeOntologiesQueryResult = Apollo.QueryResult<StudiesTreeOntologiesQuery, StudiesTreeOntologiesQueryVariables>;
export const SingleStudyInfoDocument = gql`
    query singleStudyInfo($studyId: Int!) {
  studyOverviewsList(filter: {studyId: {equalTo: $studyId}}) {
    ...StudyInfo
  }
}
    ${StudyInfoFragmentDoc}`;

/**
 * __useSingleStudyInfoQuery__
 *
 * To run a query within a React component, call `useSingleStudyInfoQuery` and pass it any options that fit your needs.
 * When your component renders, `useSingleStudyInfoQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useSingleStudyInfoQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *   },
 * });
 */
export function useSingleStudyInfoQuery(baseOptions: Apollo.QueryHookOptions<SingleStudyInfoQuery, SingleStudyInfoQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<SingleStudyInfoQuery, SingleStudyInfoQueryVariables>(SingleStudyInfoDocument, options);
      }
export function useSingleStudyInfoLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<SingleStudyInfoQuery, SingleStudyInfoQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<SingleStudyInfoQuery, SingleStudyInfoQueryVariables>(SingleStudyInfoDocument, options);
        }
export type SingleStudyInfoQueryHookResult = ReturnType<typeof useSingleStudyInfoQuery>;
export type SingleStudyInfoLazyQueryHookResult = ReturnType<typeof useSingleStudyInfoLazyQuery>;
export type SingleStudyInfoQueryResult = Apollo.QueryResult<SingleStudyInfoQuery, SingleStudyInfoQueryVariables>;
export const ReferenceStudiesOverviewDocument = gql`
    query referenceStudiesOverview($organismTaxId: String!) {
  referenceStudyOverviewsList(condition: {organismTaxId: $organismTaxId}) {
    ...ReferenceStudyInfo
  }
}
    ${ReferenceStudyInfoFragmentDoc}`;

/**
 * __useReferenceStudiesOverviewQuery__
 *
 * To run a query within a React component, call `useReferenceStudiesOverviewQuery` and pass it any options that fit your needs.
 * When your component renders, `useReferenceStudiesOverviewQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useReferenceStudiesOverviewQuery({
 *   variables: {
 *      organismTaxId: // value for 'organismTaxId'
 *   },
 * });
 */
export function useReferenceStudiesOverviewQuery(baseOptions: Apollo.QueryHookOptions<ReferenceStudiesOverviewQuery, ReferenceStudiesOverviewQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<ReferenceStudiesOverviewQuery, ReferenceStudiesOverviewQueryVariables>(ReferenceStudiesOverviewDocument, options);
      }
export function useReferenceStudiesOverviewLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<ReferenceStudiesOverviewQuery, ReferenceStudiesOverviewQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<ReferenceStudiesOverviewQuery, ReferenceStudiesOverviewQueryVariables>(ReferenceStudiesOverviewDocument, options);
        }
export type ReferenceStudiesOverviewQueryHookResult = ReturnType<typeof useReferenceStudiesOverviewQuery>;
export type ReferenceStudiesOverviewLazyQueryHookResult = ReturnType<typeof useReferenceStudiesOverviewLazyQuery>;
export type ReferenceStudiesOverviewQueryResult = Apollo.QueryResult<ReferenceStudiesOverviewQuery, ReferenceStudiesOverviewQueryVariables>;
export const GeneSpecificityStudyDocument = gql`
    query geneSpecificityStudy($studyId: Int!) {
  referenceStudyOverviewsList(filter: {studyId: {equalTo: $studyId}}) {
    ...ReferenceStudyInfo
  }
  studyAnnotationFrontendGroupsList(condition: {studyId: $studyId}) {
    annotationGroupId
    annotationValuesList {
      annotationValueId
      displayValue
      color
      sampleCount
    }
  }
}
    ${ReferenceStudyInfoFragmentDoc}`;

/**
 * __useGeneSpecificityStudyQuery__
 *
 * To run a query within a React component, call `useGeneSpecificityStudyQuery` and pass it any options that fit your needs.
 * When your component renders, `useGeneSpecificityStudyQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useGeneSpecificityStudyQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *   },
 * });
 */
export function useGeneSpecificityStudyQuery(baseOptions: Apollo.QueryHookOptions<GeneSpecificityStudyQuery, GeneSpecificityStudyQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<GeneSpecificityStudyQuery, GeneSpecificityStudyQueryVariables>(GeneSpecificityStudyDocument, options);
      }
export function useGeneSpecificityStudyLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<GeneSpecificityStudyQuery, GeneSpecificityStudyQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<GeneSpecificityStudyQuery, GeneSpecificityStudyQueryVariables>(GeneSpecificityStudyDocument, options);
        }
export type GeneSpecificityStudyQueryHookResult = ReturnType<typeof useGeneSpecificityStudyQuery>;
export type GeneSpecificityStudyLazyQueryHookResult = ReturnType<typeof useGeneSpecificityStudyLazyQuery>;
export type GeneSpecificityStudyQueryResult = Apollo.QueryResult<GeneSpecificityStudyQuery, GeneSpecificityStudyQueryVariables>;
export const StudiesWithMarkerGenesDocument = gql`
    query studiesWithMarkerGenes($omicsIds: [Int!]!) {
  differentialExpressionsList(
    filter: {omicsId: {in: $omicsIds}}
    orderBy: LOG2_FOLDCHANGE_DESC
  ) {
    ...DifferentialMarker
  }
}
    ${DifferentialMarkerFragmentDoc}`;

/**
 * __useStudiesWithMarkerGenesQuery__
 *
 * To run a query within a React component, call `useStudiesWithMarkerGenesQuery` and pass it any options that fit your needs.
 * When your component renders, `useStudiesWithMarkerGenesQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useStudiesWithMarkerGenesQuery({
 *   variables: {
 *      omicsIds: // value for 'omicsIds'
 *   },
 * });
 */
export function useStudiesWithMarkerGenesQuery(baseOptions: Apollo.QueryHookOptions<StudiesWithMarkerGenesQuery, StudiesWithMarkerGenesQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<StudiesWithMarkerGenesQuery, StudiesWithMarkerGenesQueryVariables>(StudiesWithMarkerGenesDocument, options);
      }
export function useStudiesWithMarkerGenesLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<StudiesWithMarkerGenesQuery, StudiesWithMarkerGenesQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<StudiesWithMarkerGenesQuery, StudiesWithMarkerGenesQueryVariables>(StudiesWithMarkerGenesDocument, options);
        }
export type StudiesWithMarkerGenesQueryHookResult = ReturnType<typeof useStudiesWithMarkerGenesQuery>;
export type StudiesWithMarkerGenesLazyQueryHookResult = ReturnType<typeof useStudiesWithMarkerGenesLazyQuery>;
export type StudiesWithMarkerGenesQueryResult = Apollo.QueryResult<StudiesWithMarkerGenesQuery, StudiesWithMarkerGenesQueryVariables>;
export const AllGenesDocument = gql`
    query allGenes {
  omicsBasesList(filter: {omicsType: {equalTo: GENE}}) {
    ...OmicsGene
    value: displaySymbol
    ontology: omicsType
  }
}
    ${OmicsGeneFragmentDoc}`;

/**
 * __useAllGenesQuery__
 *
 * To run a query within a React component, call `useAllGenesQuery` and pass it any options that fit your needs.
 * When your component renders, `useAllGenesQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useAllGenesQuery({
 *   variables: {
 *   },
 * });
 */
export function useAllGenesQuery(baseOptions?: Apollo.QueryHookOptions<AllGenesQuery, AllGenesQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<AllGenesQuery, AllGenesQueryVariables>(AllGenesDocument, options);
      }
export function useAllGenesLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<AllGenesQuery, AllGenesQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<AllGenesQuery, AllGenesQueryVariables>(AllGenesDocument, options);
        }
export type AllGenesQueryHookResult = ReturnType<typeof useAllGenesQuery>;
export type AllGenesLazyQueryHookResult = ReturnType<typeof useAllGenesLazyQuery>;
export type AllGenesQueryResult = Apollo.QueryResult<AllGenesQuery, AllGenesQueryVariables>;
export const StudyOmicsDocument = gql`
    query studyOmics($studyId: Int!) {
  studyOmicsList(filter: {studyId: {equalTo: $studyId}}) {
    omics {
      omicsId
      displayName
      displaySymbol
    }
  }
}
    `;

/**
 * __useStudyOmicsQuery__
 *
 * To run a query within a React component, call `useStudyOmicsQuery` and pass it any options that fit your needs.
 * When your component renders, `useStudyOmicsQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useStudyOmicsQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *   },
 * });
 */
export function useStudyOmicsQuery(baseOptions: Apollo.QueryHookOptions<StudyOmicsQuery, StudyOmicsQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<StudyOmicsQuery, StudyOmicsQueryVariables>(StudyOmicsDocument, options);
      }
export function useStudyOmicsLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<StudyOmicsQuery, StudyOmicsQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<StudyOmicsQuery, StudyOmicsQueryVariables>(StudyOmicsDocument, options);
        }
export type StudyOmicsQueryHookResult = ReturnType<typeof useStudyOmicsQuery>;
export type StudyOmicsLazyQueryHookResult = ReturnType<typeof useStudyOmicsLazyQuery>;
export type StudyOmicsQueryResult = Apollo.QueryResult<StudyOmicsQuery, StudyOmicsQueryVariables>;
export const StudyBasicsDocument = gql`
    query StudyBasics($studyId: Int!) {
  study(studyId: $studyId) {
    studyId
    studyName
    cellCount
    projections
    organismTaxId
  }
  studyLayersList(condition: {studyId: $studyId}) {
    layer
    studyLayerId
  }
}
    `;

/**
 * __useStudyBasicsQuery__
 *
 * To run a query within a React component, call `useStudyBasicsQuery` and pass it any options that fit your needs.
 * When your component renders, `useStudyBasicsQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useStudyBasicsQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *   },
 * });
 */
export function useStudyBasicsQuery(baseOptions: Apollo.QueryHookOptions<StudyBasicsQuery, StudyBasicsQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<StudyBasicsQuery, StudyBasicsQueryVariables>(StudyBasicsDocument, options);
      }
export function useStudyBasicsLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<StudyBasicsQuery, StudyBasicsQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<StudyBasicsQuery, StudyBasicsQueryVariables>(StudyBasicsDocument, options);
        }
export type StudyBasicsQueryHookResult = ReturnType<typeof useStudyBasicsQuery>;
export type StudyBasicsLazyQueryHookResult = ReturnType<typeof useStudyBasicsLazyQuery>;
export type StudyBasicsQueryResult = Apollo.QueryResult<StudyBasicsQuery, StudyBasicsQueryVariables>;
export const StudyBasics2Document = gql`
    query StudyBasics2($studyId: Int!) {
  studyAnnotationFrontendGroupsList(condition: {studyId: $studyId}) {
    annotationGroupId
    isPrimary
    ordering
    displayGroup
    differentialExpressionCalculated
    createdByUser
    currentUserIsOwner
    privateToUser
    annotationValuesList {
      annotationValueId
      displayValue
      color
      sampleCount
    }
  }
  studySampleProjectionSubsamplingTransposedsList(condition: {studyId: $studyId}) {
    projectionType
    studySampleId
    projection
    modality
  }
}
    `;

/**
 * __useStudyBasics2Query__
 *
 * To run a query within a React component, call `useStudyBasics2Query` and pass it any options that fit your needs.
 * When your component renders, `useStudyBasics2Query` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useStudyBasics2Query({
 *   variables: {
 *      studyId: // value for 'studyId'
 *   },
 * });
 */
export function useStudyBasics2Query(baseOptions: Apollo.QueryHookOptions<StudyBasics2Query, StudyBasics2QueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<StudyBasics2Query, StudyBasics2QueryVariables>(StudyBasics2Document, options);
      }
export function useStudyBasics2LazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<StudyBasics2Query, StudyBasics2QueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<StudyBasics2Query, StudyBasics2QueryVariables>(StudyBasics2Document, options);
        }
export type StudyBasics2QueryHookResult = ReturnType<typeof useStudyBasics2Query>;
export type StudyBasics2LazyQueryHookResult = ReturnType<typeof useStudyBasics2LazyQuery>;
export type StudyBasics2QueryResult = Apollo.QueryResult<StudyBasics2Query, StudyBasics2QueryVariables>;
export const StudyBasics3Document = gql`
    query StudyBasics3($studyId: Int!) {
  studySampleAnnotationSubsamplingsList(condition: {studyId: $studyId}) {
    annotationValueId
    studySampleIds
  }
  studyOmicsTransposedsList(condition: {studyId: $studyId}) {
    displayName
    displaySymbol
    omicsId
    omicsType
  }
}
    `;

/**
 * __useStudyBasics3Query__
 *
 * To run a query within a React component, call `useStudyBasics3Query` and pass it any options that fit your needs.
 * When your component renders, `useStudyBasics3Query` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useStudyBasics3Query({
 *   variables: {
 *      studyId: // value for 'studyId'
 *   },
 * });
 */
export function useStudyBasics3Query(baseOptions: Apollo.QueryHookOptions<StudyBasics3Query, StudyBasics3QueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<StudyBasics3Query, StudyBasics3QueryVariables>(StudyBasics3Document, options);
      }
export function useStudyBasics3LazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<StudyBasics3Query, StudyBasics3QueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<StudyBasics3Query, StudyBasics3QueryVariables>(StudyBasics3Document, options);
        }
export type StudyBasics3QueryHookResult = ReturnType<typeof useStudyBasics3Query>;
export type StudyBasics3LazyQueryHookResult = ReturnType<typeof useStudyBasics3LazyQuery>;
export type StudyBasics3QueryResult = Apollo.QueryResult<StudyBasics3Query, StudyBasics3QueryVariables>;
export const ExpressionByOmicsIdsDocument = gql`
    query ExpressionByOmicsIds($studyLayerId: Int!, $omicsIds: [Int!]!, $subsamplingProjection: String) {
  expressionByOmicsIdsList(
    pStudyLayerId: $studyLayerId
    pOmicsIds: $omicsIds
    pSubsamplingProjection: $subsamplingProjection
  ) {
    omicsId
    studySampleIds
    values
  }
}
    `;

/**
 * __useExpressionByOmicsIdsQuery__
 *
 * To run a query within a React component, call `useExpressionByOmicsIdsQuery` and pass it any options that fit your needs.
 * When your component renders, `useExpressionByOmicsIdsQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useExpressionByOmicsIdsQuery({
 *   variables: {
 *      studyLayerId: // value for 'studyLayerId'
 *      omicsIds: // value for 'omicsIds'
 *      subsamplingProjection: // value for 'subsamplingProjection'
 *   },
 * });
 */
export function useExpressionByOmicsIdsQuery(baseOptions: Apollo.QueryHookOptions<ExpressionByOmicsIdsQuery, ExpressionByOmicsIdsQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<ExpressionByOmicsIdsQuery, ExpressionByOmicsIdsQueryVariables>(ExpressionByOmicsIdsDocument, options);
      }
export function useExpressionByOmicsIdsLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<ExpressionByOmicsIdsQuery, ExpressionByOmicsIdsQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<ExpressionByOmicsIdsQuery, ExpressionByOmicsIdsQueryVariables>(ExpressionByOmicsIdsDocument, options);
        }
export type ExpressionByOmicsIdsQueryHookResult = ReturnType<typeof useExpressionByOmicsIdsQuery>;
export type ExpressionByOmicsIdsLazyQueryHookResult = ReturnType<typeof useExpressionByOmicsIdsLazyQuery>;
export type ExpressionByOmicsIdsQueryResult = Apollo.QueryResult<ExpressionByOmicsIdsQuery, ExpressionByOmicsIdsQueryVariables>;
export const ExpressionViolinPlotDocument = gql`
    query ExpressionViolinPlot($studyId: Int!, $studyLayerId: Int!, $omicsId: Int!, $annotationGroupId: Int!, $excludeAnnotationValueIds: [Int!]!, $annotationSecondaryGroupId: Int) {
  violinPlot(
    pStudyId: $studyId
    pStudyLayerId: $studyLayerId
    pOmicsId: $omicsId
    pAnnotationGroupId: $annotationGroupId
    pExcludeAnnotationValueIds: $excludeAnnotationValueIds
    pSecondaryAnnotationGroupId: $annotationSecondaryGroupId
  )
}
    `;

/**
 * __useExpressionViolinPlotQuery__
 *
 * To run a query within a React component, call `useExpressionViolinPlotQuery` and pass it any options that fit your needs.
 * When your component renders, `useExpressionViolinPlotQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useExpressionViolinPlotQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *      studyLayerId: // value for 'studyLayerId'
 *      omicsId: // value for 'omicsId'
 *      annotationGroupId: // value for 'annotationGroupId'
 *      excludeAnnotationValueIds: // value for 'excludeAnnotationValueIds'
 *      annotationSecondaryGroupId: // value for 'annotationSecondaryGroupId'
 *   },
 * });
 */
export function useExpressionViolinPlotQuery(baseOptions: Apollo.QueryHookOptions<ExpressionViolinPlotQuery, ExpressionViolinPlotQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<ExpressionViolinPlotQuery, ExpressionViolinPlotQueryVariables>(ExpressionViolinPlotDocument, options);
      }
export function useExpressionViolinPlotLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<ExpressionViolinPlotQuery, ExpressionViolinPlotQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<ExpressionViolinPlotQuery, ExpressionViolinPlotQueryVariables>(ExpressionViolinPlotDocument, options);
        }
export type ExpressionViolinPlotQueryHookResult = ReturnType<typeof useExpressionViolinPlotQuery>;
export type ExpressionViolinPlotLazyQueryHookResult = ReturnType<typeof useExpressionViolinPlotLazyQuery>;
export type ExpressionViolinPlotQueryResult = Apollo.QueryResult<ExpressionViolinPlotQuery, ExpressionViolinPlotQueryVariables>;
export const ExpressionGroupTableDocument = gql`
    query ExpressionGroupTable($studyId: Int!, $studyLayerId: Int!, $omicsId: Int!, $annotationGroupId: Int!, $excludeAnnotationValueIds: [Int!]!, $annotationSecondaryGroupId: Int!) {
  expressionByTwoAnnotationsList(
    pStudyId: $studyId
    pStudyLayerId: $studyLayerId
    pOmicsIds: [$omicsId]
    pAnnotationGroupId: $annotationGroupId
    pSecondAnnotationGroupId: $annotationSecondaryGroupId
    pExcludeAnnotationValueIds: $excludeAnnotationValueIds
    pDropoutsAsZero: true
  ) {
    annotationValueId
    annotationDisplayValue
    secondAnnotationValueId
    secondAnnotationDisplayValue
    valueCount
    median
    nonZeroValueCount
  }
}
    `;

/**
 * __useExpressionGroupTableQuery__
 *
 * To run a query within a React component, call `useExpressionGroupTableQuery` and pass it any options that fit your needs.
 * When your component renders, `useExpressionGroupTableQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useExpressionGroupTableQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *      studyLayerId: // value for 'studyLayerId'
 *      omicsId: // value for 'omicsId'
 *      annotationGroupId: // value for 'annotationGroupId'
 *      excludeAnnotationValueIds: // value for 'excludeAnnotationValueIds'
 *      annotationSecondaryGroupId: // value for 'annotationSecondaryGroupId'
 *   },
 * });
 */
export function useExpressionGroupTableQuery(baseOptions: Apollo.QueryHookOptions<ExpressionGroupTableQuery, ExpressionGroupTableQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<ExpressionGroupTableQuery, ExpressionGroupTableQueryVariables>(ExpressionGroupTableDocument, options);
      }
export function useExpressionGroupTableLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<ExpressionGroupTableQuery, ExpressionGroupTableQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<ExpressionGroupTableQuery, ExpressionGroupTableQueryVariables>(ExpressionGroupTableDocument, options);
        }
export type ExpressionGroupTableQueryHookResult = ReturnType<typeof useExpressionGroupTableQuery>;
export type ExpressionGroupTableLazyQueryHookResult = ReturnType<typeof useExpressionGroupTableLazyQuery>;
export type ExpressionGroupTableQueryResult = Apollo.QueryResult<ExpressionGroupTableQuery, ExpressionGroupTableQueryVariables>;
export const ExpressionTTestDocument = gql`
    query ExpressionTTest($studyId: Int!, $studyLayerId: Int!, $omicsId: Int!, $annotationGroupId: Int!, $excludeAnnotationValueIds: [Int!]!, $secondAnnotationGroupId: Int, $sample1AnnotationValueId: Int!, $sample1SecondAnnotationValueId: Int, $sample2AnnotationValueId: Int!, $sample2SecondAnnotationValueId: Int) {
  expressionTtest(
    pStudyId: $studyId
    pStudyLayerId: $studyLayerId
    pOmicsId: $omicsId
    pAnnotationGroupId: $annotationGroupId
    pSecondaryAnnotationGroupId: $secondAnnotationGroupId
    pExcludeAnnotationValueIds: $excludeAnnotationValueIds
    pSample1AnnotationValueId: $sample1AnnotationValueId
    pSample1SecondAnnotationValueId: $sample1SecondAnnotationValueId
    pSample2AnnotationValueId: $sample2AnnotationValueId
    pSample2SecondAnnotationValueId: $sample2SecondAnnotationValueId
  )
}
    `;

/**
 * __useExpressionTTestQuery__
 *
 * To run a query within a React component, call `useExpressionTTestQuery` and pass it any options that fit your needs.
 * When your component renders, `useExpressionTTestQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useExpressionTTestQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *      studyLayerId: // value for 'studyLayerId'
 *      omicsId: // value for 'omicsId'
 *      annotationGroupId: // value for 'annotationGroupId'
 *      excludeAnnotationValueIds: // value for 'excludeAnnotationValueIds'
 *      secondAnnotationGroupId: // value for 'secondAnnotationGroupId'
 *      sample1AnnotationValueId: // value for 'sample1AnnotationValueId'
 *      sample1SecondAnnotationValueId: // value for 'sample1SecondAnnotationValueId'
 *      sample2AnnotationValueId: // value for 'sample2AnnotationValueId'
 *      sample2SecondAnnotationValueId: // value for 'sample2SecondAnnotationValueId'
 *   },
 * });
 */
export function useExpressionTTestQuery(baseOptions: Apollo.QueryHookOptions<ExpressionTTestQuery, ExpressionTTestQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<ExpressionTTestQuery, ExpressionTTestQueryVariables>(ExpressionTTestDocument, options);
      }
export function useExpressionTTestLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<ExpressionTTestQuery, ExpressionTTestQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<ExpressionTTestQuery, ExpressionTTestQueryVariables>(ExpressionTTestDocument, options);
        }
export type ExpressionTTestQueryHookResult = ReturnType<typeof useExpressionTTestQuery>;
export type ExpressionTTestLazyQueryHookResult = ReturnType<typeof useExpressionTTestLazyQuery>;
export type ExpressionTTestQueryResult = Apollo.QueryResult<ExpressionTTestQuery, ExpressionTTestQueryVariables>;
export const ExpressionCorrelationTrianglePlotDocument = gql`
    query ExpressionCorrelationTrianglePlot($studyId: Int!, $studyLayerId: Int!, $omicsIds: [Int!]!, $excludeAnnotationValueIds: [Int!]!) {
  correlationTrianglePlot(
    pStudyId: $studyId
    pStudyLayerId: $studyLayerId
    pOmicsIds: $omicsIds
    pExcludeAnnotationValueIds: $excludeAnnotationValueIds
  )
}
    `;

/**
 * __useExpressionCorrelationTrianglePlotQuery__
 *
 * To run a query within a React component, call `useExpressionCorrelationTrianglePlotQuery` and pass it any options that fit your needs.
 * When your component renders, `useExpressionCorrelationTrianglePlotQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useExpressionCorrelationTrianglePlotQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *      studyLayerId: // value for 'studyLayerId'
 *      omicsIds: // value for 'omicsIds'
 *      excludeAnnotationValueIds: // value for 'excludeAnnotationValueIds'
 *   },
 * });
 */
export function useExpressionCorrelationTrianglePlotQuery(baseOptions: Apollo.QueryHookOptions<ExpressionCorrelationTrianglePlotQuery, ExpressionCorrelationTrianglePlotQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<ExpressionCorrelationTrianglePlotQuery, ExpressionCorrelationTrianglePlotQueryVariables>(ExpressionCorrelationTrianglePlotDocument, options);
      }
export function useExpressionCorrelationTrianglePlotLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<ExpressionCorrelationTrianglePlotQuery, ExpressionCorrelationTrianglePlotQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<ExpressionCorrelationTrianglePlotQuery, ExpressionCorrelationTrianglePlotQueryVariables>(ExpressionCorrelationTrianglePlotDocument, options);
        }
export type ExpressionCorrelationTrianglePlotQueryHookResult = ReturnType<typeof useExpressionCorrelationTrianglePlotQuery>;
export type ExpressionCorrelationTrianglePlotLazyQueryHookResult = ReturnType<typeof useExpressionCorrelationTrianglePlotLazyQuery>;
export type ExpressionCorrelationTrianglePlotQueryResult = Apollo.QueryResult<ExpressionCorrelationTrianglePlotQuery, ExpressionCorrelationTrianglePlotQueryVariables>;
export const AutocompleteDocument = gql`
    query autocomplete($query: String!) {
  autocompleteList(searchQuery: $query, first: 20) {
    isSynonymOfPreferredTerm
    label
    labelHighlight
    ontCode
    ontology
  }
}
    `;

/**
 * __useAutocompleteQuery__
 *
 * To run a query within a React component, call `useAutocompleteQuery` and pass it any options that fit your needs.
 * When your component renders, `useAutocompleteQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useAutocompleteQuery({
 *   variables: {
 *      query: // value for 'query'
 *   },
 * });
 */
export function useAutocompleteQuery(baseOptions: Apollo.QueryHookOptions<AutocompleteQuery, AutocompleteQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<AutocompleteQuery, AutocompleteQueryVariables>(AutocompleteDocument, options);
      }
export function useAutocompleteLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<AutocompleteQuery, AutocompleteQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<AutocompleteQuery, AutocompleteQueryVariables>(AutocompleteDocument, options);
        }
export type AutocompleteQueryHookResult = ReturnType<typeof useAutocompleteQuery>;
export type AutocompleteLazyQueryHookResult = ReturnType<typeof useAutocompleteLazyQuery>;
export type AutocompleteQueryResult = Apollo.QueryResult<AutocompleteQuery, AutocompleteQueryVariables>;
export const OmicsAutocompleteDocument = gql`
    query OmicsAutocomplete($searchQuery: String!, $omicsType: OmicsType!, $taxId: Int!) {
  omicsAutocompleteList(
    searchQuery: $searchQuery
    omicsTypeFilter: $omicsType
    taxIdFilter: $taxId
  ) {
    displaySymbol
    labelHighlight
    omicsId
    omicsType
  }
}
    `;

/**
 * __useOmicsAutocompleteQuery__
 *
 * To run a query within a React component, call `useOmicsAutocompleteQuery` and pass it any options that fit your needs.
 * When your component renders, `useOmicsAutocompleteQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useOmicsAutocompleteQuery({
 *   variables: {
 *      searchQuery: // value for 'searchQuery'
 *      omicsType: // value for 'omicsType'
 *      taxId: // value for 'taxId'
 *   },
 * });
 */
export function useOmicsAutocompleteQuery(baseOptions: Apollo.QueryHookOptions<OmicsAutocompleteQuery, OmicsAutocompleteQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<OmicsAutocompleteQuery, OmicsAutocompleteQueryVariables>(OmicsAutocompleteDocument, options);
      }
export function useOmicsAutocompleteLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<OmicsAutocompleteQuery, OmicsAutocompleteQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<OmicsAutocompleteQuery, OmicsAutocompleteQueryVariables>(OmicsAutocompleteDocument, options);
        }
export type OmicsAutocompleteQueryHookResult = ReturnType<typeof useOmicsAutocompleteQuery>;
export type OmicsAutocompleteLazyQueryHookResult = ReturnType<typeof useOmicsAutocompleteLazyQuery>;
export type OmicsAutocompleteQueryResult = Apollo.QueryResult<OmicsAutocompleteQuery, OmicsAutocompleteQueryVariables>;
export const OntologiesDocument = gql`
    query ontologies {
  ontologiesList {
    ...ontologyOverview
  }
}
    ${OntologyOverviewFragmentDoc}`;

/**
 * __useOntologiesQuery__
 *
 * To run a query within a React component, call `useOntologiesQuery` and pass it any options that fit your needs.
 * When your component renders, `useOntologiesQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useOntologiesQuery({
 *   variables: {
 *   },
 * });
 */
export function useOntologiesQuery(baseOptions?: Apollo.QueryHookOptions<OntologiesQuery, OntologiesQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<OntologiesQuery, OntologiesQueryVariables>(OntologiesDocument, options);
      }
export function useOntologiesLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<OntologiesQuery, OntologiesQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<OntologiesQuery, OntologiesQueryVariables>(OntologiesDocument, options);
        }
export type OntologiesQueryHookResult = ReturnType<typeof useOntologiesQuery>;
export type OntologiesLazyQueryHookResult = ReturnType<typeof useOntologiesLazyQuery>;
export type OntologiesQueryResult = Apollo.QueryResult<OntologiesQuery, OntologiesQueryVariables>;
export const ExpressionByAnnotationDocument = gql`
    query expressionByAnnotation($studyLayerIds: [Int!]!, $omicsIds: [Int!]!, $annotationGroupId: Int!, $excludeAnnotationValueIds: [Int!]!) {
  expressionByAnnotationList(
    pStudyLayerIds: $studyLayerIds
    pOmicsIds: $omicsIds
    pAnnotationGroupId: $annotationGroupId
    pExcludeAnnotationValueIds: $excludeAnnotationValueIds
  ) {
    ...DotPlotElement
  }
}
    ${DotPlotElementFragmentDoc}`;

/**
 * __useExpressionByAnnotationQuery__
 *
 * To run a query within a React component, call `useExpressionByAnnotationQuery` and pass it any options that fit your needs.
 * When your component renders, `useExpressionByAnnotationQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useExpressionByAnnotationQuery({
 *   variables: {
 *      studyLayerIds: // value for 'studyLayerIds'
 *      omicsIds: // value for 'omicsIds'
 *      annotationGroupId: // value for 'annotationGroupId'
 *      excludeAnnotationValueIds: // value for 'excludeAnnotationValueIds'
 *   },
 * });
 */
export function useExpressionByAnnotationQuery(baseOptions: Apollo.QueryHookOptions<ExpressionByAnnotationQuery, ExpressionByAnnotationQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<ExpressionByAnnotationQuery, ExpressionByAnnotationQueryVariables>(ExpressionByAnnotationDocument, options);
      }
export function useExpressionByAnnotationLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<ExpressionByAnnotationQuery, ExpressionByAnnotationQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<ExpressionByAnnotationQuery, ExpressionByAnnotationQueryVariables>(ExpressionByAnnotationDocument, options);
        }
export type ExpressionByAnnotationQueryHookResult = ReturnType<typeof useExpressionByAnnotationQuery>;
export type ExpressionByAnnotationLazyQueryHookResult = ReturnType<typeof useExpressionByAnnotationLazyQuery>;
export type ExpressionByAnnotationQueryResult = Apollo.QueryResult<ExpressionByAnnotationQuery, ExpressionByAnnotationQueryVariables>;
export const CellOAnnotationGroupIdDocument = gql`
    query CellOAnnotationGroupId {
  annotationGroupsList(filter: {h5AdColumn: {equalTo: "CellO_celltype"}}) {
    annotationGroupId
  }
}
    `;

/**
 * __useCellOAnnotationGroupIdQuery__
 *
 * To run a query within a React component, call `useCellOAnnotationGroupIdQuery` and pass it any options that fit your needs.
 * When your component renders, `useCellOAnnotationGroupIdQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useCellOAnnotationGroupIdQuery({
 *   variables: {
 *   },
 * });
 */
export function useCellOAnnotationGroupIdQuery(baseOptions?: Apollo.QueryHookOptions<CellOAnnotationGroupIdQuery, CellOAnnotationGroupIdQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<CellOAnnotationGroupIdQuery, CellOAnnotationGroupIdQueryVariables>(CellOAnnotationGroupIdDocument, options);
      }
export function useCellOAnnotationGroupIdLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<CellOAnnotationGroupIdQuery, CellOAnnotationGroupIdQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<CellOAnnotationGroupIdQuery, CellOAnnotationGroupIdQueryVariables>(CellOAnnotationGroupIdDocument, options);
        }
export type CellOAnnotationGroupIdQueryHookResult = ReturnType<typeof useCellOAnnotationGroupIdQuery>;
export type CellOAnnotationGroupIdLazyQueryHookResult = ReturnType<typeof useCellOAnnotationGroupIdLazyQuery>;
export type CellOAnnotationGroupIdQueryResult = Apollo.QueryResult<CellOAnnotationGroupIdQuery, CellOAnnotationGroupIdQueryVariables>;
export const HalfAVolcanoDocument = gql`
    query halfAVolcano($annotationValueId: Int!, $studyId: Int!) {
  differentialExpressionsList(
    filter: {annotationValueId: {equalTo: $annotationValueId}, studyId: {equalTo: $studyId}}
  ) {
    log2Foldchange
    pvalueAdj
  }
}
    `;

/**
 * __useHalfAVolcanoQuery__
 *
 * To run a query within a React component, call `useHalfAVolcanoQuery` and pass it any options that fit your needs.
 * When your component renders, `useHalfAVolcanoQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useHalfAVolcanoQuery({
 *   variables: {
 *      annotationValueId: // value for 'annotationValueId'
 *      studyId: // value for 'studyId'
 *   },
 * });
 */
export function useHalfAVolcanoQuery(baseOptions: Apollo.QueryHookOptions<HalfAVolcanoQuery, HalfAVolcanoQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<HalfAVolcanoQuery, HalfAVolcanoQueryVariables>(HalfAVolcanoDocument, options);
      }
export function useHalfAVolcanoLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<HalfAVolcanoQuery, HalfAVolcanoQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<HalfAVolcanoQuery, HalfAVolcanoQueryVariables>(HalfAVolcanoDocument, options);
        }
export type HalfAVolcanoQueryHookResult = ReturnType<typeof useHalfAVolcanoQuery>;
export type HalfAVolcanoLazyQueryHookResult = ReturnType<typeof useHalfAVolcanoLazyQuery>;
export type HalfAVolcanoQueryResult = Apollo.QueryResult<HalfAVolcanoQuery, HalfAVolcanoQueryVariables>;
export const AnnotationValueCoocurrenceDocument = gql`
    query annotationValueCoocurrence($studyId: Int!, $annotationGroupId1: Int!, $annotationGroupId2: Int!) {
  annotationValueCoocurrenceList(
    studyId: $studyId
    annotationGroupId1: $annotationGroupId1
    annotationGroupId2: $annotationGroupId2
  ) {
    annotationValueId1
    annotationValueId2
    occurrence
  }
}
    `;

/**
 * __useAnnotationValueCoocurrenceQuery__
 *
 * To run a query within a React component, call `useAnnotationValueCoocurrenceQuery` and pass it any options that fit your needs.
 * When your component renders, `useAnnotationValueCoocurrenceQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useAnnotationValueCoocurrenceQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *      annotationGroupId1: // value for 'annotationGroupId1'
 *      annotationGroupId2: // value for 'annotationGroupId2'
 *   },
 * });
 */
export function useAnnotationValueCoocurrenceQuery(baseOptions: Apollo.QueryHookOptions<AnnotationValueCoocurrenceQuery, AnnotationValueCoocurrenceQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<AnnotationValueCoocurrenceQuery, AnnotationValueCoocurrenceQueryVariables>(AnnotationValueCoocurrenceDocument, options);
      }
export function useAnnotationValueCoocurrenceLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<AnnotationValueCoocurrenceQuery, AnnotationValueCoocurrenceQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<AnnotationValueCoocurrenceQuery, AnnotationValueCoocurrenceQueryVariables>(AnnotationValueCoocurrenceDocument, options);
        }
export type AnnotationValueCoocurrenceQueryHookResult = ReturnType<typeof useAnnotationValueCoocurrenceQuery>;
export type AnnotationValueCoocurrenceLazyQueryHookResult = ReturnType<typeof useAnnotationValueCoocurrenceLazyQuery>;
export type AnnotationValueCoocurrenceQueryResult = Apollo.QueryResult<AnnotationValueCoocurrenceQuery, AnnotationValueCoocurrenceQueryVariables>;
export const SaveUserAnnotationDocument = gql`
    mutation SaveUserAnnotation($studyId: Int!, $annotationGroupName: String!, $selectedSampleIds: String!, $unexpressedSamplesOmicsIds: [Int!]) {
  userAnnotationDefine(
    input: {pStudyId: $studyId, pAnnotationGroupName: $annotationGroupName, pSelectedSampleIds: $selectedSampleIds, pUnexpressedSamplesOmicsIds: $unexpressedSamplesOmicsIds}
  ) {
    clientMutationId
    integer
  }
}
    `;

/**
 * __useSaveUserAnnotationMutation__
 *
 * To run a mutation, you first call `useSaveUserAnnotationMutation` within a React component and pass it any options that fit your needs.
 * When your component renders, `useSaveUserAnnotationMutation` returns a tuple that includes:
 * - A mutate function that you can call at any time to execute the mutation
 * - An object with fields that represent the current status of the mutation's execution
 *
 * @param baseOptions options that will be passed into the mutation, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options-2;
 *
 * @example
 * const [saveUserAnnotationMutation, { data, loading, error }] = useSaveUserAnnotationMutation({
 *   variables: {
 *      studyId: // value for 'studyId'
 *      annotationGroupName: // value for 'annotationGroupName'
 *      selectedSampleIds: // value for 'selectedSampleIds'
 *      unexpressedSamplesOmicsIds: // value for 'unexpressedSamplesOmicsIds'
 *   },
 * });
 */
export function useSaveUserAnnotationMutation(baseOptions?: Apollo.MutationHookOptions<SaveUserAnnotationMutation, SaveUserAnnotationMutationVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useMutation<SaveUserAnnotationMutation, SaveUserAnnotationMutationVariables>(SaveUserAnnotationDocument, options);
      }
export type SaveUserAnnotationMutationHookResult = ReturnType<typeof useSaveUserAnnotationMutation>;
export type SaveUserAnnotationMutationResult = Apollo.MutationResult<SaveUserAnnotationMutation>;
export type SaveUserAnnotationMutationOptions = Apollo.BaseMutationOptions<SaveUserAnnotationMutation, SaveUserAnnotationMutationVariables>;
export const EditUserAnnotationDocument = gql`
    mutation EditUserAnnotation($studyId: Int!, $annotationGroupId: Int!, $privateToUser: Boolean!) {
  userAnnotationEdit(
    input: {pStudyId: $studyId, pAnnotationGroupId: $annotationGroupId, pPrivateToUser: $privateToUser}
  ) {
    boolean
  }
}
    `;

/**
 * __useEditUserAnnotationMutation__
 *
 * To run a mutation, you first call `useEditUserAnnotationMutation` within a React component and pass it any options that fit your needs.
 * When your component renders, `useEditUserAnnotationMutation` returns a tuple that includes:
 * - A mutate function that you can call at any time to execute the mutation
 * - An object with fields that represent the current status of the mutation's execution
 *
 * @param baseOptions options that will be passed into the mutation, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options-2;
 *
 * @example
 * const [editUserAnnotationMutation, { data, loading, error }] = useEditUserAnnotationMutation({
 *   variables: {
 *      studyId: // value for 'studyId'
 *      annotationGroupId: // value for 'annotationGroupId'
 *      privateToUser: // value for 'privateToUser'
 *   },
 * });
 */
export function useEditUserAnnotationMutation(baseOptions?: Apollo.MutationHookOptions<EditUserAnnotationMutation, EditUserAnnotationMutationVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useMutation<EditUserAnnotationMutation, EditUserAnnotationMutationVariables>(EditUserAnnotationDocument, options);
      }
export type EditUserAnnotationMutationHookResult = ReturnType<typeof useEditUserAnnotationMutation>;
export type EditUserAnnotationMutationResult = Apollo.MutationResult<EditUserAnnotationMutation>;
export type EditUserAnnotationMutationOptions = Apollo.BaseMutationOptions<EditUserAnnotationMutation, EditUserAnnotationMutationVariables>;
export const DeleteUserAnnotationDocument = gql`
    mutation DeleteUserAnnotation($studyId: Int!, $annotationGroupId: Int!) {
  userAnnotationDelete(
    input: {pStudyId: $studyId, pAnnotationGroupId: $annotationGroupId}
  ) {
    boolean
  }
}
    `;

/**
 * __useDeleteUserAnnotationMutation__
 *
 * To run a mutation, you first call `useDeleteUserAnnotationMutation` within a React component and pass it any options that fit your needs.
 * When your component renders, `useDeleteUserAnnotationMutation` returns a tuple that includes:
 * - A mutate function that you can call at any time to execute the mutation
 * - An object with fields that represent the current status of the mutation's execution
 *
 * @param baseOptions options that will be passed into the mutation, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options-2;
 *
 * @example
 * const [deleteUserAnnotationMutation, { data, loading, error }] = useDeleteUserAnnotationMutation({
 *   variables: {
 *      studyId: // value for 'studyId'
 *      annotationGroupId: // value for 'annotationGroupId'
 *   },
 * });
 */
export function useDeleteUserAnnotationMutation(baseOptions?: Apollo.MutationHookOptions<DeleteUserAnnotationMutation, DeleteUserAnnotationMutationVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useMutation<DeleteUserAnnotationMutation, DeleteUserAnnotationMutationVariables>(DeleteUserAnnotationDocument, options);
      }
export type DeleteUserAnnotationMutationHookResult = ReturnType<typeof useDeleteUserAnnotationMutation>;
export type DeleteUserAnnotationMutationResult = Apollo.MutationResult<DeleteUserAnnotationMutation>;
export type DeleteUserAnnotationMutationOptions = Apollo.BaseMutationOptions<DeleteUserAnnotationMutation, DeleteUserAnnotationMutationVariables>;
export const StudyAdminListDocument = gql`
    query studyAdminList {
  studyAdminDetailsList {
    ...StudyAdminDetails
  }
  userStudyUploadConfigured
}
    ${StudyAdminDetailsFragmentDoc}`;

/**
 * __useStudyAdminListQuery__
 *
 * To run a query within a React component, call `useStudyAdminListQuery` and pass it any options that fit your needs.
 * When your component renders, `useStudyAdminListQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useStudyAdminListQuery({
 *   variables: {
 *   },
 * });
 */
export function useStudyAdminListQuery(baseOptions?: Apollo.QueryHookOptions<StudyAdminListQuery, StudyAdminListQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<StudyAdminListQuery, StudyAdminListQueryVariables>(StudyAdminListDocument, options);
      }
export function useStudyAdminListLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<StudyAdminListQuery, StudyAdminListQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<StudyAdminListQuery, StudyAdminListQueryVariables>(StudyAdminListDocument, options);
        }
export type StudyAdminListQueryHookResult = ReturnType<typeof useStudyAdminListQuery>;
export type StudyAdminListLazyQueryHookResult = ReturnType<typeof useStudyAdminListLazyQuery>;
export type StudyAdminListQueryResult = Apollo.QueryResult<StudyAdminListQuery, StudyAdminListQueryVariables>;
export const StudyLogsDocument = gql`
    query studyLogs($studyId: Int!) {
  studyImportLogsList(condition: {studyId: $studyId}) {
    importFile
    importLog
  }
}
    `;

/**
 * __useStudyLogsQuery__
 *
 * To run a query within a React component, call `useStudyLogsQuery` and pass it any options that fit your needs.
 * When your component renders, `useStudyLogsQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useStudyLogsQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *   },
 * });
 */
export function useStudyLogsQuery(baseOptions: Apollo.QueryHookOptions<StudyLogsQuery, StudyLogsQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<StudyLogsQuery, StudyLogsQueryVariables>(StudyLogsDocument, options);
      }
export function useStudyLogsLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<StudyLogsQuery, StudyLogsQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<StudyLogsQuery, StudyLogsQueryVariables>(StudyLogsDocument, options);
        }
export type StudyLogsQueryHookResult = ReturnType<typeof useStudyLogsQuery>;
export type StudyLogsLazyQueryHookResult = ReturnType<typeof useStudyLogsLazyQuery>;
export type StudyLogsQueryResult = Apollo.QueryResult<StudyLogsQuery, StudyLogsQueryVariables>;
export const GeneSpecificityDocument = gql`
    query GeneSpecificity($studyId: Int!, $studyLayerId: Int!, $omicsIds: [Int!]!, $annotationGroupId: Int!, $secondAnnotationGroupId: Int!, $excludeAnnotationValueIds: [Int!]!) {
  expressionByTwoAnnotationsList(
    pStudyId: $studyId
    pStudyLayerId: $studyLayerId
    pOmicsIds: $omicsIds
    pAnnotationGroupId: $annotationGroupId
    pExcludeAnnotationValueIds: $excludeAnnotationValueIds
    pSecondAnnotationGroupId: $secondAnnotationGroupId
    pDropoutsAsZero: false
  ) {
    omicsId
    annotationValueId
    annotationDisplayValue
    secondAnnotationDisplayValue
    secondAnnotationValueId
    valueCount
    exprSamplesFraction
    mean
    color
  }
  annotationGroupsList(
    filter: {annotationGroupId: {in: [$annotationGroupId, $secondAnnotationGroupId]}}
  ) {
    displayGroup
    annotationGroupId
  }
}
    `;

/**
 * __useGeneSpecificityQuery__
 *
 * To run a query within a React component, call `useGeneSpecificityQuery` and pass it any options that fit your needs.
 * When your component renders, `useGeneSpecificityQuery` returns an object from Apollo Client that contains loading, error, and data properties
 * you can use to render your UI.
 *
 * @param baseOptions options that will be passed into the query, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options;
 *
 * @example
 * const { data, loading, error } = useGeneSpecificityQuery({
 *   variables: {
 *      studyId: // value for 'studyId'
 *      studyLayerId: // value for 'studyLayerId'
 *      omicsIds: // value for 'omicsIds'
 *      annotationGroupId: // value for 'annotationGroupId'
 *      secondAnnotationGroupId: // value for 'secondAnnotationGroupId'
 *      excludeAnnotationValueIds: // value for 'excludeAnnotationValueIds'
 *   },
 * });
 */
export function useGeneSpecificityQuery(baseOptions: Apollo.QueryHookOptions<GeneSpecificityQuery, GeneSpecificityQueryVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useQuery<GeneSpecificityQuery, GeneSpecificityQueryVariables>(GeneSpecificityDocument, options);
      }
export function useGeneSpecificityLazyQuery(baseOptions?: Apollo.LazyQueryHookOptions<GeneSpecificityQuery, GeneSpecificityQueryVariables>) {
          const options = {...defaultOptions, ...baseOptions}
          return Apollo.useLazyQuery<GeneSpecificityQuery, GeneSpecificityQueryVariables>(GeneSpecificityDocument, options);
        }
export type GeneSpecificityQueryHookResult = ReturnType<typeof useGeneSpecificityQuery>;
export type GeneSpecificityLazyQueryHookResult = ReturnType<typeof useGeneSpecificityLazyQuery>;
export type GeneSpecificityQueryResult = Apollo.QueryResult<GeneSpecificityQuery, GeneSpecificityQueryVariables>;
export const StudyUpdateDocument = gql`
    mutation studyUpdate($studyId: Int!, $studyName: String!, $description: String, $readerPermissions: [String!], $adminPermissions: [String!], $tissueNcitIds: [String!], $diseaseMeshIds: [String!], $visible: Boolean!, $externalWebsite: String) {
  updateStudy(
    input: {studyId: $studyId, patch: {studyName: $studyName, description: $description, readerPermissions: $readerPermissions, adminPermissions: $adminPermissions, tissueNcitIds: $tissueNcitIds, diseaseMeshIds: $diseaseMeshIds, visible: $visible, externalWebsite: $externalWebsite}}
  ) {
    clientMutationId
  }
}
    `;

/**
 * __useStudyUpdateMutation__
 *
 * To run a mutation, you first call `useStudyUpdateMutation` within a React component and pass it any options that fit your needs.
 * When your component renders, `useStudyUpdateMutation` returns a tuple that includes:
 * - A mutate function that you can call at any time to execute the mutation
 * - An object with fields that represent the current status of the mutation's execution
 *
 * @param baseOptions options that will be passed into the mutation, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options-2;
 *
 * @example
 * const [studyUpdateMutation, { data, loading, error }] = useStudyUpdateMutation({
 *   variables: {
 *      studyId: // value for 'studyId'
 *      studyName: // value for 'studyName'
 *      description: // value for 'description'
 *      readerPermissions: // value for 'readerPermissions'
 *      adminPermissions: // value for 'adminPermissions'
 *      tissueNcitIds: // value for 'tissueNcitIds'
 *      diseaseMeshIds: // value for 'diseaseMeshIds'
 *      visible: // value for 'visible'
 *      externalWebsite: // value for 'externalWebsite'
 *   },
 * });
 */
export function useStudyUpdateMutation(baseOptions?: Apollo.MutationHookOptions<StudyUpdateMutation, StudyUpdateMutationVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useMutation<StudyUpdateMutation, StudyUpdateMutationVariables>(StudyUpdateDocument, options);
      }
export type StudyUpdateMutationHookResult = ReturnType<typeof useStudyUpdateMutation>;
export type StudyUpdateMutationResult = Apollo.MutationResult<StudyUpdateMutation>;
export type StudyUpdateMutationOptions = Apollo.BaseMutationOptions<StudyUpdateMutation, StudyUpdateMutationVariables>;
export const StudyDeleteDocument = gql`
    mutation studyDelete($studyId: Int!) {
  deleteAllStudyData(input: {pStudyId: $studyId}) {
    boolean
  }
}
    `;

/**
 * __useStudyDeleteMutation__
 *
 * To run a mutation, you first call `useStudyDeleteMutation` within a React component and pass it any options that fit your needs.
 * When your component renders, `useStudyDeleteMutation` returns a tuple that includes:
 * - A mutate function that you can call at any time to execute the mutation
 * - An object with fields that represent the current status of the mutation's execution
 *
 * @param baseOptions options that will be passed into the mutation, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options-2;
 *
 * @example
 * const [studyDeleteMutation, { data, loading, error }] = useStudyDeleteMutation({
 *   variables: {
 *      studyId: // value for 'studyId'
 *   },
 * });
 */
export function useStudyDeleteMutation(baseOptions?: Apollo.MutationHookOptions<StudyDeleteMutation, StudyDeleteMutationVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useMutation<StudyDeleteMutation, StudyDeleteMutationVariables>(StudyDeleteDocument, options);
      }
export type StudyDeleteMutationHookResult = ReturnType<typeof useStudyDeleteMutation>;
export type StudyDeleteMutationResult = Apollo.MutationResult<StudyDeleteMutation>;
export type StudyDeleteMutationOptions = Apollo.BaseMutationOptions<StudyDeleteMutation, StudyDeleteMutationVariables>;
export const CreateStudyUploadDocument = gql`
    mutation createStudyUpload($filename: String!) {
  createStudyUpload(input: {filename: $filename}) {
    json
  }
}
    `;

/**
 * __useCreateStudyUploadMutation__
 *
 * To run a mutation, you first call `useCreateStudyUploadMutation` within a React component and pass it any options that fit your needs.
 * When your component renders, `useCreateStudyUploadMutation` returns a tuple that includes:
 * - A mutate function that you can call at any time to execute the mutation
 * - An object with fields that represent the current status of the mutation's execution
 *
 * @param baseOptions options that will be passed into the mutation, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options-2;
 *
 * @example
 * const [createStudyUploadMutation, { data, loading, error }] = useCreateStudyUploadMutation({
 *   variables: {
 *      filename: // value for 'filename'
 *   },
 * });
 */
export function useCreateStudyUploadMutation(baseOptions?: Apollo.MutationHookOptions<CreateStudyUploadMutation, CreateStudyUploadMutationVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useMutation<CreateStudyUploadMutation, CreateStudyUploadMutationVariables>(CreateStudyUploadDocument, options);
      }
export type CreateStudyUploadMutationHookResult = ReturnType<typeof useCreateStudyUploadMutation>;
export type CreateStudyUploadMutationResult = Apollo.MutationResult<CreateStudyUploadMutation>;
export type CreateStudyUploadMutationOptions = Apollo.BaseMutationOptions<CreateStudyUploadMutation, CreateStudyUploadMutationVariables>;
export const StudyDefinitionUpdateDocument = gql`
    mutation studyDefinitionUpdate {
  studyDefinitionUpdate(input: {}) {
    boolean
  }
}
    `;

/**
 * __useStudyDefinitionUpdateMutation__
 *
 * To run a mutation, you first call `useStudyDefinitionUpdateMutation` within a React component and pass it any options that fit your needs.
 * When your component renders, `useStudyDefinitionUpdateMutation` returns a tuple that includes:
 * - A mutate function that you can call at any time to execute the mutation
 * - An object with fields that represent the current status of the mutation's execution
 *
 * @param baseOptions options that will be passed into the mutation, supported options are listed on: https://www.apollographql.com/docs/react/api/react-hooks/#options-2;
 *
 * @example
 * const [studyDefinitionUpdateMutation, { data, loading, error }] = useStudyDefinitionUpdateMutation({
 *   variables: {
 *   },
 * });
 */
export function useStudyDefinitionUpdateMutation(baseOptions?: Apollo.MutationHookOptions<StudyDefinitionUpdateMutation, StudyDefinitionUpdateMutationVariables>) {
        const options = {...defaultOptions, ...baseOptions}
        return Apollo.useMutation<StudyDefinitionUpdateMutation, StudyDefinitionUpdateMutationVariables>(StudyDefinitionUpdateDocument, options);
      }
export type StudyDefinitionUpdateMutationHookResult = ReturnType<typeof useStudyDefinitionUpdateMutation>;
export type StudyDefinitionUpdateMutationResult = Apollo.MutationResult<StudyDefinitionUpdateMutation>;
export type StudyDefinitionUpdateMutationOptions = Apollo.BaseMutationOptions<StudyDefinitionUpdateMutation, StudyDefinitionUpdateMutationVariables>;