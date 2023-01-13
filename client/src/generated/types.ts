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
  Datetime: any;
  JSON: any;
};

export type AnnotationGroup = Node & {
  __typename?: 'AnnotationGroup';
  annotationGroupId: Scalars['Int'];
  /** Reads and enables pagination through a set of `AnnotationValue`. */
  annotationValuesList: Array<AnnotationValue>;
  displayGroup: Scalars['String'];
  h5AdColumn: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads and enables pagination through a set of `StudyAnnotationGroupUi`. */
  studyAnnotationGroupUisList: Array<StudyAnnotationGroupUi>;
};


export type AnnotationGroupAnnotationValuesListArgs = {
  condition: InputMaybe<AnnotationValueCondition>;
  filter: InputMaybe<AnnotationValueFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<AnnotationValuesOrderBy>>;
};


export type AnnotationGroupStudyAnnotationGroupUisListArgs = {
  condition: InputMaybe<StudyAnnotationGroupUiCondition>;
  filter: InputMaybe<StudyAnnotationGroupUiFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyAnnotationGroupUisOrderBy>>;
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
  color: Maybe<Scalars['String']>;
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressionsList: Array<DifferentialExpression>;
  displayValue: Scalars['String'];
  h5AdValue: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads and enables pagination through a set of `StudySampleAnnotation`. */
  studySampleAnnotationsList: Array<StudySampleAnnotation>;
};


export type AnnotationValueDifferentialExpressionsListArgs = {
  condition: InputMaybe<DifferentialExpressionCondition>;
  filter: InputMaybe<DifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
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
  /** Checks for equality with the object’s `color` field. */
  color: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `displayValue` field. */
  displayValue: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `h5AdValue` field. */
  h5AdValue: InputMaybe<Scalars['String']>;
};

/** A filter to be used against `AnnotationValue` object types. All fields are combined with a logical ‘and.’ */
export type AnnotationValueFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<AnnotationValueFilter>>;
  /** Filter by the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<IntFilter>;
  /** Filter by the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<IntFilter>;
  /** Filter by the object’s `color` field. */
  color: InputMaybe<StringFilter>;
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
  color: InputMaybe<Scalars['String']>;
  displayValue: Scalars['String'];
  h5AdValue: Scalars['String'];
};

/** Represents an update to a `AnnotationValue`. Fields that are set will be updated. */
export type AnnotationValuePatch = {
  annotationGroupId: InputMaybe<Scalars['Int']>;
  annotationValueId: InputMaybe<Scalars['Int']>;
  color: InputMaybe<Scalars['String']>;
  displayValue: InputMaybe<Scalars['String']>;
  h5AdValue: InputMaybe<Scalars['String']>;
};

/** Methods to use when ordering `AnnotationValue`. */
export enum AnnotationValuesOrderBy {
  AnnotationGroupIdAsc = 'ANNOTATION_GROUP_ID_ASC',
  AnnotationGroupIdDesc = 'ANNOTATION_GROUP_ID_DESC',
  AnnotationValueIdAsc = 'ANNOTATION_VALUE_ID_ASC',
  AnnotationValueIdDesc = 'ANNOTATION_VALUE_ID_DESC',
  ColorAsc = 'COLOR_ASC',
  ColorDesc = 'COLOR_DESC',
  DisplayValueAsc = 'DISPLAY_VALUE_ASC',
  DisplayValueDesc = 'DISPLAY_VALUE_DESC',
  H5AdValueAsc = 'H5AD_VALUE_ASC',
  H5AdValueDesc = 'H5AD_VALUE_DESC',
  Natural = 'NATURAL',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC'
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

/** An input for mutations affecting `BoxplotValue` */
export type BoxplotValueInput = {
  median: InputMaybe<Scalars['Float']>;
  n: InputMaybe<Scalars['Int']>;
  outliers: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  q1: InputMaybe<Scalars['Float']>;
  q1Whisker: InputMaybe<Scalars['Float']>;
  q3: InputMaybe<Scalars['Float']>;
  q3Whisker: InputMaybe<Scalars['Float']>;
};

export type Concept = Node & {
  __typename?: 'Concept';
  /** Reads and enables pagination through a set of `Concept`. */
  allChildrenList: Maybe<Array<Concept>>;
  /** Reads and enables pagination through a set of `Concept`. */
  allParentsList: Maybe<Array<Concept>>;
  allParentsPathsList: Maybe<Array<Maybe<ConceptAllParentsPathsRecord>>>;
  /** Reads and enables pagination through a set of `Concept`. */
  childrenDepthLimitList: Maybe<Array<Concept>>;
  childrenPathsDepthLimitList: Maybe<Array<Maybe<ConceptChildrenPathsDepthLimitRecord>>>;
  cid: Scalars['Int'];
  /** Reads and enables pagination through a set of `ConceptHierarchy`. */
  conceptHierarchiesByCidList: Array<ConceptHierarchy>;
  /** Reads and enables pagination through a set of `ConceptHierarchy`. */
  conceptHierarchiesByParentCidList: Array<ConceptHierarchy>;
  /** Reads and enables pagination through a set of `ConceptSynonym`. */
  conceptSynonymsByCidList: Array<ConceptSynonym>;
  /** Reads and enables pagination through a set of `Concept`. */
  directChildrenList: Maybe<Array<Concept>>;
  /** Reads and enables pagination through a set of `Concept`. */
  directParentsList: Maybe<Array<Concept>>;
  label: Maybe<Scalars['String']>;
  labelTsvector: Maybe<Scalars['String']>;
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  ontCode: Maybe<Scalars['String']>;
  ontid: Maybe<Scalars['Int']>;
  /** Reads a single `Ontology` that is related to this `Concept`. */
  ontologyByOntid: Maybe<Ontology>;
  /** Reads and enables pagination through a set of `Concept`. */
  testPotentialChildrenList: Maybe<Array<Concept>>;
};


export type ConceptAllChildrenListArgs = {
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptAllParentsListArgs = {
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptAllParentsPathsListArgs = {
  filter: InputMaybe<ConceptAllParentsPathsRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptChildrenDepthLimitListArgs = {
  depthLimit: InputMaybe<Scalars['Int']>;
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptChildrenPathsDepthLimitListArgs = {
  depthLimit: InputMaybe<Scalars['Int']>;
  filter: InputMaybe<ConceptChildrenPathsDepthLimitRecordFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptConceptHierarchiesByCidListArgs = {
  condition: InputMaybe<ConceptHierarchyCondition>;
  filter: InputMaybe<ConceptHierarchyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
};


export type ConceptConceptHierarchiesByParentCidListArgs = {
  condition: InputMaybe<ConceptHierarchyCondition>;
  filter: InputMaybe<ConceptHierarchyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
};


export type ConceptConceptSynonymsByCidListArgs = {
  condition: InputMaybe<ConceptSynonymCondition>;
  filter: InputMaybe<ConceptSynonymFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ConceptSynonymsOrderBy>>;
};


export type ConceptDirectChildrenListArgs = {
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
};


export type ConceptDirectParentsListArgs = {
  filter: InputMaybe<ConceptFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
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
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Reads a single `Ontology` that is related to this `Concept`. */
  ontologyByOntid: Maybe<Ontology>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Reads a single `OmicsBase` that is related to this `DifferentialExpression`. */
  omics: Maybe<OmicsBase>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `Study` that is related to this `DifferentialExpression`. */
  study: Maybe<Study>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Reads a single `OmicsProteinAntibodyTag` that is related to this `OmicsProteinAntibodyTagGene`. */
  proteinAntibodyTag: Maybe<OmicsProteinAntibodyTag>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Reads a single `OmicsBase` that is related to this `OmicsProteinAntibodyTag`. */
  proteinAntibodyTag: Maybe<OmicsBase>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
  /** Reads a single `OmicsTranscriptionFactor` that is related to this `OmicsTranscriptionFactorGene`. */
  transcriptionFactor: Maybe<OmicsTranscriptionFactor>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Reads a single `StudySample` that is related to this `StudySampleProjection`. */
  studyStudySample: Maybe<StudySample>;
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
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedAnnotationGroupNodeId: Maybe<Scalars['ID']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  deletedAnnotationValueNodeId: Maybe<Scalars['ID']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  deletedConceptNodeId: Maybe<Scalars['ID']>;
  /** Reads a single `Ontology` that is related to this `Concept`. */
  ontologyByOntid: Maybe<Ontology>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Reads a single `OmicsBase` that is related to this `OmicsProteinAntibodyTag`. */
  proteinAntibodyTag: Maybe<OmicsBase>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  log2Foldchange: Maybe<Scalars['Float']>;
  omicsId: Maybe<Scalars['Int']>;
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
  /** Filter by the object’s `log2Foldchange` field. */
  log2Foldchange: InputMaybe<FloatFilter>;
  /** Negates the expression. */
  not: InputMaybe<DifferentialExpressionVFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
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

/** Methods to use when ordering `DifferentialExpressionV`. */
export enum DifferentialExpressionVsOrderBy {
  AnnotationValueIdAsc = 'ANNOTATION_VALUE_ID_ASC',
  AnnotationValueIdDesc = 'ANNOTATION_VALUE_ID_DESC',
  DisplayNameAsc = 'DISPLAY_NAME_ASC',
  DisplayNameDesc = 'DISPLAY_NAME_DESC',
  DisplaySymbolAsc = 'DISPLAY_SYMBOL_ASC',
  DisplaySymbolDesc = 'DISPLAY_SYMBOL_DESC',
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

export type ExpressionByAnnotationBoxplot = {
  __typename?: 'ExpressionByAnnotationBoxplot';
  annotationGroupId: Maybe<Scalars['Int']>;
  annotationValueId: Maybe<Scalars['Int']>;
  boxplot: Maybe<BoxplotValue>;
  omicsId: Maybe<Scalars['Int']>;
  studyLayerId: Maybe<Scalars['Int']>;
  values: Maybe<Array<Maybe<Scalars['Float']>>>;
};

/**
 * A condition to be used against `ExpressionByAnnotationBoxplot` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type ExpressionByAnnotationBoxplotCondition = {
  /** Checks for equality with the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `boxplot` field. */
  boxplot: InputMaybe<BoxplotValueInput>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `values` field. */
  values: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A filter to be used against `ExpressionByAnnotationBoxplot` object types. All fields are combined with a logical ‘and.’ */
export type ExpressionByAnnotationBoxplotFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<ExpressionByAnnotationBoxplotFilter>>;
  /** Filter by the object’s `annotationGroupId` field. */
  annotationGroupId: InputMaybe<IntFilter>;
  /** Filter by the object’s `annotationValueId` field. */
  annotationValueId: InputMaybe<IntFilter>;
  /** Filter by the object’s `boxplot` field. */
  boxplot: InputMaybe<BoxplotValueFilter>;
  /** Negates the expression. */
  not: InputMaybe<ExpressionByAnnotationBoxplotFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<ExpressionByAnnotationBoxplotFilter>>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId: InputMaybe<IntFilter>;
  /** Filter by the object’s `values` field. */
  values: InputMaybe<FloatListFilter>;
};

/** Methods to use when ordering `ExpressionByAnnotationBoxplot`. */
export enum ExpressionByAnnotationBoxplotsOrderBy {
  AnnotationGroupIdAsc = 'ANNOTATION_GROUP_ID_ASC',
  AnnotationGroupIdDesc = 'ANNOTATION_GROUP_ID_DESC',
  AnnotationValueIdAsc = 'ANNOTATION_VALUE_ID_ASC',
  AnnotationValueIdDesc = 'ANNOTATION_VALUE_ID_DESC',
  BoxplotAsc = 'BOXPLOT_ASC',
  BoxplotDesc = 'BOXPLOT_DESC',
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  StudyLayerIdAsc = 'STUDY_LAYER_ID_ASC',
  StudyLayerIdDesc = 'STUDY_LAYER_ID_DESC',
  ValuesAsc = 'VALUES_ASC',
  ValuesDesc = 'VALUES_DESC'
}

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
  /** Creates a single `OmicsTranscriptionFactor`. */
  createOmicsTranscriptionFactor: Maybe<CreateOmicsTranscriptionFactorPayload>;
  /** Creates a single `OmicsTranscriptionFactorGene`. */
  createOmicsTranscriptionFactorGene: Maybe<CreateOmicsTranscriptionFactorGenePayload>;
  /** Creates a single `Ontology`. */
  createOntology: Maybe<CreateOntologyPayload>;
  /** Creates a single `Study`. */
  createStudy: Maybe<CreateStudyPayload>;
  /** Creates a single `StudyAnnotationGroupUi`. */
  createStudyAnnotationGroupUi: Maybe<CreateStudyAnnotationGroupUiPayload>;
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
export type MutationCreateStudyArgs = {
  input: CreateStudyInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudyAnnotationGroupUiArgs = {
  input: CreateStudyAnnotationGroupUiInput;
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
  /** Filter by the object’s `taxId` field. */
  taxId: InputMaybe<IntFilter>;
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
  TaxIdAsc = 'TAX_ID_ASC',
  TaxIdDesc = 'TAX_ID_DESC'
}

export type OmicsBase = Node & {
  __typename?: 'OmicsBase';
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressionsByOmicsIdList: Array<DifferentialExpression>;
  displayName: Maybe<Scalars['String']>;
  displaySymbol: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads a single `OmicsGene` that is related to this `OmicsBase`. */
  omicsGeneByGeneId: Maybe<OmicsGene>;
  omicsId: Scalars['Int'];
  /** Reads a single `OmicsProteinAntibodyTag` that is related to this `OmicsBase`. */
  omicsProteinAntibodyTagByProteinAntibodyTagId: Maybe<OmicsProteinAntibodyTag>;
  /** Reads a single `OmicsTranscriptionFactor` that is related to this `OmicsBase`. */
  omicsTranscriptionFactorByOmicsId: Maybe<OmicsTranscriptionFactor>;
  omicsType: OmicsType;
  /** Reads and enables pagination through a set of `StudyOmic`. */
  studyOmicsByOmicsIdList: Array<StudyOmic>;
  taxId: Scalars['Int'];
};


export type OmicsBaseDifferentialExpressionsByOmicsIdListArgs = {
  condition: InputMaybe<DifferentialExpressionCondition>;
  filter: InputMaybe<DifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
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
  /** Checks for equality with the object’s `displaySymbol` field. */
  displaySymbol: InputMaybe<Scalars['String']>;
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
  displaySymbol: Scalars['String'];
  omicsId: InputMaybe<Scalars['Int']>;
  omicsType: OmicsType;
  taxId: Scalars['Int'];
};

/** Represents an update to a `OmicsBase`. Fields that are set will be updated. */
export type OmicsBasePatch = {
  displayName: InputMaybe<Scalars['String']>;
  displaySymbol: InputMaybe<Scalars['String']>;
  omicsId: InputMaybe<Scalars['Int']>;
  omicsType: InputMaybe<OmicsType>;
  taxId: InputMaybe<Scalars['Int']>;
};

/** Methods to use when ordering `OmicsBase`. */
export enum OmicsBasesOrderBy {
  DisplayNameAsc = 'DISPLAY_NAME_ASC',
  DisplayNameDesc = 'DISPLAY_NAME_DESC',
  DisplaySymbolAsc = 'DISPLAY_SYMBOL_ASC',
  DisplaySymbolDesc = 'DISPLAY_SYMBOL_DESC',
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
  omicsProteinAntibodyTagGenesByGeneIdList: Array<OmicsProteinAntibodyTagGene>;
  /** Reads and enables pagination through a set of `OmicsTranscriptionFactorGene`. */
  omicsTranscriptionFactorGenesByGeneIdList: Array<OmicsTranscriptionFactorGene>;
};


export type OmicsGeneOmicsProteinAntibodyTagGenesByGeneIdListArgs = {
  condition: InputMaybe<OmicsProteinAntibodyTagGeneCondition>;
  filter: InputMaybe<OmicsProteinAntibodyTagGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsProteinAntibodyTagGenesOrderBy>>;
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
  omicsProteinAntibodyTagGenesByProteinAntibodyTagIdList: Array<OmicsProteinAntibodyTagGene>;
  /** Reads a single `OmicsBase` that is related to this `OmicsProteinAntibodyTag`. */
  proteinAntibodyTag: Maybe<OmicsBase>;
  proteinAntibodyTagId: Scalars['Int'];
  taxId: Scalars['Int'];
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

export type OmicsTranscriptionFactor = Node & {
  __typename?: 'OmicsTranscriptionFactor';
  jasparMatrixId: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads a single `OmicsBase` that is related to this `OmicsTranscriptionFactor`. */
  omics: Maybe<OmicsBase>;
  omicsId: Scalars['Int'];
  /** Reads and enables pagination through a set of `OmicsTranscriptionFactorGene`. */
  omicsTranscriptionFactorGenesByTranscriptionFactorIdList: Array<OmicsTranscriptionFactorGene>;
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
  conceptsByOntidList: Array<Concept>;
  name: Maybe<Scalars['String']>;
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  ontid: Scalars['Int'];
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

export enum ProjectionType {
  Pca = 'PCA',
  Tsne = 'TSNE',
  Umap = 'UMAP'
}

/** A filter to be used against ProjectionType fields. All fields are combined with a logical ‘and.’ */
export type ProjectionTypeFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom: InputMaybe<ProjectionType>;
  /** Equal to the specified value. */
  equalTo: InputMaybe<ProjectionType>;
  /** Greater than the specified value. */
  greaterThan: InputMaybe<ProjectionType>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo: InputMaybe<ProjectionType>;
  /** Included in the specified list. */
  in: InputMaybe<Array<ProjectionType>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan: InputMaybe<ProjectionType>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo: InputMaybe<ProjectionType>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom: InputMaybe<ProjectionType>;
  /** Not equal to the specified value. */
  notEqualTo: InputMaybe<ProjectionType>;
  /** Not included in the specified list. */
  notIn: InputMaybe<Array<ProjectionType>>;
};

/** The root query type which gives access points into the data universe. */
export type Query = Node & {
  __typename?: 'Query';
  /** Reads a set of `_AllUsedOntologyId`. */
  _allUsedOntologyIdsList: Maybe<Array<_AllUsedOntologyId>>;
  _conceptHierarchyMinimumTreesImpl: Maybe<MinimumTreesResult>;
  _finalBoxplot: Maybe<BoxplotValue>;
  _semanticOrderImpl: Maybe<Array<Maybe<Scalars['Int']>>>;
  annotationGroup: Maybe<AnnotationGroup>;
  /** Reads a single `AnnotationGroup` using its globally unique `ID`. */
  annotationGroupByNodeId: Maybe<AnnotationGroup>;
  /** Reads a set of `AnnotationGroup`. */
  annotationGroupsList: Maybe<Array<AnnotationGroup>>;
  annotationValue: Maybe<AnnotationValue>;
  /** Reads a single `AnnotationValue` using its globally unique `ID`. */
  annotationValueByNodeId: Maybe<AnnotationValue>;
  /** Reads a set of `AnnotationValue`. */
  annotationValuesList: Maybe<Array<AnnotationValue>>;
  /** Reads and enables pagination through a set of `AutocompleteResult`. */
  autocompleteList: Maybe<Array<AutocompleteResult>>;
  concept: Maybe<Concept>;
  /** Reads a single `Concept` using its globally unique `ID`. */
  conceptByNodeId: Maybe<Concept>;
  conceptCidArrayToCodes: Maybe<Array<Maybe<Scalars['String']>>>;
  /** Reads a set of `ConceptHierarchy`. */
  conceptHierarchiesList: Maybe<Array<ConceptHierarchy>>;
  conceptHierarchyMinimumTreesParentsListsList: Maybe<Array<Maybe<ConceptHierarchyMinimumTreesParentsListsRecord>>>;
  /** Reads a set of `ConceptSynonym`. */
  conceptSynonymsList: Maybe<Array<ConceptSynonym>>;
  /** Reads and enables pagination through a set of `Concept`. */
  conceptsInSemanticOrderList: Maybe<Array<Concept>>;
  /** Reads a set of `Concept`. */
  conceptsList: Maybe<Array<Concept>>;
  correlationTrianglePlot: Maybe<Scalars['String']>;
  /** Reads a set of `DifferentialExpressionV`. */
  differentialExpressionVsList: Maybe<Array<DifferentialExpressionV>>;
  /** Reads a set of `DifferentialExpression`. */
  differentialExpressionsList: Maybe<Array<DifferentialExpression>>;
  /** Reads a set of `ExpressionByAnnotationBoxplot`. */
  expressionByAnnotationBoxplotsList: Maybe<Array<ExpressionByAnnotationBoxplot>>;
  /** Reads and enables pagination through a set of `ExpressionByOmic`. */
  expressionByOmicsIdsList: Maybe<Array<ExpressionByOmic>>;
  /** Fetches an object given its globally unique `ID`. */
  node: Maybe<Node>;
  /** The root query type must be a `Node` to work well with Relay 1 mutations. This just resolves to `query`. */
  nodeId: Scalars['ID'];
  /** Reads a set of `OmicsAll`. */
  omicsAllsList: Maybe<Array<OmicsAll>>;
  omicsBase: Maybe<OmicsBase>;
  /** Reads a single `OmicsBase` using its globally unique `ID`. */
  omicsBaseByNodeId: Maybe<OmicsBase>;
  /** Reads a set of `OmicsBase`. */
  omicsBasesList: Maybe<Array<OmicsBase>>;
  omicsGene: Maybe<OmicsGene>;
  /** Reads a single `OmicsGene` using its globally unique `ID`. */
  omicsGeneByNodeId: Maybe<OmicsGene>;
  /** Reads a set of `OmicsGene`. */
  omicsGenesList: Maybe<Array<OmicsGene>>;
  omicsProteinAntibodyTag: Maybe<OmicsProteinAntibodyTag>;
  /** Reads a single `OmicsProteinAntibodyTag` using its globally unique `ID`. */
  omicsProteinAntibodyTagByNodeId: Maybe<OmicsProteinAntibodyTag>;
  /** Reads a set of `OmicsProteinAntibodyTagGene`. */
  omicsProteinAntibodyTagGenesList: Maybe<Array<OmicsProteinAntibodyTagGene>>;
  /** Reads a set of `OmicsProteinAntibodyTag`. */
  omicsProteinAntibodyTagsList: Maybe<Array<OmicsProteinAntibodyTag>>;
  omicsTranscriptionFactor: Maybe<OmicsTranscriptionFactor>;
  /** Reads a single `OmicsTranscriptionFactor` using its globally unique `ID`. */
  omicsTranscriptionFactorByNodeId: Maybe<OmicsTranscriptionFactor>;
  /** Reads a set of `OmicsTranscriptionFactorGene`. */
  omicsTranscriptionFactorGenesList: Maybe<Array<OmicsTranscriptionFactorGene>>;
  /** Reads a set of `OmicsTranscriptionFactor`. */
  omicsTranscriptionFactorsList: Maybe<Array<OmicsTranscriptionFactor>>;
  ontCodesInfoList: Maybe<Array<Maybe<OntCodesInfoRecord>>>;
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
  /** Reads a set of `Study`. */
  studiesList: Maybe<Array<Study>>;
  study: Maybe<Study>;
  /** Reads a set of `StudyAnnotationGroupUi`. */
  studyAnnotationGroupUisList: Maybe<Array<StudyAnnotationGroupUi>>;
  /** Reads a single `Study` using its globally unique `ID`. */
  studyByNodeId: Maybe<Study>;
  studyLayer: Maybe<StudyLayer>;
  /** Reads a single `StudyLayer` using its globally unique `ID`. */
  studyLayerByNodeId: Maybe<StudyLayer>;
  /** Reads a set of `StudyLayer`. */
  studyLayersList: Maybe<Array<StudyLayer>>;
  /** Reads a set of `StudyOmic`. */
  studyOmicsList: Maybe<Array<StudyOmic>>;
  /** Reads a set of `StudyOmicsTransposed`. */
  studyOmicsTransposedsList: Maybe<Array<StudyOmicsTransposed>>;
  /** Reads a set of `StudyOverviewOntology`. */
  studyOverviewOntologiesList: Maybe<Array<StudyOverviewOntology>>;
  /** Reads a set of `StudyOverview`. */
  studyOverviewsList: Maybe<Array<StudyOverview>>;
  studySample: Maybe<StudySample>;
  /** Reads a set of `StudySampleAnnotationSubsampling`. */
  studySampleAnnotationSubsamplingsList: Maybe<Array<StudySampleAnnotationSubsampling>>;
  /** Reads a set of `StudySampleAnnotation`. */
  studySampleAnnotationsList: Maybe<Array<StudySampleAnnotation>>;
  /** Reads a single `StudySample` using its globally unique `ID`. */
  studySampleByNodeId: Maybe<StudySample>;
  /** Reads a set of `StudySampleProjectionSubsamplingTransposed`. */
  studySampleProjectionSubsamplingTransposedsList: Maybe<Array<StudySampleProjectionSubsamplingTransposed>>;
  /** Reads a set of `StudySampleProjection`. */
  studySampleProjectionsList: Maybe<Array<StudySampleProjection>>;
  /** Reads a set of `StudySample`. */
  studySamplesList: Maybe<Array<StudySample>>;
  /** Reads a set of `TreeOntology`. */
  treeOntologiesList: Maybe<Array<TreeOntology>>;
  violinPlot: Maybe<Scalars['String']>;
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
export type QueryAnnotationValuesListArgs = {
  condition: InputMaybe<AnnotationValueCondition>;
  filter: InputMaybe<AnnotationValueFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<AnnotationValuesOrderBy>>;
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
export type QueryConceptHierarchiesListArgs = {
  condition: InputMaybe<ConceptHierarchyCondition>;
  filter: InputMaybe<ConceptHierarchyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
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
export type QueryConceptSynonymsListArgs = {
  condition: InputMaybe<ConceptSynonymCondition>;
  filter: InputMaybe<ConceptSynonymFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ConceptSynonymsOrderBy>>;
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
  pOmicsIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pStudyLayerId: InputMaybe<Scalars['Int']>;
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
export type QueryDifferentialExpressionsListArgs = {
  condition: InputMaybe<DifferentialExpressionCondition>;
  filter: InputMaybe<DifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpressionByAnnotationBoxplotsListArgs = {
  condition: InputMaybe<ExpressionByAnnotationBoxplotCondition>;
  filter: InputMaybe<ExpressionByAnnotationBoxplotFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<ExpressionByAnnotationBoxplotsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpressionByOmicsIdsListArgs = {
  filter: InputMaybe<ExpressionByOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  pOmicsIds: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  pStudyLayerId: InputMaybe<Scalars['Int']>;
  pSubsamplingProjection: InputMaybe<ProjectionType>;
};


/** The root query type which gives access points into the data universe. */
export type QueryNodeArgs = {
  nodeId: Scalars['ID'];
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
export type QueryOmicsBaseArgs = {
  omicsId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsBaseByNodeIdArgs = {
  nodeId: Scalars['ID'];
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
export type QueryOmicsProteinAntibodyTagGenesListArgs = {
  condition: InputMaybe<OmicsProteinAntibodyTagGeneCondition>;
  filter: InputMaybe<OmicsProteinAntibodyTagGeneFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsProteinAntibodyTagGenesOrderBy>>;
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
export type QueryOmicsTranscriptionFactorArgs = {
  omicsId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsTranscriptionFactorByNodeIdArgs = {
  nodeId: Scalars['ID'];
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
export type QueryOmicsTranscriptionFactorsListArgs = {
  condition: InputMaybe<OmicsTranscriptionFactorCondition>;
  filter: InputMaybe<OmicsTranscriptionFactorFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<OmicsTranscriptionFactorsOrderBy>>;
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
export type QueryStudyLayerArgs = {
  studyLayerId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyLayerByNodeIdArgs = {
  nodeId: Scalars['ID'];
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
export type QueryStudyOmicsListArgs = {
  condition: InputMaybe<StudyOmicCondition>;
  filter: InputMaybe<StudyOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyOmicsOrderBy>>;
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
export type QueryStudyOverviewOntologiesListArgs = {
  condition: InputMaybe<StudyOverviewOntologyCondition>;
  filter: InputMaybe<StudyOverviewOntologyFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyOverviewOntologiesOrderBy>>;
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
export type QueryStudySampleAnnotationSubsamplingsListArgs = {
  condition: InputMaybe<StudySampleAnnotationSubsamplingCondition>;
  filter: InputMaybe<StudySampleAnnotationSubsamplingFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySampleAnnotationSubsamplingsOrderBy>>;
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
export type QueryStudySampleProjectionSubsamplingTransposedsListArgs = {
  condition: InputMaybe<StudySampleProjectionSubsamplingTransposedCondition>;
  filter: InputMaybe<StudySampleProjectionSubsamplingTransposedFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySampleProjectionSubsamplingTransposedsOrderBy>>;
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
export type QueryStudySamplesListArgs = {
  condition: InputMaybe<StudySampleCondition>;
  filter: InputMaybe<StudySampleFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySamplesOrderBy>>;
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
export type QueryViolinPlotArgs = {
  pAnnotationGroupId: InputMaybe<Scalars['Int']>;
  pOmicsId: InputMaybe<Scalars['Int']>;
  pStudyId: InputMaybe<Scalars['Int']>;
  pStudyLayerId: InputMaybe<Scalars['Int']>;
};

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

/** Methods to use when ordering `Study`. */
export enum StudiesOrderBy {
  AttributeValueFreqAsc = 'ATTRIBUTE_VALUE_FREQ_ASC',
  AttributeValueFreqDesc = 'ATTRIBUTE_VALUE_FREQ_DESC',
  CellCountAsc = 'CELL_COUNT_ASC',
  CellCountDesc = 'CELL_COUNT_DESC',
  CellOntologyIdsAsc = 'CELL_ONTOLOGY_IDS_ASC',
  CellOntologyIdsDesc = 'CELL_ONTOLOGY_IDS_DESC',
  ClusterColorMapAsc = 'CLUSTER_COLOR_MAP_ASC',
  ClusterColorMapDesc = 'CLUSTER_COLOR_MAP_DESC',
  ClusterHullsAsc = 'CLUSTER_HULLS_ASC',
  ClusterHullsDesc = 'CLUSTER_HULLS_DESC',
  DescriptionAsc = 'DESCRIPTION_ASC',
  DescriptionDesc = 'DESCRIPTION_DESC',
  DiseaseMeshIdsAsc = 'DISEASE_MESH_IDS_ASC',
  DiseaseMeshIdsDesc = 'DISEASE_MESH_IDS_DESC',
  H5AdfileModifiedDateAsc = 'H5ADFILE_MODIFIED_DATE_ASC',
  H5AdfileModifiedDateDesc = 'H5ADFILE_MODIFIED_DATE_DESC',
  ImportStatusAsc = 'IMPORT_STATUS_ASC',
  ImportStatusDesc = 'IMPORT_STATUS_DESC',
  ImportStatusUpdatedAsc = 'IMPORT_STATUS_UPDATED_ASC',
  ImportStatusUpdatedDesc = 'IMPORT_STATUS_UPDATED_DESC',
  Natural = 'NATURAL',
  OrganismTaxIdAsc = 'ORGANISM_TAX_ID_ASC',
  OrganismTaxIdDesc = 'ORGANISM_TAX_ID_DESC',
  PlotCoordsAsc = 'PLOT_COORDS_ASC',
  PlotCoordsDesc = 'PLOT_COORDS_DESC',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  ProjectionCellCoordsAsc = 'PROJECTION_CELL_COORDS_ASC',
  ProjectionCellCoordsDesc = 'PROJECTION_CELL_COORDS_DESC',
  ProjectionCellIndicesAsc = 'PROJECTION_CELL_INDICES_ASC',
  ProjectionCellIndicesDesc = 'PROJECTION_CELL_INDICES_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudyNameAsc = 'STUDY_NAME_ASC',
  StudyNameDesc = 'STUDY_NAME_DESC',
  TissueNcitIdsAsc = 'TISSUE_NCIT_IDS_ASC',
  TissueNcitIdsDesc = 'TISSUE_NCIT_IDS_DESC'
}

export type Study = Node & {
  __typename?: 'Study';
  attributeValueFreq: Maybe<Scalars['JSON']>;
  cellCount: Maybe<Scalars['Int']>;
  cellOntologyIds: Maybe<Array<Maybe<Scalars['String']>>>;
  clusterColorMap: Maybe<Scalars['JSON']>;
  clusterHulls: Maybe<Scalars['JSON']>;
  description: Maybe<Scalars['String']>;
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressionsList: Array<DifferentialExpression>;
  diseaseMeshIds: Maybe<Array<Maybe<Scalars['String']>>>;
  h5AdfileModifiedDate: Maybe<Scalars['Datetime']>;
  importStatus: Maybe<Scalars['String']>;
  importStatusUpdated: Maybe<Scalars['Datetime']>;
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  organismTaxId: Maybe<Scalars['String']>;
  plotCoords: Maybe<Scalars['JSON']>;
  projectionCellCoords: Maybe<Scalars['JSON']>;
  projectionCellIndices: Maybe<Scalars['JSON']>;
  /** Reads and enables pagination through a set of `StudyAnnotationGroupUi`. */
  studyAnnotationGroupUisList: Array<StudyAnnotationGroupUi>;
  studyId: Scalars['Int'];
  /** Reads and enables pagination through a set of `StudyLayer`. */
  studyLayersList: Array<StudyLayer>;
  studyName: Scalars['String'];
  /** Reads and enables pagination through a set of `StudyOmic`. */
  studyOmicsList: Array<StudyOmic>;
  /** Reads and enables pagination through a set of `StudyOmicsTransposed`. */
  studyOmicsTransposedList: Array<StudyOmicsTransposed>;
  /** Reads and enables pagination through a set of `StudySampleAnnotationSubsampling`. */
  studySampleAnnotationSubsamplingList: Array<StudySampleAnnotationSubsampling>;
  /** Reads and enables pagination through a set of `StudySampleAnnotation`. */
  studySampleAnnotationsList: Array<StudySampleAnnotation>;
  /** Reads and enables pagination through a set of `StudySampleProjectionSubsamplingTransposed`. */
  studySampleProjectionSubsamplingTransposedList: Array<StudySampleProjectionSubsamplingTransposed>;
  /** Reads and enables pagination through a set of `StudySample`. */
  studySamplesList: Array<StudySample>;
  tissueNcitIds: Maybe<Array<Maybe<Scalars['String']>>>;
};


export type StudyDifferentialExpressionsListArgs = {
  condition: InputMaybe<DifferentialExpressionCondition>;
  filter: InputMaybe<DifferentialExpressionFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


export type StudyStudyAnnotationGroupUisListArgs = {
  condition: InputMaybe<StudyAnnotationGroupUiCondition>;
  filter: InputMaybe<StudyAnnotationGroupUiFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyAnnotationGroupUisOrderBy>>;
};


export type StudyStudyLayersListArgs = {
  condition: InputMaybe<StudyLayerCondition>;
  filter: InputMaybe<StudyLayerFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyLayersOrderBy>>;
};


export type StudyStudyOmicsListArgs = {
  condition: InputMaybe<StudyOmicCondition>;
  filter: InputMaybe<StudyOmicFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyOmicsOrderBy>>;
};


export type StudyStudyOmicsTransposedListArgs = {
  condition: InputMaybe<StudyOmicsTransposedCondition>;
  filter: InputMaybe<StudyOmicsTransposedFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudyOmicsTransposedsOrderBy>>;
};


export type StudyStudySampleAnnotationSubsamplingListArgs = {
  condition: InputMaybe<StudySampleAnnotationSubsamplingCondition>;
  filter: InputMaybe<StudySampleAnnotationSubsamplingFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySampleAnnotationSubsamplingsOrderBy>>;
};


export type StudyStudySampleAnnotationsListArgs = {
  condition: InputMaybe<StudySampleAnnotationCondition>;
  filter: InputMaybe<StudySampleAnnotationFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySampleAnnotationsOrderBy>>;
};


export type StudyStudySampleProjectionSubsamplingTransposedListArgs = {
  condition: InputMaybe<StudySampleProjectionSubsamplingTransposedCondition>;
  filter: InputMaybe<StudySampleProjectionSubsamplingTransposedFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySampleProjectionSubsamplingTransposedsOrderBy>>;
};


export type StudyStudySamplesListArgs = {
  condition: InputMaybe<StudySampleCondition>;
  filter: InputMaybe<StudySampleFilter>;
  first: InputMaybe<Scalars['Int']>;
  offset: InputMaybe<Scalars['Int']>;
  orderBy: InputMaybe<Array<StudySamplesOrderBy>>;
};

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
  /** Checks for equality with the object’s `attributeValueFreq` field. */
  attributeValueFreq: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `cellCount` field. */
  cellCount: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `cellOntologyIds` field. */
  cellOntologyIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `clusterColorMap` field. */
  clusterColorMap: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `clusterHulls` field. */
  clusterHulls: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `description` field. */
  description: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `diseaseMeshIds` field. */
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `h5AdfileModifiedDate` field. */
  h5AdfileModifiedDate: InputMaybe<Scalars['Datetime']>;
  /** Checks for equality with the object’s `importStatus` field. */
  importStatus: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `importStatusUpdated` field. */
  importStatusUpdated: InputMaybe<Scalars['Datetime']>;
  /** Checks for equality with the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `plotCoords` field. */
  plotCoords: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `projectionCellCoords` field. */
  projectionCellCoords: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `projectionCellIndices` field. */
  projectionCellIndices: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyName` field. */
  studyName: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `tissueNcitIds` field. */
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** A filter to be used against `Study` object types. All fields are combined with a logical ‘and.’ */
export type StudyFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudyFilter>>;
  /** Filter by the object’s `attributeValueFreq` field. */
  attributeValueFreq: InputMaybe<JsonFilter>;
  /** Filter by the object’s `cellCount` field. */
  cellCount: InputMaybe<IntFilter>;
  /** Filter by the object’s `cellOntologyIds` field. */
  cellOntologyIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `clusterColorMap` field. */
  clusterColorMap: InputMaybe<JsonFilter>;
  /** Filter by the object’s `clusterHulls` field. */
  clusterHulls: InputMaybe<JsonFilter>;
  /** Filter by the object’s `description` field. */
  description: InputMaybe<StringFilter>;
  /** Filter by the object’s `diseaseMeshIds` field. */
  diseaseMeshIds: InputMaybe<StringListFilter>;
  /** Filter by the object’s `h5AdfileModifiedDate` field. */
  h5AdfileModifiedDate: InputMaybe<DatetimeFilter>;
  /** Filter by the object’s `importStatus` field. */
  importStatus: InputMaybe<StringFilter>;
  /** Filter by the object’s `importStatusUpdated` field. */
  importStatusUpdated: InputMaybe<DatetimeFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyFilter>>;
  /** Filter by the object’s `organismTaxId` field. */
  organismTaxId: InputMaybe<StringFilter>;
  /** Filter by the object’s `plotCoords` field. */
  plotCoords: InputMaybe<JsonFilter>;
  /** Filter by the object’s `projectionCellCoords` field. */
  projectionCellCoords: InputMaybe<JsonFilter>;
  /** Filter by the object’s `projectionCellIndices` field. */
  projectionCellIndices: InputMaybe<JsonFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyName` field. */
  studyName: InputMaybe<StringFilter>;
  /** Filter by the object’s `tissueNcitIds` field. */
  tissueNcitIds: InputMaybe<StringListFilter>;
};

/** An input for mutations affecting `Study` */
export type StudyInput = {
  attributeValueFreq: InputMaybe<Scalars['JSON']>;
  cellCount: InputMaybe<Scalars['Int']>;
  cellOntologyIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  clusterColorMap: InputMaybe<Scalars['JSON']>;
  clusterHulls: InputMaybe<Scalars['JSON']>;
  description: InputMaybe<Scalars['String']>;
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  h5AdfileModifiedDate: InputMaybe<Scalars['Datetime']>;
  importStatus: InputMaybe<Scalars['String']>;
  importStatusUpdated: InputMaybe<Scalars['Datetime']>;
  organismTaxId: InputMaybe<Scalars['String']>;
  plotCoords: InputMaybe<Scalars['JSON']>;
  projectionCellCoords: InputMaybe<Scalars['JSON']>;
  projectionCellIndices: InputMaybe<Scalars['JSON']>;
  studyId: InputMaybe<Scalars['Int']>;
  studyName: Scalars['String'];
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

export type StudyLayer = Node & {
  __typename?: 'StudyLayer';
  layer: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  omicsType: OmicsType;
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
  omicsType: OmicsType;
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
  /** Reads a single `Study` that is related to this `StudyOmicsTransposed`. */
  study: Maybe<Study>;
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
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
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
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

export type StudyOverview = {
  __typename?: 'StudyOverview';
  cellCount: Maybe<Scalars['Int']>;
  description: Maybe<Scalars['String']>;
  studyId: Maybe<Scalars['Int']>;
  studyName: Maybe<Scalars['String']>;
  /** Reads and enables pagination through a set of `StudyOverviewOntology`. */
  studyOntologyList: Array<StudyOverviewOntology>;
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
  /** Checks for equality with the object’s `description` field. */
  description: InputMaybe<Scalars['String']>;
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
  /** Filter by the object’s `description` field. */
  description: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not: InputMaybe<StudyOverviewFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudyOverviewFilter>>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyName` field. */
  studyName: InputMaybe<StringFilter>;
};

/** An input for mutations affecting `StudyOverview` */
export type StudyOverviewInput = {
  cellCount: InputMaybe<Scalars['Int']>;
  description: InputMaybe<Scalars['String']>;
  studyId: InputMaybe<Scalars['Int']>;
  studyName: InputMaybe<Scalars['String']>;
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

/** Methods to use when ordering `StudyOverview`. */
export enum StudyOverviewsOrderBy {
  CellCountAsc = 'CELL_COUNT_ASC',
  CellCountDesc = 'CELL_COUNT_DESC',
  DescriptionAsc = 'DESCRIPTION_ASC',
  DescriptionDesc = 'DESCRIPTION_DESC',
  Natural = 'NATURAL',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudyNameAsc = 'STUDY_NAME_ASC',
  StudyNameDesc = 'STUDY_NAME_DESC'
}

/** Represents an update to a `Study`. Fields that are set will be updated. */
export type StudyPatch = {
  attributeValueFreq: InputMaybe<Scalars['JSON']>;
  cellCount: InputMaybe<Scalars['Int']>;
  cellOntologyIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  clusterColorMap: InputMaybe<Scalars['JSON']>;
  clusterHulls: InputMaybe<Scalars['JSON']>;
  description: InputMaybe<Scalars['String']>;
  diseaseMeshIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  h5AdfileModifiedDate: InputMaybe<Scalars['Datetime']>;
  importStatus: InputMaybe<Scalars['String']>;
  importStatusUpdated: InputMaybe<Scalars['Datetime']>;
  organismTaxId: InputMaybe<Scalars['String']>;
  plotCoords: InputMaybe<Scalars['JSON']>;
  projectionCellCoords: InputMaybe<Scalars['JSON']>;
  projectionCellIndices: InputMaybe<Scalars['JSON']>;
  studyId: InputMaybe<Scalars['Int']>;
  studyName: InputMaybe<Scalars['String']>;
  tissueNcitIds: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

export type StudySample = Node & {
  __typename?: 'StudySample';
  h5AdObsIndex: Scalars['Int'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads a single `Study` that is related to this `StudySample`. */
  study: Maybe<Study>;
  studyId: Scalars['Int'];
  studySampleId: Scalars['Int'];
  /** Reads and enables pagination through a set of `StudySampleProjection`. */
  studySampleProjectionsByStudyIdAndStudySampleIdList: Array<StudySampleProjection>;
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
  studyId: Scalars['Int'];
  studySampleIds: Array<InputMaybe<Scalars['Int']>>;
};

export type StudySampleAnnotationSubsampling = {
  __typename?: 'StudySampleAnnotationSubsampling';
  annotationValueId: Maybe<Scalars['Int']>;
  /** Reads a single `Study` that is related to this `StudySampleAnnotationSubsampling`. */
  study: Maybe<Study>;
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
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds: InputMaybe<IntListFilter>;
};

/** Methods to use when ordering `StudySampleAnnotationSubsampling`. */
export enum StudySampleAnnotationSubsamplingsOrderBy {
  AnnotationValueIdAsc = 'ANNOTATION_VALUE_ID_ASC',
  AnnotationValueIdDesc = 'ANNOTATION_VALUE_ID_DESC',
  Natural = 'NATURAL',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC'
}

/** Methods to use when ordering `StudySampleAnnotation`. */
export enum StudySampleAnnotationsOrderBy {
  AnnotationValueIdAsc = 'ANNOTATION_VALUE_ID_ASC',
  AnnotationValueIdDesc = 'ANNOTATION_VALUE_ID_DESC',
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
  studyId: Scalars['Int'];
  studySampleId: Scalars['Int'];
};

/** Represents an update to a `StudySample`. Fields that are set will be updated. */
export type StudySamplePatch = {
  h5AdObsIndex: InputMaybe<Scalars['Int']>;
  studyId: InputMaybe<Scalars['Int']>;
  studySampleId: InputMaybe<Scalars['Int']>;
};

export type StudySampleProjection = {
  __typename?: 'StudySampleProjection';
  displaySubsampling: Scalars['Boolean'];
  projection: Array<Maybe<Scalars['Float']>>;
  projectionType: ProjectionType;
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
  /** Checks for equality with the object’s `projection` field. */
  projection: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Checks for equality with the object’s `projectionType` field. */
  projectionType: InputMaybe<ProjectionType>;
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
  /** Negates the expression. */
  not: InputMaybe<StudySampleProjectionFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudySampleProjectionFilter>>;
  /** Filter by the object’s `projection` field. */
  projection: InputMaybe<FloatListFilter>;
  /** Filter by the object’s `projectionType` field. */
  projectionType: InputMaybe<ProjectionTypeFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleId` field. */
  studySampleId: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `StudySampleProjection` */
export type StudySampleProjectionInput = {
  displaySubsampling: Scalars['Boolean'];
  projection: Array<InputMaybe<Scalars['Float']>>;
  projectionType: ProjectionType;
  studyId: Scalars['Int'];
  studySampleId: Scalars['Int'];
};

export type StudySampleProjectionSubsamplingTransposed = {
  __typename?: 'StudySampleProjectionSubsamplingTransposed';
  projection: Maybe<Array<Maybe<Scalars['Float']>>>;
  projectionType: Maybe<ProjectionType>;
  /** Reads a single `Study` that is related to this `StudySampleProjectionSubsamplingTransposed`. */
  study: Maybe<Study>;
  studyId: Maybe<Scalars['Int']>;
  studySampleId: Maybe<Array<Maybe<Scalars['Int']>>>;
};

/**
 * A condition to be used against `StudySampleProjectionSubsamplingTransposed`
 * object types. All fields are tested for equality and combined with a logical ‘and.’
 */
export type StudySampleProjectionSubsamplingTransposedCondition = {
  /** Checks for equality with the object’s `projection` field. */
  projection: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Checks for equality with the object’s `projectionType` field. */
  projectionType: InputMaybe<ProjectionType>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleId` field. */
  studySampleId: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
};

/** A filter to be used against `StudySampleProjectionSubsamplingTransposed` object types. All fields are combined with a logical ‘and.’ */
export type StudySampleProjectionSubsamplingTransposedFilter = {
  /** Checks for all expressions in this list. */
  and: InputMaybe<Array<StudySampleProjectionSubsamplingTransposedFilter>>;
  /** Negates the expression. */
  not: InputMaybe<StudySampleProjectionSubsamplingTransposedFilter>;
  /** Checks for any expressions in this list. */
  or: InputMaybe<Array<StudySampleProjectionSubsamplingTransposedFilter>>;
  /** Filter by the object’s `projection` field. */
  projection: InputMaybe<FloatListFilter>;
  /** Filter by the object’s `projectionType` field. */
  projectionType: InputMaybe<ProjectionTypeFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleId` field. */
  studySampleId: InputMaybe<IntListFilter>;
};

/** Methods to use when ordering `StudySampleProjectionSubsamplingTransposed`. */
export enum StudySampleProjectionSubsamplingTransposedsOrderBy {
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

/** Methods to use when ordering `StudySampleProjection`. */
export enum StudySampleProjectionsOrderBy {
  DisplaySubsamplingAsc = 'DISPLAY_SUBSAMPLING_ASC',
  DisplaySubsamplingDesc = 'DISPLAY_SUBSAMPLING_DESC',
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

/** Methods to use when ordering `StudySample`. */
export enum StudySamplesOrderBy {
  H5AdObsIndexAsc = 'H5AD_OBS_INDEX_ASC',
  H5AdObsIndexDesc = 'H5AD_OBS_INDEX_DESC',
  Natural = 'NATURAL',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudySampleIdAsc = 'STUDY_SAMPLE_ID_ASC',
  StudySampleIdDesc = 'STUDY_SAMPLE_ID_DESC'
}

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
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Reads a single `Ontology` that is related to this `Concept`. */
  ontologyByOntid: Maybe<Ontology>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Reads a single `OmicsBase` that is related to this `OmicsProteinAntibodyTag`. */
  proteinAntibodyTag: Maybe<OmicsBase>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query: Maybe<Query>;
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
};

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

/** Methods to use when ordering `_AllUsedOntologyId`. */
export enum _AllUsedOntologyIdsOrderBy {
  Natural = 'NATURAL',
  OntologyAsc = 'ONTOLOGY_ASC',
  OntologyDesc = 'ONTOLOGY_DESC',
  OntCodeAsc = 'ONT_CODE_ASC',
  OntCodeDesc = 'ONT_CODE_DESC'
}

export type StudyInfoFragment = { __typename?: 'StudyOverview', studyId: number, studyName: string, description: string, cellCount: number, studyOntologyList: Array<{ __typename?: 'StudyOverviewOntology', ontCodes: Array<string>, labels: Array<string>, ontology: string, parentIds: Array<string> }> };

export type TreeOntologyOverviewFragment = { __typename?: 'TreeOntology', label: string, ontCode: string, ontology: string, parentOntCodePath: Array<string> };

export type DegQueryVariables = Exact<{
  studyId: Scalars['Int'];
  annotationValueId: Scalars['Int'];
}>;


export type DegQuery = { __typename?: 'Query', differentialExpressionVsList: Array<{ __typename?: 'DifferentialExpressionV', omicsId: number, studyId: number, annotationValueId: number, displayName: string, displaySymbol: string, pvalueAdj: number, log2Foldchange: number }> };

export type StudiesQueryVariables = Exact<{ [key: string]: never; }>;


export type StudiesQuery = { __typename?: 'Query', studyOverviewsList: Array<{ __typename?: 'StudyOverview', studyId: number, studyName: string, description: string, cellCount: number, studyOntologyList: Array<{ __typename?: 'StudyOverviewOntology', ontCodes: Array<string>, labels: Array<string>, ontology: string, parentIds: Array<string> }> }>, treeOntologiesList: Array<{ __typename?: 'TreeOntology', label: string, ontCode: string, ontology: string, parentOntCodePath: Array<string> }> };

export type StudyBasicsFragment = { __typename?: 'Study', studyId: number, studyName: string, studyLayersList: Array<{ __typename?: 'StudyLayer', layer: string, studyLayerId: number }>, studyOmicsTransposedList: Array<{ __typename?: 'StudyOmicsTransposed', displayName: Array<string>, displaySymbol: Array<string>, omicsId: Array<number>, omicsType: Array<OmicsType> }>, studyAnnotationGroupUisList: Array<{ __typename?: 'StudyAnnotationGroupUi', isPrimary: boolean, ordering: number, differentialExpressionCalculated: boolean, annotationGroup: { __typename?: 'AnnotationGroup', annotationGroupId: number, displayGroup: string, annotationValuesList: Array<{ __typename?: 'AnnotationValue', annotationValueId: number, displayValue: string, color: string }> } }>, studySampleAnnotationSubsamplingList: Array<{ __typename?: 'StudySampleAnnotationSubsampling', annotationValueId: number, studySampleIds: Array<number> }>, studySampleProjectionSubsamplingTransposedList: Array<{ __typename?: 'StudySampleProjectionSubsamplingTransposed', projectionType: ProjectionType, studySampleId: Array<number>, projection: Array<number> }> };

export type StudyOmicsQueryVariables = Exact<{
  studyId: Scalars['Int'];
}>;


export type StudyOmicsQuery = { __typename?: 'Query', studyOmicsList: Array<{ __typename?: 'StudyOmic', omics: { __typename?: 'OmicsBase', omicsId: number, displayName: string, displaySymbol: string } }> };

export type StudyBasicsQueryVariables = Exact<{
  studyId: Scalars['Int'];
}>;


export type StudyBasicsQuery = { __typename?: 'Query', study: { __typename?: 'Study', studyId: number, studyName: string, studyLayersList: Array<{ __typename?: 'StudyLayer', layer: string, studyLayerId: number }>, studyOmicsTransposedList: Array<{ __typename?: 'StudyOmicsTransposed', displayName: Array<string>, displaySymbol: Array<string>, omicsId: Array<number>, omicsType: Array<OmicsType> }>, studyAnnotationGroupUisList: Array<{ __typename?: 'StudyAnnotationGroupUi', isPrimary: boolean, ordering: number, differentialExpressionCalculated: boolean, annotationGroup: { __typename?: 'AnnotationGroup', annotationGroupId: number, displayGroup: string, annotationValuesList: Array<{ __typename?: 'AnnotationValue', annotationValueId: number, displayValue: string, color: string }> } }>, studySampleAnnotationSubsamplingList: Array<{ __typename?: 'StudySampleAnnotationSubsampling', annotationValueId: number, studySampleIds: Array<number> }>, studySampleProjectionSubsamplingTransposedList: Array<{ __typename?: 'StudySampleProjectionSubsamplingTransposed', projectionType: ProjectionType, studySampleId: Array<number>, projection: Array<number> }> } };

export type ExpressionByOmicsIdsQueryVariables = Exact<{
  studyLayerId: Scalars['Int'];
  omicsIds: Array<Scalars['Int']> | Scalars['Int'];
}>;


export type ExpressionByOmicsIdsQuery = { __typename?: 'Query', expressionByOmicsIdsList: Array<{ __typename?: 'ExpressionByOmic', omicsId: number, studySampleIds: Array<number>, values: Array<number> }> };

export type ExpressionViolinPlotQueryVariables = Exact<{
  studyId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  omicsId: Scalars['Int'];
  annotationGroupId: Scalars['Int'];
}>;


export type ExpressionViolinPlotQuery = { __typename?: 'Query', violinPlot: string };

export type ExpressionCorrelationTrianglePlotQueryVariables = Exact<{
  studyLayerId: Scalars['Int'];
  omicsIds: Array<Scalars['Int']> | Scalars['Int'];
}>;


export type ExpressionCorrelationTrianglePlotQuery = { __typename?: 'Query', correlationTrianglePlot: string };

export type AutocompleteQueryVariables = Exact<{
  query: Scalars['String'];
}>;


export type AutocompleteQuery = { __typename?: 'Query', autocompleteList: Array<{ __typename?: 'AutocompleteResult', isSynonymOfPreferredTerm: string, label: string, labelHighlight: string, ontCode: string, ontology: string }> };

export type OntologyOverviewFragment = { __typename?: 'Ontology', name: string, ontid: number, nodeId: string };

export type OntologiesQueryVariables = Exact<{ [key: string]: never; }>;


export type OntologiesQuery = { __typename?: 'Query', ontologiesList: Array<{ __typename?: 'Ontology', name: string, ontid: number, nodeId: string }> };

export const StudyInfoFragmentDoc = gql`
    fragment StudyInfo on StudyOverview {
  studyId
  studyName
  description
  cellCount
  studyOntologyList {
    ontCodes
    labels
    ontology
    parentIds
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
export const StudyBasicsFragmentDoc = gql`
    fragment StudyBasics on Study {
  studyId
  studyName
  studyLayersList {
    layer
    studyLayerId
  }
  studyOmicsTransposedList {
    displayName
    displaySymbol
    omicsId
    omicsType
  }
  studyAnnotationGroupUisList {
    annotationGroup {
      annotationGroupId
      displayGroup
      annotationValuesList {
        annotationValueId
        displayValue
        color
      }
    }
    isPrimary
    ordering
    differentialExpressionCalculated
  }
  studySampleAnnotationSubsamplingList {
    annotationValueId
    studySampleIds
  }
  studySampleProjectionSubsamplingTransposedList {
    projectionType
    studySampleId
    projection
  }
}
    `;
export const OntologyOverviewFragmentDoc = gql`
    fragment ontologyOverview on Ontology {
  name
  ontid
  nodeId
}
    `;
export const DegDocument = gql`
    query deg($studyId: Int!, $annotationValueId: Int!) {
  differentialExpressionVsList(
    filter: {annotationValueId: {equalTo: $annotationValueId}, studyId: {equalTo: $studyId}}
  ) {
    omicsId
    studyId
    annotationValueId
    displayName
    displaySymbol
    pvalueAdj
    log2Foldchange
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
  treeOntologiesList {
    ...TreeOntologyOverview
  }
}
    ${StudyInfoFragmentDoc}
${TreeOntologyOverviewFragmentDoc}`;

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
    ...StudyBasics
  }
}
    ${StudyBasicsFragmentDoc}`;

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
export const ExpressionByOmicsIdsDocument = gql`
    query ExpressionByOmicsIds($studyLayerId: Int!, $omicsIds: [Int!]!) {
  expressionByOmicsIdsList(
    pStudyLayerId: $studyLayerId
    pOmicsIds: $omicsIds
    pSubsamplingProjection: UMAP
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
    query ExpressionViolinPlot($studyId: Int!, $studyLayerId: Int!, $omicsId: Int!, $annotationGroupId: Int!) {
  violinPlot(
    pStudyId: $studyId
    pStudyLayerId: $studyLayerId
    pOmicsId: $omicsId
    pAnnotationGroupId: $annotationGroupId
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
export const ExpressionCorrelationTrianglePlotDocument = gql`
    query ExpressionCorrelationTrianglePlot($studyLayerId: Int!, $omicsIds: [Int!]!) {
  correlationTrianglePlot(pStudyLayerId: $studyLayerId, pOmicsIds: $omicsIds)
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
 *      studyLayerId: // value for 'studyLayerId'
 *      omicsIds: // value for 'omicsIds'
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