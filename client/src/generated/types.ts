/* eslint-disable */
import { gql } from '@apollo/client';
import * as Apollo from '@apollo/client';
export type Maybe<T> = T | null;
export type InputMaybe<T> = Maybe<T>;
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

/** A filter to be used against Boolean fields. All fields are combined with a logical ‘and.’ */
export type BooleanFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom?: InputMaybe<Scalars['Boolean']>;
  /** Equal to the specified value. */
  equalTo?: InputMaybe<Scalars['Boolean']>;
  /** Greater than the specified value. */
  greaterThan?: InputMaybe<Scalars['Boolean']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo?: InputMaybe<Scalars['Boolean']>;
  /** Included in the specified list. */
  in?: InputMaybe<Array<Scalars['Boolean']>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull?: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan?: InputMaybe<Scalars['Boolean']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo?: InputMaybe<Scalars['Boolean']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom?: InputMaybe<Scalars['Boolean']>;
  /** Not equal to the specified value. */
  notEqualTo?: InputMaybe<Scalars['Boolean']>;
  /** Not included in the specified list. */
  notIn?: InputMaybe<Array<Scalars['Boolean']>>;
};

export type Concept = Node & {
  __typename?: 'Concept';
  cid: Scalars['Int'];
  /** Reads and enables pagination through a set of `ConceptHierarchy`. */
  conceptHierarchiesByCidList: Array<ConceptHierarchy>;
  /** Reads and enables pagination through a set of `ConceptHierarchy`. */
  conceptHierarchiesByParentCidList: Array<ConceptHierarchy>;
  /** Reads and enables pagination through a set of `ConceptSynonym`. */
  conceptSynonymsByCidList: Array<ConceptSynonym>;
  label?: Maybe<Scalars['String']>;
  labelTsvector?: Maybe<Scalars['String']>;
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  ontCode?: Maybe<Scalars['String']>;
  ontid?: Maybe<Scalars['Int']>;
  /** Reads a single `Ontology` that is related to this `Concept`. */
  ontologyByOntid?: Maybe<Ontology>;
  umlsDisease?: Maybe<Scalars['Boolean']>;
  umlsSemanticGroups?: Maybe<Array<Maybe<Scalars['String']>>>;
  umlsSemanticTypes?: Maybe<Array<Maybe<Scalars['String']>>>;
};


export type ConceptConceptHierarchiesByCidListArgs = {
  condition?: InputMaybe<ConceptHierarchyCondition>;
  filter?: InputMaybe<ConceptHierarchyFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
};


export type ConceptConceptHierarchiesByParentCidListArgs = {
  condition?: InputMaybe<ConceptHierarchyCondition>;
  filter?: InputMaybe<ConceptHierarchyFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
};


export type ConceptConceptSynonymsByCidListArgs = {
  condition?: InputMaybe<ConceptSynonymCondition>;
  filter?: InputMaybe<ConceptSynonymFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptSynonymsOrderBy>>;
};

/** A condition to be used against `Concept` object types. All fields are tested for equality and combined with a logical ‘and.’ */
export type ConceptCondition = {
  /** Checks for equality with the object’s `cid` field. */
  cid?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `label` field. */
  label?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `labelTsvector` field. */
  labelTsvector?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `ontCode` field. */
  ontCode?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `ontid` field. */
  ontid?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `umlsDisease` field. */
  umlsDisease?: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `umlsSemanticGroups` field. */
  umlsSemanticGroups?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `umlsSemanticTypes` field. */
  umlsSemanticTypes?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** A filter to be used against `Concept` object types. All fields are combined with a logical ‘and.’ */
export type ConceptFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<ConceptFilter>>;
  /** Filter by the object’s `cid` field. */
  cid?: InputMaybe<IntFilter>;
  /** Filter by the object’s `label` field. */
  label?: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not?: InputMaybe<ConceptFilter>;
  /** Filter by the object’s `ontCode` field. */
  ontCode?: InputMaybe<StringFilter>;
  /** Filter by the object’s `ontid` field. */
  ontid?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<ConceptFilter>>;
  /** Filter by the object’s `umlsDisease` field. */
  umlsDisease?: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `umlsSemanticGroups` field. */
  umlsSemanticGroups?: InputMaybe<StringListFilter>;
  /** Filter by the object’s `umlsSemanticTypes` field. */
  umlsSemanticTypes?: InputMaybe<StringListFilter>;
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
  cid?: Maybe<Scalars['Int']>;
  /** Reads a single `Concept` that is related to this `ConceptHierarchy`. */
  conceptByCid?: Maybe<Concept>;
  /** Reads a single `Concept` that is related to this `ConceptHierarchy`. */
  conceptByParentCid?: Maybe<Concept>;
  parentCid?: Maybe<Scalars['Int']>;
};

/**
 * A condition to be used against `ConceptHierarchy` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type ConceptHierarchyCondition = {
  /** Checks for equality with the object’s `cid` field. */
  cid?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `parentCid` field. */
  parentCid?: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `ConceptHierarchy` object types. All fields are combined with a logical ‘and.’ */
export type ConceptHierarchyFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<ConceptHierarchyFilter>>;
  /** Filter by the object’s `cid` field. */
  cid?: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not?: InputMaybe<ConceptHierarchyFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<ConceptHierarchyFilter>>;
  /** Filter by the object’s `parentCid` field. */
  parentCid?: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `ConceptHierarchy` */
export type ConceptHierarchyInput = {
  cid?: InputMaybe<Scalars['Int']>;
  parentCid?: InputMaybe<Scalars['Int']>;
};

/** An input for mutations affecting `Concept` */
export type ConceptInput = {
  cid?: InputMaybe<Scalars['Int']>;
  label?: InputMaybe<Scalars['String']>;
  labelTsvector?: InputMaybe<Scalars['String']>;
  ontCode?: InputMaybe<Scalars['String']>;
  ontid?: InputMaybe<Scalars['Int']>;
  umlsDisease?: InputMaybe<Scalars['Boolean']>;
  umlsSemanticGroups?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  umlsSemanticTypes?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** Represents an update to a `Concept`. Fields that are set will be updated. */
export type ConceptPatch = {
  cid?: InputMaybe<Scalars['Int']>;
  label?: InputMaybe<Scalars['String']>;
  labelTsvector?: InputMaybe<Scalars['String']>;
  ontCode?: InputMaybe<Scalars['String']>;
  ontid?: InputMaybe<Scalars['Int']>;
  umlsDisease?: InputMaybe<Scalars['Boolean']>;
  umlsSemanticGroups?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  umlsSemanticTypes?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

export type ConceptSynonym = {
  __typename?: 'ConceptSynonym';
  cid: Scalars['Int'];
  /** Reads a single `Concept` that is related to this `ConceptSynonym`. */
  conceptByCid?: Maybe<Concept>;
  synonym?: Maybe<Scalars['String']>;
  synonymTsvector?: Maybe<Scalars['String']>;
};

/**
 * A condition to be used against `ConceptSynonym` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type ConceptSynonymCondition = {
  /** Checks for equality with the object’s `cid` field. */
  cid?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `synonym` field. */
  synonym?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `synonymTsvector` field. */
  synonymTsvector?: InputMaybe<Scalars['String']>;
};

/** A filter to be used against `ConceptSynonym` object types. All fields are combined with a logical ‘and.’ */
export type ConceptSynonymFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<ConceptSynonymFilter>>;
  /** Filter by the object’s `cid` field. */
  cid?: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not?: InputMaybe<ConceptSynonymFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<ConceptSynonymFilter>>;
  /** Filter by the object’s `synonym` field. */
  synonym?: InputMaybe<StringFilter>;
};

/** An input for mutations affecting `ConceptSynonym` */
export type ConceptSynonymInput = {
  cid: Scalars['Int'];
  synonym?: InputMaybe<Scalars['String']>;
  synonymTsvector?: InputMaybe<Scalars['String']>;
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
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  UmlsDiseaseAsc = 'UMLS_DISEASE_ASC',
  UmlsDiseaseDesc = 'UMLS_DISEASE_DESC',
  UmlsSemanticGroupsAsc = 'UMLS_SEMANTIC_GROUPS_ASC',
  UmlsSemanticGroupsDesc = 'UMLS_SEMANTIC_GROUPS_DESC',
  UmlsSemanticTypesAsc = 'UMLS_SEMANTIC_TYPES_ASC',
  UmlsSemanticTypesDesc = 'UMLS_SEMANTIC_TYPES_DESC'
}

/** All input for the create `ConceptHierarchy` mutation. */
export type CreateConceptHierarchyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: Maybe<Scalars['String']>;
  /** Reads a single `Concept` that is related to this `ConceptHierarchy`. */
  conceptByCid?: Maybe<Concept>;
  /** Reads a single `Concept` that is related to this `ConceptHierarchy`. */
  conceptByParentCid?: Maybe<Concept>;
  /** The `ConceptHierarchy` that was created by this mutation. */
  conceptHierarchy?: Maybe<ConceptHierarchy>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
};

/** All input for the create `Concept` mutation. */
export type CreateConceptInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Concept` that was created by this mutation. */
  concept?: Maybe<Concept>;
  /** Reads a single `Ontology` that is related to this `Concept`. */
  ontologyByOntid?: Maybe<Ontology>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
};

/** All input for the create `ConceptSynonym` mutation. */
export type CreateConceptSynonymInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: Maybe<Scalars['String']>;
  /** Reads a single `Concept` that is related to this `ConceptSynonym`. */
  conceptByCid?: Maybe<Concept>;
  /** The `ConceptSynonym` that was created by this mutation. */
  conceptSynonym?: Maybe<ConceptSynonym>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
};

/** All input for the create `DifferentialExpression` mutation. */
export type CreateDifferentialExpressionInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `DifferentialExpression` to be created by this mutation. */
  differentialExpression: DifferentialExpressionInput;
};

/** The output of our create `DifferentialExpression` mutation. */
export type CreateDifferentialExpressionPayload = {
  __typename?: 'CreateDifferentialExpressionPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `DifferentialExpression` that was created by this mutation. */
  differentialExpression?: Maybe<DifferentialExpression>;
  /** Reads a single `Omic` that is related to this `DifferentialExpression`. */
  omics?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `SampleAnnotation` that is related to this `DifferentialExpression`. */
  sampleAnnotation?: Maybe<SampleAnnotation>;
  /** Reads a single `SampleAnnotationValue` that is related to this `DifferentialExpression`. */
  sampleAnnotationValue?: Maybe<SampleAnnotationValue>;
  /** Reads a single `Study` that is related to this `DifferentialExpression`. */
  study?: Maybe<Study>;
};

/** All input for the create `Expression1` mutation. */
export type CreateExpression1Input = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `Expression1` to be created by this mutation. */
  expression1: Expression1Input;
};

/** The output of our create `Expression1` mutation. */
export type CreateExpression1Payload = {
  __typename?: 'CreateExpression1Payload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Expression1` that was created by this mutation. */
  expression1?: Maybe<Expression1>;
  /** Reads a single `Omic` that is related to this `Expression1`. */
  omics?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `StudyLayer` that is related to this `Expression1`. */
  studyLayer?: Maybe<StudyLayer>;
};

/** All input for the create `Expression2` mutation. */
export type CreateExpression2Input = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `Expression2` to be created by this mutation. */
  expression2: Expression2Input;
};

/** The output of our create `Expression2` mutation. */
export type CreateExpression2Payload = {
  __typename?: 'CreateExpression2Payload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Expression2` that was created by this mutation. */
  expression2?: Maybe<Expression2>;
  /** Reads a single `Omic` that is related to this `Expression2`. */
  omics?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `StudyLayer` that is related to this `Expression2`. */
  studyLayer?: Maybe<StudyLayer>;
};

/** All input for the create `Expression3` mutation. */
export type CreateExpression3Input = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `Expression3` to be created by this mutation. */
  expression3: Expression3Input;
};

/** The output of our create `Expression3` mutation. */
export type CreateExpression3Payload = {
  __typename?: 'CreateExpression3Payload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Expression3` that was created by this mutation. */
  expression3?: Maybe<Expression3>;
  /** Reads a single `Omic` that is related to this `Expression3`. */
  omics?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `StudyLayer` that is related to this `Expression3`. */
  studyLayer?: Maybe<StudyLayer>;
};

/** All input for the create `Expression4` mutation. */
export type CreateExpression4Input = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `Expression4` to be created by this mutation. */
  expression4: Expression4Input;
};

/** The output of our create `Expression4` mutation. */
export type CreateExpression4Payload = {
  __typename?: 'CreateExpression4Payload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Expression4` that was created by this mutation. */
  expression4?: Maybe<Expression4>;
  /** Reads a single `Omic` that is related to this `Expression4`. */
  omics?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `StudyLayer` that is related to this `Expression4`. */
  studyLayer?: Maybe<StudyLayer>;
};

/** All input for the create `Expression5` mutation. */
export type CreateExpression5Input = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `Expression5` to be created by this mutation. */
  expression5: Expression5Input;
};

/** The output of our create `Expression5` mutation. */
export type CreateExpression5Payload = {
  __typename?: 'CreateExpression5Payload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Expression5` that was created by this mutation. */
  expression5?: Maybe<Expression5>;
  /** Reads a single `Omic` that is related to this `Expression5`. */
  omics?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `StudyLayer` that is related to this `Expression5`. */
  studyLayer?: Maybe<StudyLayer>;
};

/** All input for the create `Expression6` mutation. */
export type CreateExpression6Input = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `Expression6` to be created by this mutation. */
  expression6: Expression6Input;
};

/** The output of our create `Expression6` mutation. */
export type CreateExpression6Payload = {
  __typename?: 'CreateExpression6Payload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Expression6` that was created by this mutation. */
  expression6?: Maybe<Expression6>;
  /** Reads a single `Omic` that is related to this `Expression6`. */
  omics?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `StudyLayer` that is related to this `Expression6`. */
  studyLayer?: Maybe<StudyLayer>;
};

/** All input for the create `Expression7` mutation. */
export type CreateExpression7Input = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `Expression7` to be created by this mutation. */
  expression7: Expression7Input;
};

/** The output of our create `Expression7` mutation. */
export type CreateExpression7Payload = {
  __typename?: 'CreateExpression7Payload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Expression7` that was created by this mutation. */
  expression7?: Maybe<Expression7>;
  /** Reads a single `Omic` that is related to this `Expression7`. */
  omics?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `StudyLayer` that is related to this `Expression7`. */
  studyLayer?: Maybe<StudyLayer>;
};

/** All input for the create `Expression8` mutation. */
export type CreateExpression8Input = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `Expression8` to be created by this mutation. */
  expression8: Expression8Input;
};

/** The output of our create `Expression8` mutation. */
export type CreateExpression8Payload = {
  __typename?: 'CreateExpression8Payload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Expression8` that was created by this mutation. */
  expression8?: Maybe<Expression8>;
  /** Reads a single `Omic` that is related to this `Expression8`. */
  omics?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `StudyLayer` that is related to this `Expression8`. */
  studyLayer?: Maybe<StudyLayer>;
};

/** All input for the create `Expression9` mutation. */
export type CreateExpression9Input = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `Expression9` to be created by this mutation. */
  expression9: Expression9Input;
};

/** The output of our create `Expression9` mutation. */
export type CreateExpression9Payload = {
  __typename?: 'CreateExpression9Payload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Expression9` that was created by this mutation. */
  expression9?: Maybe<Expression9>;
  /** Reads a single `Omic` that is related to this `Expression9`. */
  omics?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `StudyLayer` that is related to this `Expression9`. */
  studyLayer?: Maybe<StudyLayer>;
};

/** All input for the create `Expression10` mutation. */
export type CreateExpression10Input = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `Expression10` to be created by this mutation. */
  expression10: Expression10Input;
};

/** The output of our create `Expression10` mutation. */
export type CreateExpression10Payload = {
  __typename?: 'CreateExpression10Payload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Expression10` that was created by this mutation. */
  expression10?: Maybe<Expression10>;
  /** Reads a single `Omic` that is related to this `Expression10`. */
  omics?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `StudyLayer` that is related to this `Expression10`. */
  studyLayer?: Maybe<StudyLayer>;
};

/** All input for the create `Omic` mutation. */
export type CreateOmicInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `Omic` to be created by this mutation. */
  omic: OmicInput;
};

/** The output of our create `Omic` mutation. */
export type CreateOmicPayload = {
  __typename?: 'CreateOmicPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Omic` that was created by this mutation. */
  omic?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
};

/** All input for the create `OmicsRegionGene` mutation. */
export type CreateOmicsRegionGeneInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: Maybe<Scalars['String']>;
  /** Reads a single `Omic` that is related to this `OmicsRegionGene`. */
  omics?: Maybe<Omic>;
  /** The `OmicsRegionGene` that was created by this mutation. */
  omicsRegionGene?: Maybe<OmicsRegionGene>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
};

/** All input for the create `Ontology` mutation. */
export type CreateOntologyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Ontology` that was created by this mutation. */
  ontology?: Maybe<Ontology>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
};

/** All input for the create `SampleAnnotation` mutation. */
export type CreateSampleAnnotationInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `SampleAnnotation` to be created by this mutation. */
  sampleAnnotation: SampleAnnotationInput;
};

/** The output of our create `SampleAnnotation` mutation. */
export type CreateSampleAnnotationPayload = {
  __typename?: 'CreateSampleAnnotationPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** The `SampleAnnotation` that was created by this mutation. */
  sampleAnnotation?: Maybe<SampleAnnotation>;
};

/** All input for the create `SampleAnnotationValue` mutation. */
export type CreateSampleAnnotationValueInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `SampleAnnotationValue` to be created by this mutation. */
  sampleAnnotationValue: SampleAnnotationValueInput;
};

/** The output of our create `SampleAnnotationValue` mutation. */
export type CreateSampleAnnotationValuePayload = {
  __typename?: 'CreateSampleAnnotationValuePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `SampleAnnotation` that is related to this `SampleAnnotationValue`. */
  sampleAnnotation?: Maybe<SampleAnnotation>;
  /** The `SampleAnnotationValue` that was created by this mutation. */
  sampleAnnotationValue?: Maybe<SampleAnnotationValue>;
};

/** All input for the create `Study` mutation. */
export type CreateStudyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `Study` to be created by this mutation. */
  study: StudyInput;
};

/** All input for the create `StudyLayer` mutation. */
export type CreateStudyLayerInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudyLayer`. */
  study?: Maybe<Study>;
  /** The `StudyLayer` that was created by this mutation. */
  studyLayer?: Maybe<StudyLayer>;
};

/** All input for the create `StudyOmic` mutation. */
export type CreateStudyOmicInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: Maybe<Scalars['String']>;
  /** Reads a single `Omic` that is related to this `StudyOmic`. */
  omics?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudyOmic`. */
  study?: Maybe<Study>;
  /** The `StudyOmic` that was created by this mutation. */
  studyOmic?: Maybe<StudyOmic>;
};

/** The output of our create `Study` mutation. */
export type CreateStudyPayload = {
  __typename?: 'CreateStudyPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** The `Study` that was created by this mutation. */
  study?: Maybe<Study>;
};

/** All input for the create `StudySampleAnnotationUi` mutation. */
export type CreateStudySampleAnnotationUiInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The `StudySampleAnnotationUi` to be created by this mutation. */
  studySampleAnnotationUi: StudySampleAnnotationUiInput;
};

/** The output of our create `StudySampleAnnotationUi` mutation. */
export type CreateStudySampleAnnotationUiPayload = {
  __typename?: 'CreateStudySampleAnnotationUiPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `SampleAnnotation` that is related to this `StudySampleAnnotationUi`. */
  sampleAnnotation?: Maybe<SampleAnnotation>;
  /** Reads a single `Study` that is related to this `StudySampleAnnotationUi`. */
  study?: Maybe<Study>;
  /** The `StudySampleAnnotationUi` that was created by this mutation. */
  studySampleAnnotationUi?: Maybe<StudySampleAnnotationUi>;
};

/** All input for the create `StudySample` mutation. */
export type CreateStudySampleInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudySample`. */
  study?: Maybe<Study>;
  /** The `StudySample` that was created by this mutation. */
  studySample?: Maybe<StudySample>;
};

/** A filter to be used against Datetime fields. All fields are combined with a logical ‘and.’ */
export type DatetimeFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom?: InputMaybe<Scalars['Datetime']>;
  /** Equal to the specified value. */
  equalTo?: InputMaybe<Scalars['Datetime']>;
  /** Greater than the specified value. */
  greaterThan?: InputMaybe<Scalars['Datetime']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo?: InputMaybe<Scalars['Datetime']>;
  /** Included in the specified list. */
  in?: InputMaybe<Array<Scalars['Datetime']>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull?: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan?: InputMaybe<Scalars['Datetime']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo?: InputMaybe<Scalars['Datetime']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom?: InputMaybe<Scalars['Datetime']>;
  /** Not equal to the specified value. */
  notEqualTo?: InputMaybe<Scalars['Datetime']>;
  /** Not included in the specified list. */
  notIn?: InputMaybe<Array<Scalars['Datetime']>>;
};

/** All input for the `deleteConceptByNodeId` mutation. */
export type DeleteConceptByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: InputMaybe<Scalars['String']>;
};

/** The output of our delete `Concept` mutation. */
export type DeleteConceptPayload = {
  __typename?: 'DeleteConceptPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Concept` that was deleted by this mutation. */
  concept?: Maybe<Concept>;
  deletedConceptNodeId?: Maybe<Scalars['ID']>;
  /** Reads a single `Ontology` that is related to this `Concept`. */
  ontologyByOntid?: Maybe<Ontology>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
};

/** All input for the `deleteOmicByNodeId` mutation. */
export type DeleteOmicByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `Omic` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteOmic` mutation. */
export type DeleteOmicInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  omicsId: Scalars['Int'];
};

/** The output of our delete `Omic` mutation. */
export type DeleteOmicPayload = {
  __typename?: 'DeleteOmicPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  deletedOmicNodeId?: Maybe<Scalars['ID']>;
  /** The `Omic` that was deleted by this mutation. */
  omic?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
};

/** All input for the `deleteOntologyByNodeId` mutation. */
export type DeleteOntologyByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `Ontology` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteOntology` mutation. */
export type DeleteOntologyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  ontid: Scalars['Int'];
};

/** The output of our delete `Ontology` mutation. */
export type DeleteOntologyPayload = {
  __typename?: 'DeleteOntologyPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  deletedOntologyNodeId?: Maybe<Scalars['ID']>;
  /** The `Ontology` that was deleted by this mutation. */
  ontology?: Maybe<Ontology>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
};

/** All input for the `deleteSampleAnnotationByNodeId` mutation. */
export type DeleteSampleAnnotationByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `SampleAnnotation` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteSampleAnnotation` mutation. */
export type DeleteSampleAnnotationInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  sampleAnnotationId: Scalars['Int'];
};

/** The output of our delete `SampleAnnotation` mutation. */
export type DeleteSampleAnnotationPayload = {
  __typename?: 'DeleteSampleAnnotationPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  deletedSampleAnnotationNodeId?: Maybe<Scalars['ID']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** The `SampleAnnotation` that was deleted by this mutation. */
  sampleAnnotation?: Maybe<SampleAnnotation>;
};

/** All input for the `deleteSampleAnnotationValueByNodeId` mutation. */
export type DeleteSampleAnnotationValueByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `SampleAnnotationValue` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteSampleAnnotationValue` mutation. */
export type DeleteSampleAnnotationValueInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  sampleAnnotationValueId: Scalars['Int'];
};

/** The output of our delete `SampleAnnotationValue` mutation. */
export type DeleteSampleAnnotationValuePayload = {
  __typename?: 'DeleteSampleAnnotationValuePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  deletedSampleAnnotationValueNodeId?: Maybe<Scalars['ID']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `SampleAnnotation` that is related to this `SampleAnnotationValue`. */
  sampleAnnotation?: Maybe<SampleAnnotation>;
  /** The `SampleAnnotationValue` that was deleted by this mutation. */
  sampleAnnotationValue?: Maybe<SampleAnnotationValue>;
};

/** All input for the `deleteStudyByNodeId` mutation. */
export type DeleteStudyByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `Study` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteStudy` mutation. */
export type DeleteStudyInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  studyId: Scalars['Int'];
};

/** All input for the `deleteStudyLayerByNodeId` mutation. */
export type DeleteStudyLayerByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `StudyLayer` to be deleted. */
  nodeId: Scalars['ID'];
};

/** All input for the `deleteStudyLayer` mutation. */
export type DeleteStudyLayerInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  studyLayerId: Scalars['Int'];
};

/** The output of our delete `StudyLayer` mutation. */
export type DeleteStudyLayerPayload = {
  __typename?: 'DeleteStudyLayerPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  deletedStudyLayerNodeId?: Maybe<Scalars['ID']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudyLayer`. */
  study?: Maybe<Study>;
  /** The `StudyLayer` that was deleted by this mutation. */
  studyLayer?: Maybe<StudyLayer>;
};

/** The output of our delete `Study` mutation. */
export type DeleteStudyPayload = {
  __typename?: 'DeleteStudyPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  deletedStudyNodeId?: Maybe<Scalars['ID']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** The `Study` that was deleted by this mutation. */
  study?: Maybe<Study>;
};

export type DifferentialExpression = {
  __typename?: 'DifferentialExpression';
  log2Foldchange?: Maybe<Scalars['Float']>;
  /** Reads a single `Omic` that is related to this `DifferentialExpression`. */
  omics?: Maybe<Omic>;
  omicsId?: Maybe<Scalars['Int']>;
  pvalue?: Maybe<Scalars['Float']>;
  pvalueAdj?: Maybe<Scalars['Float']>;
  /** Reads a single `SampleAnnotation` that is related to this `DifferentialExpression`. */
  sampleAnnotation?: Maybe<SampleAnnotation>;
  sampleAnnotationId: Scalars['Int'];
  /** Reads a single `SampleAnnotationValue` that is related to this `DifferentialExpression`. */
  sampleAnnotationValue?: Maybe<SampleAnnotationValue>;
  sampleAnnotationValueId: Scalars['Int'];
  score?: Maybe<Scalars['Float']>;
  /** Reads a single `Study` that is related to this `DifferentialExpression`. */
  study?: Maybe<Study>;
  studyId: Scalars['Int'];
};

/**
 * A condition to be used against `DifferentialExpression` object types. All fields
 * are tested for equality and combined with a logical ‘and.’
 */
export type DifferentialExpressionCondition = {
  /** Checks for equality with the object’s `log2Foldchange` field. */
  log2Foldchange?: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `pvalue` field. */
  pvalue?: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `pvalueAdj` field. */
  pvalueAdj?: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `sampleAnnotationId` field. */
  sampleAnnotationId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `sampleAnnotationValueId` field. */
  sampleAnnotationValueId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `score` field. */
  score?: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId?: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `DifferentialExpression` object types. All fields are combined with a logical ‘and.’ */
export type DifferentialExpressionFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<DifferentialExpressionFilter>>;
  /** Filter by the object’s `log2Foldchange` field. */
  log2Foldchange?: InputMaybe<FloatFilter>;
  /** Negates the expression. */
  not?: InputMaybe<DifferentialExpressionFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<DifferentialExpressionFilter>>;
  /** Filter by the object’s `pvalue` field. */
  pvalue?: InputMaybe<FloatFilter>;
  /** Filter by the object’s `pvalueAdj` field. */
  pvalueAdj?: InputMaybe<FloatFilter>;
  /** Filter by the object’s `sampleAnnotationId` field. */
  sampleAnnotationId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `sampleAnnotationValueId` field. */
  sampleAnnotationValueId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `score` field. */
  score?: InputMaybe<FloatFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId?: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `DifferentialExpression` */
export type DifferentialExpressionInput = {
  log2Foldchange?: InputMaybe<Scalars['Float']>;
  omicsId?: InputMaybe<Scalars['Int']>;
  pvalue?: InputMaybe<Scalars['Float']>;
  pvalueAdj?: InputMaybe<Scalars['Float']>;
  sampleAnnotationId: Scalars['Int'];
  sampleAnnotationValueId: Scalars['Int'];
  score?: InputMaybe<Scalars['Float']>;
  studyId: Scalars['Int'];
};

/** Methods to use when ordering `DifferentialExpression`. */
export enum DifferentialExpressionsOrderBy {
  Log2FoldchangeAsc = 'LOG2_FOLDCHANGE_ASC',
  Log2FoldchangeDesc = 'LOG2_FOLDCHANGE_DESC',
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  PvalueAdjAsc = 'PVALUE_ADJ_ASC',
  PvalueAdjDesc = 'PVALUE_ADJ_DESC',
  PvalueAsc = 'PVALUE_ASC',
  PvalueDesc = 'PVALUE_DESC',
  SampleAnnotationIdAsc = 'SAMPLE_ANNOTATION_ID_ASC',
  SampleAnnotationIdDesc = 'SAMPLE_ANNOTATION_ID_DESC',
  SampleAnnotationValueIdAsc = 'SAMPLE_ANNOTATION_VALUE_ID_ASC',
  SampleAnnotationValueIdDesc = 'SAMPLE_ANNOTATION_VALUE_ID_DESC',
  ScoreAsc = 'SCORE_ASC',
  ScoreDesc = 'SCORE_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

export type Expression1 = {
  __typename?: 'Expression1';
  /** Reads a single `Omic` that is related to this `Expression1`. */
  omics?: Maybe<Omic>;
  omicsId: Scalars['Int'];
  /** Reads a single `StudyLayer` that is related to this `Expression1`. */
  studyLayer?: Maybe<StudyLayer>;
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<Maybe<Scalars['Int']>>;
  values: Array<Maybe<Scalars['Float']>>;
};

/**
 * A condition to be used against `Expression1` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type Expression1Condition = {
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `values` field. */
  values?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A filter to be used against `Expression1` object types. All fields are combined with a logical ‘and.’ */
export type Expression1Filter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<Expression1Filter>>;
  /** Negates the expression. */
  not?: InputMaybe<Expression1Filter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<Expression1Filter>>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<IntListFilter>;
  /** Filter by the object’s `values` field. */
  values?: InputMaybe<FloatListFilter>;
};

/** An input for mutations affecting `Expression1` */
export type Expression1Input = {
  omicsId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<InputMaybe<Scalars['Int']>>;
  values: Array<InputMaybe<Scalars['Float']>>;
};

/** Methods to use when ordering `Expression1`. */
export enum Expression1sOrderBy {
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  StudyLayerIdAsc = 'STUDY_LAYER_ID_ASC',
  StudyLayerIdDesc = 'STUDY_LAYER_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC',
  ValuesAsc = 'VALUES_ASC',
  ValuesDesc = 'VALUES_DESC'
}

export type Expression2 = {
  __typename?: 'Expression2';
  /** Reads a single `Omic` that is related to this `Expression2`. */
  omics?: Maybe<Omic>;
  omicsId: Scalars['Int'];
  /** Reads a single `StudyLayer` that is related to this `Expression2`. */
  studyLayer?: Maybe<StudyLayer>;
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<Maybe<Scalars['Int']>>;
  values: Array<Maybe<Scalars['Float']>>;
};

/**
 * A condition to be used against `Expression2` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type Expression2Condition = {
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `values` field. */
  values?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A filter to be used against `Expression2` object types. All fields are combined with a logical ‘and.’ */
export type Expression2Filter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<Expression2Filter>>;
  /** Negates the expression. */
  not?: InputMaybe<Expression2Filter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<Expression2Filter>>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<IntListFilter>;
  /** Filter by the object’s `values` field. */
  values?: InputMaybe<FloatListFilter>;
};

/** An input for mutations affecting `Expression2` */
export type Expression2Input = {
  omicsId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<InputMaybe<Scalars['Int']>>;
  values: Array<InputMaybe<Scalars['Float']>>;
};

/** Methods to use when ordering `Expression2`. */
export enum Expression2sOrderBy {
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  StudyLayerIdAsc = 'STUDY_LAYER_ID_ASC',
  StudyLayerIdDesc = 'STUDY_LAYER_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC',
  ValuesAsc = 'VALUES_ASC',
  ValuesDesc = 'VALUES_DESC'
}

export type Expression3 = {
  __typename?: 'Expression3';
  /** Reads a single `Omic` that is related to this `Expression3`. */
  omics?: Maybe<Omic>;
  omicsId: Scalars['Int'];
  /** Reads a single `StudyLayer` that is related to this `Expression3`. */
  studyLayer?: Maybe<StudyLayer>;
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<Maybe<Scalars['Int']>>;
  values: Array<Maybe<Scalars['Float']>>;
};

/**
 * A condition to be used against `Expression3` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type Expression3Condition = {
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `values` field. */
  values?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A filter to be used against `Expression3` object types. All fields are combined with a logical ‘and.’ */
export type Expression3Filter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<Expression3Filter>>;
  /** Negates the expression. */
  not?: InputMaybe<Expression3Filter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<Expression3Filter>>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<IntListFilter>;
  /** Filter by the object’s `values` field. */
  values?: InputMaybe<FloatListFilter>;
};

/** An input for mutations affecting `Expression3` */
export type Expression3Input = {
  omicsId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<InputMaybe<Scalars['Int']>>;
  values: Array<InputMaybe<Scalars['Float']>>;
};

/** Methods to use when ordering `Expression3`. */
export enum Expression3sOrderBy {
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  StudyLayerIdAsc = 'STUDY_LAYER_ID_ASC',
  StudyLayerIdDesc = 'STUDY_LAYER_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC',
  ValuesAsc = 'VALUES_ASC',
  ValuesDesc = 'VALUES_DESC'
}

export type Expression4 = {
  __typename?: 'Expression4';
  /** Reads a single `Omic` that is related to this `Expression4`. */
  omics?: Maybe<Omic>;
  omicsId: Scalars['Int'];
  /** Reads a single `StudyLayer` that is related to this `Expression4`. */
  studyLayer?: Maybe<StudyLayer>;
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<Maybe<Scalars['Int']>>;
  values: Array<Maybe<Scalars['Float']>>;
};

/**
 * A condition to be used against `Expression4` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type Expression4Condition = {
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `values` field. */
  values?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A filter to be used against `Expression4` object types. All fields are combined with a logical ‘and.’ */
export type Expression4Filter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<Expression4Filter>>;
  /** Negates the expression. */
  not?: InputMaybe<Expression4Filter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<Expression4Filter>>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<IntListFilter>;
  /** Filter by the object’s `values` field. */
  values?: InputMaybe<FloatListFilter>;
};

/** An input for mutations affecting `Expression4` */
export type Expression4Input = {
  omicsId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<InputMaybe<Scalars['Int']>>;
  values: Array<InputMaybe<Scalars['Float']>>;
};

/** Methods to use when ordering `Expression4`. */
export enum Expression4sOrderBy {
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  StudyLayerIdAsc = 'STUDY_LAYER_ID_ASC',
  StudyLayerIdDesc = 'STUDY_LAYER_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC',
  ValuesAsc = 'VALUES_ASC',
  ValuesDesc = 'VALUES_DESC'
}

export type Expression5 = {
  __typename?: 'Expression5';
  /** Reads a single `Omic` that is related to this `Expression5`. */
  omics?: Maybe<Omic>;
  omicsId: Scalars['Int'];
  /** Reads a single `StudyLayer` that is related to this `Expression5`. */
  studyLayer?: Maybe<StudyLayer>;
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<Maybe<Scalars['Int']>>;
  values: Array<Maybe<Scalars['Float']>>;
};

/**
 * A condition to be used against `Expression5` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type Expression5Condition = {
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `values` field. */
  values?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A filter to be used against `Expression5` object types. All fields are combined with a logical ‘and.’ */
export type Expression5Filter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<Expression5Filter>>;
  /** Negates the expression. */
  not?: InputMaybe<Expression5Filter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<Expression5Filter>>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<IntListFilter>;
  /** Filter by the object’s `values` field. */
  values?: InputMaybe<FloatListFilter>;
};

/** An input for mutations affecting `Expression5` */
export type Expression5Input = {
  omicsId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<InputMaybe<Scalars['Int']>>;
  values: Array<InputMaybe<Scalars['Float']>>;
};

/** Methods to use when ordering `Expression5`. */
export enum Expression5sOrderBy {
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  StudyLayerIdAsc = 'STUDY_LAYER_ID_ASC',
  StudyLayerIdDesc = 'STUDY_LAYER_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC',
  ValuesAsc = 'VALUES_ASC',
  ValuesDesc = 'VALUES_DESC'
}

export type Expression6 = {
  __typename?: 'Expression6';
  /** Reads a single `Omic` that is related to this `Expression6`. */
  omics?: Maybe<Omic>;
  omicsId: Scalars['Int'];
  /** Reads a single `StudyLayer` that is related to this `Expression6`. */
  studyLayer?: Maybe<StudyLayer>;
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<Maybe<Scalars['Int']>>;
  values: Array<Maybe<Scalars['Float']>>;
};

/**
 * A condition to be used against `Expression6` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type Expression6Condition = {
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `values` field. */
  values?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A filter to be used against `Expression6` object types. All fields are combined with a logical ‘and.’ */
export type Expression6Filter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<Expression6Filter>>;
  /** Negates the expression. */
  not?: InputMaybe<Expression6Filter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<Expression6Filter>>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<IntListFilter>;
  /** Filter by the object’s `values` field. */
  values?: InputMaybe<FloatListFilter>;
};

/** An input for mutations affecting `Expression6` */
export type Expression6Input = {
  omicsId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<InputMaybe<Scalars['Int']>>;
  values: Array<InputMaybe<Scalars['Float']>>;
};

/** Methods to use when ordering `Expression6`. */
export enum Expression6sOrderBy {
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  StudyLayerIdAsc = 'STUDY_LAYER_ID_ASC',
  StudyLayerIdDesc = 'STUDY_LAYER_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC',
  ValuesAsc = 'VALUES_ASC',
  ValuesDesc = 'VALUES_DESC'
}

export type Expression7 = {
  __typename?: 'Expression7';
  /** Reads a single `Omic` that is related to this `Expression7`. */
  omics?: Maybe<Omic>;
  omicsId: Scalars['Int'];
  /** Reads a single `StudyLayer` that is related to this `Expression7`. */
  studyLayer?: Maybe<StudyLayer>;
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<Maybe<Scalars['Int']>>;
  values: Array<Maybe<Scalars['Float']>>;
};

/**
 * A condition to be used against `Expression7` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type Expression7Condition = {
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `values` field. */
  values?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A filter to be used against `Expression7` object types. All fields are combined with a logical ‘and.’ */
export type Expression7Filter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<Expression7Filter>>;
  /** Negates the expression. */
  not?: InputMaybe<Expression7Filter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<Expression7Filter>>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<IntListFilter>;
  /** Filter by the object’s `values` field. */
  values?: InputMaybe<FloatListFilter>;
};

/** An input for mutations affecting `Expression7` */
export type Expression7Input = {
  omicsId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<InputMaybe<Scalars['Int']>>;
  values: Array<InputMaybe<Scalars['Float']>>;
};

/** Methods to use when ordering `Expression7`. */
export enum Expression7sOrderBy {
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  StudyLayerIdAsc = 'STUDY_LAYER_ID_ASC',
  StudyLayerIdDesc = 'STUDY_LAYER_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC',
  ValuesAsc = 'VALUES_ASC',
  ValuesDesc = 'VALUES_DESC'
}

export type Expression8 = {
  __typename?: 'Expression8';
  /** Reads a single `Omic` that is related to this `Expression8`. */
  omics?: Maybe<Omic>;
  omicsId: Scalars['Int'];
  /** Reads a single `StudyLayer` that is related to this `Expression8`. */
  studyLayer?: Maybe<StudyLayer>;
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<Maybe<Scalars['Int']>>;
  values: Array<Maybe<Scalars['Float']>>;
};

/**
 * A condition to be used against `Expression8` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type Expression8Condition = {
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `values` field. */
  values?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A filter to be used against `Expression8` object types. All fields are combined with a logical ‘and.’ */
export type Expression8Filter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<Expression8Filter>>;
  /** Negates the expression. */
  not?: InputMaybe<Expression8Filter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<Expression8Filter>>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<IntListFilter>;
  /** Filter by the object’s `values` field. */
  values?: InputMaybe<FloatListFilter>;
};

/** An input for mutations affecting `Expression8` */
export type Expression8Input = {
  omicsId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<InputMaybe<Scalars['Int']>>;
  values: Array<InputMaybe<Scalars['Float']>>;
};

/** Methods to use when ordering `Expression8`. */
export enum Expression8sOrderBy {
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  StudyLayerIdAsc = 'STUDY_LAYER_ID_ASC',
  StudyLayerIdDesc = 'STUDY_LAYER_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC',
  ValuesAsc = 'VALUES_ASC',
  ValuesDesc = 'VALUES_DESC'
}

export type Expression9 = {
  __typename?: 'Expression9';
  /** Reads a single `Omic` that is related to this `Expression9`. */
  omics?: Maybe<Omic>;
  omicsId: Scalars['Int'];
  /** Reads a single `StudyLayer` that is related to this `Expression9`. */
  studyLayer?: Maybe<StudyLayer>;
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<Maybe<Scalars['Int']>>;
  values: Array<Maybe<Scalars['Float']>>;
};

/**
 * A condition to be used against `Expression9` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type Expression9Condition = {
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `values` field. */
  values?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A filter to be used against `Expression9` object types. All fields are combined with a logical ‘and.’ */
export type Expression9Filter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<Expression9Filter>>;
  /** Negates the expression. */
  not?: InputMaybe<Expression9Filter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<Expression9Filter>>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<IntListFilter>;
  /** Filter by the object’s `values` field. */
  values?: InputMaybe<FloatListFilter>;
};

/** An input for mutations affecting `Expression9` */
export type Expression9Input = {
  omicsId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<InputMaybe<Scalars['Int']>>;
  values: Array<InputMaybe<Scalars['Float']>>;
};

/** Methods to use when ordering `Expression9`. */
export enum Expression9sOrderBy {
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  StudyLayerIdAsc = 'STUDY_LAYER_ID_ASC',
  StudyLayerIdDesc = 'STUDY_LAYER_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC',
  ValuesAsc = 'VALUES_ASC',
  ValuesDesc = 'VALUES_DESC'
}

export type Expression10 = {
  __typename?: 'Expression10';
  /** Reads a single `Omic` that is related to this `Expression10`. */
  omics?: Maybe<Omic>;
  omicsId: Scalars['Int'];
  /** Reads a single `StudyLayer` that is related to this `Expression10`. */
  studyLayer?: Maybe<StudyLayer>;
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<Maybe<Scalars['Int']>>;
  values: Array<Maybe<Scalars['Float']>>;
};

/**
 * A condition to be used against `Expression10` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type Expression10Condition = {
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `values` field. */
  values?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A filter to be used against `Expression10` object types. All fields are combined with a logical ‘and.’ */
export type Expression10Filter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<Expression10Filter>>;
  /** Negates the expression. */
  not?: InputMaybe<Expression10Filter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<Expression10Filter>>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleIds` field. */
  studySampleIds?: InputMaybe<IntListFilter>;
  /** Filter by the object’s `values` field. */
  values?: InputMaybe<FloatListFilter>;
};

/** An input for mutations affecting `Expression10` */
export type Expression10Input = {
  omicsId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
  studySampleIds: Array<InputMaybe<Scalars['Int']>>;
  values: Array<InputMaybe<Scalars['Float']>>;
};

/** Methods to use when ordering `Expression10`. */
export enum Expression10sOrderBy {
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  StudyLayerIdAsc = 'STUDY_LAYER_ID_ASC',
  StudyLayerIdDesc = 'STUDY_LAYER_ID_DESC',
  StudySampleIdsAsc = 'STUDY_SAMPLE_IDS_ASC',
  StudySampleIdsDesc = 'STUDY_SAMPLE_IDS_DESC',
  ValuesAsc = 'VALUES_ASC',
  ValuesDesc = 'VALUES_DESC'
}

/** A filter to be used against Float fields. All fields are combined with a logical ‘and.’ */
export type FloatFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom?: InputMaybe<Scalars['Float']>;
  /** Equal to the specified value. */
  equalTo?: InputMaybe<Scalars['Float']>;
  /** Greater than the specified value. */
  greaterThan?: InputMaybe<Scalars['Float']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo?: InputMaybe<Scalars['Float']>;
  /** Included in the specified list. */
  in?: InputMaybe<Array<Scalars['Float']>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull?: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan?: InputMaybe<Scalars['Float']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo?: InputMaybe<Scalars['Float']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom?: InputMaybe<Scalars['Float']>;
  /** Not equal to the specified value. */
  notEqualTo?: InputMaybe<Scalars['Float']>;
  /** Not included in the specified list. */
  notIn?: InputMaybe<Array<Scalars['Float']>>;
};

/** A filter to be used against Float List fields. All fields are combined with a logical ‘and.’ */
export type FloatListFilter = {
  /** Any array item is equal to the specified value. */
  anyEqualTo?: InputMaybe<Scalars['Float']>;
  /** Any array item is greater than the specified value. */
  anyGreaterThan?: InputMaybe<Scalars['Float']>;
  /** Any array item is greater than or equal to the specified value. */
  anyGreaterThanOrEqualTo?: InputMaybe<Scalars['Float']>;
  /** Any array item is less than the specified value. */
  anyLessThan?: InputMaybe<Scalars['Float']>;
  /** Any array item is less than or equal to the specified value. */
  anyLessThanOrEqualTo?: InputMaybe<Scalars['Float']>;
  /** Any array item is not equal to the specified value. */
  anyNotEqualTo?: InputMaybe<Scalars['Float']>;
  /** Contained by the specified list of values. */
  containedBy?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Contains the specified list of values. */
  contains?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Equal to the specified value. */
  equalTo?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Greater than the specified value. */
  greaterThan?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull?: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Not equal to the specified value. */
  notEqualTo?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
  /** Overlaps the specified list of values. */
  overlaps?: InputMaybe<Array<InputMaybe<Scalars['Float']>>>;
};

/** A filter to be used against Int fields. All fields are combined with a logical ‘and.’ */
export type IntFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom?: InputMaybe<Scalars['Int']>;
  /** Equal to the specified value. */
  equalTo?: InputMaybe<Scalars['Int']>;
  /** Greater than the specified value. */
  greaterThan?: InputMaybe<Scalars['Int']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo?: InputMaybe<Scalars['Int']>;
  /** Included in the specified list. */
  in?: InputMaybe<Array<Scalars['Int']>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull?: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan?: InputMaybe<Scalars['Int']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo?: InputMaybe<Scalars['Int']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom?: InputMaybe<Scalars['Int']>;
  /** Not equal to the specified value. */
  notEqualTo?: InputMaybe<Scalars['Int']>;
  /** Not included in the specified list. */
  notIn?: InputMaybe<Array<Scalars['Int']>>;
};

/** A filter to be used against Int List fields. All fields are combined with a logical ‘and.’ */
export type IntListFilter = {
  /** Any array item is equal to the specified value. */
  anyEqualTo?: InputMaybe<Scalars['Int']>;
  /** Any array item is greater than the specified value. */
  anyGreaterThan?: InputMaybe<Scalars['Int']>;
  /** Any array item is greater than or equal to the specified value. */
  anyGreaterThanOrEqualTo?: InputMaybe<Scalars['Int']>;
  /** Any array item is less than the specified value. */
  anyLessThan?: InputMaybe<Scalars['Int']>;
  /** Any array item is less than or equal to the specified value. */
  anyLessThanOrEqualTo?: InputMaybe<Scalars['Int']>;
  /** Any array item is not equal to the specified value. */
  anyNotEqualTo?: InputMaybe<Scalars['Int']>;
  /** Contained by the specified list of values. */
  containedBy?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Contains the specified list of values. */
  contains?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Equal to the specified value. */
  equalTo?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Greater than the specified value. */
  greaterThan?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull?: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Not equal to the specified value. */
  notEqualTo?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Overlaps the specified list of values. */
  overlaps?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
};

/** A filter to be used against JSON fields. All fields are combined with a logical ‘and.’ */
export type JsonFilter = {
  /** Contained by the specified JSON. */
  containedBy?: InputMaybe<Scalars['JSON']>;
  /** Contains the specified JSON. */
  contains?: InputMaybe<Scalars['JSON']>;
  /** Contains all of the specified keys. */
  containsAllKeys?: InputMaybe<Array<Scalars['String']>>;
  /** Contains any of the specified keys. */
  containsAnyKeys?: InputMaybe<Array<Scalars['String']>>;
  /** Contains the specified key. */
  containsKey?: InputMaybe<Scalars['String']>;
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom?: InputMaybe<Scalars['JSON']>;
  /** Equal to the specified value. */
  equalTo?: InputMaybe<Scalars['JSON']>;
  /** Greater than the specified value. */
  greaterThan?: InputMaybe<Scalars['JSON']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo?: InputMaybe<Scalars['JSON']>;
  /** Included in the specified list. */
  in?: InputMaybe<Array<Scalars['JSON']>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull?: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan?: InputMaybe<Scalars['JSON']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo?: InputMaybe<Scalars['JSON']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom?: InputMaybe<Scalars['JSON']>;
  /** Not equal to the specified value. */
  notEqualTo?: InputMaybe<Scalars['JSON']>;
  /** Not included in the specified list. */
  notIn?: InputMaybe<Array<Scalars['JSON']>>;
};

/** The root mutation type which contains root level fields which mutate data. */
export type Mutation = {
  __typename?: 'Mutation';
  /** Creates a single `Concept`. */
  createConcept?: Maybe<CreateConceptPayload>;
  /** Creates a single `ConceptHierarchy`. */
  createConceptHierarchy?: Maybe<CreateConceptHierarchyPayload>;
  /** Creates a single `ConceptSynonym`. */
  createConceptSynonym?: Maybe<CreateConceptSynonymPayload>;
  /** Creates a single `DifferentialExpression`. */
  createDifferentialExpression?: Maybe<CreateDifferentialExpressionPayload>;
  /** Creates a single `Expression1`. */
  createExpression1?: Maybe<CreateExpression1Payload>;
  /** Creates a single `Expression2`. */
  createExpression2?: Maybe<CreateExpression2Payload>;
  /** Creates a single `Expression3`. */
  createExpression3?: Maybe<CreateExpression3Payload>;
  /** Creates a single `Expression4`. */
  createExpression4?: Maybe<CreateExpression4Payload>;
  /** Creates a single `Expression5`. */
  createExpression5?: Maybe<CreateExpression5Payload>;
  /** Creates a single `Expression6`. */
  createExpression6?: Maybe<CreateExpression6Payload>;
  /** Creates a single `Expression7`. */
  createExpression7?: Maybe<CreateExpression7Payload>;
  /** Creates a single `Expression8`. */
  createExpression8?: Maybe<CreateExpression8Payload>;
  /** Creates a single `Expression9`. */
  createExpression9?: Maybe<CreateExpression9Payload>;
  /** Creates a single `Expression10`. */
  createExpression10?: Maybe<CreateExpression10Payload>;
  /** Creates a single `Omic`. */
  createOmic?: Maybe<CreateOmicPayload>;
  /** Creates a single `OmicsRegionGene`. */
  createOmicsRegionGene?: Maybe<CreateOmicsRegionGenePayload>;
  /** Creates a single `Ontology`. */
  createOntology?: Maybe<CreateOntologyPayload>;
  /** Creates a single `SampleAnnotation`. */
  createSampleAnnotation?: Maybe<CreateSampleAnnotationPayload>;
  /** Creates a single `SampleAnnotationValue`. */
  createSampleAnnotationValue?: Maybe<CreateSampleAnnotationValuePayload>;
  /** Creates a single `Study`. */
  createStudy?: Maybe<CreateStudyPayload>;
  /** Creates a single `StudyLayer`. */
  createStudyLayer?: Maybe<CreateStudyLayerPayload>;
  /** Creates a single `StudyOmic`. */
  createStudyOmic?: Maybe<CreateStudyOmicPayload>;
  /** Creates a single `StudySample`. */
  createStudySample?: Maybe<CreateStudySamplePayload>;
  /** Creates a single `StudySampleAnnotationUi`. */
  createStudySampleAnnotationUi?: Maybe<CreateStudySampleAnnotationUiPayload>;
  /** Deletes a single `Concept` using a unique key. */
  deleteConcept?: Maybe<DeleteConceptPayload>;
  /** Deletes a single `Concept` using its globally unique id. */
  deleteConceptByNodeId?: Maybe<DeleteConceptPayload>;
  /** Deletes a single `Omic` using a unique key. */
  deleteOmic?: Maybe<DeleteOmicPayload>;
  /** Deletes a single `Omic` using its globally unique id. */
  deleteOmicByNodeId?: Maybe<DeleteOmicPayload>;
  /** Deletes a single `Ontology` using a unique key. */
  deleteOntology?: Maybe<DeleteOntologyPayload>;
  /** Deletes a single `Ontology` using its globally unique id. */
  deleteOntologyByNodeId?: Maybe<DeleteOntologyPayload>;
  /** Deletes a single `SampleAnnotation` using a unique key. */
  deleteSampleAnnotation?: Maybe<DeleteSampleAnnotationPayload>;
  /** Deletes a single `SampleAnnotation` using its globally unique id. */
  deleteSampleAnnotationByNodeId?: Maybe<DeleteSampleAnnotationPayload>;
  /** Deletes a single `SampleAnnotationValue` using a unique key. */
  deleteSampleAnnotationValue?: Maybe<DeleteSampleAnnotationValuePayload>;
  /** Deletes a single `SampleAnnotationValue` using its globally unique id. */
  deleteSampleAnnotationValueByNodeId?: Maybe<DeleteSampleAnnotationValuePayload>;
  /** Deletes a single `Study` using a unique key. */
  deleteStudy?: Maybe<DeleteStudyPayload>;
  /** Deletes a single `Study` using its globally unique id. */
  deleteStudyByNodeId?: Maybe<DeleteStudyPayload>;
  /** Deletes a single `StudyLayer` using a unique key. */
  deleteStudyLayer?: Maybe<DeleteStudyLayerPayload>;
  /** Deletes a single `StudyLayer` using its globally unique id. */
  deleteStudyLayerByNodeId?: Maybe<DeleteStudyLayerPayload>;
  /** Updates a single `Concept` using a unique key and a patch. */
  updateConcept?: Maybe<UpdateConceptPayload>;
  /** Updates a single `Concept` using its globally unique id and a patch. */
  updateConceptByNodeId?: Maybe<UpdateConceptPayload>;
  /** Updates a single `Omic` using a unique key and a patch. */
  updateOmic?: Maybe<UpdateOmicPayload>;
  /** Updates a single `Omic` using its globally unique id and a patch. */
  updateOmicByNodeId?: Maybe<UpdateOmicPayload>;
  /** Updates a single `Ontology` using a unique key and a patch. */
  updateOntology?: Maybe<UpdateOntologyPayload>;
  /** Updates a single `Ontology` using its globally unique id and a patch. */
  updateOntologyByNodeId?: Maybe<UpdateOntologyPayload>;
  /** Updates a single `SampleAnnotation` using a unique key and a patch. */
  updateSampleAnnotation?: Maybe<UpdateSampleAnnotationPayload>;
  /** Updates a single `SampleAnnotation` using its globally unique id and a patch. */
  updateSampleAnnotationByNodeId?: Maybe<UpdateSampleAnnotationPayload>;
  /** Updates a single `SampleAnnotationValue` using a unique key and a patch. */
  updateSampleAnnotationValue?: Maybe<UpdateSampleAnnotationValuePayload>;
  /** Updates a single `SampleAnnotationValue` using its globally unique id and a patch. */
  updateSampleAnnotationValueByNodeId?: Maybe<UpdateSampleAnnotationValuePayload>;
  /** Updates a single `Study` using a unique key and a patch. */
  updateStudy?: Maybe<UpdateStudyPayload>;
  /** Updates a single `Study` using its globally unique id and a patch. */
  updateStudyByNodeId?: Maybe<UpdateStudyPayload>;
  /** Updates a single `StudyLayer` using a unique key and a patch. */
  updateStudyLayer?: Maybe<UpdateStudyLayerPayload>;
  /** Updates a single `StudyLayer` using its globally unique id and a patch. */
  updateStudyLayerByNodeId?: Maybe<UpdateStudyLayerPayload>;
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
export type MutationCreateExpression1Args = {
  input: CreateExpression1Input;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateExpression2Args = {
  input: CreateExpression2Input;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateExpression3Args = {
  input: CreateExpression3Input;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateExpression4Args = {
  input: CreateExpression4Input;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateExpression5Args = {
  input: CreateExpression5Input;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateExpression6Args = {
  input: CreateExpression6Input;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateExpression7Args = {
  input: CreateExpression7Input;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateExpression8Args = {
  input: CreateExpression8Input;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateExpression9Args = {
  input: CreateExpression9Input;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateExpression10Args = {
  input: CreateExpression10Input;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateOmicArgs = {
  input: CreateOmicInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateOmicsRegionGeneArgs = {
  input: CreateOmicsRegionGeneInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateOntologyArgs = {
  input: CreateOntologyInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateSampleAnnotationArgs = {
  input: CreateSampleAnnotationInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateSampleAnnotationValueArgs = {
  input: CreateSampleAnnotationValueInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudyArgs = {
  input: CreateStudyInput;
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
export type MutationCreateStudySampleArgs = {
  input: CreateStudySampleInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationCreateStudySampleAnnotationUiArgs = {
  input: CreateStudySampleAnnotationUiInput;
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
export type MutationDeleteOmicArgs = {
  input: DeleteOmicInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteOmicByNodeIdArgs = {
  input: DeleteOmicByNodeIdInput;
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
export type MutationDeleteSampleAnnotationArgs = {
  input: DeleteSampleAnnotationInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteSampleAnnotationByNodeIdArgs = {
  input: DeleteSampleAnnotationByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteSampleAnnotationValueArgs = {
  input: DeleteSampleAnnotationValueInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationDeleteSampleAnnotationValueByNodeIdArgs = {
  input: DeleteSampleAnnotationValueByNodeIdInput;
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
export type MutationUpdateConceptArgs = {
  input: UpdateConceptInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateConceptByNodeIdArgs = {
  input: UpdateConceptByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOmicArgs = {
  input: UpdateOmicInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateOmicByNodeIdArgs = {
  input: UpdateOmicByNodeIdInput;
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
export type MutationUpdateSampleAnnotationArgs = {
  input: UpdateSampleAnnotationInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateSampleAnnotationByNodeIdArgs = {
  input: UpdateSampleAnnotationByNodeIdInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateSampleAnnotationValueArgs = {
  input: UpdateSampleAnnotationValueInput;
};


/** The root mutation type which contains root level fields which mutate data. */
export type MutationUpdateSampleAnnotationValueByNodeIdArgs = {
  input: UpdateSampleAnnotationValueByNodeIdInput;
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

/** An object with a globally unique `ID`. */
export type Node = {
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
};

export type Omic = Node & {
  __typename?: 'Omic';
  antibodyName?: Maybe<Scalars['String']>;
  antibodySymbol?: Maybe<Scalars['String']>;
  build?: Maybe<Scalars['String']>;
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressionsByOmicsIdList: Array<DifferentialExpression>;
  displayName?: Maybe<Scalars['String']>;
  displaySymbol: Scalars['String'];
  ensemblGeneId?: Maybe<Scalars['String']>;
  entrezGeneIds?: Maybe<Array<Maybe<Scalars['String']>>>;
  /** Reads and enables pagination through a set of `Expression1`. */
  expression1SByOmicsIdList: Array<Expression1>;
  /** Reads and enables pagination through a set of `Expression2`. */
  expression2SByOmicsIdList: Array<Expression2>;
  /** Reads and enables pagination through a set of `Expression3`. */
  expression3SByOmicsIdList: Array<Expression3>;
  /** Reads and enables pagination through a set of `Expression4`. */
  expression4SByOmicsIdList: Array<Expression4>;
  /** Reads and enables pagination through a set of `Expression5`. */
  expression5SByOmicsIdList: Array<Expression5>;
  /** Reads and enables pagination through a set of `Expression6`. */
  expression6SByOmicsIdList: Array<Expression6>;
  /** Reads and enables pagination through a set of `Expression7`. */
  expression7SByOmicsIdList: Array<Expression7>;
  /** Reads and enables pagination through a set of `Expression8`. */
  expression8SByOmicsIdList: Array<Expression8>;
  /** Reads and enables pagination through a set of `Expression9`. */
  expression9SByOmicsIdList: Array<Expression9>;
  /** Reads and enables pagination through a set of `Expression10`. */
  expression10SByOmicsIdList: Array<Expression10>;
  hgncSymbols?: Maybe<Array<Maybe<Scalars['String']>>>;
  jasparMatrixId?: Maybe<Scalars['String']>;
  linkedGenes: Array<Maybe<Scalars['Int']>>;
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  omicsId: Scalars['Int'];
  /** Reads and enables pagination through a set of `OmicsRegionGene`. */
  omicsRegionGenesByOmicsIdList: Array<OmicsRegionGene>;
  omicsType: Scalars['String'];
  regionChr?: Maybe<Scalars['String']>;
  regionEnd?: Maybe<Scalars['Int']>;
  regionStart?: Maybe<Scalars['Int']>;
  /** Reads and enables pagination through a set of `StudyOmic`. */
  studyOmicsByOmicsIdList: Array<StudyOmic>;
  taxId: Scalars['Int'];
};


export type OmicDifferentialExpressionsByOmicsIdListArgs = {
  condition?: InputMaybe<DifferentialExpressionCondition>;
  filter?: InputMaybe<DifferentialExpressionFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


export type OmicExpression1SByOmicsIdListArgs = {
  condition?: InputMaybe<Expression1Condition>;
  filter?: InputMaybe<Expression1Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression1sOrderBy>>;
};


export type OmicExpression2SByOmicsIdListArgs = {
  condition?: InputMaybe<Expression2Condition>;
  filter?: InputMaybe<Expression2Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression2sOrderBy>>;
};


export type OmicExpression3SByOmicsIdListArgs = {
  condition?: InputMaybe<Expression3Condition>;
  filter?: InputMaybe<Expression3Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression3sOrderBy>>;
};


export type OmicExpression4SByOmicsIdListArgs = {
  condition?: InputMaybe<Expression4Condition>;
  filter?: InputMaybe<Expression4Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression4sOrderBy>>;
};


export type OmicExpression5SByOmicsIdListArgs = {
  condition?: InputMaybe<Expression5Condition>;
  filter?: InputMaybe<Expression5Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression5sOrderBy>>;
};


export type OmicExpression6SByOmicsIdListArgs = {
  condition?: InputMaybe<Expression6Condition>;
  filter?: InputMaybe<Expression6Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression6sOrderBy>>;
};


export type OmicExpression7SByOmicsIdListArgs = {
  condition?: InputMaybe<Expression7Condition>;
  filter?: InputMaybe<Expression7Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression7sOrderBy>>;
};


export type OmicExpression8SByOmicsIdListArgs = {
  condition?: InputMaybe<Expression8Condition>;
  filter?: InputMaybe<Expression8Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression8sOrderBy>>;
};


export type OmicExpression9SByOmicsIdListArgs = {
  condition?: InputMaybe<Expression9Condition>;
  filter?: InputMaybe<Expression9Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression9sOrderBy>>;
};


export type OmicExpression10SByOmicsIdListArgs = {
  condition?: InputMaybe<Expression10Condition>;
  filter?: InputMaybe<Expression10Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression10sOrderBy>>;
};


export type OmicOmicsRegionGenesByOmicsIdListArgs = {
  condition?: InputMaybe<OmicsRegionGeneCondition>;
  filter?: InputMaybe<OmicsRegionGeneFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsRegionGenesOrderBy>>;
};


export type OmicStudyOmicsByOmicsIdListArgs = {
  condition?: InputMaybe<StudyOmicCondition>;
  filter?: InputMaybe<StudyOmicFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyOmicsOrderBy>>;
};

/** A condition to be used against `Omic` object types. All fields are tested for equality and combined with a logical ‘and.’ */
export type OmicCondition = {
  /** Checks for equality with the object’s `antibodyName` field. */
  antibodyName?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `antibodySymbol` field. */
  antibodySymbol?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `build` field. */
  build?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `displayName` field. */
  displayName?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `displaySymbol` field. */
  displaySymbol?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `ensemblGeneId` field. */
  ensemblGeneId?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `entrezGeneIds` field. */
  entrezGeneIds?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `hgncSymbols` field. */
  hgncSymbols?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Checks for equality with the object’s `jasparMatrixId` field. */
  jasparMatrixId?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `linkedGenes` field. */
  linkedGenes?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `omicsType` field. */
  omicsType?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `regionChr` field. */
  regionChr?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `regionEnd` field. */
  regionEnd?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `regionStart` field. */
  regionStart?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `taxId` field. */
  taxId?: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `Omic` object types. All fields are combined with a logical ‘and.’ */
export type OmicFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<OmicFilter>>;
  /** Filter by the object’s `antibodyName` field. */
  antibodyName?: InputMaybe<StringFilter>;
  /** Filter by the object’s `antibodySymbol` field. */
  antibodySymbol?: InputMaybe<StringFilter>;
  /** Filter by the object’s `build` field. */
  build?: InputMaybe<StringFilter>;
  /** Filter by the object’s `displayName` field. */
  displayName?: InputMaybe<StringFilter>;
  /** Filter by the object’s `displaySymbol` field. */
  displaySymbol?: InputMaybe<StringFilter>;
  /** Filter by the object’s `ensemblGeneId` field. */
  ensemblGeneId?: InputMaybe<StringFilter>;
  /** Filter by the object’s `entrezGeneIds` field. */
  entrezGeneIds?: InputMaybe<StringListFilter>;
  /** Filter by the object’s `hgncSymbols` field. */
  hgncSymbols?: InputMaybe<StringListFilter>;
  /** Filter by the object’s `jasparMatrixId` field. */
  jasparMatrixId?: InputMaybe<StringFilter>;
  /** Filter by the object’s `linkedGenes` field. */
  linkedGenes?: InputMaybe<IntListFilter>;
  /** Negates the expression. */
  not?: InputMaybe<OmicFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `omicsType` field. */
  omicsType?: InputMaybe<StringFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<OmicFilter>>;
  /** Filter by the object’s `regionChr` field. */
  regionChr?: InputMaybe<StringFilter>;
  /** Filter by the object’s `regionEnd` field. */
  regionEnd?: InputMaybe<IntFilter>;
  /** Filter by the object’s `regionStart` field. */
  regionStart?: InputMaybe<IntFilter>;
  /** Filter by the object’s `taxId` field. */
  taxId?: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `Omic` */
export type OmicInput = {
  antibodyName?: InputMaybe<Scalars['String']>;
  antibodySymbol?: InputMaybe<Scalars['String']>;
  build?: InputMaybe<Scalars['String']>;
  displayName?: InputMaybe<Scalars['String']>;
  displaySymbol: Scalars['String'];
  ensemblGeneId?: InputMaybe<Scalars['String']>;
  entrezGeneIds?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  hgncSymbols?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  jasparMatrixId?: InputMaybe<Scalars['String']>;
  linkedGenes?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  omicsId?: InputMaybe<Scalars['Int']>;
  omicsType: Scalars['String'];
  regionChr?: InputMaybe<Scalars['String']>;
  regionEnd?: InputMaybe<Scalars['Int']>;
  regionStart?: InputMaybe<Scalars['Int']>;
  taxId: Scalars['Int'];
};

/** Represents an update to a `Omic`. Fields that are set will be updated. */
export type OmicPatch = {
  antibodyName?: InputMaybe<Scalars['String']>;
  antibodySymbol?: InputMaybe<Scalars['String']>;
  build?: InputMaybe<Scalars['String']>;
  displayName?: InputMaybe<Scalars['String']>;
  displaySymbol?: InputMaybe<Scalars['String']>;
  ensemblGeneId?: InputMaybe<Scalars['String']>;
  entrezGeneIds?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  hgncSymbols?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  jasparMatrixId?: InputMaybe<Scalars['String']>;
  linkedGenes?: InputMaybe<Array<InputMaybe<Scalars['Int']>>>;
  omicsId?: InputMaybe<Scalars['Int']>;
  omicsType?: InputMaybe<Scalars['String']>;
  regionChr?: InputMaybe<Scalars['String']>;
  regionEnd?: InputMaybe<Scalars['Int']>;
  regionStart?: InputMaybe<Scalars['Int']>;
  taxId?: InputMaybe<Scalars['Int']>;
};

/** Methods to use when ordering `Omic`. */
export enum OmicsOrderBy {
  AntibodyNameAsc = 'ANTIBODY_NAME_ASC',
  AntibodyNameDesc = 'ANTIBODY_NAME_DESC',
  AntibodySymbolAsc = 'ANTIBODY_SYMBOL_ASC',
  AntibodySymbolDesc = 'ANTIBODY_SYMBOL_DESC',
  BuildAsc = 'BUILD_ASC',
  BuildDesc = 'BUILD_DESC',
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
  JasparMatrixIdAsc = 'JASPAR_MATRIX_ID_ASC',
  JasparMatrixIdDesc = 'JASPAR_MATRIX_ID_DESC',
  LinkedGenesAsc = 'LINKED_GENES_ASC',
  LinkedGenesDesc = 'LINKED_GENES_DESC',
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC',
  OmicsTypeAsc = 'OMICS_TYPE_ASC',
  OmicsTypeDesc = 'OMICS_TYPE_DESC',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  RegionChrAsc = 'REGION_CHR_ASC',
  RegionChrDesc = 'REGION_CHR_DESC',
  RegionEndAsc = 'REGION_END_ASC',
  RegionEndDesc = 'REGION_END_DESC',
  RegionStartAsc = 'REGION_START_ASC',
  RegionStartDesc = 'REGION_START_DESC',
  TaxIdAsc = 'TAX_ID_ASC',
  TaxIdDesc = 'TAX_ID_DESC'
}

export type OmicsRegionGene = {
  __typename?: 'OmicsRegionGene';
  ensemblGeneId?: Maybe<Scalars['String']>;
  evidence?: Maybe<Scalars['String']>;
  evidenceScore?: Maybe<Scalars['Float']>;
  evidenceSource?: Maybe<Scalars['String']>;
  gene?: Maybe<Scalars['String']>;
  /** Reads a single `Omic` that is related to this `OmicsRegionGene`. */
  omics?: Maybe<Omic>;
  omicsId: Scalars['Int'];
};

/**
 * A condition to be used against `OmicsRegionGene` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type OmicsRegionGeneCondition = {
  /** Checks for equality with the object’s `ensemblGeneId` field. */
  ensemblGeneId?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `evidence` field. */
  evidence?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `evidenceScore` field. */
  evidenceScore?: InputMaybe<Scalars['Float']>;
  /** Checks for equality with the object’s `evidenceSource` field. */
  evidenceSource?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `gene` field. */
  gene?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `OmicsRegionGene` object types. All fields are combined with a logical ‘and.’ */
export type OmicsRegionGeneFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<OmicsRegionGeneFilter>>;
  /** Filter by the object’s `ensemblGeneId` field. */
  ensemblGeneId?: InputMaybe<StringFilter>;
  /** Filter by the object’s `evidence` field. */
  evidence?: InputMaybe<StringFilter>;
  /** Filter by the object’s `evidenceScore` field. */
  evidenceScore?: InputMaybe<FloatFilter>;
  /** Filter by the object’s `evidenceSource` field. */
  evidenceSource?: InputMaybe<StringFilter>;
  /** Filter by the object’s `gene` field. */
  gene?: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not?: InputMaybe<OmicsRegionGeneFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<OmicsRegionGeneFilter>>;
};

/** An input for mutations affecting `OmicsRegionGene` */
export type OmicsRegionGeneInput = {
  ensemblGeneId?: InputMaybe<Scalars['String']>;
  evidence?: InputMaybe<Scalars['String']>;
  evidenceScore?: InputMaybe<Scalars['Float']>;
  evidenceSource?: InputMaybe<Scalars['String']>;
  gene?: InputMaybe<Scalars['String']>;
  omicsId: Scalars['Int'];
};

/** Methods to use when ordering `OmicsRegionGene`. */
export enum OmicsRegionGenesOrderBy {
  EnsemblGeneIdAsc = 'ENSEMBL_GENE_ID_ASC',
  EnsemblGeneIdDesc = 'ENSEMBL_GENE_ID_DESC',
  EvidenceAsc = 'EVIDENCE_ASC',
  EvidenceDesc = 'EVIDENCE_DESC',
  EvidenceScoreAsc = 'EVIDENCE_SCORE_ASC',
  EvidenceScoreDesc = 'EVIDENCE_SCORE_DESC',
  EvidenceSourceAsc = 'EVIDENCE_SOURCE_ASC',
  EvidenceSourceDesc = 'EVIDENCE_SOURCE_DESC',
  GeneAsc = 'GENE_ASC',
  GeneDesc = 'GENE_DESC',
  Natural = 'NATURAL',
  OmicsIdAsc = 'OMICS_ID_ASC',
  OmicsIdDesc = 'OMICS_ID_DESC'
}

/** Methods to use when ordering `Ontology`. */
export enum OntologiesOrderBy {
  NameAsc = 'NAME_ASC',
  NameDesc = 'NAME_DESC',
  Natural = 'NATURAL',
  OntidAsc = 'ONTID_ASC',
  OntidDesc = 'ONTID_DESC',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  UmlsSabAsc = 'UMLS_SAB_ASC',
  UmlsSabDesc = 'UMLS_SAB_DESC'
}

export type Ontology = Node & {
  __typename?: 'Ontology';
  /** Reads and enables pagination through a set of `Concept`. */
  conceptsByOntidList: Array<Concept>;
  name?: Maybe<Scalars['String']>;
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  ontid: Scalars['Int'];
  umlsSab?: Maybe<Scalars['String']>;
};


export type OntologyConceptsByOntidListArgs = {
  condition?: InputMaybe<ConceptCondition>;
  filter?: InputMaybe<ConceptFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptsOrderBy>>;
};

/**
 * A condition to be used against `Ontology` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type OntologyCondition = {
  /** Checks for equality with the object’s `name` field. */
  name?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `ontid` field. */
  ontid?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `umlsSab` field. */
  umlsSab?: InputMaybe<Scalars['String']>;
};

/** A filter to be used against `Ontology` object types. All fields are combined with a logical ‘and.’ */
export type OntologyFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<OntologyFilter>>;
  /** Filter by the object’s `name` field. */
  name?: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not?: InputMaybe<OntologyFilter>;
  /** Filter by the object’s `ontid` field. */
  ontid?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<OntologyFilter>>;
  /** Filter by the object’s `umlsSab` field. */
  umlsSab?: InputMaybe<StringFilter>;
};

/** An input for mutations affecting `Ontology` */
export type OntologyInput = {
  name?: InputMaybe<Scalars['String']>;
  ontid: Scalars['Int'];
  umlsSab?: InputMaybe<Scalars['String']>;
};

/** Represents an update to a `Ontology`. Fields that are set will be updated. */
export type OntologyPatch = {
  name?: InputMaybe<Scalars['String']>;
  ontid?: InputMaybe<Scalars['Int']>;
  umlsSab?: InputMaybe<Scalars['String']>;
};

/** The root query type which gives access points into the data universe. */
export type Query = Node & {
  __typename?: 'Query';
  concept?: Maybe<Concept>;
  /** Reads a single `Concept` using its globally unique `ID`. */
  conceptByNodeId?: Maybe<Concept>;
  /** Reads a set of `ConceptHierarchy`. */
  conceptHierarchiesList?: Maybe<Array<ConceptHierarchy>>;
  /** Reads a set of `ConceptSynonym`. */
  conceptSynonymsList?: Maybe<Array<ConceptSynonym>>;
  /** Reads a set of `Concept`. */
  conceptsList?: Maybe<Array<Concept>>;
  /** Reads a set of `DifferentialExpression`. */
  differentialExpressionsList?: Maybe<Array<DifferentialExpression>>;
  /** Reads a set of `Expression1`. */
  expression1sList?: Maybe<Array<Expression1>>;
  /** Reads a set of `Expression2`. */
  expression2sList?: Maybe<Array<Expression2>>;
  /** Reads a set of `Expression3`. */
  expression3sList?: Maybe<Array<Expression3>>;
  /** Reads a set of `Expression4`. */
  expression4sList?: Maybe<Array<Expression4>>;
  /** Reads a set of `Expression5`. */
  expression5sList?: Maybe<Array<Expression5>>;
  /** Reads a set of `Expression6`. */
  expression6sList?: Maybe<Array<Expression6>>;
  /** Reads a set of `Expression7`. */
  expression7sList?: Maybe<Array<Expression7>>;
  /** Reads a set of `Expression8`. */
  expression8sList?: Maybe<Array<Expression8>>;
  /** Reads a set of `Expression9`. */
  expression9sList?: Maybe<Array<Expression9>>;
  /** Reads a set of `Expression10`. */
  expression10sList?: Maybe<Array<Expression10>>;
  /** Fetches an object given its globally unique `ID`. */
  node?: Maybe<Node>;
  /** The root query type must be a `Node` to work well with Relay 1 mutations. This just resolves to `query`. */
  nodeId: Scalars['ID'];
  omic?: Maybe<Omic>;
  /** Reads a single `Omic` using its globally unique `ID`. */
  omicByNodeId?: Maybe<Omic>;
  /** Reads a set of `Omic`. */
  omicsList?: Maybe<Array<Omic>>;
  /** Reads a set of `OmicsRegionGene`. */
  omicsRegionGenesList?: Maybe<Array<OmicsRegionGene>>;
  /** Reads a set of `Ontology`. */
  ontologiesList?: Maybe<Array<Ontology>>;
  ontology?: Maybe<Ontology>;
  /** Reads a single `Ontology` using its globally unique `ID`. */
  ontologyByNodeId?: Maybe<Ontology>;
  /**
   * Exposes the root query type nested one level down. This is helpful for Relay 1
   * which can only query top level fields if they are in a particular form.
   */
  query: Query;
  sampleAnnotation?: Maybe<SampleAnnotation>;
  /** Reads a single `SampleAnnotation` using its globally unique `ID`. */
  sampleAnnotationByNodeId?: Maybe<SampleAnnotation>;
  sampleAnnotationValue?: Maybe<SampleAnnotationValue>;
  /** Reads a single `SampleAnnotationValue` using its globally unique `ID`. */
  sampleAnnotationValueByNodeId?: Maybe<SampleAnnotationValue>;
  /** Reads a set of `SampleAnnotationValue`. */
  sampleAnnotationValuesList?: Maybe<Array<SampleAnnotationValue>>;
  /** Reads a set of `SampleAnnotation`. */
  sampleAnnotationsList?: Maybe<Array<SampleAnnotation>>;
  /** Reads a set of `Study`. */
  studiesList?: Maybe<Array<Study>>;
  study?: Maybe<Study>;
  /** Reads a single `Study` using its globally unique `ID`. */
  studyByNodeId?: Maybe<Study>;
  studyLayer?: Maybe<StudyLayer>;
  /** Reads a single `StudyLayer` using its globally unique `ID`. */
  studyLayerByNodeId?: Maybe<StudyLayer>;
  /** Reads a set of `StudyLayer`. */
  studyLayersList?: Maybe<Array<StudyLayer>>;
  /** Reads a set of `StudyOmic`. */
  studyOmicsList?: Maybe<Array<StudyOmic>>;
  /** Reads a set of `StudySampleAnnotationUi`. */
  studySampleAnnotationUisList?: Maybe<Array<StudySampleAnnotationUi>>;
  /** Reads a set of `StudySample`. */
  studySamplesList?: Maybe<Array<StudySample>>;
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
export type QueryConceptHierarchiesListArgs = {
  condition?: InputMaybe<ConceptHierarchyCondition>;
  filter?: InputMaybe<ConceptHierarchyFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptHierarchiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptSynonymsListArgs = {
  condition?: InputMaybe<ConceptSynonymCondition>;
  filter?: InputMaybe<ConceptSynonymFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptSynonymsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryConceptsListArgs = {
  condition?: InputMaybe<ConceptCondition>;
  filter?: InputMaybe<ConceptFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<ConceptsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryDifferentialExpressionsListArgs = {
  condition?: InputMaybe<DifferentialExpressionCondition>;
  filter?: InputMaybe<DifferentialExpressionFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpression1sListArgs = {
  condition?: InputMaybe<Expression1Condition>;
  filter?: InputMaybe<Expression1Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression1sOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpression2sListArgs = {
  condition?: InputMaybe<Expression2Condition>;
  filter?: InputMaybe<Expression2Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression2sOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpression3sListArgs = {
  condition?: InputMaybe<Expression3Condition>;
  filter?: InputMaybe<Expression3Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression3sOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpression4sListArgs = {
  condition?: InputMaybe<Expression4Condition>;
  filter?: InputMaybe<Expression4Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression4sOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpression5sListArgs = {
  condition?: InputMaybe<Expression5Condition>;
  filter?: InputMaybe<Expression5Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression5sOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpression6sListArgs = {
  condition?: InputMaybe<Expression6Condition>;
  filter?: InputMaybe<Expression6Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression6sOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpression7sListArgs = {
  condition?: InputMaybe<Expression7Condition>;
  filter?: InputMaybe<Expression7Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression7sOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpression8sListArgs = {
  condition?: InputMaybe<Expression8Condition>;
  filter?: InputMaybe<Expression8Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression8sOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpression9sListArgs = {
  condition?: InputMaybe<Expression9Condition>;
  filter?: InputMaybe<Expression9Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression9sOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryExpression10sListArgs = {
  condition?: InputMaybe<Expression10Condition>;
  filter?: InputMaybe<Expression10Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression10sOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryNodeArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicArgs = {
  omicsId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsListArgs = {
  condition?: InputMaybe<OmicCondition>;
  filter?: InputMaybe<OmicFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOmicsRegionGenesListArgs = {
  condition?: InputMaybe<OmicsRegionGeneCondition>;
  filter?: InputMaybe<OmicsRegionGeneFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OmicsRegionGenesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryOntologiesListArgs = {
  condition?: InputMaybe<OntologyCondition>;
  filter?: InputMaybe<OntologyFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<OntologiesOrderBy>>;
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
export type QuerySampleAnnotationArgs = {
  sampleAnnotationId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QuerySampleAnnotationByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QuerySampleAnnotationValueArgs = {
  sampleAnnotationValueId: Scalars['Int'];
};


/** The root query type which gives access points into the data universe. */
export type QuerySampleAnnotationValueByNodeIdArgs = {
  nodeId: Scalars['ID'];
};


/** The root query type which gives access points into the data universe. */
export type QuerySampleAnnotationValuesListArgs = {
  condition?: InputMaybe<SampleAnnotationValueCondition>;
  filter?: InputMaybe<SampleAnnotationValueFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<SampleAnnotationValuesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QuerySampleAnnotationsListArgs = {
  condition?: InputMaybe<SampleAnnotationCondition>;
  filter?: InputMaybe<SampleAnnotationFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<SampleAnnotationsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudiesListArgs = {
  condition?: InputMaybe<StudyCondition>;
  filter?: InputMaybe<StudyFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudiesOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyArgs = {
  studyId: Scalars['Int'];
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
  condition?: InputMaybe<StudyLayerCondition>;
  filter?: InputMaybe<StudyLayerFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyLayersOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudyOmicsListArgs = {
  condition?: InputMaybe<StudyOmicCondition>;
  filter?: InputMaybe<StudyOmicFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyOmicsOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySampleAnnotationUisListArgs = {
  condition?: InputMaybe<StudySampleAnnotationUiCondition>;
  filter?: InputMaybe<StudySampleAnnotationUiFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySampleAnnotationUisOrderBy>>;
};


/** The root query type which gives access points into the data universe. */
export type QueryStudySamplesListArgs = {
  condition?: InputMaybe<StudySampleCondition>;
  filter?: InputMaybe<StudySampleFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySamplesOrderBy>>;
};

export type SampleAnnotation = Node & {
  __typename?: 'SampleAnnotation';
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressionsList: Array<DifferentialExpression>;
  display: Scalars['String'];
  h5AdColumn: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  sampleAnnotationId: Scalars['Int'];
  /** Reads and enables pagination through a set of `SampleAnnotationValue`. */
  sampleAnnotationValuesList: Array<SampleAnnotationValue>;
  /** Reads and enables pagination through a set of `StudySampleAnnotationUi`. */
  studySampleAnnotationUisList: Array<StudySampleAnnotationUi>;
};


export type SampleAnnotationDifferentialExpressionsListArgs = {
  condition?: InputMaybe<DifferentialExpressionCondition>;
  filter?: InputMaybe<DifferentialExpressionFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


export type SampleAnnotationSampleAnnotationValuesListArgs = {
  condition?: InputMaybe<SampleAnnotationValueCondition>;
  filter?: InputMaybe<SampleAnnotationValueFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<SampleAnnotationValuesOrderBy>>;
};


export type SampleAnnotationStudySampleAnnotationUisListArgs = {
  condition?: InputMaybe<StudySampleAnnotationUiCondition>;
  filter?: InputMaybe<StudySampleAnnotationUiFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySampleAnnotationUisOrderBy>>;
};

/**
 * A condition to be used against `SampleAnnotation` object types. All fields are
 * tested for equality and combined with a logical ‘and.’
 */
export type SampleAnnotationCondition = {
  /** Checks for equality with the object’s `display` field. */
  display?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `h5AdColumn` field. */
  h5AdColumn?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `sampleAnnotationId` field. */
  sampleAnnotationId?: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `SampleAnnotation` object types. All fields are combined with a logical ‘and.’ */
export type SampleAnnotationFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<SampleAnnotationFilter>>;
  /** Filter by the object’s `display` field. */
  display?: InputMaybe<StringFilter>;
  /** Filter by the object’s `h5AdColumn` field. */
  h5AdColumn?: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not?: InputMaybe<SampleAnnotationFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<SampleAnnotationFilter>>;
  /** Filter by the object’s `sampleAnnotationId` field. */
  sampleAnnotationId?: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `SampleAnnotation` */
export type SampleAnnotationInput = {
  display: Scalars['String'];
  h5AdColumn: Scalars['String'];
  sampleAnnotationId?: InputMaybe<Scalars['Int']>;
};

/** Represents an update to a `SampleAnnotation`. Fields that are set will be updated. */
export type SampleAnnotationPatch = {
  display?: InputMaybe<Scalars['String']>;
  h5AdColumn?: InputMaybe<Scalars['String']>;
  sampleAnnotationId?: InputMaybe<Scalars['Int']>;
};

export type SampleAnnotationValue = Node & {
  __typename?: 'SampleAnnotationValue';
  color?: Maybe<Scalars['String']>;
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressionsList: Array<DifferentialExpression>;
  display: Scalars['String'];
  h5AdValue: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads a single `SampleAnnotation` that is related to this `SampleAnnotationValue`. */
  sampleAnnotation?: Maybe<SampleAnnotation>;
  sampleAnnotationId: Scalars['Int'];
  sampleAnnotationValueId: Scalars['Int'];
};


export type SampleAnnotationValueDifferentialExpressionsListArgs = {
  condition?: InputMaybe<DifferentialExpressionCondition>;
  filter?: InputMaybe<DifferentialExpressionFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};

/**
 * A condition to be used against `SampleAnnotationValue` object types. All fields
 * are tested for equality and combined with a logical ‘and.’
 */
export type SampleAnnotationValueCondition = {
  /** Checks for equality with the object’s `color` field. */
  color?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `display` field. */
  display?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `h5AdValue` field. */
  h5AdValue?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `sampleAnnotationId` field. */
  sampleAnnotationId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `sampleAnnotationValueId` field. */
  sampleAnnotationValueId?: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `SampleAnnotationValue` object types. All fields are combined with a logical ‘and.’ */
export type SampleAnnotationValueFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<SampleAnnotationValueFilter>>;
  /** Filter by the object’s `color` field. */
  color?: InputMaybe<StringFilter>;
  /** Filter by the object’s `display` field. */
  display?: InputMaybe<StringFilter>;
  /** Filter by the object’s `h5AdValue` field. */
  h5AdValue?: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not?: InputMaybe<SampleAnnotationValueFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<SampleAnnotationValueFilter>>;
  /** Filter by the object’s `sampleAnnotationId` field. */
  sampleAnnotationId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `sampleAnnotationValueId` field. */
  sampleAnnotationValueId?: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `SampleAnnotationValue` */
export type SampleAnnotationValueInput = {
  color?: InputMaybe<Scalars['String']>;
  display: Scalars['String'];
  h5AdValue: Scalars['String'];
  sampleAnnotationId: Scalars['Int'];
  sampleAnnotationValueId?: InputMaybe<Scalars['Int']>;
};

/** Represents an update to a `SampleAnnotationValue`. Fields that are set will be updated. */
export type SampleAnnotationValuePatch = {
  color?: InputMaybe<Scalars['String']>;
  display?: InputMaybe<Scalars['String']>;
  h5AdValue?: InputMaybe<Scalars['String']>;
  sampleAnnotationId?: InputMaybe<Scalars['Int']>;
  sampleAnnotationValueId?: InputMaybe<Scalars['Int']>;
};

/** Methods to use when ordering `SampleAnnotationValue`. */
export enum SampleAnnotationValuesOrderBy {
  ColorAsc = 'COLOR_ASC',
  ColorDesc = 'COLOR_DESC',
  DisplayAsc = 'DISPLAY_ASC',
  DisplayDesc = 'DISPLAY_DESC',
  H5AdValueAsc = 'H5AD_VALUE_ASC',
  H5AdValueDesc = 'H5AD_VALUE_DESC',
  Natural = 'NATURAL',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  SampleAnnotationIdAsc = 'SAMPLE_ANNOTATION_ID_ASC',
  SampleAnnotationIdDesc = 'SAMPLE_ANNOTATION_ID_DESC',
  SampleAnnotationValueIdAsc = 'SAMPLE_ANNOTATION_VALUE_ID_ASC',
  SampleAnnotationValueIdDesc = 'SAMPLE_ANNOTATION_VALUE_ID_DESC'
}

/** Methods to use when ordering `SampleAnnotation`. */
export enum SampleAnnotationsOrderBy {
  DisplayAsc = 'DISPLAY_ASC',
  DisplayDesc = 'DISPLAY_DESC',
  H5AdColumnAsc = 'H5AD_COLUMN_ASC',
  H5AdColumnDesc = 'H5AD_COLUMN_DESC',
  Natural = 'NATURAL',
  PrimaryKeyAsc = 'PRIMARY_KEY_ASC',
  PrimaryKeyDesc = 'PRIMARY_KEY_DESC',
  SampleAnnotationIdAsc = 'SAMPLE_ANNOTATION_ID_ASC',
  SampleAnnotationIdDesc = 'SAMPLE_ANNOTATION_ID_DESC'
}

/** A filter to be used against String fields. All fields are combined with a logical ‘and.’ */
export type StringFilter = {
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom?: InputMaybe<Scalars['String']>;
  /** Not equal to the specified value, treating null like an ordinary value (case-insensitive). */
  distinctFromInsensitive?: InputMaybe<Scalars['String']>;
  /** Ends with the specified string (case-sensitive). */
  endsWith?: InputMaybe<Scalars['String']>;
  /** Ends with the specified string (case-insensitive). */
  endsWithInsensitive?: InputMaybe<Scalars['String']>;
  /** Equal to the specified value. */
  equalTo?: InputMaybe<Scalars['String']>;
  /** Equal to the specified value (case-insensitive). */
  equalToInsensitive?: InputMaybe<Scalars['String']>;
  /** Greater than the specified value. */
  greaterThan?: InputMaybe<Scalars['String']>;
  /** Greater than the specified value (case-insensitive). */
  greaterThanInsensitive?: InputMaybe<Scalars['String']>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo?: InputMaybe<Scalars['String']>;
  /** Greater than or equal to the specified value (case-insensitive). */
  greaterThanOrEqualToInsensitive?: InputMaybe<Scalars['String']>;
  /** Included in the specified list. */
  in?: InputMaybe<Array<Scalars['String']>>;
  /** Included in the specified list (case-insensitive). */
  inInsensitive?: InputMaybe<Array<Scalars['String']>>;
  /** Contains the specified string (case-sensitive). */
  includes?: InputMaybe<Scalars['String']>;
  /** Contains the specified string (case-insensitive). */
  includesInsensitive?: InputMaybe<Scalars['String']>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull?: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan?: InputMaybe<Scalars['String']>;
  /** Less than the specified value (case-insensitive). */
  lessThanInsensitive?: InputMaybe<Scalars['String']>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo?: InputMaybe<Scalars['String']>;
  /** Less than or equal to the specified value (case-insensitive). */
  lessThanOrEqualToInsensitive?: InputMaybe<Scalars['String']>;
  /** Matches the specified pattern (case-sensitive). An underscore (_) matches any single character; a percent sign (%) matches any sequence of zero or more characters. */
  like?: InputMaybe<Scalars['String']>;
  /** Matches the specified pattern (case-insensitive). An underscore (_) matches any single character; a percent sign (%) matches any sequence of zero or more characters. */
  likeInsensitive?: InputMaybe<Scalars['String']>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom?: InputMaybe<Scalars['String']>;
  /** Equal to the specified value, treating null like an ordinary value (case-insensitive). */
  notDistinctFromInsensitive?: InputMaybe<Scalars['String']>;
  /** Does not end with the specified string (case-sensitive). */
  notEndsWith?: InputMaybe<Scalars['String']>;
  /** Does not end with the specified string (case-insensitive). */
  notEndsWithInsensitive?: InputMaybe<Scalars['String']>;
  /** Not equal to the specified value. */
  notEqualTo?: InputMaybe<Scalars['String']>;
  /** Not equal to the specified value (case-insensitive). */
  notEqualToInsensitive?: InputMaybe<Scalars['String']>;
  /** Not included in the specified list. */
  notIn?: InputMaybe<Array<Scalars['String']>>;
  /** Not included in the specified list (case-insensitive). */
  notInInsensitive?: InputMaybe<Array<Scalars['String']>>;
  /** Does not contain the specified string (case-sensitive). */
  notIncludes?: InputMaybe<Scalars['String']>;
  /** Does not contain the specified string (case-insensitive). */
  notIncludesInsensitive?: InputMaybe<Scalars['String']>;
  /** Does not match the specified pattern (case-sensitive). An underscore (_) matches any single character; a percent sign (%) matches any sequence of zero or more characters. */
  notLike?: InputMaybe<Scalars['String']>;
  /** Does not match the specified pattern (case-insensitive). An underscore (_) matches any single character; a percent sign (%) matches any sequence of zero or more characters. */
  notLikeInsensitive?: InputMaybe<Scalars['String']>;
  /** Does not start with the specified string (case-sensitive). */
  notStartsWith?: InputMaybe<Scalars['String']>;
  /** Does not start with the specified string (case-insensitive). */
  notStartsWithInsensitive?: InputMaybe<Scalars['String']>;
  /** Starts with the specified string (case-sensitive). */
  startsWith?: InputMaybe<Scalars['String']>;
  /** Starts with the specified string (case-insensitive). */
  startsWithInsensitive?: InputMaybe<Scalars['String']>;
};

/** A filter to be used against String List fields. All fields are combined with a logical ‘and.’ */
export type StringListFilter = {
  /** Any array item is equal to the specified value. */
  anyEqualTo?: InputMaybe<Scalars['String']>;
  /** Any array item is greater than the specified value. */
  anyGreaterThan?: InputMaybe<Scalars['String']>;
  /** Any array item is greater than or equal to the specified value. */
  anyGreaterThanOrEqualTo?: InputMaybe<Scalars['String']>;
  /** Any array item is less than the specified value. */
  anyLessThan?: InputMaybe<Scalars['String']>;
  /** Any array item is less than or equal to the specified value. */
  anyLessThanOrEqualTo?: InputMaybe<Scalars['String']>;
  /** Any array item is not equal to the specified value. */
  anyNotEqualTo?: InputMaybe<Scalars['String']>;
  /** Contained by the specified list of values. */
  containedBy?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Contains the specified list of values. */
  contains?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Not equal to the specified value, treating null like an ordinary value. */
  distinctFrom?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Equal to the specified value. */
  equalTo?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Greater than the specified value. */
  greaterThan?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Greater than or equal to the specified value. */
  greaterThanOrEqualTo?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Is null (if `true` is specified) or is not null (if `false` is specified). */
  isNull?: InputMaybe<Scalars['Boolean']>;
  /** Less than the specified value. */
  lessThan?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Less than or equal to the specified value. */
  lessThanOrEqualTo?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Equal to the specified value, treating null like an ordinary value. */
  notDistinctFrom?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Not equal to the specified value. */
  notEqualTo?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
  /** Overlaps the specified list of values. */
  overlaps?: InputMaybe<Array<InputMaybe<Scalars['String']>>>;
};

/** Methods to use when ordering `Study`. */
export enum StudiesOrderBy {
  AttributeValueFreqAsc = 'ATTRIBUTE_VALUE_FREQ_ASC',
  AttributeValueFreqDesc = 'ATTRIBUTE_VALUE_FREQ_DESC',
  CellCountAsc = 'CELL_COUNT_ASC',
  CellCountDesc = 'CELL_COUNT_DESC',
  ClusterColorMapAsc = 'CLUSTER_COLOR_MAP_ASC',
  ClusterColorMapDesc = 'CLUSTER_COLOR_MAP_DESC',
  ClusterHullsAsc = 'CLUSTER_HULLS_ASC',
  ClusterHullsDesc = 'CLUSTER_HULLS_DESC',
  DescriptionAsc = 'DESCRIPTION_ASC',
  DescriptionDesc = 'DESCRIPTION_DESC',
  H5AdfileModifiedDateAsc = 'H5ADFILE_MODIFIED_DATE_ASC',
  H5AdfileModifiedDateDesc = 'H5ADFILE_MODIFIED_DATE_DESC',
  ImportStatusAsc = 'IMPORT_STATUS_ASC',
  ImportStatusDesc = 'IMPORT_STATUS_DESC',
  ImportStatusUpdatedAsc = 'IMPORT_STATUS_UPDATED_ASC',
  ImportStatusUpdatedDesc = 'IMPORT_STATUS_UPDATED_DESC',
  Natural = 'NATURAL',
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
  StudyNameDesc = 'STUDY_NAME_DESC'
}

export type Study = Node & {
  __typename?: 'Study';
  attributeValueFreq?: Maybe<Scalars['JSON']>;
  cellCount?: Maybe<Scalars['Int']>;
  clusterColorMap?: Maybe<Scalars['JSON']>;
  clusterHulls?: Maybe<Scalars['JSON']>;
  description?: Maybe<Scalars['String']>;
  /** Reads and enables pagination through a set of `DifferentialExpression`. */
  differentialExpressionsList: Array<DifferentialExpression>;
  h5AdfileModifiedDate?: Maybe<Scalars['Datetime']>;
  importStatus?: Maybe<Scalars['String']>;
  importStatusUpdated?: Maybe<Scalars['Datetime']>;
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  plotCoords?: Maybe<Scalars['JSON']>;
  projectionCellCoords?: Maybe<Scalars['JSON']>;
  projectionCellIndices?: Maybe<Scalars['JSON']>;
  studyId: Scalars['Int'];
  /** Reads and enables pagination through a set of `StudyLayer`. */
  studyLayersList: Array<StudyLayer>;
  studyName: Scalars['String'];
  /** Reads and enables pagination through a set of `StudyOmic`. */
  studyOmicsList: Array<StudyOmic>;
  /** Reads and enables pagination through a set of `StudySampleAnnotationUi`. */
  studySampleAnnotationUisList: Array<StudySampleAnnotationUi>;
  /** Reads and enables pagination through a set of `StudySample`. */
  studySamplesList: Array<StudySample>;
};


export type StudyDifferentialExpressionsListArgs = {
  condition?: InputMaybe<DifferentialExpressionCondition>;
  filter?: InputMaybe<DifferentialExpressionFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<DifferentialExpressionsOrderBy>>;
};


export type StudyStudyLayersListArgs = {
  condition?: InputMaybe<StudyLayerCondition>;
  filter?: InputMaybe<StudyLayerFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyLayersOrderBy>>;
};


export type StudyStudyOmicsListArgs = {
  condition?: InputMaybe<StudyOmicCondition>;
  filter?: InputMaybe<StudyOmicFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudyOmicsOrderBy>>;
};


export type StudyStudySampleAnnotationUisListArgs = {
  condition?: InputMaybe<StudySampleAnnotationUiCondition>;
  filter?: InputMaybe<StudySampleAnnotationUiFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySampleAnnotationUisOrderBy>>;
};


export type StudyStudySamplesListArgs = {
  condition?: InputMaybe<StudySampleCondition>;
  filter?: InputMaybe<StudySampleFilter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<StudySamplesOrderBy>>;
};

/** A condition to be used against `Study` object types. All fields are tested for equality and combined with a logical ‘and.’ */
export type StudyCondition = {
  /** Checks for equality with the object’s `attributeValueFreq` field. */
  attributeValueFreq?: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `cellCount` field. */
  cellCount?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `clusterColorMap` field. */
  clusterColorMap?: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `clusterHulls` field. */
  clusterHulls?: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `description` field. */
  description?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `h5AdfileModifiedDate` field. */
  h5AdfileModifiedDate?: InputMaybe<Scalars['Datetime']>;
  /** Checks for equality with the object’s `importStatus` field. */
  importStatus?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `importStatusUpdated` field. */
  importStatusUpdated?: InputMaybe<Scalars['Datetime']>;
  /** Checks for equality with the object’s `plotCoords` field. */
  plotCoords?: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `projectionCellCoords` field. */
  projectionCellCoords?: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `projectionCellIndices` field. */
  projectionCellIndices?: InputMaybe<Scalars['JSON']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyName` field. */
  studyName?: InputMaybe<Scalars['String']>;
};

/** A filter to be used against `Study` object types. All fields are combined with a logical ‘and.’ */
export type StudyFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<StudyFilter>>;
  /** Filter by the object’s `attributeValueFreq` field. */
  attributeValueFreq?: InputMaybe<JsonFilter>;
  /** Filter by the object’s `cellCount` field. */
  cellCount?: InputMaybe<IntFilter>;
  /** Filter by the object’s `clusterColorMap` field. */
  clusterColorMap?: InputMaybe<JsonFilter>;
  /** Filter by the object’s `clusterHulls` field. */
  clusterHulls?: InputMaybe<JsonFilter>;
  /** Filter by the object’s `description` field. */
  description?: InputMaybe<StringFilter>;
  /** Filter by the object’s `h5AdfileModifiedDate` field. */
  h5AdfileModifiedDate?: InputMaybe<DatetimeFilter>;
  /** Filter by the object’s `importStatus` field. */
  importStatus?: InputMaybe<StringFilter>;
  /** Filter by the object’s `importStatusUpdated` field. */
  importStatusUpdated?: InputMaybe<DatetimeFilter>;
  /** Negates the expression. */
  not?: InputMaybe<StudyFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<StudyFilter>>;
  /** Filter by the object’s `plotCoords` field. */
  plotCoords?: InputMaybe<JsonFilter>;
  /** Filter by the object’s `projectionCellCoords` field. */
  projectionCellCoords?: InputMaybe<JsonFilter>;
  /** Filter by the object’s `projectionCellIndices` field. */
  projectionCellIndices?: InputMaybe<JsonFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyName` field. */
  studyName?: InputMaybe<StringFilter>;
};

/** An input for mutations affecting `Study` */
export type StudyInput = {
  attributeValueFreq?: InputMaybe<Scalars['JSON']>;
  cellCount?: InputMaybe<Scalars['Int']>;
  clusterColorMap?: InputMaybe<Scalars['JSON']>;
  clusterHulls?: InputMaybe<Scalars['JSON']>;
  description?: InputMaybe<Scalars['String']>;
  h5AdfileModifiedDate?: InputMaybe<Scalars['Datetime']>;
  importStatus?: InputMaybe<Scalars['String']>;
  importStatusUpdated?: InputMaybe<Scalars['Datetime']>;
  plotCoords?: InputMaybe<Scalars['JSON']>;
  projectionCellCoords?: InputMaybe<Scalars['JSON']>;
  projectionCellIndices?: InputMaybe<Scalars['JSON']>;
  studyId?: InputMaybe<Scalars['Int']>;
  studyName: Scalars['String'];
};

export type StudyLayer = Node & {
  __typename?: 'StudyLayer';
  /** Reads and enables pagination through a set of `Expression1`. */
  expression1sList: Array<Expression1>;
  /** Reads and enables pagination through a set of `Expression2`. */
  expression2sList: Array<Expression2>;
  /** Reads and enables pagination through a set of `Expression3`. */
  expression3sList: Array<Expression3>;
  /** Reads and enables pagination through a set of `Expression4`. */
  expression4sList: Array<Expression4>;
  /** Reads and enables pagination through a set of `Expression5`. */
  expression5sList: Array<Expression5>;
  /** Reads and enables pagination through a set of `Expression6`. */
  expression6sList: Array<Expression6>;
  /** Reads and enables pagination through a set of `Expression7`. */
  expression7sList: Array<Expression7>;
  /** Reads and enables pagination through a set of `Expression8`. */
  expression8sList: Array<Expression8>;
  /** Reads and enables pagination through a set of `Expression9`. */
  expression9sList: Array<Expression9>;
  /** Reads and enables pagination through a set of `Expression10`. */
  expression10sList: Array<Expression10>;
  layer: Scalars['String'];
  /** A globally unique identifier. Can be used in various places throughout the system to identify this single value. */
  nodeId: Scalars['ID'];
  /** Reads a single `Study` that is related to this `StudyLayer`. */
  study?: Maybe<Study>;
  studyId: Scalars['Int'];
  studyLayerId: Scalars['Int'];
};


export type StudyLayerExpression1sListArgs = {
  condition?: InputMaybe<Expression1Condition>;
  filter?: InputMaybe<Expression1Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression1sOrderBy>>;
};


export type StudyLayerExpression2sListArgs = {
  condition?: InputMaybe<Expression2Condition>;
  filter?: InputMaybe<Expression2Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression2sOrderBy>>;
};


export type StudyLayerExpression3sListArgs = {
  condition?: InputMaybe<Expression3Condition>;
  filter?: InputMaybe<Expression3Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression3sOrderBy>>;
};


export type StudyLayerExpression4sListArgs = {
  condition?: InputMaybe<Expression4Condition>;
  filter?: InputMaybe<Expression4Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression4sOrderBy>>;
};


export type StudyLayerExpression5sListArgs = {
  condition?: InputMaybe<Expression5Condition>;
  filter?: InputMaybe<Expression5Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression5sOrderBy>>;
};


export type StudyLayerExpression6sListArgs = {
  condition?: InputMaybe<Expression6Condition>;
  filter?: InputMaybe<Expression6Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression6sOrderBy>>;
};


export type StudyLayerExpression7sListArgs = {
  condition?: InputMaybe<Expression7Condition>;
  filter?: InputMaybe<Expression7Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression7sOrderBy>>;
};


export type StudyLayerExpression8sListArgs = {
  condition?: InputMaybe<Expression8Condition>;
  filter?: InputMaybe<Expression8Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression8sOrderBy>>;
};


export type StudyLayerExpression9sListArgs = {
  condition?: InputMaybe<Expression9Condition>;
  filter?: InputMaybe<Expression9Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression9sOrderBy>>;
};


export type StudyLayerExpression10sListArgs = {
  condition?: InputMaybe<Expression10Condition>;
  filter?: InputMaybe<Expression10Filter>;
  first?: InputMaybe<Scalars['Int']>;
  offset?: InputMaybe<Scalars['Int']>;
  orderBy?: InputMaybe<Array<Expression10sOrderBy>>;
};

/**
 * A condition to be used against `StudyLayer` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type StudyLayerCondition = {
  /** Checks for equality with the object’s `layer` field. */
  layer?: InputMaybe<Scalars['String']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudyLayer` object types. All fields are combined with a logical ‘and.’ */
export type StudyLayerFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<StudyLayerFilter>>;
  /** Filter by the object’s `layer` field. */
  layer?: InputMaybe<StringFilter>;
  /** Negates the expression. */
  not?: InputMaybe<StudyLayerFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<StudyLayerFilter>>;
  /** Filter by the object’s `studyId` field. */
  studyId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyLayerId` field. */
  studyLayerId?: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `StudyLayer` */
export type StudyLayerInput = {
  layer: Scalars['String'];
  studyId: Scalars['Int'];
  studyLayerId?: InputMaybe<Scalars['Int']>;
};

/** Represents an update to a `StudyLayer`. Fields that are set will be updated. */
export type StudyLayerPatch = {
  layer?: InputMaybe<Scalars['String']>;
  studyId?: InputMaybe<Scalars['Int']>;
  studyLayerId?: InputMaybe<Scalars['Int']>;
};

/** Methods to use when ordering `StudyLayer`. */
export enum StudyLayersOrderBy {
  LayerAsc = 'LAYER_ASC',
  LayerDesc = 'LAYER_DESC',
  Natural = 'NATURAL',
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
  /** Reads a single `Omic` that is related to this `StudyOmic`. */
  omics?: Maybe<Omic>;
  omicsId: Scalars['Int'];
  regionEnd?: Maybe<Scalars['Int']>;
  regionStart?: Maybe<Scalars['Int']>;
  /** Reads a single `Study` that is related to this `StudyOmic`. */
  study?: Maybe<Study>;
  studyId: Scalars['Int'];
};

/**
 * A condition to be used against `StudyOmic` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type StudyOmicCondition = {
  /** Checks for equality with the object’s `h5AdVarIndex` field. */
  h5AdVarIndex?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `omicsId` field. */
  omicsId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `regionEnd` field. */
  regionEnd?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `regionStart` field. */
  regionStart?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId?: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudyOmic` object types. All fields are combined with a logical ‘and.’ */
export type StudyOmicFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<StudyOmicFilter>>;
  /** Filter by the object’s `h5AdVarIndex` field. */
  h5AdVarIndex?: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not?: InputMaybe<StudyOmicFilter>;
  /** Filter by the object’s `omicsId` field. */
  omicsId?: InputMaybe<IntFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<StudyOmicFilter>>;
  /** Filter by the object’s `regionEnd` field. */
  regionEnd?: InputMaybe<IntFilter>;
  /** Filter by the object’s `regionStart` field. */
  regionStart?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId?: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `StudyOmic` */
export type StudyOmicInput = {
  h5AdVarIndex: Scalars['Int'];
  omicsId: Scalars['Int'];
  regionEnd?: InputMaybe<Scalars['Int']>;
  regionStart?: InputMaybe<Scalars['Int']>;
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

/** Represents an update to a `Study`. Fields that are set will be updated. */
export type StudyPatch = {
  attributeValueFreq?: InputMaybe<Scalars['JSON']>;
  cellCount?: InputMaybe<Scalars['Int']>;
  clusterColorMap?: InputMaybe<Scalars['JSON']>;
  clusterHulls?: InputMaybe<Scalars['JSON']>;
  description?: InputMaybe<Scalars['String']>;
  h5AdfileModifiedDate?: InputMaybe<Scalars['Datetime']>;
  importStatus?: InputMaybe<Scalars['String']>;
  importStatusUpdated?: InputMaybe<Scalars['Datetime']>;
  plotCoords?: InputMaybe<Scalars['JSON']>;
  projectionCellCoords?: InputMaybe<Scalars['JSON']>;
  projectionCellIndices?: InputMaybe<Scalars['JSON']>;
  studyId?: InputMaybe<Scalars['Int']>;
  studyName?: InputMaybe<Scalars['String']>;
};

export type StudySample = {
  __typename?: 'StudySample';
  displaySubsampling: Scalars['Boolean'];
  h5AdObsIndex: Scalars['Int'];
  /** Reads a single `Study` that is related to this `StudySample`. */
  study?: Maybe<Study>;
  studyId: Scalars['Int'];
  studySampleId: Scalars['Int'];
};

export type StudySampleAnnotationUi = {
  __typename?: 'StudySampleAnnotationUi';
  differentialExpressionCalculated: Scalars['Boolean'];
  isPrimary: Scalars['Boolean'];
  ordering: Scalars['Int'];
  /** Reads a single `SampleAnnotation` that is related to this `StudySampleAnnotationUi`. */
  sampleAnnotation?: Maybe<SampleAnnotation>;
  sampleAnnotationId: Scalars['Int'];
  /** Reads a single `Study` that is related to this `StudySampleAnnotationUi`. */
  study?: Maybe<Study>;
  studyId: Scalars['Int'];
};

/**
 * A condition to be used against `StudySampleAnnotationUi` object types. All
 * fields are tested for equality and combined with a logical ‘and.’
 */
export type StudySampleAnnotationUiCondition = {
  /** Checks for equality with the object’s `differentialExpressionCalculated` field. */
  differentialExpressionCalculated?: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `isPrimary` field. */
  isPrimary?: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `ordering` field. */
  ordering?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `sampleAnnotationId` field. */
  sampleAnnotationId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId?: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudySampleAnnotationUi` object types. All fields are combined with a logical ‘and.’ */
export type StudySampleAnnotationUiFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<StudySampleAnnotationUiFilter>>;
  /** Filter by the object’s `differentialExpressionCalculated` field. */
  differentialExpressionCalculated?: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `isPrimary` field. */
  isPrimary?: InputMaybe<BooleanFilter>;
  /** Negates the expression. */
  not?: InputMaybe<StudySampleAnnotationUiFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<StudySampleAnnotationUiFilter>>;
  /** Filter by the object’s `ordering` field. */
  ordering?: InputMaybe<IntFilter>;
  /** Filter by the object’s `sampleAnnotationId` field. */
  sampleAnnotationId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studyId` field. */
  studyId?: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `StudySampleAnnotationUi` */
export type StudySampleAnnotationUiInput = {
  differentialExpressionCalculated: Scalars['Boolean'];
  isPrimary: Scalars['Boolean'];
  ordering: Scalars['Int'];
  sampleAnnotationId: Scalars['Int'];
  studyId: Scalars['Int'];
};

/** Methods to use when ordering `StudySampleAnnotationUi`. */
export enum StudySampleAnnotationUisOrderBy {
  DifferentialExpressionCalculatedAsc = 'DIFFERENTIAL_EXPRESSION_CALCULATED_ASC',
  DifferentialExpressionCalculatedDesc = 'DIFFERENTIAL_EXPRESSION_CALCULATED_DESC',
  IsPrimaryAsc = 'IS_PRIMARY_ASC',
  IsPrimaryDesc = 'IS_PRIMARY_DESC',
  Natural = 'NATURAL',
  OrderingAsc = 'ORDERING_ASC',
  OrderingDesc = 'ORDERING_DESC',
  SampleAnnotationIdAsc = 'SAMPLE_ANNOTATION_ID_ASC',
  SampleAnnotationIdDesc = 'SAMPLE_ANNOTATION_ID_DESC',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC'
}

/**
 * A condition to be used against `StudySample` object types. All fields are tested
 * for equality and combined with a logical ‘and.’
 */
export type StudySampleCondition = {
  /** Checks for equality with the object’s `displaySubsampling` field. */
  displaySubsampling?: InputMaybe<Scalars['Boolean']>;
  /** Checks for equality with the object’s `h5AdObsIndex` field. */
  h5AdObsIndex?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studyId` field. */
  studyId?: InputMaybe<Scalars['Int']>;
  /** Checks for equality with the object’s `studySampleId` field. */
  studySampleId?: InputMaybe<Scalars['Int']>;
};

/** A filter to be used against `StudySample` object types. All fields are combined with a logical ‘and.’ */
export type StudySampleFilter = {
  /** Checks for all expressions in this list. */
  and?: InputMaybe<Array<StudySampleFilter>>;
  /** Filter by the object’s `displaySubsampling` field. */
  displaySubsampling?: InputMaybe<BooleanFilter>;
  /** Filter by the object’s `h5AdObsIndex` field. */
  h5AdObsIndex?: InputMaybe<IntFilter>;
  /** Negates the expression. */
  not?: InputMaybe<StudySampleFilter>;
  /** Checks for any expressions in this list. */
  or?: InputMaybe<Array<StudySampleFilter>>;
  /** Filter by the object’s `studyId` field. */
  studyId?: InputMaybe<IntFilter>;
  /** Filter by the object’s `studySampleId` field. */
  studySampleId?: InputMaybe<IntFilter>;
};

/** An input for mutations affecting `StudySample` */
export type StudySampleInput = {
  displaySubsampling: Scalars['Boolean'];
  h5AdObsIndex: Scalars['Int'];
  studyId: Scalars['Int'];
  studySampleId: Scalars['Int'];
};

/** Methods to use when ordering `StudySample`. */
export enum StudySamplesOrderBy {
  DisplaySubsamplingAsc = 'DISPLAY_SUBSAMPLING_ASC',
  DisplaySubsamplingDesc = 'DISPLAY_SUBSAMPLING_DESC',
  H5AdObsIndexAsc = 'H5AD_OBS_INDEX_ASC',
  H5AdObsIndexDesc = 'H5AD_OBS_INDEX_DESC',
  Natural = 'NATURAL',
  StudyIdAsc = 'STUDY_ID_ASC',
  StudyIdDesc = 'STUDY_ID_DESC',
  StudySampleIdAsc = 'STUDY_SAMPLE_ID_ASC',
  StudySampleIdDesc = 'STUDY_SAMPLE_ID_DESC'
}

/** All input for the `updateConceptByNodeId` mutation. */
export type UpdateConceptByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Concept` that was updated by this mutation. */
  concept?: Maybe<Concept>;
  /** Reads a single `Ontology` that is related to this `Concept`. */
  ontologyByOntid?: Maybe<Ontology>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
};

/** All input for the `updateOmicByNodeId` mutation. */
export type UpdateOmicByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `Omic` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `Omic` being updated. */
  patch: OmicPatch;
};

/** All input for the `updateOmic` mutation. */
export type UpdateOmicInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  omicsId: Scalars['Int'];
  /** An object where the defined keys will be set on the `Omic` being updated. */
  patch: OmicPatch;
};

/** The output of our update `Omic` mutation. */
export type UpdateOmicPayload = {
  __typename?: 'UpdateOmicPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Omic` that was updated by this mutation. */
  omic?: Maybe<Omic>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
};

/** All input for the `updateOntologyByNodeId` mutation. */
export type UpdateOntologyByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: Maybe<Scalars['String']>;
  /** The `Ontology` that was updated by this mutation. */
  ontology?: Maybe<Ontology>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
};

/** All input for the `updateSampleAnnotationByNodeId` mutation. */
export type UpdateSampleAnnotationByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `SampleAnnotation` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `SampleAnnotation` being updated. */
  patch: SampleAnnotationPatch;
};

/** All input for the `updateSampleAnnotation` mutation. */
export type UpdateSampleAnnotationInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** An object where the defined keys will be set on the `SampleAnnotation` being updated. */
  patch: SampleAnnotationPatch;
  sampleAnnotationId: Scalars['Int'];
};

/** The output of our update `SampleAnnotation` mutation. */
export type UpdateSampleAnnotationPayload = {
  __typename?: 'UpdateSampleAnnotationPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** The `SampleAnnotation` that was updated by this mutation. */
  sampleAnnotation?: Maybe<SampleAnnotation>;
};

/** All input for the `updateSampleAnnotationValueByNodeId` mutation. */
export type UpdateSampleAnnotationValueByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** The globally unique `ID` which will identify a single `SampleAnnotationValue` to be updated. */
  nodeId: Scalars['ID'];
  /** An object where the defined keys will be set on the `SampleAnnotationValue` being updated. */
  patch: SampleAnnotationValuePatch;
};

/** All input for the `updateSampleAnnotationValue` mutation. */
export type UpdateSampleAnnotationValueInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
  /** An object where the defined keys will be set on the `SampleAnnotationValue` being updated. */
  patch: SampleAnnotationValuePatch;
  sampleAnnotationValueId: Scalars['Int'];
};

/** The output of our update `SampleAnnotationValue` mutation. */
export type UpdateSampleAnnotationValuePayload = {
  __typename?: 'UpdateSampleAnnotationValuePayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `SampleAnnotation` that is related to this `SampleAnnotationValue`. */
  sampleAnnotation?: Maybe<SampleAnnotation>;
  /** The `SampleAnnotationValue` that was updated by this mutation. */
  sampleAnnotationValue?: Maybe<SampleAnnotationValue>;
};

/** All input for the `updateStudyByNodeId` mutation. */
export type UpdateStudyByNodeIdInput = {
  /**
   * An arbitrary string value with no semantic meaning. Will be included in the
   * payload verbatim. May be used to track mutations by the client.
   */
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: InputMaybe<Scalars['String']>;
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
  clientMutationId?: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** Reads a single `Study` that is related to this `StudyLayer`. */
  study?: Maybe<Study>;
  /** The `StudyLayer` that was updated by this mutation. */
  studyLayer?: Maybe<StudyLayer>;
};

/** The output of our update `Study` mutation. */
export type UpdateStudyPayload = {
  __typename?: 'UpdateStudyPayload';
  /**
   * The exact same `clientMutationId` that was provided in the mutation input,
   * unchanged and unused. May be used by a client to track mutations.
   */
  clientMutationId?: Maybe<Scalars['String']>;
  /** Our root query field type. Allows us to run any query from our mutation payload. */
  query?: Maybe<Query>;
  /** The `Study` that was updated by this mutation. */
  study?: Maybe<Study>;
};

export type StudiesQueryVariables = Exact<{ [key: string]: never; }>;


export type StudiesQuery = { __typename?: 'Query', studiesList?: Array<{ __typename?: 'Study', attributeValueFreq?: any | null, cellCount?: number | null, clusterColorMap?: any | null, clusterHulls?: any | null, description?: string | null }> | null };


export const StudiesDocument = gql`
    query studies {
  studiesList {
    attributeValueFreq
    cellCount
    clusterColorMap
    clusterHulls
    description
  }
}
    `;

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