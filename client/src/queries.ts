import gql from 'graphql-tag';

gql`
  query correlatedgenes($studyId: Int!, $omicsId: Int!) {
    getCorrelatedGenesList(studyId: $studyId, omicsId: $omicsId) {
      displayName
      displaySymbol
      omicsId
      r
    }
  }

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

  fragment TreeOntologyOverview on TreeOntology {
    label
    ontCode
    ontology
    parentOntCodePath
  }

  query deg($filter: DifferentialExpressionVFilter!) {
    differentialExpressionVsList(filter: $filter) {
      omicsId
      studyId
      annotationValueId
      otherAnnotationValueId
      omicsType
      displayName
      displaySymbol
      pvalueAdj
      log2Foldchange
      linkedGenes
    }
  }

  query studies {
    studyOverviewsList {
      ...StudyInfo
    }
  }

  query studiesTreeOntologies {
    treeOntologiesList {
      ...TreeOntologyOverview
    }
  }

  query singleStudyInfo($studyId: Int!) {
    studyOverviewsList(filter: { studyId: { equalTo: $studyId } }) {
      ...StudyInfo
    }
  }

  query referenceStudiesOverview($organismTaxId: String!) {
    referenceStudyOverviewsList(condition: { organismTaxId: $organismTaxId }) {
      ...ReferenceStudyInfo
    }
  }

  query geneSpecificityStudy($studyId: Int!) {
    referenceStudyOverviewsList(filter: { studyId: { equalTo: $studyId } }) {
      ...ReferenceStudyInfo
    }
    studyAnnotationFrontendGroupsList(condition: { studyId: $studyId }) {
      annotationGroupId
      annotationValuesList {
        annotationValueId
        displayValue
        color
        sampleCount
      }
    }
  }

  fragment AnnotationGrp on StudyAnnotationFrontendGroup {
    annotationGroupId
    isPrimary
    ordering
    displayGroup
    differentialExpressionCalculated
    pairwiseDifferentialExpressionCalculated
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

  query studiesWithMarkerGenes($omicsIds: [Int!]!) {
    differentialExpressionsList(filter: { omicsId: { in: $omicsIds }, otherAnnotationValueId: { isNull: true } }, orderBy: LOG2_FOLDCHANGE_DESC) {
      ...DifferentialMarker
    }
  }

  fragment OmicsGene on OmicsBase {
    displayName
    displaySymbol
    omicsId
    taxId
  }

  query allGenes {
    omicsBasesList(filter: { omicsType: { equalTo: GENE } }) {
      ...OmicsGene
      value: displaySymbol
      ontology: omicsType
    }
  }

  query studyOmics($studyId: Int!) {
    studyOmicsList(filter: { studyId: { equalTo: $studyId } }) {
      omics {
        omicsId
        displayName
        displaySymbol
      }
    }
  }

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

  query StudyBasics($studyId: Int!) {
    study(studyId: $studyId) {
      studyId
      studyName
      cellCount
      projections
      organismTaxId
    }
    studyLayersList(condition: { studyId: $studyId }) {
      layer
      studyLayerId
    }
  }

  query StudyBasics2($studyId: Int!) {
    studyAnnotationFrontendGroupsList(condition: { studyId: $studyId }) {
      annotationGroupId
      isPrimary
      ordering
      displayGroup
      differentialExpressionCalculated
      pairwiseDifferentialExpressionCalculated
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
    studySampleProjectionSubsamplingTransposedsList(condition: { studyId: $studyId }) {
      projectionType
      studySampleId
      projection
      modality
    }
  }
  query StudyBasics3($studyId: Int!) {
    studySampleAnnotationSubsamplingsList(condition: { studyId: $studyId }) {
      annotationValueId
      studySampleIds
    }

    studyOmicsTransposedsList(condition: { studyId: $studyId }) {
      displayName
      displaySymbol
      omicsId
      omicsType
    }
  }

  query ExpressionByOmicsIds($studyLayerId: Int!, $omicsIds: [Int!]!, $subsamplingProjection: String) {
    expressionByOmicsIdsList(pStudyLayerId: $studyLayerId, pOmicsIds: $omicsIds, pSubsamplingProjection: $subsamplingProjection) {
      omicsId
      studySampleIds
      values
    }
  }

  query ExpressionViolinPlot(
    $studyId: Int!
    $studyLayerId: Int!
    $omicsId: Int!
    $annotationGroupId: Int!
    $excludeAnnotationValueIds: [Int!]!
    $annotationSecondaryGroupId: Int
  ) {
    violinPlot(
      pStudyId: $studyId
      pStudyLayerId: $studyLayerId
      pOmicsId: $omicsId
      pAnnotationGroupId: $annotationGroupId
      pExcludeAnnotationValueIds: $excludeAnnotationValueIds
      pSecondaryAnnotationGroupId: $annotationSecondaryGroupId
    )
  }

  query ExpressionGroupTable(
    $studyId: Int!
    $studyLayerId: Int!
    $omicsId: Int!
    $annotationGroupId: Int!
    $excludeAnnotationValueIds: [Int!]!
    $annotationSecondaryGroupId: Int!
  ) {
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

  query ExpressionTTest(
    $studyId: Int!
    $studyLayerId: Int!
    $omicsId: Int!
    $annotationGroupId: Int!
    $excludeAnnotationValueIds: [Int!]!
    $secondAnnotationGroupId: Int
    $sample1AnnotationValueId: Int!
    $sample1SecondAnnotationValueId: Int
    $sample2AnnotationValueId: Int!
    $sample2SecondAnnotationValueId: Int
  ) {
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

  query ExpressionCorrelationTrianglePlot($studyId: Int!, $studyLayerId: Int!, $omicsIds: [Int!]!, $excludeAnnotationValueIds: [Int!]!) {
    correlationTrianglePlot(pStudyId: $studyId, pStudyLayerId: $studyLayerId, pOmicsIds: $omicsIds, pExcludeAnnotationValueIds: $excludeAnnotationValueIds)
  }

  query autocomplete($query: String!) {
    autocompleteList(searchQuery: $query, first: 20) {
      isSynonymOfPreferredTerm
      label
      labelHighlight
      ontCode
      ontology
    }
  }

  query OmicsAutocomplete($searchQuery: String!, $omicsType: OmicsType!, $taxId: Int!) {
    omicsAutocompleteList(searchQuery: $searchQuery, omicsTypeFilter: $omicsType, taxIdFilter: $taxId) {
      displaySymbol
      labelHighlight
      omicsId
      omicsType
    }
  }

  fragment ontologyOverview on Ontology {
    name
    ontid
    nodeId
  }
  query ontologies {
    ontologiesList {
      ...ontologyOverview
    }
  }

  fragment DotPlotElement on ExpressionByAnnotation {
    studyLayerId
    omicsId
    annotationValueId
    annotationDisplayValue
    q3
    median
    exprSamplesFraction
  }

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

  query CellOAnnotationGroupId {
    annotationGroupsList(filter: { h5AdColumn: { equalTo: "CellO_celltype" } }) {
      annotationGroupId
    }
  }

  query halfAVolcano($annotationValueId: Int!, $studyId: Int!) {
    differentialExpressionsList(
      filter: { annotationValueId: { equalTo: $annotationValueId }, studyId: { equalTo: $studyId }, otherAnnotationValueId: { isNull: true } }
    ) {
      log2Foldchange
      pvalueAdj
    }
  }

  query annotationValueCoocurrence($studyId: Int!, $annotationGroupId1: Int!, $annotationGroupId2: Int!) {
    annotationValueCoocurrenceList(studyId: $studyId, annotationGroupId1: $annotationGroupId1, annotationGroupId2: $annotationGroupId2) {
      annotationValueId1
      annotationValueId2
      occurrence
    }
  }

  mutation SaveUserAnnotation($studyId: Int!, $annotationGroupName: String!, $selectedSampleIds: String!, $unexpressedSamplesOmicsIds: [Int!]) {
    userAnnotationDefine(
      input: {
        pStudyId: $studyId
        pAnnotationGroupName: $annotationGroupName
        pSelectedSampleIds: $selectedSampleIds
        pUnexpressedSamplesOmicsIds: $unexpressedSamplesOmicsIds
      }
    ) {
      clientMutationId
      integer
    }
  }

  mutation EditUserAnnotation($studyId: Int!, $annotationGroupId: Int!, $privateToUser: Boolean!) {
    userAnnotationEdit(input: { pStudyId: $studyId, pAnnotationGroupId: $annotationGroupId, pPrivateToUser: $privateToUser }) {
      boolean
    }
  }

  mutation DeleteUserAnnotation($studyId: Int!, $annotationGroupId: Int!) {
    userAnnotationDelete(input: { pStudyId: $studyId, pAnnotationGroupId: $annotationGroupId }) {
      boolean
    }
  }

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

  query studyAdminList {
    studyAdminDetailsList {
      ...StudyAdminDetails
    }
    userStudyUploadConfigured
  }

  query studyLogs($studyId: Int!) {
    studyImportLogsList(condition: { studyId: $studyId }) {
      importFile
      importLog
    }
  }

  query GeneSpecificity(
    $studyId: Int!
    $studyLayerId: Int!
    $omicsIds: [Int!]!
    $annotationGroupId: Int!
    $secondAnnotationGroupId: Int!
    $excludeAnnotationValueIds: [Int!]!
  ) {
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
    annotationGroupsList(filter: { annotationGroupId: { in: [$annotationGroupId, $secondAnnotationGroupId] } }) {
      displayGroup
      annotationGroupId
    }
  }

  mutation studyUpdate(
    $studyId: Int!
    $studyName: String!
    $description: String
    $readerPermissions: [String!]
    $adminPermissions: [String!]
    $tissueNcitIds: [String!]
    $diseaseMeshIds: [String!]
    $visible: Boolean!
    $externalWebsite: String
  ) {
    updateStudy(
      input: {
        studyId: $studyId
        patch: {
          studyName: $studyName
          description: $description
          readerPermissions: $readerPermissions
          adminPermissions: $adminPermissions
          tissueNcitIds: $tissueNcitIds
          diseaseMeshIds: $diseaseMeshIds
          visible: $visible
          externalWebsite: $externalWebsite
        }
      }
    ) {
      clientMutationId
    }
  }

  mutation studyDelete($studyId: Int!) {
    deleteAllStudyData(input: { pStudyId: $studyId }) {
      boolean
    }
  }

  mutation createStudyUpload($filename: String!) {
    createStudyUpload(input: { filename: $filename }) {
      json
    }
  }

  mutation studyDefinitionUpdate {
    studyDefinitionUpdate(input: {}) {
      boolean
    }
  }
`;
