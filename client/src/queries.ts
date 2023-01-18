import gql from 'graphql-tag';

gql`
fragment StudyInfo on StudyOverview   {
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

fragment TreeOntologyOverview on TreeOntology {
    label
    ontCode
    ontology
    parentOntCodePath
}

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

query studies {
    studyOverviewsList {
        ...StudyInfo
    }
    treeOntologiesList {
        ...TreeOntologyOverview
    }
}

fragment AnnotationGrp on StudyAnnotationFrontendGroup {
    annotationGroupId
    isPrimary
    ordering
    displayGroup
    differentialExpressionCalculated
    annotationValuesList {
      annotationValueId
      displayValue
      color
      sampleCount
    }
} 

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
    annotationGroupsList {
        ...AnnotationGrp
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

fragment DifferentialMarker on DifferentialExpression {
    annotationValueId
    log2Foldchange
    pvalueAdj
    score
    study{
      studyName
      studyId
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
  differentialExpressionsList(filter: {omicsId: {in: $omicsIds}} orderBy: LOG2_FOLDCHANGE_DESC) {
    ...DifferentialMarker
  }
}

fragment OmicsGene on OmicsBase{
    displayName
    displaySymbol
    omicsId
    taxId
}

query allGenes {
  omicsBasesList(filter: {omicsType: {equalTo: GENE}}) {
    ...OmicsGene
    value:displaySymbol
    ontology:omicsType
  }
}

query studyOmics($studyId: Int!) {
    studyOmicsList(filter: {studyId: {equalTo: $studyId}}) {
        omics{
            omicsId
            displayName
            displaySymbol
        }
    }
}

query StudyBasics($studyId: Int!) {
  study(studyId: $studyId) {
     ...StudyBasics
  }
}


query ExpressionByOmicsIds($studyLayerId: Int!, $omicsIds: [Int!]!) {
  expressionByOmicsIdsList(pStudyLayerId:$studyLayerId, pOmicsIds:$omicsIds, pSubsamplingProjection:UMAP) {
    omicsId
    studySampleIds
    values
  }
}

query ExpressionViolinPlot($studyId: Int!, $studyLayerId: Int!, $omicsId: Int!, $annotationGroupId: Int!, $excludeAnnotationValueIds: [Int!]!) {
  violinPlot(pStudyId:$studyId, pStudyLayerId:$studyLayerId, pOmicsId: $omicsId, pAnnotationGroupId: $annotationGroupId, pExcludeAnnotationValueIds: $excludeAnnotationValueIds)
}

query ExpressionCorrelationTrianglePlot($studyId: Int!, $studyLayerId: Int!, $omicsIds: [Int!]!, $excludeAnnotationValueIds: [Int!]!) {
  correlationTrianglePlot(pStudyId:$studyId, pStudyLayerId:$studyLayerId, pOmicsIds: $omicsIds, pExcludeAnnotationValueIds: $excludeAnnotationValueIds)
}

query autocomplete($query:String!) {
  autocompleteList(searchQuery:$query, first: 20) {
    isSynonymOfPreferredTerm
    label
    labelHighlight
    ontCode
    ontology
  }
}
fragment ontologyOverview on Ontology {
    name
    ontid
    nodeId
}
query ontologies {
    ontologiesList{
        ...ontologyOverview
    }
}

query expressionByCelltype($omicsIds: [Int!]!) {
  expressionByCelltypesList(
    filter: { omicsId: { in: $omicsIds } }
    # , studyId: { in: [1, 2] }
  ) {
    celltype
    studyId
    omicsId
    q3
    exprCellsFraction
  }
}
query halfAVolcano($annotationValueId: Int!, $studyId:Int!) {
  differentialExpressionsList(
    filter: {annotationValueId: {equalTo: $annotationValueId}, studyId: {equalTo: $studyId}}
  ) {
    log2Foldchange
    pvalueAdj
  }
}


`;