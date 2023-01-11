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

fragment StudyBasics on Study {
    studyId
    studyName
    studyLayersList {
      layer
      studyLayerId
    }
    studyOmicsTransposedList {
      omicsId
      omicsType
      displayName
      displaySymbol
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
    studySampleAnnotationsList {
      # subsampling?
      studySampleIds
      annotationValueId
    }
    studySampleProjectionSubsamplingTransposedList {
      projectionType
      studySampleId
      projection
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
  expressionByOmicsIdsList(pStudyLayerId:$studyLayerId, pOmicsIds:$omicsIds, pSubsamplingProjection:null) {
    omicsId
    studySampleIds
    values
  }
}

query ExpressionViolinPlot($studyId: Int! $studyLayerId: Int!, $omicsId: Int!, $annotationGroupId: Int!) {
  violinPlot(pStudyId:$studyId, pStudyLayerId:$studyLayerId, pOmicsId: $omicsId, pAnnotationGroupId: $annotationGroupId)
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

`;