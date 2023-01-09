import gql from 'graphql-tag';

gql`
fragment StudyOverview on Study {
    studyId
    studyName
    cellCount
    description
    organismTaxId
    tissueNcitIds
    diseaseMeshIds
}
fragment TreeTissueOverview on TreeTissue {
    cid
    label
    ontCode
    parentCids
    parentOntCodePath
}

fragment TreeDiseaseOverview on TreeDisease {
    cid
    label
    ontCode
    parentCids
    parentOntCodePath
}

query deg($annotationValueId: Int!) {
  differentialExpressionVsList(
    filter: {annotationValueId: {equalTo: $annotationValueId}}
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
  studiesList {
    ...StudyOverview
  }
  treeTissuesList {
    ...TreeTissueOverview
  }
  treeDiseasesList {
    ...TreeDiseaseOverview
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

query autocomplete($query:String!) {
  autocompleteList(searchQuery:$query, first: 20) {
    isSynonymOfPreferredTerm
    label
    labelHighlight
    ontCode
    ontology
  }
}

`;