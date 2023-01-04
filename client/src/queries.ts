import gql from 'graphql-tag';

gql`
query studies {
  studiesList {
    attributeValueFreq
    cellCount
    clusterColorMap
    clusterHulls
    description
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
    studySampleProjectionSubsamplingTransposedList {
      projectionType
      studySampleId
      projection
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

`;