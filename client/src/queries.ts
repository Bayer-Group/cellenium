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
    studyOmicsList {
      omics {
        omicsId
        displaySymbol
        displayName
      }
    }
    studySampleAnnotationUisList {
      annotation {
        annotationValuesList {
          annotationValueId
          displayValue
          color
        }
        annotationId
        displayGroup
      }
      isPrimary
      ordering
      differentialExpressionCalculated
    }
}

query StudyBasics($studyId: Int!) {
  study(studyId: $studyId) {
     ...StudyBasics
  }
}

`;