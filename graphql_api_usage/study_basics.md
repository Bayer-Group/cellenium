# Study Basics

Once a study is selected, API queries are needed to populate a basic UI and keep some master data (e.g. list
of annotations, list of genes) for as long as the study is used.

Besides gene expression, we're storing measurements of protein expression or transcription factors - hence the
umbrella term "omics" instead of gene. Each study references the "omics" elements that it contains data for.
Omics elements can reference each other, e.g. a transcription factor can reference the relevant genes.

Layers could be the normalized data layer vs. an additional imputed values layer, and a study supports multiple
layers. The studyLayerId is required for various other queries to retrieve expression data.

```gql
query StudyBasics {
  study(studyId: 1) {
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
}
```
