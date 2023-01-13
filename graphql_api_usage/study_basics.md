# Study Basics

Once a study is selected, API queries are needed to populate a basic UI and keep some master data (e.g. list
of annotations, list of genes) for as long as the study is used.

Besides gene expression, we're storing measurements of protein expression or transcription factors - hence the
umbrella term "omics" instead of gene. Each study references the "omics" elements that it contains data for.
Omics elements can reference each other, e.g. a transcription factor can reference the relevant genes.

Layers could be the normalized data layer vs. an additional imputed values layer, and a study supports multiple
layers. The studyLayerId is required for various other queries to retrieve expression data.

Sample annotations and sample UMAP projection x/y coordinates are also retrieved in the beginning. Those are
used to make scatterplots of samples colored by annotation, or colored by their expression of a requested gene.

We also include a fairly large study in the data examples, with 880k cells. The web browser is still able to
download and display so many data points, but it feels sluggish. A density based subsampling of the samples
is done in the study preparation step, according to point density in UMAP space, to bring huge sample sizes down
to 50k samples by filtering out overlapping points. Sample annotations and expression values are retrieved
in the "subsampling" GraphQL nodes below. If studies have less than 50k samples, the subsampling is just a no-op
and the sample size is not reduced.

```gql
query StudyBasics {
  study(studyId: 1) {
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
}
```
