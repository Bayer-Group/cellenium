# Expression Plots

### All expression data for a particular gene

This query is used to populate the singe gene scatter plot, where you see each cell
in e.g. a UMAP projection and the cells are colored by expression level. As there can be
"hotspots" in the projection with many overlapping cells in a small area of the plot,
density-based subsampling can be used to retrieve a fraction of the cells' expression values.

Coexpression plots could be done as well, but the subsampling parameter should be turned off
when matching multiple genes, as not to miss the sometimes few common expressed cells between two genes.

```gql
query SingleGeneScatterPlots {
  expressionByOmicsIdsList(
    pOmicsIds: [1, 2]
    pStudyLayerId: 1
    pSubsamplingProjection:UMAP
  ) {
    omicsId
    studySampleIds
    values
  }
}
```

### Expression Data For Box Plot

The study samples can be grouped by a selected sample annotation, and the boxstats (box, median, whiskers, outliers),
are calculated per gene (omics ID) and per annotation group.

As an alternative, the list of values can be retrieved.

This is just a GraphQL example, the query isn't included in the cellenium frontend.

```gql
query ExpressionBoxplotByAnnotation {
  expressionByAnnotationList(
    pStudyLayerIds: [1],
    pOmicsIds: [1, 2],
    pAnnotationGroupId: 1,
    pExcludeAnnotationValueIds: []
  ) {
    omicsId
    annotationValueId
    boxplotParams {
      n
      q1Whisker
      q1
      median
      q3
      q3Whisker
      outliers
    }
    exprSamplesFraction
  }
}

```

### Expression Data For Cross-Study Cell Type Bubble Plot

The cell type vs study bubble chart shows a gene's expression in multiple studies, grouped by the
standardized cell type annotation (using the CellO cell type prediction for consistent annotation).

This is basically the query above, with a fixed `pAnnotationGroupId` (the CellO cell type annotation group)
and multiple `pStudyLayerIds`, to cover all the studies.

### Violin plot (server side generated)

To avoid the retrieval of many data points in the web browser, plots can be done in database stored functions
and transmitted to the client in rasterized form (base64 encoded, graphQL doesn't support binary results).

```gql
query ExpressionViolinPlot($studyId: Int! $studyLayerId: Int!, $omicsId: Int!, $annotationGroupId: Int!) {
  violinPlot(pStudyId:$studyId, pStudyLayerId:$studyLayerId, pOmicsId: $omicsId, pAnnotationGroupId: $annotationGroupId)
}

{
  "studyId": 1,
  "studyLayerId": 1,
  "omicsId": 8356,
  "annotationGroupId": 1
}
```
