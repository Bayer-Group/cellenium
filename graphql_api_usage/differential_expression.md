# Differentially Expressioned Genes

Based on a particular cell annotation (e.g. a subset of cells has the xyz cell type "annotation value"
in the "cell type" annotation group), differentially expressed genes for cells with the particular
annotation vs. all the other cells can be retrieved.

The differentially expressed genes are pre-calculated for all groups during study import.

An annotation value is specific to an annotation group, so the group ID isn't relevant in the query.

```gql
query DifferentialExpression {
  differentialExpressionsList(
    filter: {studyId: {equalTo: 1}, annotationValueId: {equalTo: 8}}
  ) {
    omicsId
    pvalue
    pvalueAdj
    score
    log2Foldchange
  }
}
```
