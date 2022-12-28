# API for Study Overview

### List of Studies, Ontology Tree contents

In studyOverviewsList, the ontology Ids and ParentIds can be used to filter the list of studies for
selected ontology terms (a selected parent term must match studies which have a more specialized term).

The treeTissuesList etc. contain all information to construct an ontology tree. The API filters the
ontology for actual study annotations and their parent concepts, so that the tree contains concepts
which yield search results. As the ontologies are not trees but networks (i.e. a node can have
multiple parents), the API performs an additional step to force the nodes into a tree structure
(settling on parents with most relevance only). We also tried an ontology browser UI that allowed
to browser multiple parents for a "current node", and we found that this UI was too confusing.

```gql
query OverviewPage {
  studyOverviewsList {
    studyId
    studyName
    description
    tissueLabels
    tissueNcitIds
    tissueParentIds
    diseaseLabels
    diseaseMeshIds
    diseaseParentIds
  }
  treeTissuesList {
    cid
    label
    ontCode
    parentCids
    parentOntCodePath
  }
  treeDiseasesList {
    cid
    label
    ontCode
    parentCids
    parentOntCodePath
  }
}
```

### Autocomplete Search

In the search bar, users can enter disease, tissue, cell type concept labels and synonyms. The API
presents a couple of matching concepts for the typed characters. For selected concepts, the UI can
keep the Ids and filter the study list (using Ids and ParentIds).


```gql
query Autocomplete {
  autocompleteList(searchQuery: "pan") {
    ontology
    ontCode
    label
    labelHighlight
    isSynonymOfPreferredTerm
  }
}
```

TODO: we can also restrict the search results on actual study annotation data, as in the tree structures above.
