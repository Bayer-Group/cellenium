import { TreeOntologyOverviewFragment } from "../generated/types";
import { OntologyItem } from "../model";

export function ontology2Color(ontology: string) {
  console.log({ ontology });
  let color;
  switch (ontology) {
    case "taxonomy":
      color = "blue"; //theme.colors.yellow[5];
      break;
    case "NCIT":
      color = "pink"; //theme.colors.red[5];
      break;
    case "MeSH":
      color = "teal"; //theme.colors.violet[5];
      break;
    case "GENE":
      color = "lime";
      break;
    case "CO":
      color = "orange";
      break;
    default:
      color = "gray"; //theme.colors.gray[5];
      break;
  }
  return color;
}

export function generateOntologyTrees(
  nodeList: TreeOntologyOverviewFragment[],
) {
  // @ts-ignore
  // const ontologies = [...new Set(nodeList.map(item => item.ontology))];
  const ontologyItemMap = new Map<string, OntologyItem>();
  const ontologyMap = new Map<string, OntologyItem>();
  // generate the nodes =
  const allNodes: OntologyItem[] = nodeList.map((nd) => {
    return {
      id: nd.ontCode,
      unique_id: `${nd.ontology}_${nd.ontCode}`,
      label: nd.label,
      parent_unique_id: nd.parentOntCodePath
        ? `${nd.ontology}_${nd.parentOntCodePath[0]}`
        : undefined,
      ontology: nd.ontology,
      children: [],
    };
  });

  // setup the hash
  allNodes.map((nd: OntologyItem) => {
    ontologyItemMap.set(nd.unique_id, nd);
  });

  // fill the children
  allNodes.map((nd: OntologyItem) => {
    if (nd.parent_unique_id) {
      let parent = ontologyItemMap.get(nd.parent_unique_id);
      if (parent && parent.children) parent.children.push(nd);
    }
  });

  // now generate the Map ontology label --> root node of ontology. for each ontology there is just one root node!
  allNodes
    .filter((nd) => nd.parent_unique_id === undefined)
    .forEach((nd) => ontologyMap.set(nd.ontology, nd));
  return ontologyMap;
}
