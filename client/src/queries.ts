import {gql} from '@apollo/client';

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

`;