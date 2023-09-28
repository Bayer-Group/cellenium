import { useEffect, useMemo, useState } from 'react';
import { Group, Loader } from '@mantine/core';
import { StudyInfoFragment, useStudiesQuery } from '../../generated/types';
import { OntologyItem } from '../../model';
import { generateOntologyTrees } from '../../utils/helper';
import { OfferingItem, SearchBar } from './SearchBar';

function metadataValues(study: StudyInfoFragment) {
  if (study.metadata) {
    return Object.keys(study.metadata).map((key) => study.metadata[key]);
  }
  return [];
}

export function StudySearchBar({ onStudyListUpdate }: { onStudyListUpdate: (studies: StudyInfoFragment[]) => void }) {
  const { data, loading } = useStudiesQuery();
  const allStudies = useMemo(
    () =>
      data?.studyOverviewsList &&
      data.studyOverviewsList.map((study) => ({
        ...study,
        allOntCodes: [study.studyOntologyList.map((ont) => ont.ontCodes), study.studyOntologyList.map((ont) => ont.parentIds)]
          .flat(2)
          .filter((ontCode) => !!ontCode),
        allSearchableStrings: [study.studyName, study.description, ...metadataValues(study)]
          .filter((s) => s?.length > 0)
          .map((s) => `${s}`.toLocaleLowerCase()),
      })),
    [data],
  );

  const ontologyTrees: Map<string, OntologyItem> | undefined = useMemo(
    () => data?.treeOntologiesList && generateOntologyTrees(data.treeOntologiesList),
    [data],
  );
  const [filters, setFilters] = useState<OfferingItem[]>([]);

  const filteredStudies = useMemo(() => {
    if (!allStudies) {
      return undefined;
    }
    const groupedFilters = filters.reduce((entryMap, e) => entryMap.set(e.ontology, [...(entryMap.get(e.ontology) || []), e.ontcode || e.value]), new Map());
    let keepStudies = allStudies;
    groupedFilters.forEach((values: string[], ontology: string) => {
      if (ontology === 'FREETEXT') {
        const valuesLower = values.map((v) => v.toLocaleLowerCase());
        keepStudies = keepStudies.filter((s) => s.allSearchableStrings.find((str) => valuesLower.find((v) => str.includes(v))));
      } else {
        // all 'real' ontology annotations
        keepStudies = keepStudies.filter((s) => s.allOntCodes.find((o) => values.indexOf(o) !== -1));
      }
    });
    return keepStudies;
  }, [allStudies, filters]);

  useEffect(() => onStudyListUpdate(filteredStudies || []), [filteredStudies]);

  if (loading) {
    return (
      <Group position="center" mt="5rem">
        <Loader variant="dots" />
      </Group>
    );
  }
  return <SearchBar ontologies={ontologyTrees} onSearchFiltersUpdate={setFilters} />;
}
