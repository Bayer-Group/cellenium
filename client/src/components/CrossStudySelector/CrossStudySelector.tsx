import { Center, Loader, Stack, Text } from '@mantine/core';
import { useEffect, useState } from 'react';
import { StudyCard } from '../StudyCard/StudyCard';
import { ReferenceStudyOverview, StudyInfoFragment, useReferenceStudiesOverviewQuery } from '../../generated/types';

export function CrossStudySelector({ selectStudy, taxId }: { selectStudy: (studyId: number) => void; taxId: string }) {
  const [allStudies, setReferenceStudies] = useState<ReferenceStudyOverview[]>([]);
  const { data, loading } = useReferenceStudiesOverviewQuery({ variables: { organismTaxId: taxId } });
  useEffect(() => {
    if (data) {
      setReferenceStudies((data?.referenceStudyOverviewsList as ReferenceStudyOverview[]) || []);
    }
  }, [data, setReferenceStudies]);

  return (
    <Stack w="100%" h="100%" style={{ overflowY: 'scroll' }} pos="relative" pt="xl">
      {loading && (
        <Stack w="100%" h="100%" align="center" justify="center">
          <Loader />
        </Stack>
      )}
      {allStudies.length === 0 && !loading && (
        <Center h="100%" w="100%">
          <Text>No studies found</Text>
        </Center>
      )}
      <Stack px="sm">
        {allStudies.map((study) => (
          <StudyCard
            study={study as StudyInfoFragment}
            key={study.studyId}
            onClick={() => {
              selectStudy(study.studyId);
            }}
          />
        ))}
      </Stack>
    </Stack>
  );
}
