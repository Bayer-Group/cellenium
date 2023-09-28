import { useRecoilState, useRecoilValue } from 'recoil';
import { Stack, Text } from '@mantine/core';
import { selectedAnnotationFilterState, studyState } from '../../atoms';

export function AnnotationFilter() {
  const [selectedAnnotationFilter] = useRecoilState(selectedAnnotationFilterState);
  const study = useRecoilValue(studyState);
  return (
    <Stack spacing={1}>
      {study &&
        study.annotationValueMap &&
        selectedAnnotationFilter.map((av: number) => {
          const annotationValue = study.annotationValueMap.get(av)?.displayValue;
          return (
            <Text size="xs" key={`${annotationValue}-annotation`}>
              {annotationValue}
            </Text>
          );
        })}
    </Stack>
  );
}
