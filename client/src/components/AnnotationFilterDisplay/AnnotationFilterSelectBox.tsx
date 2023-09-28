import { useCallback, useMemo } from 'react';
import { MultiSelect, SelectItem } from '@mantine/core';
import { useRecoilState, useRecoilValue } from 'recoil';
import { selectedAnnotationFilterState, studyState } from '../../atoms';

export function AnnotationFilterSelectBox() {
  const study = useRecoilValue(studyState);
  const [selectedAnnotationFilter, setSelectedAnnotationFilter] = useRecoilState(selectedAnnotationFilterState);
  const annotations: SelectItem[] = useMemo(() => {
    const anns: SelectItem[] = [];
    if (study) {
      study.annotationGroupMap.forEach((value) => {
        const groupName = value.displayGroup;
        value.annotationValuesList.map((av) =>
          anns.push({
            value: av.annotationValueId.toString(),
            label: av.displayValue,
            group: groupName,
          }),
        );
      });
    }
    return anns;
  }, [study]);

  const onChange = useCallback(
    (value: string[]) => {
      setSelectedAnnotationFilter(value.map((v: string) => parseInt(v, 10)));
    },
    [setSelectedAnnotationFilter],
  );

  return (
    <MultiSelect
      style={{ width: 210, maxWidth: 210, minWidth: 210 }}
      value={selectedAnnotationFilter.map((s) => s.toString())}
      onChange={onChange}
      searchable
      placeholder="Select"
      data={annotations}
    />
  );
}
