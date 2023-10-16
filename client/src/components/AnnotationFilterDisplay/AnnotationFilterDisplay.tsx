import { ActionIcon, Group, MultiSelect, SelectItem, Stack, Text } from '@mantine/core';
import { IconTrash } from '@tabler/icons-react';
import { useRecoilState, useRecoilValue } from 'recoil';
import { useCallback, useEffect, useMemo } from 'react';
import { selectedAnnotationFilterState, studyState } from '../../atoms';
import { SwitchIcon } from '../../assets/icons/Icons';

export function AnnotationFilterDisplay() {
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

  useEffect(() => {
    setSelectedAnnotationFilter([]);
  }, [setSelectedAnnotationFilter]);

  const invertSelection = useCallback(() => {
    if (study) {
      const groups: number[] = [];
      const anns: number[] = [];
      study.annotationGroupMap.forEach((g) => {
        g.annotationValuesList.forEach((ann) => {
          if (selectedAnnotationFilter.includes(ann.annotationValueId)) {
            groups.push(g.annotationGroupId);
          }
        });
      });
      study.annotationGroupMap.forEach((g) => {
        g.annotationValuesList.forEach((ann) => {
          if (!selectedAnnotationFilter.includes(ann.annotationValueId) && groups.includes(g.annotationGroupId)) {
            anns.push(ann.annotationValueId);
          }
        });
      });
      setSelectedAnnotationFilter(anns);
    }
  }, [selectedAnnotationFilter, setSelectedAnnotationFilter, study]);

  const deleteSelectedAnnotations = useCallback(() => {
    setSelectedAnnotationFilter([]);
  }, [setSelectedAnnotationFilter]);

  return (
    <Stack spacing={0}>
      <Group align="center" position="apart" noWrap>
        <Text size="xs">Filter out cells with annotation(s)</Text>
        <Group align="center" noWrap>
          <ActionIcon size="xs" aria-label="Invert Selection" onClick={invertSelection}>
            <SwitchIcon />
          </ActionIcon>
          <ActionIcon size="xs" onClick={deleteSelectedAnnotations}>
            <IconTrash />
          </ActionIcon>
        </Group>
      </Group>
      <Group>
        <MultiSelect
          w="100%"
          value={selectedAnnotationFilter.map((s) => s.toString())}
          onChange={onChange}
          searchable
          placeholder="Select"
          data={annotations}
        />
      </Group>
    </Stack>
  );
}
