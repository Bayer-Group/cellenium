import { ActionIcon, Group, MultiSelect, SelectItem, Stack, Text } from '@mantine/core';
import { IconTrash } from '@tabler/icons-react';
import { useRecoilState, useRecoilValue } from 'recoil';
import { useCallback, useEffect, useMemo } from 'react';
import { selectedAnnotationFilterState, studyState } from '../../atoms';

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
      const anns: number[] = [];
      study.annotationGroupMap.forEach((g) => {
        g.annotationValuesList.forEach((ann) => {
          if (!selectedAnnotationFilter.includes(ann.annotationValueId)) {
            anns.push(ann.annotationValueId);
          }
        });
      });
      setSelectedAnnotationFilter(anns);
    }
  }, [selectedAnnotationFilter, setSelectedAnnotationFilter, study]);

  return (
    <Stack spacing={0}>
      <Group align="center" position="apart" noWrap>
        <Text size="xs">Filter out cells with annotation(s)</Text>
        <Group align="center" noWrap>
          <ActionIcon size="xs" aria-label="Invert Selection" onClick={invertSelection}>
            <svg xmlns="http://www.w3.org/2000/svg" height="1em" viewBox="0 0 576 512" style={{ fill: 'gray' }}>
              <path
                color="currentColor"
                d="M272 416c17.7 0 32-14.3 32-32s-14.3-32-32-32H160c-17.7 0-32-14.3-32-32V192h32c12.9 0 24.6-7.8 29.6-19.8s2.2-25.7-6.9-34.9l-64-64c-12.5-12.5-32.8-12.5-45.3 0l-64 64c-9.2 9.2-11.9 22.9-6.9 34.9s16.6 19.8 29.6 19.8l32 0 0 128c0 53 43 96 96 96H272zM304 96c-17.7 0-32 14.3-32 32s14.3 32 32 32l112 0c17.7 0 32 14.3 32 32l0 128H416c-12.9 0-24.6 7.8-29.6 19.8s-2.2 25.7 6.9 34.9l64 64c12.5 12.5 32.8 12.5 45.3 0l64-64c9.2-9.2 11.9-22.9 6.9-34.9s-16.6-19.8-29.6-19.8l-32 0V192c0-53-43-96-96-96L304 96z"
              />
            </svg>
          </ActionIcon>
          <ActionIcon
            size="xs"
            onClick={() => {
              setSelectedAnnotationFilter([]);
            }}
          >
            <IconTrash />
          </ActionIcon>
        </Group>
      </Group>
      <Group>
        <MultiSelect
          style={{ width: 210, maxWidth: 210, minWidth: 210 }}
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
