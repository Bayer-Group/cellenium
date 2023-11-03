import React, { forwardRef, useCallback, useMemo, useState } from 'react';
import { Group, Select, Text } from '@mantine/core';
import { RecoilState, useRecoilState, useRecoilValue, useSetRecoilState } from 'recoil';
import { IconCalculator } from '@tabler/icons-react';
import { annotationGroupIdState, annotationSecondaryGroupIdState, selectedAnnotationState, studyState } from '../../atoms';

interface ItemProps extends React.ComponentPropsWithoutRef<'div'> {
  value: string;
  label: string;
  differentialExpressionCalculated: boolean;
}

const SelectItem = forwardRef<HTMLDivElement, ItemProps>(function SelectItem({ label, differentialExpressionCalculated, ...others }: ItemProps, ref) {
  return (
    <div ref={ref} {...others}>
      <Group position="apart" align="center" noWrap>
        <Text>{label}</Text>
        {differentialExpressionCalculated && (
          <Group title="DEG calculated" align="center">
            <IconCalculator color="gray" size={20} />
          </Group>
        )}
      </Group>
    </div>
  );
});

function AnnotationGroupSelector({
  label,
  annotationGroupState,
  hideAnnotationGroup,
  includeDeselectEntry,
  onSelect,
}: {
  label?: string;
  annotationGroupState: RecoilState<number | undefined>;
  hideAnnotationGroup?: number;
  includeDeselectEntry?: boolean;
  onSelect?: (annotationGroupId: number | undefined) => void;
}) {
  const study = useRecoilValue(studyState);
  const [annotationGroupId, setAnnotationGroupId] = useRecoilState(annotationGroupState);
  const [value, setValue] = useState<string | undefined>();
  const annotations: ItemProps[] = useMemo(() => {
    const anns: ItemProps[] = [];
    if (study) {
      if (includeDeselectEntry) {
        anns.push({
          value: 'NONE',
          label: '(none)',
          differentialExpressionCalculated: false,
        });
      }
      study.annotationGroupMap.forEach((v, key) => {
        if (key !== hideAnnotationGroup) {
          anns.push({
            value: key.toString(),
            label: v.displayGroup,
            differentialExpressionCalculated: v.differentialExpressionCalculated,
          });
        }
      });
    }
    setValue(annotationGroupId?.toString() || 'NONE');
    return anns;
  }, [study, annotationGroupId, includeDeselectEntry, hideAnnotationGroup]);

  const update = useCallback(
    (v: string | null) => {
      if (v) {
        const intVal = v === 'NONE' ? undefined : parseInt(v, 10);
        setValue(v);
        setAnnotationGroupId(intVal);
        onSelect && onSelect(intVal);
      }
    },
    [onSelect, setAnnotationGroupId],
  );

  return (
    <Select
      w="100%"
      value={value}
      onChange={update}
      label={label || 'Select annotation group'}
      labelProps={{ size: 'xs' }}
      itemComponent={SelectItem}
      placeholder="Pick one"
      transitionProps={{
        duration: 80,
        timingFunction: 'ease',
      }}
      data={annotations}
    />
  );
}

export function AnnotationGroupSelectBox() {
  const setSelectedAnnotation = useSetRecoilState(selectedAnnotationState);
  const [annotationSecondaryGroupId, setAnnotationSecondaryGroupId] = useRecoilState(annotationSecondaryGroupIdState);
  return (
    <AnnotationGroupSelector
      annotationGroupState={annotationGroupIdState}
      onSelect={(annotationGroupId) => {
        setSelectedAnnotation(0);
        if (annotationGroupId === annotationSecondaryGroupId) {
          setAnnotationSecondaryGroupId(undefined);
        }
      }}
    />
  );
}

export function AnnotationSecondGroupSelectBox() {
  const annotationGroupId = useRecoilValue(annotationGroupIdState);
  return (
    <AnnotationGroupSelector
      label="Group by second annotation"
      annotationGroupState={annotationSecondaryGroupIdState}
      hideAnnotationGroup={annotationGroupId}
      includeDeselectEntry
    />
  );
}
