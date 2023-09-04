import React, { forwardRef, useMemo, useState } from 'react';
import { Group, Select, Text } from '@mantine/core';
import { RecoilState, useRecoilState, useRecoilValue, useSetRecoilState } from 'recoil';
import { annotationGroupIdState, annotationSecondaryGroupIdState, selectedAnnotationState, studyState } from '../../atoms';
import { IconCalculator } from '@tabler/icons-react';

interface ItemProps extends React.ComponentPropsWithoutRef<'div'> {
  value: string;
  label: string;
  differentialExpressionCalculated: boolean;
}

const SelectItem = forwardRef<HTMLDivElement, ItemProps>(({ label, differentialExpressionCalculated, ...others }: ItemProps, ref) => (
  <div ref={ref} {...others}>
    <Group position={'apart'} align={'center'} noWrap>
      <Text>{label}</Text>
      {differentialExpressionCalculated && (
        <Group title="DEG calculated" align={'center'}>
          <IconCalculator color={'gray'} size={20} />
        </Group>
      )}
    </Group>
  </div>
));

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
      study.annotationGroupMap.forEach((value, key) => {
        if (key !== hideAnnotationGroup) {
          anns.push({
            value: key.toString(),
            label: value.displayGroup,
            differentialExpressionCalculated: value.differentialExpressionCalculated,
          });
        }
      });
    }
    setValue(annotationGroupId?.toString() || 'NONE');
    return anns;
  }, [study, annotationGroupId]);

  function update(value: string | null) {
    if (value) {
      const intVal = value === null || value === 'NONE' ? undefined : parseInt(value);
      setValue(value);
      setAnnotationGroupId(intVal);
      onSelect && onSelect(intVal);
    }
  }

  return (
    <Select
      style={{ maxWidth: 210, width: 210, minWidth: 210 }}
      value={value}
      onChange={(value) => update(value)}
      label={label || 'Select annotation group'}
      labelProps={{ size: 'xs' }}
      itemComponent={SelectItem}
      placeholder="Pick one"
      transitionProps={{
        duration: 80,
        timingFunction: 'ease',
      }}
      data={annotations as any}
    />
  );
}

function AnnotationGroupSelectBox() {
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
function AnnotationSecondGroupSelectBox() {
  const annotationGroupId = useRecoilValue(annotationGroupIdState);
  return (
    <AnnotationGroupSelector
      label={'Group by second annotation'}
      annotationGroupState={annotationSecondaryGroupIdState}
      hideAnnotationGroup={annotationGroupId}
      includeDeselectEntry
    />
  );
}

export { AnnotationGroupSelectBox, AnnotationSecondGroupSelectBox };
