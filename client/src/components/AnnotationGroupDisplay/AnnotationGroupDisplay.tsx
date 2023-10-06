import { Stack } from '@mantine/core';
import { useRecoilState, useRecoilValue } from 'recoil';
import { Annotation } from '../Annotation/Annotation';
import { annotationGroupIdState, highlightAnnotationState, studyState } from '../../atoms';

export function AnnotationGroupDisplay({ disableSelection }: { disableSelection?: boolean }) {
  const [, setHighlightAnnotation] = useRecoilState(highlightAnnotationState);
  const annotationGroupId = useRecoilValue(annotationGroupIdState);

  const study = useRecoilValue(studyState);

  if (!study || !annotationGroupId) {
    return null;
  }
  const annotations = study.annotationGroupMap.get(annotationGroupId)?.annotationValuesList;
  const isSelectable = disableSelection ? false : (study.annotationGroupMap.get(annotationGroupId)?.differentialExpressionCalculated as boolean);
  return (
    <Stack spacing={2} onMouseLeave={() => setHighlightAnnotation(0)} style={{ maxWidth: 205 }}>
      {annotations !== undefined &&
        annotations.map((annot) => {
          return (
            <Annotation
              key={annot.annotationValueId}
              label={annot.displayValue}
              sampleCount={annot.sampleCount}
              color={annot.color}
              annotationId={annot.annotationValueId}
              isSelectable={isSelectable}
              sampleCountPercentage={((annot.sampleCount / study.cellCount) * 100).toFixed(2)}
            />
          );
        })}
    </Stack>
  );
}
