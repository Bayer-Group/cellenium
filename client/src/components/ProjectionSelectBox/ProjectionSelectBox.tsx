import { Select, Text } from '@mantine/core';
import { useRecoilState, useRecoilValue } from 'recoil';
import { useCallback, useMemo } from 'react';
import { selectedProjectionState, studyState } from '../../atoms';

function niceLabel(value: string) {
  return value.replace('umap', 'UMAP').replace('tsne', 't-SNE').replace('pca', 'PCA');
}

export function ProjectionSelectBox() {
  const study = useRecoilValue(studyState);
  const [projection, setProjection] = useRecoilState(selectedProjectionState);
  const options = useMemo(
    () =>
      Array.from(study?.projections || []).map((p) => ({
        value: p,
        label: niceLabel(p),
      })),
    [study?.projections],
  );

  const onChange = useCallback(
    (p: string | null) => {
      setProjection(p || '');
    },
    [setProjection],
  );

  if (options.length === 1) {
    return <Text size="xs">Projection: {options[0].label}</Text>;
  }

  return <Select w="100%" labelProps={{ size: 'xs' }} label="Select projection" value={projection} onChange={onChange} data={options} />;
}
