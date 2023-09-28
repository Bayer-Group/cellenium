import { ColorSwatch, createStyles, Grid, Group, Text } from '@mantine/core';
import { useRecoilState } from 'recoil';
import { useCallback } from 'react';
import { highlightAnnotationState, selectedAnnotationState, selectedGenesState } from '../../atoms';

const useStyles = createStyles((theme) => ({
  hovered: {
    backgroundColor: theme.colors.gray[3],
    borderRadius: theme.radius.xs,
  },
  clicked: {
    backgroundColor: theme.colors.blue[1],
    borderRadius: theme.radius.xs,
  },
}));

export function Annotation({
  label,
  color,
  sampleCount,
  annotationId,
  isSelectable = false,
}: {
  label: string;
  color: string;
  sampleCount: number;
  annotationId: number;
  isSelectable: boolean;
}) {
  const { classes, cx } = useStyles();
  const [highlight, setHighlight] = useRecoilState(highlightAnnotationState);
  const [selected, setSelected] = useRecoilState(selectedAnnotationState);
  const [, setSelectedGenes] = useRecoilState(selectedGenesState);
  const annotationIsSelected = selected === annotationId && isSelectable;
  const showBold = annotationIsSelected ? 800 : 'md';

  const onClick = useCallback(() => {
    if (!isSelectable) return;
    if (highlight === selected) {
      setSelected(0);
    } else {
      setSelectedGenes([]);
      setSelected(annotationId);
    }
  }, [annotationId, highlight, isSelectable, selected, setSelected, setSelectedGenes]);

  return (
    <Grid
      columns={12}
      pl={10}
      gutter={0}
      sx={{ cursor: 'pointer' }}
      justify="space-between"
      align="center"
      onMouseOver={() => setHighlight(annotationId)}
      onClick={onClick}
      className={cx({
        [classes.hovered]: annotationId === highlight,
        [classes.clicked]: annotationIsSelected,
      })}
    >
      <Grid.Col span={7}>
        <Group pr={2} spacing={2}>
          <Text title={label} size="xs" weight={showBold} lineClamp={1}>
            {label}
          </Text>
        </Group>
      </Grid.Col>
      <Grid.Col span={4} style={{ textAlign: 'right' }}>
        {sampleCount ? (
          <Text size="xs" weight={showBold} lineClamp={1}>
            ({sampleCount})
          </Text>
        ) : null}
      </Grid.Col>
      <Grid.Col span={1} pl={5}>
        <div>
          <ColorSwatch key={color} color={color} size={annotationIsSelected ? 12 : 15} style={{ border: annotationIsSelected ? '2px solid black' : '' }} />
        </div>
      </Grid.Col>
    </Grid>
  );
}
