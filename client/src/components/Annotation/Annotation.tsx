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
  cursor: {
    cursor: 'pointer',
  },
  border: {
    border: '2px solid black',
  },
  label: {
    whiteSpace: 'nowrap',
  },
}));

export function Annotation({
  label,
  color,
  sampleCount,
  sampleCountPercentage,
  annotationId,
  isSelectable = false,
}: {
  label: string;
  color: string;
  sampleCount: number;
  sampleCountPercentage: string;
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
      columns={15}
      pl="sm"
      gutter={0}
      justify="space-between"
      align="center"
      onMouseOver={() => setHighlight(annotationId)}
      onClick={onClick}
      className={cx(
        {
          [classes.hovered]: annotationId === highlight,
          [classes.clicked]: annotationIsSelected,
        },
        classes.cursor,
      )}
    >
      <Grid.Col span={7}>
        <Group pr={2} spacing={2}>
          <Text title={label} size="xs" weight={showBold} lineClamp={1} className={classes.label}>
            {label}
          </Text>
        </Group>
      </Grid.Col>
      <Grid.Col span={6}>
        {sampleCount ? (
          <Text size="xs" weight={showBold} lineClamp={1} align="right">
            ({sampleCountPercentage ? `${sampleCountPercentage}%` : null})
          </Text>
        ) : null}
      </Grid.Col>
      <Grid.Col span={2} pl={5}>
        <ColorSwatch key={color} color={color} size={annotationIsSelected ? 12 : 15} className={annotationIsSelected ? classes.border : undefined} />
      </Grid.Col>
    </Grid>
  );
}
