import { ColorSwatch, createStyles, Group, Text } from '@mantine/core';
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
    outline: '2px solid black',
  },
  spacing_on_hover: {
    letterSpacing: 0,
  },
  nowrap: {
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
    <Group
      w="100%"
      onMouseOver={() => setHighlight(annotationId)}
      onClick={onClick}
      className={cx(
        {
          [classes.hovered]: annotationId === highlight,
          [classes.clicked]: annotationIsSelected,
        },
        classes.cursor,
      )}
      position="apart"
      noWrap
      spacing="xs"
      px="0.25rem"
    >
      <Text title={label} size="xs" weight={showBold} w="7.5rem" truncate classNames={cx(classes.nowrap, showBold ? classes.spacing_on_hover : undefined)}>
        {label}
      </Text>

      <Group position="right" noWrap spacing="xs">
        {sampleCount || sampleCountPercentage ? (
          <Text
            size="xs"
            title={sampleCountPercentage ? `${sampleCountPercentage}%` : sampleCount.toString(10)}
            weight={showBold}
            lineClamp={1}
            align="right"
            classNames={cx(classes.nowrap, showBold ? classes.spacing_on_hover : undefined)}
          >
            ({sampleCountPercentage ? `${sampleCountPercentage}%` : sampleCount.toString(10)})
          </Text>
        ) : null}
        <ColorSwatch key={color} color={color} size={annotationIsSelected ? 12 : 14} className={annotationIsSelected ? classes.border : undefined} />
      </Group>
    </Group>
  );
}
