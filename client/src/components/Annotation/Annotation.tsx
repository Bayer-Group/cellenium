import { ColorSwatch, createStyles, Group, Text } from '@mantine/core';
import { useRecoilState, useSetRecoilState } from 'recoil';
import { useCallback } from 'react';
import { highlightAnnotationState, selectedAnnotationState, selectedDEGComparisonAnnotationState, selectedGenesState } from '../../atoms';

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
    userSelect: 'none',
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
  enableDEGComparisonSelection = false,
}: {
  label: string;
  color: string;
  sampleCount: number;
  sampleCountPercentage: string;
  annotationId: number;
  isSelectable: boolean;
  enableDEGComparisonSelection: boolean;
}) {
  const { classes, cx } = useStyles();
  const [highlight, setHighlight] = useRecoilState(highlightAnnotationState);
  const [selected, setSelected] = useRecoilState(selectedAnnotationState);
  const [selectedDEGComparison, setSelectedDEGComparison] = useRecoilState(selectedDEGComparisonAnnotationState);
  const setSelectedGenes = useSetRecoilState(selectedGenesState);
  const annotationIsSelected = (selected === annotationId || selectedDEGComparison === annotationId) && isSelectable;
  const showBold = annotationIsSelected ? 800 : 'md';

  const onClick = useCallback(
    (event: MouseEvent) => {
      if (!isSelectable) return;
      if (highlight === selected) {
        setSelected(0);
      } else {
        setSelectedGenes([]);
        if (event.shiftKey && selected && enableDEGComparisonSelection) {
          setSelectedDEGComparison(annotationId);
        } else {
          setSelected(annotationId);
          setSelectedDEGComparison(0);
        }
      }
    },
    [annotationId, highlight, isSelectable, selected, setSelected, setSelectedDEGComparison, setSelectedGenes],
  );

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
