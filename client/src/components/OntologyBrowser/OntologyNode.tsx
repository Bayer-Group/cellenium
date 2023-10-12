import { ActionIcon, createStyles, Group, Text, useMantineTheme } from '@mantine/core';
import { IconChevronDown, IconChevronRight } from '@tabler/icons-react';
import { useMemo } from 'react';
import { OntologyItem } from '../../model';
import { StudyOverview } from '../../generated/types';

const useStyles = createStyles((theme) => ({
  main: {
    paddingLeft: '2px',
    paddingRight: '2px',
    '&:hover': {
      backgroundColor: theme.colors.gray[1],
      borderRadius: '2px',
    },
  },
  cursor: {
    cursor: 'pointer',
  },
}));

export function OntologyNode({
  item,
  selected,
  hasChildren,
  level,
  onToggle,
  handleAddOntologyItem,
  studies,
}: {
  item: OntologyItem;
  selected: boolean;
  hasChildren: boolean | undefined;
  level: number;
  onToggle: () => void;
  handleAddOntologyItem: (item: OntologyItem) => void;
  studies?: StudyOverview & { allOntCodes: string[] }[];
}) {
  const theme = useMantineTheme();
  const { classes } = useStyles();

  const numberOfStudies = useMemo(() => {
    return studies?.filter((study) => study.allOntCodes.includes(item.id)).length;
  }, [item.id, studies]);

  return (
    <Group pl={`${level * 10}px`} spacing={0}>
      {hasChildren && (
        <ActionIcon size="xs" variant="subtle" onClick={onToggle}>
          {selected ? <IconChevronDown color={theme.colors.dark[9]} /> : <IconChevronRight color={theme.colors.dark[9]} />}
        </ActionIcon>
      )}
      <Text onClick={() => handleAddOntologyItem(item)} className={classes.main} size="xs" classNames={classes.cursor}>
        {item.label} ({numberOfStudies})
      </Text>
    </Group>
  );
}
