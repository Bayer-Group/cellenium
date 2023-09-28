import { ActionIcon, createStyles, Group, Text, useMantineTheme } from '@mantine/core';
import { IconChevronDown, IconChevronRight } from '@tabler/icons-react';
import { OntologyItem } from '../../model';

const useStyles = createStyles((theme) => ({
  main: {
    paddingLeft: '2px',
    paddingRight: '2px',
    '&:hover': {
      backgroundColor: theme.colors.gray[1],
      borderRadius: '2px',
    },
  },
}));

export function OntologyNode({
  item,
  selected,
  hasChildren,
  level,
  onToggle,
  handleAddOntologyItem,
}: {
  item: OntologyItem;
  selected: boolean;
  hasChildren: boolean | undefined;
  level: number;
  onToggle: () => void;
  handleAddOntologyItem: (item: OntologyItem) => void;
}) {
  const theme = useMantineTheme();
  const { classes } = useStyles();
  return (
    <Group pl={`${level * 10}px`} spacing={0}>
      {hasChildren && (
        <ActionIcon size="xs" variant="subtle" onClick={onToggle}>
          {selected ? <IconChevronDown color={theme.colors.dark[9]} /> : <IconChevronRight color={theme.colors.dark[9]} />}
        </ActionIcon>
      )}
      <Text onClick={() => handleAddOntologyItem(item)} className={classes.main} size="xs" style={{ cursor: 'pointer' }}>
        {item.label}
      </Text>
    </Group>
  );
}
