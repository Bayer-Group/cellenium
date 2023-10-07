import { Card, createStyles, Group, Stack, Text } from '@mantine/core';
import { Link } from 'react-router-dom';
import { DifferentialMarkerFragment } from '../../generated/types';
import { InlineFoldChangePlot } from '../InlineFoldChangePlot/InlineFoldChangePlot';

const useStyles = createStyles((theme) => ({
  main: {
    backgroundColor: 'white',
    cursor: 'pointer',
    '&:hover': {
      backgroundColor: theme.colors.gray[1],
    },
  },
}));

export function MarkerCard({ data }: { data: DifferentialMarkerFragment }) {
  const { classes } = useStyles();
  const newStudyUrl = `/study/${data.study.studyId}?page=CellMarkerAnalysis&annotationGroupId=${data.annotationValue.annotationGroup.annotationGroupId}&annotationValueId=${data.annotationValueId}&omicsId=${data.omics.omicsId}`;
  return (
    <Card w="100%" shadow="sm" p="lg" radius="md" withBorder component={Link} to={newStudyUrl}>
      <Card.Section className={classes.main} withBorder inheritPadding py="xs">
        <Group position="left" spacing="xs" noWrap pr={20}>
          <Stack spacing={4} justify="flex-end" align="flex-start" style={{ minWidth: 90 }}>
            <Text span fw={100} lineClamp={1} size="xs">
              Annotation group
            </Text>
            <Text fw={100} size="xs">
              Annotation
            </Text>
          </Stack>
          <Stack spacing={4} justify="flex-start">
            <Text fw={700} lineClamp={1} size="xs">
              {data.annotationValue.annotationGroup.displayGroup}
            </Text>
            <Text fw={700} lineClamp={1} size="xs">
              {data.annotationValue.displayValue}
            </Text>
          </Stack>
        </Group>
      </Card.Section>
      <Stack pt={10}>
        <Group position="left" spacing="xs" noWrap>
          <Stack align="flex-start" w="6rem">
            <Text fw={100} size="xs">
              Gene
            </Text>
            <Text size="xs" fw={100} lineClamp={1} w="6rem">
              p-value (adj.)
            </Text>
            <Text size="xs" fw={100} w="6rem">
              log2FC
            </Text>
          </Stack>
          <Stack>
            <Text size="xs" fw={700} w="6rem">
              {data.omics.displaySymbol}
            </Text>
            <Text size="xs" fw={700} w="6rem">
              {data.pvalueAdj.toExponential(2)}
            </Text>
            <Text size="xs" fw={700} w="6rem">
              {data.log2Foldchange.toFixed(2)}
            </Text>
          </Stack>
          <Group h="100%" align="center" position="center" w="100%">
            <InlineFoldChangePlot studyId={data.study.studyId} annotationValueId={data.annotationValueId} pval={data.pvalueAdj} log2fc={data.log2Foldchange} />
          </Group>
        </Group>
      </Stack>
    </Card>
  );
}
