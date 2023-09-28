import { Card, Container, createStyles, Group, Stack, Text } from '@mantine/core';
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
    <Card pr={0} shadow="sm" p="lg" radius="md" withBorder component={Link} to={newStudyUrl}>
      <Card.Section className={classes.main} withBorder inheritPadding py="xs">
        <Group position="left" spacing="xs" noWrap pr={20}>
          <Stack spacing={4} justify="flex-end" align="flex-start" style={{ minWidth: 90 }}>
            <Text span fw={100} lineClamp={1} size="xs">
              Cell annotation
            </Text>
            <Text fw={100} size="xs">
              Study
            </Text>
          </Stack>
          <Stack spacing={4} justify="flex-start">
            <Text fw={700} lineClamp={1} size="xs">{`${data.annotationValue.displayValue} (${data.annotationValue.annotationGroup.displayGroup})`}</Text>
            <Text fw={700} lineClamp={1} size="xs">
              {data.study.studyName}
            </Text>
          </Stack>
        </Group>
      </Card.Section>
      <Stack pt={10}>
        <Group position="left" spacing="xs" noWrap>
          <Stack align="flex-start">
            <Text fw={100} size="xs">
              Gene
            </Text>
            <Text size="xs" fw={100} lineClamp={1}>
              p-value (adj.)
            </Text>
            <Text size="xs" fw={100}>
              log2FC
            </Text>
          </Stack>
          <Stack>
            <Text size="xs" fw={700}>
              {data.omics.displaySymbol}
            </Text>
            <Text size="xs" fw={700}>
              {data.pvalueAdj.toExponential(2)}
            </Text>
            <Text size="xs" fw={700}>
              {data.log2Foldchange.toFixed(2)}
            </Text>
          </Stack>
          <Container>
            <InlineFoldChangePlot studyId={data.study.studyId} annotationValueId={data.annotationValueId} pval={data.pvalueAdj} log2fc={data.log2Foldchange} />
          </Container>
        </Group>
      </Stack>
    </Card>
  );
}
