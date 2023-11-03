import { createStyles, Group, Loader, Modal, Stack, Text } from '@mantine/core';
import { Prism } from '@mantine/prism';
import { StudyAdminDetailsFragment, useStudyLogsQuery } from '../../generated/types';

const useStyles = createStyles(() => ({
  height: {
    maxHeight: '70vh',
  },
}));

export function StudyLogModal({ opened, reset, study }: { opened: boolean; reset: () => void; study: StudyAdminDetailsFragment | undefined }) {
  const { classes } = useStyles();
  const { data, loading } = useStudyLogsQuery({
    variables: { studyId: study?.studyId ?? -1 },
    skip: study === undefined || study.studyId === undefined,
  });

  return (
    <Modal opened={opened} onClose={reset} size="xl">
      <Stack>
        <Text weight="bold" size="xl">
          Study Import Log {study?.studyName}
        </Text>
        {loading && (
          <Group position="center">
            <Loader color="blue" />
          </Group>
        )}
        {data && (
          <Prism language="bash" copyLabel="Copy log to clipboard" copiedLabel="Log copied to clipboard" className={classes.height}>
            {data.studyImportLogsList.map((log) => log.importLog).join('\n')}
          </Prism>
        )}
      </Stack>
    </Modal>
  );
}
