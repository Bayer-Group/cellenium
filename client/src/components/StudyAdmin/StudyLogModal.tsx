import { Group, Loader, Modal, Stack, Text } from '@mantine/core';
import { Prism } from '@mantine/prism';
import { StudyAdminDetailsFragment, useStudyLogsQuery } from '../../generated/types';

export function StudyLogModal({ opened, reset, study }: { opened: boolean; reset: () => void; study: StudyAdminDetailsFragment | undefined }) {
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
            <Loader />
          </Group>
        )}
        {data && (
          <Prism language="bash" copyLabel="Copy log to clipboard" copiedLabel="Log copied to clipboard" style={{ maxHeight: '70vh' }}>
            {data.studyImportLogsList.map((log) => log.importLog).join('\n')}
          </Prism>
        )}
      </Stack>
    </Modal>
  );
}
