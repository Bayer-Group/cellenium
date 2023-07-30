import { StudyAdminDetailsFragment, useStudyDeleteMutation } from '../../generated/types.ts';
import { useCallback, useState } from 'react';
import { showNotification } from '@mantine/notifications';
import { Button, Group, Modal, Stack, Text, TextInput } from '@mantine/core';

export function DeleteStudyModal({ opened, reset, study }: { opened: boolean; reset: () => void; study: StudyAdminDetailsFragment | undefined }) {
  const [confirmationInput, setConfirmationInput] = useState('');
  const [deleteStudyMutation, { loading }] = useStudyDeleteMutation();

  const deleteStudy = useCallback(() => {
    if (study === undefined) {
      return;
    }
    // TODO: delete s3 file
    deleteStudyMutation({
      variables: {
        studyId: study.studyId,
      },
    })
      .then(() => {
        setConfirmationInput('');
        reset();
      })
      .catch((reason) => {
        showNotification({
          title: 'Could not delete study',
          message: reason.message,
          color: 'red',
        });
      });
  }, [study, deleteStudyMutation, reset]);

  return (
    <Modal opened={opened} onClose={reset}>
      <Stack>
        <Text weight="bold" size="xl">
          Delete Study
        </Text>
        <Text>
          Are you sure you want to delete the study <b>{study?.studyName}</b>?
        </Text>
        <Text>
          To confirm deletion, type <i>permanently delete</i> in the text input field.
        </Text>
        <TextInput value={confirmationInput} onChange={(event) => setConfirmationInput(event.currentTarget.value)} placeholder="permanently delete" />
        <Group position="right">
          <Button color="red" loading={loading} disabled={confirmationInput !== 'permanently delete' || loading} onClick={deleteStudy}>
            Delete
          </Button>
        </Group>
      </Stack>
    </Modal>
  );
}
