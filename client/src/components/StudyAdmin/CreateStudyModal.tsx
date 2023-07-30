import { useForm } from '@mantine/form';
import { useCreateStudyUploadMutation } from '../../generated/types.ts';
import { useCallback } from 'react';
import { showNotification } from '@mantine/notifications';
import { Button, Group, Modal, Select, Stack, Text, TextInput } from '@mantine/core';
import { Form } from 'react-router-dom';
import { Prism } from '@mantine/prism';

export function CreateStudyModal({ opened, reset }: { opened: boolean; reset: () => void }) {
  const form = useForm({
    initialValues: {
      studyName: '',
      filetype: '.h5ad',
    },
    validate: (values) => {
      const errors: Record<string, string> = {};
      if (values.studyName === '') {
        errors.studyName = 'Study name is required';
      }
      return errors;
    },
  });

  const [createStudyUploadMutation, { data: studyUploadData, loading, error }] = useCreateStudyUploadMutation();

  const createStudy = useCallback(() => {
    createStudyUploadMutation({
      variables: {
        studyName: form.values.studyName,
        filetype: form.values.filetype,
      },
    }).catch((reason: any) => {
      showNotification({
        title: 'Could not create study',
        message: reason.message,
        color: 'red',
      });
      form.reset();
      reset();
    });
  }, [form, createStudyUploadMutation, reset]);

  return (
    <Modal
      opened={opened}
      onClose={() => {
        reset();
        form.reset();
      }}
      size="xl"
    >
      <Stack>
        <Text weight="bold" size="xl">
          Create Study
        </Text>
        <Form>
          <TextInput label="Study Title" {...form.getInputProps('studyName')} disabled={loading || studyUploadData !== undefined || error !== undefined} />
          <Select
            data={['.h5ad', '.h5mu']}
            label="Filetype"
            placeholder="Select a filetype"
            {...form.getInputProps('filetype')}
            disabled={loading || studyUploadData !== undefined || error !== undefined}
          />
        </Form>
        {studyUploadData !== undefined && (
          <>
            <Text>Please upload your study file using the following curl command (replace &lt;PATH_TO_YOUR_FILE&gt; with the path to your file):</Text>
            <Prism language="bash" copyLabel="Command code to clipboard" copiedLabel="Command copied to clipboard">
              {`curl -v ${Object.entries(studyUploadData?.createStudyUpload.json['fields'])
                .map(([key, value]) => '-F ' + key + '=' + value)
                .join(' ')} -F file=@<PATH_TO_YOUR_FILE> ${studyUploadData?.createStudyUpload.json['url']}`}
            </Prism>
          </>
        )}

        <Group position="right">
          <Button
            color="blue"
            onClick={createStudy}
            loading={loading}
            disabled={!form.isValid() || loading || studyUploadData !== undefined || error !== undefined}
          >
            Create
          </Button>
        </Group>
      </Stack>
    </Modal>
  );
}
