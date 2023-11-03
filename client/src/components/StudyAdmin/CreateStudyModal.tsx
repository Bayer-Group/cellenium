import { useForm } from '@mantine/form';
import { useCallback } from 'react';
import { showNotification } from '@mantine/notifications';
import { Button, Group, Modal, Stack, Text, TextInput } from '@mantine/core';
import { Form } from 'react-router-dom';
import { Prism } from '@mantine/prism';
import { useCreateStudyUploadMutation } from '../../generated/types';

export function CreateStudyModal({ opened, reset }: { opened: boolean; reset: () => void }) {
  const form = useForm({
    initialValues: {
      filename: '',
    },
    validate: (values) => {
      const errors: Record<string, string> = {};
      if (values.filename === '') {
        errors.filename = 'Filename is required';
      } else if (!values.filename.endsWith('.h5ad') && !values.filename.endsWith('.h5mu')) {
        errors.filename = 'Only .h5ad or .h5mu files are accepted';
      }
      return errors;
    },
  });

  const [createStudyUploadMutation, { data: studyUploadData, loading, error }] = useCreateStudyUploadMutation();

  const createStudy = useCallback(() => {
    createStudyUploadMutation({
      variables: {
        filename: form.values.filename,
      },
    }).catch((reason: { message: string }) => {
      showNotification({
        title: 'Could not create study',
        message: reason.message,
        color: 'red',
      });
      form.reset();
      reset();
    });
  }, [form, createStudyUploadMutation, reset]);

  const modalClose = useCallback(() => {
    reset();
    form.reset();
  }, [form, reset]);

  return (
    <Modal opened={opened} onClose={modalClose} size="xl">
      <Stack>
        <Text weight="bold" size="xl">
          Create Study
        </Text>
        <Text>
          Cellenium can import studies in h5ad or h5mu file format. Enter your local file name to generate create a cellenium study placeholder (to be filled
          with data from the file) and a curl command for uploading your file to cellenium.
        </Text>
        <Form>
          <TextInput
            label="Filename"
            {...form.getInputProps('filename')}
            placeholder="my_study.h5ad"
            disabled={loading || studyUploadData !== undefined || error !== undefined}
          />
        </Form>
        {studyUploadData !== undefined && (
          <>
            <Text>
              Please upload your study file using the following curl command. (The local filename is found after the @ symbol in the proposed curl command, feel
              free to modify the actual local name and path.)
            </Text>
            <Prism language="bash" copyLabel="Command code to clipboard" copiedLabel="Command copied to clipboard">
              {`curl -v ${Object.entries(studyUploadData?.createStudyUpload.json.fields)
                .map(([key, value]) => `-F ${key}=${value}`)
                .join(' ')} -F file=@${form.values.filename} ${studyUploadData?.createStudyUpload.json.url}`}
            </Prism>
            <Text>You can close this dialog. Once uploaded, data will be processed automatically. Refresh the study list to observe the current status.</Text>
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
