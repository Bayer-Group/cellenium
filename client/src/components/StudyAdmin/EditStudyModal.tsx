import { useForm } from '@mantine/form';
import { showNotification } from '@mantine/notifications';
import { useEffect } from 'react';
import { Button, Checkbox, Group, Modal, Stack, Text, Textarea, TextInput } from '@mantine/core';
import { Form } from 'react-router-dom';
import { InputMaybe, StudyAdminDetailsFragment, useStudyDefinitionUpdateMutation, useStudyUpdateMutation } from '../../generated/types';

export function EditStudyModal({ opened, reset, study }: { opened: boolean; reset: () => void; study: StudyAdminDetailsFragment | undefined }) {
  const [studyUpdateMutation, { loading: studyUpdateLoading }] = useStudyUpdateMutation();
  const [studyDefinitionUpdateMutation, { loading: studyUpdateLoading2 }] = useStudyDefinitionUpdateMutation();

  const form = useForm({
    initialValues: {
      studyName: '',
      description: '',
      readerPermissions: '',
      adminPermissions: '',
      tissueNcitIds: '',
      diseaseMeshIds: '',
      visible: false,
      externalWebsite: '',
    },
    validate: {},
  });

  const submit = () => {
    if (study) {
      const splitToArray = (s: string) => {
        const a = s.split(';').map((si) => si.trim());
        if (a.length === 1 && a[0] === '') {
          return null;
        }
        return a;
      };

      studyUpdateMutation({
        variables: {
          studyId: study.studyId,
          studyName: form.values.studyName,
          description: form.values.description,
          readerPermissions: splitToArray(form.values.readerPermissions) as InputMaybe<string[]>,
          adminPermissions: splitToArray(form.values.adminPermissions) as InputMaybe<string[]>,
          tissueNcitIds: splitToArray(form.values.tissueNcitIds) as InputMaybe<string[]>,
          diseaseMeshIds: splitToArray(form.values.diseaseMeshIds) as InputMaybe<string[]>,
          visible: form.values.visible,
          externalWebsite: form.values.externalWebsite,
        },
      })
        .then(() => {
          studyDefinitionUpdateMutation().then(() => {
            reset();
          });
        })
        .catch((reason) => {
          showNotification({
            title: 'Could not save study changes',
            message: reason.message,
            color: 'red',
          });
        });
    }
  };

  useEffect(() => {
    form.setValues({
      studyName: study?.studyName || '',
      description: study?.description || '',
      readerPermissions: (study?.readerPermissions || []).join('; '),
      adminPermissions: (study?.adminPermissions || []).join('; '),
      tissueNcitIds: (study?.tissueNcitIds || []).join('; '),
      diseaseMeshIds: (study?.diseaseMeshIds || []).join('; '),
      visible: study?.visible || false,
      externalWebsite: study?.externalWebsite || '',
    });
  }, [study]);

  return (
    <Modal
      opened={opened}
      onClose={() => {
        reset();
      }}
      size="80vw"
    >
      <Stack>
        <Text weight="bold" size="xl">
          Edit Study
        </Text>
        <Form>
          <TextInput label="Title" {...form.getInputProps('studyName')} />
          <Textarea label="Description" {...form.getInputProps('description')} />
          <TextInput label="Reader Permissions, separate multiple groups / usernames with ;" {...form.getInputProps('readerPermissions')} />
          <TextInput label="Admin Permissions, separate multiple groups / usernames with ;" {...form.getInputProps('adminPermissions')} />
          <Checkbox mt="md" label="Study is visible" {...form.getInputProps('visible', { type: 'checkbox' })} />
          <TextInput label="Tissue NCIT IDs, separate multiple with ;" {...form.getInputProps('tissueNcitIds')} />
          <TextInput
            label="Disease MeSH IDs, separate multiple with ; and use the pseudo-ID HEALTHY to indicate 'not diseased'"
            {...form.getInputProps('diseaseMeshIds')}
          />
          <TextInput label="External Website" {...form.getInputProps('externalWebsite')} />
          <Group position="right" mt="md">
            <Button disabled={!study?.adminPermissionGranted} onClick={submit} loading={studyUpdateLoading || studyUpdateLoading2}>
              Save Changes
            </Button>
          </Group>
        </Form>
      </Stack>
    </Modal>
  );
}
