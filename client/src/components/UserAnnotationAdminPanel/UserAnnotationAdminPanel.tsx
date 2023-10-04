import { useState } from 'react';
import { Button, Divider, Stack, Switch, Text } from '@mantine/core';
import { useRecoilValue, useSetRecoilState } from 'recoil';
import { annotationGroupIdState, studyReloadHelperState, studyState } from '../../atoms.ts';
import { useDeleteUserAnnotationMutation, useEditUserAnnotationMutation } from '../../generated/types.ts';

export function UserAnnotationAdminPanel() {
  const study = useRecoilValue(studyState);
  const annotationGroupId = useRecoilValue(annotationGroupIdState);
  const annotationGroup = study?.annotationGroupMap.get(annotationGroupId || -1);
  const [editUserAnnotationMutation] = useEditUserAnnotationMutation();
  const [deleteUserAnnotationMutation, { loading: deleteLoading }] = useDeleteUserAnnotationMutation();
  const setStudyReloadHelper = useSetRecoilState(studyReloadHelperState);
  const [reallyDelete, setReallyDelete] = useState(false);
  return (
    <>
      <Divider size="xs" label="Manage user annotation" />
      {annotationGroup?.createdByUser && (
        <Stack justify="flex-start" align={'flex-start'}>
          {!annotationGroup.currentUserIsOwner && <Text>Annotation was created by {annotationGroup.createdByUser}</Text>}
          {annotationGroup.currentUserIsOwner && (
            <>
              <Switch
                size="xs"
                checked={!annotationGroup.privateToUser}
                onChange={(event) =>
                  editUserAnnotationMutation({
                    variables: {
                      studyId: study?.studyId || -1,
                      annotationGroupId: annotationGroupId || -1,
                      privateToUser: !event.currentTarget.checked,
                    },
                  }).then(() => {
                    setStudyReloadHelper((prev) => prev + 1);
                  })
                }
                label={'visible to all users'}
              />

              <Button
                color={reallyDelete ? 'red' : ''}
                variant={reallyDelete ? 'filled' : 'light'}
                size="xs"
                loading={deleteLoading}
                onClick={() => {
                  if (!reallyDelete) {
                    setReallyDelete(true);
                  } else {
                    deleteUserAnnotationMutation({
                      variables: {
                        studyId: study?.studyId || -1,
                        annotationGroupId: annotationGroupId || -1,
                      },
                    }).then(() => {
                      setStudyReloadHelper((prev) => prev + 1);
                    });
                  }
                }}
              >
                {reallyDelete ? 'Really? Click again to delete' : 'Delete annotation group'}
              </Button>
            </>
          )}
        </Stack>
      )}
    </>
  );
}
