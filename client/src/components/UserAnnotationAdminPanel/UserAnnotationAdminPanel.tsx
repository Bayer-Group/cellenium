import { ChangeEventHandler, useCallback, useState } from 'react';
import { Button, Divider, Stack, Switch, Text } from '@mantine/core';
import { useRecoilValue, useSetRecoilState } from 'recoil';
import { annotationGroupIdState, studyReloadHelperState, studyState } from '../../atoms';
import { useDeleteUserAnnotationMutation, useEditUserAnnotationMutation } from '../../generated/types';

export function UserAnnotationAdminPanel() {
  const study = useRecoilValue(studyState);
  const annotationGroupId = useRecoilValue(annotationGroupIdState);
  const annotationGroup = study?.annotationGroupMap.get(annotationGroupId || -1);
  const [editUserAnnotationMutation] = useEditUserAnnotationMutation();
  const [deleteUserAnnotationMutation, { loading: deleteLoading }] = useDeleteUserAnnotationMutation();
  const setStudyReloadHelper = useSetRecoilState(studyReloadHelperState);
  const [reallyDelete, setReallyDelete] = useState(false);

  const deleteAnnotation = useCallback(() => {
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
  }, [annotationGroupId, deleteUserAnnotationMutation, reallyDelete, setStudyReloadHelper, study?.studyId]);

  const onChange: ChangeEventHandler<HTMLInputElement> = useCallback(
    (event) =>
      editUserAnnotationMutation({
        variables: {
          studyId: study?.studyId || -1,
          annotationGroupId: annotationGroupId || -1,
          privateToUser: !event.currentTarget.checked,
        },
      }).then(() => {
        setStudyReloadHelper((prev) => prev + 1);
      }),
    [annotationGroupId, editUserAnnotationMutation, setStudyReloadHelper, study?.studyId],
  );

  return (
    <>
      <Divider size="xs" label="Manage user annotation" />
      {annotationGroup?.createdByUser && (
        <Stack justify="flex-start" align="flex-start">
          {!annotationGroup.currentUserIsOwner && <Text>Annotation was created by {annotationGroup.createdByUser}</Text>}
          {annotationGroup.currentUserIsOwner && (
            <>
              <Switch size="xs" checked={!annotationGroup.privateToUser} onChange={onChange} label="visible to all users" />
              <Button
                color={reallyDelete ? 'red' : ''}
                variant={reallyDelete ? 'filled' : 'light'}
                size="xs"
                loading={deleteLoading}
                onClick={deleteAnnotation}
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
