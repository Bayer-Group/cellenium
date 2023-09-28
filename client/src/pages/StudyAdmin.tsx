import { useCallback, useState } from 'react';
import { NavBar } from '../components';
import { Button, Container, Group, Loader, Space, Stack, Text } from '@mantine/core';
import { IconArticle, IconDotsVertical, IconPlus, IconX } from '@tabler/icons-react';
import { StudyAdminDetailsFragment, useStudyAdminListQuery } from '../generated/types';
import DataTable from 'react-data-table-component';
import { CreateStudyModal } from '../components/StudyAdmin/CreateStudyModal.tsx';
import { DeleteStudyModal } from '../components/StudyAdmin/DeleteStudyModal.tsx';
import { StudyLogModal } from '../components/StudyAdmin/StudyLogModal.tsx';
import { EditStudyModal } from '../components/StudyAdmin/EditStudyModal.tsx';

export default function StudyAdmin() {
  const [newStudyModalOpen, setNewStudyModalOpen] = useState(false);

  const { data, loading, refetch, error } = useStudyAdminListQuery();
  const [selectedEditStudy, setSelectedEditStudy] = useState<StudyAdminDetailsFragment | undefined>(undefined);
  const [selectedDeleteStudy, setSelectedDeleteStudy] = useState<StudyAdminDetailsFragment | undefined>(undefined);
  const [selectedLogStudy, setSelectedLogStudy] = useState<StudyAdminDetailsFragment | undefined>(undefined);

  const resetNewStudyModal = useCallback(() => {
    setNewStudyModalOpen(false);
    void refetch();
  }, [setNewStudyModalOpen, refetch]);

  const resetDeleteModal = useCallback(() => {
    setSelectedDeleteStudy(undefined);
    void refetch();
  }, [setSelectedDeleteStudy, refetch]);

  const resetEditModal = useCallback(() => {
    setSelectedEditStudy(undefined);
    void refetch();
  }, [setSelectedEditStudy, refetch]);

  const columns = [
    {
      name: 'ID',
      selector: (row: StudyAdminDetailsFragment) => row.studyId,
      sortable: true,
      width: '10%',
    },
    {
      name: 'Title',
      selector: (row: StudyAdminDetailsFragment) => row.studyName,
      sortable: true,
    },
    {
      name: 'Filename',
      selector: (row: StudyAdminDetailsFragment) => row.filename,
      width: '20%',
    },
    {
      name: 'Your Role',
      selector: (row: StudyAdminDetailsFragment) => (row.adminPermissionGranted ? 'Admin' : row.readerPermissionGranted ? 'View' : 'No Access'),
      sortable: true,
      width: '10%',
    },
    {
      name: 'Import Status',
      selector: (row: StudyAdminDetailsFragment) =>
        row.importStarted
          ? row.importFailed
            ? 'Failed'
            : row.importFinished
            ? 'Imported'
            : 'Unknown Error'
          : row.importFailed
          ? 'Failed'
          : row.importFinished
          ? 'Imported'
          : 'Not Started',
      sortable: true,
      width: '130px',
    },
    {
      name: '',
      width: '50px',
      cell: (row: StudyAdminDetailsFragment) =>
        !row.hasImportLog ? null : (
          <span title={'View import logs'}>
            <IconArticle onClick={() => setSelectedLogStudy(row)} style={{ cursor: 'pointer' }} />
          </span>
        ),
    },
    {
      name: '',
      width: '50px',
      cell: (row: StudyAdminDetailsFragment) =>
        !row.adminPermissionGranted ? null : <IconDotsVertical onClick={() => setSelectedEditStudy(row)} style={{ cursor: 'pointer' }} />,
    },
    {
      name: '',
      width: '50px',
      cell: (row: StudyAdminDetailsFragment) =>
        !row.adminPermissionGranted ? null : <IconX onClick={() => setSelectedDeleteStudy(row)} style={{ cursor: 'pointer' }} />,
    },
  ];

  return (
    <Container fluid={true}>
      <EditStudyModal opened={selectedEditStudy !== undefined} reset={resetEditModal} study={selectedEditStudy} />
      {newStudyModalOpen && <CreateStudyModal opened={newStudyModalOpen} reset={resetNewStudyModal} />}
      <DeleteStudyModal study={selectedDeleteStudy} reset={resetDeleteModal} opened={selectedDeleteStudy !== undefined} />
      <StudyLogModal opened={selectedLogStudy !== undefined} study={selectedLogStudy} reset={() => setSelectedLogStudy(undefined)} />
      <NavBar />
      <Space h="xl" />
      <Stack px="md">
        <Group position="right">
          {(data?.userStudyUploadConfigured || true) && (
            <Button onClick={() => setNewStudyModalOpen(true)}>
              <Group spacing="xs">
                <IconPlus />
                <span>New Study</span>
              </Group>
            </Button>
          )}
        </Group>

        {loading && (
          <Group position="center">
            <Loader variant="dots" color="blue" size="md" />
          </Group>
        )}
        {!loading && data ? (
          <DataTable data={data?.studyAdminDetailsList || []} columns={columns} defaultSortFieldId={1} defaultSortAsc={false} style={{ width: '100%' }} />
        ) : null}
        {error && (
          <Group position="center">
            <Text weight="bold">An Error occurred</Text>
          </Group>
        )}
      </Stack>
    </Container>
  );
}
