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

export function StudyAdmin() {
  const [newStudyModalOpen, setNewStudyModalOpen] = useState(false);

  const { data, loading, refetch, error } = useStudyAdminListQuery();
  const [selectedEditStudy, setSelectedEditStudy] = useState<StudyAdminDetailsFragment | undefined>(undefined);
  const [selectedDeleteStudy, setSelectedDeleteStudy] = useState<StudyAdminDetailsFragment | undefined>(undefined);
  const [selectedLogStudy, setSelectedLogStudy] = useState<StudyAdminDetailsFragment | undefined>(undefined);

  const resetDeleteModal = useCallback(() => {
    setSelectedDeleteStudy(undefined);
    void refetch();
  }, []);

  const resetEditModal = useCallback(() => {
    setSelectedEditStudy(undefined);
    void refetch();
  }, []);

  const columns = [
    {
      name: 'ID',
      selector: (row: StudyAdminDetailsFragment) => row.studyId,
      sortable: true,
    },
    {
      name: 'Title',
      selector: (row: StudyAdminDetailsFragment) => row.studyName,
      sortable: true,
    },
    {
      name: 'Filename',
      selector: (row: StudyAdminDetailsFragment) => row.filename,
      sortable: true,
    },
    {
      name: 'Your Role',
      selector: (row: StudyAdminDetailsFragment) => (row.adminPermissionGranted ? 'Admin' : row.readerPermissionGranted ? 'View' : 'No Access'),
      sortable: true,
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
          : row.importFinished
          ? 'Imported'
          : row.importFailed
          ? 'Failed'
          : 'Not Started',
      sortable: true,
    },
    {
      name: '',
      width: '50px',
      cell: (row: StudyAdminDetailsFragment) =>
        !row.hasImportLog ? null : <IconArticle onClick={() => setSelectedLogStudy(row)} style={{ cursor: 'pointer' }} />,
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

      <CreateStudyModal
        opened={newStudyModalOpen}
        reset={() => {
          setNewStudyModalOpen(false);
          void refetch();
        }}
      />
      <DeleteStudyModal study={selectedDeleteStudy} reset={resetDeleteModal} opened={selectedDeleteStudy !== undefined} />
      <StudyLogModal opened={selectedLogStudy !== undefined} study={selectedLogStudy} reset={() => setSelectedLogStudy(undefined)} />
      <NavBar />
      <Space h="xl" />
      <Stack px="md">
        <Group position="right">
          {data?.userStudyUploadConfigured && (
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
            <Loader variant="dots" color="gray" size="xl" />
          </Group>
        )}
        {!loading && data ? <DataTable data={data?.studyAdminDetailsList || []} columns={columns} /> : null}
        {error && (
          <Group position="center">
            <Text weight="bold">An Error occurred</Text>
          </Group>
        )}
      </Stack>
    </Container>
  );
}
