import { useCallback, useMemo, useState } from 'react';
import { Button, createStyles, Group, Loader, Space, Stack, Text } from '@mantine/core';
import { IconArticle, IconDotsVertical, IconPlus, IconX } from '@tabler/icons-react';
import DataTable from 'react-data-table-component';
import { StudyAdminDetailsFragment, useStudyAdminListQuery } from '../generated/types';
import { CreateStudyModal } from '../components/StudyAdmin/CreateStudyModal';
import { DeleteStudyModal } from '../components/StudyAdmin/DeleteStudyModal';
import { StudyLogModal } from '../components/StudyAdmin/StudyLogModal';
import { EditStudyModal } from '../components/StudyAdmin/EditStudyModal';
import { NavBarProvider } from '../components/NavBar/NavBar';

const useStyles = createStyles(() => ({
  fullW: {
    width: '100%',
  },
  cursor: {
    cursor: 'pointer',
  },
}));

function StudyImportLogCell({ row, setSelectedLogStudy }: { row: StudyAdminDetailsFragment; setSelectedLogStudy: (row: StudyAdminDetailsFragment) => void }) {
  const { classes } = useStyles();
  const selectStudyLog = useCallback(() => {
    setSelectedLogStudy(row);
  }, [row, setSelectedLogStudy]);

  return !row.hasImportLog ? null : <IconArticle onClick={selectStudyLog} className={classes.cursor} />;
}

function StudyImportAdminCell({
  row,
  setSelectedEditStudy,
}: {
  row: StudyAdminDetailsFragment;
  setSelectedEditStudy: (row: StudyAdminDetailsFragment) => void;
}) {
  const { classes } = useStyles();
  const selectStudyEdit = useCallback(() => {
    setSelectedEditStudy(row);
  }, [row, setSelectedEditStudy]);

  return !row.adminPermissionGranted ? null : <IconDotsVertical onClick={selectStudyEdit} className={classes.cursor} />;
}

function StudyImportDeleteCell({
  row,
  setSelectedDeleteStudy,
}: {
  row: StudyAdminDetailsFragment;
  setSelectedDeleteStudy: (row: StudyAdminDetailsFragment) => void;
}) {
  const { classes } = useStyles();
  const selectStudyDelete = useCallback(() => {
    setSelectedDeleteStudy(row);
  }, [row, setSelectedDeleteStudy]);

  return !row.adminPermissionGranted ? null : <IconX onClick={selectStudyDelete} className={classes.cursor} />;
}

export default function StudyAdmin() {
  const { classes } = useStyles();
  const [newStudyModalOpen, setNewStudyModalOpen] = useState(false);

  const { data, loading, refetch, error } = useStudyAdminListQuery();
  const [selectedEditStudy, setSelectedEditStudy] = useState<StudyAdminDetailsFragment | undefined>(undefined);
  const [selectedDeleteStudy, setSelectedDeleteStudy] = useState<StudyAdminDetailsFragment | undefined>(undefined);
  const [selectedLogStudy, setSelectedLogStudy] = useState<StudyAdminDetailsFragment | undefined>(undefined);

  const resetNewStudyModal = useCallback(async () => {
    setNewStudyModalOpen(false);
    await refetch();
  }, [setNewStudyModalOpen, refetch]);

  const resetDeleteModal = useCallback(async () => {
    setSelectedDeleteStudy(undefined);
    await refetch();
  }, [setSelectedDeleteStudy, refetch]);

  const resetEditModal = useCallback(async () => {
    setSelectedEditStudy(undefined);
    await refetch();
  }, [setSelectedEditStudy, refetch]);

  const StudyImportLogCellFn = useCallback((row: StudyAdminDetailsFragment) => {
    return <StudyImportLogCell row={row} setSelectedLogStudy={setSelectedLogStudy} />;
  }, []);

  const StudyImportAdminCellFn = useCallback((row: StudyAdminDetailsFragment) => {
    return <StudyImportAdminCell row={row} setSelectedEditStudy={setSelectedEditStudy} />;
  }, []);

  const StudyImportDeleteCellFn = useCallback((row: StudyAdminDetailsFragment) => {
    return <StudyImportDeleteCell row={row} setSelectedDeleteStudy={setSelectedDeleteStudy} />;
  }, []);

  const columns = useMemo(() => {
    return [
      {
        name: 'ID',
        selector: (row: StudyAdminDetailsFragment) => row.studyId,
        sortable: true,
        width: '15%',
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
              : 'Importing...'
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
        cell: StudyImportLogCellFn,
      },
      {
        name: '',
        width: '50px',
        cell: StudyImportAdminCellFn,
      },
      {
        name: '',
        width: '50px',
        cell: StudyImportDeleteCellFn,
      },
    ];
  }, [StudyImportAdminCellFn, StudyImportDeleteCellFn, StudyImportLogCellFn]);

  const openNewStudyModal = useCallback(() => {
    setNewStudyModalOpen(!newStudyModalOpen);
  }, [newStudyModalOpen]);

  const resetStudyLogModal = useCallback(() => {
    setSelectedLogStudy(undefined);
  }, []);

  return (
    <NavBarProvider>
      <EditStudyModal opened={selectedEditStudy !== undefined} reset={resetEditModal} study={selectedEditStudy} />
      <CreateStudyModal opened={newStudyModalOpen} reset={resetNewStudyModal} />
      <DeleteStudyModal study={selectedDeleteStudy} reset={resetDeleteModal} opened={selectedDeleteStudy !== undefined} />
      <StudyLogModal opened={selectedLogStudy !== undefined} study={selectedLogStudy} reset={resetStudyLogModal} />
      <Space h="xl" />
      <Stack px="md">
        <Group position="right">
          {(data?.userStudyUploadConfigured || true) && (
            <Button onClick={openNewStudyModal}>
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
          <DataTable data={data?.studyAdminDetailsList || []} columns={columns} defaultSortFieldId={1} defaultSortAsc={false} className={classes.fullW} />
        ) : null}
        {error && (
          <Group position="center">
            <Text weight="bold">An Error occurred</Text>
          </Group>
        )}
      </Stack>
    </NavBarProvider>
  );
}
