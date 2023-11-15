import { createStyles, Modal, Stack, Text } from '@mantine/core';
import { useRecoilValue } from 'recoil';
import { useCallback, useState } from 'react';
import { StudyCard } from '../StudyCard/StudyCard';
import { StudyInfoFragment, useSingleStudyInfoQuery } from '../../generated/types';
import { studyState } from '../../atoms';

export function StudyInfoModal({ opened, onClose, study }: { opened: boolean; onClose: () => void; study?: StudyInfoFragment }) {
  return (
    <Modal opened={opened} onClose={onClose} size="xl">
      {study && <StudyCard study={study} detailed />}
    </Modal>
  );
}

const useStyles = createStyles(() => ({
  cursor: {
    cursor: 'pointer',
  },
}));

export function StudyTitle() {
  const { classes } = useStyles();
  const study = useRecoilValue(studyState);
  const [modalOpen, setModalOpen] = useState(false);
  const { data: singleStudyInfoList } = useSingleStudyInfoQuery({
    variables: {
      studyId: study?.studyId || -1,
    },
    skip: !study,
  });

  const modalClick = useCallback(() => {
    setModalOpen(!modalOpen);
  }, [modalOpen]);

  return (
    <Stack>
      <Text weight="bold" truncate="end" size="md" title="click for complete study information" onClick={modalClick} className={classes.cursor}>
        {study && study.studyName}
      </Text>
      {singleStudyInfoList && <StudyInfoModal opened={modalOpen} onClose={modalClick} study={singleStudyInfoList?.studyOverviewsList[0]} />}
    </Stack>
  );
}
