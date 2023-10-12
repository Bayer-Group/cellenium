import { Modal, Stack, Text } from '@mantine/core';
import { useRecoilValue } from 'recoil';
import { useCallback, useState } from 'react';
import { StudyCard } from '../StudyCard/StudyCard';
import { StudyInfoFragment, useSingleStudyInfoQuery } from '../../generated/types';
import { studyState } from '../../atoms';

function StudyInfoModal({ opened, onClose, study }: { opened: boolean; onClose: () => void; study: StudyInfoFragment }) {
  return (
    <Modal
      opened={opened}
      onClose={onClose}
      title={
        <Text weight="bold" size="md">
          Study information
        </Text>
      }
      size="xl"
    >
      <StudyCard study={study} detailed />
    </Modal>
  );
}

export function StudyTitle() {
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
      <Text weight="bold" truncate="end" size="md" title="click for complete study information" onClick={modalClick} style={{ cursor: 'pointer' }}>
        {study && study.studyName}
      </Text>
      {singleStudyInfoList && <StudyInfoModal opened={modalOpen} onClose={modalClick} study={singleStudyInfoList?.studyOverviewsList[0]} />}
    </Stack>
  );
}
