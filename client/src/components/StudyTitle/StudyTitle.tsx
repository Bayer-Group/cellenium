import { Modal, Text } from '@mantine/core';
import { StudyInfoFragment, useSingleStudyInfoQuery } from '../../generated/types';
import { useRecoilValue } from 'recoil';
import { studyState } from '../../atoms.ts';
import { useState } from 'react';
import { StudyCard } from '../StudyCard/StudyCard';

function StudyInfoModal({ opened, onClose, study }: { opened: boolean; onClose: () => void; study: StudyInfoFragment }) {
  return (
    <Modal opened={opened} onClose={onClose} title={'Study information'} size="xl">
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

  return (
    <>
      <Text size="s" title={'click for complete study information'} onClick={() => setModalOpen(true)}>
        {study && study.studyName}
      </Text>
      {singleStudyInfoList && <StudyInfoModal opened={modalOpen} onClose={() => setModalOpen(false)} study={singleStudyInfoList?.studyOverviewsList[0]} />}
    </>
  );
}
