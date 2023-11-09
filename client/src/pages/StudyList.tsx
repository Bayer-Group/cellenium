import { Stack } from '@mantine/core';
import { useRecoilState } from 'recoil';
import { StudySearchBar } from '../components/SearchBar/StudySearchBar';
import { StudyCard } from '../components/StudyCard/StudyCard';
import { StudySearchList } from '../atoms';

function StudyList() {
  const [studyList] = useRecoilState(StudySearchList);

  return (
    <Stack p="md" spacing={0}>
      <StudySearchBar initialFocus />
      <Stack w="100%" pb="md" mt="1rem">
        {studyList.map((study) => (
          <StudyCard study={study} key={study.studyId} />
        ))}
      </Stack>
    </Stack>
  );
}

export default StudyList;
