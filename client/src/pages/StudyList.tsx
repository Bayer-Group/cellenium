import { useState } from 'react';
import { Container, Space, Stack } from '@mantine/core';
import { StudyInfoFragment } from '../generated/types';
import { StudySearchBar } from '../components/SearchBar/StudySearchBar';
import { StudyCard } from '../components/StudyCard/StudyCard';
import { NavBarProvider } from '../components/NavBar/NavBar';

function StudyList() {
  const [studyList, setStudyList] = useState<StudyInfoFragment[]>([]);

  return (
    <NavBarProvider scrollable>
      <Space h="xs" />
      <Container w="100%" size="xl" p={0}>
        <StudySearchBar onStudyListUpdate={setStudyList} />
      </Container>
      <Container w="100%" size="xl">
        <Stack w="100%" pb="md">
          {studyList.map((study) => (
            <StudyCard study={study} key={study.studyId} />
          ))}
        </Stack>
      </Container>
    </NavBarProvider>
  );
}

export default StudyList;
