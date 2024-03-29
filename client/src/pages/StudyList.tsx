import { useState } from 'react';
import { Stack } from '@mantine/core';
import { StudyInfoFragment } from '../generated/types';
import { StudySearchBar } from '../components/SearchBar/StudySearchBar';
import { StudyCard } from '../components/StudyCard/StudyCard';
import { NavBarProvider } from '../components/NavBar/NavBar';

function StudyList() {
  const [studyList, setStudyList] = useState<StudyInfoFragment[]>([]);

  return (
    <NavBarProvider scrollable>
      <Stack p="md" spacing={0}>
        <StudySearchBar onStudyListUpdate={setStudyList} />
        <Stack w="100%" pb="md" mt="1rem">
          {studyList.map((study) => (
            <StudyCard study={study} key={study.studyId} />
          ))}
        </Stack>
      </Stack>
    </NavBarProvider>
  );
}

export default StudyList;
