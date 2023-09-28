import { useState } from 'react';
import { Container, Grid, Space } from '@mantine/core';
import { StudyInfoFragment } from '../generated/types';
import { StudySearchBar } from '../components/SearchBar/StudySearchBar';
import { NavBar } from '../components/NavBar/NavBar';
import { StudyCard } from '../components/StudyCard/StudyCard';

function StudyList() {
  const [studyList, setStudyList] = useState<StudyInfoFragment[]>([]);

  return (
    <Container fluid>
      <NavBar />
      <Space h="xl" />
      <Container size="xl" style={{ paddingBottom: '2rem' }}>
        <StudySearchBar onStudyListUpdate={setStudyList} />
      </Container>
      <Container size="xl">
        <Grid>
          {studyList.map((study) => (
            <Grid.Col span={12} key={study.studyId}>
              <StudyCard study={study} />
            </Grid.Col>
          ))}
        </Grid>
      </Container>
    </Container>
  );
}

export default StudyList;
