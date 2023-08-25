import { useMemo, useState } from 'react';
import { Center, Group, Text } from '@mantine/core';
import { LeftSidePanel } from '../components';
import { SankeyPlot } from '../components/';
import { SankeyAnnotationGroupSelector } from '../components/SankeyAnnotationGroupSelector/SankeyAnnotationGroupSelector';
import { useRecoilValue } from 'recoil';
import { studyState } from '../atoms';
import { SelectBoxItem } from '../model';
import { showNotification } from '@mantine/notifications';

const AnnotationComparison = () => {
  const study = useRecoilValue(studyState);

  const [value1, setValue1] = useState<string | undefined>();
  const [value2, setValue2] = useState<string | undefined>();

  function handleChange1(value: string) {
    if (value !== value2) setValue1(value);
    else {
      showNotification({
        title: 'Please choose two different annotation groups. ',
        message: 'Nice try!',
        color: 'red',
        autoClose: 2500,
      });
    }
  }

  const annotations: SelectBoxItem[] = useMemo(() => {
    const anns: SelectBoxItem[] = [];
    if (study) {
      study.annotationGroupMap.forEach((value, key) => {
        anns.push({
          value: key.toString(),
          label: value.displayGroup,
        });
      });
    }
    if (anns.length > 1) {
      setValue1(anns[0].value);
      setValue2(anns[1].value);
    }
    return anns;
  }, [study]);
  return (
    <Group position={'apart'} noWrap>
      <LeftSidePanel>
        <SankeyAnnotationGroupSelector annotationGroups={annotations} handleChange1={handleChange1} value1={value1} handleChange2={setValue2} value2={value2} />
      </LeftSidePanel>
      <main
        style={{
          height: '100vh',
          overflowY: 'scroll',
          flexGrow: 1,
          paddingTop: 60,
        }}
      >
        {study && annotations.length >= 2 && value1 && value2 && (
          <SankeyPlot
            annotationValues1={study.annotationGroupMap.get(parseInt(value1))?.annotationValuesList}
            annotationValues2={study.annotationGroupMap.get(parseInt(value2))?.annotationValuesList}
            annotationGroupId1={parseInt(value1)}
            annotationGroupId2={parseInt(value2)}
            studyId={study.studyId}
          />
        )}
        {annotations.length < 2 && (
          <Center style={{ height: '100%', width: '100%' }}>
            <Text>My friend. Please think twice: Does a comparison of annotation groups make sense when there is only 1???</Text>
          </Center>
        )}
      </main>
    </Group>
  );
};

export default AnnotationComparison;
