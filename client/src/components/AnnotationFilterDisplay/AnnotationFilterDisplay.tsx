import { ActionIcon, Group, Stack, Text } from '@mantine/core';
import AnnotationFilterSelectBox from './AnnotationFilterSelectBox';
import { IconTrash } from '@tabler/icons-react';
import { useRecoilState } from 'recoil';
import { selectedAnnotationFilterState } from '../../atoms';

const AnnotationFilterDisplay = () => {
  const [, setSelectedAnnotationFilter] = useRecoilState(selectedAnnotationFilterState);
  return (
    <Stack spacing={0}>
      <Group align={'center'} position={'apart'} noWrap={true}>
        <Text size={'xs'}>Filter out cells with annotation(s)</Text>
        <ActionIcon
          size={'xs'}
          onClick={() => {
            setSelectedAnnotationFilter([]);
          }}
        >
          <IconTrash />
        </ActionIcon>
      </Group>
      <Group>
        <AnnotationFilterSelectBox />
      </Group>
    </Stack>
  );
};

export { AnnotationFilterDisplay };
