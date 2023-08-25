import { Group } from '@mantine/core';
import { LeftSidePanel, RightSidePanel } from '../components';

const UserAnnotation = () => {
  return (
    <Group position={'apart'}>
      <LeftSidePanel></LeftSidePanel>

      <RightSidePanel></RightSidePanel>
    </Group>
  );
};

export default UserAnnotation;
