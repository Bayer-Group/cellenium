import { Group } from '@mantine/core';
import { RightSidePanel } from '../components/RightSidePanel/RightSidePanel';
import { LeftSidePanel } from '../components/LeftSidePanel/LeftSidePanel';

function UserAnnotation() {
  return (
    <Group position="apart">
      <LeftSidePanel />
      <RightSidePanel />
    </Group>
  );
}

export default UserAnnotation;
