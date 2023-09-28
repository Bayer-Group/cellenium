import { Center, Loader } from '@mantine/core';

function GlobalLoading() {
  return (
    <Center style={{ width: '100vw', height: '100vh' }}>
      <Loader size="xl" color="gray" variant="dots" />
    </Center>
  );
}

export { GlobalLoading };
