import { Loader, Stack } from '@mantine/core';

function GlobalLoading() {
  return (
    <Stack style={{ width: '100%', height: '100%' }} align="center" justify="center">
      <Loader size="xl" color="blue" variant="dots" />
    </Stack>
  );
}

export { GlobalLoading };
