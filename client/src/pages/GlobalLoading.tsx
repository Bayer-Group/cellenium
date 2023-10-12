import { Loader, Stack } from '@mantine/core';

export function GlobalLoading() {
  return (
    <Stack w="100%" h="100%" align="center" justify="center">
      <Loader size="xl" color="blue" variant="dots" />
    </Stack>
  );
}
