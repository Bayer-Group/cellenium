import { useRouteError } from 'react-router-dom';
import { Stack, Title, Text } from '@mantine/core';
import { NavBarProvider } from '../components/NavBar/NavBar';

export default function ErrorPage() {
  const error: { statusText: string; message: string } = useRouteError() as { statusText: string; message: string };
  console.error(error);

  return (
    <NavBarProvider>
      <Stack w="100vw" h="100vh" align="center" pt="md">
        <Title>Oops!</Title>
        <Text>Sorry, an unexpected error has occurred.</Text>
        <Text>
          <i>{error.statusText || error.message}</i>
        </Text>
      </Stack>
    </NavBarProvider>
  );
}
