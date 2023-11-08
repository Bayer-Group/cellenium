import React, { lazy } from 'react';
import ReactDOM from 'react-dom/client';
import './index.css';
import { MantineProvider, Stack } from '@mantine/core';
import './assets/fonts/Exo-Bold.ttf';
import { createBrowserRouter, RouterProvider } from 'react-router-dom';
import { RecoilRoot } from 'recoil';
import { Notifications } from '@mantine/notifications';
import { ModalsProvider } from '@mantine/modals';
import { ApolloProvider } from '@apollo/client';
import ErrorPage from './pages/ErrorPage';
import { GlobalLoading } from './pages/GlobalLoading';
import reportWebVitals from './reportWebVitals';
import { Docs } from './pages/Docs';
import { apolloClient } from './client';
import { NavBarProvider } from './components/NavBar/NavBar';

const StudyList = lazy(() => import('./pages/StudyList'));
const MarkerGeneSearch = lazy(() => import('./pages/MarkerGeneSearch'));
const CrossStudySearch = lazy(() => import('./pages/CrossStudySearch'));
const StudyPage = lazy(() => import('./pages/StudyPage'));
const StudyAdmin = lazy(() => import('./pages/StudyAdmin'));

const router = createBrowserRouter([
  {
    path: '/',
    element: (
      <NavBarProvider scrollable>
        <StudyList />
      </NavBarProvider>
    ),
    errorElement: <ErrorPage />,
  },
  {
    path: '/markergene',
    element: (
      <NavBarProvider scrollable>
        <MarkerGeneSearch />
      </NavBarProvider>
    ),
    errorElement: <ErrorPage />,
  },
  {
    path: '/crossstudy',
    element: (
      <NavBarProvider scrollable>
        <CrossStudySearch />
      </NavBarProvider>
    ),
    errorElement: <ErrorPage />,
  },
  {
    path: '/study/:studyId',
    element: <StudyPage />,
  },
  {
    path: '/analysis/:studyId',
    element: <StudyPage />,
  },
  {
    path: 'study-admin',
    element: (
      <NavBarProvider>
        <StudyAdmin />
      </NavBarProvider>
    ),
  },
  {
    path: 'manageStudyImport',
    element: (
      <NavBarProvider>
        <StudyAdmin />
      </NavBarProvider>
    ),
  },
  {
    path: 'docs',
    element: (
      <NavBarProvider>
        <Docs />
      </NavBarProvider>
    ),
  },
]);

const root = ReactDOM.createRoot(document.getElementById('root') as HTMLElement);

root.render(
  <MantineProvider
    theme={{
      fontFamily: 'Rubik, sans-serif',
      headings: {
        fontFamily: 'Exo-bold, sans-serif',
      },
    }}
  >
    <Stack w="100vw" pos="relative" spacing={0} h="100vh" style={{ overflow: 'hidden' }}>
      <ApolloProvider client={apolloClient}>
        <RecoilRoot>
          <React.Suspense fallback={<GlobalLoading />}>
            <Notifications />
            <ModalsProvider>
              <RouterProvider router={router} />
            </ModalsProvider>
          </React.Suspense>
        </RecoilRoot>
      </ApolloProvider>
    </Stack>
  </MantineProvider>,
);

// If you want to start measuring performance in your app, pass a function
// to log results (for example: reportWebVitals(console.log))
// or send to an analytics endpoint. Learn more: https://bit.ly/CRA-vitals
reportWebVitals();

export { apolloClient };
