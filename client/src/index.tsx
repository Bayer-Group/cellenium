import React, { lazy } from 'react';
import ReactDOM from 'react-dom/client';
import './index.css';
import reportWebVitals from './reportWebVitals';
import { ApolloClient, ApolloProvider, HttpLink, InMemoryCache, NormalizedCacheObject } from '@apollo/client';
import { Container, MantineProvider } from '@mantine/core';
import './fonts/Exo-Bold.ttf';
import { createBrowserRouter, RouterProvider } from 'react-router-dom';
import { GlobalLoading } from './pages/GlobalLoading';
import { RecoilRoot } from 'recoil';
import ErrorPage from './pages/ErrorPage';
import { Notifications } from '@mantine/notifications';
import { ModalsProvider } from '@mantine/modals';
import { Docs } from './pages/Docs.tsx';

export const apolloCache = new InMemoryCache();

const link = new HttpLink({
  uri: '/postgraphile/',
});

const apolloClient: ApolloClient<NormalizedCacheObject> = new ApolloClient({
  cache: apolloCache,
  link,
});

const StudyList = lazy(() => import('./pages/StudyList'));
const MarkerGeneSearch = lazy(() => import('./pages/MarkerGeneSearch'));
const CrossStudySearch = lazy(() => import('./pages/CrossStudySearch'));
const StudyPage = lazy(() => import('./pages/StudyPage'));
const StudyAdmin = lazy(() => import('./pages/StudyAdmin'));
const OntologySandbox = lazy(() => import('./pages/OntologySandbox'));

const router = createBrowserRouter([
  {
    path: '/',
    element: <StudyList />,
    errorElement: <ErrorPage />,
  },
  {
    path: '/markergene',
    element: <MarkerGeneSearch />,
    errorElement: <ErrorPage />,
  },
  {
    path: '/crossstudy',
    element: <CrossStudySearch />,
    errorElement: <ErrorPage />,
  },
  {
    path: '/study/:studyId',
    element: <StudyPage />,
  },
  {
    path: 'study-admin',
    element: <StudyAdmin />,
  },
  {
    path: 'ontology-sandbox',
    element: <OntologySandbox />,
  },
  {
    path: 'docs',
    element: <Docs />,
  },
]);

const root = ReactDOM.createRoot(document.getElementById('root') as HTMLElement);

root.render(
  <ApolloProvider client={apolloClient}>
    <MantineProvider
      theme={{
        fontFamily: 'Rubik, sans-serif',
        headings: {
          fontFamily: 'Exo-bold, sans-serif',
        },
      }}
    >
      <RecoilRoot>
        <React.Suspense fallback={<GlobalLoading />}>
          <Notifications />
          <ModalsProvider>
            <Container style={{ padding: 0 }} fluid={true}>
              <RouterProvider router={router} />
            </Container>
          </ModalsProvider>
        </React.Suspense>
      </RecoilRoot>
    </MantineProvider>
  </ApolloProvider>,
);

// If you want to start measuring performance in your app, pass a function
// to log results (for example: reportWebVitals(console.log))
// or send to an analytics endpoint. Learn more: https://bit.ly/CRA-vitals
reportWebVitals();

export { apolloClient };
