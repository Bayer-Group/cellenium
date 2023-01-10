import React from 'react';
import ReactDOM from 'react-dom/client';
import './index.css';
import reportWebVitals from './reportWebVitals';
import {ApolloClient, ApolloProvider, HttpLink, InMemoryCache, NormalizedCacheObject} from '@apollo/client';
import {MantineProvider} from '@mantine/core';
import './fonts/Exo-Bold.ttf';
import {createBrowserRouter, RouterProvider,} from "react-router-dom";
import SearchResults from "./pages/SearchResults";
import {GlobalLoading} from "./pages/GlobalLoading";
import {RecoilRoot} from "recoil";
import ErrorPage from "./pages/ErrorPage";
import {NotificationsProvider} from '@mantine/notifications';
import {ModalsProvider} from "@mantine/modals";
import OntologySandbox from "./pages/OntologySandbox";
import {StudyPage} from "./pages/StudyPage";

export const apolloCache = new InMemoryCache();

const link = new HttpLink({
    uri: '/postgraphile/'
});

const apolloClient: ApolloClient<NormalizedCacheObject> = new ApolloClient({
    cache: apolloCache,
    link
});

const router = createBrowserRouter([
    {
        path: "/",
        element: <SearchResults/>,
        errorElement: <ErrorPage/>,
    },
    {
        path: "/study/:studyId",
        element: <StudyPage/>,
    },
    {
        path: "ontologySandbox",
        element: <OntologySandbox/>
    }
]);
const root = ReactDOM.createRoot(
    document.getElementById('root') as HTMLElement
);

root.render(
    <ApolloProvider client={apolloClient}>
        <MantineProvider theme={{
            fontFamily: 'Rubik, sans-serif',
            headings: {
                fontFamily: 'Exo-bold, sans-serif'
            }
        }}>
            <RecoilRoot>
                <React.Suspense fallback={<GlobalLoading/>}>
                    <NotificationsProvider zIndex={10000000}>
                        <ModalsProvider>
                            <RouterProvider router={router}/>
                        </ModalsProvider>
                    </NotificationsProvider>
                </React.Suspense>
            </RecoilRoot>
        </MantineProvider>
    </ApolloProvider>
)
;

// If you want to start measuring performance in your app, pass a function
// to log results (for example: reportWebVitals(console.log))
// or send to an analytics endpoint. Learn more: https://bit.ly/CRA-vitals
reportWebVitals();

export {apolloClient};
