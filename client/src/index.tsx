import React from 'react';
import ReactDOM from 'react-dom/client';
import './index.css';
import reportWebVitals from './reportWebVitals';
import {ApolloClient, ApolloProvider, HttpLink, InMemoryCache, NormalizedCacheObject} from '@apollo/client';
import {MantineProvider} from '@mantine/core';
import './fonts/Exo-Bold.ttf';
import {createBrowserRouter, RouterProvider,} from "react-router-dom";
import SearchResults from "./pages/SearchResults";
import CellMarkerAnalysis from "./pages/CellMarkerAnalysis";
import AnnotationsInUmapScatterplotTestPage from "./pages/AnnotationsInUmapScatterplotTestPage";
import {GlobalLoading} from "./pages/GlobalLoading";
import {RecoilRoot} from "recoil";
import SingleGeneExpressionInUmapScatterplotTestPage from "./pages/SingleGeneExpressionInUmapScatterplotTestPage";
import ExpressionAnalysis from "./pages/ExpressionAnalysis";

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
    },
    {
        path: "/cellmarkeranalysis",
        element: <CellMarkerAnalysis/>,
    },
    {
        path: "/expressionanalysis",
        element: <ExpressionAnalysis/>,
    },
    {
        path: "/AnnotationsInUmapScatterplotTestPage",
        element: <AnnotationsInUmapScatterplotTestPage/>,
    },
    {
        path: "/SingleGeneExpressionInUmapScatterplotTestPage",
        element: <SingleGeneExpressionInUmapScatterplotTestPage/>,
    },
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
                    <RouterProvider router={router}/>
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
