import { ApolloClient, HttpLink, InMemoryCache, NormalizedCacheObject } from '@apollo/client';

export const apolloCache = new InMemoryCache();

const link = new HttpLink({
  uri: '/postgraphile/',
});

export const apolloClient: ApolloClient<NormalizedCacheObject> = new ApolloClient({
  cache: apolloCache,
  link,
});
