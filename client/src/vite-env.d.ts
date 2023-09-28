/// <reference types="vite/client" />

declare module '*.md' {
  import React from 'react';

  const ReactComponent: React.VFC;
  export default ReactComponent;
  export const attributes = Record<string, any>;
}
