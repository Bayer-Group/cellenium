import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import svgr from 'vite-plugin-svgr';
import Markdown from 'vite-plugin-react-markdown';

export default defineConfig(() => ({
  server: {
    open: true,
    proxy: {
      '^/postgraphile/.*': {
        target: 'http://localhost:5000',
      },
    },
  },
  build: {
    outDir: 'build',
    sourcemap: true,
  },
  plugins: [
    Markdown(),
    react({
      include: [/\.tsx$/, /\.md$/], // <-- add .md
    }),
    svgr({ svgrOptions: { icon: true } }),
  ],
}));
