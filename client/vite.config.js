import {defineConfig} from "vite";
import react from "@vitejs/plugin-react";
import svgr from "vite-plugin-svgr";
const fetch = (...args) => import('node-fetch').then(({default: fetch}) => fetch(...args));

export default defineConfig(() => ({
    server: {
        open: true,
        proxy: {
            "^/postgraphile/.*": {
                target: "http://localhost:4000"
            },
        },
    },
    build: {
        outDir: "build",
        sourcemap: true,
    },
    plugins: [react(), svgr({svgrOptions: {icon: true}})],
}));
