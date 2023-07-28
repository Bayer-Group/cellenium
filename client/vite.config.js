import {defineConfig} from "vite";
import react from "@vitejs/plugin-react";
import svgr from "vite-plugin-svgr";

export default defineConfig(() => ({
    server: {
        open: true,
        proxy: {
            "^/postgraphile/.*": {
                target: "http://localhost:6000"
            },
        },
    },
    build: {
        outDir: "build",
        sourcemap: true,
    },
    plugins: [react(), svgr({svgrOptions: {icon: true}})],
}));
