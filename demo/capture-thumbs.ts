// Thumbnail Capture — How It Works
//
// Thumbnails are captured from within the running demo app itself.
// On localhost, an "Update Thumbnails" button appears in the sidebar footer.
//
// How to update thumbnails:
// 1. Start the dev server: bun run dev
// 2. Open http://localhost:3000
// 3. Click the "Update Thumbnails" button in the sidebar
// 4. The app will navigate to each demo, wait for render, capture the canvas,
//    and POST the JPEG to the dev server which saves it to public/thumbs/
//
// The button is only visible on localhost — it won't appear on GitHub Pages.
//
// Technical details:
// - Canvas has preserveDrawingBuffer: true (set in three-viewer.ts)
// - Each capture waits 2 seconds for WASM + Three.js to fully render
// - Images are saved as 320x240 JPEG at 85% quality
// - The dev server's POST /api/save-thumb endpoint handles file saving
//
// To add a new demo's thumbnail:
// 1. Add the route name to DEMO_ROUTES array in src/main.ts
// 2. Click "Update Thumbnails" to regenerate all
