import { readFileSync, existsSync, statSync } from "fs";
import { join, extname } from "path";

const PORT = parseInt(process.env.PORT || "3000");
const ROOT = import.meta.dir;

const MIME_TYPES: Record<string, string> = {
  ".html": "text/html",
  ".css": "text/css",
  ".js": "application/javascript",
  ".ts": "application/javascript",
  ".wasm": "application/wasm",
  ".json": "application/json",
  ".svg": "image/svg+xml",
  ".png": "image/png",
  ".ico": "image/x-icon",
};

const transpiler = new Bun.Transpiler({ loader: "ts", target: "browser" });

function tryServe(pathname: string): { content: string | Uint8Array; mime: string } | null {
  // index.html for root
  if (pathname === "/" || pathname === "/index.html") {
    const p = join(ROOT, "index.html");
    if (existsSync(p)) return { content: readFileSync(p, "utf-8"), mime: "text/html" };
  }

  // TypeScript files from src/ — transpile on the fly
  if (pathname.startsWith("/src/") && pathname.endsWith(".ts")) {
    const p = join(ROOT, pathname);
    if (existsSync(p)) {
      const code = transpiler.transformSync(readFileSync(p, "utf-8"));
      return { content: code, mime: "application/javascript" };
    }
  }

  // Any file relative to demo root (handles /styles/*, /public/pkg/*, etc.)
  const filePath = join(ROOT, pathname);
  if (existsSync(filePath) && statSync(filePath).isFile()) {
    const ext = extname(filePath);
    const mime = MIME_TYPES[ext] || "application/octet-stream";
    if (ext === ".wasm") {
      return { content: new Uint8Array(readFileSync(filePath)), mime };
    }
    return { content: readFileSync(filePath, "utf-8"), mime };
  }

  return null;
}

const server = Bun.serve({
  port: PORT,
  fetch(req) {
    const url = new URL(req.url);
    const resolved = tryServe(url.pathname);
    if (resolved) {
      return new Response(resolved.content, {
        headers: {
          "Content-Type": resolved.mime,
          "Access-Control-Allow-Origin": "*",
          "Cache-Control": "no-cache",
        },
      });
    }

    // SPA fallback: serve index.html for unresolved routes (client-side routing)
    const indexPath = join(ROOT, "index.html");
    if (existsSync(indexPath)) {
      return new Response(readFileSync(indexPath, "utf-8"), {
        headers: { "Content-Type": "text/html", "Cache-Control": "no-cache" },
      });
    }

    return new Response("Not found", { status: 404 });
  },
});

console.log(`Manifold Demo running at http://localhost:${server.port}`);
