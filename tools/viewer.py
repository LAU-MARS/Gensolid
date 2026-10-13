#!/usr/bin/env python3
"""Local-only geometry viewer: static UI + bounded Rust scene subprocesses."""
import argparse
import json
import mimetypes
import pathlib
import subprocess
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

ROOT = pathlib.Path(__file__).resolve().parents[1]
GATE = threading.BoundedSemaphore(2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--no-build", action="store_true")
    args = parser.parse_args()
    if not args.no_build:
        subprocess.run(["cargo", "build", "-p", "gensolid-viewer"], cwd=ROOT, check=True)
    binary = ROOT / "target" / "debug" / ("gensolid-viewer.exe" if __import__("os").name == "nt" else "gensolid-viewer")
    if not binary.is_file():
        parser.error(f"build the viewer first: {binary}")

    class Handler(BaseHTTPRequestHandler):
        def reply(self, status, data, mime="application/json; charset=utf-8"):
            self.send_response(status)
            self.send_header("Content-Type", mime)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Content-Security-Policy", "default-src 'self'; style-src 'self'; script-src 'self'; connect-src 'self'; img-src 'self' data:; frame-ancestors 'none'")
            self.end_headers()
            self.wfile.write(data)

        def error(self, status, message):
            self.reply(status, json.dumps({"status": "error", "message": message}).encode())

        def do_GET(self):
            name = {"/": "index.html", "/app.js": "app.js", "/style.css": "style.css"}.get(self.path)
            if name is None:
                self.error(404, "not found")
                return
            self.reply(200, (ROOT / "viewer" / name).read_bytes(), (mimetypes.guess_type(name)[0] or "text/plain") + "; charset=utf-8")

        def do_POST(self):
            if self.path != "/api/scene":
                self.error(404, "not found")
                return
            # Require JSON and same-origin browser requests; no permissive CORS.
            if self.headers.get("Content-Type", "").split(";")[0] != "application/json":
                self.error(415, "JSON required")
                return
            origin = self.headers.get("Origin")
            if origin and origin not in (f"http://127.0.0.1:{args.port}", f"http://localhost:{args.port}"):
                self.error(403, "same-origin requests only")
                return
            try:
                size = int(self.headers.get("Content-Length", "0"))
            except ValueError:
                size = 0
            if not 0 < size <= 65536:
                self.error(413, "request must contain 1..65536 bytes")
                return
            data = self.rfile.read(size)
            if not GATE.acquire(blocking=False):
                self.error(429, "two computations already running; retry shortly")
                return
            try:
                run = subprocess.run([str(binary)], input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ROOT, timeout=90)
                result = json.loads(run.stdout)
                self.reply(200 if run.returncode == 0 else 422, json.dumps(result, ensure_ascii=False).encode("utf-8"))
            except subprocess.TimeoutExpired:
                self.error(504, "geometry computation exceeded 90 seconds")
            except (OSError, ValueError) as error:
                self.error(500, str(error))
            finally:
                GATE.release()

    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"Gensolid viewer: http://127.0.0.1:{args.port}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
