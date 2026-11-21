#!/bin/bash
# Open the Riemann–Fermat Atlas Gallery through a local web server on 127.0.0.1 (Firefox shows the favicon only for http pages). Serves the repository root so the image links resolve.
# Usage: ./serve_gallery.sh [port=8769]   (loopback only; stop with: pkill -f "http.server 8769")
cd "$(dirname "$0")/.."; PORT=${1:-8769}
curl -s -o /dev/null "http://127.0.0.1:$PORT/image_gallery/favicon.png" || { nohup python3 -m http.server "$PORT" --bind 127.0.0.1 >/dev/null 2>&1 & sleep 1; }
xdg-open "http://127.0.0.1:$PORT/image_gallery/Riemann_Fermat_Atlas_Gallery.html" >/dev/null 2>&1
