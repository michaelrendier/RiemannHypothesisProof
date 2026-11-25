#!/bin/bash
# Open the JavaScript visualizations through a local web server on 127.0.0.1 (loopback only). Usage: ./serve_visualizations.sh [port=8770] [page=webgl_aizawa_vs_zeta_uns.html]
cd "$(dirname "$0")/javascript"; PORT=${1:-8770}; PAGE=${2:-index.html}
curl -s -o /dev/null "http://127.0.0.1:$PORT/favicon.png" || { nohup python3 -m http.server "$PORT" --bind 127.0.0.1 >/dev/null 2>&1 & sleep 1; }
xdg-open "http://127.0.0.1:$PORT/$PAGE" >/dev/null 2>&1
