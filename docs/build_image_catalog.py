#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
"""Build the image-catalog PDF from docs/IMAGE_CATALOG.md: renders each figure (SVG → PNG), inserts it under its heading, pandoc → HTML → headless Chrome → PDF.
Run from the repo root:  python3 docs/build_image_catalog.py [OUT.pdf]   (needs pandoc, rsvg-convert, google-chrome)"""
import pathlib, re, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "docs" / "RHP_Image_Catalog.pdf"
CSS = """
@page { size: A4; margin: 14mm 13mm; }
body { font-family: 'DejaVu Sans', sans-serif; font-size: 9.2pt; line-height: 1.38; color: #1b1b1f; max-width: none; }
h1 { font-size: 17pt; margin: 0 0 4px; } h2 { font-size: 12.5pt; margin: 16px 0 4px; border-bottom: 1px solid #999; padding-bottom: 2px; page-break-after: avoid; }
code { font-family: 'DejaVu Sans Mono', monospace; font-size: 8.4pt; background: #f1f1f4; padding: 0 2px; }
img { max-width: 100%; display: block; margin: 5px 0 7px; border: 1px solid #ddd; }
li { margin-bottom: 2px; } p, li { orphans: 3; widows: 3; } hr { border: 0; border-top: 1px solid #bbb; margin: 12px 0; }
"""


def main():
    tmp = pathlib.Path(tempfile.mkdtemp(prefix="rhp_catalog_"))
    out_lines = []
    for line in (ROOT / "docs" / "IMAGE_CATALOG.md").read_text().splitlines():
        if line.startswith("- ") and out_lines and out_lines[-1].strip() and not out_lines[-1].startswith("- "):
            out_lines.append("")                      # pandoc needs a blank line before a list
        out_lines.append(line)
        m = re.match(r"## \d+\. `([^`]+)`", line)
        if m:
            src = ROOT / m.group(1)
            if src.suffix == ".svg" and src.exists():
                png = tmp / (src.stem + ".png")
                subprocess.run(["rsvg-convert", "-w", "1500", str(src), "-o", str(png)], check=True)
                out_lines += ["", f"![]({png})", ""]
            elif src.suffix == ".png" and src.exists():
                out_lines += ["", f"![]({src})", ""]
    md = tmp / "catalog.md"; md.write_text("\n".join(out_lines))
    css = tmp / "catalog.css"; css.write_text(CSS)
    html = tmp / "catalog.html"
    subprocess.run(["pandoc", "-f", "markdown-yaml_metadata_block", "-t", "html5", "-s", "--embed-resources", "--css", str(css), "--metadata", "pagetitle=RiemannHypothesisProof image catalog", "-o", str(html), str(md)], check=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["google-chrome", "--headless", "--no-sandbox", "--disable-gpu", f"--print-to-pdf={OUT}", "--no-pdf-header-footer", str(html)], check=True, stderr=subprocess.DEVNULL)
    print(OUT)


if __name__ == "__main__":
    main()
