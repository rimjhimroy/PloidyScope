#!/usr/bin/env python3
"""Export SVG logo assets to PNG.

Usage:
    python scripts/export_logo_png.py

The script looks for one of these rasterization backends:
1. rsvg-convert
2. inkscape
3. magick
4. convert
5. cairosvg in the active Python environment

If none are available, it exits with a clear error message describing the
supported options.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ASSETS_DIR = PROJECT_ROOT / "assets"
LOGO_SVG = ASSETS_DIR / "ploidyscope-logo.svg"
ICON_SVG = ASSETS_DIR / "ploidyscope-icon.svg"
LOGO_PNG = ASSETS_DIR / "ploidyscope-logo.png"
ICON_PNG = ASSETS_DIR / "ploidyscope-icon.png"


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def export_with_rsvg() -> bool:
    executable = shutil.which("rsvg-convert")
    if executable is None:
        return False
    run([executable, "-w", "1200", "-h", "360", "-o", str(LOGO_PNG), str(LOGO_SVG)])
    run([executable, "-w", "256", "-h", "256", "-o", str(ICON_PNG), str(ICON_SVG)])
    return True


def export_with_inkscape() -> bool:
    executable = shutil.which("inkscape")
    if executable is None:
        return False
    run([executable, str(LOGO_SVG), "--export-filename", str(LOGO_PNG), "--export-width", "1200", "--export-height", "360"])
    run([executable, str(ICON_SVG), "--export-filename", str(ICON_PNG), "--export-width", "256", "--export-height", "256"])
    return True


def export_with_imagemagick(command_name: str) -> bool:
    executable = shutil.which(command_name)
    if executable is None:
        return False
    run([executable, str(LOGO_SVG), "-resize", "1200x360", str(LOGO_PNG)])
    run([executable, str(ICON_SVG), "-resize", "256x256", str(ICON_PNG)])
    return True


def export_with_cairosvg() -> bool:
    try:
        import cairosvg
    except ModuleNotFoundError:
        return False

    cairosvg.svg2png(url=str(LOGO_SVG), write_to=str(LOGO_PNG), output_width=1200, output_height=360)
    cairosvg.svg2png(url=str(ICON_SVG), write_to=str(ICON_PNG), output_width=256, output_height=256)
    return True


def main() -> int:
    for exporter in (
        export_with_rsvg,
        export_with_inkscape,
        lambda: export_with_imagemagick("magick"),
        lambda: export_with_imagemagick("convert"),
        export_with_cairosvg,
    ):
        if exporter():
            print(f"Wrote {LOGO_PNG}")
            print(f"Wrote {ICON_PNG}")
            return 0

    print(
        "No SVG rasterizer found. Install one of: rsvg-convert, inkscape, magick, convert, or Python package cairosvg.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())