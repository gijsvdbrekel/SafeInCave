#!/usr/bin/env python3
"""
render_domain_overview.py — shows the full 3D computational domain (the
450x450x660 m salt block) with each cavern void embedded inside, to answer the
reviewer question of whether the full domain is meshed without symmetry cuts.

Style follows the thesis schematic (Figure 3.7): a plain wireframe box, with
the cavern rendered inside in the same smooth salt-shaded style as the
individual cavern renders. The box is sized to a fixed, shared bounding cube
(large enough to contain the tallest cavern plus a small margin) rather than
the true 450x450x660 m domain, so the caverns fill most of the box and remain
clearly legible; the true domain scale is not the point of this figure, only
that each geometry is embedded in an explicit 3D box with no symmetry cuts.

Renders all six cavern geometries and assembles them into the same 3x2
montage layout as cavern_comparison.png (three per row, two rows), with dark
print-style labels above each panel.

Off-screen (headless / PuTTY) safe. Reads only the grids, so it runs anywhere.

Usage:
    python render_domain_overview.py
"""

import os
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont

import render_cavern_shapes as base   # reuse the cavern-wall loader + style

pv.OFF_SCREEN = True

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GRIDS_ROOT = base.GRIDS_ROOT
OUT_DIR = os.path.join(SCRIPT_DIR, "cavern_renders")

# (display label, grid folder) -- same six shapes as render_cavern_shapes.py
CAVERNS = base.CAVERNS

# Shared local box: a box large enough to contain the tallest cavern (the
# tilted cavern, ~249 m) plus a margin, centered on each cavern's own center.
# Taller than wide (height margin > width margin) to echo the thesis-style
# box proportions (Figure 3.7) rather than a perfect cube. This is
# deliberately much smaller than the true 450x450x660 m domain, so the
# caverns fill most of the box and stay legible in a small panel.
_TALLEST_EXTENT = 249.3   # m; tallest single-dimension cavern extent (tilt)
BOX_WIDTH = _TALLEST_EXTENT * 1.15
BOX_HEIGHT = _TALLEST_EXTENT * 1.40

EDGE_COLOR = "#1a2530"
EDGE_TUBE_RADIUS_FRAC = 0.006   # fraction of the box diagonal

WINDOW = (1500, 1500)
SUPERSAMPLE = 2
CAMERA_AZIMUTH = 35.0
CAMERA_ELEVATION = 18.0
FRAME_MARGIN = 1.12   # extra headroom around the box when fitting the camera


def _load_font(size):
    for path in ["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                 "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"]:
        if os.path.isfile(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def render_one(msh_path, out_path):
    """Wireframe local box + cavern inside, for a single geometry. The box is
    centered on the cavern's own center and sized to the shared BOX_WIDTH x
    BOX_WIDTH x BOX_HEIGHT, so every panel uses an identical box but is
    individually centered on its own cavern."""
    cavern = base.prepare_surface(base.load_cavern_surface(msh_path, base.CAVERN_TAG))
    cx, cy, cz = cavern.center
    hw, hh = BOX_WIDTH / 2, BOX_HEIGHT / 2
    local_bounds = (cx - hw, cx + hw, cy - hw, cy + hw, cz - hh, cz + hh)

    box = pv.Box(bounds=local_bounds)
    box_edges = box.extract_all_edges()
    diag = np.sqrt(2 * BOX_WIDTH ** 2 + BOX_HEIGHT ** 2)
    box_tubes = box_edges.tube(radius=EDGE_TUBE_RADIUS_FRAC * diag)

    p = pv.Plotter(off_screen=True, window_size=WINDOW)
    p.set_background("white")
    p.add_mesh(box_tubes, color=EDGE_COLOR, smooth_shading=True, ambient=0.4, diffuse=0.6)
    base._add_cavern(p, cavern)

    p.enable_parallel_projection()
    p.camera_position = "yz"
    p.camera.azimuth = CAMERA_AZIMUTH
    p.camera.elevation = CAMERA_ELEVATION
    # Fit the camera directly to the (padded) box bounds so the box's true
    # anisotropic, taller-than-wide shape is respected on screen -- rather
    # than estimating the zoom from the 3D body diagonal, which conflates
    # the height and width and washes out any height increase.
    hw_pad, hh_pad = hw * FRAME_MARGIN, hh * FRAME_MARGIN
    padded_bounds = (cx - hw_pad, cx + hw_pad, cy - hw_pad, cy + hw_pad,
                     cz - hh_pad, cz + hh_pad)
    p.reset_camera(bounds=padded_bounds)
    # Note: combining enable_anti_aliasing("ssaa") with screenshot(scale=2)
    # triggers a VTK aspect-ratio distortion bug when the camera is fit via
    # reset_camera() on an anisotropic box. The scale=2 supersampling below
    # already gives clean high-resolution output, so ssaa is omitted here.

    p.screenshot(out_path, transparent_background=True, scale=SUPERSAMPLE)
    p.close()
    print(f"[SAVED] {out_path}")


def montage_print(labels, png_paths, out_path, ncols=3):
    """Same 3x2 (six-shape) layout as render_cavern_shapes.render_comparison,
    but with plain dark labels for print rather than the white slide style."""
    panels = []
    for p in png_paths:
        im = Image.open(p).convert("RGBA")
        bb = im.getbbox()
        panels.append(im.crop(bb) if bb else im)

    pw = max(im.width for im in panels)
    ph = max(im.height for im in panels)
    pad = int(0.06 * pw)
    label_h = int(0.14 * ph)
    cell_w, cell_h = pw + 2 * pad, ph + label_h + pad

    scratch = ImageDraw.Draw(Image.new("RGBA", (10, 10)))
    fs = int(0.5 * label_h)
    while fs > 12:
        f = _load_font(fs)
        if max(scratch.textbbox((0, 0), l, font=f)[2] for l in labels) <= 0.94 * cell_w:
            break
        fs -= max(2, fs // 18)
    font = _load_font(fs)

    nrows = (len(panels) + ncols - 1) // ncols
    m = Image.new("RGBA", (cell_w * ncols, cell_h * nrows), (0, 0, 0, 0))
    d = ImageDraw.Draw(m)
    for i, (im, label) in enumerate(zip(panels, labels)):
        r, c = divmod(i, ncols)
        cx = c * cell_w + cell_w // 2
        m.alpha_composite(im, (cx - im.width // 2, r * cell_h + label_h))
        tw = d.textbbox((0, 0), label, font=font)[2]
        d.text((cx - tw // 2, r * cell_h + int(0.12 * label_h)),
               label, fill=(35, 35, 35, 255), font=font)
    m.save(out_path)
    print(f"[SAVED] {out_path}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    labels, png_paths = [], []
    for label, folder in CAVERNS:
        msh_path = os.path.join(GRIDS_ROOT, folder, "geom.msh")
        if not os.path.isfile(msh_path):
            print(f"[SKIP] {label}: {msh_path} not found")
            continue
        key = label.lower().replace(" ", "_").replace("-", "")
        out_path = os.path.join(OUT_DIR, f"domain_overview_{key}.png")
        render_one(msh_path, out_path)
        labels.append(label)
        png_paths.append(out_path)

    montage_print(labels, png_paths,
                  os.path.join(OUT_DIR, "domain_overview_comparison.png"), ncols=3)
    print("\n[DONE]")


if __name__ == "__main__":
    main()
