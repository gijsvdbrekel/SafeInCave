#!/usr/bin/env python3
"""
render_cavern_shapes.py — clean, high-resolution 3D renders of the cavern
geometries for the final-defense slides.

Replaces the low-resolution Gmsh mesh-viewer screenshots (cluttered with the
surrounding tetrahedral background mesh) with smooth-shaded "sculptures" of each
cavern, rendered from a single consistent camera/scale so the shapes can be
compared directly.

Reads the cavern wall surface (physical group "Cavern", tag 29) straight from
each grid's geom.msh — no simulation results required, so this runs anywhere the
grids are available.

Off-screen (headless / PuTTY) safe: uses VTK off-screen rendering and writes PNGs.

Outputs (into ./cavern_renders/):
  cavern_<key>.png        one render per shape (transparent background)
  cavern_comparison.png   all shapes in a single row on a white background

Usage:
    python render_cavern_shapes.py
"""

import os
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

import numpy as np
import pyvista as pv
import gmsh
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image, ImageDraw, ImageFont

pv.OFF_SCREEN = True

# Light "salt" top fading to a medium slate at depth — gives depth cueing
# without the near-black bottom of a full grayscale ramp.
SALT_CMAP = LinearSegmentedColormap.from_list(
    "salt_depth", ["#eef2f6", "#c4ceda", "#8b9bad", "#5d6c80"])

# =============================================================================
# CONFIGURATION
# =============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GRIDS_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "grids"))
OUT_DIR = os.path.join(SCRIPT_DIR, "cavern_renders")

CAVERN_TAG = 29   # physical-group tag of the cavern wall surface

# (display label, grid folder) — the six thesis shapes, in presentation order.
CAVERNS = [
    ("Regular",              "cavern_regular_1200_3D"),
    ("Direct-circulation",   "cavern_directcirculation_1200_3D"),
    ("Reversed-circulation", "cavern_reversedcirculation_1200_3D"),
    ("String-failure",       "cavern_tubefailure_1200_3D"),
    ("Tilt",                 "cavern_tilted_1200_3D"),
    ("Fast-leached",         "cavern_fastleached_1200_3D"),
]

# Rendering style.
SURFACE_COLOR = "#cdd7e2"      # cool light grey-blue "salt sculpture"
EDGE_COLOR = None              # set to a hex string to overlay wireframe; None = smooth
BACKGROUND = "white"
WINDOW = (1100, 1500)          # per-shape render size (px)
SUPERSAMPLE = 2               # anti-aliasing scale factor
SHADE_BY_DEPTH = True          # subtle vertical colour gradient for depth cueing

# Single shared camera, applied to every shape so scale/angle are identical.
CAMERA_AZIMUTH = 35.0          # degrees around the vertical axis
CAMERA_ELEVATION = 12.0        # degrees above horizontal


# =============================================================================
# GEOMETRY
# =============================================================================

def load_cavern_surface(msh_path, cavern_tag=CAVERN_TAG):
    """Return the cavern wall surface (physical tag) as a PyVista PolyData."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.open(msh_path)
        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.asarray(coords).reshape(-1, 3)
        index_of = {int(t): i for i, t in enumerate(node_tags)}

        ents = gmsh.model.getEntitiesForPhysicalGroup(2, cavern_tag)
        tri_blocks = []
        for e in ents:
            etypes, _etags, enodes = gmsh.model.mesh.getElements(2, e)
            for et, en in zip(etypes, enodes):
                if et == 2:  # 3-node triangle
                    tri_blocks.append(np.asarray(en, dtype=np.int64).reshape(-1, 3))
        if not tri_blocks:
            raise RuntimeError(f"No cavern triangles (tag {cavern_tag}) in {msh_path}")
        tris = np.vstack(tri_blocks)
    finally:
        gmsh.finalize()

    tri_idx = np.vectorize(index_of.get)(tris)
    faces = np.hstack([np.full((len(tri_idx), 1), 3, dtype=np.int64), tri_idx]).ravel()
    surf = pv.PolyData(coords, faces)
    surf = surf.clean()
    return surf


def prepare_surface(surf):
    """Smooth normals + optional depth scalar for nicer shading."""
    surf = surf.compute_normals(cell_normals=False, point_normals=True,
                                auto_orient_normals=True)
    if SHADE_BY_DEPTH:
        z = surf.points[:, 2]
        surf["depth"] = z.max() - z   # 0 at top, increasing downward
    return surf


# =============================================================================
# RENDERING
# =============================================================================

def _set_shared_camera(plotter, bounds_list):
    """Parallel projection with a single scale derived from the tallest cavern,
    so every shape is drawn at the same physical scale."""
    max_extent = 0.0
    for b in bounds_list:
        dx, dy, dz = b[1] - b[0], b[3] - b[2], b[5] - b[4]
        max_extent = max(max_extent, dx, dy, dz)
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = 0.55 * max_extent


def _add_cavern(plotter, surf):
    kwargs = dict(
        color=SURFACE_COLOR,
        smooth_shading=True,
        specular=0.3,
        specular_power=15,
        ambient=0.30,
        diffuse=0.65,
    )
    if SHADE_BY_DEPTH:
        kwargs.pop("color")
        plotter.add_mesh(surf, scalars="depth", cmap=SALT_CMAP,
                         show_scalar_bar=False, clim=[surf["depth"].min(),
                                                      surf["depth"].max()],
                         **kwargs)
    else:
        plotter.add_mesh(surf, **kwargs)
    if EDGE_COLOR:
        plotter.add_mesh(surf, style="wireframe", color=EDGE_COLOR,
                         line_width=0.5, opacity=0.25)


def render_single(label, surf, out_path, shared_scale):
    p = pv.Plotter(off_screen=True, window_size=WINDOW)
    p.set_background(BACKGROUND)
    _add_cavern(p, surf)
    p.enable_parallel_projection()
    # Center the camera on the cavern and apply the shared scale.
    cx, cy, cz = surf.center
    p.camera.focal_point = (cx, cy, cz)
    p.camera_position = "yz"
    p.camera.azimuth = CAMERA_AZIMUTH
    p.camera.elevation = CAMERA_ELEVATION
    p.camera.parallel_scale = shared_scale
    p.enable_anti_aliasing("ssaa")
    p.screenshot(out_path, transparent_background=True, scale=SUPERSAMPLE)
    p.close()
    print(f"[SAVED] {out_path}")


def _load_font(size):
    for path in ["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                 "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"]:
        if os.path.isfile(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def render_comparison(labels, png_paths, out_path, ncols=3):
    """Composite the individual transparent renders into one grid (default 3
    columns -> 3x2 for six shapes: three on top, three on the bottom), with the
    shape name above each cavern.
    Transparent background, so the montage matches the individual renders.

    Each panel is cropped to its content (alpha bounding box) so the name sits
    directly above its cavern instead of floating over empty margin."""
    panels = []
    for p in png_paths:
        im = Image.open(p).convert("RGBA")
        bb = im.getbbox()          # tight bounds of non-transparent content
        panels.append(im.crop(bb) if bb else im)

    pw = max(im.width for im in panels)
    ph = max(im.height for im in panels)
    pad = int(0.06 * pw)
    label_h = int(0.16 * ph)
    cell_w = pw + 2 * pad
    cell_h = ph + label_h + pad

    # Auto-shrink the label font so the longest name fits the panel width.
    _scratch = ImageDraw.Draw(Image.new("RGBA", (10, 10)))
    font_size = int(0.52 * label_h)
    while font_size > 12:
        f = _load_font(font_size)
        st = max(2, font_size // 16)
        widest = max(_scratch.textbbox((0, 0), lab, font=f, stroke_width=st)[2]
                     for lab in labels)
        if widest <= 0.94 * cell_w:
            break
        font_size -= max(2, font_size // 18)
    font = _load_font(font_size)
    stroke = max(2, font_size // 16)   # heavier strokes -> easier to read

    nrows = (len(panels) + ncols - 1) // ncols
    montage = Image.new("RGBA", (cell_w * ncols, cell_h * nrows), (0, 0, 0, 0))
    draw = ImageDraw.Draw(montage)
    for i, (im, label) in enumerate(zip(panels, labels)):
        r, c = divmod(i, ncols)
        cx = c * cell_w + cell_w // 2
        montage.alpha_composite(im, (cx - im.width // 2, r * cell_h + label_h))
        bbox = draw.textbbox((0, 0), label, font=font, stroke_width=stroke)
        tw = bbox[2] - bbox[0]
        # White text with a thin dark outline so it stays legible on any slide.
        draw.text((cx - tw // 2, r * cell_h + int(0.14 * label_h)),
                  label, fill=(255, 255, 255, 255), font=font,
                  stroke_width=stroke, stroke_fill=(35, 35, 35, 255))
    montage.save(out_path)   # keep alpha -> transparent background
    print(f"[SAVED] {out_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    labels, surfs, bounds = [], [], []
    for label, folder in CAVERNS:
        msh = os.path.join(GRIDS_ROOT, folder, "geom.msh")
        if not os.path.isfile(msh):
            print(f"[SKIP] {label}: {msh} not found")
            continue
        surf = prepare_surface(load_cavern_surface(msh))
        labels.append(label)
        surfs.append(surf)
        bounds.append(surf.bounds)
        print(f"[LOADED] {label:22s} {surf.n_points:6d} pts  "
              f"z-extent {surf.bounds[5]-surf.bounds[4]:.0f} m")

    if not surfs:
        print("[ERROR] No cavern surfaces loaded.")
        return

    # Shared parallel scale from the tallest shape.
    max_extent = max(max(b[1]-b[0], b[3]-b[2], b[5]-b[4]) for b in bounds)
    shared_scale = 0.60 * max_extent

    png_paths = []
    for label, surf in zip(labels, surfs):
        key = label.lower().replace(" ", "_").replace("-", "")
        path = os.path.join(OUT_DIR, f"cavern_{key}.png")
        render_single(label, surf, path, shared_scale)
        png_paths.append(path)

    render_comparison(labels, png_paths,
                      os.path.join(OUT_DIR, "cavern_comparison.png"))
    print("\n[DONE]")


if __name__ == "__main__":
    main()
