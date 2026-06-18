#!/usr/bin/env python3
"""
render_heterogeneous_caverns.py — clean 3D renders of the heterogeneous caverns
(regular cavern crossed by a dipping anhydrite interlayer) for the final-defense
slides.

Renders the two interlayer cases together, in the same style as the cavern-shape
comparison:
    Heterogeneous below   (interlayer below the cavern centre)
    Heterogeneous above   (interlayer above the cavern centre)

Each render shows the cavern (opaque, salt grey) together with the anhydrite
interlayer (semi-transparent slab, clipped to a local box around the cavern so
the slab reads as a band crossing the cavern rather than a domain-wide wall).
Same camera / scale for both, plus a 2-panel named montage.

Off-screen / headless (PuTTY) safe. Reads only the grids, so it runs anywhere.

Usage:
    python render_heterogeneous_caverns.py
"""

import os
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

import numpy as np
import pyvista as pv
import gmsh

import render_cavern_shapes as base   # reuse surface loader, montage, fonts

pv.OFF_SCREEN = True

# =============================================================================
# CONFIGURATION
# =============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GRIDS_ROOT = base.GRIDS_ROOT
OUT_DIR = os.path.join(SCRIPT_DIR, "cavern_renders")

# (display label, grid folder)
CASES = [
    ("Heterogeneous below", "cavern_spike_lower_1200_3D"),
    ("Heterogeneous above", "cavern_spike_upper_1200_3D"),
]

CAVERN_TAG = 29
INTERLAYER_TAG = 32        # "Interlayer_1" volume region

CAV_COLOR = "#cdd7e2"      # salt grey (matches the shape renders)
IL_COLOR = "#c08a3e"       # anhydrite tan
IL_OPACITY = 0.55
IL_BOX_MARGIN = 45.0       # how far (m) the shown interlayer extends past the cavern

WINDOW = (1100, 1500)
SUPERSAMPLE = 2
CAMERA_AZIMUTH = 35.0
CAMERA_ELEVATION = 12.0


# =============================================================================
# GEOMETRY
# =============================================================================

def load_volume_region_surface(msh_path, tag):
    """Boundary surface of a physical volume region (tetrahedra) as PolyData."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.open(msh_path)
        ntags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.asarray(coords).reshape(-1, 3)
        index_of = {int(t): i for i, t in enumerate(ntags)}
        blocks = []
        for e in gmsh.model.getEntitiesForPhysicalGroup(3, tag):
            etypes, _et, enodes = gmsh.model.mesh.getElements(3, e)
            for et, en in zip(etypes, enodes):
                if et == 4:
                    blocks.append(np.asarray(en, dtype=np.int64).reshape(-1, 4))
        tets = np.vstack(blocks)
    finally:
        gmsh.finalize()
    ti = np.vectorize(index_of.get)(tets)
    cells = np.hstack([np.full((len(ti), 1), 4, dtype=np.int64), ti]).ravel()
    ctypes = np.full(len(ti), pv.CellType.TETRA, dtype=np.uint8)
    return pv.UnstructuredGrid(cells, ctypes, coords)


def clip_to_cavern_box(grid, cav_surf, margin=IL_BOX_MARGIN):
    """Keep only the cells whose centroid lies within the cavern x-y footprint
    expanded by `margin`, so the interlayer shows as a local band."""
    b = cav_surf.bounds
    xmin, xmax = b[0] - margin, b[1] + margin
    ymin, ymax = b[2] - margin, b[3] + margin
    c = grid.cell_centers().points
    mask = ((c[:, 0] >= xmin) & (c[:, 0] <= xmax) &
            (c[:, 1] >= ymin) & (c[:, 1] <= ymax))
    return grid.extract_cells(np.where(mask)[0]).extract_surface()


# =============================================================================
# RENDERING
# =============================================================================

def render_one(cav, il, out_path, shared_scale):
    p = pv.Plotter(off_screen=True, window_size=WINDOW)
    p.set_background("white")
    p.add_mesh(il, color=IL_COLOR, opacity=IL_OPACITY, smooth_shading=True,
               specular=0.2, ambient=0.35, diffuse=0.7)
    p.add_mesh(cav, color=CAV_COLOR, smooth_shading=True, specular=0.3,
               specular_power=15, ambient=0.30, diffuse=0.65)
    p.enable_parallel_projection()
    p.camera.focal_point = cav.center
    p.camera_position = "yz"
    p.camera.azimuth = CAMERA_AZIMUTH
    p.camera.elevation = CAMERA_ELEVATION
    p.camera.parallel_scale = shared_scale
    p.enable_anti_aliasing("ssaa")
    p.screenshot(out_path, transparent_background=True, scale=SUPERSAMPLE)
    p.close()
    print(f"[SAVED] {out_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    loaded = []
    for label, folder in CASES:
        msh = os.path.join(GRIDS_ROOT, folder, "geom.msh")
        if not os.path.isfile(msh):
            print(f"[SKIP] {label}: {msh} not found")
            continue
        cav = base.load_cavern_surface(msh, CAVERN_TAG)
        il_vol = load_volume_region_surface(msh, INTERLAYER_TAG)
        il = clip_to_cavern_box(il_vol, cav)
        loaded.append((label, cav, il))
        print(f"[LOADED] {label:22s} cavern {cav.n_points} pts, interlayer band "
              f"{il.n_points} pts")

    if not loaded:
        print("[ERROR] No heterogeneous cases loaded.")
        return

    # Shared scale to include the cavern + the clipped interlayer band.
    max_ext = 0.0
    for _, cav, il in loaded:
        for m in (cav, il):
            b = m.bounds
            max_ext = max(max_ext, b[1]-b[0], b[3]-b[2], b[5]-b[4])
    shared_scale = 0.60 * max_ext

    png_paths, labels = [], []
    for label, cav, il in loaded:
        key = label.lower().replace(" ", "_").replace("-", "")
        path = os.path.join(OUT_DIR, f"cavern_{key}.png")
        render_one(cav, il, path, shared_scale)
        png_paths.append(path)
        labels.append(label)

    base.render_comparison(labels, png_paths,
                           os.path.join(OUT_DIR, "heterogeneous_comparison.png"),
                           ncols=2)
    print("\n[DONE]")


if __name__ == "__main__":
    main()
