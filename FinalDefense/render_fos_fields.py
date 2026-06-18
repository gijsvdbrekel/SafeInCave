#!/usr/bin/env python3
"""
render_fos_fields.py — high-resolution Factor-of-Safety renders for the
final-defense slides.

Replaces the small, low-res ParaView FOS screenshots with consistent high-res
renders, and adds the views you cannot get from a single screenshot:

VIEW options (FOS_VIEW):
  "surface"  cavern wall coloured by FOS (outside view).
  "section"  cavern wall clipped at the mid-plane — look inside the cavern.
  "salt"     vertical cross-section through the salt around the cavern, coloured
             by FOS, so you see how FOS develops into the wall with distance
             (the "distance to cavern wall" view from paraview.py). Limit the
             depth shown with SALT_MAX_DIST_M.

TIME options (TIME_MODE):
  "fos_min"        minimum FOS over the whole operation (worst case per cell).
  "min_pressure"   FOS at the timestep of minimum cavern pressure.
  "timestep"       FOS at TIMESTEP_INDEX (negative indexes from the end).
  "sequence"       render every TIME_STRIDE-th timestep to frames and assemble
                   an animation (GIF always; MP4 if ffmpeg is present) — shows
                   how the FOS field develops over time.

Two separate colour bars, split at FOS = 1: a red bar for the dilatant cells
(FOS < 1) and a blue bar for the safe cells (FOS > 1).

FOS is the thesis definition: De Vries (2005) boundary with the per-cell Lode
angle from the stress tensor, FOS = q_dilatancy(p, psi) / q.

Off-screen / headless (PuTTY) safe. Edit CASES + the options below; the script
works for any cavern geometry / pressure scheme that has operation output.

Usage:
    python render_fos_fields.py
"""

import os
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

import json
import numpy as np
import pyvista as pv
import gmsh
import meshio
from scipy.spatial import cKDTree
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image

import safeincave.PostProcessingTools as post

pv.OFF_SCREEN = True

# =============================================================================
# CONFIGURATION
# =============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_ROOT = os.path.normpath(os.path.join(
    SCRIPT_DIR, "..", "examples", "mechanics", "nobian", "Simulation", "output"))
OUT_DIR = os.path.join(SCRIPT_DIR, "fos_renders")

# (display label, case folder relative to OUTPUT_ROOT or absolute).
CASES = [
    ("Regular", "case_leaching_linear_industry(21)_365days_SB_MD_regular1200"),
]

FOS_VIEW = "salt"          # "surface" | "section" | "salt"
TIME_MODE = "min_pressure" # "fos_min" | "min_pressure" | "timestep" | "sequence"
TIMESTEP_INDEX = -1        # used when TIME_MODE == "timestep"
TIME_STRIDE = 4            # frame stride for TIME_MODE == "sequence" (and fos_min)

CAVERN_TAG = 29
CLIP_AXIS = "y"            # cross-section plane normal for "section"/"salt"
SALT_MAX_DIST_M = 30.0     # for "salt": only show salt within this distance of wall

# FOS colour scale, split at 1 with two separate colour bars.
#   safe  (FOS >= 1): light blue at 1.01  ->  dark blue at 5
#   dilatant (FOS < 1): red at 1  ->  purple at 0.2
FOS_DIL_CLIM = (0.2, 1.0)
FOS_SAFE_CLIM = (1.0, 5.0)
DIL_CMAP = LinearSegmentedColormap.from_list(
    "fos_dil", ["#6a0dad", "#b5179e", "#e5383b"])     # 0.2 purple -> 1.0 red
SAFE_CMAP = LinearSegmentedColormap.from_list(
    "fos_safe", ["#9ecae1", "#4292c6", "#08306b"])    # 1.0 light blue -> 5 dark blue

WINDOW = (1200, 1500)
SUPERSAMPLE = 2
BACKGROUND = "white"
ANIM_FPS = 8

# Camera. Surface view uses a 3D angle; section/salt look onto the cut plane.
CAMERA_AZIMUTH = 35.0
CAMERA_ELEVATION = 12.0

MPA = 1.0e6
DAY_S = 86400.0


# =============================================================================
# FOS COMPUTATION  (thesis definition)
# =============================================================================

def _psi_from_tensor33(sig_Pa, compression_positive=True):
    sig = 0.5 * (sig_Pa + np.swapaxes(sig_Pa, -1, -2))
    sig_eff = -sig if compression_positive else sig
    vals = np.linalg.eigvalsh(sig_eff)
    s1, s2, s3 = vals[:, 2], vals[:, 1], vals[:, 0]
    mean = (s1 + s2 + s3) / 3.0
    s1d, s2d, s3d = s1 - mean, s2 - mean, s3 - mean
    J2 = (1.0/6.0) * ((s1d - s2d)**2 + (s2d - s3d)**2 + (s3d - s1d)**2)
    J3 = s1d * s2d * s3d
    J2_safe = np.maximum(J2, 1e-30)
    x = np.clip((3.0*np.sqrt(3.0)/2.0) * (J3 / (J2_safe**1.5)), -1.0, 1.0)
    return (1.0/3.0) * np.arccos(x) - np.pi/6.0


def _q_dil_devries(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    I1 = 3.0 * np.asarray(p_MPa, float)
    denom = np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi)
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)
    num = D1 * (np.abs(I1) / sigma_ref) ** m + T0
    return np.sqrt(3.0) * (num / denom)


def fos_one_step(p_Pa, q_Pa, sig33, q_tol_MPa=1e-3):
    p_MPa = -p_Pa / MPA
    q_MPa = q_Pa / MPA
    psi = _psi_from_tensor33(sig33, compression_positive=True)
    q_dil = _q_dil_devries(p_MPa, psi)
    fos = np.full(len(p_MPa), np.inf)
    mask = q_MPa >= q_tol_MPa
    fos[mask] = q_dil[mask] / q_MPa[mask]
    return np.clip(fos, 0.0, 1e6)


class CaseData:
    """Loads and holds a case's operation FOS inputs + pressure schedule."""

    def __init__(self, case_dir):
        op = os.path.join(case_dir, "operation")
        self.centroids, t_p, self.p_el = post.read_cell_scalar(
            os.path.join(op, "p_elems", "p_elems.xdmf"))
        _, t_q, self.q_el = post.read_cell_scalar(os.path.join(op, "q_elems", "q_elems.xdmf"))
        _, t_s, self.sig = post.read_cell_tensor(os.path.join(op, "sig", "sig.xdmf"))
        self.nt = min(self.p_el.shape[0], self.q_el.shape[0], self.sig.shape[0])
        self.t_days = np.asarray(t_p[:self.nt], float) / DAY_S
        self.geom_msh = os.path.join(op, "mesh", "geom.msh")
        self.case_dir = case_dir
        self._load_pressure()

    def _load_pressure(self):
        pj = os.path.join(self.case_dir, "pressure_schedule.json")
        self.t_pres = self.p_pres = None
        if os.path.isfile(pj):
            d = json.load(open(pj))
            if "t_hours" in d and "p_MPa" in d:
                self.t_pres = np.asarray(d["t_hours"], float) / 24.0
                self.p_pres = np.asarray(d["p_MPa"], float)

    def fos_at(self, it):
        return fos_one_step(self.p_el[it], self.q_el[it], self.sig[it])

    def fos_min(self, stride=1):
        fmin = np.full(self.p_el.shape[1], np.inf)
        for it in range(0, self.nt, stride):
            np.minimum(fmin, self.fos_at(it), out=fmin)
        return np.clip(fmin, 0.0, 1e6)

    def min_pressure_step(self):
        """Index (into the operation timeline) of minimum cavern pressure."""
        if self.t_pres is None:
            # fall back to the step with the lowest mean stress proxy
            return int(np.argmin([self.p_el[it].mean() for it in range(self.nt)]))
        p_at_t = np.interp(self.t_days, self.t_pres, self.p_pres)
        return int(np.argmin(p_at_t))


# =============================================================================
# GEOMETRY
# =============================================================================

def load_cavern_surface(msh_path, cavern_tag=CAVERN_TAG):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.open(msh_path)
        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.asarray(coords).reshape(-1, 3)
        index_of = {int(t): i for i, t in enumerate(node_tags)}
        ents = gmsh.model.getEntitiesForPhysicalGroup(2, cavern_tag)
        blocks = []
        for e in ents:
            etypes, _et, enodes = gmsh.model.mesh.getElements(2, e)
            for et, en in zip(etypes, enodes):
                if et == 2:
                    blocks.append(np.asarray(en, dtype=np.int64).reshape(-1, 3))
        tris = np.vstack(blocks)
    finally:
        gmsh.finalize()
    tri_idx = np.vectorize(index_of.get)(tris)
    faces = np.hstack([np.full((len(tri_idx), 1), 3, dtype=np.int64), tri_idx]).ravel()
    return pv.PolyData(coords, faces).clean()


def load_volume(msh_path):
    """Full salt volume as a PyVista UnstructuredGrid (tetrahedra)."""
    m = meshio.read(msh_path)
    cells = m.cells_dict if hasattr(m, "cells_dict") else m.cells
    tets = np.asarray(cells["tetra"], dtype=np.int64)
    n = tets.shape[0]
    cell_arr = np.hstack([np.full((n, 1), 4, dtype=np.int64), tets]).ravel()
    celltypes = np.full(n, pv.CellType.TETRA, dtype=np.uint8)
    grid = pv.UnstructuredGrid(cell_arr, celltypes, np.asarray(m.points, float))
    return grid


def nearest_map(target_pts, source_pts):
    return cKDTree(source_pts).query(target_pts)[1]


# =============================================================================
# RENDERING
# =============================================================================

def _camera(plotter, mesh, view):
    b = mesh.bounds
    extent = max(b[1]-b[0], b[3]-b[2], b[5]-b[4])
    scale = 0.58 * extent
    cx, cy, cz = mesh.center
    plotter.enable_parallel_projection()
    # Angled 3/4 view for all views. For the cut-away views ("section", "salt")
    # the near half is removed, so this angle looks into the opened interior.
    plotter.camera.focal_point = (cx, cy, cz)
    plotter.camera_position = "yz"
    plotter.camera.azimuth = CAMERA_AZIMUTH
    plotter.camera.elevation = CAMERA_ELEVATION
    plotter.camera.parallel_scale = scale


def _add_split_fos(plotter, mesh, lighting=True):
    """Add the mesh split at FOS=1 into a red (dilatant) and blue (safe) part,
    each with its own colour bar."""
    fos = mesh.cell_data["FOS"]
    has_dil = np.any(fos < 1.0)
    has_safe = np.any(fos >= 1.0)
    shade = dict(smooth_shading=True, specular=0.15, ambient=0.35, diffuse=0.7) \
        if lighting else dict(lighting=False)

    if has_safe:
        safe = mesh.threshold(1.0, scalars="FOS", invert=False)
        plotter.add_mesh(safe, scalars="FOS", cmap=SAFE_CMAP, clim=FOS_SAFE_CLIM,
                         show_scalar_bar=False, **shade)
        plotter.add_scalar_bar(title="FOS  (safe)", n_labels=5, fmt="%.1f",
                               title_font_size=30, label_font_size=24, color="black",
                               vertical=True, position_x=0.86, position_y=0.50,
                               width=0.09, height=0.42)
    if has_dil:
        dil = mesh.threshold(1.0, scalars="FOS", invert=True)
        plotter.add_mesh(dil, scalars="FOS", cmap=DIL_CMAP, clim=FOS_DIL_CLIM,
                         show_scalar_bar=False, **shade)
        plotter.add_scalar_bar(title="FOS  (dilatant)", n_labels=3, fmt="%.1f",
                               title_font_size=30, label_font_size=24, color="black",
                               vertical=True, position_x=0.86, position_y=0.05,
                               width=0.09, height=0.40)


def build_view_mesh(view, surf, vol, fos_vol, surf_parent, dist_to_wall, clip_origin):
    """Return the mesh to render for the requested view, with cell_data['FOS'].

    surface : full cavern wall.
    section : cavern wall only, clipped at the mid-plane so the 3/4 view looks
              into the interior wall (cut-away, no salt).
    salt    : the near-wall salt shell (cells within SALT_MAX_DIST_M of the wall)
              clipped at the mid-plane and viewed at an angle — so you see the
              cavern wall AND the salt beside it together, like the ParaView
              cut-away: the cut face shows FOS in the salt, and looking into the
              opening shows the cavern's interior wall (red where dilatant).
    """
    if view in ("surface", "section"):
        s = surf.copy()
        s.cell_data["FOS"] = fos_vol[surf_parent]
        if view == "section":
            # Remove the half nearest the angled camera (+x side) so the view
            # looks into the opened interior wall (cut-away).
            s = s.clip(normal="x", origin=clip_origin, invert=True)
        return s
    elif view == "salt":
        v = vol.copy()
        v.cell_data["FOS"] = fos_vol
        v.cell_data["dist"] = dist_to_wall
        if SALT_MAX_DIST_M:
            # Keep only the salt shell within SALT_MAX_DIST_M of the cavern wall.
            v = v.threshold(SALT_MAX_DIST_M, scalars="dist", invert=True)
        # Cut the shell in half and keep the far side, revealing the cavern's
        # interior wall framed by the salt cross-section.
        v = v.clip(normal="x", origin=clip_origin, invert=True)
        return v
    raise ValueError(f"unknown FOS_VIEW {view}")


def render_frame(view_mesh, view, out_path, title=None):
    p = pv.Plotter(off_screen=True, window_size=WINDOW)
    p.set_background(BACKGROUND)
    _add_split_fos(p, view_mesh, lighting=True)
    if title:
        p.add_text(title, position="upper_left", font_size=14, color="black")
    _camera(p, view_mesh, view)
    p.enable_anti_aliasing("ssaa")
    p.screenshot(out_path, scale=SUPERSAMPLE)
    p.close()


def assemble_animation(frame_paths, stem):
    frames = [Image.open(f).convert("RGB") for f in frame_paths]
    gif = os.path.join(OUT_DIR, stem + ".gif")
    frames[0].save(gif, save_all=True, append_images=frames[1:],
                   duration=int(1000 / ANIM_FPS), loop=0)
    print(f"[SAVED] {gif}")
    # MP4 if ffmpeg present
    try:
        import shutil, subprocess, tempfile
        if shutil.which("ffmpeg"):
            mp4 = os.path.join(OUT_DIR, stem + ".mp4")
            with tempfile.TemporaryDirectory() as td:
                for i, im in enumerate(frames):
                    im.save(os.path.join(td, f"f{i:04d}.png"))
                subprocess.run(["ffmpeg", "-y", "-framerate", str(ANIM_FPS),
                                "-i", os.path.join(td, "f%04d.png"),
                                "-pix_fmt", "yuv420p", "-vf",
                                "pad=ceil(iw/2)*2:ceil(ih/2)*2", mp4],
                               check=True, capture_output=True)
            print(f"[SAVED] {mp4}")
        else:
            print("[INFO] ffmpeg not found — MP4 skipped (GIF written).")
    except Exception as e:
        print(f"[INFO] MP4 step skipped: {e}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    for label, case in CASES:
        case_dir = case if os.path.isabs(case) else os.path.join(OUTPUT_ROOT, case)
        if not os.path.isdir(case_dir):
            print(f"[SKIP] {label}: {case_dir} not found")
            continue
        key = label.lower().replace(" ", "_").replace("-", "")
        print(f"[CASE] {label}  view={FOS_VIEW}  time={TIME_MODE}")

        data = CaseData(case_dir)
        surf = load_cavern_surface(data.geom_msh)
        vol = load_volume(data.geom_msh)

        # FOS lives on the volume cells; map to surface faces / build dist field.
        vol_centroids = vol.cell_centers().points
        cell_of_volcell = nearest_map(vol_centroids, data.centroids)   # vol cell -> data cell
        surf_parent_data = nearest_map(surf.cell_centers().points, data.centroids)
        # cavern wall node coords for distance-to-wall
        wall_pts = surf.points
        dist_to_wall = cKDTree(wall_pts).query(vol_centroids)[0]

        clip_origin = (surf.center[0], surf.center[1], surf.center[2])

        if TIME_MODE == "sequence":
            steps = list(range(0, data.nt, TIME_STRIDE))
            frame_paths = []
            tmp = os.path.join(OUT_DIR, f"_frames_{key}_{FOS_VIEW}")
            os.makedirs(tmp, exist_ok=True)
            for fi, it in enumerate(steps):
                fos_cells = data.fos_at(it)
                vm = build_view_mesh(FOS_VIEW, surf, vol, fos_cells[cell_of_volcell],
                                     surf_parent_data, dist_to_wall, clip_origin)
                title = f"{label}   t = {data.t_days[it]:.0f} d"
                fp = os.path.join(tmp, f"frame_{fi:04d}.png")
                render_frame(vm, FOS_VIEW, fp, title=title)
                frame_paths.append(fp)
                print(f"   frame {fi+1}/{len(steps)}  (t={data.t_days[it]:.0f} d)")
            assemble_animation(frame_paths, f"fos_{key}_{FOS_VIEW}_sequence")
        else:
            if TIME_MODE == "fos_min":
                fos_cells = data.fos_min(stride=TIME_STRIDE)
                tag = "fosmin"
            elif TIME_MODE == "min_pressure":
                it = data.min_pressure_step()
                fos_cells = data.fos_at(it)
                tag = f"minP_t{data.t_days[it]:.0f}d"
                print(f"   min-pressure timestep: it={it}  t={data.t_days[it]:.0f} d")
            elif TIME_MODE == "timestep":
                it = TIMESTEP_INDEX % data.nt
                fos_cells = data.fos_at(it)
                tag = f"t{data.t_days[it]:.0f}d"
            else:
                raise ValueError(f"unknown TIME_MODE {TIME_MODE}")

            vm = build_view_mesh(FOS_VIEW, surf, vol, fos_cells[cell_of_volcell],
                                 surf_parent_data, dist_to_wall, clip_origin)
            n_dil = int(np.sum(vm.cell_data["FOS"] < 1.0))
            out = os.path.join(OUT_DIR, f"fos_{key}_{FOS_VIEW}_{tag}.png")
            render_frame(vm, FOS_VIEW, out,
                         title=f"{label}  ({FOS_VIEW})")
            print(f"[SAVED] {out}   (cells FOS<1 shown: {n_dil})")

    print("\n[DONE]")


if __name__ == "__main__":
    main()
