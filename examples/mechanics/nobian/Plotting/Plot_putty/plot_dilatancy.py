#!/usr/bin/env python3
import os
import json
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import meshio
from mpi4py import MPI
import dolfinx as dfx

import safeincave.PostProcessingTools as post
from case_index import detect_layout_and_collect_cases, filter_cases, one_case_per_cavern_label


# -------------------------
# Global constants
# -------------------------
hour = 3600.0
day = 24.0 * hour
MPa = 1e6

# =============================================================================
# USER CONFIGURATION - Edit this section to customize the plot
# =============================================================================

# --- Output folder containing simulation results ---
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/nobian/Simulation/output"

# --- Case selection filters ---
# Set any filter to None to include all values for that parameter
SELECT = {
    "caverns": None,                           # e.g. ["Regular", "Tilted"] or None for all
    "pressure": None,                          # "sinus", "linear", "irregular", "csv_profile", or None
    "scenario": None,                          # e.g. ["full", "desai_only"] or None
    "n_cycles": None,                          # e.g. 10, 50, or None for any
    "operation_days": None,                    # e.g. 365, or None for any
    "case_contains": None,                     # substring filter on case name, or None
    "one_case_per_cavern": True,               # only one case per cavern shape
    "t_step": "last",                          # "last" or specific timestep index (int)
}

# --- Color coding by cavern shape ---
CAVERN_COLORS = {
    "Asymmetric":    "#1f77b4",   # blue
    "Irregular":     "#ff7f0e",   # orange
    "IrregularFine": "#d62728",   # red
    "Multichamber":  "#2ca02c",   # green
    "Regular":       "#9467bd",   # purple
    "Teardrop":      "#8c564b",   # brown
    "Tilt":          "#e377c2",   # pink
    "Fastleached":   "#17becf",   # cyan
    "Tubingfailure": "#bcbd22",   # olive
}

# --- Linestyle coding by scenario ---
SCENARIO_LINESTYLES = {
    "disloc_old_only":  "-",
    "disloc_new_only":  "--",
    "desai_only":       "-.",
    "full_minus_desai": ":",
    "full":             (0, (3, 1, 1, 1)),  # dash-dot-dot
    None:               "-",
}

# --- Ordering for legend ---
CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt", "Fastleached", "Tubingfailure", "IrregularFine"]
SCENARIO_ORDER = ["disloc_old_only", "disloc_new_only", "desai_only", "full_minus_desai", "full", None]

# --- Profile extraction settings (IMPORTANT for 3D meshes) ---
PROFILE = {
    # slice axis for building a 2D profile out of 3D cavern surface points
    # "y" usually works well for caverns centered around y=0
    "slice_axis": "y",         # "x" or "y"
    "slice_center": "auto",    # "auto" or float (e.g. 0.0)
    "slice_tol_m": "auto",     # "auto" or float in meters (e.g. 1.0)
    "max_points": 800,         # downsample to keep plots light
}

# --- Other settings ---
CAVERN_PHYS_TAG = 29          # Physical tag for cavern boundary in mesh
OUTDIR = os.path.join(ROOT, "_plots_dashboard")
DPI = 220
SHOW = False                  # Show plot interactively after saving

# =============================================================================
# END OF USER CONFIGURATION
# =============================================================================


def get_cavern_color(cavern_label):
    """Get color for a cavern label, with fallback."""
    return CAVERN_COLORS.get(cavern_label, "#333333")


def get_scenario_linestyle(scenario):
    """Get linestyle for a scenario, with fallback."""
    return SCENARIO_LINESTYLES.get(scenario, "-")


def print_config_summary():
    """Print configuration summary at startup."""
    print("=" * 60)
    print("DILATANCY DASHBOARD PLOT CONFIGURATION")
    print("=" * 60)
    print(f"  ROOT:              {ROOT}")
    print(f"  Time step:         {SELECT.get('t_step', 'last')}")
    print(f"  One per cavern:    {SELECT.get('one_case_per_cavern', True)}")
    print(f"  Caverns:           {SELECT.get('caverns', 'all')}")
    print(f"  Pressure:          {SELECT.get('pressure', 'all')}")
    print(f"  Scenario:          {SELECT.get('scenario', 'all')}")
    print(f"  Contains:          {SELECT.get('case_contains', 'any')}")
    print(f"  Profile slice:     axis={PROFILE['slice_axis']}, center={PROFILE['slice_center']}")
    print("=" * 60)


# -------------------------
# layout helpers
# -------------------------
def apply_theme(fig, axes):
    fig.patch.set_facecolor("white")
    for ax in axes:
        if ax is None:
            continue
        ax.grid(True, alpha=0.25)
        ax.set_axisbelow(True)

def create_layout():
    fig = plt.figure(figsize=(16, 9))
    fig.subplots_adjust(top=0.95, bottom=0.08, left=0.06, right=0.98, hspace=0.44, wspace=0.64)
    gs = GridSpec(18, 19, figure=fig)

    ax_info = fig.add_subplot(gs[0:2, 0:19])
    ax_shape = fig.add_subplot(gs[3:12, 0:4])

    ax00 = fig.add_subplot(gs[3:7, 5:9])
    ax01 = fig.add_subplot(gs[3:7, 10:14])
    ax02 = fig.add_subplot(gs[3:7, 15:])

    ax10 = fig.add_subplot(gs[8:12, 5:9])
    ax11 = fig.add_subplot(gs[8:12, 10:14])
    ax12 = fig.add_subplot(gs[8:12, 15:])

    ax_pressure = fig.add_subplot(gs[14:, 0:7])
    ax_empty1 = fig.add_subplot(gs[14:, 8:13])
    ax_empty2 = fig.add_subplot(gs[14:, 14:19])
    ax_empty1.axis("off")
    ax_empty2.axis("off")

    apply_theme(fig, [ax_info, ax_shape, ax00, ax01, ax02, ax10, ax11, ax12, ax_pressure])
    return fig, ax_info, ax_shape, [ax00, ax01, ax02, ax10, ax11, ax12], ax_pressure


def path_pressure_json(case_folder):
    return os.path.join(case_folder, "pressure_schedule.json")

def read_pressure_schedule(case_path: str):
    pjson = path_pressure_json(case_path)
    if not os.path.isfile(pjson):
        return None, None
    with open(pjson, "r") as f:
        data = json.load(f)

    if "t_hours" in data and "p_MPa" in data:
        return np.asarray(data["t_hours"], float), np.asarray(data["p_MPa"], float)
    if "t_values_s" in data and "p_values_Pa" in data:
        return np.asarray(data["t_values_s"], float) / hour, np.asarray(data["p_values_Pa"], float) / MPa
    if "t_values" in data and "p_values" in data:
        return np.asarray(data["t_values"], float) / hour, np.asarray(data["p_values"], float) / MPa
    return None, None


# -------------------------
# 3D cavern surface -> 2D profile helpers
# -------------------------
def _axis_index(ax: str) -> int:
    ax = str(ax).lower()
    if ax == "x":
        return 0
    if ax == "y":
        return 1
    if ax == "z":
        return 2
    raise ValueError("slice_axis must be 'x','y', or 'z'")

def read_cavern_surface_vertices_from_msh(geom_msh_path: str, cavern_tag: int = 29):
    """
    Returns unique vertex indices on cavern boundary facets (tag=cavern_tag),
    plus the coordinate array X (Nnodes,3).
    """
    mesh, cell_tags, facet_tags = dfx.io.gmshio.read_from_msh(geom_msh_path, MPI.COMM_SELF, 0)
    if facet_tags is None:
        raise RuntimeError("No facet_tags found in geom.msh (cannot locate cavern surface).")

    dim = mesh.topology.dim
    fdim = dim - 1
    cavern_facets = facet_tags.find(cavern_tag)
    if cavern_facets.size == 0:
        raise RuntimeError(f"No facets found with physical tag {cavern_tag} in {geom_msh_path}")

    mesh.topology.create_connectivity(fdim, 0)
    f2v = mesh.topology.connectivity(fdim, 0)

    verts = []
    for f in cavern_facets:
        vs = f2v.links(int(f))
        for v in vs:
            verts.append(int(v))
    verts = np.unique(np.asarray(verts, dtype=int))
    X = mesh.geometry.x.copy()
    return X, verts

def build_profile_indices_from_surface_points(X: np.ndarray, surface_verts: np.ndarray,
                                             *, slice_axis="y", slice_center="auto", slice_tol_m="auto",
                                             max_points=800):
    """
    Pick a thin slice through the cavern surface point cloud.
    Returns selected vertex indices (into X).
    """
    ai = _axis_index(slice_axis)
    pts = X[surface_verts]

    # choose center
    if slice_center == "auto":
        c0 = float(np.median(pts[:, ai]))
    else:
        c0 = float(slice_center)

    # choose tolerance
    if slice_tol_m == "auto":
        # start with 1% of span on that axis, minimum 0.25m
        span = float(np.nanmax(pts[:, ai]) - np.nanmin(pts[:, ai]))
        tol = max(0.25, 0.01 * span)
    else:
        tol = float(slice_tol_m)

    # enlarge tol if too few points
    for _ in range(8):
        mask = np.abs(pts[:, ai] - c0) <= tol
        pick = surface_verts[mask]
        if pick.size >= 60:
            break
        tol *= 1.8

    if pick.size < 10:
        raise RuntimeError(
            f"Could not build a profile slice: only {pick.size} points found. "
            f"Try PROFILE['slice_axis']='x' or set PROFILE['slice_center'] to a better value."
        )

    # downsample (keep extremes in z)
    if pick.size > max_points:
        # keep points sorted by z and take evenly spaced
        z = X[pick, 2]
        order = np.argsort(z)
        pick = pick[order]
        idx = np.linspace(0, pick.size - 1, max_points).astype(int)
        pick = pick[idx]

    return pick, {"slice_axis": slice_axis, "slice_center": c0, "slice_tol_m": tol, "n_points": int(pick.size)}


# -------------------------
# WallProfile + probes (3D-proof)
# -------------------------
class WallProfileData:
    """
    Provides a 2D-ish wall profile in x-z using either:
      - line elements (if present), OR
      - cavern surface slice (3D) if no line elements exist.
    """
    def __init__(self, operation_folder, scale=1.0):
        u_xdmf = os.path.join(operation_folder, "u", "u.xdmf")
        geom_msh = os.path.join(operation_folder, "mesh", "geom.msh")

        points_xdmf, self.time_list, u_field = post.read_node_vector(u_xdmf)  # (Nt,N,3)
        self.points_xdmf = points_xdmf
        self.u_field = u_field
        self.scale = float(scale)

        # --- attempt A: line cells (2D meshes) ---
        wall_idx = None
        try:
            reader_msh = meshio.read(geom_msh)
            points_msh = reader_msh.points

            if hasattr(reader_msh, "cells_dict") and "line" in reader_msh.cells_dict:
                wall_idx = np.unique(reader_msh.cells_dict["line"].flatten())
            else:
                for cell_block in reader_msh.cells:
                    if getattr(cell_block, "type", None) == "line":
                        wall_idx = np.unique(cell_block.data.flatten())
                        break

            if wall_idx is not None and wall_idx.size > 0:
                mapping = post.build_mapping(points_msh, points_xdmf)
                wall_idx = np.asarray([mapping[int(i)] for i in wall_idx], dtype=int)
                self.profile_note = "profile from mesh 'line' cells"
                self.profile_meta = {}
            else:
                wall_idx = None
        except Exception:
            wall_idx = None

        # --- attempt B: 3D cavern surface slice ---
        if wall_idx is None:
            X, surface_verts = read_cavern_surface_vertices_from_msh(geom_msh, cavern_tag=CAVERN_PHYS_TAG)
            pick_msh, meta = build_profile_indices_from_surface_points(
                X, surface_verts,
                slice_axis=PROFILE["slice_axis"],
                slice_center=PROFILE["slice_center"],
                slice_tol_m=PROFILE["slice_tol_m"],
                max_points=PROFILE["max_points"],
            )
            mapping = post.build_mapping(X, points_xdmf)
            wall_idx = np.asarray([mapping[int(i)] for i in pick_msh], dtype=int)
            self.profile_note = "profile from cavern surface slice (3D)"
            self.profile_meta = meta

        # store
        self.wall_idx = wall_idx
        self.wall_points = points_xdmf[wall_idx]             # (Nw,3)
        self.wall_u = u_field[:, wall_idx, :]                # (Nt,Nw,3)

        # sort by z for a clean profile curve
        order = np.argsort(self.wall_points[:, 2])
        self.wall_points = self.wall_points[order]
        self.wall_u = self.wall_u[:, order, :]
        self.wall_idx = self.wall_idx[order]

    def get_wall_coordinates(self, time_step: int):
        return self.wall_points + self.scale * self.wall_u[time_step, :]


def auto_generate_probes(wall_profile: WallProfileData, n_bend_probes: int = 2, min_gap=10):
    # Use x-z curve behavior; points are 3D but we probe based on x-z curvature
    pts = wall_profile.wall_points
    x = pts[:, 0]
    z = pts[:, 2]

    idx_bottom = int(np.argmin(z))
    idx_top = int(np.argmax(z))
    z_mid = 0.5 * (z[idx_bottom] + z[idx_top])
    idx_mid = int(np.argmin(np.abs(z - z_mid)))

    base = [idx_top, idx_mid, idx_bottom]

    dx = np.gradient(x)
    dz = np.gradient(z)
    ddx = np.gradient(dx)
    ddz = np.gradient(dz)

    denom = (dx * dx + dz * dz) ** 1.5
    denom[denom == 0.0] = np.inf
    curvature = np.abs(dx * ddz - dz * ddx) / denom
    curvature[np.isnan(curvature)] = 0.0

    idx_all = np.argsort(curvature)[::-1]

    def too_close(i, existing):
        return any(abs(int(i) - int(e)) < min_gap for e in existing)

    bend = []
    for i in idx_all:
        if len(bend) >= int(n_bend_probes):
            break
        if int(i) in base:
            continue
        if too_close(i, base + bend):
            continue
        bend.append(int(i))

    pick = sorted(base + bend, key=lambda k: z[k], reverse=True)
    return pts[pick]


# -------------------------
# Stress paths
# -------------------------
class StressPath:
    def __init__(self, operation_folder, probes):
        self.probes = probes
        self.points, self.time_list, self.p_elems = post.read_cell_scalar(os.path.join(operation_folder, "p_elems", "p_elems.xdmf"))
        _, _, self.q_elems = post.read_cell_scalar(os.path.join(operation_folder, "q_elems", "q_elems.xdmf"))

        self.p_probes = np.zeros((probes.shape[0], self.time_list.size))
        self.q_probes = np.zeros((probes.shape[0], self.time_list.size))
        for i, probe in enumerate(probes):
            idx = post.find_closest_point(probe, self.points)
            self.p_probes[i, :] = -self.p_elems[:, idx] / MPa
            self.q_probes[i, :] = self.q_elems[:, idx] / MPa


def plot_shape(ax, wall: WallProfileData, t_step: int, probes, cavern_label=None):
    w0 = wall.get_wall_coordinates(0)
    wt = wall.get_wall_coordinates(t_step)

    color = get_cavern_color(cavern_label) if cavern_label else "#1f77b4"
    ax.plot(w0[:, 0], w0[:, 2], "-", linewidth=2, color="#888888", alpha=0.7, label="Initial")
    ax.plot(wt[:, 0], wt[:, 2], "-", linewidth=2, color=color, label=f"t_step={t_step}")

    for i, pr in enumerate(probes):
        ax.scatter(pr[0], pr[2], s=60, edgecolors="black")

    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.axis("equal")
    ax.legend(fontsize=8, frameon=True)

    # annotate extraction mode
    txt = wall.profile_note
    if wall.profile_meta:
        m = wall.profile_meta
        txt += f"\n(slice {m['slice_axis']}={m['slice_center']:.3f} ± {m['slice_tol_m']:.3f} m, N={m['n_points']})"
    ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=7, va="bottom")


def plot_pressure(ax, case_path: str, t_step: int, cavern_label=None):
    tH, pMPa = read_pressure_schedule(case_path)
    if tH is None:
        ax.text(0.5, 0.5, "No pressure_schedule.json", ha="center", va="center", transform=ax.transAxes)
        return
    color = get_cavern_color(cavern_label) if cavern_label else "#1f77b4"
    ax.plot(tH, pMPa, linewidth=2, color=color)
    idx = min(int(t_step), len(tH) - 1)
    ax.scatter([tH[idx]], [pMPa[idx]], s=50, color=color, edgecolors="black", zorder=10)
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Pressure (MPa)")


def plot_stress(ax, sp: StressPath, i_probe: int, t_step: int, cavern_label=None):
    color = get_cavern_color(cavern_label) if cavern_label else "#1f77b4"
    ax.plot(sp.p_probes[i_probe], sp.q_probes[i_probe], linewidth=2, color=color)
    ax.scatter([sp.p_probes[i_probe, t_step]], [sp.q_probes[i_probe, t_step]], s=40, color=color, edgecolors="black", zorder=10)
    ax.set_xlabel("p (MPa)")
    ax.set_ylabel("q (MPa)")
    ax.grid(True, alpha=0.25)


def render_dashboard(case_meta: dict):
    case_path = case_meta["case_path"]
    op = os.path.join(case_path, "operation")

    fig, ax_info, ax_shape, stress_axes, ax_pressure = create_layout()

    wall = WallProfileData(op, scale=1.0)
    Nt = len(wall.time_list)

    tsel = SELECT.get("t_step", "last")
    if isinstance(tsel, str) and tsel.lower() == "last":
        t_step = Nt - 1
    else:
        t_step = int(tsel)
        t_step = max(0, min(t_step, Nt - 1))

    probes = auto_generate_probes(wall, n_bend_probes=2)
    sp = StressPath(op, probes)

    ax_info.axis("off")
    ax_info.text(
        0.01, 0.65,
        f"Case: {case_meta['case_name']} | cavern={case_meta.get('cavern_label')} | "
        f"pressure={case_meta.get('pressure_scenario')} | scenario={case_meta.get('scenario_preset')} | "
        f"t_step={t_step}/{Nt-1}",
        fontsize=11
    )

    plot_shape(ax_shape, wall, t_step, probes, cavern_label=case_meta.get("cavern_label"))
    plot_pressure(ax_pressure, case_path, t_step, cavern_label=case_meta.get("cavern_label"))

    n_use = min(len(stress_axes), probes.shape[0])
    for i in range(n_use):
        plot_stress(stress_axes[i], sp, i, t_step, cavern_label=case_meta.get("cavern_label"))
        stress_axes[i].set_title(f"Probe {i}", fontsize=9)
    for j in range(n_use, len(stress_axes)):
        stress_axes[j].axis("off")

    os.makedirs(OUTDIR, exist_ok=True)
    outname = f"dashboard_{case_meta.get('cavern_label','Unknown')}_{case_meta['case_name']}.png"
    outpath = os.path.join(OUTDIR, outname)
    fig.savefig(outpath, dpi=DPI, bbox_inches="tight")
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close(fig)


def main():
    print_config_summary()

    all_cases = detect_layout_and_collect_cases(ROOT)
    selected = filter_cases(all_cases, SELECT)
    if not selected:
        print("[INFO] Found cases (unfiltered):")
        for m in all_cases[:60]:
            print(" -", m["case_name"], "| pressure:", m.get("pressure_scenario"),
                  "| scenario:", m.get("scenario_preset"), "| cavern:", m.get("cavern_label"))
        raise RuntimeError(f"No cases match SELECT={SELECT}")

    if SELECT.get("one_case_per_cavern", True):
        selected = one_case_per_cavern_label(selected)

    print(f"[INFO] Processing {len(selected)} case(s)...")

    # order by preferred cavern order
    by_label = {c["cavern_label"]: c for c in selected}
    for lab in CAVERN_ORDER:
        if lab in by_label:
            render_dashboard(by_label[lab])

    # extras
    extras = [k for k in by_label.keys() if k not in CAVERN_ORDER]
    for lab in sorted(extras):
        render_dashboard(by_label[lab])


if __name__ == "__main__":
    main()

