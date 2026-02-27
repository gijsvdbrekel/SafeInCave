"""
plot_results.py — Combined plotting script for salt cavern simulations.

Merges the functionality of:
  1. plot_convergence.py      — volume convergence + pressure schedule
  2. plot_stress_state.py     — p-q stress paths at 5 probes + dilatancy boundaries
  3. plot_FOS.py              — full-field FOS: global min + roof/mid/floor stats + XDMF
  4. plot_fracture_propagation.py — FOS at radial distances from wall + dilating zone depth

Usage:
  python plot_results.py
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import meshio
import xml.etree.ElementTree as ET

import dolfinx as dfx
from mpi4py import MPI

import safeincave.PostProcessingTools as post
from case_index import detect_layout_and_collect_cases, filter_cases


# =============================================================================
# USER CONFIGURATION
# =============================================================================

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", "Simulation", "output"))

# ── CASE SELECTION ─────────────────────────────────────────────────────────────
# SELECT: Filter which cases to plot (set to None to include all)
#
# Available filters:
#   "caverns"        - Cavern shapes to include: (e.g. "regular1200", "tilted1200","directcirculation600", "reversedcirculation1200")                                                   "fastleached1200", "tubefailure1200").
#   "pressure"       - Pressure scenario: "industry", "transport", "power_generation", "csv", or None
#   "scenario"       - Material scenario(s): "A_SIC", "A_MD", "B_SIC", "B_MD" or None
#   "n_cycles"       - Number of cycles (int) or None
#   "operation_days" - Operation duration (int) or None
#   "case_contains"  - Substring match in case name or None

SELECT = {
    "caverns": ["regular1200"],
    "pressure": ["industry", "power_generation", "csv", "transport"],
    "scenario": None,
    "n_cycles": None,
    "operation_days": None,
    "case_contains": None,
}

# ── PLOT MODE ─────────────────────────────────────────────────────────────────
# PLOT_MODE controls how figures are grouped:
#
#   "compare_shapes"    - One figure per plot type, all cavern shapes overlaid.
#                         Use when running all shapes with the same settings
#                         (same pressure, scenario, cycles).
#                         Color = cavern shape, linestyle = scenario.
#
#   "compare_scenarios" - One figure PER CAVERN SHAPE, with different
#                         scenarios / pressures / cycles overlaid on it.
#                         Use when comparing what-if settings for the same cavern.
#                         Color = scenario, linestyle = scenario.
#
#   "compare_pressures" - One figure PER CAVERN SHAPE, with different
#                         pressure schemes overlaid on it.
#                         Use when comparing pressure schedules for the same
#                         cavern and scenario (e.g. industry vs power_generation).
#                         Color = pressure scheme, linestyle = pressure scheme.

PLOT_MODE = "compare_pressures"    # "compare_shapes", "compare_scenarios", or "compare_pressures"

FIGURES = {
    "convergence": False,         # Figure 1: volume convergence
    "stress_state": False,        # Figure 2: p-q stress paths
    "fos": True,                 # Figure 3: FOS over time
    "fracture_propagation": False,# Figure 4: dilatancy zone analysis
}

# Stress state options
# Available: "ratigan_027", "ratigan_018", "spiers", "devries_comp", "devries_ext"
SHOW_DILATANCY = ["ratigan_027", "spiers", "devries_comp", "devries_ext"]

# FOS plot tuning (rolling quantile bands)
FOS_MAX_POINTS = 1200       # downsample for cleaner plots
FOS_BAND_WINDOW = 31        # rolling window size (odd)
FOS_SHOW_BAND = True        # show P10-P90 band + median
PROBE_ORDER = ["top", "quarter", "mid", "threequarter", "bottom"]

# Fracture propagation options
RADIAL_DISTANCES = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]

# FOS (Factor of Safety) = q_boundary / q
#   FOS < 1.0 means dilating (unsafe)
#   FOS > 1.0 means safe

FOS_THRESHOLD = 1.0

# FOS options
WRITE_FOS_XDMF = False
XDMF_SUBDIR = os.path.join("operation", "_fos_outputs")
XDMF_FILENAME = "fos.xdmf"

OUT_DIR = os.path.join(ROOT, "_figures")
SHOW = True
DPI = 180
CAVERN_PHYS_TAG = 29

# =============================================================================
# END OF USER CONFIGURATION
# =============================================================================


# =============================================================================
# 1. CONSTANTS & STYLE DICTS
# =============================================================================

MPA = 1e6
HOUR = 3600.0
DAY = 24.0 * HOUR

CAVERN_COLORS = {
    "Asymmetric":          "#1f77b4",
    "Direct-circulation":  "#ff7f0e",
    "IrregularFine":       "#d62728",
    "Regular":             "#9467bd",
    "Reversed-circulation":"#8c564b",
    "Tilt":                "#e377c2",
    "Fast-leached":        "#17becf",
    "Tube-failure":        "#bcbd22",
}

SCENARIO_COLORS = {
    # Current: Scenario A/B × SafeInCave/Munson-Dawson
    "A_SIC":            "#1f77b4",   # Scenario A, SafeInCave   (blue)
    "A_MD":             "#ff7f0e",   # Scenario A, Munson-Dawson (orange)
    "B_SIC":            "#2ca02c",   # Scenario B, SafeInCave   (green)
    "B_MD":             "#d62728",   # Scenario B, Munson-Dawson (red)
    # Legacy names (kept for backward compatibility)
    "disloc_old_only":  "#1f77b4",
    "disloc_new_only":  "#ff7f0e",
    "desai_only":       "#2ca02c",
    "full_minus_desai": "#2ca02c",
    "full":             "#d62728",
    "full_minus_ps":    "#d62728",
    "md_only":          "#1f77b4",
    "md_steady_only":   "#bcbd22",
    "full_md":          "#e377c2",
    "interlayer":       "#7f7f7f",
    "nointerlayer":     "#8c564b",
    None:               "#333333",
}

SCENARIO_LINESTYLES = {
    # Current
    "A_SIC":            "-",
    "A_MD":             "--",
    "B_SIC":            "-",
    "B_MD":             "--",
    # Legacy
    "disloc_old_only":  "-",
    "disloc_new_only":  "--",
    "desai_only":       "-.",
    "full_minus_desai": "-",
    "full":             "-",
    "full_minus_ps":    "-",
    "md_only":          "-",
    "md_steady_only":   "--",
    "full_md":          "-",
    "interlayer":       "-",
    "nointerlayer":     "--",
    None:               "-",
}

PRESSURE_COLORS = {
    "industry":         "#1f77b4",   # blue
    "transport":        "#ff7f0e",   # orange
    "power_generation": "#2ca02c",   # green
    "csv":              "#d62728",   # red
    # Legacy
    "sinus":            "#9467bd",   # purple
    "irregular":        "#8c564b",   # brown
    "linear":           "#e377c2",   # pink
    None:               "#333333",
}

PRESSURE_LINESTYLES = {
    "industry":         "-",
    "transport":        "--",
    "power_generation": "-.",
    "csv":              ":",
    # Legacy
    "sinus":            "-",
    "irregular":        "--",
    "linear":           "-.",
    None:               "-",
}

PRESSURE_LABELS = {
    "industry":         "Industry",
    "transport":        "Transport",
    "power_generation": "Power generation",
    "csv":              "CSV profile",
    "sinus":            "Sinusoidal",
    "irregular":        "Irregular",
    "linear":           "Linear",
}

PROBE_COLORS = {
    "top":          "#e41a1c",
    "quarter":      "#377eb8",
    "mid":          "#4daf4a",
    "threequarter": "#984ea3",
    "bottom":       "#ff7f00",
}

CAVERN_ORDER = [
    "Asymmetric", "Direct-circulation", "Regular", "Reversed-circulation",
    "Tilt", "Fast-leached", "Tube-failure", "IrregularFine",
]
SCENARIO_ORDER = [
    # Current
    "A_SIC", "A_MD", "B_SIC", "B_MD",
    # Legacy
    "disloc_old_only", "disloc_new_only", "desai_only", "full_minus_desai",
    "full", "full_minus_ps", "md_only", "md_steady_only", "full_md",
    "interlayer", "nointerlayer", None,
]


# =============================================================================
# 2. PATH HELPERS
# =============================================================================

def path_u_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.xdmf")

def path_u_h5(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.h5")

def path_geom_msh(case_folder):
    return os.path.join(case_folder, "operation", "mesh", "geom.msh")

def path_p_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "p_elems", "p_elems.xdmf")

def path_q_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "q_elems", "q_elems.xdmf")

def path_sig_xdmf(case_folder):
    candidates = [
        os.path.join(case_folder, "operation", "sig", "sig.xdmf"),
        os.path.join(case_folder, "operation", "sig.xdmf"),
    ]
    for p in candidates:
        if os.path.isfile(p):
            return p
    return candidates[0]  # default path even if missing

def path_pressure_json(case_folder):
    return os.path.join(case_folder, "pressure_schedule.json")

def is_valid_hdf5(h5_path: str) -> bool:
    try:
        with open(h5_path, "rb") as f:
            return f.read(8) == b"\x89HDF\r\n\x1a\n"
    except Exception:
        return False

def read_pressure_schedule(case_folder: str):
    pjson = path_pressure_json(case_folder)
    if not os.path.isfile(pjson):
        return None, None
    with open(pjson, "r") as f:
        data = json.load(f)
    if "t_hours" in data and "p_MPa" in data:
        t = np.asarray(data["t_hours"], float)
        p = np.asarray(data["p_MPa"], float)
        return t, p
    if "t_values_s" in data and "p_values_Pa" in data:
        t = np.asarray(data["t_values_s"], float) / HOUR
        p = np.asarray(data["p_values_Pa"], float) / MPA
        return t, p
    if "t_values" in data and "p_values" in data:
        t = np.asarray(data["t_values"], float) / HOUR
        p = np.asarray(data["p_values"], float) / MPA
        return t, p
    return None, None


# =============================================================================
# 3. CASE HELPERS
# =============================================================================

def get_case_color_and_style(cavern_label, scenario_preset, pressure_scenario=None, mode=None):
    """Return (color, linestyle) based on the active PLOT_MODE."""
    if mode is None:
        mode = PLOT_MODE.lower().replace(" ", "_")
    if mode == "compare_pressures":
        color = PRESSURE_COLORS.get(pressure_scenario, "#333333")
        linestyle = PRESSURE_LINESTYLES.get(pressure_scenario, "-")
    elif mode == "compare_scenarios":
        if scenario_preset is not None:
            color = SCENARIO_COLORS.get(scenario_preset, "#333333")
        else:
            color = CAVERN_COLORS.get(cavern_label, "#333333")
        linestyle = SCENARIO_LINESTYLES.get(scenario_preset, "-")
    else:  # compare_shapes
        color = CAVERN_COLORS.get(cavern_label, "#333333")
        linestyle = SCENARIO_LINESTYLES.get(scenario_preset, "-")
    return color, linestyle


def get_case_label(meta, mode=None):
    """Return a human-readable label for a case based on the active PLOT_MODE."""
    if mode is None:
        mode = PLOT_MODE.lower().replace(" ", "_")
    cav = meta.get("cavern_label", "")
    sc = meta.get("scenario_preset")
    ps = meta.get("pressure_scenario")
    if mode == "compare_pressures":
        return PRESSURE_LABELS.get(ps, ps or "unknown")
    elif mode == "compare_scenarios":
        return f"{sc}" if sc is not None else meta.get("case_name", "")
    else:  # compare_shapes
        return f"{cav} | {sc}" if sc is not None else f"{cav}"


def get_series_key(meta):
    """Return a unique series key for a case (used in stress state dicts)."""
    return (meta.get("cavern_label"), meta.get("scenario_preset"), meta.get("pressure_scenario"))


def pick_one_case_per_series(cases_meta):
    """In compare_shapes mode, keep one case per (cavern, scenario, pressure) combo."""
    out = {}
    for m in cases_meta:
        key = (m.get("cavern_label"), m.get("scenario_preset"), m.get("pressure_scenario"))
        if key not in out:
            out[key] = m
    return list(out.values())


def group_cases_by_cavern(cases):
    """Group cases by cavern_label. Returns dict: cavern_label -> [cases]."""
    groups = {}
    for c in cases:
        cav = c.get("cavern_label", "unknown")
        groups.setdefault(cav, []).append(c)
    return groups


def case_has_convergence_files(case_path):
    return (
        os.path.isfile(path_u_xdmf(case_path)) and
        os.path.isfile(path_u_h5(case_path)) and
        os.path.isfile(path_geom_msh(case_path)) and
        is_valid_hdf5(path_u_h5(case_path))
    )

def case_has_stress_files(case_path):
    return all(os.path.isfile(p) for p in [
        path_u_xdmf(case_path),
        path_geom_msh(case_path),
        path_p_xdmf(case_path),
        path_q_xdmf(case_path),
    ])

def case_has_fos_files(case_path):
    sig_ok = os.path.isfile(path_sig_xdmf(case_path))
    return sig_ok and all(os.path.isfile(p) for p in [
        path_geom_msh(case_path),
        path_p_xdmf(case_path),
        path_q_xdmf(case_path),
    ])

def case_has_fracture_files(case_path):
    return all(os.path.isfile(p) for p in [
        path_u_xdmf(case_path),
        path_geom_msh(case_path),
        path_sig_xdmf(case_path),
        path_p_xdmf(case_path),
        path_q_xdmf(case_path),
    ])


# =============================================================================
# 4. GEOMETRY HELPERS
# =============================================================================

def get_wall_indices_from_msh(msh_path):
    msh = meshio.read(msh_path)
    wall_idx = None

    if hasattr(msh, "cells_dict") and "line" in msh.cells_dict:
        wall_idx = np.unique(np.asarray(msh.cells_dict["line"]).reshape(-1))

    if wall_idx is None and isinstance(getattr(msh, "cells", None), dict):
        if "line" in msh.cells:
            wall_idx = np.unique(np.asarray(msh.cells["line"]).reshape(-1))

    if wall_idx is None:
        for cb in msh.cells:
            if getattr(cb, "type", None) == "line":
                wall_idx = np.unique(np.asarray(cb.data).reshape(-1))
                break

    if wall_idx is None:
        raise ValueError(f"No 'line' cells found in {msh_path}")

    return msh.points, wall_idx


def load_wall_points(case_folder):
    u_xdmf = path_u_xdmf(case_folder)
    msh_path = path_geom_msh(case_folder)

    points, _, _ = post.read_node_vector(u_xdmf)
    points_msh, wall_idx_msh = get_wall_indices_from_msh(msh_path)

    mapping = post.build_mapping(points_msh, points)
    wall_idx = np.array([mapping[i] for i in wall_idx_msh], dtype=int)

    wall_points = points[wall_idx]
    order = np.argsort(wall_points[:, 2])
    return wall_points[order]


def auto_generate_probes_from_wall_points(wall_points_sorted_z):
    pts = wall_points_sorted_z
    z = pts[:, 2]
    z_min, z_max = z.min(), z.max()

    fractions = {"bottom": 0.0, "threequarter": 0.25, "mid": 0.5, "quarter": 0.75, "top": 1.0}

    probes = {}
    for name, frac in fractions.items():
        z_target = z_min + frac * (z_max - z_min)
        idx = int(np.argmin(np.abs(z - z_target)))
        probes[name] = pts[idx]

    return probes


def auto_generate_probes_with_index(wall_points_sorted_z):
    """Like auto_generate_probes but also returns wall_idx (for fracture propagation)."""
    pts = wall_points_sorted_z
    z = pts[:, 2]
    z_min, z_max = z.min(), z.max()

    fractions = {"bottom": 0.0, "threequarter": 0.25, "mid": 0.5, "quarter": 0.75, "top": 1.0}

    probes = {}
    for name, frac in fractions.items():
        z_target = z_min + frac * (z_max - z_min)
        idx = int(np.argmin(np.abs(z - z_target)))
        probes[name] = {"wall_idx": idx, "wall_point": pts[idx]}

    return probes


def compute_outward_normal_2d(wall_points):
    n_pts = len(wall_points)
    normals = np.zeros((n_pts, 3))

    for i in range(n_pts):
        if i == 0:
            tangent = wall_points[1] - wall_points[0]
        elif i == n_pts - 1:
            tangent = wall_points[-1] - wall_points[-2]
        else:
            tangent = wall_points[i+1] - wall_points[i-1]

        tangent[1] = 0
        t_len = np.sqrt(tangent[0]**2 + tangent[2]**2)
        if t_len > 1e-10:
            tangent = tangent / t_len

        normal = np.array([tangent[2], 0.0, -tangent[0]])

        if normal[0] < 0:
            normal = -normal

        normals[i] = normal

    return normals


# =============================================================================
# 5. FIGURE 1 FUNCTIONS — Convergence
# =============================================================================

def extract_cavern_facets_from_msh_dolfinx(geom_msh_path: str, cavern_tag: int = 29):
    mesh, cell_tags, facet_tags = dfx.io.gmshio.read_from_msh(geom_msh_path, MPI.COMM_WORLD, 0)
    if facet_tags is None:
        raise RuntimeError("No facet tags were read from the .msh")
    fdim = mesh.topology.dim - 1
    cavern_facets = facet_tags.find(cavern_tag)
    if cavern_facets.size == 0:
        raise RuntimeError(f"No facets found with physical tag {cavern_tag}")
    mesh.topology.create_connectivity(fdim, 0)
    f2v = mesh.topology.connectivity(fdim, 0)
    facets = []
    for f in cavern_facets:
        verts = f2v.links(int(f))
        if len(verts) == 3:
            facets.append(verts)
    facets = np.asarray(facets, dtype=int)
    X = mesh.geometry.x.copy()
    return X, facets


def _area_vectors(points_xyz: np.ndarray, tri: np.ndarray):
    p0 = points_xyz[tri[:, 0]]
    p1 = points_xyz[tri[:, 1]]
    p2 = points_xyz[tri[:, 2]]
    return 0.5 * np.cross(p1 - p0, p2 - p0)


def orient_area_vectors_outward(points: np.ndarray, tris: np.ndarray, area_vecs: np.ndarray):
    center = np.mean(points[np.unique(tris.reshape(-1))], axis=0)
    tri_centroids = (points[tris[:, 0]] + points[tris[:, 1]] + points[tris[:, 2]]) / 3.0
    v = tri_centroids - center
    s = np.einsum("ij,ij->i", area_vecs, v)
    flip = s < 0.0
    area_vecs2 = area_vecs.copy()
    area_vecs2[flip] *= -1.0
    return area_vecs2


def compute_convergence_3d_percent(case_path: str, cavern_phys_tag: int = 29):
    geom_msh = path_geom_msh(case_path)
    u_xdmf = path_u_xdmf(case_path)

    pts_msh, tris_msh = extract_cavern_facets_from_msh_dolfinx(geom_msh, cavern_tag=cavern_phys_tag)

    area_vecs = _area_vectors(pts_msh, tris_msh)
    area_vecs = orient_area_vectors_outward(pts_msh, tris_msh, area_vecs)

    tri_centroids = (pts_msh[tris_msh[:, 0]] + pts_msh[tris_msh[:, 1]] + pts_msh[tris_msh[:, 2]]) / 3.0
    V0 = (1.0 / 3.0) * np.sum(np.einsum("ij,ij->i", tri_centroids, area_vecs))
    V0 = float(abs(V0))
    if not np.isfinite(V0) or V0 <= 0.0:
        raise RuntimeError(f"Invalid V0={V0}")

    points_xdmf, time_list, u_field = post.read_node_vector(u_xdmf)

    mapping = post.build_mapping(pts_msh, points_xdmf)
    tri_xdmf = np.vectorize(mapping.__getitem__)(tris_msh).astype(int)
    v0, v1, v2 = tri_xdmf[:, 0], tri_xdmf[:, 1], tri_xdmf[:, 2]

    u_ref = u_field[0]
    Nt = u_field.shape[0]
    dV = np.zeros(Nt, dtype=float)
    for k in range(Nt):
        u0 = u_field[k, v0, :] - u_ref[v0, :]
        u1 = u_field[k, v1, :] - u_ref[v1, :]
        u2 = u_field[k, v2, :] - u_ref[v2, :]
        uc = (u0 + u1 + u2) / 3.0
        dV[k] = np.sum(np.einsum("ij,ij->i", uc, area_vecs))

    conv_pct = -100.0 * (dV / V0)
    t_days = np.asarray(time_list, float) / DAY
    conv_pct[0] = 0.0
    return t_days, conv_pct


# =============================================================================
# 6. FIGURE 2 FUNCTIONS — Stress state
# =============================================================================

def read_stress_paths(case_folder, probes_dict):
    p_path = path_p_xdmf(case_folder)
    q_path = path_q_xdmf(case_folder)

    points_p, time_list, p_elems = post.read_cell_scalar(p_path)
    _, time_list2, q_elems = post.read_cell_scalar(q_path)

    n = min(len(time_list), len(time_list2))
    time_list = time_list[:n]
    p_elems = p_elems[:n]
    q_elems = q_elems[:n]

    out = {}
    for key, probe_xyz in probes_dict.items():
        idx = post.find_closest_point(probe_xyz, points_p)
        p = -p_elems[:, idx] / MPA
        q =  q_elems[:, idx] / MPA
        out[key] = (p, q)

    return out


def plot_dilatancy_boundaries(ax, show_boundaries=None, p_min=0.01, p_max=40.0, npts=500):
    if show_boundaries is None:
        show_boundaries = ["ratigan_027", "ratigan_018", "spiers", "devries_comp", "devries_ext"]

    p = np.linspace(p_min, p_max, npts)
    I1 = 3.0 * p

    def q_from_sqrtJ2(sqrtJ2):
        return np.sqrt(3.0) * sqrtJ2

    boundary_styles = {
        "ratigan_027":  {"color": "#7570b3", "linestyle": "--", "linewidth": 1.3, "alpha": 0.85},
        "ratigan_018":  {"color": "#7570b3", "linestyle": ":",  "linewidth": 1.3, "alpha": 0.85},
        "spiers":       {"color": "#66a61e", "linestyle": "-.", "linewidth": 1.3, "alpha": 0.85},
        "devries_comp": {"color": "#e7298a", "linestyle": "-",  "linewidth": 1.5, "alpha": 0.90},
        "devries_ext":  {"color": "#e7298a", "linestyle": "--", "linewidth": 1.5, "alpha": 0.90},
    }

    if "ratigan_027" in show_boundaries:
        D = 0.27
        sqrtJ2 = D * I1
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), label="Ratigan 1991 (D=0.27)",
                **boundary_styles["ratigan_027"])

    if "ratigan_018" in show_boundaries:
        D = 0.18
        sqrtJ2 = D * I1
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), label="Ratigan 1991 (D=0.18)",
                **boundary_styles["ratigan_018"])

    if "spiers" in show_boundaries:
        D_sp, b_sp = 0.27, 1.9
        sqrtJ2 = D_sp * I1 + b_sp
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), label="Spiers 1988",
                **boundary_styles["spiers"])

    def devries_q(I1_MPa, psi_rad, D1=0.683, D2=0.512, T0=1.50, m=0.75, sigma_ref=1.0):
        sgn = np.sign(I1_MPa)
        sgn[sgn == 0.0] = 1.0
        denom = (np.sqrt(3.0) * np.cos(psi_rad) - D2 * np.sin(psi_rad))
        sqrtJ2_ = D1 * ((I1_MPa / (sgn * sigma_ref)) ** m) / denom + T0
        return np.sqrt(3.0) * sqrtJ2_

    psi_comp = -np.pi / 6.0
    psi_ext = np.pi / 6.0

    if "devries_comp" in show_boundaries:
        ax.plot(p, devries_q(I1, psi_comp), label="De Vries 2005 (comp)",
                **boundary_styles["devries_comp"])

    if "devries_ext" in show_boundaries:
        ax.plot(p, devries_q(I1, psi_ext), label="De Vries 2005 (ext)",
                **boundary_styles["devries_ext"])


# =============================================================================
# 6b. FOS PLOT HELPERS (rolling bands, downsampling)
# =============================================================================

def _downsample_xy(x, y, max_points=1200):
    """Thin x, y arrays to at most max_points evenly-spaced samples."""
    n = len(x)
    if n <= max_points:
        return np.asarray(x), np.asarray(y)
    idx = np.linspace(0, n - 1, max_points).astype(int)
    idx[0] = 0
    idx[-1] = n - 1
    return np.asarray(x)[idx], np.asarray(y)[idx]


def _rolling_quantile_band(y, window=31, qlo=0.1, qhi=0.9):
    """Rolling quantile band + median for oscillatory series."""
    y = np.asarray(y)
    n = len(y)
    lo = np.empty(n)
    hi = np.empty(n)
    med = np.empty(n)
    half = window // 2
    for i in range(n):
        a = max(0, i - half)
        b = min(n, i + half + 1)
        w = y[a:b]
        lo[i] = np.quantile(w, qlo)
        hi[i] = np.quantile(w, qhi)
        med[i] = np.quantile(w, 0.5)
    return lo, med, hi


def _compute_fos_probes_for_case(c):
    """Compute FOS time series at 5 wall probes + global field minimum.

    Returns (t_days, dict[probe_name -> FOS_array]).
    The dict includes a special key "global_min" with the per-timestep
    minimum FOS across ALL cells in the mesh.
    """
    folder = c["case_path"]
    p_path = path_p_xdmf(folder)
    q_path = path_q_xdmf(folder)
    sig_path = path_sig_xdmf(folder)

    centroids_p, time_list_p, p_elems = post.read_cell_scalar(p_path)
    _, time_list_q, q_elems = post.read_cell_scalar(q_path)
    centroids_sig, time_list_sig, sig33 = post.read_cell_tensor(sig_path)

    n_times = min(len(time_list_p), len(time_list_q), len(time_list_sig))
    time_days = np.array(time_list_p[:n_times]) / DAY
    p_elems = p_elems[:n_times]
    q_elems = q_elems[:n_times]
    sig33 = sig33[:n_times]

    p_MPa = -p_elems / MPA
    q_MPa = q_elems / MPA
    sig33_MPa = -sig33 / MPA

    # --- Global field FOS (all cells) ---
    FOS_field = compute_FOS_field(p_elems, q_elems, sig33, compression_positive=True, q_tol_MPa=1e-3)
    fos_global_min = np.nanmin(np.where(np.isfinite(FOS_field), FOS_field, np.nan), axis=1)

    # --- Probe-based FOS ---
    wall_points = load_wall_points(folder)
    probes = auto_generate_probes_from_wall_points(wall_points)

    fos_by_probe = {}
    for probe_name, probe_xyz in probes.items():
        idx_p = post.find_closest_point(probe_xyz, centroids_p)
        idx_sig = post.find_closest_point(probe_xyz, centroids_sig)

        p_t = p_MPa[:, idx_p]
        q_t = q_MPa[:, idx_p]
        sig_point6 = tensor33_to_voigt6(sig33_MPa[:, idx_sig, :, :])
        psi_t = _psi_from_voigt6(sig_point6)

        q_boundary = q_dil_devries(p_t, psi_t)
        q_safe = np.where(q_t < 1e-3, 1e-3, q_t)
        fos = np.where(q_t < 1e-3, 100.0, q_boundary / q_safe)
        fos_by_probe[probe_name] = fos

    fos_by_probe["global_min"] = fos_global_min

    return time_days, fos_by_probe


# =============================================================================
# 7. FIGURE 3 FUNCTIONS — FOS field
# =============================================================================

def load_p_q(case_folder):
    p_path = path_p_xdmf(case_folder)
    q_path = path_q_xdmf(case_folder)
    if not os.path.isfile(p_path):
        raise FileNotFoundError(f"Missing p_elems for {case_folder}")
    if not os.path.isfile(q_path):
        raise FileNotFoundError(f"Missing q_elems for {case_folder}")
    _, t_p, p_elems = post.read_cell_scalar(p_path)
    _, t_q, q_elems = post.read_cell_scalar(q_path)
    n = min(len(t_p), len(t_q))
    return np.asarray(t_p[:n], float), np.asarray(p_elems[:n], float), np.asarray(q_elems[:n], float)


def load_sig(case_folder):
    sig_path = path_sig_xdmf(case_folder)
    if not os.path.isfile(sig_path):
        raise FileNotFoundError(f"Missing sig for {case_folder}")
    _, t, sig_vals = post.read_cell_tensor(sig_path)
    return np.asarray(t, float), np.asarray(sig_vals, float), sig_path


def _psi_from_tensor33(sig_Pa, compression_positive=True):
    """Compute Lode angle from (nt, nc, 3, 3) stress tensor via eigenvalues."""
    sig = 0.5 * (sig_Pa + np.swapaxes(sig_Pa, -1, -2))
    sig_eff = -sig if compression_positive else sig
    vals = np.linalg.eigvalsh(sig_eff)
    s1, s2, s3 = vals[:, 2], vals[:, 1], vals[:, 0]
    I1 = s1 + s2 + s3
    mean = I1 / 3.0
    s1d, s2d, s3d = s1 - mean, s2 - mean, s3 - mean
    J2 = (1.0/6.0) * ((s1d - s2d)**2 + (s2d - s3d)**2 + (s3d - s1d)**2)
    J3 = s1d * s2d * s3d
    J2_safe = np.maximum(J2, 1e-30)
    x = (3.0*np.sqrt(3.0)/2.0) * (J3 / (J2_safe**1.5))
    x = np.clip(x, -1.0, 1.0)
    theta = (1.0/3.0) * np.arccos(x)
    psi = theta - np.pi/6.0
    return psi


def q_dil_rd(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    p = np.asarray(p_MPa, dtype=float)
    psi = np.asarray(psi, dtype=float)
    I1 = 3.0 * p
    denom = (np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi))
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)
    num = D1 * (np.abs(I1) / sigma_ref) ** m + T0
    sqrtJ2_dil = num / denom
    return np.sqrt(3.0) * sqrtJ2_dil


def compute_FOS_field(p_elems, q_elems, sig_vals, *, compression_positive=True, q_tol_MPa=1e-3):
    """Compute FOS for full field using stored p/q and psi from tensor."""
    nt, nc = p_elems.shape[0], p_elems.shape[1]
    FOS = np.full((nt, nc), np.inf, dtype=float)

    for it in range(nt):
        p_MPa = -p_elems[it] / MPA
        q_MPa = q_elems[it] / MPA

        psi = _psi_from_tensor33(sig_vals[it], compression_positive=compression_positive)

        q_dil = q_dil_rd(p_MPa, psi)

        mask = q_MPa >= float(q_tol_MPa)
        FOS[it, mask] = q_dil[mask] / q_MPa[mask]

    FOS[~np.isfinite(FOS)] = np.inf
    return np.clip(FOS, 0.0, 1e6)


def extract_cavern_wall_cells_slices(geom_msh_path, cavern_phys_tag=29, top_fraction=0.2, bottom_fraction=0.2):
    mesh, cell_tags, facet_tags = dfx.io.gmshio.read_from_msh(geom_msh_path, MPI.COMM_SELF, 0)
    if facet_tags is None:
        raise RuntimeError("No facet tags present in geom.msh")
    dim = mesh.topology.dim
    fdim = dim - 1
    facets = facet_tags.find(cavern_phys_tag)
    if facets.size == 0:
        raise RuntimeError(f"No facets with tag {cavern_phys_tag}")
    mesh.topology.create_connectivity(fdim, 0)
    f2v = mesh.topology.connectivity(fdim, 0)
    X = mesh.geometry.x
    zc = np.empty(facets.size, dtype=float)
    for i, f in enumerate(facets):
        vs = f2v.links(int(f))
        zc[i] = X[vs, 2].mean()
    z_min, z_max = float(zc.min()), float(zc.max())
    H = z_max - z_min
    z_top_thr = z_max - float(top_fraction) * H
    z_bot_thr = z_min + float(bottom_fraction) * H
    roof_mask  = zc >= z_top_thr
    floor_mask = zc <= z_bot_thr
    mid_mask   = ~(roof_mask | floor_mask)
    mesh.topology.create_connectivity(fdim, dim)
    f2c = mesh.topology.connectivity(fdim, dim)
    def cells_from_mask(mask):
        subset = facets[mask]
        cells = set()
        for f in subset:
            for c in f2c.links(int(f)):
                cells.add(int(c))
        return np.array(sorted(cells), dtype=int)
    return {"roof": cells_from_mask(roof_mask), "mid": cells_from_mask(mid_mask), "floor": cells_from_mask(floor_mask)}


def fos_series_stats(FOS, cell_idx):
    Nt = FOS.shape[0]
    out = {"min": np.full(Nt, np.nan), "mean": np.full(Nt, np.nan), "p05": np.full(Nt, np.nan)}
    if cell_idx.size == 0:
        return out
    vals = FOS[:, cell_idx].copy()
    vals[~np.isfinite(vals)] = np.nan
    out["min"]  = np.nanmin(vals, axis=1)
    out["mean"] = np.nanmean(vals, axis=1)
    out["p05"]  = np.nanpercentile(vals, 5.0, axis=1)
    return out


# XDMF mesh helpers

def _infer_xdmf_mesh_grid_name(xdmf_path: str) -> str:
    fallbacks = ["Grid", "mesh", "Mesh", "domain", "Domain"]
    try:
        tree = ET.parse(xdmf_path)
        root = tree.getroot()
        for grid in root.iter():
            if grid.tag.lower().endswith("grid"):
                name = grid.attrib.get("Name") or grid.attrib.get("name")
                if name:
                    return name
    except Exception:
        pass
    return fallbacks[0]


def _read_mesh_from_sig_xdmf(sig_xdmf_path: str):
    grid_name = _infer_xdmf_mesh_grid_name(sig_xdmf_path)
    candidates = [grid_name, "Grid", "mesh", "Mesh", "domain", "Domain"]
    tried = []
    with dfx.io.XDMFFile(MPI.COMM_SELF, sig_xdmf_path, "r") as xf:
        last_err = None
        for name in candidates:
            if name in tried:
                continue
            tried.append(name)
            try:
                mesh = xf.read_mesh(name=name)
                return mesh, name
            except Exception as e:
                last_err = e
                continue
    raise RuntimeError(
        f"Could not read mesh from {sig_xdmf_path}. Tried Grid names: {tried}. "
        f"Last error: {last_err}"
    )


def write_fos_xdmf_from_sig(case_folder: str, sig_xdmf_path: str, time_list: np.ndarray, FOS: np.ndarray):
    out_dir = os.path.join(case_folder, XDMF_SUBDIR)
    os.makedirs(out_dir, exist_ok=True)
    out_xdmf = os.path.join(out_dir, XDMF_FILENAME)

    mesh, used_name = _read_mesh_from_sig_xdmf(sig_xdmf_path)

    V0 = dfx.fem.functionspace(mesh, ("DG", 0))
    fos_fun = dfx.fem.Function(V0)
    fos_fun.name = "FoS"

    ncells = mesh.topology.index_map(mesh.topology.dim).size_local
    if FOS.shape[1] != ncells:
        raise RuntimeError(
            f"[FoS XDMF] Cell count mismatch for case={case_folder}\n"
            f"  FOS has ncells={FOS.shape[1]}\n"
            f"  mesh (from sig.xdmf grid '{used_name}') has ncells={ncells}\n"
            "This means you are NOT aligned with the mesh used to write sig.xdmf."
        )

    FOS_write = np.asarray(FOS, float).copy()
    FOS_write[~np.isfinite(FOS_write)] = np.nan

    with dfx.io.XDMFFile(MPI.COMM_SELF, out_xdmf, "w") as xw:
        xw.write_mesh(mesh)
        for it, t in enumerate(time_list):
            fos_fun.x.array[:] = FOS_write[it, :]
            xw.write_function(fos_fun, t=float(t))

    print("[SAVED XDMF]", out_xdmf, "(+ .h5)")


# =============================================================================
# 8. FIGURE 4 FUNCTIONS — Fracture propagation
# =============================================================================

def _psi_from_voigt6(sig_voigt):
    """Compute Lode angle from (..., 6) Voigt notation [s11,s22,s33,s12,s13,s23]."""
    s11 = sig_voigt[..., 0]
    s22 = sig_voigt[..., 1]
    s33 = sig_voigt[..., 2]
    s12 = sig_voigt[..., 3]
    s13 = sig_voigt[..., 4]
    s23 = sig_voigt[..., 5]

    p = (s11 + s22 + s33) / 3.0
    dev11 = s11 - p
    dev22 = s22 - p
    dev33 = s33 - p

    J2 = 0.5 * (dev11**2 + dev22**2 + dev33**2) + s12**2 + s13**2 + s23**2
    J2 = np.maximum(J2, 1e-30)

    J3 = (dev11 * (dev22 * dev33 - s23**2)
          - s12 * (s12 * dev33 - s23 * s13)
          + s13 * (s12 * s23 - dev22 * s13))

    sqrtJ2 = np.sqrt(J2)
    arg = (3.0 * np.sqrt(3.0) / 2.0) * J3 / (sqrtJ2**3)
    arg = np.clip(arg, -1.0, 1.0)
    psi = -(1.0 / 3.0) * np.arcsin(arg)
    return psi


def q_dil_devries(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    I1 = 3.0 * p_MPa
    denom = np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi)
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)
    num = D1 * (np.abs(I1) / sigma_ref) ** m + T0
    sqrtJ2_dil = num / denom
    return np.sqrt(3.0) * sqrtJ2_dil


def tensor33_to_voigt6(sig33):
    s11 = sig33[..., 0, 0]
    s22 = sig33[..., 1, 1]
    s33 = sig33[..., 2, 2]
    s12 = sig33[..., 0, 1]
    s13 = sig33[..., 0, 2]
    s23 = sig33[..., 1, 2]
    return np.stack([s11, s22, s33, s12, s13, s23], axis=-1)


def compute_FOS_point(p_mpa, q_mpa, psi=None, q_tol_MPa=1e-3):
    """Compute FOS at individual points (for fracture propagation)."""
    if psi is None:
        psi = -np.pi / 6.0
    q_boundary = q_dil_devries(p_mpa, psi)
    q_safe = np.where(q_mpa < q_tol_MPa, q_tol_MPa, q_mpa)
    FOS = q_boundary / q_safe
    FOS = np.where(q_mpa < q_tol_MPa, 100.0, FOS)
    return FOS


def generate_radial_sample_points(wall_points, probes, normals, radial_distances):
    sample_points = {}
    for probe_name, probe_info in probes.items():
        wall_idx = probe_info["wall_idx"]
        wall_pt = probe_info["wall_point"]
        normal = normals[wall_idx]
        points_at_distances = []
        for dist in radial_distances:
            pt = wall_pt + dist * normal
            points_at_distances.append((dist, pt))
        sample_points[probe_name] = points_at_distances
    return sample_points


def read_stress_at_points(case_folder, sample_points):
    p_path = path_p_xdmf(case_folder)
    q_path = path_q_xdmf(case_folder)

    centroids_p, time_list_p, p_elems = post.read_cell_scalar(p_path)
    centroids_q, time_list_q, q_elems = post.read_cell_scalar(q_path)

    sig_path = path_sig_xdmf(case_folder)
    centroids_sig, time_list_sig, sig33 = post.read_cell_tensor(sig_path)

    n_times = min(len(time_list_p), len(time_list_q), len(time_list_sig))
    time_list = time_list_p[:n_times]
    p_elems = p_elems[:n_times]
    q_elems = q_elems[:n_times]
    sig33 = sig33[:n_times]

    time_days = np.array(time_list) / DAY

    p_elems_MPa = -p_elems / MPA
    q_elems_MPa = q_elems / MPA
    sig33_MPa = -sig33 / MPA

    stress_data = {}
    for probe_name, dist_points in sample_points.items():
        stress_data[probe_name] = {}
        for dist, xyz in dist_points:
            idx_p = post.find_closest_point(xyz, centroids_p)
            idx_sig = post.find_closest_point(xyz, centroids_sig)

            p = p_elems_MPa[:, idx_p]
            q = q_elems_MPa[:, idx_p]

            sig_point33 = sig33_MPa[:, idx_sig, :, :]
            sig_point6 = tensor33_to_voigt6(sig_point33)
            psi = _psi_from_voigt6(sig_point6)

            stress_data[probe_name][dist] = (p, q, psi)

    return time_days, stress_data


def compute_FOS_data(time_days, stress_data):
    fos_data = {}
    for probe_name, dist_data in stress_data.items():
        fos_data[probe_name] = {}
        for dist, (p, q, psi) in dist_data.items():
            fos = compute_FOS_point(p, q, psi)
            fos_data[probe_name][dist] = fos
    return fos_data


# =============================================================================
# 9. FIGURE 1 PLOTTING — Convergence
# =============================================================================

def plot_convergence_combined(cases):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 7), sharex=True,
                                   gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12})
    fig.suptitle(f"3D cavern volume convergence | pressure={SELECT.get('pressure')}")

    for c in cases:
        cav = c.get("cavern_label")
        sc = c.get("scenario_preset")
        ps = c.get("pressure_scenario")
        col, ls = get_case_color_and_style(cav, sc, ps)
        label = get_case_label(c)
        try:
            t_days, conv = compute_convergence_3d_percent(c["case_path"], cavern_phys_tag=CAVERN_PHYS_TAG)
        except Exception as e:
            print(f"[SKIP] {cav}/{sc}: {e}")
            continue
        ax1.plot(t_days, conv, linewidth=2.0, color=col, linestyle=ls, alpha=0.95, label=label)

    ax1.set_ylabel("Convergence (ΔV/V0) (%)")
    ax1.grid(True, alpha=0.3)

    h, l = ax1.get_legend_handles_labels()
    uniq = {}
    for hh, ll in zip(h, l):
        if ll not in uniq:
            uniq[ll] = hh
    if uniq:
        ax1.legend(uniq.values(), uniq.keys(), loc="best", fontsize=9, frameon=True)

    # Plot pressure schedule(s) — overlay all distinct schedules in compare_pressures
    _plotted_pressures = set()
    for c in cases:
        ps = c.get("pressure_scenario")
        if ps in _plotted_pressures:
            continue
        _plotted_pressures.add(ps)
        tH, pMPa = read_pressure_schedule(c["case_path"])
        if tH is not None:
            pcol, pls = get_case_color_and_style(c.get("cavern_label"), c.get("scenario_preset"), ps)
            ax2.plot(tH / 24.0, pMPa, linewidth=1.7, color=pcol, linestyle=pls,
                     label=PRESSURE_LABELS.get(ps, ps))
    if not _plotted_pressures:
        ax2.text(0.5, 0.5, "No pressure_schedule.json found.", ha="center", va="center", transform=ax2.transAxes)
    if len(_plotted_pressures) > 1:
        ax2.legend(fontsize=8, frameon=True, loc="best")
    ax2.set_ylabel("Pressure (MPa)")
    ax2.set_xlabel("Time (days)")
    ax2.grid(True, alpha=0.3)

    outname = f"convergence_combined_pressure={SELECT.get('pressure')}_scenario={SELECT.get('scenario')}.png"
    outpath = os.path.join(OUT_DIR, outname.replace(" ", ""))
    fig.savefig(outpath, dpi=DPI)
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close(fig)


def plot_convergence_separate(cases):
    for c in cases:
        cav = c.get("cavern_label")
        sc = c.get("scenario_preset")
        ps = c.get("pressure_scenario")
        col, ls = get_case_color_and_style(cav, sc, ps)
        label = get_case_label(c)

        try:
            t_days, conv = compute_convergence_3d_percent(c["case_path"], cavern_phys_tag=CAVERN_PHYS_TAG)
        except Exception as e:
            print(f"[SKIP] {cav}/{sc}: {e}")
            continue

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True,
                                       gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12})
        fig.suptitle(f"Convergence: {label}")

        ax1.plot(t_days, conv, linewidth=2.0, color=col, alpha=0.95, label=label)
        ax1.set_ylabel("Convergence (ΔV/V0) (%)")
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc="best", fontsize=9, frameon=True)

        tH, pMPa = read_pressure_schedule(c["case_path"])
        if tH is None:
            ax2.text(0.5, 0.5, "No pressure_schedule.json found.", ha="center", va="center", transform=ax2.transAxes)
        else:
            ax2.plot(tH / 24.0, pMPa, linewidth=1.7)
        ax2.set_ylabel("Pressure (MPa)")
        ax2.set_xlabel("Time (days)")
        ax2.grid(True, alpha=0.3)

        safe_name = c.get("case_name", "unknown").replace(" ", "_")
        outname = f"convergence_{safe_name}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


def plot_convergence_per_cavern(cases):
    """One convergence figure per cavern, scenarios/pressures overlaid."""
    groups = group_cases_by_cavern(cases)
    for cav_label, cav_cases in groups.items():
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 7), sharex=True,
                                       gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12})
        fig.suptitle(f"3D cavern volume convergence — {cav_label}")

        for c in cav_cases:
            sc = c.get("scenario_preset")
            ps = c.get("pressure_scenario")
            col, ls = get_case_color_and_style(cav_label, sc, ps)
            label = get_case_label(c)
            try:
                t_days, conv = compute_convergence_3d_percent(c["case_path"], cavern_phys_tag=CAVERN_PHYS_TAG)
            except Exception as e:
                print(f"[SKIP] {cav_label}/{sc}: {e}")
                continue
            ax1.plot(t_days, conv, linewidth=2.0, color=col, linestyle=ls, alpha=0.95, label=label)

        ax1.set_ylabel("Convergence (ΔV/V0) (%)")
        ax1.grid(True, alpha=0.3)
        h, l = ax1.get_legend_handles_labels()
        uniq = {}
        for hh, ll in zip(h, l):
            if ll not in uniq:
                uniq[ll] = hh
        if uniq:
            ax1.legend(uniq.values(), uniq.keys(), loc="best", fontsize=9, frameon=True)

        # Plot pressure schedule(s) — overlay all distinct schedules
        _plotted_pressures = set()
        for c in cav_cases:
            ps = c.get("pressure_scenario")
            if ps in _plotted_pressures:
                continue
            _plotted_pressures.add(ps)
            tH, pMPa = read_pressure_schedule(c["case_path"])
            if tH is not None:
                pcol, pls = get_case_color_and_style(cav_label, c.get("scenario_preset"), ps)
                ax2.plot(tH / 24.0, pMPa, linewidth=1.7, color=pcol, linestyle=pls,
                         label=PRESSURE_LABELS.get(ps, ps))
        if not _plotted_pressures:
            ax2.text(0.5, 0.5, "No pressure_schedule.json found.", ha="center", va="center", transform=ax2.transAxes)
        if len(_plotted_pressures) > 1:
            ax2.legend(fontsize=8, frameon=True, loc="best")
        ax2.set_ylabel("Pressure (MPa)")
        ax2.set_xlabel("Time (days)")
        ax2.grid(True, alpha=0.3)

        safe_cav = cav_label.replace(" ", "_")
        outname = f"convergence_{safe_cav}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


# =============================================================================
# 10. FIGURE 2 PLOTTING — Stress state
# =============================================================================

def plot_stress_combined(cases, stress_by_series):
    probe_types = ["top", "quarter", "mid", "threequarter", "bottom"]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    for i, ptype in enumerate(probe_types):
        ax = axes[i]
        plot_dilatancy_boundaries(ax, show_boundaries=SHOW_DILATANCY)

        for (cav, sc, ps), d in stress_by_series.items():
            p, q = d[ptype]
            color, linestyle = get_case_color_and_style(cav, sc, ps)
            label = get_case_label({"cavern_label": cav, "scenario_preset": sc, "pressure_scenario": ps})

            ax.plot(p, q, linewidth=2.0, color=color, linestyle=linestyle, label=label)
            ax.scatter(p[-1], q[-1], s=30, edgecolors="black", linewidths=0.6,
                       color=color, zorder=5)

        ax.set_title(f"p-q stress path: {ptype}")
        ax.set_xlabel("Mean stress p (MPa)")
        ax.set_ylabel("Von Mises q (MPa)")
        ax.grid(True, alpha=0.3)

    axp = axes[5]
    _plotted_pressures = set()
    for c in cases:
        ps = c.get("pressure_scenario")
        if ps in _plotted_pressures:
            continue
        _plotted_pressures.add(ps)
        tH, pMPa = read_pressure_schedule(c["case_path"])
        if tH is not None:
            pcol, pls = get_case_color_and_style(c.get("cavern_label"), c.get("scenario_preset"), ps)
            axp.plot(tH / 24.0, pMPa, linewidth=2.0, color=pcol, linestyle=pls,
                     label=PRESSURE_LABELS.get(ps, ps))
    if not _plotted_pressures:
        axp.text(0.5, 0.5, "No pressure_schedule.json found.",
                 ha="center", va="center", transform=axp.transAxes)
        axp.axis("off")
    else:
        axp.set_title("Pressure schedule")
        axp.set_xlabel("Time (days)")
        axp.set_ylabel("Pressure (MPa)")
        axp.grid(True, alpha=0.3)
        if len(_plotted_pressures) > 1:
            axp.legend(fontsize=8, frameon=True)

    handles, labels = axes[0].get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h
    fig.legend(uniq.values(), uniq.keys(), loc="upper center", bbox_to_anchor=(0.5, 0.96),
               ncol=min(5, len(uniq)), frameon=True, fontsize=9)

    fig.tight_layout(rect=[0, 0, 1, 0.92])

    outname = f"pq_paths_combined_pressure={SELECT.get('pressure')}_scenario={SELECT.get('scenario')}.png"
    outpath = os.path.join(OUT_DIR, outname.replace(" ", "").replace("'", "").replace("[", "").replace("]", ""))
    fig.savefig(outpath, dpi=DPI)
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close(fig)


def plot_stress_separate(cases, stress_by_series):
    probe_types = ["top", "quarter", "mid", "threequarter", "bottom"]

    for (cav, sc, ps), d in stress_by_series.items():
        fig, axes = plt.subplots(2, 3, figsize=(14, 8))
        axes = axes.flatten()

        case_label = get_case_label({"cavern_label": cav, "scenario_preset": sc, "pressure_scenario": ps})
        color, linestyle = get_case_color_and_style(cav, sc, ps)

        for i, ptype in enumerate(probe_types):
            ax = axes[i]
            plot_dilatancy_boundaries(ax, show_boundaries=SHOW_DILATANCY)

            p, q = d[ptype]
            ax.plot(p, q, linewidth=2.5, color=color, linestyle=linestyle, label=case_label)
            ax.scatter(p[-1], q[-1], s=40, edgecolors="black", linewidths=0.8,
                       color=color, zorder=5)

            ax.set_title(f"p-q stress path: {ptype}")
            ax.set_xlabel("Mean stress p (MPa)")
            ax.set_ylabel("Von Mises q (MPa)")
            ax.grid(True, alpha=0.3)
            if i == 0:
                ax.legend(loc="upper left", fontsize=9, frameon=True)

        axp = axes[5]
        case_path = None
        for m in cases:
            if m.get("cavern_label") == cav and m.get("scenario_preset") == sc and m.get("pressure_scenario") == ps:
                case_path = m["case_path"]
                break

        if case_path:
            tH, pMPa = read_pressure_schedule(case_path)
            if tH is not None:
                axp.plot(tH / 24.0, pMPa, linewidth=2.0, color=color)
                axp.set_title("Pressure schedule")
                axp.set_xlabel("Time (days)")
                axp.set_ylabel("Pressure (MPa)")
                axp.grid(True, alpha=0.3)
            else:
                axp.text(0.5, 0.5, "No pressure data", ha="center", va="center", transform=axp.transAxes)
                axp.axis("off")
        else:
            axp.axis("off")

        fig.suptitle(f"Stress State: {case_label}", fontsize=12, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])

        safe_name = f"{cav}_{sc}_{ps}".replace(" ", "_").replace("/", "_").replace("None", "")
        outname = f"pq_paths_{safe_name}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


def plot_stress_per_cavern(cases):
    """One stress-state figure per cavern, scenarios/pressures overlaid."""
    probe_types = ["top", "quarter", "mid", "threequarter", "bottom"]
    groups = group_cases_by_cavern(cases)

    for cav_label, cav_cases in groups.items():
        # Build stress paths for all cases of this cavern
        stress_by_series = {}
        series_meta = {}
        for m in cav_cases:
            folder = m["case_path"]
            try:
                wall_points = load_wall_points(folder)
                probes = auto_generate_probes_from_wall_points(wall_points)
                series_key = get_series_key(m)
                stress_by_series[series_key] = read_stress_paths(folder, probes)
                series_meta[series_key] = m
            except Exception as e:
                print(f"[WARN] Skipping stress state for {m.get('case_name')}: {e}")
        if not stress_by_series:
            print(f"[STRESS STATE] No valid data for {cav_label}")
            continue

        fig, axes = plt.subplots(2, 3, figsize=(16, 9))
        axes = axes.flatten()

        for i, ptype in enumerate(probe_types):
            ax = axes[i]
            plot_dilatancy_boundaries(ax, show_boundaries=SHOW_DILATANCY)

            for (cav, sc, ps), d in stress_by_series.items():
                p, q = d[ptype]
                color, linestyle = get_case_color_and_style(cav, sc, ps)
                label = get_case_label({"cavern_label": cav, "scenario_preset": sc, "pressure_scenario": ps})

                ax.plot(p, q, linewidth=2.0, color=color, linestyle=linestyle, label=label)
                ax.scatter(p[-1], q[-1], s=30, edgecolors="black", linewidths=0.6,
                           color=color, zorder=5)

            ax.set_title(f"p-q stress path: {ptype}")
            ax.set_xlabel("Mean stress p (MPa)")
            ax.set_ylabel("Von Mises q (MPa)")
            ax.grid(True, alpha=0.3)

        # Pressure subplot — overlay all distinct pressure schedules
        axp = axes[5]
        _plotted_pressures = set()
        for c in cav_cases:
            ps = c.get("pressure_scenario")
            if ps in _plotted_pressures:
                continue
            _plotted_pressures.add(ps)
            tH, pMPa = read_pressure_schedule(c["case_path"])
            if tH is not None:
                pcol, pls = get_case_color_and_style(cav_label, c.get("scenario_preset"), ps)
                axp.plot(tH / 24.0, pMPa, linewidth=2.0, color=pcol, linestyle=pls,
                         label=PRESSURE_LABELS.get(ps, ps))
        if not _plotted_pressures:
            axp.text(0.5, 0.5, "No pressure_schedule.json found.",
                     ha="center", va="center", transform=axp.transAxes)
            axp.axis("off")
        else:
            axp.set_title("Pressure schedule")
            axp.set_xlabel("Time (days)")
            axp.set_ylabel("Pressure (MPa)")
            axp.grid(True, alpha=0.3)
            if len(_plotted_pressures) > 1:
                axp.legend(fontsize=8, frameon=True)

        handles, labels = axes[0].get_legend_handles_labels()
        uniq = {}
        for h, l in zip(handles, labels):
            if l not in uniq:
                uniq[l] = h
        fig.legend(uniq.values(), uniq.keys(), loc="upper center", bbox_to_anchor=(0.5, 0.96),
                   ncol=min(5, len(uniq)), frameon=True, fontsize=9)

        fig.suptitle(f"Stress State — {cav_label}", fontsize=12, fontweight="bold", y=0.99)
        fig.tight_layout(rect=[0, 0, 1, 0.92])

        safe_cav = cav_label.replace(" ", "_")
        outname = f"pq_paths_{safe_cav}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


# =============================================================================
# 11. FIGURE 3 PLOTTING — FOS
# =============================================================================

def _compute_fos_for_case(c):
    """Load data and compute FOS for a single case. Returns dict with all needed arrays."""
    time_list_pq, p_elems, q_elems = load_p_q(c["case_path"])
    time_list_sig, sig_vals, sig_path = load_sig(c["case_path"])

    n_times = min(len(time_list_pq), len(time_list_sig))
    time_list = time_list_pq[:n_times]
    p_elems = p_elems[:n_times]
    q_elems = q_elems[:n_times]
    sig_vals = sig_vals[:n_times]

    FOS = compute_FOS_field(p_elems, q_elems, sig_vals, compression_positive=True, q_tol_MPa=1e-3)

    if WRITE_FOS_XDMF:
        write_fos_xdmf_from_sig(c["case_path"], sig_path, time_list, FOS)

    fos_min_global = np.nanmin(np.where(np.isfinite(FOS), FOS, np.nan), axis=1)

    geom_msh = path_geom_msh(c["case_path"])
    slices = extract_cavern_wall_cells_slices(geom_msh, cavern_phys_tag=CAVERN_PHYS_TAG,
                                              top_fraction=0.20, bottom_fraction=0.20)

    roof_stats  = fos_series_stats(FOS, slices["roof"])
    mid_stats   = fos_series_stats(FOS, slices["mid"])
    floor_stats = fos_series_stats(FOS, slices["floor"])

    t_days = time_list / DAY

    return {
        "t_days": t_days,
        "fos_min_global": fos_min_global,
        "roof_stats": roof_stats,
        "mid_stats": mid_stats,
        "floor_stats": floor_stats,
    }


def _plot_fos_on_axes(axes, fos_by_case, title_prefix=""):
    """
    Plot FOS on a 2x3 grid: 5 probe subplots + 1 min-FOS summary.

    fos_by_case: list of (label, color, linestyle, t_days, fos_by_probe)
    """
    for i, ptype in enumerate(PROBE_ORDER):
        ax = axes[i]
        for label, col, ls, t_days, fos_probes in fos_by_case:
            if ptype not in fos_probes:
                continue
            fos = np.asarray(fos_probes[ptype])
            tx, fy = _downsample_xy(t_days, fos, max_points=FOS_MAX_POINTS)

            if FOS_SHOW_BAND:
                lo, med, hi = _rolling_quantile_band(fy, window=FOS_BAND_WINDOW)
                ax.fill_between(tx, lo, hi, color=col, alpha=0.12, linewidth=0)
                ax.plot(tx, med, linewidth=1.8, linestyle=ls, color=col, label=label)
            else:
                ax.plot(tx, fy, linewidth=1.6, linestyle=ls, color=col, alpha=0.9, label=label)

        ax.axhline(y=1.0, color='red', linestyle=':', linewidth=1.5, alpha=0.8)
        ax.set_title(f"FOS: {ptype}")
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Factor of Safety")
        ax.grid(True, alpha=0.3)

    # --- Summary subplot: global field minimum FOS ---
    ax_sum = axes[5]
    for label, col, ls, t_days, fos_probes in fos_by_case:
        if "global_min" in fos_probes:
            min_fos = fos_probes["global_min"]
        else:
            # Fallback: min across probes
            all_fos = np.array([fos_probes[p] for p in PROBE_ORDER if p in fos_probes])
            min_fos = np.min(all_fos, axis=0)
        tx, my = _downsample_xy(t_days, min_fos, max_points=FOS_MAX_POINTS)

        if FOS_SHOW_BAND:
            lo, med, hi = _rolling_quantile_band(my, window=FOS_BAND_WINDOW)
            ax_sum.fill_between(tx, lo, hi, color=col, alpha=0.12, linewidth=0)
            ax_sum.plot(tx, med, linewidth=2.0, linestyle=ls, color=col, label=label)
        else:
            ax_sum.plot(tx, my, linewidth=2.0, linestyle=ls, color=col, label=label)

    ax_sum.axhline(y=1.0, color='red', linestyle=':', linewidth=1.5, alpha=0.8)
    ax_sum.set_title("Min FOS (global field)")
    ax_sum.set_xlabel("Time (days)")
    ax_sum.set_ylabel("Factor of Safety")
    ax_sum.grid(True, alpha=0.3)
    ax_sum.legend(fontsize=8, frameon=True, loc="best")


def plot_fos_combined(cases):
    """Single 2x3 figure with all cases overlaid (compare_shapes mode)."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    fos_by_case = []
    for c in cases:
        cav = c.get("cavern_label")
        sc = c.get("scenario_preset")
        ps = c.get("pressure_scenario")
        col, ls = get_case_color_and_style(cav, sc, ps)
        label = get_case_label(c)

        try:
            t_days, fos_probes = _compute_fos_probes_for_case(c)
            min_val = min(np.min(v) for v in fos_probes.values())
            print(f"    [FOS] {label}: min FOS = {min_val:.3f}")
            fos_by_case.append((label, col, ls, t_days, fos_probes))
        except Exception as e:
            print(f"[SKIP] {cav}/{sc}: {e}")
            continue

    if not fos_by_case:
        plt.close(fig)
        return

    _plot_fos_on_axes(axes, fos_by_case)
    fig.suptitle(f"Factor of Safety | pressure={SELECT.get('pressure')}", fontsize=12, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    outname = f"fos_combined_pressure={SELECT.get('pressure')}_scenario={SELECT.get('scenario')}.png"
    outpath = os.path.join(OUT_DIR, outname.replace(" ", ""))
    fig.savefig(outpath, dpi=DPI)
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close(fig)


def plot_fos_separate(cases):
    """One 2x3 figure per case."""
    for c in cases:
        cav = c.get("cavern_label")
        sc = c.get("scenario_preset")
        ps = c.get("pressure_scenario")
        col, _ = get_case_color_and_style(cav, sc, ps)
        label = get_case_label(c)

        try:
            t_days, fos_probes = _compute_fos_probes_for_case(c)
        except Exception as e:
            print(f"[SKIP] {cav}/{sc}: {e}")
            continue

        fig, axes = plt.subplots(2, 3, figsize=(16, 9))
        axes = axes.flatten()
        _plot_fos_on_axes(axes, [(label, col, "-", t_days, fos_probes)])
        fig.suptitle(f"Factor of Safety: {label}", fontsize=12, fontweight='bold')
        fig.tight_layout(rect=[0, 0, 1, 0.95])

        safe_name = c.get("case_name", "unknown").replace(" ", "_")
        outname = f"fos_{safe_name}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


def plot_fos_per_cavern(cases):
    """One 2x3 FOS figure per cavern, scenarios/pressures overlaid."""
    groups = group_cases_by_cavern(cases)
    for cav_label, cav_cases in groups.items():
        fig, axes = plt.subplots(2, 3, figsize=(16, 9))
        axes = axes.flatten()

        fos_by_case = []
        for c in cav_cases:
            sc = c.get("scenario_preset")
            ps = c.get("pressure_scenario")
            col, ls = get_case_color_and_style(cav_label, sc, ps)
            label = get_case_label(c)

            try:
                t_days, fos_probes = _compute_fos_probes_for_case(c)
                min_val = min(np.min(v) for v in fos_probes.values())
                print(f"    [FOS] {label}: min FOS = {min_val:.3f}")
                fos_by_case.append((label, col, ls, t_days, fos_probes))
            except Exception as e:
                print(f"[SKIP] {cav_label}/{sc}: {e}")
                continue

        if not fos_by_case:
            plt.close(fig)
            continue

        _plot_fos_on_axes(axes, fos_by_case)
        fig.suptitle(f"Factor of Safety — {cav_label}", fontsize=12, fontweight='bold')
        fig.tight_layout(rect=[0, 0, 1, 0.95])

        safe_cav = cav_label.replace(" ", "_")
        outname = f"fos_{safe_cav}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


# =============================================================================
# 12. FIGURE 4 PLOTTING — Fracture propagation (always per-case)
# =============================================================================

def plot_cavern_with_probes(ax, wall_points, sample_points):
    ax.plot(wall_points[:, 0], wall_points[:, 2], 'k-', linewidth=2, label='Cavern wall')

    for probe_name, dist_points in sample_points.items():
        color = PROBE_COLORS.get(probe_name, 'gray')

        xs = [pt[0] for _, pt in dist_points]
        zs = [pt[2] for _, pt in dist_points]

        ax.plot(xs, zs, '-', color=color, linewidth=1.5, alpha=0.7)

        for i, (dist, pt) in enumerate(dist_points):
            marker = 'o' if dist == 0 else 's'
            size = 80 if dist == 0 else 40
            ax.scatter(pt[0], pt[2], c=color, s=size, marker=marker,
                      edgecolors='black', linewidths=0.5, zorder=5)

        wall_pt = dist_points[0][1]
        ax.annotate(probe_name, (wall_pt[0], wall_pt[2]),
                   textcoords="offset points", xytext=(10, 0),
                   fontsize=9, color=color, fontweight='bold')

    ax.set_xlabel('Radial distance (m)')
    ax.set_ylabel('Height z (m)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_title('Cavern profile with sample points')


def plot_propagation_depth(ax, time_days, fos_data, radial_distances, threshold=1.0):
    for probe_name in fos_data.keys():
        color = PROBE_COLORS.get(probe_name, 'gray')

        max_depth = np.zeros(len(time_days))
        for t_idx in range(len(time_days)):
            depth = 0.0
            for dist in radial_distances:
                fos = fos_data[probe_name][dist][t_idx]
                if fos < threshold:
                    depth = dist
            max_depth[t_idx] = depth

        ax.plot(time_days, max_depth, linewidth=2, color=color, label=probe_name)

    ax.set_xlabel('Time (days)')
    ax.set_ylabel('FOS<1 penetration depth (m)')
    ax.set_title('Connected dilating zone size')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=9)


def plot_fracture_propagation(case_meta):
    """Create fracture propagation figure for a single case (always per-case)."""
    folder = case_meta["case_path"]
    cav = case_meta.get("cavern_label")
    sc = case_meta.get("scenario_preset")
    ps = case_meta.get("pressure_scenario")
    ps_label = PRESSURE_LABELS.get(ps, ps or "unknown")
    case_label = f"{cav} | {sc} | {ps_label}" if sc else f"{cav} | {ps_label}"

    print(f"\n[PROCESSING fracture propagation] {cav} | {sc} | {ps}")

    try:
        wall_points = load_wall_points(folder)
        normals = compute_outward_normal_2d(wall_points)
        probes = auto_generate_probes_with_index(wall_points)

        sample_points = generate_radial_sample_points(
            wall_points, probes, normals, RADIAL_DISTANCES
        )

        time_days, stress_data = read_stress_at_points(folder, sample_points)
        fos_data = compute_FOS_data(time_days, stress_data)

        fig, (ax_cavern, ax_depth) = plt.subplots(1, 2, figsize=(16, 8))

        plot_cavern_with_probes(ax_cavern, wall_points, sample_points)
        plot_propagation_depth(ax_depth, time_days, fos_data, RADIAL_DISTANCES, FOS_THRESHOLD)

        fig.suptitle(f'Dilatancy Zone Analysis: {case_label}\n'
                     f'(De Vries criterion, radial distances: {RADIAL_DISTANCES} m)',
                     fontsize=12, fontweight='bold')
        fig.tight_layout(rect=[0, 0, 1, 0.93])

        safe_name = f"{cav}_{sc}_{ps}".replace(" ", "_").replace("/", "_").replace("None", "none")
        outname = f"fracture_propagation_FOS_{safe_name}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print(f"[SAVED] {outpath}")

        if SHOW:
            plt.show()
        plt.close(fig)

    except Exception as e:
        print(f"[ERROR] Failed to process fracture propagation: {e}")
        import traceback
        traceback.print_exc()


# =============================================================================
# 13. MAIN
# =============================================================================

def main():
    if MPI.COMM_WORLD.rank != 0:
        return

    mode = PLOT_MODE.lower().replace(" ", "_")

    print("=" * 70)
    print("COMBINED PLOT SCRIPT — plot_results.py")
    print("=" * 70)
    print(f"  ROOT:        {ROOT}")
    print(f"  PLOT_MODE:   {PLOT_MODE}")
    print(f"  Caverns:     {SELECT.get('caverns', 'all')}")
    print(f"  Pressure:    {SELECT.get('pressure', 'all')}")
    print(f"  Scenario:    {SELECT.get('scenario', 'all')}")
    print(f"  Contains:    {SELECT.get('case_contains', 'any')}")
    print(f"  Figures:     {[k for k,v in FIGURES.items() if v]}")
    print("=" * 70)

    os.makedirs(OUT_DIR, exist_ok=True)

    all_cases = detect_layout_and_collect_cases(ROOT)
    cases = filter_cases(all_cases, SELECT)
    if not cases:
        print("[DEBUG] Found (first 20):")
        for m in all_cases[:20]:
            print(" -", m.get("cavern_label"), m.get("scenario_preset"),
                  m.get("pressure_scenario"), m.get("case_name"))
        raise RuntimeError(f"No cases matched SELECT={SELECT}")

    # In compare_shapes mode, deduplicate so there is one case per series
    if mode == "compare_shapes":
        cases = pick_one_case_per_series(cases)

    print(f"[INFO] {len(cases)} case(s) after filtering:")
    for m in cases:
        print(f"  - {m.get('cavern_label')} | {m.get('scenario_preset')} | {m.get('case_name')}")

    # --- Figure 1: Convergence ---
    if FIGURES.get("convergence"):
        conv_cases = [c for c in cases if case_has_convergence_files(c["case_path"])]
        if conv_cases:
            print(f"\n[CONVERGENCE] {len(conv_cases)} case(s) with required files")
            if mode == "compare_shapes":
                plot_convergence_combined(conv_cases)
            elif mode in ("compare_scenarios", "compare_pressures"):
                plot_convergence_per_cavern(conv_cases)
            else:
                print(f"[WARN] Unknown PLOT_MODE '{PLOT_MODE}', using compare_shapes")
                plot_convergence_combined(conv_cases)
        else:
            print("[CONVERGENCE] No cases with required files (u.xdmf + geom.msh)")

    # --- Figure 2: Stress state ---
    if FIGURES.get("stress_state"):
        stress_cases = [c for c in cases if case_has_stress_files(c["case_path"])]
        if stress_cases:
            print(f"\n[STRESS STATE] {len(stress_cases)} case(s) with required files")
            if mode == "compare_shapes":
                stress_by_series = {}
                for m in stress_cases:
                    folder = m["case_path"]
                    try:
                        wall_points = load_wall_points(folder)
                        probes = auto_generate_probes_from_wall_points(wall_points)
                        series_key = get_series_key(m)
                        stress_by_series[series_key] = read_stress_paths(folder, probes)
                    except Exception as e:
                        print(f"[WARN] Skipping stress state for {m.get('case_name')}: {e}")
                if stress_by_series:
                    plot_stress_combined(stress_cases, stress_by_series)
                else:
                    print("[STRESS STATE] No valid stress data could be loaded.")
            elif mode in ("compare_scenarios", "compare_pressures"):
                plot_stress_per_cavern(stress_cases)
            else:
                print(f"[WARN] Unknown PLOT_MODE '{PLOT_MODE}', skipping stress state")
        else:
            print("[STRESS STATE] No cases with required files (u.xdmf + geom.msh + p/q_elems)")

    # --- Figure 3: FOS ---
    if FIGURES.get("fos"):
        fos_cases = [c for c in cases if case_has_fos_files(c["case_path"])]
        if fos_cases:
            print(f"\n[FOS] {len(fos_cases)} case(s) with required files")
            if mode == "compare_shapes":
                plot_fos_combined(fos_cases)
            elif mode in ("compare_scenarios", "compare_pressures"):
                plot_fos_per_cavern(fos_cases)
            else:
                print(f"[WARN] Unknown PLOT_MODE '{PLOT_MODE}', using compare_shapes")
                plot_fos_combined(fos_cases)
        else:
            print("[FOS] No cases with required files (geom.msh + p/q_elems + sig.xdmf)")

    # --- Figure 4: Fracture propagation (always per-case) ---
    if FIGURES.get("fracture_propagation"):
        frac_cases = [c for c in cases if case_has_fracture_files(c["case_path"])]
        if frac_cases:
            print(f"\n[FRACTURE PROPAGATION] {len(frac_cases)} case(s) with required files")
            for case_meta in frac_cases:
                plot_fracture_propagation(case_meta)
        else:
            print("[FRACTURE PROPAGATION] No cases with required files (u.xdmf + geom.msh + p/q_elems + sig.xdmf)")

    print("\n[DONE]")


if __name__ == "__main__":
    main()
