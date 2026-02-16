"""
plot_pressure_swing_comparison.py — Compare results across pressure swing amplitudes.

Reads all case_swing_*bar_* folders and produces comparison figures showing
how cavern stability metrics change with daily pressure swing rate (bar/day).

Supports multiple cavern shapes: produces figures per cavern, or set
CAVERN to filter for a specific shape.

Figures per cavern:
  1. Convergence — volume convergence vs time + pressure schedules
  2. Stress state — p-q stress paths at 5 wall probes with dilatancy boundaries
  3. Factor of Safety — FOS vs time at 5 wall probes + min-FOS summary
"""

import os
import re
import json
import numpy as np
import matplotlib.pyplot as plt
import meshio
import dolfinx as dfx
from mpi4py import MPI

import safeincave.PostProcessingTools as post

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                           USER CONFIGURATION                                  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/nobian/Simulation/output"

# CAVERN: Filter for a specific cavern shape, or None to auto-detect all.
#   Use the key as it appears in the folder name (e.g. "regular1200", "tilted1200",
#   "asymmetric1200", "directcirculation1200", "reversedcirculation1200",
#   "fastleached1200", "tubefailure1200").
#   Set to None to process all cavern shapes found.
CAVERN = None

# SWING_BARS: Filter for specific swing values, or None to auto-detect all.
# SWING_BARS = [10, 12, 14, 16, 18, 20]
SWING_BARS = None

# Dilatancy boundaries to show on p-q plots
SHOW_DILATANCY = ["ratigan_027", "spiers", "devries_comp", "devries_ext"]

OUT_DIR = os.path.join(ROOT, "_figures")
SHOW = False
DPI = 180

CAVERN_PHYS_TAG = 29

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                        END OF USER CONFIGURATION                              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

MPA = 1e6
DAY = 24.0 * 3600.0
HOUR = 3600.0

PROBE_ORDER = ["top", "quarter", "mid", "threequarter", "bottom"]


# ══════════════════════════════════════════════════════════════════════════════
# SWING COLOR MAP
# ══════════════════════════════════════════════════════════════════════════════

def build_swing_colors(swing_bars):
    """Build a color map for swing values. Uses a perceptually-ordered colormap."""
    n = len(swing_bars)
    cmap = plt.cm.plasma
    return {b: cmap(i / max(n - 1, 1)) for i, b in enumerate(sorted(swing_bars))}


# ══════════════════════════════════════════════════════════════════════════════
# SWING LINESTYLE HELPER  (NEW)
# ══════════════════════════════════════════════════════════════════════════════

def _swing_linestyle(bar):
    """Make 30 bar/day dashed, others solid."""
    return "--" if bar == 30 else "-"


# ══════════════════════════════════════════════════════════════════════════════
# CAVERN LABEL HELPER
# ══════════════════════════════════════════════════════════════════════════════

def _nice_cavern_label(cavern_key):
    """Convert cavern key (e.g. 'regular1200') to display label."""
    low = (cavern_key or "").lower()
    if low.startswith("asymmetric"):
        return "Asymmetric"
    if low.startswith("directcirculation"):
        return "Direct-circulation"
    if low.startswith("reversedcirculation"):
        return "Reversed-circulation"
    if low.startswith("tilt"):
        return "Tilt"
    if low.startswith("regular"):
        return "Regular"
    if low.startswith("fastleached"):
        return "Fast-leached"
    if low.startswith("tubefailure"):
        return "Tube-failure"
    return cavern_key or "Unknown"


# ══════════════════════════════════════════════════════════════════════════════
# CASE DISCOVERY
# ══════════════════════════════════════════════════════════════════════════════

def discover_swing_cases(root):
    """
    Find all case_swing_*bar_* folders, extract swing amplitude and cavern key.

    Folder pattern: case_swing_<N>bar_<cavern_key>
    Example:        case_swing_10bar_regular1200
    """
    pattern = re.compile(r"case_swing_(\d+)bar_(.+)")
    cases = []

    if not os.path.isdir(root):
        raise FileNotFoundError(f"ROOT not found: {root}")

    for name in sorted(os.listdir(root)):
        m = pattern.match(name)
        if m and os.path.isdir(os.path.join(root, name)):
            swing_bar = int(m.group(1))
            cavern_key = m.group(2)
            cases.append({
                "swing_bar": swing_bar,
                "swing_mpa": swing_bar / 10.0,
                "cavern_key": cavern_key,
                "cavern_label": _nice_cavern_label(cavern_key),
                "case_path": os.path.join(root, name),
                "case_name": name,
            })

    return sorted(cases, key=lambda c: (c["cavern_key"], c["swing_bar"]))


# ══════════════════════════════════════════════════════════════════════════════
# PATH HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def path_u_xdmf(cp):
    return os.path.join(cp, "operation", "u", "u.xdmf")

def path_geom_msh(cp):
    return os.path.join(cp, "operation", "mesh", "geom.msh")

def path_sig_xdmf(cp):
    return os.path.join(cp, "operation", "sig", "sig.xdmf")

def path_p_xdmf(cp):
    return os.path.join(cp, "operation", "p_elems", "p_elems.xdmf")

def path_q_xdmf(cp):
    return os.path.join(cp, "operation", "q_elems", "q_elems.xdmf")


def has_convergence_files(cp):
    return os.path.isfile(path_u_xdmf(cp)) and os.path.isfile(path_geom_msh(cp))

def has_stress_files(cp):
    return all(os.path.isfile(p) for p in [path_p_xdmf(cp), path_q_xdmf(cp), path_geom_msh(cp)])

def has_fos_files(cp):
    return all(os.path.isfile(p) for p in [path_p_xdmf(cp), path_q_xdmf(cp), path_sig_xdmf(cp), path_geom_msh(cp)])


def read_pressure_schedule(case_folder):
    """Returns (t_days, p_MPa) or (None, None)."""
    pjson = os.path.join(case_folder, "pressure_schedule.json")
    if not os.path.isfile(pjson):
        return None, None
    with open(pjson, "r") as f:
        data = json.load(f)
    if "t_hours" in data and "p_MPa" in data:
        t = np.asarray(data["t_hours"], float)
        p = np.asarray(data["p_MPa"], float)
        return t / 24.0, p
    if "t_values_s" in data and "p_values_Pa" in data:
        t = np.asarray(data["t_values_s"], float)
        p = np.asarray(data["p_values_Pa"], float)
        return t / DAY, p / MPA
    return None, None


# ══════════════════════════════════════════════════════════════════════════════
# GEOMETRY HELPERS
# ══════════════════════════════════════════════════════════════════════════════

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


def auto_generate_probes(wall_points_sorted_z):
    """Generate probe locations at evenly-spaced heights along the wall."""
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


# ══════════════════════════════════════════════════════════════════════════════
# CONVERGENCE COMPUTATION
# ══════════════════════════════════════════════════════════════════════════════

def _area_vectors(points_xyz, tri):
    p0 = points_xyz[tri[:, 0]]
    p1 = points_xyz[tri[:, 1]]
    p2 = points_xyz[tri[:, 2]]
    return 0.5 * np.cross(p1 - p0, p2 - p0)


def _orient_area_vectors_outward(points, tris, area_vecs):
    center = np.mean(points[np.unique(tris.reshape(-1))], axis=0)
    tri_centroids = (points[tris[:, 0]] + points[tris[:, 1]] + points[tris[:, 2]]) / 3.0
    v = tri_centroids - center
    s = np.einsum("ij,ij->i", area_vecs, v)
    flip = s < 0.0
    area_vecs2 = area_vecs.copy()
    area_vecs2[flip] *= -1.0
    return area_vecs2


def compute_convergence(case_path, cavern_tag=29):
    geom_msh = path_geom_msh(case_path)
    u_xdmf = path_u_xdmf(case_path)

    mesh, cell_tags, facet_tags = dfx.io.gmshio.read_from_msh(geom_msh, MPI.COMM_WORLD, 0)
    fdim = mesh.topology.dim - 1
    cavern_facets = facet_tags.find(cavern_tag)
    mesh.topology.create_connectivity(fdim, 0)
    f2v = mesh.topology.connectivity(fdim, 0)

    facets = []
    for f in cavern_facets:
        verts = f2v.links(int(f))
        if len(verts) == 3:
            facets.append(verts)
    tris = np.asarray(facets, dtype=int)
    pts_msh = mesh.geometry.x.copy()

    area_vecs = _area_vectors(pts_msh, tris)
    area_vecs = _orient_area_vectors_outward(pts_msh, tris, area_vecs)

    tri_centroids = (pts_msh[tris[:, 0]] + pts_msh[tris[:, 1]] + pts_msh[tris[:, 2]]) / 3.0
    V0 = abs((1.0 / 3.0) * np.sum(np.einsum("ij,ij->i", tri_centroids, area_vecs)))

    points_xdmf, time_list, u_field = post.read_node_vector(u_xdmf)
    mapping = post.build_mapping(pts_msh, points_xdmf)
    tri_xdmf = np.vectorize(mapping.__getitem__)(tris).astype(int)
    v0, v1, v2 = tri_xdmf[:, 0], tri_xdmf[:, 1], tri_xdmf[:, 2]

    u_ref = u_field[0]
    Nt = u_field.shape[0]
    dV = np.zeros(Nt, float)
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


# ══════════════════════════════════════════════════════════════════════════════
# STRESS PATH READING
# ══════════════════════════════════════════════════════════════════════════════

def read_stress_paths(case_folder, probes_dict):
    """
    Read p-q stress paths at probe locations.

    Returns: dict[probe_name] -> (p_array_MPa, q_array_MPa)
    """
    p_path = path_p_xdmf(case_folder)
    q_path = path_q_xdmf(case_folder)

    points_p, time_list, p_elems = post.read_cell_scalar(p_path)
    _, time_list2, q_elems = post.read_cell_scalar(q_path)

    n = min(len(time_list), len(time_list2))
    p_elems = p_elems[:n]
    q_elems = q_elems[:n]

    out = {}
    for key, probe_xyz in probes_dict.items():
        idx = post.find_closest_point(probe_xyz, points_p)
        p = -p_elems[:, idx] / MPA
        q = q_elems[:, idx] / MPA
        out[key] = (p, q)

    return out


# ══════════════════════════════════════════════════════════════════════════════
# FOS COMPUTATION
# ══════════════════════════════════════════════════════════════════════════════

def compute_psi_from_sigma(sig_voigt):
    s11, s22, s33 = sig_voigt[..., 0], sig_voigt[..., 1], sig_voigt[..., 2]
    s12, s13, s23 = sig_voigt[..., 3], sig_voigt[..., 4], sig_voigt[..., 5]
    p = (s11 + s22 + s33) / 3.0
    dev11, dev22, dev33 = s11 - p, s22 - p, s33 - p
    J2 = 0.5 * (dev11**2 + dev22**2 + dev33**2) + s12**2 + s13**2 + s23**2
    J2 = np.maximum(J2, 1e-30)
    J3 = (dev11 * (dev22 * dev33 - s23**2)
          - s12 * (s12 * dev33 - s23 * s13)
          + s13 * (s12 * s23 - dev22 * s13))
    sqrtJ2 = np.sqrt(J2)
    arg = np.clip((3.0 * np.sqrt(3.0) / 2.0) * J3 / (sqrtJ2**3), -1.0, 1.0)
    return (1.0 / 3.0) * np.arcsin(arg)


def q_dil_devries(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    I1 = 3.0 * p_MPa
    denom = np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi)
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)
    num = D1 * (np.abs(I1) / sigma_ref) ** m + T0
    return np.sqrt(3.0) * num / denom


def tensor33_to_voigt6(sig33):
    return np.stack([sig33[..., 0, 0], sig33[..., 1, 1], sig33[..., 2, 2],
                     sig33[..., 0, 1], sig33[..., 0, 2], sig33[..., 1, 2]], axis=-1)


def compute_FOS_per_probe(case_folder, probes_dict):
    """
    Compute FOS time series at each probe location.

    Returns: t_days, dict[probe_name] -> FOS_array
    """
    p_path = path_p_xdmf(case_folder)
    q_path = path_q_xdmf(case_folder)
    sig_path = path_sig_xdmf(case_folder)

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

    fos_by_probe = {}
    for probe_name, probe_xyz in probes_dict.items():
        idx_p = post.find_closest_point(probe_xyz, centroids_p)
        idx_sig = post.find_closest_point(probe_xyz, centroids_sig)

        p_t = p_MPa[:, idx_p]
        q_t = q_MPa[:, idx_p]
        sig_point6 = tensor33_to_voigt6(sig33_MPa[:, idx_sig, :, :])
        psi_t = compute_psi_from_sigma(sig_point6)

        q_boundary = q_dil_devries(p_t, psi_t)
        q_safe = np.where(q_t < 1e-3, 1e-3, q_t)
        fos = np.where(q_t < 1e-3, 100.0, q_boundary / q_safe)
        fos_by_probe[probe_name] = fos

    return time_days, fos_by_probe


# ══════════════════════════════════════════════════════════════════════════════
# DILATANCY BOUNDARIES (from plot_stress_state.py)
# ══════════════════════════════════════════════════════════════════════════════

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
        sqrtJ2 = 0.27 * I1
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), label="Ratigan 1991 (D=0.27)", **boundary_styles["ratigan_027"])

    if "ratigan_018" in show_boundaries:
        sqrtJ2 = 0.18 * I1
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), label="Ratigan 1991 (D=0.18)", **boundary_styles["ratigan_018"])

    if "spiers" in show_boundaries:
        sqrtJ2 = 0.27 * I1 + 1.9
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), label="Spiers 1988", **boundary_styles["spiers"])

    def devries_q(I1_MPa, psi_rad, D1=0.683, D2=0.512, T0=1.50, m=0.75, sigma_ref=1.0):
        sgn = np.sign(I1_MPa)
        sgn[sgn == 0.0] = 1.0
        denom = (np.sqrt(3.0) * np.cos(psi_rad) - D2 * np.sin(psi_rad))
        sqrtJ2_ = D1 * ((I1_MPa / (sgn * sigma_ref)) ** m) / denom + T0
        return np.sqrt(3.0) * sqrtJ2_

    if "devries_comp" in show_boundaries:
        ax.plot(p, devries_q(I1, -np.pi / 6.0), label="De Vries 2005 (comp)", **boundary_styles["devries_comp"])

    if "devries_ext" in show_boundaries:
        ax.plot(p, devries_q(I1, np.pi / 6.0), label="De Vries 2005 (ext)", **boundary_styles["devries_ext"])


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: CONVERGENCE
# ══════════════════════════════════════════════════════════════════════════════

def plot_fig_convergence(cavern_key, cavern_label, cases, swing_colors):
    """
    Two-panel figure: top = convergence lines, bottom = pressure schedules.
    One line per swing value, same colors.
    """
    fig, (ax_conv, ax_pres) = plt.subplots(
        2, 1, figsize=(13, 7), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12})

    for c in cases:
        bar = c["swing_bar"]
        color = swing_colors[bar]
        cp = c["case_path"]

        if not has_convergence_files(cp):
            continue

        try:
            t, conv = compute_convergence(cp, cavern_tag=CAVERN_PHYS_TAG)
            ax_conv.plot(
                t, conv,
                linewidth=2,
                linestyle=_swing_linestyle(bar),   # NEW
                color=color,
                alpha=0.9,
                label=f"{bar} bar/day"
            )
            print(f"    [CONV] {bar} bar: {conv[-1]:.4f}% at t={t[-1]:.1f} d")
        except Exception as e:
            print(f"    [WARN] Convergence {bar} bar: {e}")

        t_d, p_mpa = read_pressure_schedule(cp)
        if t_d is not None:
            ax_pres.plot(
                t_d, p_mpa,
                linewidth=1.2,
                linestyle=_swing_linestyle(bar),   # NEW
                color=color,
                alpha=0.7
            )

    ax_conv.set_ylabel("Convergence (dV/V0) (%)")
    ax_conv.grid(True, alpha=0.3)
    ax_conv.legend(fontsize=9, frameon=True, loc="best")

    ax_pres.set_ylabel("Pressure (MPa)")
    ax_pres.set_xlabel("Time (days)")
    ax_pres.grid(True, alpha=0.3)

    fig.suptitle(f"Volume Convergence — {cavern_label}", fontsize=13, fontweight='bold')

    safe_key = cavern_key.replace(" ", "_")
    outpath = os.path.join(OUT_DIR, f"swing_convergence_{safe_key}.png")
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight')
    print(f"  [SAVED] {outpath}")

    if SHOW:
        plt.show()
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: STRESS STATE (p-q paths)
# ══════════════════════════════════════════════════════════════════════════════

def plot_fig_stress_state(cavern_key, cavern_label, cases, swing_colors):
    """
    2x3 figure like plot_stress_state.py: 5 probe subplots + pressure schedule.
    Each subplot overlays p-q paths for all swing values with dilatancy boundaries.
    """
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    stress_by_swing = {}
    for c in cases:
        bar = c["swing_bar"]
        cp = c["case_path"]
        if not has_stress_files(cp):
            continue
        try:
            wall_points = load_wall_points(cp)
            probes = auto_generate_probes(wall_points)
            stress_by_swing[bar] = read_stress_paths(cp, probes)
            print(f"    [STRESS] {bar} bar: loaded")
        except Exception as e:
            print(f"    [WARN] Stress {bar} bar: {e}")

    if not stress_by_swing:
        plt.close(fig)
        print("  [SKIP] No stress data available")
        return

    for i, ptype in enumerate(PROBE_ORDER):
        ax = axes[i]
        plot_dilatancy_boundaries(ax, show_boundaries=SHOW_DILATANCY)

        for bar in sorted(stress_by_swing.keys()):
            d = stress_by_swing[bar]
            if ptype not in d:
                continue
            p, q = d[ptype]
            color = swing_colors[bar]
            ax.plot(
                p, q,
                linewidth=2.0,
                linestyle=_swing_linestyle(bar),   # NEW
                color=color,
                label=f"{bar} bar/day"
            )
            ax.scatter(p[-1], q[-1], s=30, edgecolors="black", linewidths=0.6,
                       color=color, zorder=5)

        ax.set_title(f"p-q: {ptype}")
        ax.set_xlabel("Mean stress p (MPa)")
        ax.set_ylabel("Von Mises q (MPa)")
        ax.grid(True, alpha=0.3)

    axp = axes[5]
    for c in cases:
        bar = c["swing_bar"]
        t_d, p_mpa = read_pressure_schedule(c["case_path"])
        if t_d is not None:
            axp.plot(
                t_d, p_mpa,
                linewidth=1.2,
                linestyle=_swing_linestyle(bar),   # NEW
                color=swing_colors[bar],
                alpha=0.7
            )
    axp.set_title("Pressure schedule")
    axp.set_xlabel("Time (days)")
    axp.set_ylabel("Pressure (MPa)")
    axp.grid(True, alpha=0.3)

    handles, labels = axes[0].get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h
    fig.legend(uniq.values(), uniq.keys(), loc="upper center", bbox_to_anchor=(0.5, 0.96),
               ncol=min(6, len(uniq)), frameon=True, fontsize=9)

    fig.suptitle(f"Stress State — {cavern_label}", fontsize=13, fontweight='bold', y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.92])

    safe_key = cavern_key.replace(" ", "_")
    outpath = os.path.join(OUT_DIR, f"swing_stress_state_{safe_key}.png")
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight')
    print(f"  [SAVED] {outpath}")

    if SHOW:
        plt.show()
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: FACTOR OF SAFETY
# ══════════════════════════════════════════════════════════════════════════════

def plot_fig_fos(cavern_key, cavern_label, cases, swing_colors):
    """
    2x3 figure: 5 probe subplots showing FOS vs time for each swing,
    + 1 summary subplot showing minimum FOS across all probes.
    """
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    fos_by_swing = {}
    for c in cases:
        bar = c["swing_bar"]
        cp = c["case_path"]
        if not has_fos_files(cp):
            continue
        try:
            wall_points = load_wall_points(cp)
            probes = auto_generate_probes(wall_points)
            t_days, fos_probes = compute_FOS_per_probe(cp, probes)
            fos_by_swing[bar] = (t_days, fos_probes)
            min_val = min(np.min(v) for v in fos_probes.values())
            print(f"    [FOS] {bar} bar: min FOS = {min_val:.3f}")
        except Exception as e:
            print(f"    [WARN] FOS {bar} bar: {e}")

    if not fos_by_swing:
        plt.close(fig)
        print("  [SKIP] No FOS data available")
        return

    for i, ptype in enumerate(PROBE_ORDER):
        ax = axes[i]

        for bar in sorted(fos_by_swing.keys()):
            t_days, fos_probes = fos_by_swing[bar]
            if ptype not in fos_probes:
                continue
            fos = fos_probes[ptype]
            color = swing_colors[bar]
            ax.plot(
                t_days, fos,
                linewidth=1.8,
                linestyle=_swing_linestyle(bar),   # NEW
                color=color,
                label=f"{bar} bar/day"
            )

        ax.axhline(y=1.0, color='red', linestyle=':', linewidth=1.5, alpha=0.8)
        ax.set_title(f"FOS: {ptype}")
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Factor of Safety")
        ax.grid(True, alpha=0.3)

    ax_sum = axes[5]
    for bar in sorted(fos_by_swing.keys()):
        t_days, fos_probes = fos_by_swing[bar]
        all_fos = np.array([fos_probes[p] for p in PROBE_ORDER if p in fos_probes])
        min_fos = np.min(all_fos, axis=0)
        color = swing_colors[bar]
        ax_sum.plot(
            t_days, min_fos,
            linewidth=2,
            linestyle=_swing_linestyle(bar),   # NEW
            color=color,
            label=f"{bar} bar/day"
        )

    ax_sum.axhline(y=1.0, color='red', linestyle=':', linewidth=1.5, alpha=0.8)
    ax_sum.set_title("Min FOS (all probes)")
    ax_sum.set_xlabel("Time (days)")
    ax_sum.set_ylabel("Factor of Safety")
    ax_sum.grid(True, alpha=0.3)

    handles, labels = axes[0].get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h
    fig.legend(uniq.values(), uniq.keys(), loc="upper center", bbox_to_anchor=(0.5, 0.96),
               ncol=min(6, len(uniq)), frameon=True, fontsize=9)

    fig.suptitle(f"Factor of Safety — {cavern_label}", fontsize=13, fontweight='bold', y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.92])

    safe_key = cavern_key.replace(" ", "_")
    outpath = os.path.join(OUT_DIR, f"swing_fos_{safe_key}.png")
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight')
    print(f"  [SAVED] {outpath}")

    if SHOW:
        plt.show()
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    all_swing = discover_swing_cases(ROOT)

    if CAVERN is not None:
        cavern_low = CAVERN.lower()
        all_swing = [c for c in all_swing if c["cavern_key"].lower() == cavern_low]

    if SWING_BARS is not None:
        all_swing = [c for c in all_swing if c["swing_bar"] in SWING_BARS]

    if not all_swing:
        print(f"[ERROR] No swing cases found in {ROOT}")
        print("  Expected folders like: case_swing_10bar_regular1200/")
        if CAVERN is not None:
            print(f"  (filtered for CAVERN = '{CAVERN}')")
        return

    cavern_groups = {}
    for c in all_swing:
        key = c["cavern_key"]
        if key not in cavern_groups:
            cavern_groups[key] = []
        cavern_groups[key].append(c)

    all_bars = sorted(set(c["swing_bar"] for c in all_swing))
    swing_colors = build_swing_colors(all_bars)

    print("=" * 70)
    print("PRESSURE SWING COMPARISON")
    print("=" * 70)
    for ckey in sorted(cavern_groups):
        label = _nice_cavern_label(ckey)
        swings = [c["swing_bar"] for c in cavern_groups[ckey]]
        print(f"  {label:25s}  swings: {swings} bar/day")
    print(f"  Color map:  {all_bars} bar/day")
    print("=" * 70)

    for ckey in sorted(cavern_groups):
        cases = cavern_groups[ckey]
        cavern_label = cases[0]["cavern_label"]

        print(f"\n{'=' * 70}")
        print(f"  CAVERN: {cavern_label} ({ckey})")
        print(f"{'=' * 70}")

        print("\n  --- Figure 1: Convergence ---")
        plot_fig_convergence(ckey, cavern_label, cases, swing_colors)

        print("\n  --- Figure 2: Stress state ---")
        plot_fig_stress_state(ckey, cavern_label, cases, swing_colors)

        print("\n  --- Figure 3: Factor of Safety ---")
        plot_fig_fos(ckey, cavern_label, cases, swing_colors)

    print(f"\n{'=' * 70}")
    print("[DONE]")


if __name__ == "__main__":
    main()