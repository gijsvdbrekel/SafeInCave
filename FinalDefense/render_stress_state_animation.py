#!/usr/bin/env python3
"""
render_stress_state_animation.py — stress-state figures for the final-defense
slides.

Produces, for one or more simulations:

1. A STATIC p-q figure (stress_state_static_*.png) overlaying the full stress
   path of every case, with the start (open ring) and end (filled dot) of each
   path clearly marked.

2. An ANIMATED explainer (GIF always, MP4 when ffmpeg is present) with three
   synchronised panels over the operation timeline:
     [ 3D cavern + probe ]   [ p-q stress path ]   [ cavern pressure(t) ]
   The p-q path grows from its start point to the current point while a marker
   sweeps the pressure schedule — so the audience sees the stress state move as
   the pressure changes. Multiple cases are overlaid (one colour each).

Only the De Vries (2005) dilatancy boundary is drawn. Reads each case's
operation output (p_elems / q_elems), pressure schedule and wall geometry.
Off-screen / headless (PuTTY) safe.

Usage:
    python render_stress_state_animation.py
Edit CASES + the options below; works for any geometry / pressure scheme.
"""

import os
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
import json
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.animation as animation

import meshio
import pyvista as pv
import gmsh
import safeincave.PostProcessingTools as post

pv.OFF_SCREEN = True

# =============================================================================
# CONFIGURATION
# =============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_ROOT = os.path.normpath(os.path.join(
    SCRIPT_DIR, "..", "examples", "mechanics", "nobian", "Simulation", "output"))

# ---- Product 1: CASE COMPARISON (overlay several simulations at one probe) ----
# Cases to plot: (label, case folder [rel. to OUTPUT_ROOT or absolute], probe).
# List several to overlay them; the animation shows each cavern side by side
# (each with its own probe), and the p-q / pressure panels overlay all cases.
CASES = [
    ("Industry",  "case_leaching_linear_industry(21)_365days_SB_MD_regular1200", "mid"),
    # ("Transport", "case_leaching_linear_transport(180)_365days_SB_MD_regular1200", "mid"),
]

MAKE_STATIC = True           # static p-q figure with start/end markers
MAKE_ANIMATION = True        # animated explainer (caverns + p-q + pressure)

# ---- Product 2: MULTI-PROBE (one panel per probe, cases overlaid) --------------
# A grid with one p-q panel per probe. Every case in MULTIPROBE_CASES is overlaid
# in each probe panel (coloured per case), so you compare e.g. all probes of the
# fast-leached cavern vs. all probes of the string-failure cavern, etc. Like the
# plot_results.py stress-state grid, but cleaner. Static + animated.
MULTIPROBE_CASES = [
    ("Fast-leached",   "case_leaching_linear_industry(21)_365days_SB_MD_fastleached1200"),
    # ("String-failure", "case_..._tubefailure1200"),
    # ("Regular",        "case_leaching_linear_industry(21)_365days_SB_MD_regular1200"),
]
PROBES = ["top", "quarter", "mid", "threequarter", "bottom"]
MAKE_MULTIPROBE_STATIC = True
MAKE_MULTIPROBE_ANIM = True

OUT_DIR = os.path.join(SCRIPT_DIR, "stress_state_anim")

N_FRAMES_MAX = 180           # subsample the timeline to at most this many frames
FPS = 18
GIF_DPI = 110
MP4_DPI = 140

# Only the De Vries (2005) boundary is drawn.
SHOW_DILATANCY = ["devries_comp", "devries_ext"]

# Per-case colours (used when several cases are overlaid).
CASE_COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf"]

# 3D cavern thumbnail (left panel of the animation).
CAV_COLOR = "#cdd7e2"
CAV_CAM_AZIMUTH = 35.0
CAV_CAM_ELEVATION = 12.0

# Depth conversion for the cavern panel (1200k grids: model top at 908 m depth).
OVERBURDEN_TOP_DEPTH_M = 908.0
DOMAIN_HEIGHT_M = 660.0

PROBE_FRACTIONS = {"bottom": 0.0, "threequarter": 0.25, "mid": 0.5,
                   "quarter": 0.75, "top": 1.0}
PROBE_COLORS = {"top": "#e41a1c", "quarter": "#377eb8", "mid": "#4daf4a",
                "threequarter": "#984ea3", "bottom": "#ff7f00"}

MPA = 1.0e6
DAY_S = 86400.0
HOUR = 3600.0


# =============================================================================
# DATA LOADING
# =============================================================================

def _z_to_depth(z):
    return OVERBURDEN_TOP_DEPTH_M + (DOMAIN_HEIGHT_M - np.asarray(z, float))


def load_wall_points(case_dir):
    op = os.path.join(case_dir, "operation")
    m = meshio.read(os.path.join(op, "mesh", "geom.msh"))
    cells = m.cells_dict if hasattr(m, "cells_dict") else m.cells
    line = np.unique(np.asarray(cells["line"]).reshape(-1))
    upts, _, _ = post.read_node_vector(os.path.join(op, "u", "u.xdmf"))
    mapping = post.build_mapping(m.points, upts)
    wall_idx = np.array([mapping[i] for i in line], dtype=int)
    wp = upts[wall_idx]
    return wp[np.argsort(wp[:, 2])]


def probes_from_wall(wall_points):
    z = wall_points[:, 2]
    z_min, z_max = z.min(), z.max()
    probes = {}
    for name, frac in PROBE_FRACTIONS.items():
        zt = z_min + frac * (z_max - z_min)
        probes[name] = wall_points[int(np.argmin(np.abs(z - zt)))]
    return probes


def load_stress_path(case_dir, probe_xyz):
    op = os.path.join(case_dir, "operation")
    pts_p, t_p, p_el = post.read_cell_scalar(os.path.join(op, "p_elems", "p_elems.xdmf"))
    _, t_q, q_el = post.read_cell_scalar(os.path.join(op, "q_elems", "q_elems.xdmf"))
    n = min(len(t_p), len(t_q))
    idx = post.find_closest_point(probe_xyz, pts_p)
    t_days = np.asarray(t_p[:n], float) / DAY_S
    p = -p_el[:n, idx] / MPA      # compression-positive
    q = q_el[:n, idx] / MPA
    return t_days, p, q


def load_pressure(case_dir):
    d = json.load(open(os.path.join(case_dir, "pressure_schedule.json")))
    if "t_hours" in d and "p_MPa" in d:
        return np.asarray(d["t_hours"], float) / 24.0, np.asarray(d["p_MPa"], float)
    if "t_values_s" in d and "p_values_Pa" in d:
        return np.asarray(d["t_values_s"], float) / DAY_S, np.asarray(d["p_values_Pa"], float) / MPA
    raise RuntimeError("No usable pressure schedule fields")


# =============================================================================
# DILATANCY BOUNDARIES
# =============================================================================

def boundary_curves(p_lo, p_hi, npts=400):
    p = np.linspace(max(p_lo, 0.01), p_hi, npts)
    I1 = 3.0 * p
    out = {}

    def q_from(sq):
        return np.sqrt(3.0) * sq

    if "ratigan_027" in SHOW_DILATANCY:
        out["Ratigan 1991 (D=0.27)"] = (p, q_from(0.27 * I1),
                                        dict(color="#000000", ls="--", lw=1.4, alpha=0.8))
    if "spiers" in SHOW_DILATANCY:
        out["Spiers 1988"] = (p, q_from(0.27 * I1 + 1.9),
                              dict(color="#555555", ls="-.", lw=1.4, alpha=0.8))

    def devries(psi):
        sgn = np.sign(I1); sgn[sgn == 0] = 1.0
        denom = np.sqrt(3.0) * np.cos(psi) - 0.512 * np.sin(psi)
        return q_from(0.683 * ((I1 / sgn) ** 0.75) / denom + 1.5)

    if "devries_comp" in SHOW_DILATANCY:
        out["De Vries 2005 (comp)"] = (p, devries(np.pi / 6.0),
                                       dict(color="#888888", ls="-", lw=1.6, alpha=1.0))
    if "devries_ext" in SHOW_DILATANCY:
        out["De Vries 2005 (ext)"] = (p, devries(-np.pi / 6.0),
                                      dict(color="#000000", ls="--", lw=2.6, alpha=1.0))
    return out


# =============================================================================
# 3D CAVERN THUMBNAIL (left panel of the animation)
# =============================================================================

def load_cavern_surface(geom_msh, tag=29):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.open(geom_msh)
        ntags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.asarray(coords).reshape(-1, 3)
        index_of = {int(t): i for i, t in enumerate(ntags)}
        ents = gmsh.model.getEntitiesForPhysicalGroup(2, tag)
        blocks = []
        for e in ents:
            etypes, _et, enodes = gmsh.model.mesh.getElements(2, e)
            for et, en in zip(etypes, enodes):
                if et == 2:
                    blocks.append(np.asarray(en, dtype=np.int64).reshape(-1, 3))
        tris = np.vstack(blocks)
    finally:
        gmsh.finalize()
    ti = np.vectorize(index_of.get)(tris)
    faces = np.hstack([np.full((len(ti), 1), 3, dtype=np.int64), ti]).ravel()
    return pv.PolyData(coords, faces).clean()


def render_cavern_thumb(geom_msh, probe_xyz, out_png, pcolor):
    surf = load_cavern_surface(geom_msh)
    surf = surf.compute_normals(cell_normals=False, point_normals=True,
                                auto_orient_normals=True)
    b = surf.bounds
    ext = max(b[1]-b[0], b[3]-b[2], b[5]-b[4])
    cx, cy, cz = surf.center
    # Probe marker placed just proud of the wall, along the outward radial dir.
    radial = np.array([probe_xyz[0]-cx, probe_xyz[1]-cy, 0.0], float)
    n = np.linalg.norm(radial)
    radial = radial / n if n > 1e-9 else np.array([1.0, 0.0, 0.0])
    marker = np.array(probe_xyz, float) + radial * 0.02 * ext

    p = pv.Plotter(off_screen=True, window_size=(800, 1100))
    p.set_background("white")
    p.add_mesh(surf, color=CAV_COLOR, smooth_shading=True, specular=0.3,
               specular_power=15, ambient=0.30, diffuse=0.65)
    p.add_mesh(pv.Sphere(radius=0.05*ext, center=marker), color=pcolor,
               smooth_shading=True)
    p.enable_parallel_projection()
    p.camera.focal_point = (cx, cy, cz)
    p.camera_position = "yz"
    p.camera.azimuth = CAV_CAM_AZIMUTH
    p.camera.elevation = CAV_CAM_ELEVATION
    p.camera.parallel_scale = 0.60 * ext
    p.enable_anti_aliasing("ssaa")
    p.screenshot(out_png, transparent_background=True, scale=2)
    p.close()
    return out_png


# =============================================================================
# CASE LOADING
# =============================================================================

def load_case(label, case_rel, probe_name, color):
    case_dir = case_rel if os.path.isabs(case_rel) else os.path.join(OUTPUT_ROOT, case_rel)
    if not os.path.isdir(case_dir):
        raise SystemExit(f"[ERROR] case not found: {case_dir}")
    wall = load_wall_points(case_dir)
    probes = probes_from_wall(wall)
    probe_xyz = probes[probe_name]
    t_days, p_path, q_path = load_stress_path(case_dir, probe_xyz)
    t_pres, p_pres = load_pressure(case_dir)
    p_at_t = np.interp(t_days, t_pres, p_pres)
    return dict(label=label, case_dir=case_dir, probe=probe_name, color=color,
                probe_xyz=probe_xyz,
                t_days=t_days, p_path=p_path, q_path=q_path, p_at_t=p_at_t,
                geom=os.path.join(case_dir, "operation", "mesh", "geom.msh"))


def _pq_limits(cases):
    p_lo = max(0.0, min(c["p_path"].min() for c in cases) - 2.0)
    p_hi = max(c["p_path"].max() for c in cases) + 4.0
    q_hi = 1.30 * max(c["q_path"].max() for c in cases)
    return p_lo, p_hi, q_hi


# =============================================================================
# STATIC FIGURE  (full paths, start + end marked)
# =============================================================================

def plot_static(cases, out_path):
    fig, ax = plt.subplots(figsize=(9.5, 8.0))
    p_lo, p_hi, q_hi = _pq_limits(cases)
    for label, (pp, qq, st) in boundary_curves(p_lo, p_hi).items():
        ax.plot(pp, qq, label=label, **st)
    for c in cases:
        ax.plot(c["p_path"], c["q_path"], color=c["color"], lw=1.6, alpha=0.85,
                label=c["label"], zorder=3)
        ax.scatter([c["p_path"][0]], [c["q_path"][0]], s=150, facecolors="white",
                   edgecolors=c["color"], linewidths=2.4, zorder=5)
        ax.scatter([c["p_path"][-1]], [c["q_path"][-1]], s=150, color=c["color"],
                   edgecolors="black", linewidths=1.4, zorder=6)
    # Proxy handles explaining the start/end markers.
    ax.scatter([], [], s=150, facecolors="white", edgecolors="black",
               linewidths=2.0, label="start")
    ax.scatter([], [], s=150, facecolors="black", edgecolors="black",
               label="end")
    ax.set_xlim(p_lo, p_hi)
    ax.set_ylim(0.0, q_hi)
    ax.set_xlabel("Mean stress p (MPa)", fontsize=15)
    ax.set_ylabel("Differential stress q (MPa)", fontsize=15)
    ax.set_title("Stress paths over time", fontsize=16, fontweight="bold")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=10, loc="upper left", framealpha=0.92)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[SAVED] {out_path}")


# =============================================================================
# ANIMATION  (3 panels: cavern + p-q + pressure; multi-case overlay)
# =============================================================================

def build_animation(cases, cavern_pngs):
    n = len(cases)
    fig = plt.figure(figsize=(16 + 1.6 * n, 6.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.6 * n + 0.3, 1.25, 1.55],
                          wspace=0.24)
    # --- Cavern panel(s): one 3D thumbnail per case, side by side. ---
    cav_gs = gs[0].subgridspec(1, n, wspace=0.04)
    for i, c in enumerate(cases):
        axc = fig.add_subplot(cav_gs[0, i])
        axc.imshow(mpimg.imread(cavern_pngs[i]))
        axc.axis("off")
        axc.set_title(f"{c['label']}\n(probe: {c['probe']})", fontsize=12,
                      fontweight="bold", color=c["color"])
    ax_pq = fig.add_subplot(gs[1])
    ax_pr = fig.add_subplot(gs[2])

    # --- p-q panel ---
    p_lo, p_hi, q_hi = _pq_limits(cases)
    for label, (pp, qq, st) in boundary_curves(p_lo, p_hi).items():
        ax_pq.plot(pp, qq, label=label, **st)
    ax_pq.set_xlim(p_lo, p_hi)
    ax_pq.set_ylim(0.0, q_hi)
    ax_pq.set_xlabel("Mean stress p (MPa)", fontsize=13)
    ax_pq.set_ylabel("Differential stress q (MPa)", fontsize=13)
    ax_pq.set_title("Stress path", fontsize=15, fontweight="bold")
    ax_pq.grid(True, alpha=0.25)

    # --- pressure panel ---
    t_max = max(c["t_days"].max() for c in cases)
    pr_lo = min(c["p_at_t"].min() for c in cases) - 1.0
    pr_hi = max(c["p_at_t"].max() for c in cases) + 1.0
    ax_pr.set_xlim(0.0, t_max)
    ax_pr.set_ylim(pr_lo, pr_hi)
    ax_pr.set_xlabel("Time (days)", fontsize=13)
    ax_pr.set_ylabel("Cavern pressure (MPa)", fontsize=13)
    ax_pr.set_title("Cavern pressure", fontsize=15, fontweight="bold")
    ax_pr.grid(True, alpha=0.25)

    art = []
    for c in cases:
        col = c["color"]
        ax_pq.plot(c["p_path"][0], c["q_path"][0], "o", ms=11, mfc="white",
                   mec=col, mew=2.0, zorder=5)              # start ring
        (path_line,) = ax_pq.plot([], [], color=col, lw=2.0, zorder=4)
        (cur_pt,) = ax_pq.plot([], [], "o", color=col, ms=11, mec="black",
                               mew=1.1, zorder=6)
        ax_pr.plot(c["t_days"], c["p_at_t"], color=col, lw=1.0, alpha=0.30, zorder=1)
        (pr_line,) = ax_pr.plot([], [], color=col, lw=2.0, zorder=2,
                                label=c["label"])
        (pr_pt,) = ax_pr.plot([], [], "o", color=col, ms=9, mec="white",
                              mew=1.0, zorder=3)
        art.append((c, path_line, cur_pt, pr_line, pr_pt))

    ax_pq.legend(fontsize=9, loc="upper left", framealpha=0.9)
    ax_pr.legend(fontsize=10, loc="upper right", framealpha=0.9)
    time_txt = fig.text(0.5, 0.975, "", ha="center", va="top", fontsize=15,
                        fontweight="bold")

    n_frames = N_FRAMES_MAX

    def update(fi):
        prog = fi / (n_frames - 1)
        changed = [time_txt]
        for (c, path_line, cur_pt, pr_line, pr_pt) in art:
            k = int(round(prog * (len(c["t_days"]) - 1)))
            path_line.set_data(c["p_path"][:k+1], c["q_path"][:k+1])
            cur_pt.set_data([c["p_path"][k]], [c["q_path"][k]])
            pr_line.set_data(c["t_days"][:k+1], c["p_at_t"][:k+1])
            pr_pt.set_data([c["t_days"][k]], [c["p_at_t"][k]])
            changed += [path_line, cur_pt, pr_line, pr_pt]
        k0 = int(round(prog * (len(cases[0]["t_days"]) - 1)))
        c0 = cases[0]
        time_txt.set_text(f"t = {c0['t_days'][k0]:6.1f} days      "
                          f"cavern pressure = {c0['p_at_t'][k0]:4.1f} MPa")
        return changed

    anim = animation.FuncAnimation(fig, update, frames=n_frames,
                                   interval=1000 / FPS, blit=False)
    return fig, anim


# =============================================================================
# MULTI-PROBE  (one p-q panel per probe; cases overlaid in each panel)
# =============================================================================

def load_mp_case(label, case_rel, color, probe_names):
    case_dir = case_rel if os.path.isabs(case_rel) else os.path.join(OUTPUT_ROOT, case_rel)
    if not os.path.isdir(case_dir):
        raise SystemExit(f"[ERROR] case not found: {case_dir}")
    op = os.path.join(case_dir, "operation")
    probes = probes_from_wall(load_wall_points(case_dir))
    pts_p, t_p, p_el = post.read_cell_scalar(os.path.join(op, "p_elems", "p_elems.xdmf"))
    _, t_q, q_el = post.read_cell_scalar(os.path.join(op, "q_elems", "q_elems.xdmf"))
    n = min(len(t_p), len(t_q))
    t_days = np.asarray(t_p[:n], float) / DAY_S
    paths = {}
    for name in probe_names:
        idx = post.find_closest_point(probes[name], pts_p)
        paths[name] = (-p_el[:n, idx] / MPA, q_el[:n, idx] / MPA)
    t_pres, p_pres = load_pressure(case_dir)
    return dict(label=label, color=color, paths=paths, t_days=t_days,
                p_at_t=np.interp(t_days, t_pres, p_pres))


def _mp_limits(mps, probe_names):
    allp = [v[0] for mp in mps for v in mp["paths"].values()]
    allq = [v[1] for mp in mps for v in mp["paths"].values()]
    p_lo = max(0.0, min(p.min() for p in allp) - 2.0)
    p_hi = max(p.max() for p in allp) + 4.0
    q_hi = 1.30 * max(q.max() for q in allq)
    return p_lo, p_hi, q_hi


def _mp_grid(mps, probe_names):
    """Figure + (probe-name -> axis) for a 2x3 grid: 5 probes + pressure."""
    p_lo, p_hi, q_hi = _mp_limits(mps, probe_names)
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    probe_axes = {}
    for ax, name in zip(axes, probe_names):
        for _lab, (pp, qq, st) in boundary_curves(p_lo, p_hi).items():
            ax.plot(pp, qq, **st)
        ax.set_xlim(p_lo, p_hi)
        ax.set_ylim(0.0, q_hi)
        ax.set_title(name.capitalize(), fontsize=15, fontweight="bold", color="#333")
        ax.set_xlabel("Mean stress p (MPa)", fontsize=12)
        ax.set_ylabel("Differential stress q (MPa)", fontsize=12)
        ax.grid(True, alpha=0.25)
        probe_axes[name] = ax
    ax_pr = axes[len(probe_names)] if len(probe_names) < 6 else None
    # Title on top, then one legend entry per case below it.
    fig.suptitle("Stress paths per probe", fontsize=18, fontweight="bold", y=0.99)
    handles = [plt.Line2D([], [], color=mp["color"], lw=3, label=mp["label"])
               for mp in mps]
    fig.legend(handles=handles, loc="upper center", ncol=min(len(mps), 5),
               fontsize=12, frameon=True, bbox_to_anchor=(0.5, 0.955))
    return fig, axes, probe_axes, ax_pr


def plot_static_multiprobe(mps, probe_names, out_path):
    fig, axes, probe_axes, ax_pr = _mp_grid(mps, probe_names)
    for name, ax in probe_axes.items():
        for mp in mps:
            pp, qq = mp["paths"][name]
            col = mp["color"]
            ax.plot(pp, qq, color=col, lw=1.5, alpha=0.85, zorder=3)
            ax.scatter([pp[0]], [qq[0]], s=110, facecolors="white", edgecolors=col,
                       linewidths=2.2, zorder=5)
            ax.scatter([pp[-1]], [qq[-1]], s=110, color=col, edgecolors="black",
                       linewidths=1.3, zorder=6)
    if ax_pr is not None:
        for mp in mps:
            ax_pr.plot(mp["t_days"], mp["p_at_t"], color=mp["color"], lw=1.4)
        ax_pr.set_title("Cavern pressure", fontsize=15, fontweight="bold")
        ax_pr.set_xlabel("Time (days)", fontsize=12)
        ax_pr.set_ylabel("Cavern pressure (MPa)", fontsize=12)
        ax_pr.grid(True, alpha=0.25)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"[SAVED] {out_path}  (start=ring, end=dot)")


def build_animation_multiprobe(mps, probe_names):
    fig, axes, probe_axes, ax_pr = _mp_grid(mps, probe_names)
    arts = []   # (mp, {name: (line, pt)})
    for mp in mps:
        col = mp["color"]
        per_probe = {}
        for name, ax in probe_axes.items():
            pp, qq = mp["paths"][name]
            ax.plot(pp[0], qq[0], "o", ms=9, mfc="white", mec=col, mew=2.0, zorder=5)
            (line,) = ax.plot([], [], color=col, lw=1.8, zorder=4)
            (pt,) = ax.plot([], [], "o", color=col, ms=10, mec="black", mew=1.0, zorder=6)
            per_probe[name] = (line, pt)
        arts.append((mp, per_probe))

    pr_arts = []
    if ax_pr is not None:
        for mp in mps:
            ax_pr.plot(mp["t_days"], mp["p_at_t"], color=mp["color"], lw=1.0, alpha=0.30)
            (prl,) = ax_pr.plot([], [], color=mp["color"], lw=2.0)
            (prp,) = ax_pr.plot([], [], "o", color=mp["color"], ms=9, mec="white", mew=1.0)
            pr_arts.append((mp, prl, prp))
        t_max = max(mp["t_days"].max() for mp in mps)
        ax_pr.set_xlim(0.0, t_max)
        ax_pr.set_ylim(min(mp["p_at_t"].min() for mp in mps) - 1,
                       max(mp["p_at_t"].max() for mp in mps) + 1)
        ax_pr.set_title("Cavern pressure", fontsize=15, fontweight="bold")
        ax_pr.set_xlabel("Time (days)", fontsize=12)
        ax_pr.set_ylabel("Cavern pressure (MPa)", fontsize=12)
        ax_pr.grid(True, alpha=0.25)

    time_txt = fig.text(0.5, 0.915, "", ha="center", va="top", fontsize=14,
                        fontweight="bold", color="#222")
    n_frames = N_FRAMES_MAX

    def update(fi):
        prog = fi / (n_frames - 1)
        changed = [time_txt]
        for mp, per_probe in arts:
            k = int(round(prog * (len(mp["t_days"]) - 1)))
            for name, (line, pt) in per_probe.items():
                pp, qq = mp["paths"][name]
                line.set_data(pp[:k+1], qq[:k+1])
                pt.set_data([pp[k]], [qq[k]])
                changed += [line, pt]
        for mp, prl, prp in pr_arts:
            k = int(round(prog * (len(mp["t_days"]) - 1)))
            prl.set_data(mp["t_days"][:k+1], mp["p_at_t"][:k+1])
            prp.set_data([mp["t_days"][k]], [mp["p_at_t"][k]])
            changed += [prl, prp]
        k0 = int(round(prog * (len(mps[0]["t_days"]) - 1)))
        time_txt.set_text(f"t = {mps[0]['t_days'][k0]:6.1f} days")
        return changed

    fig.tight_layout(rect=[0, 0, 1, 0.88])
    anim = animation.FuncAnimation(fig, update, frames=n_frames,
                                   interval=1000 / FPS, blit=False)
    return fig, anim


def _save_anim(fig, anim, stem):
    gif = os.path.join(OUT_DIR, stem + ".gif")
    anim.save(gif, writer=animation.PillowWriter(fps=FPS), dpi=GIF_DPI)
    print(f"[SAVED] {gif}")
    if animation.writers.is_available("ffmpeg"):
        mp4 = os.path.join(OUT_DIR, stem + ".mp4")
        anim.save(mp4, writer=animation.FFMpegWriter(fps=FPS, bitrate=4000), dpi=MP4_DPI)
        print(f"[SAVED] {mp4}")
    else:
        print("[INFO] ffmpeg not available - MP4 skipped (GIF written).")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---- Product 1: case comparison (overlay) ----
    if MAKE_STATIC or MAKE_ANIMATION:
        cases = [load_case(lbl, rel, probe, CASE_COLORS[i % len(CASE_COLORS)])
                 for i, (lbl, rel, probe) in enumerate(CASES)]
        for c in cases:
            print(f"[CASE] {c['label']:14s} probe={c['probe']}  nt={len(c['t_days'])}")
        tag = "_".join(c["label"].lower().replace(" ", "") for c in cases)

        if MAKE_STATIC:
            plot_static(cases, os.path.join(OUT_DIR, f"stress_state_static_{tag}.png"))

        if MAKE_ANIMATION:
            thumbs = []
            for c in cases:
                t = os.path.join(OUT_DIR, f"_cavern_{c['label'].lower().replace(' ', '')}.png")
                render_cavern_thumb(c["geom"], c["probe_xyz"], t, c["color"])
                thumbs.append(t)
            fig, anim = build_animation(cases, thumbs)
            _save_anim(fig, anim, f"stress_state_anim_{tag}")

    # ---- Product 2: multi-probe (one panel per probe; cases overlaid) ----
    if MAKE_MULTIPROBE_STATIC or MAKE_MULTIPROBE_ANIM:
        mps = [load_mp_case(lbl, rel, CASE_COLORS[i % len(CASE_COLORS)], PROBES)
               for i, (lbl, rel) in enumerate(MULTIPROBE_CASES)]
        for mp in mps:
            print(f"[MULTIPROBE] {mp['label']:16s} probes={PROBES}  nt={len(mp['t_days'])}")
        mtag = "_".join(mp["label"].lower().replace(" ", "").replace("-", "") for mp in mps)
        if MAKE_MULTIPROBE_STATIC:
            plot_static_multiprobe(mps, PROBES,
                                   os.path.join(OUT_DIR, f"stress_multiprobe_static_{mtag}.png"))
        if MAKE_MULTIPROBE_ANIM:
            fig, anim = build_animation_multiprobe(mps, PROBES)
            _save_anim(fig, anim, f"stress_multiprobe_anim_{mtag}")

    print("[DONE]")


if __name__ == "__main__":
    main()
