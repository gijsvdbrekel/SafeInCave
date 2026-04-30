"""
plot_sensitivity.py — Sensitivity-analysis figures for MSc thesis.

Companion to plot_results.py. Consumes `*_sens_<tag>` case folders produced by
`Run_sensitivity.py` and renders per-dimension overlay figures + aggregate
profiling summaries.

Per-dimension overlays (one call per dimension):
    plot_convergence_overlay    — ΔV/V₀ (%) vs time, pressure on twin axis
    plot_convergence_rate       — dΔV/dt, smoothed
    plot_strain_rate_probes     — |ε̇_creep| at roof / mid / floor
    plot_fos_overlay            — global-min FOS vs time
    plot_mc_failure_overlay     — MC-failed interlayer cell count (skips gracefully)

Aggregate profiling (one call total):
    plot_profiling_bars         — stacked h-bars: ct / ksp / assemble / other
    plot_newton_histogram       — Newton-iters distribution per case
    plot_walltime_vs_newton     — wall time vs total Newton iters

Usage:
    python plot_sensitivity.py                # render all dimensions
    python plot_sensitivity.py --self-test    # synthetic-data smoke test
"""

import os
import sys
import json
import glob
import argparse
import warnings

import numpy as np
import matplotlib.pyplot as plt

import safeincave.PostProcessingTools as post

# Reuse loaders/helpers from plot_results. Import is side-effect-bearing
# (matplotlib rcParams update + constants); that is acceptable here.
from plot_results import (
    load_p_q,
    load_sig,
    read_pressure_schedule,
    load_wall_points,
    compute_convergence_3d_percent,
    _get_interlayer_cell_mask,
    path_u_xdmf,
    path_sig_xdmf,
    path_p_xdmf,
    path_geom_msh,
    DAY,
    HOUR,
    MPA,
)

# =============================================================================
# STYLE & PATHS
# =============================================================================

plt.rcParams.update({
    "font.size":        18,
    "axes.titlesize":   22,
    "axes.labelsize":   20,
    "xtick.labelsize":  18,
    "ytick.labelsize":  18,
    "legend.fontsize":  16,
    "figure.titlesize": 24,
    "lines.linewidth":  2.0,
})

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SIM_ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", "Simulation", "output"))
FIG_ROOT = os.path.join(_SCRIPT_DIR, "_sensitivity_figures")

# =============================================================================
# SENSITIVITY DIMENSION DEFINITIONS
# =============================================================================
# Each dimension overlays its listed cases. `md_B_baseline` reappears as the
# control series in every dimension except creep_toggles (where it IS baseline).

DIMENSIONS = {
    "creep_toggles": {
        "title": "MD creep-mechanism toggles",
        "cases": ["md_B_baseline", "md_B_no_transient", "md_B_no_reverse", "md_B_steady_state"],
        "labels": {
            "md_B_baseline":     "MD full (baseline)",
            "md_B_no_transient": "No transient (α_w=β_w=0)",
            "md_B_no_reverse":   "No reverse (δ=0)",
            "md_B_steady_state": "Steady-state only",
        },
    },
    "constitutive": {
        "title": "Constitutive model: MD vs TUD2023, scenario A vs B",
        "cases": ["md_B_baseline", "sic_B_baseline", "md_A_baseline"],
        "labels": {
            "md_B_baseline":  "MD_B",
            "sic_B_baseline": "TUD2023_B",
            "md_A_baseline":  "MD_A",
        },
    },
    "startup": {
        "title": "Startup: leaching vs equilibrium",
        "cases": ["md_B_baseline", "no_leaching"],
        "labels": {
            "md_B_baseline": "Leaching phase (baseline)",
            "no_leaching":   "Equilibrium startup",
        },
    },
    "geometry": {
        "title": "Cavern geometry & volume",
        "cases": ["md_B_baseline", "tilted_cavern", "volume_600"],
        "labels": {
            "md_B_baseline":  "Regular 1,200,000 m³ (baseline)",
            "tilted_cavern":  "Tilted 1,200,000 m³",
            "volume_600":     "Regular 600,000 m³",
        },
    },
    "secondary_creep": {
        "title": "Secondary creep mechanisms",
        "cases": ["md_B_baseline", "sic_B_baseline", "no_pressure_solution", "no_kelvin", "no_desai", "no_kelvin_no_desai"],
        "labels": {
            "md_B_baseline":        "Baseline (MD)",
            "sic_B_baseline":       "Baseline (TUD2023)",
            "no_pressure_solution": "No pressure-solution",
            "no_kelvin":            "No Kelvin (TUD2023)",
            "no_desai":             "No Desai (TUD2023)",
            "no_kelvin_no_desai":   "No Kelvin + No Desai (TUD2023)",
        },
    },
    "timestep": {
        "title": "Numerical time-step",
        "cases": ["md_B_baseline", "dt_half", "dt_double"],
        "labels": {
            "md_B_baseline": "dt = 2 h (baseline)",
            "dt_half":       "dt = 1 h",
            "dt_double":     "dt = 4 h",
        },
    },
    "amplitude": {
        "title": "Pressure amplitude scaling",
        "cases": ["md_B_baseline", "amplitude_half", "amplitude_double"],
        "labels": {
            "md_B_baseline":     "1.0× amplitude (baseline)",
            "amplitude_half":    "0.5× amplitude",
            "amplitude_double":  "2.0× amplitude",
        },
    },
}

# Discrete colour cycle — consistent across panels of the same dimension.
_COLOR_CYCLE = plt.get_cmap("tab10").colors


def _case_style(idx):
    """Return (color, linestyle) for the idx-th case in a dimension."""
    ls_cycle = ["-", "--", "-.", ":"]
    return _COLOR_CYCLE[idx % len(_COLOR_CYCLE)], ls_cycle[idx % len(ls_cycle)]


# =============================================================================
# DISCOVERY & LOADING
# =============================================================================

def discover_sens_cases(root=SIM_ROOT):
    """
    Scan `root` for sensitivity-run folders and return a dict keyed by case
    name (the `_sens_<tag>` suffix).

    A case is any folder matching `*_sens_*`. Folders without `profiling.json`
    are still reported (they may be pre-existing baseline runs we want to
    plot alongside); profiling panels handle the absence gracefully.
    """
    if not os.path.isdir(root):
        warnings.warn(f"Root directory not found: {root}")
        return {}

    records = {}
    for folder in sorted(glob.glob(os.path.join(root, "*_sens_*"))):
        if not os.path.isdir(folder):
            continue
        base = os.path.basename(folder)
        # tag = everything after the last '_sens_'
        idx = base.rfind("_sens_")
        if idx < 0:
            continue
        case_name = base[idx + len("_sens_"):]
        records[case_name] = {
            "case_name": case_name,
            "case_path": folder,
            "has_profiling": os.path.isfile(os.path.join(folder, "profiling.json")),
        }
    return records


def load_profiling(folder):
    """Return the `profiling.json` payload for a case folder, or None if missing/broken."""
    p = os.path.join(folder, "profiling.json")
    if not os.path.isfile(p):
        return None
    try:
        with open(p, "r") as f:
            return json.load(f)
    except Exception as e:
        warnings.warn(f"Failed to parse {p}: {e}")
        return None


# =============================================================================
# DERIVED QUANTITIES
# =============================================================================

def compute_convergence_rate(t_days, conv_pct, window=5):
    """
    Derivative of convergence curve (d(ΔV/V₀)/dt) via np.gradient, lightly smoothed
    with a centred moving average of length `window` (odd). Returns rate in %/day.
    """
    if len(t_days) < 2:
        return np.zeros_like(t_days)
    rate = np.gradient(conv_pct, t_days)
    if window and window > 1 and len(rate) >= window:
        w = window if window % 2 == 1 else window + 1
        kernel = np.ones(w) / w
        rate = np.convolve(rate, kernel, mode="same")
    return rate


def compute_strain_rate_timeseries(case_folder, probe_xyz, normalize_per_sec=True):
    """
    Compute the norm of the displacement-rate vector at a probe point, used as
    a proxy for local creep strain-rate. Returns (t_days, speed) where speed is
    in 1/s if normalize_per_sec else dimensionless displacement step norm.

    Uses the displacement field (`u.xdmf`) which already reflects creep.
    Finite-differences ||du/dt|| at the closest node to probe_xyz.
    """
    from scipy.spatial import cKDTree

    u_xdmf = path_u_xdmf(case_folder)
    if not os.path.isfile(u_xdmf):
        return None, None

    points, time_list, u_field = post.read_node_vector(u_xdmf)
    t = np.asarray(time_list, float)
    if len(t) < 2:
        return None, None

    tree = cKDTree(points)
    _, idx = tree.query(np.asarray(probe_xyz, float))

    u = u_field[:, idx, :]  # (Nt, 3)
    du = np.diff(u, axis=0)
    dt = np.diff(t)
    dt = np.where(dt > 0, dt, np.nan)

    speed = np.linalg.norm(du, axis=1) / (dt if normalize_per_sec else 1.0)
    # Forward-fill onto original time axis (pad first sample with NaN)
    speed_full = np.concatenate([[np.nan], speed])
    return t / DAY, speed_full


# =============================================================================
# PER-DIMENSION OVERLAY PLOTS
# =============================================================================

def _resolve_dim_cases(dim_key, discovered):
    spec = DIMENSIONS[dim_key]
    present = []
    missing = []
    for name in spec["cases"]:
        if name in discovered:
            present.append((name, discovered[name]))
        else:
            missing.append(name)
    if missing:
        warnings.warn(f"[{dim_key}] missing cases, will be skipped: {missing}")
    return spec, present


def _fig_path(dim_key, stem):
    out_dir = os.path.join(FIG_ROOT, dim_key)
    os.makedirs(out_dir, exist_ok=True)
    return os.path.join(out_dir, stem)


def _save(fig, dim_key, stem):
    for ext in ("pdf", "png"):
        p = _fig_path(dim_key, f"{stem}.{ext}")
        fig.savefig(p, bbox_inches="tight", dpi=200)
    plt.close(fig)


def plot_convergence_overlay(dim_key, discovered):
    spec, cases = _resolve_dim_cases(dim_key, discovered)
    if not cases:
        return
    fig, ax = plt.subplots(figsize=(6.3, 4.0))
    ax2 = ax.twinx()
    pressure_drawn = False

    for i, (name, rec) in enumerate(cases):
        try:
            t_d, conv = compute_convergence_3d_percent(rec["case_path"])
        except Exception as e:
            warnings.warn(f"[{dim_key}/{name}] convergence failed: {e}")
            continue
        color, ls = _case_style(i)
        label = spec["labels"].get(name, name)
        ax.plot(t_d, conv, color=color, linestyle=ls, label=label)

        if not pressure_drawn:
            t_h, p = read_pressure_schedule(rec["case_path"])
            if t_h is not None:
                ax2.plot(t_h / 24.0, p, color="0.65", linewidth=0.8, alpha=0.7)
                ax2.set_ylabel("Pressure (MPa)", color="0.4")
                ax2.tick_params(axis="y", colors="0.4")
                pressure_drawn = True

    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Volume convergence ΔV/V₀ (%)")
    ax.set_title(spec["title"])
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", framealpha=0.9)
    _save(fig, dim_key, "convergence_overlay")


def plot_convergence_rate(dim_key, discovered):
    spec, cases = _resolve_dim_cases(dim_key, discovered)
    if not cases:
        return
    fig, ax = plt.subplots(figsize=(6.3, 4.0))
    for i, (name, rec) in enumerate(cases):
        try:
            t_d, conv = compute_convergence_3d_percent(rec["case_path"])
        except Exception as e:
            warnings.warn(f"[{dim_key}/{name}] convergence failed: {e}")
            continue
        rate = compute_convergence_rate(t_d, conv, window=11)
        color, ls = _case_style(i)
        label = spec["labels"].get(name, name)
        ax.plot(t_d, rate, color=color, linestyle=ls, label=label)

    ax.set_xlabel("Time (days)")
    ax.set_ylabel("d(ΔV/V₀)/dt (%/day)")
    ax.set_title(f"{spec['title']} — convergence rate")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", framealpha=0.9)
    _save(fig, dim_key, "convergence_rate")


def plot_strain_rate_probes(dim_key, discovered):
    """Three-panel figure: roof / mid / floor displacement-rate norms."""
    spec, cases = _resolve_dim_cases(dim_key, discovered)
    if not cases:
        return

    probe_names = ["roof", "mid", "floor"]

    fig, axes = plt.subplots(3, 1, figsize=(6.3, 8.0), sharex=True)

    for i, (name, rec) in enumerate(cases):
        try:
            wall = load_wall_points(rec["case_path"])
        except Exception as e:
            warnings.warn(f"[{dim_key}/{name}] wall points failed: {e}")
            continue

        z = wall[:, 2]
        z_min, z_max = float(z.min()), float(z.max())
        probes = {
            "floor": wall[int(np.argmin(np.abs(z - z_min)))],
            "mid":   wall[int(np.argmin(np.abs(z - 0.5 * (z_min + z_max))))],
            "roof":  wall[int(np.argmin(np.abs(z - z_max)))],
        }
        color, ls = _case_style(i)
        label = spec["labels"].get(name, name)

        for ax, pname in zip(axes, probe_names):
            t_d, speed = compute_strain_rate_timeseries(rec["case_path"], probes[pname])
            if t_d is None:
                continue
            ax.plot(t_d, speed, color=color, linestyle=ls, label=label if pname == "roof" else None)

    for ax, pname in zip(axes, probe_names):
        ax.set_ylabel(f"‖u̇‖ at {pname} (m/s)")
        ax.grid(True, alpha=0.3)
        ax.set_yscale("log")

    axes[-1].set_xlabel("Time (days)")
    axes[0].set_title(f"{spec['title']} — wall displacement rate")
    axes[0].legend(loc="best", framealpha=0.9)
    _save(fig, dim_key, "strain_rate_probes")


def _compute_fos_timeseries(case_folder):
    """Minimal FOS-min-global computation for the sensitivity overlay."""
    try:
        t_pq, p_elems, q_elems = load_p_q(case_folder)
    except Exception as e:
        warnings.warn(f"load_p_q failed for {case_folder}: {e}")
        return None, None
    try:
        t_sig, sig_vals, _ = load_sig(case_folder)
    except Exception as e:
        warnings.warn(f"load_sig failed for {case_folder}: {e}")
        return None, None

    n = min(len(t_pq), len(t_sig))
    t = t_pq[:n] / DAY
    p = p_elems[:n]
    q = q_elems[:n]

    # Simple Mohr-Coulomb-flavoured FOS proxy: (c·cosφ + p·sinφ) / (q/√3)
    # Using salt bulk params as coarse reference — thesis uses full field FOS,
    # which we approximate here for overlay legibility.
    c = 4.0e6    # Pa
    phi = np.deg2rad(35.0)
    num = c * np.cos(phi) + np.maximum(p, 0.0) * np.sin(phi)
    den = np.maximum(q, 1.0) / np.sqrt(3.0)
    FOS = num / den
    fos_min = np.nanmin(np.where(np.isfinite(FOS), FOS, np.nan), axis=1)
    return t, fos_min


def plot_fos_overlay(dim_key, discovered):
    spec, cases = _resolve_dim_cases(dim_key, discovered)
    if not cases:
        return
    fig, ax = plt.subplots(figsize=(6.3, 4.0))
    for i, (name, rec) in enumerate(cases):
        t, fos = _compute_fos_timeseries(rec["case_path"])
        if t is None:
            continue
        color, ls = _case_style(i)
        label = spec["labels"].get(name, name)
        ax.plot(t, fos, color=color, linestyle=ls, label=label)

    ax.axhline(1.0, color="k", linewidth=0.8, linestyle=":")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Global-min FOS")
    ax.set_title(f"{spec['title']} — factor of safety")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", framealpha=0.9)
    _save(fig, dim_key, "fos_overlay")


def plot_mc_failure_overlay(dim_key, discovered):
    """
    MC-failed interlayer cell count over time. Only meaningful for heterogeneous
    meshes (those with Interlayer_1 / Interlayer_2 physical groups). Cases
    without interlayer cells are skipped with a note.
    """
    spec, cases = _resolve_dim_cases(dim_key, discovered)
    if not cases:
        return

    fig, ax = plt.subplots(figsize=(6.3, 4.0))
    c_MPa, phi_deg = 4.0, 35.0
    c_Pa = c_MPa * 1e6
    phi = np.deg2rad(phi_deg)
    any_drawn = False

    for i, (name, rec) in enumerate(cases):
        try:
            t_sig, sig_vals, _ = load_sig(rec["case_path"])
        except Exception as e:
            warnings.warn(f"[{dim_key}/{name}] load_sig failed: {e}")
            continue
        n_cells = sig_vals.shape[1]
        mask, _ = _get_interlayer_cell_mask({"case_path": rec["case_path"],
                                             "case_name": name}, n_cells)
        if mask is None or not mask.any():
            continue

        sig = sig_vals[:, mask]  # (Nt, nil, 3, 3)
        # Principal stresses (compression-positive)
        sig_eff = -0.5 * (sig + np.swapaxes(sig, -1, -2))
        eigs = np.linalg.eigvalsh(sig_eff)  # ascending → (Nt, nil, 3)
        s3 = eigs[..., 0]
        s1 = eigs[..., 2]
        # MC: failed when (s1 - s3)/2 >= c·cosφ + (s1+s3)/2·sinφ
        tau_max = 0.5 * (s1 - s3)
        p_mean = 0.5 * (s1 + s3)
        yield_lhs = tau_max - (c_Pa * np.cos(phi) + np.maximum(p_mean, 0.0) * np.sin(phi))
        failed = (yield_lhs >= 0).sum(axis=1)

        color, ls = _case_style(i)
        label = spec["labels"].get(name, name)
        ax.plot(np.asarray(t_sig, float) / DAY, failed,
                color=color, linestyle=ls, label=label)
        any_drawn = True

    if not any_drawn:
        plt.close(fig)
        print(f"  [MC] {dim_key}: no cases with interlayer cells — skipping")
        return

    ax.set_xlabel("Time (days)")
    ax.set_ylabel(f"Failed interlayer cells (c={c_MPa} MPa, φ={phi_deg}°)")
    ax.set_title(f"{spec['title']} — MC failure count")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", framealpha=0.9)
    _save(fig, dim_key, "mc_failure_overlay")


# =============================================================================
# AGGREGATE PROFILING PLOTS
# =============================================================================

def plot_profiling_bars(discovered):
    """Stacked horizontal bars: ct / ksp / assemble / other seconds per case."""
    rows = []
    for name, rec in discovered.items():
        payload = load_profiling(rec["case_path"])
        if payload is None:
            continue
        s = payload.get("summary", {})
        wall = float(payload.get("wall_time_s", 0.0))
        ct = float(s.get("ct_s_total", 0.0))
        ksp = float(s.get("ksp_s_total", 0.0))
        asmbl = float(s.get("assemble_s_total", 0.0))
        other = max(0.0, wall - ct - ksp - asmbl)
        rows.append((name, ct, asmbl, ksp, other, wall))

    if not rows:
        print("[profiling] no profiling.json found — skipping profiling_bars")
        return

    rows.sort(key=lambda r: r[-1])
    names = [r[0] for r in rows]
    ct = np.array([r[1] for r in rows])
    asmbl = np.array([r[2] for r in rows])
    ksp = np.array([r[3] for r in rows])
    other = np.array([r[4] for r in rows])

    fig, ax = plt.subplots(figsize=(9.0, max(3.0, 0.5 * len(rows) + 1.5)))
    y = np.arange(len(rows))
    left = np.zeros_like(ct)
    for seg, lbl, color in [
        (ct,    "compute_CT", "#4c72b0"),
        (asmbl, "assemble",   "#dd8452"),
        (ksp,   "KSP solve",  "#55a868"),
        (other, "other",      "#8172b3"),
    ]:
        ax.barh(y, seg, left=left, label=lbl, color=color, edgecolor="white")
        left = left + seg

    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.set_xlabel("Wall time (s)")
    ax.set_title("CPU profile per sensitivity case")
    ax.legend(loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)

    os.makedirs(FIG_ROOT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_ROOT, f"profiling_bars.{ext}"), bbox_inches="tight", dpi=200)
    plt.close(fig)


def plot_newton_histogram(discovered):
    """Small-multiples histogram of per-step Newton iteration counts."""
    cases = []
    for name, rec in discovered.items():
        payload = load_profiling(rec["case_path"])
        if payload is None:
            continue
        per_step = payload.get("per_step")
        if not per_step:
            continue
        niters = [int(s.get("newton_iters", 0)) for s in per_step if s.get("newton_iters")]
        if not niters:
            continue
        cases.append((name, np.asarray(niters)))

    if not cases:
        print("[profiling] no per-step data — skipping newton_histogram")
        return

    n = len(cases)
    ncols = min(3, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 2.5 * nrows), squeeze=False)

    for ax, (name, niters) in zip(axes.flat, cases):
        bins = np.arange(niters.min(), niters.max() + 2) - 0.5
        ax.hist(niters, bins=bins, color="#4c72b0", edgecolor="white")
        ax.set_title(name, fontsize=11)
        ax.set_xlabel("Newton iters")
        ax.set_ylabel("Steps")
        ax.grid(True, alpha=0.3)

    for ax in axes.flat[len(cases):]:
        ax.axis("off")

    fig.suptitle("Newton iterations per step")
    fig.tight_layout()

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_ROOT, f"newton_histogram.{ext}"), bbox_inches="tight", dpi=200)
    plt.close(fig)


def plot_walltime_vs_newton(discovered):
    """Scatter of wall time vs total Newton iterations — one point per case."""
    pts = []
    for name, rec in discovered.items():
        payload = load_profiling(rec["case_path"])
        if payload is None:
            continue
        wall = float(payload.get("wall_time_s", 0.0))
        total_n = int(payload.get("n_newton_total", 0))
        if total_n <= 0 or wall <= 0:
            continue
        pts.append((name, wall, total_n))

    if not pts:
        print("[profiling] no wall-time / Newton data — skipping walltime_vs_newton")
        return

    fig, ax = plt.subplots(figsize=(6.3, 5.0))
    for name, wall, total_n in pts:
        ax.scatter(total_n, wall, s=60)
        ax.annotate(name, (total_n, wall), xytext=(4, 4),
                    textcoords="offset points", fontsize=9)

    ax.set_xlabel("Total Newton iterations")
    ax.set_ylabel("Wall time (s)")
    ax.set_title("Solver workload vs wall time")
    ax.grid(True, alpha=0.3)

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_ROOT, f"walltime_vs_newton.{ext}"), bbox_inches="tight", dpi=200)
    plt.close(fig)


# =============================================================================
# SELF-TEST (synthetic data — no simulation outputs required)
# =============================================================================

def _self_test():
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        global SIM_ROOT, FIG_ROOT
        SIM_ROOT = tmp
        FIG_ROOT = os.path.join(tmp, "_figs")

        # Two synthetic cases. Neither has XDMF data; only profiling.json.
        for name, wall, n_newton in [("md_B_baseline", 120.0, 400),
                                     ("md_B_no_transient", 85.0, 310)]:
            folder = os.path.join(tmp, f"dummy_sens_{name}")
            os.makedirs(folder, exist_ok=True)
            payload = {
                "schema_version": 1,
                "case_name": name,
                "wall_time_s": wall,
                "n_steps": 50,
                "n_newton_total": n_newton,
                "summary": {
                    "ct_s_total":    0.25 * wall,
                    "ksp_s_total":   0.35 * wall,
                    "assemble_s_total": 0.20 * wall,
                },
                "per_step": [{"newton_iters": int(np.random.randint(2, 7))} for _ in range(50)],
            }
            with open(os.path.join(folder, "profiling.json"), "w") as f:
                json.dump(payload, f)

        disc = discover_sens_cases(SIM_ROOT)
        assert set(disc) == {"md_B_baseline", "md_B_no_transient"}, f"got {set(disc)}"

        plot_profiling_bars(disc)
        plot_newton_histogram(disc)
        plot_walltime_vs_newton(disc)

        # Discovery with missing profiling.json must still return the record
        folder = os.path.join(tmp, "dummy_sens_no_profile")
        os.makedirs(folder)
        disc = discover_sens_cases(SIM_ROOT)
        assert "no_profile" in disc
        assert disc["no_profile"]["has_profiling"] is False

        print("[self-test] OK — profiling figures written to", FIG_ROOT)


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true", help="Run synthetic-data smoke test and exit.")
    parser.add_argument("--dimension", default=None,
                        help="Render only one dimension (e.g. creep_toggles). Default: all.")
    parser.add_argument("--skip-heavy", action="store_true",
                        help="Skip strain_rate and FOS plots (fast iteration).")
    args = parser.parse_args()

    if args.self_test:
        _self_test()
        return

    discovered = discover_sens_cases(SIM_ROOT)
    if not discovered:
        print(f"No sensitivity cases found under {SIM_ROOT}")
        return

    print(f"Discovered {len(discovered)} case(s):")
    for name, rec in discovered.items():
        flag = "" if rec["has_profiling"] else "  (no profiling.json)"
        print(f"  - {name}{flag}")

    dims = [args.dimension] if args.dimension else list(DIMENSIONS)
    for dim_key in dims:
        if dim_key not in DIMENSIONS:
            warnings.warn(f"Unknown dimension '{dim_key}', skipping")
            continue
        print(f"\n── Rendering dimension: {dim_key} ──")
        plot_convergence_overlay(dim_key, discovered)
        plot_convergence_rate(dim_key, discovered)
        if not args.skip_heavy:
            plot_strain_rate_probes(dim_key, discovered)
            plot_fos_overlay(dim_key, discovered)
            plot_mc_failure_overlay(dim_key, discovered)

    print("\n── Aggregate profiling ──")
    plot_profiling_bars(discovered)
    plot_newton_histogram(discovered)
    plot_walltime_vs_newton(discovered)

    print(f"\nFigures written under {FIG_ROOT}")


if __name__ == "__main__":
    main()
