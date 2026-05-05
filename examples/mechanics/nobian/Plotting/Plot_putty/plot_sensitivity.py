"""
plot_sensitivity.py — Sensitivity-analysis figures for MSc thesis.

Companion to plot_results.py. Consumes `*_sens_<tag>` case folders produced by
`Run_sensitivity.py` and renders per-dimension overlay figures + aggregate
profiling summaries.

Per-dimension overlays (one call per dimension):
    plot_convergence_overlay      — ΔV/V₀ (%) vs time, pressure on twin axis
    plot_convergence_rate         — dΔV/dt, smoothed
    plot_total_strain_rate_probe  — equivalent strain rate from eps_tot at mid-wall
    plot_vp_strain_rate_probe     — equivalent strain rate from eps_vp at mid-wall
    plot_fos_overlay              — global-min FOS vs time
    plot_mc_failure_overlay       — MC-failed interlayer cell count (skips gracefully)

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
    q_dil_devries,
    _psi_from_voigt6,
    tensor33_to_voigt6,
)

# =============================================================================
# STYLE & PATHS
# =============================================================================

plt.rcParams.update({
    "font.size":        14,
    "axes.titlesize":   16,
    "axes.labelsize":   14,
    "xtick.labelsize":  12,
    "ytick.labelsize":  12,
    "legend.fontsize":  12,
    "figure.titlesize": 18,
    "lines.linewidth":  2.0,
})

# =============================================================================
# CANONICAL CASE LABELS (used across every figure)
# =============================================================================
# Map raw case names → human-readable labels used in legends/axes throughout.
CASE_LABELS = {
    "md_B_baseline":        "MunsonDawson_baseline",
    "md_B_no_transient":    "MunsonDawson_noHardening",
    "md_B_no_reverse":      "MunsonDawson_noRecovery",
    "md_B_steady_state":    "MunsonDawson_SteadyState",
    "sic_B_baseline":       "TUD2023_baseline",
    "no_kelvin":            "TUD2023_noKelvin",
    "no_desai":             "TUD2023_noDesai",
    "no_kelvin_no_desai":   "TUD2023_SteadyState",
    "no_pressure_solution": "no_pressure_solution",
    "md_A_baseline":        "MunsonDawson_A",
    "no_leaching":          "no_leaching",
    "tilted_cavern":        "tilted_cavern",
    "volume_600":           "volume_600",
    "dt_half":              "dt_half",
    "dt_double":            "dt_double",
    "amplitude_half":       "amplitude_half",
    "amplitude_double":     "amplitude_double",
}


def case_label(name):
    """Return the canonical display label for a raw case name."""
    return CASE_LABELS.get(name, name)

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
        "labels": {k: case_label(k) for k in
                   ["md_B_baseline", "md_B_no_transient", "md_B_no_reverse", "md_B_steady_state"]},
    },
    "constitutive": {
        "title": "Constitutive model: MD vs TUD2023, scenario A vs B",
        # Order matters: TUD2023 (solid) drawn first so MD (dashed) sits on top.
        "cases": ["sic_B_baseline", "md_A_baseline", "md_B_baseline"],
        "labels": {
            "md_B_baseline":  "MD",
            "sic_B_baseline": "TUD2023",
            "md_A_baseline":  case_label("md_A_baseline"),
        },
    },
    "startup": {
        "title": "Startup: leaching vs equilibrium",
        "cases": ["md_B_baseline", "no_leaching"],
        "labels": {
            "md_B_baseline": case_label("md_B_baseline"),
            "no_leaching":   "Equilibrium startup",
        },
    },
    "geometry": {
        "title": "Cavern geometry & volume",
        "cases": ["md_B_baseline", "tilted_cavern", "volume_600"],
        "labels": {
            "md_B_baseline":  f"{case_label('md_B_baseline')} (Regular 1,200,000 m³)",
            "tilted_cavern":  "Tilted 1,200,000 m³",
            "volume_600":     "Regular 600,000 m³",
        },
    },
    "secondary_creep": {
        "title": "Secondary creep mechanisms",
        # Baseline MD removed per user request — TUD2023 baseline is the reference.
        "cases": ["sic_B_baseline", "no_pressure_solution", "no_kelvin", "no_desai", "no_kelvin_no_desai"],
        "labels": {k: case_label(k) for k in
                   ["sic_B_baseline", "no_pressure_solution", "no_kelvin", "no_desai", "no_kelvin_no_desai"]},
    },
    "timestep": {
        "title": "Numerical time-step",
        "cases": ["md_B_baseline", "dt_half", "dt_double"],
        "labels": {
            "md_B_baseline": f"{case_label('md_B_baseline')} (dt = 2 h)",
            "dt_half":       "dt = 1 h",
            "dt_double":     "dt = 4 h",
        },
    },
    "amplitude": {
        "title": "Pressure amplitude scaling",
        "cases": ["md_B_baseline", "amplitude_half", "amplitude_double"],
        "labels": {
            "md_B_baseline":     f"{case_label('md_B_baseline')} (1.0× amplitude)",
            "amplitude_half":    "0.5× amplitude",
            "amplitude_double":  "2.0× amplitude",
        },
    },
}

# Per-dimension explicit linestyle override:
# secondary_creep — keep TUD2023_baseline solid, all others dashed/dotted variants
# so subtle differences become readable.
# constitutive — TUD2023 solid, MD dashed for direct visual comparison.
DIMENSION_LINESTYLES = {
    "secondary_creep": {
        "sic_B_baseline":       "-",     # baseline solid
        "no_pressure_solution": "--",
        "no_kelvin":            "-.",
        "no_desai":             ":",
        "no_kelvin_no_desai":   (0, (3, 1, 1, 1, 1, 1)),  # dash-dot-dot
    },
    "constitutive": {
        "sic_B_baseline": "-",
        "md_B_baseline":  "--",
    },
}

# Per-dimension explicit colour override.
# constitutive — TUD2023 green, MD red.
DIMENSION_COLORS = {
    "constitutive": {
        "sic_B_baseline": "tab:green",
        "md_B_baseline":  "tab:red",
    },
}

# Discrete colour cycle — consistent across panels of the same dimension.
_COLOR_CYCLE = plt.get_cmap("tab10").colors


def _case_style(idx, dim_key=None, case_name=None):
    """Return (color, linestyle) for a case.

    If `dim_key` has an explicit override in DIMENSION_LINESTYLES /
    DIMENSION_COLORS and `case_name` is present, use the prescribed value.
    Otherwise cycle.
    """
    ls_cycle = ["-", "--", "-.", ":"]
    color = _COLOR_CYCLE[idx % len(_COLOR_CYCLE)]
    color_override = DIMENSION_COLORS.get(dim_key, {}) if dim_key else {}
    if case_name in color_override:
        color = color_override[case_name]
    ls_override = DIMENSION_LINESTYLES.get(dim_key, {}) if dim_key else {}
    if case_name in ls_override:
        return color, ls_override[case_name]
    return color, ls_cycle[idx % len(ls_cycle)]


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


def _path_strain_xdmf(case_folder, field_name):
    """Path to a DG0 strain tensor XDMF (`eps_tot` or `eps_vp`)."""
    candidates = [
        os.path.join(case_folder, "operation", field_name, f"{field_name}.xdmf"),
        os.path.join(case_folder, "operation", f"{field_name}.xdmf"),
    ]
    for p in candidates:
        if os.path.isfile(p):
            return p
    return None


def compute_eq_strain_rate_timeseries(case_folder, probe_xyz, field_name="eps_tot"):
    """
    Equivalent (von Mises) strain rate at the cell closest to `probe_xyz`,
    computed by finite-differencing the equivalent strain stored in the
    DG0 tensor field `<field_name>` (e.g. `eps_tot`, `eps_vp`).

    Equivalent strain:
        ε_eq = sqrt( (2/3) · ε_dev : ε_dev )

    where ε_dev = ε − (1/3) tr(ε) I. The rate is dε_eq / dt in 1/s.

    Returns (t_days, eps_rate). Both arrays have length Nt; the first sample
    is NaN (no backward difference available).
    """
    xdmf = _path_strain_xdmf(case_folder, field_name)
    if xdmf is None:
        return None, None

    centroids, time_list, eps = post.read_cell_tensor(xdmf)
    t = np.asarray(time_list, float)
    if len(t) < 2:
        return None, None

    idx = post.find_closest_point(np.asarray(probe_xyz, float), centroids)
    e = eps[:, idx, :, :]                              # (Nt, 3, 3)
    tr = e[:, 0, 0] + e[:, 1, 1] + e[:, 2, 2]
    dev = e.copy()
    dev[:, 0, 0] -= tr / 3.0
    dev[:, 1, 1] -= tr / 3.0
    dev[:, 2, 2] -= tr / 3.0
    eps_eq = np.sqrt((2.0 / 3.0) * np.einsum("tij,tij->t", dev, dev))

    deps = np.diff(eps_eq)
    dt = np.diff(t)
    dt = np.where(dt > 0, dt, np.nan)
    rate = deps / dt
    rate_full = np.concatenate([[np.nan], rate])
    return t / DAY, rate_full


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


def _legend_on_top(ax, ncol=None):
    """Place the legend above the axes, centred."""
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    n = ncol if ncol is not None else min(len(handles), 3)
    ax.legend(handles, labels, loc="lower center",
              bbox_to_anchor=(0.5, 1.02), ncol=n,
              frameon=False, handlelength=2.5)


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
        color, ls = _case_style(i, dim_key=dim_key, case_name=name)
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
    ax.set_ylabel("ΔV/V₀ (%)")
    ax.grid(True, alpha=0.3)
    _legend_on_top(ax)
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
        color, ls = _case_style(i, dim_key=dim_key, case_name=name)
        label = spec["labels"].get(name, name)
        ax.plot(t_d, rate, color=color, linestyle=ls, label=label)

    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Convergence rate (%/day)")
    ax.grid(True, alpha=0.3)
    _legend_on_top(ax)
    _save(fig, dim_key, "convergence_rate")


def _plot_eq_strain_rate(dim_key, discovered, field_name, ylabel, stem):
    """Mid-wall equivalent strain rate from a DG0 strain field."""
    spec, cases = _resolve_dim_cases(dim_key, discovered)
    if not cases:
        return

    fig, ax = plt.subplots(figsize=(6.3, 4.0))
    any_drawn = False

    for i, (name, rec) in enumerate(cases):
        try:
            wall = load_wall_points(rec["case_path"])
        except Exception as e:
            warnings.warn(f"[{dim_key}/{name}] wall points failed: {e}")
            continue

        z = wall[:, 2]
        z_min, z_max = float(z.min()), float(z.max())
        mid_probe = wall[int(np.argmin(np.abs(z - 0.5 * (z_min + z_max))))]

        t_d, rate = compute_eq_strain_rate_timeseries(
            rec["case_path"], mid_probe, field_name=field_name
        )
        if t_d is None:
            warnings.warn(f"[{dim_key}/{name}] {field_name} not available — skipping")
            continue
        color, ls = _case_style(i, dim_key=dim_key, case_name=name)
        label = spec["labels"].get(name, name)
        # Plot magnitude on log scale
        ax.plot(t_d, np.abs(rate), color=color, linestyle=ls, label=label)
        any_drawn = True

    if not any_drawn:
        plt.close(fig)
        return

    ax.set_xlabel("Time (days)")
    ax.set_ylabel(ylabel)
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    _legend_on_top(ax)
    _save(fig, dim_key, stem)


def plot_total_strain_rate_probe(dim_key, discovered):
    """Mid-wall total equivalent strain rate (from `eps_tot`)."""
    _plot_eq_strain_rate(
        dim_key, discovered,
        field_name="eps_tot",
        ylabel="Mid-wall total strain rate (1/s)",
        stem="strain_rate_total",
    )


def plot_vp_strain_rate_probe(dim_key, discovered):
    """Mid-wall viscoplastic equivalent strain rate (from `eps_vp`)."""
    _plot_eq_strain_rate(
        dim_key, discovered,
        field_name="eps_vp",
        ylabel="Mid-wall viscoplastic strain rate (1/s)",
        stem="strain_rate_vp",
    )


def _compute_fos_timeseries(case_folder):
    """De Vries (2005) field-minimum FOS time series.

    Mirrors the pipeline used in plot_results.py so the sensitivity overlay
    sits on the same FoS scale as the rest of the thesis.
    """
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
    t_days = t_pq[:n] / DAY

    p_MPa = -p_elems[:n] / MPA              # compression positive
    q_MPa = q_elems[:n] / MPA
    sig33_MPa = -sig_vals[:n] / MPA         # compression positive

    sig_voigt = tensor33_to_voigt6(sig33_MPa)
    psi = _psi_from_voigt6(sig_voigt)
    q_boundary = q_dil_devries(p_MPa, psi)

    q_safe = np.where(q_MPa < 1e-3, 1e-3, q_MPa)
    FOS = q_boundary / q_safe
    FOS = np.where(q_MPa < 1e-3, 100.0, FOS)

    fos_min = np.nanmin(np.where(np.isfinite(FOS), FOS, np.nan), axis=1)
    return t_days, fos_min


def plot_fos_overlay(dim_key, discovered):
    spec, cases = _resolve_dim_cases(dim_key, discovered)
    if not cases:
        return
    fig, ax = plt.subplots(figsize=(6.3, 4.0))
    for i, (name, rec) in enumerate(cases):
        t, fos = _compute_fos_timeseries(rec["case_path"])
        if t is None:
            continue
        color, ls = _case_style(i, dim_key=dim_key, case_name=name)
        label = spec["labels"].get(name, name)
        ax.plot(t, fos, color=color, linestyle=ls, label=label)

    ax.axhline(1.0, color="k", linewidth=0.8, linestyle=":")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Global-min FOS")
    ax.grid(True, alpha=0.3)
    _legend_on_top(ax)
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

        color, ls = _case_style(i, dim_key=dim_key, case_name=name)
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
    ax.grid(True, alpha=0.3)
    _legend_on_top(ax)
    _save(fig, dim_key, "mc_failure_overlay")


# =============================================================================
# AGGREGATE PROFILING PLOTS
# =============================================================================

def plot_profiling_bars(discovered):
    """Stacked horizontal bars: ct / ksp / assemble / other CPU hours per case.

    Time is reported in hours (same accumulated CPU time recorded by the
    profiler — this is wall time spent inside the wrapped methods).
    """
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
    names = [case_label(r[0]) for r in rows]
    HR = 3600.0
    ct = np.array([r[1] for r in rows]) / HR
    asmbl = np.array([r[2] for r in rows]) / HR
    ksp = np.array([r[3] for r in rows]) / HR
    other = np.array([r[4] for r in rows]) / HR

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
    ax.set_xlabel("CPU time (h)")
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
        ax.set_title(case_label(name), fontsize=10)
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
    """Scatter of CPU time (hours) vs total Newton iterations — one point per case.

    Each case is coloured and labelled in a legend on top of the plot, so points
    that overlap remain identifiable without text running off the figure.
    Note: this is the wall-clock time of the run (single-process), which equals
    CPU time for a serial solve.
    """
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

    HR = 3600.0
    fig, ax = plt.subplots(figsize=(7.5, 5.0))

    n_pts = len(pts)
    cmap = plt.get_cmap("tab10") if n_pts <= 10 else plt.get_cmap("tab20")
    for i, (name, wall, total_n) in enumerate(pts):
        ax.scatter(total_n, wall / HR, s=80,
                   color=cmap(i % cmap.N),
                   label=case_label(name),
                   edgecolor="black", linewidth=0.5)

    ax.set_xlabel("Total Newton iterations")
    ax.set_ylabel("CPU time (h)")
    ax.grid(True, alpha=0.3)
    _legend_on_top(ax, ncol=min(n_pts, 3))

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
            plot_total_strain_rate_probe(dim_key, discovered)
            plot_vp_strain_rate_probe(dim_key, discovered)
            plot_fos_overlay(dim_key, discovered)
            plot_mc_failure_overlay(dim_key, discovered)

    print("\n── Aggregate profiling ──")
    plot_profiling_bars(discovered)
    plot_newton_histogram(discovered)
    plot_walltime_vs_newton(discovered)

    print(f"\nFigures written under {FIG_ROOT}")


if __name__ == "__main__":
    main()
