#!/usr/bin/env python3
"""
Generate two thesis-ready figures for the methodology chapter:

1. pressure_overview.png
   All three operational scenarios overlaid on a single full-timeline plot
   (leaching → debrining → transition → operation), with annotated phases.

2. pressure_ops_zoom.png
   1×3 panel showing a zoomed view of the operational stage for each
   scenario, with consistent axis scaling to compare cycling frequency,
   pressure swing, and loading smoothness.

Parameters are aligned with Run.py defaults.

Usage:
    python plot_demand_profiles.py
"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt

from pressure_scheme import (
    build_leaching_schedule_mpa,
    build_sinus_operation_mpa,
    build_csv_operation_mpa,
    apply_fade_in,
    _sample_times_hours,
    DAY_H,
)

# =============================================================================
# PARAMETERS  (aligned with Run.py)
# =============================================================================
P_LITHOSTATIC    = 26.0       # MPa
LEACHING_END_FR  = 0.40
LEACHING_DAYS    = 91.0
LEACHING_DT_H    = 12.0
DEBRINING_DAYS   = 30.0       # Run.py default
RAMP_UP_HOURS    = 336.0      # 2-week cosine fade-in
OPERATION_DAYS   = 365.0
DT_HOURS         = 0.5

P_LEACH_END = LEACHING_END_FR * P_LITHOSTATIC   # 10.4 MPa
TOTAL_DAYS  = LEACHING_DAYS + DEBRINING_DAYS + OPERATION_DAYS  # 486

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Scenario-specific (from Run.py)
IND_N_CYCLES    = 36
IND_P_AMPL      = 5.0         # P_AMPLITUDE_MPA
IND_P_MEAN      = P_LEACH_END + IND_P_AMPL   # 15.4

TRA_N_CYCLES    = 180
TRA_P_HIGH      = P_LEACH_END + 5.0    # 15.4
TRA_P_LOW       = P_LEACH_END + 1.5    # 11.9

PG_N_EVENTS     = 10
PG_P_BASE       = P_LEACH_END + 10.0   # 20.4
PG_RECOVERY_TAU = 300.0                 # hours
PG_P_MIN        = 6.2                   # absolute floor

CSV_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "drukprofiel_zoutcaverne_2035_8760u.csv")


# =============================================================================
# INIT PHASE  (shared by all scenarios)
# =============================================================================

def build_init_phase():
    t_leach, p_leach = build_leaching_schedule_mpa(
        leaching_days=LEACHING_DAYS,
        leaching_dt_hours=LEACHING_DT_H,
        p_start_mpa=P_LITHOSTATIC,
        p_end_mpa=P_LEACH_END,
        mode="linear",
        stepped_n_steps=6,
    )
    t_deb = _sample_times_hours(DEBRINING_DAYS * DAY_H, DT_HOURS)
    p_deb = np.full_like(t_deb, P_LEACH_END)
    return t_leach, p_leach, t_deb, p_deb


def stitch(t_leach, p_leach, t_deb, p_deb, t_op, p_op):
    t_all = list(t_leach)
    p_all = list(p_leach)
    cursor = float(t_leach[-1])

    t_all.extend((t_deb[1:] + cursor).tolist())
    p_all.extend(p_deb[1:].tolist())
    cursor += float(t_deb[-1])

    t_all.extend((t_op[1:] + cursor).tolist())
    p_all.extend(p_op[1:].tolist())

    return np.asarray(t_all), np.asarray(p_all)


# =============================================================================
# SCENARIO BUILDERS
# =============================================================================

def build_industry(op_days):
    t_op, p_op = build_sinus_operation_mpa(
        operation_days=op_days, dt_hours=DT_HOURS,
        schedule_mode="stretch", n_cycles=IND_N_CYCLES,
        p_mean_mpa=IND_P_MEAN, p_ampl_mpa=IND_P_AMPL,
    )
    apply_fade_in(t_op, p_op, p_start_mpa=P_LEACH_END,
                  fade_in_hours=RAMP_UP_HOURS)
    return t_op, p_op


def build_transport(op_days):
    total_h = op_days * 24.0
    t_h = np.arange(0, total_h + DT_HOURS, DT_HOURS)
    # 48-h trapezoidal cycle (from Run.py)
    cycle_h = total_h / TRA_N_CYCLES
    base_frac = np.array([0.0, 8/48, 12/48, 28/48, 32/48, 1.0])
    base_p = np.array([TRA_P_HIGH, TRA_P_HIGH, TRA_P_LOW,
                       TRA_P_LOW, TRA_P_HIGH, TRA_P_HIGH])
    base_t = base_frac * cycle_h
    t_mod = np.mod(t_h, cycle_h)
    p = np.interp(t_mod, base_t, base_p)
    apply_fade_in(t_h, p, p_start_mpa=P_LEACH_END,
                  fade_in_hours=RAMP_UP_HOURS)
    return t_h, p


def build_power_generation(op_days):
    total_h = op_days * 24.0
    t_h = np.arange(0, total_h + DT_HOURS, DT_HOURS)
    p = np.full_like(t_h, PG_P_BASE)

    rng = np.random.RandomState(42)
    centres = np.linspace(1.0, op_days - 1.0, PG_N_EVENTS)
    centres += rng.uniform(-0.8, 0.8, size=PG_N_EVENTS)

    for day_c in centres:
        t_start = day_c * 24.0
        duration = rng.uniform(2.0, 5.0)
        depth = rng.uniform(6.5, 12.5)
        for i, t in enumerate(t_h):
            if t < t_start:
                continue
            dt_ev = t - t_start
            if dt_ev < 0.5:
                drop = depth * (dt_ev / 0.5)
            elif dt_ev < 0.5 + duration:
                drop = depth
            else:
                t_since = dt_ev - (0.5 + duration)
                drop = depth * np.exp(-t_since / PG_RECOVERY_TAU)
                if drop < 0.05:
                    break
            p[i] = min(p[i], PG_P_BASE - drop)

    p = np.maximum(p, PG_P_MIN)
    apply_fade_in(t_h, p, p_start_mpa=P_LEACH_END,
                  fade_in_hours=RAMP_UP_HOURS)
    return t_h, p


CSV_LEACH_END_FR = 0.25
CSV_P_LEACH_END  = CSV_LEACH_END_FR * P_LITHOSTATIC   # 6.5 MPa


def build_2035_projection(op_days):
    t_op, p_op = build_csv_operation_mpa(
        operation_days=op_days, dt_hours=DT_HOURS,
        schedule_mode="stretch", n_cycles=1,
        csv_file=CSV_PATH, rescale=False,
        rescale_min=0, rescale_max=0,
    )
    # Shift so minimum aligns with the 2035-specific leach end pressure
    shift = CSV_P_LEACH_END - float(np.min(p_op))
    p_op += shift
    apply_fade_in(t_op, p_op, p_start_mpa=CSV_P_LEACH_END,
                  fade_in_hours=RAMP_UP_HOURS)
    return t_op, p_op


def build_2035_init_phase():
    """Separate init phase for the 2035 projection (leaching end fraction 0.25)."""
    t_leach, p_leach = build_leaching_schedule_mpa(
        leaching_days=LEACHING_DAYS,
        leaching_dt_hours=LEACHING_DT_H,
        p_start_mpa=P_LITHOSTATIC,
        p_end_mpa=CSV_P_LEACH_END,
        mode="linear",
        stepped_n_steps=6,
    )
    t_deb = _sample_times_hours(DEBRINING_DAYS * DAY_H, DT_HOURS)
    p_deb = np.full_like(t_deb, CSV_P_LEACH_END)
    return t_leach, p_leach, t_deb, p_deb


# =============================================================================
# PLOTTING
# =============================================================================
# Colours consistent across both figures
COLORS = {
    "2035_projection":  "#2c3e50",
    "industry":         "#1f77b4",
    "power_generation": "#d62728",
    "transport":        "#2ca02c",
}
LABELS = {
    "2035_projection":  "2035 projection",
    "industry":         "Industry",
    "power_generation": "Power generation",
    "transport":        "Transport",
}

OP_START_DAY = LEACHING_DAYS + DEBRINING_DAYS          # 121


def _annotate_phases(ax, y_frac=0.96):
    """Add phase labels and vertical dividers."""
    ax.axvline(LEACHING_DAYS, color="grey", ls="--", lw=0.8, alpha=0.6)
    ax.axvline(OP_START_DAY, color="grey", ls="--", lw=0.8, alpha=0.6)
    # transition end (fade-in)
    trans_end = OP_START_DAY + RAMP_UP_HOURS / 24.0
    ax.axvline(trans_end, color="grey", ls=":", lw=0.7, alpha=0.5)

    yl, yh = ax.get_ylim()
    yt = yl + y_frac * (yh - yl)
    kw = dict(ha="center", va="top", fontsize=8.5, color="grey")
    ax.text(LEACHING_DAYS / 2,              yt, "Leaching",   **kw)
    ax.text(LEACHING_DAYS + DEBRINING_DAYS / 2, yt, "Debrining",  **kw)
    ax.text(OP_START_DAY + RAMP_UP_HOURS / 24 / 2, yt, "Transition", **kw)
    ax.text((trans_end + TOTAL_DAYS) / 2,    yt, "Operation",  **kw)


def plot_overview(scenarios):
    """Figure 1: all scenarios overlaid, full timeline."""
    fig, ax = plt.subplots(figsize=(14, 5))

    for key, (t_h, p) in scenarios.items():
        ax.plot(t_h / 24.0, p, lw=1.2, color=COLORS[key],
                label=LABELS[key], alpha=0.85)

    ax.set_xlim(0, TOTAL_DAYS)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Cavern pressure (MPa)")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.25)
    _annotate_phases(ax)

    fig.tight_layout()
    path = os.path.join(OUT_DIR, "pressure_overview.png")
    fig.savefig(path, dpi=250)
    print(f"[SAVED] {path}")
    plt.close(fig)


def plot_ops_zoom(scenarios):
    """Figure 2: 1×3 zoomed operational panels with consistent y-axis."""
    order = ["industry", "transport", "power_generation"]

    # Determine shared y-limits across all operational stages
    y_lo, y_hi = np.inf, -np.inf
    for key in order:
        t_h, p = scenarios[key]
        mask = t_h / 24.0 >= OP_START_DAY
        p_op = p[mask]
        y_lo = min(y_lo, p_op.min())
        y_hi = max(y_hi, p_op.max())
    pad = 0.06 * (y_hi - y_lo)
    y_lo -= pad
    y_hi += pad

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.2), sharey=True)

    panel_labels = ["(a)", "(b)", "(c)"]

    for ax, key, pl in zip(axes, order, panel_labels):
        t_h, p = scenarios[key]
        t_days = t_h / 24.0
        mask = t_days >= OP_START_DAY
        t_op = t_days[mask] - OP_START_DAY
        p_op = p[mask]

        ax.plot(t_op, p_op, lw=0.9, color=COLORS[key])
        ax.set_xlim(0, OPERATION_DAYS)
        ax.set_ylim(y_lo, y_hi)
        ax.set_xlabel("Operation time (days)")
        ax.set_title(f"{pl}  {LABELS[key]}", fontsize=10)
        ax.grid(True, alpha=0.25)

    axes[0].set_ylabel("Cavern pressure (MPa)")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "pressure_ops_zoom.png")
    fig.savefig(path, dpi=250)
    print(f"[SAVED] {path}")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    t_leach, p_leach, t_deb, p_deb = build_init_phase()
    t_leach_csv, p_leach_csv, t_deb_csv, p_deb_csv = build_2035_init_phase()

    scenarios = {}

    # 2035 projection uses its own init (leaching end fraction 0.25)
    t_op, p_op = build_2035_projection(OPERATION_DAYS)
    t_all, p_all = stitch(t_leach_csv, p_leach_csv, t_deb_csv, p_deb_csv,
                          t_op, p_op)
    scenarios["2035_projection"] = (t_all, p_all)

    # Other scenarios share the standard init (leaching end fraction 0.40)
    for name, builder in [("industry",         build_industry),
                           ("transport",        build_transport),
                           ("power_generation", build_power_generation)]:
        t_op, p_op = builder(OPERATION_DAYS)
        t_all, p_all = stitch(t_leach, p_leach, t_deb, p_deb, t_op, p_op)
        scenarios[name] = (t_all, p_all)

    plot_overview(scenarios)
    plot_ops_zoom(scenarios)
    print("\n[DONE] 2 figures saved.")


if __name__ == "__main__":
    main()
