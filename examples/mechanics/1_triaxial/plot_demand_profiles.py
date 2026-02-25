#!/usr/bin/env python3
"""
Generate 4 pressure scheme figures for sector-specific demand profiles.
All include linear leaching + debrining init phase.

  1. Combined 2035 projected (CSV) — full year
  2. Industry — sinusoidal, ~12 cycles/year, full year
  3. Power generation — 10 abrupt withdrawal cycles over 30-day operation window
  4. Transport — regular cycling, ~15 cycles over 30-day operation window
"""

import os
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
# SHARED SETTINGS
# =============================================================================
P_LITHOSTATIC = 26.0
LEACHING_END_FRACTION = 0.40
LEACHING_DAYS = 91.0
LEACHING_DT_HOURS = 12.0
DEBRINING_DAYS = 10.0
RAMP_UP_HOURS = 336.0          # 2-week fade-in
DT_HOURS = 0.5

P_LEACH_END = LEACHING_END_FRACTION * P_LITHOSTATIC   # 10.4 MPa

OUT_DIR = os.path.dirname(__file__)


# =============================================================================
# INIT PHASE
# =============================================================================

def build_init_phase():
    """Build linear leaching + debrining plateau."""
    t_leach, p_leach = build_leaching_schedule_mpa(
        leaching_days=LEACHING_DAYS,
        leaching_dt_hours=LEACHING_DT_HOURS,
        p_start_mpa=P_LITHOSTATIC,
        p_end_mpa=P_LEACH_END,
        mode="linear",
        stepped_n_steps=6,
    )
    t_deb = _sample_times_hours(DEBRINING_DAYS * DAY_H, DT_HOURS)
    p_deb = np.full_like(t_deb, P_LEACH_END)
    return t_leach, p_leach, t_deb, p_deb


def stitch(t_leach, p_leach, t_deb, p_deb, t_op, p_op):
    """Stitch leaching + debrining + operation into one timeline."""
    t_all = list(t_leach)
    p_all = list(p_leach)
    cursor = float(t_leach[-1])

    t_deb_shift = t_deb + cursor
    t_all.extend(t_deb_shift[1:].tolist())
    p_all.extend(p_deb[1:].tolist())
    cursor = float(t_deb_shift[-1])

    t_op_shift = t_op + cursor
    t_all.extend(t_op_shift[1:].tolist())
    p_all.extend(p_op[1:].tolist())

    return np.asarray(t_all), np.asarray(p_all)


# =============================================================================
# SCENARIO BUILDERS
# =============================================================================

def build_combined_2035(operation_days):
    """CSV-based 2035 projected profile."""
    csv_path = os.path.join(os.path.dirname(__file__), "drukprofiel_zoutcaverne_2035_8760u.csv")
    t_op, p_op = build_csv_operation_mpa(
        operation_days=operation_days,
        dt_hours=DT_HOURS,
        schedule_mode="stretch",
        n_cycles=1,
        csv_file=csv_path,
        rescale=False,
        rescale_min=0, rescale_max=0,
    )
    shift = P_LEACH_END - float(np.min(p_op))
    p_op += shift
    apply_fade_in(t_op, p_op, p_start_mpa=P_LEACH_END, fade_in_hours=RAMP_UP_HOURS)
    return t_op, p_op


def build_industry(operation_days):
    """Sinusoidal, ~12 cycles over the operation period."""
    p_mean = P_LEACH_END + 4.0   # 14.4 MPa
    p_ampl = 3.5                  # swing: 10.9 — 17.9 MPa
    t_op, p_op = build_sinus_operation_mpa(
        operation_days=operation_days,
        dt_hours=DT_HOURS,
        schedule_mode="stretch",
        n_cycles=12,
        p_mean_mpa=p_mean,
        p_ampl_mpa=p_ampl,
    )
    apply_fade_in(t_op, p_op, p_start_mpa=P_LEACH_END, fade_in_hours=RAMP_UP_HOURS)
    return t_op, p_op


def build_power_generation(operation_days):
    """10 abrupt withdrawal events with gradual re-pressurisation (tau = 48 h)."""
    total_h = operation_days * 24.0
    t_h = np.arange(0, total_h + DT_HOURS, DT_HOURS)
    n = len(t_h)

    p_base = P_LEACH_END + 7.0   # 17.4 MPa resting pressure
    p = np.full(n, p_base)

    # 10 withdrawal events spread over 30 days, semi-random timing
    rng = np.random.RandomState(42)
    n_events = 10
    recovery_tau = 48.0            # gradual re-pressurisation (hours)
    # Space events roughly evenly, then jitter
    event_centers_days = np.linspace(1.0, operation_days - 1.0, n_events)
    event_centers_days += rng.uniform(-0.8, 0.8, size=n_events)

    for day_c in event_centers_days:
        t_start = day_c * 24.0
        duration = rng.uniform(2.0, 5.0)      # sustained low: 2-5 hours
        depth = rng.uniform(3.5, 6.5)          # drop: 3.5-6.5 MPa

        for i, t in enumerate(t_h):
            if t < t_start:
                continue
            dt_event = t - t_start

            if dt_event < 0.5:
                # Sharp drop (30 min)
                drop = depth * (dt_event / 0.5)
            elif dt_event < 0.5 + duration:
                # Sustained low
                drop = depth
            else:
                # Gradual exponential recovery (tau = 48 hours)
                t_since = dt_event - (0.5 + duration)
                drop = depth * np.exp(-t_since / recovery_tau)
                if drop < 0.05:
                    break

            p[i] = min(p[i], p_base - drop)

    apply_fade_in(t_h, p, p_start_mpa=P_LEACH_END, fade_in_hours=min(RAMP_UP_HOURS, total_h * 0.3))
    return t_h, p


def build_transport(operation_days):
    """Regular trapezoidal cycling, ~15 cycles over 30 days (1 cycle per 2 days)."""
    total_h = operation_days * 24.0
    t_h = np.arange(0, total_h + DT_HOURS, DT_HOURS)

    p_high = P_LEACH_END + 5.0   # 15.4 MPa
    p_low  = P_LEACH_END + 1.0   # 11.4 MPa

    # Cycle period = 2 days = 48 hours
    # Pattern within one 48h cycle:
    #   0-8h:   plateau at high (nighttime, full)
    #   8-12h:  ramp down (morning withdrawal)
    #   12-28h: plateau at low (daytime demand)
    #   28-32h: ramp up (evening injection)
    #   32-48h: plateau at high (nighttime)
    cycle_h = 48.0
    base_t = np.array([0.0, 8.0, 12.0, 28.0, 32.0, 48.0])
    base_p = np.array([p_high, p_high, p_low, p_low, p_high, p_high])

    t_mod = np.mod(t_h, cycle_h)
    p = np.interp(t_mod, base_t, base_p)

    apply_fade_in(t_h, p, p_start_mpa=P_LEACH_END, fade_in_hours=min(RAMP_UP_HOURS, total_h * 0.3))
    return t_h, p


# =============================================================================
# PLOTTING
# =============================================================================

def plot_profile(t_h, p_mpa, total_days, title, filename, color, zoom_operation=False):
    """Plot pressure profile with phase markers and optional operation zoom."""
    t_days = t_h / 24.0
    op_start = LEACHING_DAYS + DEBRINING_DAYS

    if zoom_operation:
        # Two-panel layout: full timeline (top), zoomed operation (bottom)
        fig, (ax, ax_zoom) = plt.subplots(2, 1, figsize=(13, 8),
                                           gridspec_kw={"height_ratios": [1, 1.3], "hspace": 0.30})
    else:
        fig, ax = plt.subplots(figsize=(13, 5))

    # --- Main panel: full timeline ---
    ax.plot(t_days, p_mpa, linewidth=1.2, color=color)
    ax.axvline(LEACHING_DAYS, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axvline(op_start, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)

    y_lo, y_hi = ax.get_ylim()
    y_top = y_hi - 0.03 * (y_hi - y_lo)
    ax.text(LEACHING_DAYS / 2, y_top, "Leaching", ha="center", va="top", fontsize=9, color="gray")
    ax.text(LEACHING_DAYS + DEBRINING_DAYS / 2, y_top, "Debr.", ha="center", va="top", fontsize=9, color="gray")
    op_center = (op_start + total_days) / 2
    ax.text(op_center, y_top, "Operation", ha="center", va="top", fontsize=9, color="gray")

    ax.set_ylabel("Cavern pressure (MPa)")
    ax.set_title(title)
    ax.set_xlim(0, total_days)
    ax.grid(True, alpha=0.3)
    if not zoom_operation:
        ax.set_xlabel("Time (days)")

    # --- Zoom panel: operation phase only ---
    if zoom_operation:
        op_mask = t_days >= op_start
        t_op_days = t_days[op_mask] - op_start
        p_op = p_mpa[op_mask]
        ax_zoom.plot(t_op_days, p_op, linewidth=1.5, color=color)
        ax_zoom.set_xlabel("Operation time (days)")
        ax_zoom.set_ylabel("Cavern pressure (MPa)")
        ax_zoom.set_xlim(0, total_days - op_start)
        ax_zoom.grid(True, alpha=0.3)
        ax_zoom.set_title("Operation phase (zoomed)", fontsize=10)

    fig.tight_layout()
    outpath = os.path.join(OUT_DIR, filename)
    fig.savefig(outpath, dpi=200)
    print(f"[SAVED] {outpath}")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    t_leach, p_leach, t_deb, p_deb = build_init_phase()

    # --- 1. Combined 2035: full year (365 days) ---
    total_1 = 365.0
    op_1 = total_1 - LEACHING_DAYS - DEBRINING_DAYS
    t_op, p_op = build_combined_2035(op_1)
    t_all, p_all = stitch(t_leach, p_leach, t_deb, p_deb, t_op, p_op)
    plot_profile(t_all, p_all, total_1,
                 "Combined 2035 projected demand — all sectors",
                 "pressure_combined_2035.png", "#2c3e50")

    # --- 2. Industry: full year (365 days), 12 sinusoidal cycles ---
    total_2 = 365.0
    op_2 = total_2 - LEACHING_DAYS - DEBRINING_DAYS
    t_op, p_op = build_industry(op_2)
    t_all, p_all = stitch(t_leach, p_leach, t_deb, p_deb, t_op, p_op)
    plot_profile(t_all, p_all, total_2,
                 "Industry — steady baseload, sinusoidal supply variation",
                 "pressure_industry.png", "#1f77b4")

    # --- 3. Power generation: leaching + debrining + 30 days operation ---
    op_3 = 30.0
    total_3 = LEACHING_DAYS + DEBRINING_DAYS + op_3
    t_op, p_op = build_power_generation(op_3)
    t_all, p_all = stitch(t_leach, p_leach, t_deb, p_deb, t_op, p_op)
    plot_profile(t_all, p_all, total_3,
                 "Power generation — abrupt withdrawal peaks, gradual re-pressurisation",
                 "pressure_power_generation.png", "#d62728", zoom_operation=True)

    # --- 4. Transport: leaching + debrining + 30 days operation ---
    op_4 = 30.0
    total_4 = LEACHING_DAYS + DEBRINING_DAYS + op_4
    t_op, p_op = build_transport(op_4)
    t_all, p_all = stitch(t_leach, p_leach, t_deb, p_deb, t_op, p_op)
    plot_profile(t_all, p_all, total_4,
                 "Transport — regular daily refuelling cycle",
                 "pressure_transport.png", "#2ca02c", zoom_operation=True)

    print("\n[DONE] All 4 profiles saved.")


if __name__ == "__main__":
    main()
