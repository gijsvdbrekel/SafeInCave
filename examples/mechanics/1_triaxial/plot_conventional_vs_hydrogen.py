#!/usr/bin/env python3
"""
Two thesis figures contrasting conventional natural-gas storage operation with
the projected 2035 hydrogen storage operation in a Zuidwending salt cavern.

1. conventional_pressure.png
   The conventional gas-storage cavern pressure over one year: a seasonal
   storage cycle (fill in summer, withdraw in winter) with occasional weekly
   variation. This is a *representative / schematic* profile — there is no
   measured Zuidwending conventional-gas pressure record in the repository.

2. conventional_vs_hydrogen_pressure.png
   The same conventional profile with the 2035 hydrogen projection
   (drukprofiel_zoutcaverne_2035_8760u.csv) overlaid, to show how much more
   dynamic hydrogen operation is over the same allowable pressure band.

Both profiles share the same 6-20 MPa allowable band, so the only difference
shown is the cycling pattern, not the pressure range.

Usage:
    python plot_conventional_vs_hydrogen.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION
# =============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(SCRIPT_DIR, "drukprofiel_zoutcaverne_2035_8760u.csv")
OUT_DIR = SCRIPT_DIR

HOURS_PER_YEAR = 8760
DAY_H = 24.0

# Allowable pressure band (MPa), shared by both profiles.
P_MIN, P_MAX = 6.0, 20.0

# Conventional seasonal cycle.
SEASONAL_MEAN = 13.0     # MPa
SEASONAL_AMPL = 5.5      # MPa  (-> ~7.5 to 18.5 MPa, leaving room for ripples)
SEASONAL_PHASE = 0.15    # year fraction of the seasonal minimum (~late Feb)

# Occasional weekly variation: (centre day-of-year, gaussian width days, amplitude MPa).
# The weekly ripple is only active inside these windows; elsewhere the curve is
# the smooth seasonal cycle.
WEEKLY_WINDOWS = [
    (18,  12, 1.2),   # mid-January winter demand
    (46,  10, 1.0),   # late-February cold spell
    (300, 14, 1.1),   # late-October / autumn balancing
]

COLOR_CONV = "#1f77b4"   # conventional gas (blue)
COLOR_H2 = "#d62728"     # 2035 hydrogen (red)

# Month ticks for the time axis.
MONTH_DAYS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


# =============================================================================
# DATA
# =============================================================================

def load_h2_projection():
    """Read the 2035 hydrogen pressure projection (hour, MPa)."""
    hours, mpa = [], []
    with open(CSV_PATH) as f:
        f.readline()  # header
        for line in f:
            parts = line.strip().split(";")
            if len(parts) < 3:
                continue
            hours.append(float(parts[0]))
            mpa.append(float(parts[2].replace(",", ".")))
    t_h = np.asarray(hours) - 1.0          # 1-based hour -> 0-based
    return t_h, np.asarray(mpa)


def build_conventional(t_h):
    """Representative conventional gas-storage pressure profile.

    One seasonal cycle (minimum in late winter, maximum in late summer/autumn)
    with occasional weekly variation superimposed inside the WEEKLY_WINDOWS.
    Schematic, not measured data.
    """
    yf = t_h / HOURS_PER_YEAR              # year fraction 0..1
    p = SEASONAL_MEAN - SEASONAL_AMPL * np.cos(2 * np.pi * (yf - SEASONAL_PHASE))

    weekly = np.sin(2 * np.pi * t_h / (7 * DAY_H))
    envelope = np.zeros_like(t_h)
    for centre_day, width_day, amp in WEEKLY_WINDOWS:
        c = centre_day * DAY_H
        w = width_day * DAY_H
        envelope += amp * np.exp(-0.5 * ((t_h - c) / w) ** 2)

    p = p + weekly * envelope
    return np.clip(p, P_MIN, P_MAX)


# =============================================================================
# PLOTTING
# =============================================================================

def _style_time_axis(ax):
    ax.set_xlim(0, 365)
    ax.set_xticks(MONTH_DAYS)
    ax.set_xticklabels(MONTH_LABELS)
    ax.set_ylim(P_MIN - 0.6, P_MAX + 0.6)
    ax.set_xlabel("Time (months)", fontsize=18)
    ax.set_ylabel("Cavern pressure (MPa)", fontsize=18)
    ax.tick_params(labelsize=14)
    ax.grid(True, alpha=0.25)


def plot_conventional(t_days_conv, p_conv):
    fig, ax = plt.subplots(figsize=(16, 6))
    ax.plot(t_days_conv, p_conv, lw=2.2, color=COLOR_CONV,
            label="Conventional gas storage")
    _style_time_axis(ax)
    ax.legend(loc="upper right", fontsize=16, framealpha=0.9)
    ax.set_title("Conventional gas-storage cavern pressure (Zuidwending)",
                 fontsize=19, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "conventional_pressure.png")
    fig.savefig(path, dpi=300)
    print(f"[SAVED] {path}")
    plt.close(fig)


def plot_overlay(t_days_conv, p_conv, t_days_h2, p_h2):
    fig, ax = plt.subplots(figsize=(16, 6))
    # Hydrogen first (thin, semi-transparent) so its rapid cycling reads as a
    # dense band; conventional drawn on top so the seasonal line stays clear.
    ax.plot(t_days_h2, p_h2, lw=0.7, color=COLOR_H2, alpha=0.65,
            label="2035 hydrogen projection", zorder=2)
    ax.plot(t_days_conv, p_conv, lw=2.4, color=COLOR_CONV,
            label="Conventional gas storage", zorder=3)
    _style_time_axis(ax)
    ax.legend(loc="upper right", fontsize=16, framealpha=0.9)
    ax.set_title("Conventional gas vs. 2035 hydrogen storage operation "
                 "(Zuidwending, same pressure band)",
                 fontsize=19, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "conventional_vs_hydrogen_pressure.png")
    fig.savefig(path, dpi=300)
    print(f"[SAVED] {path}")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    t_h_h2, p_h2 = load_h2_projection()
    t_h_conv = np.arange(0, HOURS_PER_YEAR, 1.0)
    p_conv = build_conventional(t_h_conv)

    t_days_conv = t_h_conv / DAY_H
    t_days_h2 = t_h_h2 / DAY_H

    plot_conventional(t_days_conv, p_conv)
    plot_overlay(t_days_conv, p_conv, t_days_h2, p_h2)
    print("\n[DONE] 2 figures saved.")


if __name__ == "__main__":
    main()
