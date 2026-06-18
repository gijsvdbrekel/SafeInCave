#!/usr/bin/env python3
"""
Two thesis figures for the 2035 hydrogen storage operation of a Zuidwending
salt cavern.

1. projection_2035.png
   The 2035 hydrogen pressure projection alone, drawn as a single thick line
   (same line weight as the conventional-storage figure).

2. projection_2035_with_scenarios.png
   The 2035 projection with the three operational demand scenarios — industry,
   transport, and power generation — overlaid on top, with a legend.

The 2035 projection is read from drukprofiel_zoutcaverne_2035_8760u.csv; the
three scenarios are generated with the same builders used by
plot_demand_profiles.py, so the cycling patterns match the rest of the thesis.
All profiles span the same one-year operation window.

Usage:
    python plot_hydrogen_scenarios.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from plot_demand_profiles import (
    build_industry,
    build_transport,
    build_power_generation,
    OPERATION_DAYS,
    COLORS,
    LABELS,
)

# =============================================================================
# CONFIGURATION
# =============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(SCRIPT_DIR, "drukprofiel_zoutcaverne_2035_8760u.csv")
OUT_DIR = SCRIPT_DIR

DAY_H = 24.0

# 2035 projection drawn as a single thick line (matches the conventional-storage
# figure's line weight).
LW_2035 = 2.6
# Per-scenario line weight / transparency. The high-frequency scenarios
# (transport: ~180 cycles/yr, industry: 36) would otherwise fill their pressure
# band as a solid block over a full year, so they are drawn thin and translucent
# — reading as a banded envelope — while the thick dark 2035 base shows through.
SCENARIO_STYLE = {
    "industry":         {"lw": 0.7, "alpha": 0.60},
    "transport":        {"lw": 0.4, "alpha": 0.40},
    "power_generation": {"lw": 1.1, "alpha": 0.75},
}

# Month ticks for the time axis.
MONTH_DAYS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


# =============================================================================
# DATA
# =============================================================================

def load_h2_projection():
    """Read the 2035 hydrogen pressure projection (time in days, MPa)."""
    hours, mpa = [], []
    with open(CSV_PATH) as f:
        f.readline()  # header
        for line in f:
            parts = line.strip().split(";")
            if len(parts) < 3:
                continue
            hours.append(float(parts[0]))
            mpa.append(float(parts[2].replace(",", ".")))
    t_days = (np.asarray(hours) - 1.0) / DAY_H   # 1-based hour -> 0-based days
    return t_days, np.asarray(mpa)


def load_scenarios():
    """Operation-phase profiles (time in days, MPa) for the three scenarios."""
    out = {}
    for name, builder in [("industry", build_industry),
                          ("transport", build_transport),
                          ("power_generation", build_power_generation)]:
        t_h, p = builder(OPERATION_DAYS)
        out[name] = (np.asarray(t_h) / DAY_H, np.asarray(p))
    return out


# =============================================================================
# PLOTTING
# =============================================================================

def _style_time_axis(ax, y_lo, y_hi):
    ax.set_xlim(0, 365)
    ax.set_xticks(MONTH_DAYS)
    ax.set_xticklabels(MONTH_LABELS)
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel("Time (months)", fontsize=18)
    ax.set_ylabel("Cavern pressure (MPa)", fontsize=18)
    ax.tick_params(labelsize=14)
    ax.grid(True, alpha=0.25)


def plot_2035_only(t_days, p):
    fig, ax = plt.subplots(figsize=(16, 6))
    ax.plot(t_days, p, lw=LW_2035, color=COLORS["2035_projection"],
            label=LABELS["2035_projection"])
    _style_time_axis(ax, p.min() - 0.8, p.max() + 0.8)
    ax.legend(loc="upper right", fontsize=16, framealpha=0.9)
    ax.set_title("2035 hydrogen storage pressure projection (Zuidwending)",
                 fontsize=19, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "projection_2035.png")
    fig.savefig(path, dpi=300)
    print(f"[SAVED] {path}")
    plt.close(fig)


def plot_2035_with_scenarios(t_days, p, scenarios):
    fig, ax = plt.subplots(figsize=(16, 6))

    # Track global pressure extent for the y-limits.
    y_lo, y_hi = p.min(), p.max()

    # 2035 projection as the thick base layer.
    ax.plot(t_days, p, lw=LW_2035, color=COLORS["2035_projection"],
            label=LABELS["2035_projection"], zorder=2)

    # Operational scenarios overlaid on top (translucent so the 2035 base shows).
    for name in ["industry", "transport", "power_generation"]:
        t_s, p_s = scenarios[name]
        st = SCENARIO_STYLE[name]
        ax.plot(t_s, p_s, lw=st["lw"], color=COLORS[name],
                alpha=st["alpha"], label=LABELS[name], zorder=3)
        y_lo = min(y_lo, p_s.min())
        y_hi = max(y_hi, p_s.max())

    _style_time_axis(ax, y_lo - 0.8, y_hi + 0.8)
    # Thick, opaque legend handles so the thin/translucent scenario lines remain
    # identifiable in the legend.
    leg = ax.legend(loc="upper right", fontsize=15, framealpha=0.9, ncol=2)
    for line in leg.get_lines():
        line.set_linewidth(2.6)
        line.set_alpha(1.0)
    ax.set_title("2035 hydrogen projection with operational demand scenarios "
                 "(Zuidwending)", fontsize=19, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "projection_2035_with_scenarios.png")
    fig.savefig(path, dpi=300)
    print(f"[SAVED] {path}")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    t_days, p = load_h2_projection()
    scenarios = load_scenarios()

    plot_2035_only(t_days, p)
    plot_2035_with_scenarios(t_days, p, scenarios)
    print("\n[DONE] 2 figures saved.")


if __name__ == "__main__":
    main()
