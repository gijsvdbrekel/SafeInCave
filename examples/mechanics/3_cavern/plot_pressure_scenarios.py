import numpy as np
import matplotlib.pyplot as plt

import safeincave.Utils as ut

# Import the schedule builders and DAY_H from your main script
# (rename "Test" if your main file has another name)
from Test import (
    build_linear_schedule_multi,
    build_irregular_schedule_multi,
    build_sinus_schedule_multi,
    DAY_H,
)


# --- minimal "time controller" so the builders work --------------------------
class DummyTC:
    """
    Minimal stand-in for SafeInCave's TimeController.
    Only needs dt and t_final in seconds.
    """
    def __init__(self, dt_hours, total_hours):
        self.dt = dt_hours * ut.hour
        self.t_final = total_hours * ut.hour


# === USER SETTINGS ==========================================================

OPERATION_DAYS = 10.0      # total duration in days
DT_HOURS       = 0.1       # time step in hours for sampling

MODE = "repeat"            # "repeat" or "stretch" for multi-day behaviour

# Linear scenario base pattern (per day, in hours and MPa)
base_times_h_linear       = [0.0, 2.0, 14.0, 16.0, 24.0]
base_pressures_MPa_linear = [10.0, 7.0,  7.0, 10.0, 10.0]

# Irregular scenario base pattern (per day, in hours and MPa)
base_waypoints_h_irreg   = [0,  1.0,  2.0,  3.2,  4.0,  5.0,  6.4,
                            7.1, 9.0, 11.5, 13.0, 16.0, 18.0, 21.0, 24.0]
base_pressures_MPa_irreg = [9.0,12.0, 8.5, 11.8, 7.6, 10.2, 8.8,
                            11.4, 9.3, 10.7, 8.9, 11.6, 9.5, 10.2, 11.0]

# Sinus scenario parameters (absolute values in Pa inside the builder)
p_mean_MPa = 10.0     # mean pressure
p_ampl_MPa = 3.0      # amplitude


# === BUILD SCHEDULES ========================================================

tc = DummyTC(DT_HOURS, OPERATION_DAYS * DAY_H)

# Linear
t_lin, p_lin = build_linear_schedule_multi(
    tc,
    base_times_h_linear,
    base_pressures_MPa_linear,
    days=OPERATION_DAYS,
    mode=MODE,
    resample_at_dt=True,
)

# Sinus
t_sin, p_sin = build_sinus_schedule_multi(
    tc,
    p_mean=p_mean_MPa * ut.MPa,
    p_ampl=p_ampl_MPa * ut.MPa,
    days=OPERATION_DAYS,
    mode=MODE,
    daily_period_hours=24.0,
    clamp_min=0.0,
    clamp_max=None,
)

# Irregular
t_irreg, p_irreg = build_irregular_schedule_multi(
    tc,
    base_waypoints_h=base_waypoints_h_irreg,
    base_pressures_MPa=base_pressures_MPa_irreg,
    days=OPERATION_DAYS,
    mode=MODE,
    smooth=0.25,
    clamp_min=0.0,
    clamp_max=None,
    resample_at_dt=True,
)


# === CONVERT TO HOURS & MPa FOR PLOTTING ====================================

def to_hours(t_vals):
    return np.array(t_vals) / ut.hour

def to_MPa(p_vals):
    return np.array(p_vals) / ut.MPa

t_lin_h   = to_hours(t_lin)
p_lin_MPa = to_MPa(p_lin)

t_sin_h   = to_hours(t_sin)
p_sin_MPa = to_MPa(p_sin)

t_irreg_h   = to_hours(t_irreg)
p_irreg_MPa = to_MPa(p_irreg)


# === PLOT ===================================================================

fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
fig.subplots_adjust(hspace=0.25, bottom=0.10)

# Linear
axes[0].plot(t_lin_h, p_lin_MPa, "-", lw=2)
axes[0].set_ylabel("Pressure (MPa)")
axes[0].set_title(f"Linear schedule – {OPERATION_DAYS:.0f} days, mode='{MODE}'")

# Sinus
axes[1].plot(t_sin_h, p_sin_MPa, "-", lw=2)
axes[1].set_ylabel("Pressure (MPa)")
axes[1].set_title(f"Sinusoidal schedule – mean={p_mean_MPa} MPa, ampl={p_ampl_MPa} MPa")

# Irregular
axes[2].plot(t_irreg_h, p_irreg_MPa, "-", lw=2)
axes[2].set_xlabel("Time (hours)")
axes[2].set_ylabel("Pressure (MPa)")
axes[2].set_title(f"Irregular schedule – {OPERATION_DAYS:.0f} days, mode='{MODE}'")

for ax in axes:
    ax.grid(True, alpha=0.3)

plt.show()
