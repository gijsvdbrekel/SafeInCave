#!/usr/bin/env python3
"""
pressure_scheme_plot_user_config.py

Plot a SafeInCave-style pressure scheme over a 1-year timeline (365 days by default),
with USER OPTIONS at the TOP of the script (no CLI needed).

Includes dashed vertical lines:
- end of leaching phase (if USE_LEACHING)
- end of debrining phase (if USE_LEACHING and DEBRINING_DAYS > 0)

Ubuntu usage:
    python3 pressure_scheme_plot_user_config.py
It will show a window and also save a PNG by default (configurable below).
"""

from __future__ import annotations

import csv
import math
import os
from typing import List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# USER CONFIGURATION (EDIT HERE)
# =============================================================================

# --- OUTPUT ---
SHOW_PLOT = True                    # show interactive window
SAVE_PLOT = False
SAVE_PATH = "pressure_scheme.png"   # saved figure path (if SAVE_PLOT)

# --- TOTAL TIMELINE ---
TOTAL_DAYS = 365.0                  # total plotted time INCLUDING init + debrining + operation

# --- INITIALIZATION MODE ---
# If True: leaching ramp from lithostatic -> leach_end
# If False: equilibrium plateau at constant pressure
USE_LEACHING = True

# --- LITHOSTATIC (MPa) ---
# In your full sim this comes from geometry + density; here you set it directly.
P_LITHOSTATIC_MPA = 26.0

# --- LEACHING PHASE SETTINGS (only used if USE_LEACHING = True) ---
LEACHING_MODE = "stepped"           # "linear" or "stepped"
LEACHING_DAYS = 91.0
LEACHING_DT_HOURS = 12.0
STEPPED_N_STEPS = 6                 # only for stepped
LEACHING_END_FRACTION = 0.40        # p_leach_end = fraction * p_lithostatic

# --- DEBRINING PHASE (only if USE_LEACHING = True) ---
# After leaching: constant plateau at p_leach_end for DEBRINING_DAYS
DEBRINING_DAYS = 30.0               # set 0 to skip

# --- EQUILIBRIUM PHASE (only used if USE_LEACHING = False) ---
EQUILIBRIUM_HOURS = 10.0
EQUILIBRIUM_DT_HOURS = 0.5
P_EQUILIBRIUM_MPA = 15.0

# --- TRANSITION / RAMP-UP ---
# Fade-in at start of OPERATION to avoid abrupt swings:
# blends from p_start (p_leach_end or p_equilibrium) into operational schedule
RAMP_UP_HOURS = 336 # 2 weeks                # set 0 to disable (not recommended)

# --- OPERATION PHASE (fills the rest of TOTAL_DAYS) ---
PRESSURE_SCENARIO = "sinus"         # "sinus", "linear", "irregular", "csv"
SCHEDULE_MODE = "stretch"           # "stretch", "repeat"  (and "direct" only for csv)
N_CYCLES = 5                       # only meaningful for stretch
DT_HOURS = 2.0                      # sampling time step for OPERATION plot

# --- SINUS SETTINGS ---
# If USE_LEACHING=True: mean becomes (p_leach_end + P_AMPLITUDE_MPA), like your sim script
# If USE_LEACHING=False: mean is P_MEAN_MPA
P_MEAN_MPA = 15.0
P_AMPLITUDE_MPA = 6.5

# --- LINEAR SETTINGS ---
# If USE_LEACHING=True: pmin = p_leach_end, pmax = p_leach_end + PRESSURE_SWING_MPA
# If USE_LEACHING=False: pmin = P_MIN_MPA, pmax = P_MAX_MPA
P_MIN_MPA = 6.0
P_MAX_MPA = 20.0
PRESSURE_SWING_MPA = 10.0

# --- IRREGULAR SETTINGS ---
IRREGULAR_SMOOTH = 0.25
IRREGULAR_WAYPOINTS_H = [0, 1.0, 2.0, 3.2, 4.0, 5.0, 6.4, 7.1, 9.0, 11.5, 13.0, 16.0, 18.0, 21.0, 24.0]
IRREGULAR_PRESSURES_MPA = [15.0, 12.0, 8.5, 11.8, 7.6, 10.2, 8.8, 11.4, 9.3, 10.7, 8.9, 11.6, 9.5, 10.2, 11.0]
# If USE_LEACHING=True, these will be shifted so their minimum equals p_leach_end (same idea as your sim)

# --- CSV SETTINGS ---
CSV_FILE_PATH = "drukprofiel_zoutcaverne_2035_8760u.csv"
CSV_SHIFT_TO_LEACH_END = True       # if USE_LEACHING: shift csv so min equals p_leach_end
CSV_RESCALE = False                 # only used when USE_LEACHING=False (mirrors your sim logic)
CSV_RESCALE_MIN_MPA = 6.0
CSV_RESCALE_MAX_MPA = 20.0

# =============================================================================
# END USER CONFIGURATION
# =============================================================================


DAY_H = 24.0


# ---------------------------
# Helpers
# ---------------------------

def _sample_times_hours(t_end_h: float, dt_h: float) -> np.ndarray:
    if dt_h <= 0:
        raise ValueError("dt_hours must be > 0")
    n = int(math.floor(t_end_h / dt_h))
    t = np.array([k * dt_h for k in range(n + 1)], dtype=float)
    if abs(t[-1] - t_end_h) > 1e-12:
        t = np.append(t, t_end_h)
    return t


def _parse_float_auto(s: str) -> float:
    s = s.strip()
    if not s:
        return np.nan
    if "," in s and "." not in s:
        s = s.replace(",", ".")
    if "," in s and "." in s:
        s = s.replace(".", "").replace(",", ".")
    try:
        return float(s)
    except Exception:
        return np.nan


def read_pressure_csv(csv_file: str) -> np.ndarray:
    if not os.path.isfile(csv_file):
        raise FileNotFoundError(f"CSV not found: {csv_file}")

    with open(csv_file, "r", newline="", encoding="utf-8") as f:
        sample = f.read(2048)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=";,\t")
            delim = dialect.delimiter
        except Exception:
            delim = ";"

    rows: List[List[str]] = []
    with open(csv_file, "r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=delim)
        for r in reader:
            if len(r) > 0:
                rows.append(r)

    if len(rows) < 2:
        raise ValueError("CSV has too few rows")

    header = [c.strip() for c in rows[0]]
    data = rows[1:]
    header_low = [h.lower() for h in header]

    idx_mpa = idx_bar = None
    for i, h in enumerate(header_low):
        if h == "druk_mpa" or h.endswith("druk_mpa"):
            idx_mpa = i
        if h == "druk_bar" or h.endswith("druk_bar"):
            idx_bar = i

    pressures_mpa: List[float] = []
    if idx_mpa is not None:
        for r in data:
            if idx_mpa < len(r):
                pressures_mpa.append(_parse_float_auto(r[idx_mpa]))
    elif idx_bar is not None:
        for r in data:
            if idx_bar < len(r):
                pressures_mpa.append(_parse_float_auto(r[idx_bar]) / 10.0)
    else:
        ncols = len(header)
        best_i, best_count = None, -1
        for i in range(ncols):
            vals = [_parse_float_auto(r[i]) for r in data if i < len(r)]
            count = int(np.sum(np.isfinite(vals)))
            if count > best_count:
                best_count, best_i = count, i
        if best_i is None or best_count < 2:
            raise ValueError("Could not find a numeric pressure column in CSV")
        for r in data:
            if best_i < len(r):
                pressures_mpa.append(_parse_float_auto(r[best_i]))

    arr = np.asarray(pressures_mpa, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < 2:
        raise ValueError("Parsed pressure series has <2 numeric values")
    return arr


def rescale_pressure_profile(pressures_mpa: np.ndarray, new_min: float, new_max: float) -> np.ndarray:
    old_min = float(pressures_mpa.min())
    old_max = float(pressures_mpa.max())
    if old_max - old_min < 1e-9:
        return np.full_like(pressures_mpa, (new_min + new_max) / 2.0)
    frac = (pressures_mpa - old_min) / (old_max - old_min)
    return new_min + frac * (new_max - new_min)


def apply_fade_in(t_h: np.ndarray, p_mpa: np.ndarray, *, p_start_mpa: float, fade_in_hours: float) -> None:
    if fade_in_hours is None or fade_in_hours <= 0.0:
        return
    T = float(fade_in_hours)
    for i, t in enumerate(t_h):
        if t >= T:
            break
        alpha = 0.5 * (1.0 - math.cos(math.pi * t / T))
        p_mpa[i] = (1.0 - alpha) * p_start_mpa + alpha * p_mpa[i]


# ---------------------------
# Phase builders (MPa)
# ---------------------------

def build_leaching_schedule_mpa(
    *,
    leaching_days: float,
    leaching_dt_hours: float,
    p_start_mpa: float,
    p_end_mpa: float,
    mode: str,
    stepped_n_steps: int,
) -> Tuple[np.ndarray, np.ndarray]:
    t_end_h = float(leaching_days) * DAY_H
    t_h = _sample_times_hours(t_end_h, float(leaching_dt_hours))

    if mode == "linear":
        frac = np.where(t_end_h > 0, t_h / t_end_h, 1.0)
        p = p_start_mpa + frac * (p_end_mpa - p_start_mpa)
        return t_h, p

    if mode == "stepped":
        n_steps = max(2, int(stepped_n_steps))
        step_duration = t_end_h / n_steps
        step_vals = np.linspace(p_start_mpa, p_end_mpa, n_steps + 1)
        p = np.zeros_like(t_h)
        for i, t in enumerate(t_h):
            step_idx = min(int(t / step_duration), n_steps - 1)
            p[i] = step_vals[step_idx] if t < t_end_h else p_end_mpa
        return t_h, p

    raise ValueError("LEACHING_MODE must be 'linear' or 'stepped'")


def build_equilibrium_schedule_mpa(*, equilibrium_hours: float, equilibrium_dt_hours: float, p_eq_mpa: float) -> Tuple[np.ndarray, np.ndarray]:
    t_end_h = float(equilibrium_hours)
    t_h = _sample_times_hours(t_end_h, float(equilibrium_dt_hours))
    p = np.full_like(t_h, float(p_eq_mpa))
    return t_h, p


def build_sinus_operation_mpa(
    *,
    operation_days: float,
    dt_hours: float,
    schedule_mode: str,
    n_cycles: int,
    p_mean_mpa: float,
    p_ampl_mpa: float,
) -> Tuple[np.ndarray, np.ndarray]:
    total_h = float(operation_days) * DAY_H
    t_h = _sample_times_hours(total_h, float(dt_hours))

    if schedule_mode == "repeat":
        T_h = 24.0
    elif schedule_mode == "stretch":
        T_h = total_h / float(max(1, int(n_cycles)))
    else:
        raise ValueError("SCHEDULE_MODE must be 'repeat' or 'stretch' for sinus")

    two_pi_over_T = 2.0 * math.pi / T_h if T_h > 0 else 0.0
    p = p_mean_mpa + p_ampl_mpa * np.sin(two_pi_over_T * t_h)
    return t_h, p


def build_linear_operation_mpa(
    *,
    operation_days: float,
    dt_hours: float,
    schedule_mode: str,
    n_cycles: int,
    p_min_mpa: float,
    p_max_mpa: float,
) -> Tuple[np.ndarray, np.ndarray]:
    total_h = float(operation_days) * DAY_H
    t_h = _sample_times_hours(total_h, float(dt_hours))

    # daily trapezoid like your sim
    base_t = np.asarray([0.0, 2.0, 14.0, 16.0, 24.0], dtype=float)
    base_p = np.asarray([p_min_mpa, p_max_mpa, p_max_mpa, p_min_mpa, p_min_mpa], dtype=float)

    if schedule_mode == "repeat":
        t_in = np.mod(t_h, 24.0)
        p = np.interp(t_in, base_t, base_p)
        return t_h, p

    if schedule_mode == "stretch":
        cycle_len = total_h / float(max(1, int(n_cycles)))
        scaled_t = base_t * (cycle_len / 24.0)
        t_in = np.mod(t_h, cycle_len)
        p = np.interp(t_in, scaled_t, base_p)
        return t_h, p

    raise ValueError("SCHEDULE_MODE must be 'repeat' or 'stretch' for linear")


def _cardinal_segment(p0, p1, p2, p3, u, tension):
    m1 = (1 - tension) * 0.5 * (p2 - p0)
    m2 = (1 - tension) * 0.5 * (p3 - p1)
    u2 = u * u
    u3 = u2 * u
    h00 = 2*u3 - 3*u2 + 1
    h10 = u3 - 2*u2 + u
    h01 = -2*u3 + 3*u2
    h11 = u3 - u2
    return h00*p1 + h10*m1 + h01*p2 + h11*m2


def _cardinal_interp(ts, ps, t, tension):
    if t <= ts[0]:
        return ps[0]
    if t >= ts[-1]:
        return ps[-1]
    i = np.searchsorted(ts, t) - 1
    t0 = ts[i]
    t1 = ts[i+1]
    i_m1 = max(i-1, 0)
    i_p2 = min(i+2, len(ts)-1)
    u = (t - t0) / (t1 - t0)
    return _cardinal_segment(ps[i_m1], ps[i], ps[i+1], ps[i_p2], u, tension)


def build_irregular_operation_mpa(
    *,
    operation_days: float,
    dt_hours: float,
    schedule_mode: str,
    n_cycles: int,
    base_waypoints_h: List[float],
    base_pressures_mpa: List[float],
    smooth: float,
) -> Tuple[np.ndarray, np.ndarray]:
    base_t = np.asarray(base_waypoints_h, dtype=float)
    base_p = np.asarray(base_pressures_mpa, dtype=float)
    if base_t.size != base_p.size or base_t.size < 2:
        raise ValueError("Irregular waypoints and pressures must match and have >=2 points")

    total_h = float(operation_days) * DAY_H

    # Build multi-day knots
    if schedule_mode == "repeat":
        t_list, p_list = [], []
        n_days = int(math.ceil(operation_days))
        for d in range(n_days):
            off = d * 24.0
            for i in range(base_t.size):
                if d > 0 and i == 0:
                    continue
                tt = off + base_t[i]
                if tt > total_h + 1e-12:
                    break
                t_list.append(tt)
                p_list.append(base_p[i])
        knots_t = np.asarray(t_list, dtype=float)
        knots_p = np.asarray(p_list, dtype=float)

    elif schedule_mode == "stretch":
        cycles = max(1, int(n_cycles))
        base_start = float(base_t[0])
        base_end = float(base_t[-1])
        base_dur = base_end - base_start
        if base_dur <= 0:
            raise ValueError("IRREGULAR_WAYPOINTS_H must span a positive duration")

        cycle_len = total_h / float(cycles)
        scale = cycle_len / base_dur

        t_list, p_list = [], []
        for k in range(cycles):
            off = k * cycle_len
            for i in range(base_t.size):
                if k > 0 and i == 0:
                    continue
                t_scaled = off + (base_t[i] - base_start) * scale
                t_list.append(t_scaled)
                p_list.append(base_p[i])
        knots_t = np.asarray(t_list, dtype=float)
        knots_p = np.asarray(p_list, dtype=float)

    else:
        raise ValueError("SCHEDULE_MODE must be 'repeat' or 'stretch' for irregular")

    # ensure endpoints
    if knots_t[0] > 0.0:
        knots_t = np.insert(knots_t, 0, 0.0)
        knots_p = np.insert(knots_p, 0, knots_p[0])
    if knots_t[-1] < total_h:
        knots_t = np.append(knots_t, total_h)
        knots_p = np.append(knots_p, knots_p[-1])

    # resample at dt_hours
    t_h = _sample_times_hours(total_h, float(dt_hours))
    tau = float(np.clip(smooth, 0.0, 1.0))
    p = np.zeros_like(t_h)
    for i, t in enumerate(t_h):
        p[i] = _cardinal_interp(knots_t, knots_p, t, tau)
    return t_h, p


def build_csv_operation_mpa(
    *,
    operation_days: float,
    dt_hours: float,
    schedule_mode: str,
    n_cycles: int,
    csv_file: str,
    rescale: bool,
    rescale_min: float,
    rescale_max: float,
) -> Tuple[np.ndarray, np.ndarray]:
    pressures = read_pressure_csv(csv_file)  # hourly MPa
    if rescale:
        pressures = rescale_pressure_profile(pressures, rescale_min, rescale_max)

    csv_hours = int(pressures.size)
    total_h = float(operation_days) * DAY_H

    if schedule_mode == "direct":
        knots_t = np.arange(0.0, total_h + 1e-12, 1.0)
        idx = (knots_t % csv_hours).astype(int)
        knots_p = pressures[idx]

    elif schedule_mode == "repeat":
        n_rep = int(np.ceil(total_h / float(csv_hours)))
        t_list, p_list = [], []
        for r in range(n_rep):
            off = r * csv_hours
            for i in range(csv_hours):
                if r > 0 and i == 0:
                    continue
                t = off + i
                if t > total_h:
                    break
                t_list.append(t)
                p_list.append(pressures[i])
        knots_t = np.asarray(t_list, dtype=float)
        knots_p = np.asarray(p_list, dtype=float)

    elif schedule_mode == "stretch":
        cycles = max(1, int(n_cycles))
        cycle_len = total_h / float(cycles)
        scale = cycle_len / float(csv_hours)

        t_list, p_list = [], []
        for k in range(cycles):
            off = k * cycle_len
            for i in range(csv_hours):
                if k > 0 and i == 0:
                    continue
                t_list.append(off + i * scale)
                p_list.append(pressures[i])
        knots_t = np.asarray(t_list, dtype=float)
        knots_p = np.asarray(p_list, dtype=float)

    else:
        raise ValueError("CSV SCHEDULE_MODE must be 'direct', 'repeat', or 'stretch'")

    # ensure endpoints
    if knots_t[0] > 0.0:
        knots_t = np.insert(knots_t, 0, 0.0)
        knots_p = np.insert(knots_p, 0, knots_p[0])
    if knots_t[-1] < total_h:
        knots_t = np.append(knots_t, total_h)
        knots_p = np.append(knots_p, knots_p[-1])

    # resample at dt
    t_h = _sample_times_hours(total_h, float(dt_hours))
    p = np.interp(t_h, knots_t, knots_p)
    return t_h, p


# ---------------------------
# Build full schedule
# ---------------------------

def build_full_schedule_one_year() -> Tuple[np.ndarray, np.ndarray, Optional[float], Optional[float]]:
    # Validate user settings (be strict, otherwise your plot can silently lie)
    if TOTAL_DAYS <= 0:
        raise ValueError("TOTAL_DAYS must be > 0")
    if USE_LEACHING:
        if LEACHING_DAYS <= 0:
            raise ValueError("LEACHING_DAYS must be > 0 when USE_LEACHING=True")
        if LEACHING_DT_HOURS <= 0:
            raise ValueError("LEACHING_DT_HOURS must be > 0")
        if LEACHING_MODE not in ("linear", "stepped"):
            raise ValueError("LEACHING_MODE must be 'linear' or 'stepped'")
        if not (0.0 < LEACHING_END_FRACTION < 1.0):
            raise ValueError("LEACHING_END_FRACTION must be in (0,1)")
        if DEBRINING_DAYS < 0:
            raise ValueError("DEBRINING_DAYS must be >= 0")
    else:
        if EQUILIBRIUM_HOURS <= 0:
            raise ValueError("EQUILIBRIUM_HOURS must be > 0 when USE_LEACHING=False")
        if EQUILIBRIUM_DT_HOURS <= 0:
            raise ValueError("EQUILIBRIUM_DT_HOURS must be > 0")

    if PRESSURE_SCENARIO not in ("sinus", "linear", "irregular", "csv"):
        raise ValueError("PRESSURE_SCENARIO must be one of: sinus, linear, irregular, csv")
    if PRESSURE_SCENARIO != "csv" and SCHEDULE_MODE == "direct":
        raise ValueError("SCHEDULE_MODE='direct' is only allowed for PRESSURE_SCENARIO='csv'")
    if PRESSURE_SCENARIO == "csv" and SCHEDULE_MODE not in ("stretch", "repeat", "direct"):
        raise ValueError("For csv: SCHEDULE_MODE must be stretch/repeat/direct")
    if PRESSURE_SCENARIO != "csv" and SCHEDULE_MODE not in ("stretch", "repeat"):
        raise ValueError("For non-csv: SCHEDULE_MODE must be stretch/repeat")
    if DT_HOURS <= 0:
        raise ValueError("DT_HOURS must be > 0")
    if N_CYCLES <= 0:
        raise ValueError("N_CYCLES must be > 0")

    # --- Initialization ---
    leaching_end_day: Optional[float] = None
    debrining_end_day: Optional[float] = None

    if USE_LEACHING:
        p_leach_end_mpa = LEACHING_END_FRACTION * P_LITHOSTATIC_MPA

        t_init_h, p_init_mpa = build_leaching_schedule_mpa(
            leaching_days=LEACHING_DAYS,
            leaching_dt_hours=LEACHING_DT_HOURS,
            p_start_mpa=P_LITHOSTATIC_MPA,
            p_end_mpa=p_leach_end_mpa,
            mode=LEACHING_MODE,
            stepped_n_steps=STEPPED_N_STEPS,
        )
        leaching_end_day = LEACHING_DAYS

        # Debrining plateau sampled at DT_HOURS to keep plot consistent
        if DEBRINING_DAYS > 0:
            t_deb_h = _sample_times_hours(DEBRINING_DAYS * DAY_H, DT_HOURS)
            p_deb_mpa = np.full_like(t_deb_h, p_leach_end_mpa)
            debrining_end_day = LEACHING_DAYS + DEBRINING_DAYS
        else:
            t_deb_h = np.array([], dtype=float)
            p_deb_mpa = np.array([], dtype=float)

        p_anchor_mpa = p_leach_end_mpa
        used_days = LEACHING_DAYS + DEBRINING_DAYS

    else:
        t_init_h, p_init_mpa = build_equilibrium_schedule_mpa(
            equilibrium_hours=EQUILIBRIUM_HOURS,
            equilibrium_dt_hours=EQUILIBRIUM_DT_HOURS,
            p_eq_mpa=P_EQUILIBRIUM_MPA,
        )
        t_deb_h = np.array([], dtype=float)
        p_deb_mpa = np.array([], dtype=float)
        p_anchor_mpa = P_EQUILIBRIUM_MPA
        used_days = EQUILIBRIUM_HOURS / 24.0

    # --- Operation fills remainder of TOTAL_DAYS ---
    operation_days = TOTAL_DAYS - used_days
    if operation_days <= 0:
        raise ValueError(
            f"Init (+debrining) uses {used_days:.2f} days, leaving operation_days={operation_days:.2f} (<=0).\n"
            f"Fix by increasing TOTAL_DAYS or decreasing init/debrining durations."
        )

    # Build operation schedule (relative time)
    if PRESSURE_SCENARIO == "sinus":
        if USE_LEACHING:
            p_mean_use = p_anchor_mpa + P_AMPLITUDE_MPA
        else:
            p_mean_use = P_MEAN_MPA
        t_op_h, p_op_mpa = build_sinus_operation_mpa(
            operation_days=operation_days,
            dt_hours=DT_HOURS,
            schedule_mode=SCHEDULE_MODE,
            n_cycles=N_CYCLES,
            p_mean_mpa=p_mean_use,
            p_ampl_mpa=P_AMPLITUDE_MPA,
        )

    elif PRESSURE_SCENARIO == "linear":
        if USE_LEACHING:
            pmin = p_anchor_mpa
            pmax = p_anchor_mpa + PRESSURE_SWING_MPA
        else:
            if P_MIN_MPA >= P_MAX_MPA:
                raise ValueError("For linear without leaching: P_MIN_MPA must be < P_MAX_MPA")
            pmin = P_MIN_MPA
            pmax = P_MAX_MPA
        t_op_h, p_op_mpa = build_linear_operation_mpa(
            operation_days=operation_days,
            dt_hours=DT_HOURS,
            schedule_mode=SCHEDULE_MODE,
            n_cycles=N_CYCLES,
            p_min_mpa=pmin,
            p_max_mpa=pmax,
        )

    elif PRESSURE_SCENARIO == "irregular":
        if len(IRREGULAR_WAYPOINTS_H) != len(IRREGULAR_PRESSURES_MPA):
            raise ValueError("IRREGULAR_WAYPOINTS_H and IRREGULAR_PRESSURES_MPA must have same length")

        base_p = list(IRREGULAR_PRESSURES_MPA)
        if USE_LEACHING:
            shift = p_anchor_mpa - float(min(base_p))
            base_p = [p + shift for p in base_p]

        t_op_h, p_op_mpa = build_irregular_operation_mpa(
            operation_days=operation_days,
            dt_hours=DT_HOURS,
            schedule_mode=SCHEDULE_MODE,
            n_cycles=N_CYCLES,
            base_waypoints_h=IRREGULAR_WAYPOINTS_H,
            base_pressures_mpa=base_p,
            smooth=IRREGULAR_SMOOTH,
        )

    elif PRESSURE_SCENARIO == "csv":
        if not os.path.isfile(CSV_FILE_PATH):
            raise FileNotFoundError(f"CSV_FILE_PATH not found: {CSV_FILE_PATH}")

        t_op_h, p_op_mpa = build_csv_operation_mpa(
            operation_days=operation_days,
            dt_hours=DT_HOURS,
            schedule_mode=SCHEDULE_MODE,
            n_cycles=N_CYCLES,
            csv_file=CSV_FILE_PATH,
            rescale=(CSV_RESCALE if not USE_LEACHING else False),
            rescale_min=CSV_RESCALE_MIN_MPA,
            rescale_max=CSV_RESCALE_MAX_MPA,
        )

        if USE_LEACHING and CSV_SHIFT_TO_LEACH_END:
            shift = p_anchor_mpa - float(np.min(p_op_mpa))
            p_op_mpa = p_op_mpa + shift

    else:
        raise ValueError("Unknown PRESSURE_SCENARIO")

    # Fade-in start of operation
    apply_fade_in(t_op_h, p_op_mpa, p_start_mpa=p_anchor_mpa, fade_in_hours=RAMP_UP_HOURS)

    # --- Stitch timeline: init + debrining + operation ---
    t_all_h: List[float] = []
    p_all_mpa: List[float] = []

    # init
    t_all_h.extend(t_init_h.tolist())
    p_all_mpa.extend(p_init_mpa.tolist())
    t_cursor_h = float(t_init_h[-1])

    # debrining
    if USE_LEACHING and DEBRINING_DAYS > 0:
        t_deb_shift = t_deb_h + t_cursor_h
        t_all_h.extend(t_deb_shift[1:].tolist())
        p_all_mpa.extend(p_deb_mpa[1:].tolist())
        t_cursor_h = float(t_deb_shift[-1])

    # operation
    t_op_shift = t_op_h + t_cursor_h
    t_all_h.extend(t_op_shift[1:].tolist())
    p_all_mpa.extend(p_op_mpa[1:].tolist())

    t_all_h_arr = np.asarray(t_all_h, dtype=float)
    p_all_mpa_arr = np.asarray(p_all_mpa, dtype=float)

    # force exact end time
    t_all_h_arr[-1] = TOTAL_DAYS * DAY_H

    return t_all_h_arr, p_all_mpa_arr, leaching_end_day, debrining_end_day


# ---------------------------
# Plot
# ---------------------------

def plot_schedule(t_h: np.ndarray, p_mpa: np.ndarray, leach_end_day: Optional[float], debr_end_day: Optional[float]) -> None:
    t_days = t_h / 24.0

    plt.figure(figsize=(12, 7))
    plt.plot(t_days, p_mpa, linewidth=2)

    # REQUIRED: dashed vertical lines for end of leaching and end of debrining
    if leach_end_day is not None:
        plt.axvline(leach_end_day, linestyle="--", linewidth=1)
    if debr_end_day is not None:
        plt.axvline(debr_end_day, linestyle="--", linewidth=1)

    plt.xlabel("Time (days)")
    plt.ylabel("Cavern pressure (MPa)")
    plt.title(f"Pressure scheme over {TOTAL_DAYS:.0f} days (init + debrining + operation)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if SAVE_PLOT:
        plt.savefig(SAVE_PATH, dpi=200)
        print(f"[SAVED] {SAVE_PATH}")

    if SHOW_PLOT:
        plt.show()
    else:
        plt.close()


def main() -> None:
    t_h, p_mpa, leach_end_day, debr_end_day = build_full_schedule_one_year()
    plot_schedule(t_h, p_mpa, leach_end_day, debr_end_day)


if __name__ == "__main__":
    main()
