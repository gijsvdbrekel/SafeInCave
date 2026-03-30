"""
Quick test script — lightweight version of Run.py for rapid NaN diagnosis.

Runs a SHORT simulation (leaching + a few operation cycles) to verify that
the Desai model handles a given cavern/pressure/material combination without
producing NaN.  Mirrors the Run.py configuration style so that if it works
here it will also work in a full Run.py run.

Usage:
    python3 quick_test.py                                          # uses defaults below
    python3 quick_test.py fastleached 1200 B SIC power_generation  # full override
    python3 quick_test.py fastleached 1200 B MD power_generation   # MD comparison
    python3 quick_test.py tubefailure 600 A                        # partial override
"""

import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import dolfinx as do
import os
import sys
import torch as to
import math
import numpy as np
import csv


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                         USER CONFIGURATION                                  ║
# ╠══════════════════════════════════════════════════════════════════════════════╣
# ║  Change these settings just like in Run.py. The script runs a short         ║
# ║  leaching + operation phase to check for NaN / convergence issues.          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── CAVERN SELECTION ────────────────────────────────────────────────────────
CAVERN_TYPE = "regular"     # regular, tilted, directcirculation, asymmetric,
                                 # reversedcirculation, fastleached, tubefailure
CAVERN_SIZE = 1200               # 600 or 1200

# ── MATERIAL MODEL ──────────────────────────────────────────────────────────
MATERIAL_SCENARIO = "A"          # "A" or "B"
USE_MUNSON_DAWSON = False         # False = SafeInCave (SIC), True = Munson-Dawson (MD)

# ── PRESSURE SCENARIO ───────────────────────────────────────────────────────
PRESSURE_SCENARIO = "csv"   # "constant", "industry", "transport",
                                          # "power_generation", "csv"

# ── LEACHING PHASE ──────────────────────────────────────────────────────────
LEACHING_DAYS = 91
LEACHING_DT_HOURS = 12
LEACHING_END_FRACTION = 0.30

# ── OPERATION PHASE (kept short for quick testing) ──────────────────────────
OPERATION_DAYS = 30              # Short! Just enough for 1-2 pressure cycles
N_CYCLES = 2                     # Number of cycles to fit in OPERATION_DAYS
dt_hours = 2.0                   # Base time step (hours)

# ── INDUSTRY SETTINGS ──────────────────────────────────────────────────────
P_AMPLITUDE_MPA = 5.0

# ── TRANSPORT SETTINGS ─────────────────────────────────────────────────────
P_HIGH_OFFSET_MPA = 5.0
P_LOW_OFFSET_MPA  = 1.5

# ── POWER GENERATION SETTINGS ──────────────────────────────────────────────
N_EVENTS = 2
P_BASE_OFFSET_MPA = 10.0
RECOVERY_TAU_HOURS = 200.0
P_MIN_MPA = 8.5

# ── VARIABLE TIME-STEPPING ─────────────────────────────────────────────────
USE_VARIABLE_DT = True
DT_FINE_HOURS   = 0.1
DT_COARSE_HOURS = 2.0
MAX_DP_MPA      = 0.2

# ── CSV SETTINGS ────────────────────────────────────────────────────────────
CSV_FILE_PATH = "drukprofiel_zoutcaverne_2035_8760u.csv"
RESCALE_PRESSURE = True
RESCALE_MIN_MPA = 6.0
RESCALE_MAX_MPA = 20.0

# ── FADE-IN / DEBRINING ────────────────────────────────────────────────────
RAMP_UP_HOURS = 336               # Shorter than Run.py (336h) for quick test
DEBRINING_DAYS = 0               # Skip debrining to save time

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                      END OF USER CONFIGURATION                               ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── CLI overrides ───────────────────────────────────────────────────────────
if len(sys.argv) > 1:
    CAVERN_TYPE = sys.argv[1]
if len(sys.argv) > 2:
    CAVERN_SIZE = int(sys.argv[2])
if len(sys.argv) > 3:
    MATERIAL_SCENARIO = sys.argv[3]
if len(sys.argv) > 4:
    arg4 = sys.argv[4].upper()
    if arg4 in ("SIC", "FALSE"):
        USE_MUNSON_DAWSON = False
    elif arg4 in ("MD", "TRUE"):
        USE_MUNSON_DAWSON = True
    else:
        PRESSURE_SCENARIO = sys.argv[4]
if len(sys.argv) > 5:
    PRESSURE_SCENARIO = sys.argv[5]

# ══════════════════════════════════════════════════════════════════════════════
# LOOKUP TABLES (from Run.py)
# ══════════════════════════════════════════════════════════════════════════════

Z_SURFACE = 660.0
salt_density = 2200
g = -9.81

Z_MAX_BY_CAVERN = {
    "regular600": 315.26, "tilted600": 345.67, "directcirculation600": 319.86,
    "asymmetric600": 338.89, "reversedcirculation600": 353.15,
    "regular1200": 393.21, "tilted1200": 430.78, "directcirculation1200": 402.21,
    "asymmetric1200": 422.76, "reversedcirculation1200": 445.06,
    "fastleached600": 378.19, "fastleached1200": 400.08,
    "tubefailure600": 420.90, "tubefailure1200": 444.53,
}
CAVERN_HEIGHT_BY_TYPE = {
    "regular600": 150.0, "tilted600": 160.0, "directcirculation600": 145.0,
    "asymmetric600": 155.0, "reversedcirculation600": 180.0,
    "regular1200": 200.0, "tilted1200": 215.0, "directcirculation1200": 195.0,
    "asymmetric1200": 210.0, "reversedcirculation1200": 240.0,
    "fastleached600": 170.0, "fastleached1200": 215.0,
    "tubefailure600": 182.0, "tubefailure1200": 230.0,
}
P_REF_BY_SIZE = {600: 17.5, 1200: 19.3}
GRID_FOLDERS = {
    "regular600": "cavern_regular_600_3D", "tilted600": "cavern_tilted_600_3D",
    "directcirculation600": "cavern_directcirculation_600_3D",
    "asymmetric600": "cavern_asymmetric_600_3D",
    "reversedcirculation600": "cavern_reversedcirculation_600_3D",
    "regular1200": "cavern_regular_1200_3D", "tilted1200": "cavern_tilted_1200_3D",
    "directcirculation1200": "cavern_directcirculation_1200_3D",
    "asymmetric1200": "cavern_asymmetric_1200_3D",
    "reversedcirculation1200": "cavern_reversedcirculation_1200_3D",
    "fastleached600": "cavern_fastleached_600_3D", "fastleached1200": "cavern_fastleached_1200_3D",
    "tubefailure600": "cavern_tubefailure_600_3D", "tubefailure1200": "cavern_tubefailure_1200_3D",
}

DAY_H = 24.0


# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS (copied from Run.py to stay in sync)
# ══════════════════════════════════════════════════════════════════════════════

def _sample_at_dt(tc, t_end=None):
    t_end = tc.t_final if t_end is None else t_end
    n_steps = int(math.floor(t_end / tc.dt))
    t_vals = [k * tc.dt for k in range(n_steps + 1)]
    if abs(t_vals[-1] - t_end) > 1e-12:
        t_vals.append(t_end)
    return t_vals


def _cardinal_segment(p0, p1, p2, p3, u, tension):
    m1 = (1 - tension) * 0.5 * (p2 - p0)
    m2 = (1 - tension) * 0.5 * (p3 - p1)
    u2 = u * u; u3 = u2 * u
    h00 = 2*u3 - 3*u2 + 1
    h10 = u3 - 2*u2 + u
    h01 = -2*u3 + 3*u2
    h11 = u3 - u2
    return h00*p1 + h10*m1 + h01*p2 + h11*m2


def _cardinal_interp(ts, ps, t, tension):
    if t <= ts[0]: return ps[0]
    if t >= ts[-1]: return ps[-1]
    i = np.searchsorted(ts, t) - 1
    t0, t1 = ts[i], ts[i+1]
    i_m1 = max(i-1, 0)
    i_p2 = min(i+2, len(ts)-1)
    u = (t - t0) / (t1 - t0)
    return _cardinal_segment(ps[i_m1], ps[i], ps[i+1], ps[i_p2], u, tension)


def _repeat_hours(times_h, days):
    times_h = list(map(float, times_h))
    out = []
    for d in range(int(days)):
        off = d * DAY_H
        for i, t in enumerate(times_h):
            if d > 0 and i == 0 and abs(t - 0.0) < 1e-12 and abs(times_h[-1] - DAY_H) < 1e-12:
                continue
            out.append(off + t)
    return out


def build_sinus_pressure_schedule(tc, *, p_mean, p_ampl, period_hours,
                                  phase_hours=0.0, clamp_min=None, clamp_max=None):
    period = period_hours * ut.hour
    phase = phase_hours * ut.hour
    n_steps = int(math.floor(tc.t_final / tc.dt))
    t_vals = [k * tc.dt for k in range(n_steps + 1)]
    if abs(t_vals[-1] - tc.t_final) > 1e-12:
        t_vals.append(tc.t_final)
    two_pi_over_T = (2.0 * math.pi / period) if period > 0.0 else 0.0
    p_vals = []
    for t in t_vals:
        p = p_mean if period <= 0.0 else p_mean + p_ampl * math.sin(two_pi_over_T * (t - phase))
        if clamp_min is not None: p = max(p, clamp_min)
        if clamp_max is not None: p = min(p, clamp_max)
        p_vals.append(p)
    return t_vals, p_vals


def build_sinus_schedule_multi(tc, *, p_mean, p_ampl, days, mode,
                                daily_period_hours=24.0, total_cycles=1,
                                clamp_min=0.0, clamp_max=None):
    total_hours = days * DAY_H
    if mode == "stretch":
        total_cycles = max(1, int(total_cycles))
        T_hours = total_hours / float(total_cycles)
    else:
        T_hours = daily_period_hours
    return build_sinus_pressure_schedule(
        tc, p_mean=p_mean, p_ampl=p_ampl,
        period_hours=T_hours, phase_hours=0.0,
        clamp_min=clamp_min, clamp_max=clamp_max)


def build_linear_schedule_multi(tc, times_h, pressures_MPa, *, days, mode,
                                resample_at_dt=True, total_cycles=None):
    times_h = list(map(float, times_h))
    pressures_MPa = list(map(float, pressures_MPa))
    total_h = float(days) * DAY_H
    if mode == "repeat":
        t_h = _repeat_hours(times_h, days)
        p_h = []
        for d in range(int(days)):
            start = 0 if d == 0 else 1
            p_h.extend(pressures_MPa[start:])
    else:
        if total_cycles is None: total_cycles = 1
        total_cycles = max(1, int(total_cycles))
        base_start, base_end = times_h[0], times_h[-1]
        base_duration = base_end - base_start
        cycle_duration = total_h / float(total_cycles)
        scale = cycle_duration / base_duration
        t_h, p_h = [], []
        for k in range(total_cycles):
            offset = k * cycle_duration
            for i, t in enumerate(times_h):
                if k > 0 and i == 0: continue
                t_h.append(offset + (t - base_start) * scale)
                p_h.append(pressures_MPa[i])
    if t_h[0] > 0.0: t_h.insert(0, 0.0); p_h.insert(0, p_h[0])
    if t_h[-1] < total_h: t_h.append(total_h); p_h.append(p_h[-1])
    knots_t = np.array([h * ut.hour for h in t_h], dtype=float)
    knots_p = np.array([p * ut.MPa for p in p_h], dtype=float)
    if resample_at_dt:
        t_vals = _sample_at_dt(tc)
        p_vals = np.interp(t_vals, knots_t, knots_p).tolist()
    else:
        t_vals = knots_t.tolist(); p_vals = knots_p.tolist()
    return t_vals, p_vals


def build_power_generation_schedule(tc, *, p_base_pa, n_events, operation_days,
                                    recovery_tau_hours=48.0, p_min_pa=None, seed=42):
    t_vals_s = _sample_at_dt(tc)
    t_h = np.array(t_vals_s) / ut.hour
    p_base_mpa = p_base_pa / ut.MPa
    p_min_mpa = p_min_pa / ut.MPa if p_min_pa is not None else None
    p_mpa = np.full(len(t_h), p_base_mpa)
    rng = np.random.RandomState(seed)
    event_centers_days = np.linspace(1.0, operation_days - 1.0, max(1, n_events))
    event_centers_days = event_centers_days + rng.uniform(-0.8, 0.8, size=n_events)
    tau = max(0.1, float(recovery_tau_hours))
    for day_c in event_centers_days:
        t_start_h = day_c * 24.0
        duration = rng.uniform(2.0, 5.0)
        depth = rng.uniform(6.5, 12.5)
        for i, t in enumerate(t_h):
            if t < t_start_h: continue
            dt_ev = t - t_start_h
            if dt_ev < 0.5:       drop = depth * (dt_ev / 0.5)
            elif dt_ev < 0.5 + duration: drop = depth
            else:
                drop = depth * math.exp(-(dt_ev - 0.5 - duration) / tau)
                if drop < 0.05: break
            p_mpa[i] = min(p_mpa[i], p_base_mpa - drop)
    if p_min_mpa is not None:
        p_mpa = np.maximum(p_mpa, p_min_mpa)
    return t_vals_s, (p_mpa * ut.MPa).tolist()


def _parse_float_auto(s):
    s = s.strip()
    if not s: return np.nan
    if "," in s and "." not in s: s = s.replace(",", ".")
    if "," in s and "." in s: s = s.replace(".", "").replace(",", ".")
    try: return float(s)
    except Exception: return np.nan


def read_pressure_csv(csv_file):
    if not os.path.isfile(csv_file):
        raise FileNotFoundError(f"CSV not found: {csv_file}")
    with open(csv_file, "r", newline="", encoding="utf-8") as f:
        sample = f.read(2048); f.seek(0)
        try: delim = csv.Sniffer().sniff(sample, delimiters=";,\t").delimiter
        except Exception: delim = ";"
    rows = []
    with open(csv_file, "r", newline="", encoding="utf-8") as f:
        for r in csv.reader(f, delimiter=delim):
            if len(r) > 0: rows.append(r)
    if len(rows) < 2: raise ValueError("CSV has too few rows")
    header_low = [c.strip().lower() for c in rows[0]]
    data = rows[1:]
    idx_mpa = idx_bar = None
    for i, h in enumerate(header_low):
        if h == "druk_mpa" or h.endswith("druk_mpa"): idx_mpa = i
        if h == "druk_bar" or h.endswith("druk_bar"): idx_bar = i
    pressures_mpa = []
    if idx_mpa is not None:
        for r in data:
            if idx_mpa < len(r): pressures_mpa.append(_parse_float_auto(r[idx_mpa]))
    elif idx_bar is not None:
        for r in data:
            if idx_bar < len(r): pressures_mpa.append(_parse_float_auto(r[idx_bar]) / 10.0)
    else:
        ncols = len(header_low); best_i, best_count = None, -1
        for i in range(ncols):
            vals = [_parse_float_auto(r[i]) for r in data if i < len(r)]
            count = int(np.sum(np.isfinite(vals)))
            if count > best_count: best_count, best_i = count, i
        for r in data:
            if best_i < len(r): pressures_mpa.append(_parse_float_auto(r[best_i]))
    pressures_mpa = np.asarray(pressures_mpa, dtype=float)
    pressures_mpa = pressures_mpa[np.isfinite(pressures_mpa)]
    return pressures_mpa


def rescale_pressure_profile(pressures_mpa, new_min, new_max):
    old_min, old_max = pressures_mpa.min(), pressures_mpa.max()
    if old_max - old_min < 1e-9:
        return np.full_like(pressures_mpa, (new_min + new_max) / 2.0)
    return new_min + (pressures_mpa - old_min) / (old_max - old_min) * (new_max - new_min)


def build_csv_pressure_schedule(tc, csv_file, *, days, mode, total_cycles=1,
                                rescale=False, rescale_min=None, rescale_max=None,
                                resample_at_dt=True):
    pressures_mpa = read_pressure_csv(csv_file)
    csv_hours = int(pressures_mpa.size)
    print(f"[CSV] Loaded '{os.path.basename(csv_file)}': {csv_hours} hourly values, "
          f"range [{pressures_mpa.min():.2f}, {pressures_mpa.max():.2f}] MPa")
    if rescale and rescale_min is not None and rescale_max is not None:
        pressures_mpa = rescale_pressure_profile(pressures_mpa, rescale_min, rescale_max)
        print(f"[CSV] Rescaled to [{pressures_mpa.min():.2f}, {pressures_mpa.max():.2f}] MPa")
    total_hours = float(days) * 24.0
    if mode == "stretch":
        total_cycles = max(1, int(total_cycles))
        cycle_duration_h = total_hours / float(total_cycles)
        scale = cycle_duration_h / float(csv_hours)
        times_list, pres_list = [], []
        for k in range(total_cycles):
            off = k * cycle_duration_h
            for i in range(csv_hours):
                if k > 0 and i == 0: continue
                times_list.append(off + i * scale)
                pres_list.append(pressures_mpa[i])
        times_hours = np.asarray(times_list, dtype=float)
        pressures_MPa = np.asarray(pres_list, dtype=float)
    else:
        raise ValueError("quick_test only supports mode='stretch' for CSV")
    times_s = times_hours * 3600.0
    if times_s[0] > 0.0:
        times_s = np.insert(times_s, 0, 0.0)
        pressures_MPa = np.insert(pressures_MPa, 0, pressures_MPa[0])
    if times_s[-1] < tc.t_final:
        times_s = np.append(times_s, tc.t_final)
        pressures_MPa = np.append(pressures_MPa, pressures_MPa[-1])
    if resample_at_dt:
        t_vals = _sample_at_dt(tc)
        p_vals_MPa = np.interp(t_vals, times_s, pressures_MPa)
    else:
        t_vals = times_s.tolist()
        p_vals_MPa = pressures_MPa
    p_vals = [float(p) * ut.MPa for p in p_vals_MPa]
    return t_vals, p_vals


def apply_fade_in(t_pressure, p_pressure, *, p_start_pa, fade_in_hours):
    if fade_in_hours <= 0.0: return
    fade_in_s = fade_in_hours * 3600.0
    for i in range(len(t_pressure)):
        t = t_pressure[i]
        if t >= fade_in_s: break
        alpha = 0.5 * (1.0 - math.cos(math.pi * t / fade_in_s))
        p_pressure[i] = (1.0 - alpha) * p_start_pa + alpha * p_pressure[i]


def prepend_debrining(t_pressure, p_pressure, *, p_leach_end_pa, debrining_days):
    debrining_s = debrining_days * 24.0 * 3600.0
    if debrining_s <= 0.0: return t_pressure, p_pressure
    t_pre = [0.0, debrining_s]
    p_pre = [p_leach_end_pa, p_leach_end_pa]
    t_shifted = [t + debrining_s for t in t_pressure[1:]]
    p_shifted = list(p_pressure[1:])
    return t_pre + t_shifted, p_pre + p_shifted


class TimeControllerFromList:
    def __init__(self, time_list_seconds):
        self.time_list = np.asarray(time_list_seconds, dtype=float)
        self.t_initial = float(self.time_list[0])
        self.t_final = float(self.time_list[-1])
        self.t = float(self.time_list[0])
        self.step_counter = 0
        self.dt = float(self.time_list[1] - self.time_list[0])
        self.time_unit = "hour"
        self.time_conversion = ut.hour

    def keep_looping(self):
        return self.step_counter < (self.time_list.size - 1)

    def advance_time(self):
        self.step_counter += 1
        t_prev = self.t
        self.t = float(self.time_list[self.step_counter])
        self.dt = self.t - t_prev


def build_time_list_by_dp_limit(t_final_s, p_of_t, *, dt_min_s, dt_max_s, dp_max_pa):
    t = 0.0; times = [0.0]; p_prev = float(p_of_t(0.0))
    max_steps = int(np.ceil(t_final_s / dt_min_s)) + 50
    for _ in range(max_steps):
        if t >= t_final_s - 1e-12: break
        dt = dt_max_s
        while True:
            t_try = min(t + dt, t_final_s)
            p_try = float(p_of_t(t_try))
            if abs(p_try - p_prev) <= dp_max_pa or dt <= dt_min_s + 1e-12:
                t = t_try; p_prev = p_try; times.append(t); break
            dt *= 0.5
            if dt < dt_min_s: dt = dt_min_s
    if abs(times[-1] - t_final_s) > 1e-9: times.append(t_final_s)
    return times


# ══════════════════════════════════════════════════════════════════════════════
# NaN CHECKER
# ══════════════════════════════════════════════════════════════════════════════

def check_nan(label, *tensors):
    """Return True if any tensor contains NaN or Inf."""
    found = False
    for i, t in enumerate(tensors):
        if t is None: continue
        arr = t if isinstance(t, to.Tensor) else to.tensor(t)
        n_nan = int(to.isnan(arr).sum().item())
        n_inf = int(to.isinf(arr).sum().item())
        if n_nan > 0 or n_inf > 0:
            print(f"  *** {label}[{i}]: {n_nan} NaN, {n_inf} Inf / {arr.numel()} values")
            found = True
    return found


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    cavern_key = f"{CAVERN_TYPE}{CAVERN_SIZE}"
    if cavern_key not in GRID_FOLDERS:
        print(f"Unknown cavern: {cavern_key}. Valid: {list(GRID_FOLDERS.keys())}")
        sys.exit(1)

    z_max = Z_MAX_BY_CAVERN[cavern_key]
    h_cav = CAVERN_HEIGHT_BY_TYPE[cavern_key]
    z_center = z_max - h_cav / 2
    p_ref = P_REF_BY_SIZE[CAVERN_SIZE] * ut.MPa

    depth = Z_SURFACE - z_center
    p_lithostatic = p_ref + salt_density * abs(g) * depth
    p_lithostatic_mpa = p_lithostatic / ut.MPa
    p_leach_end_mpa = LEACHING_END_FRACTION * p_lithostatic_mpa
    p_leach_end = p_leach_end_mpa * ut.MPa

    sec_per_year = 365.25 * 24 * 3600
    gas_density = 0.089
    side_burden = p_ref
    over_burden = p_ref
    g_vec = [0.0, 0.0, g]

    model_tag = "MD" if USE_MUNSON_DAWSON else "SIC"

    print("=" * 70)
    print(f"QUICK TEST — {CAVERN_TYPE} {CAVERN_SIZE}k | Scenario {MATERIAL_SCENARIO} {model_tag} | {PRESSURE_SCENARIO}")
    print("=" * 70)
    print(f"  z_max={z_max:.1f}m  z_center={z_center:.1f}m")
    print(f"  P_litho={p_lithostatic_mpa:.2f} MPa  P_leach_end={p_leach_end_mpa:.2f} MPa")
    print(f"  Operation: {OPERATION_DAYS} days, {N_CYCLES} cycles, dt={dt_hours}h")
    print("=" * 70)

    # ── Load grid ────────────────────────────────────────────────────────────
    grid_path = os.path.join("..", "..", "..", "..", "grids", GRID_FOLDERS[cavern_key])
    grid = sf.GridHandlerGMSH("geom", grid_path)
    print(f"[GRID] {grid.mesh.topology.index_map(3).size_local} elements")

    # ── Momentum equation ────────────────────────────────────────────────────
    mom_eq = sf.LinearMomentum(grid, theta=0.5)
    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("cg")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-10, max_it=100)
    mom_eq.set_solver(mom_solver)

    # ── Material ─────────────────────────────────────────────────────────────
    mat = sf.Material(mom_eq.n_elems)
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    E0 = 20.425 * ut.GPa * to.ones(mom_eq.n_elems)
    nu0 = 0.25 * to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")
    mat.add_to_elastic(spring_0)

    if MATERIAL_SCENARIO == "A":
        if not USE_MUNSON_DAWSON:
            # ── SIC Scenario A: CCC Zuidwending ──
            eta = 2.5e5 * ut.GPa * to.ones(mom_eq.n_elems)
            E1  = 42.0 * ut.GPa * to.ones(mom_eq.n_elems)
            nu1 = 0.32 * to.ones(mom_eq.n_elems)
            kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")
            ndc = 4.6
            A_dc = (40.0 * (1e-6)**ndc / sec_per_year) * to.ones(mom_eq.n_elems, dtype=to.float64)
            Q_dc = (6495.0 * 8.32) * to.ones(mom_eq.n_elems)
            n_dc = ndc * to.ones(mom_eq.n_elems)
            creep_0 = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation")
            mat.add_to_non_elastic(kelvin)
            mat.add_to_non_elastic(creep_0)
        else:
            # ── MD Scenario A ──
            nmd   = 4.99
            A_md  = (18.31 * (1e-6)**nmd / sec_per_year) * to.ones(mom_eq.n_elems, dtype=to.float64)
            Q_md  = (6356.0 * 8.32) * to.ones(mom_eq.n_elems)
            n_md  = nmd  * to.ones(mom_eq.n_elems)
            K0_md = 7.0e-7 * to.ones(mom_eq.n_elems)
            c_md  = 9.02e-3 * to.ones(mom_eq.n_elems)
            m_md  = 3.0    * to.ones(mom_eq.n_elems)
            aw_md = -13.2  * to.ones(mom_eq.n_elems)
            bw_md = -7.738 * to.ones(mom_eq.n_elems)
            d_md  = 0.58   * to.ones(mom_eq.n_elems)
            mu_md = E0 / (2.0 * (1.0 + nu0))
            md = sf.MunsonDawsonCreep(A=A_md, Q=Q_md, n=n_md,
                                      K0=K0_md, c=c_md, m=m_md,
                                      alpha_w=aw_md, beta_w=bw_md, delta=d_md,
                                      mu=mu_md, name="munson_dawson")
            mat.add_to_non_elastic(md)
    elif MATERIAL_SCENARIO == "B":
        if not USE_MUNSON_DAWSON:
            # ── SIC Scenario B: calibrated against cyclic triaxial data ──
            eta = 1.0e15 * to.ones(mom_eq.n_elems)
            E1  = 1.0e11 * to.ones(mom_eq.n_elems)
            nu1 = 0.25 * to.ones(mom_eq.n_elems)
            kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")
            ndc = 5.6897
            A_dc = (25.92 * (1e-6)**ndc / sec_per_year) * to.ones(mom_eq.n_elems, dtype=to.float64)
            Q_dc = (6495.0 * 8.32) * to.ones(mom_eq.n_elems)
            n_dc = ndc * to.ones(mom_eq.n_elems)
            creep_0 = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation")
            mat.add_to_non_elastic(kelvin)
            mat.add_to_non_elastic(creep_0)
        else:
            # ── MD Scenario B ──
            nmd   = 5.6897
            A_md  = (17.28 * (1e-6)**nmd / sec_per_year) * to.ones(mom_eq.n_elems, dtype=to.float64)
            Q_md  = (6495.0 * 8.32) * to.ones(mom_eq.n_elems)
            n_md  = nmd  * to.ones(mom_eq.n_elems)
            K0_md = 2253.87  * to.ones(mom_eq.n_elems)
            c_md  = 9.02e-3  * to.ones(mom_eq.n_elems)
            m_md  = 2.466    * to.ones(mom_eq.n_elems)
            aw_md = 179.70   * to.ones(mom_eq.n_elems)
            bw_md = 60.00    * to.ones(mom_eq.n_elems)
            d_md  = 299.95   * to.ones(mom_eq.n_elems)
            mu_md = E0 / (2.0 * (1.0 + nu0))
            md = sf.MunsonDawsonCreep(A=A_md, Q=Q_md, n=n_md,
                                      K0=K0_md, c=c_md, m=m_md,
                                      alpha_w=aw_md, beta_w=bw_md, delta=d_md,
                                      mu=mu_md, name="munson_dawson")
            mat.add_to_non_elastic(md)
    else:
        print(f"Unknown MATERIAL_SCENARIO: {MATERIAL_SCENARIO}")
        sys.exit(1)

    # Pressure-solution creep (always added, same for SIC and MD)
    A_ps = (14176.0 * 1e-9 / 1e6 / sec_per_year) * to.ones(mom_eq.n_elems)
    d_ps = 5.25e-3 * to.ones(mom_eq.n_elems)
    Q_ps = (3252.0 * 8.32) * to.ones(mom_eq.n_elems)
    creep_pressure = sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "creep_pressure")
    mat.add_to_non_elastic(creep_pressure)
    mom_eq.set_material(mat)
    mom_eq.build_body_force(g_vec)

    T0_field = 298 * to.ones(mom_eq.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # ══════════════════════════════════════════════════════════════════════════
    # LEACHING PHASE
    # ══════════════════════════════════════════════════════════════════════════
    print(f"\n[LEACHING] {LEACHING_DAYS} days, dt={LEACHING_DT_HOURS}h, "
          f"{p_lithostatic_mpa:.1f} -> {p_leach_end_mpa:.1f} MPa")

    tc_leach = sf.TimeController(
        dt=LEACHING_DT_HOURS, initial_time=0.0,
        final_time=LEACHING_DAYS * 24.0, time_unit="hour")

    n_steps_leach = int(math.floor(tc_leach.t_final / tc_leach.dt))
    t_leach = [k * tc_leach.dt for k in range(n_steps_leach + 1)]
    if abs(t_leach[-1] - tc_leach.t_final) > 1e-12:
        t_leach.append(tc_leach.t_final)
    p_leach = [p_lithostatic + (p_leach_end - p_lithostatic) * (t / tc_leach.t_final)
               for t in t_leach]

    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_leach.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_leach.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_leach.t_final])
    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0,
                              [side_burden, side_burden], [0.0, tc_leach.t_final], g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0,
                               [side_burden, side_burden], [0.0, tc_leach.t_final], g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                             [over_burden, over_burden], [0.0, tc_leach.t_final], g=g_vec[2])
    bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, z_max,
                                p_leach, t_leach, g=g_vec[2])

    bc_handler = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
        bc_handler.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_handler)

    output_folder = os.path.join("output", f"quick_test_{cavern_key}_S{MATERIAL_SCENARIO}_{model_tag}_{PRESSURE_SCENARIO}")
    os.makedirs(output_folder, exist_ok=True)
    output_leach = sf.SaveFields(mom_eq)
    output_leach.set_output_folder(os.path.join(output_folder, "leaching"))
    output_leach.add_output_field("sig", "Stress (Pa)")

    sim_leach = sf.Simulator_M(mom_eq, tc_leach, [output_leach], True)
    sim_leach.run()
    print("[LEACHING] Complete.")

    # ── Post-leaching stress diagnostics ─────────────────────────────────────
    stress_np = mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3))
    stress_to = ut.numpy2torch(stress_np)

    s = -stress_to / ut.MPa
    I1 = s[:, 0, 0] + s[:, 1, 1] + s[:, 2, 2]
    I2 = (s[:, 0, 0]*s[:, 1, 1] + s[:, 1, 1]*s[:, 2, 2] + s[:, 0, 0]*s[:, 2, 2]
          - s[:, 0, 1]**2 - s[:, 0, 2]**2 - s[:, 1, 2]**2)
    J2 = (1/3)*I1**2 - I2
    I1_star = I1 + 5.0

    print(f"\n[POST-LEACHING STRESS]")
    print(f"  I1     : min={I1.min().item():.3f}  max={I1.max().item():.3f} MPa")
    print(f"  I1_star: min={I1_star.min().item():.3f}  max={I1_star.max().item():.3f} MPa")
    print(f"  J2     : min={J2.min().item():.6f}  max={J2.max().item():.3f} MPa^2")
    print(f"  Elements with I1 < 5 MPa: {int((I1 < 5).sum())}/{mom_eq.n_elems}")
    print(f"  Elements with I1 < 0 MPa: {int((I1 < 0).sum())}/{mom_eq.n_elems}")

    # ══════════════════════════════════════════════════════════════════════════
    # DESAI INITIALIZATION (SIC only — MD has no Desai element)
    # ══════════════════════════════════════════════════════════════════════════
    desai = None
    if not USE_MUNSON_DAWSON:
        print(f"\n[DESAI] Initializing (Scenario {MATERIAL_SCENARIO})...")

        if MATERIAL_SCENARIO == "A":
            mu_1    = 6.89e-12 * to.ones(mom_eq.n_elems)
            N_1     = 3.0 * to.ones(mom_eq.n_elems)
            a_1     = 1.80e-5 * to.ones(mom_eq.n_elems)
            eta_vp  = 0.82 * to.ones(mom_eq.n_elems)
            alpha_0_val = 2.0e-3
        elif MATERIAL_SCENARIO == "B":
            mu_1    = 1.016e-15 * to.ones(mom_eq.n_elems)
            N_1     = 4.2515 * to.ones(mom_eq.n_elems)
            a_1     = 1.101e-6 * to.ones(mom_eq.n_elems)
            eta_vp  = 1.7902 * to.ones(mom_eq.n_elems)
            alpha_0_val = 1.781e-3

        alpha_0 = alpha_0_val * to.ones(mom_eq.n_elems)
        n_desai = 3.0 * to.ones(mom_eq.n_elems)
        beta_1  = 0.0048 * to.ones(mom_eq.n_elems)
        beta    = 0.995 * to.ones(mom_eq.n_elems)
        m_desai = -0.5 * to.ones(mom_eq.n_elems)
        gamma   = 0.095 * to.ones(mom_eq.n_elems)
        sigma_t = 5.0 * to.ones(mom_eq.n_elems)

        desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta_vp, n_desai, beta_1, beta,
                                     m_desai, gamma, sigma_t, alpha_0, "desai")
        desai.compute_initial_hardening(stress_to, Fvp_0=0.0)

        has_nan = check_nan("alpha_0", desai.alpha_0)
        has_nan |= check_nan("alpha", desai.alpha)
        has_nan |= check_nan("Fvp", desai.Fvp)

        n_disabled = len(desai.ind_desai_disabled) if hasattr(desai, 'ind_desai_disabled') else 0
        print(f"  alpha_0 : min={desai.alpha_0.min().item():.3e}  max={desai.alpha_0.max().item():.3e}")
        print(f"  Fvp     : min={desai.Fvp.min().item():.3e}  max={desai.Fvp.max().item():.3e}")
        print(f"  Yielding (Fvp>0): {int((desai.Fvp > 0).sum())}/{mom_eq.n_elems}")
        print(f"  Disabled (alpha_0 clamped): {n_disabled}/{mom_eq.n_elems}")

        if has_nan:
            print("\n*** FAIL: NaN after Desai init ***")
            sys.exit(1)

        mat.add_to_non_elastic(desai)
        mom_eq.set_material(mat)
    else:
        print(f"\n[MD] Munson-Dawson model — no Desai initialization needed.")

    # ══════════════════════════════════════════════════════════════════════════
    # OPERATION PHASE — build pressure schedule (same logic as Run.py)
    # ══════════════════════════════════════════════════════════════════════════
    tc_cycling = sf.TimeController(
        dt=dt_hours, initial_time=0.0,
        final_time=OPERATION_DAYS * 24.0, time_unit="hour")

    if PRESSURE_SCENARIO == "constant":
        t_pressure = [0.0, tc_cycling.t_final]
        p_pressure = [p_leach_end, p_leach_end]

    elif PRESSURE_SCENARIO == "industry":
        p_mean = (p_leach_end_mpa + P_AMPLITUDE_MPA) * ut.MPa
        p_ampl = P_AMPLITUDE_MPA * ut.MPa
        t_pressure, p_pressure = build_sinus_schedule_multi(
            tc_cycling, p_mean=p_mean, p_ampl=p_ampl,
            days=OPERATION_DAYS, mode="stretch",
            total_cycles=N_CYCLES)

    elif PRESSURE_SCENARIO == "transport":
        p_high = p_leach_end_mpa + P_HIGH_OFFSET_MPA
        p_low  = p_leach_end_mpa + P_LOW_OFFSET_MPA
        base_times_h       = [0.0, 8.0, 12.0, 28.0, 32.0, 48.0]
        base_pressures_MPa = [p_high, p_high, p_low, p_low, p_high, p_high]
        t_pressure, p_pressure = build_linear_schedule_multi(
            tc_cycling, base_times_h, base_pressures_MPa,
            days=OPERATION_DAYS, mode="stretch",
            total_cycles=N_CYCLES)

    elif PRESSURE_SCENARIO == "power_generation":
        p_base_pa = (p_leach_end_mpa + P_BASE_OFFSET_MPA) * ut.MPa
        t_pressure, p_pressure = build_power_generation_schedule(
            tc_cycling, p_base_pa=p_base_pa,
            n_events=N_EVENTS, operation_days=OPERATION_DAYS,
            recovery_tau_hours=RECOVERY_TAU_HOURS,
            p_min_pa=P_MIN_MPA * ut.MPa)

    elif PRESSURE_SCENARIO == "csv":
        t_pressure, p_pressure = build_csv_pressure_schedule(
            tc_cycling, csv_file=CSV_FILE_PATH,
            days=OPERATION_DAYS, mode="stretch",
            total_cycles=N_CYCLES,
            rescale=RESCALE_PRESSURE,
            rescale_min=RESCALE_MIN_MPA,
            rescale_max=RESCALE_MAX_MPA)
    else:
        print(f"Unknown PRESSURE_SCENARIO: {PRESSURE_SCENARIO}")
        sys.exit(1)

    # Apply fade-in
    if RAMP_UP_HOURS > 0:
        apply_fade_in(t_pressure, p_pressure,
                      p_start_pa=p_leach_end, fade_in_hours=RAMP_UP_HOURS)

    # Prepend debrining
    extra_hours = 0.0
    if DEBRINING_DAYS > 0:
        extra_hours = DEBRINING_DAYS * 24.0
        t_pressure, p_pressure = prepend_debrining(
            t_pressure, p_pressure,
            p_leach_end_pa=p_leach_end, debrining_days=DEBRINING_DAYS)

    total_operation_hours = OPERATION_DAYS * 24.0 + extra_hours

    # Variable dt for power_generation / csv
    if PRESSURE_SCENARIO in ("power_generation", "csv") and USE_VARIABLE_DT:
        t_arr = np.array(t_pressure, dtype=float)
        p_arr = np.array(p_pressure, dtype=float)
        p_of_t = lambda t: float(np.interp(t, t_arr, p_arr))
        time_list = build_time_list_by_dp_limit(
            total_operation_hours * ut.hour, p_of_t,
            dt_min_s=DT_FINE_HOURS * ut.hour,
            dt_max_s=DT_COARSE_HOURS * ut.hour,
            dp_max_pa=MAX_DP_MPA * ut.MPa)
        tc_operation = TimeControllerFromList(time_list)
        n_steps = len(time_list) - 1
        print(f"[VARIABLE DT] {n_steps} steps")
    else:
        tc_operation = sf.TimeController(
            dt=dt_hours, initial_time=0.0,
            final_time=total_operation_hours, time_unit="hour")

    p_min_mpa = min(p_pressure) / ut.MPa
    p_max_mpa = max(p_pressure) / ut.MPa
    print(f"\n[OPERATION] P range: [{p_min_mpa:.2f}, {p_max_mpa:.2f}] MPa, "
          f"duration: {total_operation_hours:.0f}h")

    # ── Operation BCs ────────────────────────────────────────────────────────
    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0,
                              [side_burden, side_burden], [0.0, tc_operation.t_final], g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0,
                               [side_burden, side_burden], [0.0, tc_operation.t_final], g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                             [over_burden, over_burden], [0.0, tc_operation.t_final], g=g_vec[2])
    bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, z_max,
                                p_pressure, t_pressure, g=g_vec[2])

    bc_handler_op = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
        bc_handler_op.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_handler_op)

    output_op = sf.SaveFields(mom_eq)
    output_op.set_output_folder(os.path.join(output_folder, "operation"))
    output_op.add_output_field("sig", "Stress (Pa)")

    sim_op = sf.Simulator_M(mom_eq, tc_operation, [output_op], True)
    sim_op.run()

    # ── Final NaN check ──────────────────────────────────────────────────────
    stress_final = ut.numpy2torch(mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3)))
    has_nan = check_nan("final_stress", stress_final)
    if desai is not None:
        has_nan |= check_nan("final_alpha", desai.alpha)
        has_nan |= check_nan("final_Fvp", desai.Fvp)

    print("\n" + "=" * 70)
    if has_nan:
        print(f"*** FAIL: NaN detected after operation phase ***")
        print(f"  {CAVERN_TYPE} {CAVERN_SIZE}k | S{MATERIAL_SCENARIO} {model_tag} | {PRESSURE_SCENARIO}")
        sys.exit(1)
    else:
        print(f"*** PASS: No NaN ***")
        print(f"  {CAVERN_TYPE} {CAVERN_SIZE}k | S{MATERIAL_SCENARIO} {model_tag} | {PRESSURE_SCENARIO}")
        if desai is not None:
            print(f"  alpha: min={desai.alpha.min().item():.3e}  max={desai.alpha.max().item():.3e}")
            print(f"  Fvp  : min={desai.Fvp.min().item():.3e}  max={desai.Fvp.max().item():.3e}")
    print("=" * 70)


if __name__ == "__main__":
    main()
