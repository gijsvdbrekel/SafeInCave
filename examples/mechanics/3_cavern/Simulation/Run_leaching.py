import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import dolfinx as do
import os
import torch as to
import json
import math
import numpy as np
import csv


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                           USER CONFIGURATION                                  ║
# ╠══════════════════════════════════════════════════════════════════════════════╣
# ║  Modify the settings below to configure your simulation.                      ║
# ║  This script includes a LEACHING PHASE before operation.                      ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── CAVERN SELECTION ───────────────────────────────────────────────────────────
# CAVERN_SHAPE: Choose one of:
#   "regular"      - Standard cylindrical cavern
#   "tilted"       - Tilted/inclined cavern
#   "teardrop"     - Reversed teardrop-shaped cavern
#   "asymmetric"   - Asymmetric cavern geometry
#   "irregular"    - Irregular/complex shape (bulb)
#   "multichamber" - Multi-chamber cavern

CAVERN_SHAPE = "regular"

# CAVERN_SIZE: Choose one of:
#   600  - 600,000 m³ volume
#   1200 - 1,200,000 m³ volume

CAVERN_SIZE = 600

# ── LEACHING PHASE SETTINGS ──────────────────────────────────────────────────────
# LEACHING_MODE: How pressure decreases during leaching:
#   "linear"  - Linear decrease from lithostatic to operational pressure
#   "stepped" - Stepped decrease with plateaus (more realistic)

LEACHING_MODE = "stepped"

# LEACHING_DAYS: Duration of leaching phase in days (default ~0.25 year = 91 days)
LEACHING_DAYS = 91

# LEACHING_DT_HOURS: Time step during leaching (coarser than operation)
LEACHING_DT_HOURS = 12

# STEPPED_N_STEPS: Number of pressure steps for "stepped" mode (only used if LEACHING_MODE = "stepped")
#   Each step holds pressure constant for LEACHING_DAYS / STEPPED_N_STEPS days
STEPPED_N_STEPS = 6

# LEACHING_END_FRACTION: Fraction of lithostatic pressure at end of leaching
#   The leaching phase ends when pressure reaches this fraction of the initial lithostatic pressure.
#   All operational pressure schemes will start from this pressure (p_leach_end).
#   Example: 0.30 means leaching ends at 30% of lithostatic pressure.
LEACHING_END_FRACTION = 0.40

# ── PRESSURE SCENARIO ──────────────────────────────────────────────────────────
# PRESSURE_SCENARIO: Choose one of:
#   "sinus"     - Sinusoidal pressure variation
#   "linear"    - Piecewise linear pressure profile
#   "irregular" - Irregular/spline-smoothed pressure profile
#   "csv"       - Load pressure profile from CSV file (e.g., real operational data)

PRESSURE_SCENARIO = "csv"

# ── SINUS SETTINGS (only used when PRESSURE_SCENARIO = "sinus") ────────────────
# P_AMPLITUDE_MPA: Amplitude (MPa) - half the peak-to-peak range
#   The operational pressure will oscillate around a mean that is derived from leaching:
#   - P_MIN (operational) = p_leach_end (= LEACHING_END_FRACTION * p_lithostatic)
#   - P_MEAN (operational) = p_leach_end + P_AMPLITUDE_MPA
#   - P_MAX (operational) = p_leach_end + 2 * P_AMPLITUDE_MPA
#   Example: If p_leach_end=8 MPa and amplitude=6.5 MPa → range [8, 21] MPa

P_AMPLITUDE_MPA = 6.5

# ── LINEAR SETTINGS (only used when PRESSURE_SCENARIO = "linear") ──────────────
# PRESSURE_SWING_MPA: Pressure swing (MPa) - difference between max and min
#   The operational pressure will cycle between:
#   - P_MIN (operational) = p_leach_end (= LEACHING_END_FRACTION * p_lithostatic)
#   - P_MAX (operational) = p_leach_end + PRESSURE_SWING_MPA
#   Example: If p_leach_end=8 MPa and swing=6.4 MPa → range [8, 14.4] MPa

PRESSURE_SWING_MPA = 10

# ── SCHEDULE SETTINGS ──────────────────────────────────────────────────────────
# SCHEDULE_MODE: How to distribute cycles over the simulation period:
#   "stretch" - N_CYCLES spread evenly over OPERATION_DAYS
#   "repeat"  - Daily pattern repeated each day (N_CYCLES ignored)
#   "direct"  - (CSV only) Use hourly CSV values directly, periodic over CSV length

SCHEDULE_MODE = "stretch"

# OPERATION_DAYS: Total simulation duration in days (operation phase only)
OPERATION_DAYS = 365 * 3

# N_CYCLES: Number of pressure cycles (only used with "stretch" mode)
N_CYCLES = 5

# ── TIME STEP ──────────────────────────────────────────────────────────────────
# dt_hours: Time step size in hours (for operation phase)
dt_hours = 2

# ── CSV SETTINGS (only used when PRESSURE_SCENARIO = "csv") ────────────────────
# CSV_FILE_PATH: Path to the pressure profile CSV file
# Expected columns: Druk_MPa or Druk_bar (hourly values)
CSV_FILE_PATH = "drukprofiel_zoutcaverne_2035_8760u.csv"

# CSV_SHIFT_TO_LEACH_END: Whether to shift CSV profile so its minimum equals p_leach_end
#   - True (default): Shifts entire CSV profile so minimum pressure = p_leach_end
#     This preserves the volume swing (pressure range) while starting from leaching end pressure.
#   - False: Use original CSV values (may not match leaching end pressure)
CSV_SHIFT_TO_LEACH_END = True

# RAMP_HOURS: Startup ramp duration (hours) to smoothly transition at start of operation.
#   Helps avoid numerical instabilities (NaNs) if there's a pressure jump.
#   Set to 0 to disable ramping.
RAMP_HOURS = 24.0

# ── MATERIAL MODEL ─────────────────────────────────────────────────────────────
# USE_DESAI: Enable Desai viscoplastic model during operation phase
#   Note: Desai is NEVER enabled during leaching phase
USE_DESAI = True

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                        END OF USER CONFIGURATION                              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


# ══════════════════════════════════════════════════════════════════════════════
# AUTOMATIC CONFIGURATION (do not modify below unless you know what you're doing)
# ══════════════════════════════════════════════════════════════════════════════

# Valid options for validation
VALID_SHAPES = ["regular", "tilted", "teardrop", "asymmetric", "irregular", "multichamber"]
VALID_SIZES = [600, 1200]
VALID_SCENARIOS = ["sinus", "linear", "irregular", "csv"]
VALID_MODES_STANDARD = ["stretch", "repeat"]
VALID_MODES_CSV = ["stretch", "repeat", "direct"]
VALID_LEACHING_MODES = ["linear", "stepped"]

# Cavern z_max values (top of cavern elevation in meters, from domain bottom at z=0)
Z_MAX_BY_CAVERN = {
    "regular600": 315.26,
    "tilted600": 345.67,
    "teardrop600": 353.15,
    "asymmetric600": 338.89,
    "irregular600": 319.86,
    "multichamber600": 334.14,
    "regular1200": 393.21,
    "tilted1200": 430.78,
    "teardrop1200": 445.06,
    "asymmetric1200": 422.76,
    "irregular1200": 402.21,
    "multichamber1200": 420.82,
}

# Estimated cavern heights (meters) for center calculation
# These are approximate values based on typical cavern geometries
CAVERN_HEIGHT_BY_TYPE = {
    "regular600": 150.0,
    "tilted600": 160.0,
    "teardrop600": 180.0,
    "asymmetric600": 155.0,
    "irregular600": 145.0,
    "multichamber600": 140.0,
    "regular1200": 200.0,
    "tilted1200": 215.0,
    "teardrop1200": 240.0,
    "asymmetric1200": 210.0,
    "irregular1200": 195.0,
    "multichamber1200": 190.0,
}

# Reference pressure by cavern size (for overburden/sideburden at surface z=660m)
P_REF_BY_SIZE = {
    600: 17.5,   # MPa
    1200: 19.3,  # MPa
}

# Domain dimensions
Z_SURFACE = 660.0  # meters (top of domain)


def validate_configuration():
    """Validate user configuration and return derived values."""
    errors = []

    if CAVERN_SHAPE not in VALID_SHAPES:
        errors.append(f"CAVERN_SHAPE '{CAVERN_SHAPE}' invalid. Choose from: {VALID_SHAPES}")
    if CAVERN_SIZE not in VALID_SIZES:
        errors.append(f"CAVERN_SIZE '{CAVERN_SIZE}' invalid. Choose from: {VALID_SIZES}")
    if PRESSURE_SCENARIO not in VALID_SCENARIOS:
        errors.append(f"PRESSURE_SCENARIO '{PRESSURE_SCENARIO}' invalid. Choose from: {VALID_SCENARIOS}")
    if LEACHING_MODE not in VALID_LEACHING_MODES:
        errors.append(f"LEACHING_MODE '{LEACHING_MODE}' invalid. Choose from: {VALID_LEACHING_MODES}")

    # Validate SCHEDULE_MODE based on scenario
    if PRESSURE_SCENARIO == "csv":
        if SCHEDULE_MODE not in VALID_MODES_CSV:
            errors.append(f"SCHEDULE_MODE '{SCHEDULE_MODE}' invalid for CSV. Choose from: {VALID_MODES_CSV}")
    else:
        if SCHEDULE_MODE not in VALID_MODES_STANDARD:
            errors.append(f"SCHEDULE_MODE '{SCHEDULE_MODE}' invalid. Choose from: {VALID_MODES_STANDARD}")

    if OPERATION_DAYS <= 0:
        errors.append(f"OPERATION_DAYS must be positive, got {OPERATION_DAYS}")
    if LEACHING_DAYS <= 0:
        errors.append(f"LEACHING_DAYS must be positive, got {LEACHING_DAYS}")
    if N_CYCLES <= 0:
        errors.append(f"N_CYCLES must be positive, got {N_CYCLES}")
    if dt_hours <= 0:
        errors.append(f"dt_hours must be positive, got {dt_hours}")
    if LEACHING_DT_HOURS <= 0:
        errors.append(f"LEACHING_DT_HOURS must be positive, got {LEACHING_DT_HOURS}")
    if STEPPED_N_STEPS < 2:
        errors.append(f"STEPPED_N_STEPS must be at least 2, got {STEPPED_N_STEPS}")
    if LEACHING_END_FRACTION <= 0.0 or LEACHING_END_FRACTION >= 1.0:
        errors.append(f"LEACHING_END_FRACTION must be between 0 and 1 (exclusive), got {LEACHING_END_FRACTION}")

    # Pressure scenario validation
    if PRESSURE_SCENARIO == "sinus":
        if P_AMPLITUDE_MPA <= 0:
            errors.append(f"P_AMPLITUDE_MPA must be positive, got {P_AMPLITUDE_MPA}")
    elif PRESSURE_SCENARIO == "linear":
        if PRESSURE_SWING_MPA <= 0:
            errors.append(f"PRESSURE_SWING_MPA must be positive, got {PRESSURE_SWING_MPA}")
    elif PRESSURE_SCENARIO == "csv":
        if not os.path.isfile(CSV_FILE_PATH):
            errors.append(f"CSV_FILE_PATH not found: {CSV_FILE_PATH}")

    if errors:
        raise ValueError("Configuration errors:\n  - " + "\n  - ".join(errors))

    # Derived values
    cavern_type = f"{CAVERN_SHAPE}{CAVERN_SIZE}"
    grid_folder = f"cavern_{CAVERN_SHAPE}_{CAVERN_SIZE}_3D"
    z_max = Z_MAX_BY_CAVERN[cavern_type]
    cavern_height = CAVERN_HEIGHT_BY_TYPE[cavern_type]
    z_center = z_max - cavern_height / 2.0  # Cavern center elevation
    p_ref_mpa = P_REF_BY_SIZE[CAVERN_SIZE]

    return {
        "cavern_type": cavern_type,
        "grid_folder": grid_folder,
        "z_max": z_max,
        "z_center": z_center,
        "cavern_height": cavern_height,
        "p_ref_mpa": p_ref_mpa,
    }


def compute_lithostatic_pressure(z_center, p_ref_mpa, rho_salt, g):
    """
    Compute lithostatic pressure at cavern center.

    Parameters
    ----------
    z_center : float
        Cavern center elevation (m) from domain bottom
    p_ref_mpa : float
        Reference pressure at surface (MPa)
    rho_salt : float
        Salt density (kg/m³)
    g : float
        Gravitational acceleration (m/s², typically negative)

    Returns
    -------
    p_lithostatic : float
        Lithostatic pressure in Pa
    """
    # Depth below surface
    depth = Z_SURFACE - z_center

    # Lithostatic pressure = surface pressure + salt column pressure
    p_lithostatic = p_ref_mpa * ut.MPa + rho_salt * abs(g) * depth

    return p_lithostatic


def build_leaching_pressure_schedule(tc, *, p_start_pa, p_end_pa, mode, n_steps=6):
    """
    Build pressure schedule for leaching phase.

    Parameters
    ----------
    tc : TimeController
        Time controller for leaching phase
    p_start_pa : float
        Starting pressure (lithostatic) in Pa
    p_end_pa : float
        Ending pressure (operational starting) in Pa
    mode : str
        "linear" for linear decrease, "stepped" for stepped decrease
    n_steps : int
        Number of steps for stepped mode

    Returns
    -------
    t_vals : list
        Time values in seconds
    p_vals : list
        Pressure values in Pa
    """
    n_time_steps = int(math.floor(tc.t_final / tc.dt))
    t_vals = [k * tc.dt for k in range(n_time_steps + 1)]
    if abs(t_vals[-1] - tc.t_final) > 1e-12:
        t_vals.append(tc.t_final)

    if mode == "linear":
        # Linear interpolation from start to end
        p_vals = []
        for t in t_vals:
            frac = t / tc.t_final if tc.t_final > 0 else 1.0
            p = p_start_pa + frac * (p_end_pa - p_start_pa)
            p_vals.append(p)

    elif mode == "stepped":
        # Stepped decrease with plateaus
        step_duration = tc.t_final / n_steps
        p_step_values = np.linspace(p_start_pa, p_end_pa, n_steps + 1)

        p_vals = []
        for t in t_vals:
            step_idx = min(int(t / step_duration), n_steps - 1)
            # Hold at step value, then drop to next at step boundary
            if t < tc.t_final:
                p = p_step_values[step_idx]
            else:
                p = p_end_pa
            p_vals.append(p)

    else:
        raise ValueError(f"Unknown leaching mode: {mode}")

    return t_vals, p_vals


# ══════════════════════════════════════════════════════════════════════════════
# CSV PRESSURE PROFILE FUNCTIONS (same as Run.py)
# ══════════════════════════════════════════════════════════════════════════════

def _parse_float_auto(s: str) -> float:
    """Parse floats with either '.' or ',' as decimal separator."""
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


def read_pressure_csv(csv_file: str):
    """
    Read CSV and return pressure array in MPa.
    Supported columns (case-insensitive): Druk_MPa, Druk_bar (converted via /10).
    """
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

    rows = []
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

    pressures_mpa = []
    if idx_mpa is not None:
        for r in data:
            if idx_mpa < len(r):
                pressures_mpa.append(_parse_float_auto(r[idx_mpa]))
    elif idx_bar is not None:
        for r in data:
            if idx_bar < len(r):
                pressures_mpa.append(_parse_float_auto(r[idx_bar]) / 10.0)
    else:
        # Fallback: find first numeric column
        ncols = len(header)
        best_i, best_count = None, -1
        for i in range(ncols):
            vals = [_parse_float_auto(r[i]) for r in data if i < len(r)]
            count = np.sum(np.isfinite(vals))
            if count > best_count:
                best_count, best_i = count, i
        if best_i is None or best_count < 2:
            raise ValueError("Could not find a numeric pressure column in CSV")
        for r in data:
            if best_i < len(r):
                pressures_mpa.append(_parse_float_auto(r[best_i]))

    pressures_mpa = np.asarray(pressures_mpa, dtype=float)
    pressures_mpa = pressures_mpa[np.isfinite(pressures_mpa)]
    if pressures_mpa.size < 2:
        raise ValueError("Parsed pressure series has <2 numeric values")

    return pressures_mpa


def rescale_pressure_profile(pressures_mpa, new_min, new_max):
    """
    Rescale pressure profile from its original range to [new_min, new_max].
    Preserves relative volume (fraction of capacity) at each time step.
    """
    old_min = pressures_mpa.min()
    old_max = pressures_mpa.max()

    if old_max - old_min < 1e-9:
        return np.full_like(pressures_mpa, (new_min + new_max) / 2.0)

    fraction = (pressures_mpa - old_min) / (old_max - old_min)
    return new_min + fraction * (new_max - new_min)


def build_csv_pressure_schedule(tc, csv_file, *, days, mode, total_cycles=1,
                                rescale=False, rescale_min=None, rescale_max=None,
                                resample_at_dt=True):
    """Build pressure schedule from CSV file."""
    pressures_mpa = read_pressure_csv(csv_file)
    csv_hours = int(pressures_mpa.size)

    if MPI.COMM_WORLD.rank == 0:
        print(f"[CSV] Loaded '{os.path.basename(csv_file)}' with {csv_hours} hourly values")
        print(f"[CSV] Original range: [{pressures_mpa.min():.2f}, {pressures_mpa.max():.2f}] MPa")

    if rescale and rescale_min is not None and rescale_max is not None:
        pressures_mpa = rescale_pressure_profile(pressures_mpa, rescale_min, rescale_max)
        if MPI.COMM_WORLD.rank == 0:
            print(f"[CSV] Rescaled to: [{pressures_mpa.min():.2f}, {pressures_mpa.max():.2f}] MPa")

    total_hours = float(days) * 24.0

    if mode == "direct":
        sim_hours = np.arange(0.0, total_hours + 1e-12, 1.0)
        idx = (sim_hours % csv_hours).astype(int)
        times_hours = sim_hours
        pressures_MPa = pressures_mpa[idx]

    elif mode == "stretch":
        total_cycles = max(1, int(total_cycles))
        cycle_duration_h = total_hours / float(total_cycles)
        scale = cycle_duration_h / float(csv_hours)

        times_list, pres_list = [], []
        for k in range(total_cycles):
            off = k * cycle_duration_h
            for i in range(csv_hours):
                if k > 0 and i == 0:
                    continue
                times_list.append(off + i * scale)
                pres_list.append(pressures_mpa[i])
        times_hours = np.asarray(times_list, dtype=float)
        pressures_MPa = np.asarray(pres_list, dtype=float)

    elif mode == "repeat":
        n_rep = int(np.ceil(total_hours / float(csv_hours)))
        times_list, pres_list = [], []
        for r in range(n_rep):
            off = r * csv_hours
            for i in range(csv_hours):
                if r > 0 and i == 0:
                    continue
                t = off + i
                if t > total_hours:
                    break
                times_list.append(t)
                pres_list.append(pressures_mpa[i])
        times_hours = np.asarray(times_list, dtype=float)
        pressures_MPa = np.asarray(pres_list, dtype=float)

    else:
        raise ValueError("mode must be 'direct', 'stretch', or 'repeat'")

    times_s = times_hours * 3600.0

    if times_s[0] > 0.0:
        times_s = np.insert(times_s, 0, 0.0)
        pressures_MPa = np.insert(pressures_MPa, 0, pressures_MPa[0])
    if times_s[-1] < tc.t_final:
        times_s = np.append(times_s, tc.t_final)
        pressures_MPa = np.append(pressures_MPa, pressures_MPa[-1])

    if resample_at_dt:
        n_steps = int(math.floor(tc.t_final / tc.dt))
        t_vals = [k * tc.dt for k in range(n_steps + 1)]
        if abs(t_vals[-1] - tc.t_final) > 1e-12:
            t_vals.append(tc.t_final)
        p_vals_MPa = np.interp(t_vals, times_s, pressures_MPa)
    else:
        t_vals = times_s.tolist()
        p_vals_MPa = pressures_MPa.tolist()

    p_vals = [float(p) * ut.MPa for p in p_vals_MPa]

    if MPI.COMM_WORLD.rank == 0:
        print(f"[CSV] Schedule built: {len(t_vals)} points, mode={mode}")

    return t_vals, p_vals


def apply_startup_ramp(t_pressure, p_pressure, *, p_start_pa, ramp_hours, dt_hours):
    """
    Replace first part of schedule with a linear ramp from p_start_pa to the existing schedule.
    Operates in-place on p_pressure list.
    """
    if ramp_hours is None or ramp_hours <= 0.0:
        p_pressure[0] = p_start_pa
        return

    ramp_steps = max(1, int(round(float(ramp_hours) / float(dt_hours))))
    ramp_steps = min(ramp_steps, len(p_pressure) - 1)

    p_target = p_pressure[ramp_steps]
    p_pressure[0] = p_start_pa
    for k in range(1, ramp_steps + 1):
        a = k / float(ramp_steps)
        p_pressure[k] = (1.0 - a) * p_start_pa + a * p_target


# ══════════════════════════════════════════════════════════════════════════════
# PRESSURE SCHEDULE BUILDERS (same as Run.py)
# ══════════════════════════════════════════════════════════════════════════════

DAY_H = 24.0


def _sample_at_dt(tc, t_end=None):
    t_end = tc.t_final if t_end is None else t_end
    n_steps = int(math.floor(t_end / tc.dt))
    t_vals = [k * tc.dt for k in range(n_steps + 1)]
    if abs(t_vals[-1] - t_end) > 1e-12:
        t_vals.append(t_end)
    return t_vals


def build_sinus_pressure_schedule(tc, *, p_mean, p_ampl, period_hours, phase_hours=0.0,
                                  clamp_min=None, clamp_max=None):
    """Sinus schedule sampled at simulation time steps."""
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
        if clamp_min is not None:
            p = max(p, clamp_min)
        if clamp_max is not None:
            p = min(p, clamp_max)
        p_vals.append(p)

    return t_vals, p_vals


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


def build_irregular_pressure_schedule(tc, times_hours, pressures_MPa, *, smooth=0.3,
                                      clamp_min=None, clamp_max=None, resample_at_dt=True):
    if len(times_hours) != len(pressures_MPa) or len(times_hours) < 2:
        raise ValueError("Provide at least two waypoints with matching lengths.")

    knots_t = np.array([h * ut.hour for h in times_hours], dtype=float)
    knots_p = np.array([p * ut.MPa for p in pressures_MPa], dtype=float)

    if knots_t[0] > 0.0:
        knots_t = np.insert(knots_t, 0, 0.0)
        knots_p = np.insert(knots_p, 0, knots_p[0])
    if knots_t[-1] < tc.t_final:
        knots_t = np.append(knots_t, tc.t_final)
        knots_p = np.append(knots_p, knots_p[-1])

    if resample_at_dt:
        n_steps = int(math.floor(tc.t_final / tc.dt))
        t_vals = [k * tc.dt for k in range(n_steps + 1)]
        if abs(t_vals[-1] - tc.t_final) > 1e-12:
            t_vals.append(tc.t_final)
    else:
        t_vals = knots_t.tolist()

    p_vals = []
    if smooth is None:
        for t in t_vals:
            p = np.interp(t, knots_t, knots_p)
            if clamp_min is not None:
                p = max(p, clamp_min)
            if clamp_max is not None:
                p = min(p, clamp_max)
            p_vals.append(p)
    else:
        tau = float(np.clip(smooth, 0.0, 1.0))
        for t in t_vals:
            p = _cardinal_interp(knots_t, knots_p, t, tau)
            if clamp_min is not None:
                p = max(p, clamp_min)
            if clamp_max is not None:
                p = min(p, clamp_max)
            p_vals.append(p)

    return t_vals, p_vals


def _stretch_hours(times_h, factor_days):
    return [t * float(factor_days) for t in times_h]


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


def build_linear_schedule_multi(tc, times_h, pressures_MPa, *, days, mode,
                                resample_at_dt=True, total_cycles=None):
    if mode not in ("repeat", "stretch"):
        raise ValueError("mode must be 'repeat' or 'stretch'")

    times_h = list(map(float, times_h))
    pressures_MPa = list(map(float, pressures_MPa))
    if len(times_h) != len(pressures_MPa) or len(times_h) < 2:
        raise ValueError("Provide at least two waypoints with matching lengths.")

    total_h = float(days) * DAY_H

    if mode == "repeat":
        t_h = _repeat_hours(times_h, days)
        p_h = []
        for d in range(int(days)):
            start = 0 if d == 0 else 1
            p_h.extend(pressures_MPa[start:])

    else:  # mode == "stretch"
        if total_cycles is None:
            total_cycles = 1
        total_cycles = max(1, int(total_cycles))

        base_start = times_h[0]
        base_end = times_h[-1]
        base_duration = base_end - base_start
        if base_duration <= 0.0:
            raise ValueError("times_h must span a positive duration.")

        cycle_duration = total_h / float(total_cycles)
        scale = cycle_duration / base_duration

        t_h = []
        p_h = []
        for k in range(total_cycles):
            offset = k * cycle_duration
            for i, t in enumerate(times_h):
                if k > 0 and i == 0:
                    continue
                t_scaled = offset + (t - base_start) * scale
                t_h.append(t_scaled)
                p_h.append(pressures_MPa[i])

    if t_h[0] > 0.0:
        t_h.insert(0, 0.0)
        p_h.insert(0, p_h[0])
    if t_h[-1] < total_h:
        t_h.append(total_h)
        p_h.append(p_h[-1])

    knots_t = np.array([h * ut.hour for h in t_h], dtype=float)
    knots_p = np.array([p * ut.MPa for p in p_h], dtype=float)

    if resample_at_dt:
        t_vals = _sample_at_dt(tc)
        p_vals = np.interp(t_vals, knots_t, knots_p).tolist()
    else:
        t_vals = knots_t.tolist()
        p_vals = knots_p.tolist()

    return t_vals, p_vals


def build_irregular_schedule_multi(tc, *, base_waypoints_h, base_pressures_MPa,
                                   days, mode, smooth=0.25,
                                   clamp_min=0.0, clamp_max=None,
                                   resample_at_dt=True, total_cycles=None):
    times_h = np.asarray(base_waypoints_h, dtype=float)
    pressures = np.asarray(base_pressures_MPa, dtype=float)

    if len(times_h) != len(pressures):
        raise ValueError("base_waypoints_h and base_pressures_MPa must have same length")

    total_hours = days * DAY_H

    if mode == "repeat":
        times_h_multi = _repeat_hours(times_h, days)
        pressures_multi = []
        for d in range(int(days)):
            start = 0 if d == 0 else 1
            pressures_multi.extend(pressures[start:])

    elif mode == "stretch":
        if total_cycles is None:
            total_cycles = 1
        total_cycles = max(1, int(total_cycles))

        base_start = times_h[0]
        base_end = times_h[-1]
        base_duration = base_end - base_start
        if base_duration <= 0.0:
            raise ValueError("base_waypoints_h must span a positive duration.")

        cycle_duration = total_hours / float(total_cycles)
        scale = cycle_duration / base_duration

        times_h_multi = []
        pressures_multi = []

        for k in range(total_cycles):
            offset = k * cycle_duration
            for i, t in enumerate(times_h):
                if k > 0 and i == 0:
                    continue
                t_scaled = offset + (t - base_start) * scale
                times_h_multi.append(t_scaled)
                pressures_multi.append(pressures[i])

    else:
        raise ValueError("mode must be 'repeat' or 'stretch'")

    times_h_multi = np.asarray(times_h_multi, dtype=float)
    pressures_multi = np.asarray(pressures_multi, dtype=float)

    return build_irregular_pressure_schedule(
        tc,
        times_hours=times_h_multi.tolist(),
        pressures_MPa=pressures_multi.tolist(),
        smooth=smooth,
        clamp_min=clamp_min, clamp_max=clamp_max,
        resample_at_dt=resample_at_dt
    )


def build_sinus_schedule_multi(tc, *, p_mean, p_ampl, days, mode,
                               daily_period_hours=24.0, total_cycles=1,
                               clamp_min=0.0, clamp_max=None):
    total_hours = days * DAY_H

    if mode == "repeat":
        T_hours = daily_period_hours
    elif mode == "stretch":
        total_cycles = max(1, int(total_cycles))
        T_hours = total_hours / float(total_cycles)
    else:
        raise ValueError("mode must be 'repeat' or 'stretch'")

    return build_sinus_pressure_schedule(
        tc, p_mean=p_mean, p_ampl=p_ampl,
        period_hours=T_hours, phase_hours=0.0,
        clamp_min=clamp_min, clamp_max=clamp_max
    )


# ══════════════════════════════════════════════════════════════════════════════
# MOMENTUM EQUATION CLASSES
# ══════════════════════════════════════════════════════════════════════════════

class LinearMomentumMod(sf.LinearMomentum):
    def __init__(self, grid, theta):
        super().__init__(grid, theta)
        self.expect_vp_state = False

    def initialize(self) -> None:
        self.C.x.array[:] = to.flatten(self.mat.C)
        self.Fvp = do.fem.Function(self.DG0_1)
        self.alpha = do.fem.Function(self.DG0_1)
        self.eps_vp = do.fem.Function(self.DG0_3x3)

    def run_after_solve(self):
        if not hasattr(self, "eps_vp"):
            return

        elems = getattr(self.mat, "elems_ne", None)
        if not elems:
            return
        st = elems[-1]

        if hasattr(st, "eps_ne_k"):
            self.eps_vp.x.array[:] = to.flatten(st.eps_ne_k)

        if self.expect_vp_state:
            if not (hasattr(st, "Fvp") and hasattr(st, "alpha")):
                if MPI.COMM_WORLD.rank == 0:
                    print("[WARN] Expected Fvp/alpha but missing.")
                return
            self.Fvp.x.array[:] = st.Fvp
            self.alpha.x.array[:] = st.alpha
        else:
            if hasattr(st, "Fvp") and hasattr(st, "alpha"):
                self.Fvp.x.array[:] = st.Fvp
                self.alpha.x.array[:] = st.alpha


class SparseSaveFields(sf.SaveFields):
    """SaveFields that only writes every `interval`-th call after t=0."""
    def __init__(self, mom_eq, interval: int):
        super().__init__(mom_eq)
        self.interval = max(1, int(interval))
        self._counter = 0

    def save_fields(self, t):
        if t == 0:
            return super().save_fields(t)

        self._counter += 1
        if self._counter % self.interval == 0:
            return super().save_fields(t)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    # Validate configuration and get derived values
    config = validate_configuration()

    # Material constants (needed for lithostatic calculation)
    salt_density = 2200  # kg/m³
    g = -9.81  # m/s²

    # Compute lithostatic pressure at cavern center
    p_lithostatic = compute_lithostatic_pressure(
        config["z_center"],
        config["p_ref_mpa"],
        salt_density,
        g
    )
    p_lithostatic_mpa = p_lithostatic / ut.MPa

    # Compute leaching end pressure (fraction of lithostatic)
    p_leach_end_mpa = LEACHING_END_FRACTION * p_lithostatic_mpa
    p_leach_end = p_leach_end_mpa * ut.MPa

    # Derive operational pressures from leaching end pressure
    # All pressure schemes start at p_leach_end as their minimum
    if PRESSURE_SCENARIO == "sinus":
        p_op_min_mpa = p_leach_end_mpa
        p_op_mean_mpa = p_leach_end_mpa + P_AMPLITUDE_MPA
        p_op_max_mpa = p_leach_end_mpa + 2 * P_AMPLITUDE_MPA
    elif PRESSURE_SCENARIO == "linear":
        p_op_min_mpa = p_leach_end_mpa
        p_op_max_mpa = p_leach_end_mpa + PRESSURE_SWING_MPA
    elif PRESSURE_SCENARIO == "csv":
        p_op_min_mpa = p_leach_end_mpa  # CSV will be shifted to start here
    else:  # irregular
        p_op_min_mpa = p_leach_end_mpa  # Will be shifted to start here

    # Print configuration summary
    if MPI.COMM_WORLD.rank == 0:
        print("=" * 70)
        print("SIMULATION CONFIGURATION (WITH LEACHING PHASE)")
        print("=" * 70)
        print(f"  Cavern:           {config['cavern_type']} ({CAVERN_SHAPE}, {CAVERN_SIZE},000 m³)")
        print(f"  Cavern z_max:     {config['z_max']:.2f} m")
        print(f"  Cavern z_center:  {config['z_center']:.2f} m")
        print(f"  Cavern height:    {config['cavern_height']:.2f} m (estimated)")
        print("-" * 70)
        print("  LEACHING PHASE:")
        print(f"    Mode:           {LEACHING_MODE}")
        print(f"    Duration:       {LEACHING_DAYS} days ({LEACHING_DAYS/365:.2f} years)")
        print(f"    dt:             {LEACHING_DT_HOURS} hours")
        if LEACHING_MODE == "stepped":
            print(f"    N steps:        {STEPPED_N_STEPS}")
        print(f"    P_start:        {p_lithostatic_mpa:.2f} MPa (100% lithostatic)")
        print(f"    P_end:          {p_leach_end_mpa:.2f} MPa ({LEACHING_END_FRACTION*100:.0f}% lithostatic)")
        print(f"    Delta P:        {p_lithostatic_mpa - p_leach_end_mpa:.2f} MPa")
        print("-" * 70)
        print("  OPERATION PHASE:")
        print(f"    Scenario:       {PRESSURE_SCENARIO}")
        if PRESSURE_SCENARIO == "sinus":
            print(f"    P_min:          {p_op_min_mpa:.2f} MPa (= p_leach_end)")
            print(f"    P_mean:         {p_op_mean_mpa:.2f} MPa (derived)")
            print(f"    P_max:          {p_op_max_mpa:.2f} MPa (derived)")
            print(f"    Amplitude:      {P_AMPLITUDE_MPA} MPa")
        elif PRESSURE_SCENARIO == "linear":
            print(f"    P_min:          {p_op_min_mpa:.2f} MPa (= p_leach_end)")
            print(f"    P_max:          {p_op_max_mpa:.2f} MPa (derived)")
            print(f"    Swing:          {PRESSURE_SWING_MPA} MPa")
        elif PRESSURE_SCENARIO == "csv":
            print(f"    P_min:          {p_op_min_mpa:.2f} MPa (= p_leach_end)")
            print(f"    CSV file:       {CSV_FILE_PATH}")
            print(f"    Shift to p_end: {'enabled' if CSV_SHIFT_TO_LEACH_END else 'disabled'}")
            print(f"    Ramp:           {RAMP_HOURS} hours")
        elif PRESSURE_SCENARIO == "irregular":
            print(f"    P_min:          {p_op_min_mpa:.2f} MPa (= p_leach_end, profile shifted)")
        print(f"    Mode:           {SCHEDULE_MODE}")
        print(f"    Days:           {OPERATION_DAYS}")
        print(f"    Cycles:         {N_CYCLES}")
        print(f"    dt:             {dt_hours} hours")
        print(f"    Desai:          {'enabled' if USE_DESAI else 'disabled'}")
        print("-" * 70)
        print(f"  Grid:             {config['grid_folder']}")
        print("=" * 70)

    # Load grid
    grid_path = os.path.join("..", "..", "..", "..", "grids", config["grid_folder"])
    grid = sf.GridHandlerGMSH("geom", grid_path)

    # Get derived values
    z_max = config["z_max"]
    p_ref = config["p_ref_mpa"] * ut.MPa

    # Output folder
    output_folder = os.path.join(
        "output",
        f"case_leaching_{LEACHING_MODE}_{PRESSURE_SCENARIO}({N_CYCLES})_{OPERATION_DAYS}days_{config['cavern_type']}"
    )

    side_burden = p_ref
    over_burden = p_ref

    # Define momentum equation
    mom_eq = LinearMomentumMod(grid, theta=0.5)

    # Define solver
    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("cg")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-10, max_it=100)
    mom_eq.set_solver(mom_solver)

    # Define material properties
    mat = sf.Material(mom_eq.n_elems)

    # Set material density
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    # Elastic
    E0 = 20.425 * ut.GPa * to.ones(mom_eq.n_elems)
    nu0 = 0.25 * to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")

    # Kelvin-Voigt
    eta = 105e11 * to.ones(mom_eq.n_elems)
    E1 = 10 * ut.GPa * to.ones(mom_eq.n_elems)
    nu1 = 0.25 * to.ones(mom_eq.n_elems)
    kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")

    sec_per_year = 365.25 * 24 * 3600

    # Dislocation creep
    ndc = 4.6
    A_dc = (40.0 * (1e-6)**ndc / sec_per_year) * to.ones(mom_eq.n_elems)
    Q_dc = (6495.0 * 8.32) * to.ones(mom_eq.n_elems)
    n_dc = ndc * to.ones(mom_eq.n_elems)
    creep_0 = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation")

    # Pressure-solution creep
    A_ps = (14176.0 * 1e-9 / 1e6 / sec_per_year) * to.ones(mom_eq.n_elems)
    d_ps = 5.25e-3 * to.ones(mom_eq.n_elems)
    Q_ps = (3252.0 * 8.32) * to.ones(mom_eq.n_elems)
    creep_pressure = sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "creep_pressure")

    # Constitutive model (NO Desai during leaching)
    mat.add_to_elastic(spring_0)
    mat.add_to_non_elastic(kelvin)
    mat.add_to_non_elastic(creep_0)
    mat.add_to_non_elastic(creep_pressure)

    mom_eq.set_material(mat)

    # Body forces
    g_vec = [0.0, 0.0, g]
    mom_eq.build_body_force(g_vec)

    # Initial temperature
    T0_field = 298 * to.ones(mom_eq.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # ===== LEACHING PHASE =====
    tc_leaching = sf.TimeController(
        dt=LEACHING_DT_HOURS,
        initial_time=0.0,
        final_time=LEACHING_DAYS * 24.0,
        time_unit="hour"
    )

    # Build leaching pressure schedule
    t_leaching, p_leaching = build_leaching_pressure_schedule(
        tc_leaching,
        p_start_pa=p_lithostatic,
        p_end_pa=p_leach_end,
        mode=LEACHING_MODE,
        n_steps=STEPPED_N_STEPS
    )

    gas_density = 0.089  # Hydrogen density in kg/m³

    # Leaching BCs
    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_leaching.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_leaching.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_leaching.t_final])

    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0,
                              [side_burden, side_burden],
                              [0.0, tc_leaching.t_final],
                              g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0,
                               [side_burden, side_burden],
                               [0.0, tc_leaching.t_final],
                               g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                             [over_burden, over_burden],
                             [0.0, tc_leaching.t_final],
                             g=g_vec[2])

    bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, z_max,
                                p_leaching,
                                t_leaching,
                                g=g_vec[2])

    bc_leaching = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
        bc_leaching.add_boundary_condition(bc)

    mom_eq.set_boundary_conditions(bc_leaching)

    output_folder_leaching = os.path.join(output_folder, "leaching")
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[LEACHING] Output: {output_folder_leaching}")

    # Sparse saving during leaching (every 6th step ~ every 3 days with 12h dt)
    leaching_save_interval = max(1, int(72 / LEACHING_DT_HOURS))  # Save every ~3 days
    output_mom_leaching = SparseSaveFields(mom_eq, interval=leaching_save_interval)
    output_mom_leaching.set_output_folder(output_folder_leaching)
    output_mom_leaching.add_output_field("u", "Displacement (m)")
    output_mom_leaching.add_output_field("eps_tot", "Total strain (-)")
    output_mom_leaching.add_output_field("sig", "Stress (Pa)")
    output_mom_leaching.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom_leaching.add_output_field("q_elems", "Von Mises stress (Pa)")
    outputs_leaching = [output_mom_leaching]

    # Save leaching pressure schedule
    os.makedirs(output_folder, exist_ok=True)
    leaching_data = {
        "phase": "leaching",
        "cavern_type": config["cavern_type"],
        "leaching_mode": LEACHING_MODE,
        "leaching_days": LEACHING_DAYS,
        "dt_hours": LEACHING_DT_HOURS,
        "stepped_n_steps": STEPPED_N_STEPS if LEACHING_MODE == "stepped" else None,
        "leaching_end_fraction": LEACHING_END_FRACTION,
        "p_lithostatic_mpa": p_lithostatic_mpa,
        "p_leach_end_mpa": p_leach_end_mpa,
        "z_center_m": config["z_center"],
        "t_values_s": [float(t) for t in t_leaching],
        "p_values_Pa": [float(p) for p in p_leaching],
        "t_hours": [float(t / ut.hour) for t in t_leaching],
        "p_MPa": [float(p / ut.MPa) for p in p_leaching],
    }
    with open(os.path.join(output_folder, "leaching_schedule.json"), 'w') as f:
        json.dump(leaching_data, f, indent=2)

    if MPI.COMM_WORLD.rank == 0:
        print(f"[LEACHING] Running for {LEACHING_DAYS} days ({len(t_leaching)} steps)...")

    sim_leaching = sf.Simulator_M(mom_eq, tc_leaching, outputs_leaching, True)
    sim_leaching.run()

    if MPI.COMM_WORLD.rank == 0:
        print("[LEACHING] Complete.")

    # ===== OPERATION PHASE =====
    if USE_DESAI:
        mu_1 = 5.3665857009859815e-11 * to.ones(mom_eq.n_elems)
        N_1 = 3.1 * to.ones(mom_eq.n_elems)
        n = 3.0 * to.ones(mom_eq.n_elems)
        a_1 = 1.965018496922832e-05 * to.ones(mom_eq.n_elems)
        eta_vp = 0.8275682807874163 * to.ones(mom_eq.n_elems)
        beta_1 = 0.0048 * to.ones(mom_eq.n_elems)
        beta = 0.995 * to.ones(mom_eq.n_elems)
        m = -0.5 * to.ones(mom_eq.n_elems)
        gamma = 0.095 * to.ones(mom_eq.n_elems)
        alpha_0 = 0.0022 * to.ones(mom_eq.n_elems)
        sigma_t = 5.0 * to.ones(mom_eq.n_elems)
        desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta_vp, n, beta_1, beta, m, gamma, sigma_t, alpha_0, "desai")

        stress_to = ut.numpy2torch(mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3)))
        desai.compute_initial_hardening(stress_to, Fvp_0=0.0)

        mat.add_to_non_elastic(desai)
        mom_eq.set_material(mat)
        mom_eq.expect_vp_state = True
    else:
        mom_eq.expect_vp_state = False

    tc_operation = sf.TimeController(
        dt=dt_hours,
        initial_time=0.0,
        final_time=OPERATION_DAYS * 24.0,
        time_unit="hour"
    )

    if PRESSURE_SCENARIO == "linear":
        # Linear cycle using derived pressures from leaching end
        base_times_h = [0.0, 2.0, 14.0, 16.0, 24.0]
        base_pressures_MPa = [p_op_min_mpa, p_op_max_mpa, p_op_max_mpa, p_op_min_mpa, p_op_min_mpa]

        t_pressure, p_pressure = build_linear_schedule_multi(
            tc_operation,
            base_times_h, base_pressures_MPa,
            days=OPERATION_DAYS,
            mode=SCHEDULE_MODE,
            resample_at_dt=True,
            total_cycles=N_CYCLES,
        )

    elif PRESSURE_SCENARIO == "sinus":
        # Sinus cycle using derived pressures from leaching end
        p_mean = p_op_mean_mpa * ut.MPa
        p_ampl = P_AMPLITUDE_MPA * ut.MPa
        t_pressure, p_pressure = build_sinus_schedule_multi(
            tc_operation,
            p_mean=p_mean, p_ampl=p_ampl,
            days=OPERATION_DAYS, mode=SCHEDULE_MODE,
            daily_period_hours=24.0,
            total_cycles=N_CYCLES,
            clamp_min=None, clamp_max=None
        )

    elif PRESSURE_SCENARIO == "irregular":
        # Irregular cycle: shift base profile so minimum = p_leach_end
        base_waypoints_h = [0, 1.0, 2.0, 3.2, 4.0, 5.0, 6.4, 7.1, 9.0, 11.5,
                           13.0, 16.0, 18.0, 21.0, 24.0]
        # Original base pressures (will be shifted)
        base_pressures_orig = [15.0, 12.0, 8.5, 11.8, 7.6, 10.2, 8.8, 11.4,
                               9.3, 10.7, 8.9, 11.6, 9.5, 10.2, 11.0]
        # Shift so minimum equals p_leach_end
        orig_min = min(base_pressures_orig)
        shift = p_op_min_mpa - orig_min
        base_pressures_MPa = [p + shift for p in base_pressures_orig]

        t_pressure, p_pressure = build_irregular_schedule_multi(
            tc_operation,
            base_waypoints_h=base_waypoints_h,
            base_pressures_MPa=base_pressures_MPa,
            days=OPERATION_DAYS, mode=SCHEDULE_MODE,
            smooth=0.25, clamp_min=None, clamp_max=None,
            resample_at_dt=True,
            total_cycles=N_CYCLES,
        )

    elif PRESSURE_SCENARIO == "csv":
        t_pressure, p_pressure = build_csv_pressure_schedule(
            tc_operation,
            csv_file=CSV_FILE_PATH,
            days=OPERATION_DAYS,
            mode=SCHEDULE_MODE,
            total_cycles=N_CYCLES,
            rescale=False,  # We'll shift instead of rescale
            rescale_min=None,
            rescale_max=None,
            resample_at_dt=True
        )

        # Shift CSV profile so minimum equals p_leach_end (preserves volume swing)
        if CSV_SHIFT_TO_LEACH_END:
            p_array = np.array(p_pressure)
            csv_min_pa = p_array.min()
            shift_pa = p_leach_end - csv_min_pa
            p_pressure = [p + shift_pa for p in p_pressure]
            if MPI.COMM_WORLD.rank == 0:
                new_min = min(p_pressure) / ut.MPa
                new_max = max(p_pressure) / ut.MPa
                print(f"[CSV] Shifted profile: new range [{new_min:.2f}, {new_max:.2f}] MPa")

        # Apply startup ramp (for convergence safety)
        apply_startup_ramp(
            t_pressure, p_pressure,
            p_start_pa=p_leach_end,
            ramp_hours=RAMP_HOURS,
            dt_hours=dt_hours
        )

    else:
        raise ValueError(f"Unknown PRESSURE_SCENARIO: {PRESSURE_SCENARIO}")

    # Operation BCs
    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_operation.t_final])

    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0,
                              [side_burden, side_burden],
                              [0.0, tc_operation.t_final],
                              g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0,
                               [side_burden, side_burden],
                               [0.0, tc_operation.t_final],
                               g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                             [over_burden, over_burden],
                             [0.0, tc_operation.t_final],
                             g=g_vec[2])
    bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, z_max,
                                p_pressure,
                                t_pressure,
                                g=g_vec[2])

    bc_operation = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
        bc_operation.add_boundary_condition(bc)

    mom_eq.set_boundary_conditions(bc_operation)

    output_folder_operation = os.path.join(output_folder, "operation")
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[OPERATION] Output: {output_folder_operation}")

    # Save operation pressure schedule
    pressure_data = {
        "phase": "operation",
        "cavern_type": config["cavern_type"],
        "cavern_shape": CAVERN_SHAPE,
        "cavern_size_m3": CAVERN_SIZE * 1000,
        "scenario": PRESSURE_SCENARIO,
        "pressure_scenario": PRESSURE_SCENARIO,
        "p_leach_end_mpa": p_leach_end_mpa,
        "leaching_end_fraction": LEACHING_END_FRACTION,
        "mode": SCHEDULE_MODE,
        "n_cycles": N_CYCLES,
        "operation_days": OPERATION_DAYS,
        "dt_hours": dt_hours,
        "use_desai": USE_DESAI,
        "leaching_mode": LEACHING_MODE,
        "leaching_days": LEACHING_DAYS,
        "p_lithostatic_mpa": p_lithostatic_mpa,
        "units": {"t_raw": "s", "p_raw": "Pa", "t": "hour", "p": "MPa"},
        "t_values_s": [float(t) for t in t_pressure],
        "p_values_Pa": [float(p) for p in p_pressure],
        "t_hours": [float(t / ut.hour) for t in t_pressure],
        "p_MPa": [float(p / ut.MPa) for p in p_pressure],
    }

    if PRESSURE_SCENARIO == "sinus":
        pressure_data["p_min_mpa"] = p_op_min_mpa
        pressure_data["p_mean_mpa"] = p_op_mean_mpa
        pressure_data["p_max_mpa"] = p_op_max_mpa
        pressure_data["p_amplitude_mpa"] = P_AMPLITUDE_MPA
    elif PRESSURE_SCENARIO == "linear":
        pressure_data["p_min_mpa"] = p_op_min_mpa
        pressure_data["p_max_mpa"] = p_op_max_mpa
        pressure_data["pressure_swing_mpa"] = PRESSURE_SWING_MPA
    elif PRESSURE_SCENARIO == "csv":
        pressure_data["csv_file"] = os.path.basename(CSV_FILE_PATH)
        pressure_data["csv_shift_to_leach_end"] = CSV_SHIFT_TO_LEACH_END
        pressure_data["ramp_hours"] = RAMP_HOURS
        pressure_data["p_min_mpa"] = min(p_pressure) / ut.MPa
        pressure_data["p_max_mpa"] = max(p_pressure) / ut.MPa
    elif PRESSURE_SCENARIO == "irregular":
        pressure_data["p_min_mpa"] = min(p_pressure) / ut.MPa
        pressure_data["p_max_mpa"] = max(p_pressure) / ut.MPa

    with open(os.path.join(output_folder, "pressure_schedule.json"), 'w') as f:
        json.dump(pressure_data, f, indent=2)

    # Operation outputs – use sparse saver (every 15th step)
    output_mom_op = SparseSaveFields(mom_eq, interval=15)
    output_mom_op.set_output_folder(output_folder_operation)
    output_mom_op.add_output_field("u", "Displacement (m)")
    output_mom_op.add_output_field("eps_tot", "Total strain (-)")
    output_mom_op.add_output_field("eps_vp", "Viscoplastic strain (-)")
    output_mom_op.add_output_field("alpha", "Hardening parameter (-)")
    output_mom_op.add_output_field("Fvp", "Yield function (-)")
    output_mom_op.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom_op.add_output_field("q_elems", "Von Mises stress (Pa)")
    output_mom_op.add_output_field("sig", "Stress (Pa)")
    outputs_op = [output_mom_op]

    if MPI.COMM_WORLD.rank == 0:
        print(f"[OPERATION] Running for {OPERATION_DAYS} days...")

    sim_op = sf.Simulator_M(mom_eq, tc_operation, outputs_op, False)
    sim_op.run()

    if MPI.COMM_WORLD.rank == 0:
        print("[OPERATION] Complete.")
        print("=" * 70)
        print("SIMULATION FINISHED")
        print(f"  Total simulated time: {LEACHING_DAYS + OPERATION_DAYS} days")
        print(f"  Output folder: {output_folder}")
        print("=" * 70)


if __name__ == '__main__':
    main()
