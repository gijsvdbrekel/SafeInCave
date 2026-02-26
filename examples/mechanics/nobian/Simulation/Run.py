import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
import safeincave.HeatBC as heatBC
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
# ║  This script supports homogeneous salt caverns with optional leaching phase.  ║
# ║  Combines functionality of Run.py and Run_leaching.py.                        ║
# ║                                                                               ║
# ║  INITIALIZATION MODE:                                                         ║
# ║  - USE_LEACHING = True:  Leaching phase (pressure ramp from lithostatic)      ║
# ║  - USE_LEACHING = False: Equilibrium phase (short constant pressure phase)    ║
# ║                                                                               ║
# ║  CAVERN SELECTION:                                                            ║
# ║  - Standard shapes: regular, tilted, directcirculation, asymmetric,           ║
# ║                     reversedcirculation, fastleached, tubefailure             ║
# ║                     (sizes: 600k or 1200k m³)                                ║
# ║  - Zuidwending A5: ~1.000.000k m³, real depth 1140-1510m                            ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── INITIALIZATION MODE ───────────────────────────────────────────────────────
# USE_LEACHING: Choose initialization strategy before operation phase:
#   True  - Leaching phase: pressure decreases from lithostatic to operational
#   False - Equilibrium phase: short phase at constant equilibrium pressure
USE_LEACHING = True

# EQUILIBRIUM_HOURS: Duration of equilibrium phase (only used if USE_LEACHING = False)
EQUILIBRIUM_HOURS = 10.0

# EQUILIBRIUM_DT_HOURS: Time step during equilibrium (only used if USE_LEACHING = False)
EQUILIBRIUM_DT_HOURS = 0.5

# ── CAVERN SELECTION ───────────────────────────────────────────────────────────
# CAVERN_TYPE: Choose one of:
#   Standard shapes (with CAVERN_SIZE 600 or 1200):
#     "regular"              - Standard cylindrical cavern
#     "tilted"               - Tilted/inclined cavern
#     "directcirculation"    - Direct-circulation leached cavern
#     "asymmetric"           - Asymmetric cavern geometry
#     "reversedcirculation"  - Reversed-circulation leached cavern
#     "fastleached"          - Fast-leached rough barrel shape
#     "tubefailure"          - Tube-failure (exaggerated multi-chamber)
#   Special caverns (CAVERN_SIZE is ignored):
#     "A5"           - Zuidwending cavern A5 (~1.000.000k m³, depth 1140-1510m)

CAVERN_TYPE = "regular"

# CAVERN_SIZE: Volume in thousands of m³ (ignored for A5)
#   600  - 600,000 m³ volume
#   1200 - 1,200,000 m³ volume
CAVERN_SIZE = 1200

# ── LEACHING PHASE SETTINGS ──────────────────────────────────────────────────────
# LEACHING_MODE: How pressure decreases during leaching:
#   "linear"  - Linear decrease from lithostatic to operational pressure
#   "stepped" - Stepped decrease with plateaus (more realistic)
LEACHING_MODE = "linear"

# LEACHING_DAYS: Duration of leaching phase in days
LEACHING_DAYS = 91

# LEACHING_DT_HOURS: Time step during leaching
LEACHING_DT_HOURS = 12

# STEPPED_N_STEPS: Number of pressure steps for "stepped" mode
STEPPED_N_STEPS = 6

# LEACHING_END_FRACTION: Fraction of lithostatic pressure for operational minimum
#   When USE_LEACHING = True: Leaching ends at this fraction
#   When USE_LEACHING = False: Equilibrium pressure for industry/transport derived from this
LEACHING_END_FRACTION = 0.40

# ── DEBRINING PHASE SETTINGS (only used when USE_LEACHING = True) ────────────
# After leaching, the cavern undergoes debrining where brine is displaced.
# During debrining, pressure stays constant at the leaching end pressure.

# DEBRINING_DAYS: Duration of debrining phase in days (constant pressure)
#   Set to 0 to skip debrining phase
DEBRINING_DAYS = 30

# ── TRANSITION SETTINGS ──────────────────────────────────────────────────────
# RAMP_UP_HOURS: Duration of smooth fade-in at the start of the operational phase.
#   The operational pressure scheme starts dampened and gradually reaches full
#   amplitude over this period, using a smooth cosine transition.
#   Prevents abrupt pressure changes after leaching/equilibrium.
#   Applies to both leaching and equilibrium initialization modes.
#   Set to 0 to disable (not recommended).
RAMP_UP_HOURS = 336 # 2 weeks

# ── PRESSURE SCENARIO ──────────────────────────────────────────────────────────
# PRESSURE_SCENARIO: Choose one of:
#   "industry"         - Sinusoidal supply variation (~12 cycles/year)
#   "transport"        - Trapezoidal 2-day cycle (nighttime high → daytime low)
#   "power_generation" - Abrupt withdrawal events with gradual re-pressurisation
#   "csv"              - Load pressure profile from CSV file
PRESSURE_SCENARIO = "power_generation"

# ── INDUSTRY SETTINGS (only used when PRESSURE_SCENARIO = "industry") ──────────
# Sinusoidal schedule. With leaching: oscillates around p_leach_end + P_AMPLITUDE_MPA.
# Without leaching: uses P_MEAN_MPA as centre.
P_MEAN_MPA = 15.0          # Mean pressure (only used if USE_LEACHING = False)
P_AMPLITUDE_MPA = 5      # Half peak-to-peak amplitude (MPa)

# ── TRANSPORT SETTINGS (only used when PRESSURE_SCENARIO = "transport") ─────────
# Two-day trapezoidal cycle: 8 h at high → ramp down → 16 h at low → ramp up → 8 h at high.
# Pressures are offsets above p_leach_end (or absolute if USE_LEACHING = False).
P_HIGH_OFFSET_MPA = 5.0    # High-pressure offset above p_leach_end (MPa)
P_LOW_OFFSET_MPA  = 1.5    # Low-pressure offset above p_leach_end (MPa)

# ── POWER GENERATION SETTINGS (only used when PRESSURE_SCENARIO = "power_generation") ──
# N_EVENTS abrupt withdrawal events over the operation period, each with a sharp
# 30-min drop, sustained low, and exponential re-pressurisation.
N_EVENTS = 20              # Number of withdrawal events
P_BASE_OFFSET_MPA = 10.0    # Resting pressure offset above p_leach_end (MPa)
RECOVERY_TAU_HOURS = 48.0   # Time constant for exponential re-pressurisation (hours)
                            #   4 h  = fast recovery (~95% in 12 h)
                            #  24 h  = gradual (~63% in 1 day, ~95% in 3 days)
                            #  48 h  = slow   (~63% in 2 days, ~95% in 6 days)
P_MIN_MPA = 10.0            # Absolute minimum cavern pressure (MPa)
                            # Prevents event stacking from driving pressure too low

# ── VARIABLE TIME-STEPPING (only used when PRESSURE_SCENARIO = "power_generation") ──
# Adaptive dt: small steps during sharp pressure drops, coarse steps elsewhere.
# MAX_DP_MPA controls refinement: steps are halved until |Δp| ≤ MAX_DP_MPA per step.
USE_VARIABLE_DT = True             # True = adaptive dt for power_generation
DT_FINE_HOURS   = 0.2              # Minimum dt (12 min) during sharp events
DT_COARSE_HOURS = 2.0              # Maximum dt away from events (same as dt_hours)
MAX_DP_MPA      = 0.2              # Max pressure change per step (MPa)

# ── SCHEDULE SETTINGS ──────────────────────────────────────────────────────────
# SCHEDULE_MODE: How to distribute cycles (used by "industry" and "transport"):
#   "stretch" - N_CYCLES spread evenly over OPERATION_DAYS
#   "repeat"  - Cycle pattern repeated with its native period
#   "direct"  - (CSV only) Use hourly CSV values directly
SCHEDULE_MODE = "stretch"

# OPERATION_DAYS: Total simulation duration in days (operation phase only)
OPERATION_DAYS = 365

# N_CYCLES: Number of pressure cycles (industry: sinusoidal; transport: 2-day cycles)
N_CYCLES = 180

# ── TIME STEP ──────────────────────────────────────────────────────────────────
dt_hours = 2

# ── CSV SETTINGS (only used when PRESSURE_SCENARIO = "csv") ────────────────────
CSV_FILE_PATH = "drukprofiel_zoutcaverne_2035_8760u.csv"
CSV_SHIFT_TO_LEACH_END = True  # Shift CSV so minimum = p_leach_end
P_EQUILIBRIUM_MPA = 15.0       # Equilibrium pressure (only used if USE_LEACHING = False)
RESCALE_PRESSURE = False          # Rescale CSV pressures to fit within a specific range 
RESCALE_MIN_MPA = 6.0
RESCALE_MAX_MPA = 20.0
RAMP_HOURS = 24.0

# ── MATERIAL MODEL ─────────────────────────────────────────────────────────────
USE_SCENARIO_B    = False    # False = Scenario A (CCC Zuidwending + Herminio calibration), True = Scenario B (calibrated)
USE_MUNSON_DAWSON = False    # False = SafeInCave model (Kelvin+Desai), True = Munson-Dawson model

# ── THERMAL MODEL ─────────────────────────────────────────────────────────────
USE_THERMAL = False
THERMAL_EXPANSION_COEFF = 40e-6
THERMAL_CONDUCTIVITY = 5.2
SPECIFIC_HEAT_CAPACITY = 837.0
T_SURFACE_C = 15.0
GEOTHERMAL_GRADIENT = 30.0
H_CONVECTION = 10.0
T_GAS_AMPLITUDE_C = 20.0

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                        END OF USER CONFIGURATION                              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


# ══════════════════════════════════════════════════════════════════════════════
# AUTOMATIC CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

VALID_SHAPES = ["regular", "tilted", "directcirculation", "asymmetric", "reversedcirculation", "fastleached", "tubefailure"]
VALID_SPECIAL_CAVERNS = ["A5"]
VALID_SIZES = [600, 1200]
VALID_SCENARIOS = ["industry", "transport", "power_generation", "csv"]
VALID_MODES_STANDARD = ["stretch", "repeat"]
VALID_MODES_CSV = ["stretch", "repeat", "direct"]
VALID_LEACHING_MODES = ["linear", "stepped"]

# Cavern z_max values (top of cavern elevation in meters)
Z_MAX_BY_CAVERN = {
    "regular600": 315.26,
    "tilted600": 345.67,
    "directcirculation600": 319.86,
    "asymmetric600": 338.89,
    "reversedcirculation600": 353.15,
    "regular1200": 393.21,
    "tilted1200": 430.78,
    "directcirculation1200": 402.21,
    "asymmetric1200": 422.76,
    "reversedcirculation1200": 445.06,
    "fastleached600": 378.19,
    "fastleached1200": 400.08,
    "tubefailure600": 420.90,
    "tubefailure1200": 444.53,
    # Zuidwending A5
    "A5": 515.0,
}

# Cavern heights (estimated)
CAVERN_HEIGHT_BY_TYPE = {
    "regular600": 150.0,
    "tilted600": 160.0,
    "directcirculation600": 145.0,
    "asymmetric600": 155.0,
    "reversedcirculation600": 180.0,
    "regular1200": 200.0,
    "tilted1200": 215.0,
    "directcirculation1200": 195.0,
    "asymmetric1200": 210.0,
    "reversedcirculation1200": 240.0,
    "fastleached600": 170.0,
    "fastleached1200": 215.0,
    "tubefailure600": 182.0,
    "tubefailure1200": 230.0,
    "A5": 370.0,
}

# Reference pressure at model top (z=660m) by cavern
# Standard caverns use size-based reference, A5 has its own
P_REF_BY_SIZE = {
    600: 17.5,   # MPa
    1200: 19.3,  # MPa
}
P_REF_BY_CAVERN = {
    "A5": 21.18,  # MPa (200m sand + 795m salt overburden)
}

# Grid folder mapping
GRID_FOLDERS = {
    "regular600": "cavern_regular_600_3D",
    "tilted600": "cavern_tilted_600_3D",
    "directcirculation600": "cavern_directcirculation_600_3D",
    "asymmetric600": "cavern_asymmetric_600_3D",
    "reversedcirculation600": "cavern_reversedcirculation_600_3D",
    "regular1200": "cavern_regular_1200_3D",
    "tilted1200": "cavern_tilted_1200_3D",
    "directcirculation1200": "cavern_directcirculation_1200_3D",
    "asymmetric1200": "cavern_asymmetric_1200_3D",
    "reversedcirculation1200": "cavern_reversedcirculation_1200_3D",
    "fastleached600": "cavern_fastleached_600_3D",
    "fastleached1200": "cavern_fastleached_1200_3D",
    "tubefailure600": "cavern_tubefailure_600_3D",
    "tubefailure1200": "cavern_tubefailure_1200_3D",
    "A5": "cavern_A5_3D",
}

Z_SURFACE = 660.0  # meters (top of domain)


def validate_configuration():
    """Validate user configuration and return derived values."""
    errors = []

    # Determine if this is a special cavern or standard shape+size
    is_special = CAVERN_TYPE in VALID_SPECIAL_CAVERNS

    if not is_special:
        if CAVERN_TYPE not in VALID_SHAPES:
            errors.append(f"CAVERN_TYPE '{CAVERN_TYPE}' invalid. Choose from: {VALID_SHAPES + VALID_SPECIAL_CAVERNS}")
        if CAVERN_SIZE not in VALID_SIZES:
            errors.append(f"CAVERN_SIZE '{CAVERN_SIZE}' invalid. Choose from: {VALID_SIZES}")

    if PRESSURE_SCENARIO not in VALID_SCENARIOS:
        errors.append(f"PRESSURE_SCENARIO '{PRESSURE_SCENARIO}' invalid. Choose from: {VALID_SCENARIOS}")

    if PRESSURE_SCENARIO == "csv":
        if SCHEDULE_MODE not in VALID_MODES_CSV:
            errors.append(f"SCHEDULE_MODE '{SCHEDULE_MODE}' invalid for CSV. Choose from: {VALID_MODES_CSV}")
    else:
        if SCHEDULE_MODE not in VALID_MODES_STANDARD:
            errors.append(f"SCHEDULE_MODE '{SCHEDULE_MODE}' invalid. Choose from: {VALID_MODES_STANDARD}")

    if OPERATION_DAYS <= 0:
        errors.append(f"OPERATION_DAYS must be positive, got {OPERATION_DAYS}")
    if N_CYCLES <= 0:
        errors.append(f"N_CYCLES must be positive, got {N_CYCLES}")
    if dt_hours <= 0:
        errors.append(f"dt_hours must be positive, got {dt_hours}")

    if USE_LEACHING:
        if LEACHING_MODE not in VALID_LEACHING_MODES:
            errors.append(f"LEACHING_MODE '{LEACHING_MODE}' invalid. Choose from: {VALID_LEACHING_MODES}")
        if LEACHING_DAYS <= 0:
            errors.append(f"LEACHING_DAYS must be positive, got {LEACHING_DAYS}")
        if LEACHING_DT_HOURS <= 0:
            errors.append(f"LEACHING_DT_HOURS must be positive, got {LEACHING_DT_HOURS}")
        if STEPPED_N_STEPS < 2:
            errors.append(f"STEPPED_N_STEPS must be at least 2, got {STEPPED_N_STEPS}")
        if LEACHING_END_FRACTION <= 0.0 or LEACHING_END_FRACTION >= 1.0:
            errors.append(f"LEACHING_END_FRACTION must be between 0 and 1 (exclusive), got {LEACHING_END_FRACTION}")
        if DEBRINING_DAYS < 0:
            errors.append(f"DEBRINING_DAYS must be non-negative, got {DEBRINING_DAYS}")
        if RAMP_UP_HOURS < 0:
            errors.append(f"RAMP_UP_HOURS must be non-negative, got {RAMP_UP_HOURS}")
    else:
        if EQUILIBRIUM_HOURS <= 0:
            errors.append(f"EQUILIBRIUM_HOURS must be positive, got {EQUILIBRIUM_HOURS}")
        if EQUILIBRIUM_DT_HOURS <= 0:
            errors.append(f"EQUILIBRIUM_DT_HOURS must be positive, got {EQUILIBRIUM_DT_HOURS}")

    if PRESSURE_SCENARIO == "industry":
        if P_AMPLITUDE_MPA <= 0:
            errors.append(f"P_AMPLITUDE_MPA must be positive, got {P_AMPLITUDE_MPA}")
    elif PRESSURE_SCENARIO == "transport":
        if P_HIGH_OFFSET_MPA <= P_LOW_OFFSET_MPA:
            errors.append(f"P_HIGH_OFFSET_MPA ({P_HIGH_OFFSET_MPA}) must be > P_LOW_OFFSET_MPA ({P_LOW_OFFSET_MPA})")
    elif PRESSURE_SCENARIO == "power_generation":
        if N_EVENTS < 1:
            errors.append(f"N_EVENTS must be >= 1, got {N_EVENTS}")
    elif PRESSURE_SCENARIO == "csv":
        if not os.path.isfile(CSV_FILE_PATH):
            errors.append(f"CSV_FILE_PATH not found: {CSV_FILE_PATH}")
        if not USE_LEACHING and RESCALE_PRESSURE and RESCALE_MIN_MPA >= RESCALE_MAX_MPA:
            errors.append(f"RESCALE_MIN_MPA ({RESCALE_MIN_MPA}) must be < RESCALE_MAX_MPA ({RESCALE_MAX_MPA})")

    if errors:
        raise ValueError("Configuration errors:\n  - " + "\n  - ".join(errors))

    # Derived values
    if is_special:
        cavern_key = CAVERN_TYPE
        grid_folder = GRID_FOLDERS[cavern_key]
        z_max = Z_MAX_BY_CAVERN[cavern_key]
        cavern_height = CAVERN_HEIGHT_BY_TYPE[cavern_key]
        p_ref_mpa = P_REF_BY_CAVERN.get(cavern_key, 17.5)
        cavern_label = f"{CAVERN_TYPE} (~965k m³)"
    else:
        cavern_key = f"{CAVERN_TYPE}{CAVERN_SIZE}"
        grid_folder = GRID_FOLDERS[cavern_key]
        z_max = Z_MAX_BY_CAVERN[cavern_key]
        cavern_height = CAVERN_HEIGHT_BY_TYPE[cavern_key]
        p_ref_mpa = P_REF_BY_SIZE[CAVERN_SIZE]
        cavern_label = f"{CAVERN_TYPE} ({CAVERN_SIZE}k m³)"

    z_center = z_max - cavern_height / 2.0

    return {
        "cavern_key": cavern_key,
        "cavern_label": cavern_label,
        "grid_folder": grid_folder,
        "z_max": z_max,
        "z_center": z_center,
        "cavern_height": cavern_height,
        "p_ref_mpa": p_ref_mpa,
        "is_special": is_special,
    }


def compute_lithostatic_pressure(z_center, p_ref_mpa, rho_salt, g):
    """Compute lithostatic pressure at cavern center."""
    depth = Z_SURFACE - z_center
    p_lithostatic = p_ref_mpa * ut.MPa + rho_salt * abs(g) * depth
    return p_lithostatic


def build_leaching_pressure_schedule(tc, *, p_start_pa, p_end_pa, mode, n_steps=6):
    """Build pressure schedule for leaching phase."""
    n_time_steps = int(math.floor(tc.t_final / tc.dt))
    t_vals = [k * tc.dt for k in range(n_time_steps + 1)]
    if abs(t_vals[-1] - tc.t_final) > 1e-12:
        t_vals.append(tc.t_final)

    if mode == "linear":
        p_vals = []
        for t in t_vals:
            frac = t / tc.t_final if tc.t_final > 0 else 1.0
            p = p_start_pa + frac * (p_end_pa - p_start_pa)
            p_vals.append(p)

    elif mode == "stepped":
        step_duration = tc.t_final / n_steps
        p_step_values = np.linspace(p_start_pa, p_end_pa, n_steps + 1)

        p_vals = []
        for t in t_vals:
            step_idx = min(int(t / step_duration), n_steps - 1)
            if t < tc.t_final:
                p = p_step_values[step_idx]
            else:
                p = p_end_pa
            p_vals.append(p)

    else:
        raise ValueError(f"Unknown leaching mode: {mode}")

    return t_vals, p_vals


def prepend_debrining(t_pressure, p_pressure, *, p_leach_end_pa, debrining_days):
    """Prepend debrining plateau to operational schedule.

    After leaching, pressure stays constant during debrining before the
    operational phase begins.
    """
    debrining_s = debrining_days * 24.0 * 3600.0

    if debrining_s <= 0.0:
        return t_pressure, p_pressure

    t_pre = [0.0, debrining_s]
    p_pre = [p_leach_end_pa, p_leach_end_pa]

    # Shift operation schedule and skip its first point (debrining end replaces it)
    t_shifted = [t + debrining_s for t in t_pressure[1:]]
    p_shifted = list(p_pressure[1:])

    return t_pre + t_shifted, p_pre + p_shifted


def apply_fade_in(t_pressure, p_pressure, *, p_start_pa, fade_in_hours):
    """Apply smooth fade-in to the beginning of a pressure schedule.

    Blends from p_start_pa to the actual schedule values using a smooth
    cosine transition over fade_in_hours. This dampens initial pressure
    swings and prevents abrupt pressure changes.

    p(t) = (1 - alpha) * p_start + alpha * p_operational(t)
    where alpha = 0.5 * (1 - cos(pi * t / T_fade))
    """
    if fade_in_hours <= 0.0:
        return

    fade_in_s = fade_in_hours * 3600.0

    for i in range(len(t_pressure)):
        t = t_pressure[i]
        if t >= fade_in_s:
            break
        alpha = 0.5 * (1.0 - math.cos(math.pi * t / fade_in_s))
        p_pressure[i] = (1.0 - alpha) * p_start_pa + alpha * p_pressure[i]


# ══════════════════════════════════════════════════════════════════════════════
# CSV PRESSURE PROFILE FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

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


def read_pressure_csv(csv_file: str):
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
    old_min = pressures_mpa.min()
    old_max = pressures_mpa.max()

    if old_max - old_min < 1e-9:
        return np.full_like(pressures_mpa, (new_min + new_max) / 2.0)

    fraction = (pressures_mpa - old_min) / (old_max - old_min)
    return new_min + fraction * (new_max - new_min)


def build_csv_pressure_schedule(tc, csv_file, *, days, mode, total_cycles=1,
                                rescale=False, rescale_min=None, rescale_max=None,
                                resample_at_dt=True):
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
# PRESSURE SCHEDULE BUILDERS
# ══════════════════════════════════════════════════════════════════════════════

DAY_H = 24.0


def _sample_at_dt(tc, t_end=None):
    t_end = tc.t_final if t_end is None else t_end
    n_steps = int(math.floor(t_end / tc.dt))
    t_vals = [k * tc.dt for k in range(n_steps + 1)]
    if abs(t_vals[-1] - t_end) > 1e-12:
        t_vals.append(t_end)
    return t_vals


# ── Variable-dt time controller (for power_generation) ────────────────────────

class TimeControllerFromList:
    """
    Time controller that steps through a pre-computed list of times.
    Compatible with sf.Simulator_M (exposes t, dt, t_final, time_unit,
    time_conversion, step_counter, keep_looping, advance_time).
    """
    def __init__(self, time_list_seconds):
        self.time_list = np.asarray(time_list_seconds, dtype=float)
        if self.time_list.ndim != 1 or self.time_list.size < 2:
            raise ValueError("time_list_seconds must have at least 2 entries.")
        if not np.all(np.diff(self.time_list) > 0):
            raise ValueError("time_list_seconds must be strictly increasing.")

        self.t_initial = float(self.time_list[0])
        self.t_final = float(self.time_list[-1])
        self.t = float(self.time_list[0])
        self.step_counter = 0
        self.dt = float(self.time_list[1] - self.time_list[0])
        # Required by Simulator_M / ScreenPrinter
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
    """
    Build a variable time grid so that |p(t+dt) - p(t)| <= dp_max_pa,
    with dt clamped in [dt_min_s, dt_max_s]. All arguments in seconds/Pascals.
    """
    t = 0.0
    times = [0.0]
    p_prev = float(p_of_t(0.0))

    max_steps = int(np.ceil(t_final_s / dt_min_s)) + 50
    for _ in range(max_steps):
        if t >= t_final_s - 1e-12:
            break

        dt = dt_max_s
        while True:
            t_try = min(t + dt, t_final_s)
            p_try = float(p_of_t(t_try))
            if abs(p_try - p_prev) <= dp_max_pa or dt <= dt_min_s + 1e-12:
                t = t_try
                p_prev = p_try
                times.append(t)
                break
            dt *= 0.5
            if dt < dt_min_s:
                dt = dt_min_s

    if abs(times[-1] - t_final_s) > 1e-9:
        times.append(t_final_s)
    return times


def build_sinus_pressure_schedule(tc, *, p_mean, p_ampl, period_hours, phase_hours=0.0,
                                  clamp_min=None, clamp_max=None):
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


def build_power_generation_schedule(tc, *, p_base_pa, n_events, operation_days,
                                    recovery_tau_hours=48.0, p_min_pa=None, seed=42):
    """
    N_EVENTS abrupt withdrawal events spread over the operation period.
    Each event: sharp 30-min drop → sustained low (2–5 h) → exponential recovery.
    recovery_tau_hours controls how gradually pressure returns to base.
    p_min_pa: minimum allowed pressure (Pa) — prevents event stacking from
              driving pressure unrealistically low.
    Reproducible via seed.  Returns (t_vals_s, p_vals_Pa).
    """
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
        duration = rng.uniform(2.0, 5.0)   # sustained low: 2–5 hours
        depth    = rng.uniform(6.5, 12.5)   # pressure drop: 6.5–12.5 MPa
        for i, t in enumerate(t_h):
            if t < t_start_h:
                continue
            dt_ev = t - t_start_h
            if dt_ev < 0.5:
                drop = depth * (dt_ev / 0.5)
            elif dt_ev < 0.5 + duration:
                drop = depth
            else:
                drop = depth * math.exp(-(dt_ev - 0.5 - duration) / tau)
                if drop < 0.05:
                    break
            p_mpa[i] = min(p_mpa[i], p_base_mpa - drop)

    # Clamp to minimum pressure to prevent event stacking from going too low
    if p_min_mpa is not None:
        p_mpa = np.maximum(p_mpa, p_min_mpa)

    p_vals = (p_mpa * ut.MPa).tolist()
    return t_vals_s, p_vals


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
    config = validate_configuration()

    salt_density = 2200
    g = -9.81

    # Compute lithostatic pressure at cavern center
    p_lithostatic = compute_lithostatic_pressure(
        config["z_center"],
        config["p_ref_mpa"],
        salt_density,
        g
    )
    p_lithostatic_mpa = p_lithostatic / ut.MPa

    # Determine pressures based on initialization mode
    if USE_LEACHING:
        p_leach_end_mpa = LEACHING_END_FRACTION * p_lithostatic_mpa
        p_leach_end = p_leach_end_mpa * ut.MPa

        if PRESSURE_SCENARIO == "industry":
            p_op_min_mpa = p_leach_end_mpa
            p_op_mean_mpa = p_leach_end_mpa + P_AMPLITUDE_MPA
            p_op_max_mpa = p_leach_end_mpa + 2 * P_AMPLITUDE_MPA
        elif PRESSURE_SCENARIO == "transport":
            p_op_min_mpa = p_leach_end_mpa + P_LOW_OFFSET_MPA
            p_op_max_mpa = p_leach_end_mpa + P_HIGH_OFFSET_MPA
        else:
            p_op_min_mpa = p_leach_end_mpa
    else:
        # Equilibrium mode: use direct pressure settings
        if PRESSURE_SCENARIO == "industry":
            p_gas_mpa = P_MEAN_MPA
        elif PRESSURE_SCENARIO == "transport":
            p_gas_mpa = P_MEAN_MPA
        elif PRESSURE_SCENARIO == "csv":
            p_gas_mpa = P_EQUILIBRIUM_MPA
        else:
            p_gas_mpa = P_EQUILIBRIUM_MPA
        p_gas = p_gas_mpa * ut.MPa

    # Print configuration
    if MPI.COMM_WORLD.rank == 0:
        print("=" * 70)
        print(f"SIMULATION CONFIGURATION ({'LEACHING' if USE_LEACHING else 'EQUILIBRIUM'} MODE)")
        print("=" * 70)
        print(f"  Cavern:           {config['cavern_label']}")
        print(f"  Cavern z_max:     {config['z_max']:.2f} m")
        print(f"  Cavern z_center:  {config['z_center']:.2f} m")
        print(f"  P_ref (at z=660): {config['p_ref_mpa']:.2f} MPa")
        print(f"  P_lithostatic:    {p_lithostatic_mpa:.2f} MPa (at cavern center)")
        print("-" * 70)

        if USE_LEACHING:
            print("  LEACHING PHASE:")
            print(f"    Mode:           {LEACHING_MODE}")
            print(f"    Duration:       {LEACHING_DAYS} days")
            print(f"    P_start:        {p_lithostatic_mpa:.2f} MPa (lithostatic)")
            print(f"    P_end:          {p_leach_end_mpa:.2f} MPa ({LEACHING_END_FRACTION*100:.0f}%)")
            if DEBRINING_DAYS > 0 or RAMP_UP_HOURS > 0:
                print("  DEBRINING/RAMP-UP:")
                print(f"    Debrining:      {DEBRINING_DAYS} days at {p_leach_end_mpa:.2f} MPa")
                print(f"    Ramp-up:        {RAMP_UP_HOURS} hours")
        else:
            print("  EQUILIBRIUM PHASE:")
            print(f"    Duration:       {EQUILIBRIUM_HOURS} hours")
            print(f"    P_gas:          {p_gas_mpa:.2f} MPa")

        print("-" * 70)
        print("  OPERATION PHASE:")
        print(f"    Scenario:       {PRESSURE_SCENARIO}")
        print(f"    Mode:           {SCHEDULE_MODE}")
        print(f"    Days:           {OPERATION_DAYS}")
        print(f"    Cycles:         {N_CYCLES}")
        print(f"    dt:             {dt_hours} hours")
        print(f"    Material:       Scenario {'B' if USE_SCENARIO_B else 'A'} / {'Munson-Dawson' if USE_MUNSON_DAWSON else 'SafeInCave'}")
        print(f"    Thermal:        {'enabled' if USE_THERMAL else 'disabled'}")
        print("-" * 70)
        print(f"  Grid:             {config['grid_folder']}")
        print("=" * 70)

    # Load grid
    grid_path = os.path.join("..", "..", "..", "..", "grids", config["grid_folder"])
    grid = sf.GridHandlerGMSH("geom", grid_path)

    z_max = config["z_max"]
    p_ref = config["p_ref_mpa"] * ut.MPa

    # Output folder
    _scen_tag  = "B" if USE_SCENARIO_B else "A"
    _model_tag = "MD" if USE_MUNSON_DAWSON else "SIC"
    if USE_LEACHING:
        output_folder = os.path.join(
            "output",
            f"case_leaching_{LEACHING_MODE}_{PRESSURE_SCENARIO}({N_CYCLES})_{OPERATION_DAYS}days_S{_scen_tag}_{_model_tag}_{config['cavern_key']}"
        )
    else:
        output_folder = os.path.join(
            "output",
            f"case_{PRESSURE_SCENARIO}({N_CYCLES})_{OPERATION_DAYS}days_S{_scen_tag}_{_model_tag}_{config['cavern_key']}"
        )

    side_burden = p_ref
    over_burden = p_ref

    # Define momentum equation
    mom_eq = LinearMomentumMod(grid, theta=0.5)

    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("cg")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-10, max_it=100)
    mom_eq.set_solver(mom_solver)

    # Material
    mat = sf.Material(mom_eq.n_elems)
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    sec_per_year = 365.25 * 24 * 3600

    # -- Elastic spring (same for both scenarios) --
    E0 = 20.425 * ut.GPa * to.ones(mom_eq.n_elems)
    nu0 = 0.25 * to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")
    mat.add_to_elastic(spring_0)

    if not USE_SCENARIO_B:
        # ── Scenario A: CCC Zuidwending ───────────────────────────────────────
        if not USE_MUNSON_DAWSON:
            # SafeInCave model: Kelvin + DislocationCreep
            eta = 2.5e5 * to.ones(mom_eq.n_elems)
            E1  = 42.0 * ut.GPa * to.ones(mom_eq.n_elems)
            nu1 = 0.32 * to.ones(mom_eq.n_elems)
            kelvin  = sf.Viscoelastic(eta, E1, nu1, "kelvin")
            ndc = 4.6
            A_dc = (40.0 * (1e-6)**ndc / sec_per_year) * to.ones(mom_eq.n_elems)
            Q_dc = (6495.0 * 8.32) * to.ones(mom_eq.n_elems)
            n_dc = ndc * to.ones(mom_eq.n_elems)
            creep_0 = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation")
            mat.add_to_non_elastic(kelvin)
            mat.add_to_non_elastic(creep_0)
        else:
            # Munson-Dawson model (contains its own steady-state disloc)
            nmd   = 4.99
            A_md  = (18.31 * (1e-6)**nmd / sec_per_year) * to.ones(mom_eq.n_elems)
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
    else:
        # ── Scenario B: calibrated against TCC-1 ──────────────────────────────
        if not USE_MUNSON_DAWSON:
            # SafeInCave model: Kelvin + DislocationCreep
            eta = 5.0e12 * to.ones(mom_eq.n_elems)
            E1  = 1.5 * ut.GPa * to.ones(mom_eq.n_elems)
            nu1 = 0.25 * to.ones(mom_eq.n_elems)
            kelvin  = sf.Viscoelastic(eta, E1, nu1, "kelvin")
            ndc = 5.0
            A_dc = (40.0 * (1e-6)**ndc / sec_per_year) * to.ones(mom_eq.n_elems)
            Q_dc = (6252.0 * 8.32) * to.ones(mom_eq.n_elems)
            n_dc = ndc * to.ones(mom_eq.n_elems)
            creep_0 = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation")
            mat.add_to_non_elastic(kelvin)
            mat.add_to_non_elastic(creep_0)
        else:
            # Munson-Dawson model (shared A/n/Q with SIC Scenario B)
            nmd   = 5.0
            A_md  = (40.0 * (1e-6)**nmd / sec_per_year) * to.ones(mom_eq.n_elems)
            Q_md  = (6252.0 * 8.32) * to.ones(mom_eq.n_elems)
            n_md  = nmd  * to.ones(mom_eq.n_elems)
            K0_md = 0.60   * to.ones(mom_eq.n_elems)
            c_md  = 9.02e-3 * to.ones(mom_eq.n_elems)
            m_md  = 1.1    * to.ones(mom_eq.n_elems)
            aw_md = -17.0  * to.ones(mom_eq.n_elems)
            bw_md = -7.738 * to.ones(mom_eq.n_elems)
            d_md  = 0.25   * to.ones(mom_eq.n_elems)
            mu_md = E0 / (2.0 * (1.0 + nu0))
            md = sf.MunsonDawsonCreep(A=A_md, Q=Q_md, n=n_md,
                                      K0=K0_md, c=c_md, m=m_md,
                                      alpha_w=aw_md, beta_w=bw_md, delta=d_md,
                                      mu=mu_md, name="munson_dawson")
            mat.add_to_non_elastic(md)

    # -- Pressure-solution creep (same for both scenarios) --
    A_ps = (14176.0 * 1e-9 / 1e6 / sec_per_year) * to.ones(mom_eq.n_elems)
    d_ps = 5.25e-3 * to.ones(mom_eq.n_elems)
    Q_ps = (3252.0 * 8.32) * to.ones(mom_eq.n_elems)
    creep_pressure = sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "creep_pressure")
    mat.add_to_non_elastic(creep_pressure)

    if USE_THERMAL:
        alpha_th = THERMAL_EXPANSION_COEFF * to.ones(mom_eq.n_elems, dtype=to.float64)
        thermo = sf.Thermoelastic(alpha_th, "thermo")
        mat.add_to_thermoelastic(thermo)

    mom_eq.set_material(mat)

    g_vec = [0.0, 0.0, g]
    mom_eq.build_body_force(g_vec)

    # Temperature
    T_surface_K = T_SURFACE_C + 273.15
    dTdz = GEOTHERMAL_GRADIENT / 1000.0

    if USE_THERMAL:
        cell_centroids = grid.mesh.geometry.x[grid.mesh.topology.connectivity(3, 0).array].reshape(-1, 4, 3).mean(axis=1)
        z_coords = cell_centroids[:, 2]
        T0_field = to.tensor(T_surface_K + dTdz * (Z_SURFACE - z_coords), dtype=to.float64)
    else:
        T0_field = 298 * to.ones(mom_eq.n_elems)

    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    T_cavern_init_K = T_surface_K + dTdz * (Z_SURFACE - z_max)

    gas_density = 0.089

    # ===== INITIALIZATION PHASE (LEACHING OR EQUILIBRIUM) =====
    if USE_LEACHING:
        tc_init = sf.TimeController(
            dt=LEACHING_DT_HOURS,
            initial_time=0.0,
            final_time=LEACHING_DAYS * 24.0,
            time_unit="hour"
        )

        t_init, p_init = build_leaching_pressure_schedule(
            tc_init,
            p_start_pa=p_lithostatic,
            p_end_pa=p_leach_end,
            mode=LEACHING_MODE,
            n_steps=STEPPED_N_STEPS
        )
        init_phase_name = "leaching"
        save_interval_init = max(1, int(72 / LEACHING_DT_HOURS))
    else:
        tc_init = sf.TimeController(
            dt=EQUILIBRIUM_DT_HOURS,
            initial_time=0.0,
            final_time=EQUILIBRIUM_HOURS,
            time_unit="hour"
        )
        t_init = [0.0, tc_init.t_final]
        p_init = [p_gas, p_gas]
        init_phase_name = "equilibrium"
        save_interval_init = 1

    # Initialization BCs
    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_init.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_init.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_init.t_final])

    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0,
                              [side_burden, side_burden],
                              [0.0, tc_init.t_final], g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0,
                               [side_burden, side_burden],
                               [0.0, tc_init.t_final], g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                             [over_burden, over_burden],
                             [0.0, tc_init.t_final], g=g_vec[2])
    bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, z_max,
                                p_init, t_init, g=g_vec[2])

    bc_init = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
        bc_init.add_boundary_condition(bc)

    mom_eq.set_boundary_conditions(bc_init)

    output_folder_init = os.path.join(output_folder, init_phase_name)
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[{init_phase_name.upper()}] Output: {output_folder_init}")

    output_mom_init = SparseSaveFields(mom_eq, interval=save_interval_init)
    output_mom_init.set_output_folder(output_folder_init)
    output_mom_init.add_output_field("u", "Displacement (m)")
    output_mom_init.add_output_field("eps_tot", "Total strain (-)")
    output_mom_init.add_output_field("sig", "Stress (Pa)")
    output_mom_init.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom_init.add_output_field("q_elems", "Von Mises stress (Pa)")

    os.makedirs(output_folder, exist_ok=True)

    sim_init = sf.Simulator_M(mom_eq, tc_init, [output_mom_init], True)
    sim_init.run()

    if MPI.COMM_WORLD.rank == 0:
        print(f"[{init_phase_name.upper()}] Complete.")

    # ===== OPERATION PHASE =====
    if not USE_MUNSON_DAWSON:
        # SafeInCave model: add Desai viscoplastic element (scenario-dependent params)
        if not USE_SCENARIO_B:
            # Scenario A (CCC Zuidwending)
            mu_1    = 6.89e-12 * to.ones(mom_eq.n_elems)
            N_1     = 3.0 * to.ones(mom_eq.n_elems)
            a_1     = 1.80e-5 * to.ones(mom_eq.n_elems)
            eta_vp  = 0.82 * to.ones(mom_eq.n_elems)
            alpha_0 = 2.0e-3 * to.ones(mom_eq.n_elems)
        else:
            # Scenario B (calibrated)
            mu_1    = 1.0e-15 * to.ones(mom_eq.n_elems)
            N_1     = 2.0 * to.ones(mom_eq.n_elems)
            a_1     = 1.0e-3 * to.ones(mom_eq.n_elems)
            eta_vp  = 1.2 * to.ones(mom_eq.n_elems)
            alpha_0 = 5.0e-3 * to.ones(mom_eq.n_elems)
        # Fixed Desai shape parameters (both scenarios)
        n_desai  = 3.0 * to.ones(mom_eq.n_elems)
        beta_1   = 0.0048 * to.ones(mom_eq.n_elems)
        beta     = 0.995 * to.ones(mom_eq.n_elems)
        m_desai  = -0.5 * to.ones(mom_eq.n_elems)
        gamma    = 0.095 * to.ones(mom_eq.n_elems)
        sigma_t  = 5.0 * to.ones(mom_eq.n_elems)
        desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta_vp, n_desai, beta_1, beta,
                                     m_desai, gamma, sigma_t, alpha_0, "desai")

        stress_to = ut.numpy2torch(mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3)))
        desai.compute_initial_hardening(stress_to, Fvp_0=0.0)

        mat.add_to_non_elastic(desai)
        mom_eq.set_material(mat)
        mom_eq.expect_vp_state = True
    else:
        mom_eq.expect_vp_state = False

    # Build cycling schedule using a temporary tc for just the operational cycling period
    tc_cycling = sf.TimeController(
        dt=dt_hours,
        initial_time=0.0,
        final_time=OPERATION_DAYS * 24.0,
        time_unit="hour"
    )

    # Build operation pressure schedule
    if PRESSURE_SCENARIO == "industry":
        if USE_LEACHING:
            p_mean = p_op_mean_mpa * ut.MPa
        else:
            p_mean = P_MEAN_MPA * ut.MPa
        p_ampl = P_AMPLITUDE_MPA * ut.MPa

        t_pressure, p_pressure = build_sinus_schedule_multi(
            tc_cycling,
            p_mean=p_mean, p_ampl=p_ampl,
            days=OPERATION_DAYS, mode=SCHEDULE_MODE,
            daily_period_hours=24.0,
            total_cycles=N_CYCLES,
            clamp_min=None, clamp_max=None
        )

    elif PRESSURE_SCENARIO == "transport":
        if USE_LEACHING:
            p_high = p_leach_end_mpa + P_HIGH_OFFSET_MPA
            p_low  = p_leach_end_mpa + P_LOW_OFFSET_MPA
        else:
            p_high = P_MEAN_MPA + P_HIGH_OFFSET_MPA
            p_low  = P_MEAN_MPA + P_LOW_OFFSET_MPA
        # Two-day trapezoidal cycle: 8 h high → ramp → 16 h low → ramp → 8 h high
        base_times_h    = [0.0, 8.0, 12.0, 28.0, 32.0, 48.0]
        base_pressures_MPa = [p_high, p_high, p_low, p_low, p_high, p_high]

        t_pressure, p_pressure = build_linear_schedule_multi(
            tc_cycling,
            base_times_h, base_pressures_MPa,
            days=OPERATION_DAYS,
            mode=SCHEDULE_MODE,
            resample_at_dt=True,
            total_cycles=N_CYCLES,
        )

    elif PRESSURE_SCENARIO == "power_generation":
        if USE_LEACHING:
            p_base_pa = (p_leach_end_mpa + P_BASE_OFFSET_MPA) * ut.MPa
        else:
            p_base_pa = (P_EQUILIBRIUM_MPA + P_BASE_OFFSET_MPA) * ut.MPa

        t_pressure, p_pressure = build_power_generation_schedule(
            tc_cycling,
            p_base_pa=p_base_pa,
            n_events=N_EVENTS,
            operation_days=OPERATION_DAYS,
            recovery_tau_hours=RECOVERY_TAU_HOURS,
            p_min_pa=P_MIN_MPA * ut.MPa,
        )

    elif PRESSURE_SCENARIO == "csv":
        t_pressure, p_pressure = build_csv_pressure_schedule(
            tc_cycling,
            csv_file=CSV_FILE_PATH,
            days=OPERATION_DAYS,
            mode=SCHEDULE_MODE,
            total_cycles=N_CYCLES,
            rescale=RESCALE_PRESSURE if not USE_LEACHING else False,
            rescale_min=RESCALE_MIN_MPA,
            rescale_max=RESCALE_MAX_MPA,
            resample_at_dt=True
        )

        if USE_LEACHING and CSV_SHIFT_TO_LEACH_END:
            p_array = np.array(p_pressure)
            csv_min_pa = p_array.min()
            shift_pa = p_leach_end - csv_min_pa
            p_pressure = [p + shift_pa for p in p_pressure]
            if MPI.COMM_WORLD.rank == 0:
                new_min = min(p_pressure) / ut.MPa
                new_max = max(p_pressure) / ut.MPa
                print(f"[CSV] Shifted profile: new range [{new_min:.2f}, {new_max:.2f}] MPa")

        # Only apply CSV startup ramp if fade-in is not active
        if RAMP_UP_HOURS <= 0:
            apply_startup_ramp(
                t_pressure, p_pressure,
                p_start_pa=p_leach_end if USE_LEACHING else p_gas,
                ramp_hours=RAMP_HOURS,
                dt_hours=dt_hours
            )

    else:
        raise ValueError(f"Unknown PRESSURE_SCENARIO: {PRESSURE_SCENARIO}")

    # Apply smooth fade-in to dampen initial pressure swings
    if RAMP_UP_HOURS > 0:
        p_fade_start = p_leach_end if USE_LEACHING else p_gas
        apply_fade_in(t_pressure, p_pressure,
                      p_start_pa=p_fade_start,
                      fade_in_hours=RAMP_UP_HOURS)
        if MPI.COMM_WORLD.rank == 0:
            print(f"[TRANSITION] Fade-in: {RAMP_UP_HOURS:.1f} hours from {p_fade_start/ut.MPa:.2f} MPa")

    # Prepend debrining plateau when using leaching
    extra_hours = 0.0
    if USE_LEACHING and DEBRINING_DAYS > 0:
        extra_hours = DEBRINING_DAYS * 24.0
        t_pressure, p_pressure = prepend_debrining(
            t_pressure, p_pressure,
            p_leach_end_pa=p_leach_end,
            debrining_days=DEBRINING_DAYS
        )
        if MPI.COMM_WORLD.rank == 0:
            print(f"[TRANSITION] Debrining: {DEBRINING_DAYS} days at {p_leach_end_mpa:.2f} MPa")

    # Create tc_operation with full duration (including debrining)
    total_operation_hours = OPERATION_DAYS * 24.0 + extra_hours

    if PRESSURE_SCENARIO == "power_generation" and USE_VARIABLE_DT:
        # Adaptive time-stepping: fine dt during sharp events, coarse elsewhere
        t_arr = np.array(t_pressure, dtype=float)
        p_arr = np.array(p_pressure, dtype=float)
        p_of_t = lambda t: float(np.interp(t, t_arr, p_arr))

        time_list = build_time_list_by_dp_limit(
            total_operation_hours * ut.hour,
            p_of_t,
            dt_min_s=DT_FINE_HOURS * ut.hour,
            dt_max_s=DT_COARSE_HOURS * ut.hour,
            dp_max_pa=MAX_DP_MPA * ut.MPa,
        )
        tc_operation = TimeControllerFromList(time_list)
        if MPI.COMM_WORLD.rank == 0:
            n_steps = len(time_list) - 1
            n_coarse = int(total_operation_hours / DT_COARSE_HOURS)
            print(f"[VARIABLE DT] {n_steps} steps (vs {n_coarse} at fixed dt={DT_COARSE_HOURS}h)")
            dts = np.diff(time_list) / ut.hour
            print(f"[VARIABLE DT] dt range: {dts.min():.3f} – {dts.max():.3f} hours")
    else:
        tc_operation = sf.TimeController(
            dt=dt_hours,
            initial_time=0.0,
            final_time=total_operation_hours,
            time_unit="hour"
        )

    # Operation BCs
    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_operation.t_final])

    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0,
                              [side_burden, side_burden],
                              [0.0, tc_operation.t_final], g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0,
                               [side_burden, side_burden],
                               [0.0, tc_operation.t_final], g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                             [over_burden, over_burden],
                             [0.0, tc_operation.t_final], g=g_vec[2])
    bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, z_max,
                                p_pressure, t_pressure, g=g_vec[2])

    bc_operation = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
        bc_operation.add_boundary_condition(bc)

    mom_eq.set_boundary_conditions(bc_operation)

    output_folder_operation = os.path.join(output_folder, "operation")
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[OPERATION] Output: {output_folder_operation}")

    # Save pressure schedule
    pressure_data = {
        "cavern_key": config["cavern_key"],
        "cavern_label": config["cavern_label"],
        "use_leaching": USE_LEACHING,
        "debrining_days": DEBRINING_DAYS if USE_LEACHING else 0,
        "ramp_up_hours": RAMP_UP_HOURS if USE_LEACHING else 0,
        "scenario": PRESSURE_SCENARIO,
        "mode": SCHEDULE_MODE,
        "n_cycles": N_CYCLES,
        "operation_days": OPERATION_DAYS,
        "dt_hours": dt_hours,
        "scenario": "B" if USE_SCENARIO_B else "A",
        "model": "munson_dawson" if USE_MUNSON_DAWSON else "safeincave",
        "p_lithostatic_mpa": p_lithostatic_mpa,
        "units": {"t_raw": "s", "p_raw": "Pa", "t": "hour", "p": "MPa"},
        "t_values_s": [float(t) for t in t_pressure],
        "p_values_Pa": [float(p) for p in p_pressure],
        "t_hours": [float(t / ut.hour) for t in t_pressure],
        "p_MPa": [float(p / ut.MPa) for p in p_pressure],
    }

    with open(os.path.join(output_folder, "pressure_schedule.json"), 'w') as f:
        json.dump(pressure_data, f, indent=2)

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

    if USE_THERMAL:
        heat_eq = sf.HeatDiffusion(grid)

        heat_solver = PETSc.KSP().create(grid.mesh.comm)
        heat_solver.setType("cg")
        heat_solver.getPC().setType("asm")
        heat_solver.setTolerances(rtol=1e-12, max_it=100)
        heat_eq.set_solver(heat_solver)

        cp = SPECIFIC_HEAT_CAPACITY * to.ones(heat_eq.n_elems, dtype=to.float64)
        mat.set_specific_heat_capacity(cp)

        k_th = THERMAL_CONDUCTIVITY * to.ones(heat_eq.n_elems, dtype=to.float64)
        mat.set_thermal_conductivity(k_th)

        heat_eq.set_material(mat)

        node_coords = grid.mesh.geometry.x
        z_nodes = node_coords[:, 2]
        T0_nodes = T_surface_K + dTdz * (Z_SURFACE - z_nodes)
        T0_nodes_tensor = to.tensor(T0_nodes, dtype=to.float64)
        heat_eq.set_initial_T(T0_nodes_tensor)

        heat_time_values = [tc_operation.t_initial, tc_operation.t_final]
        nt_heat = len(heat_time_values)

        heat_bc_handler = heatBC.BcHandler(heat_eq)

        bc_heat_west = heatBC.NeumannBC("West", nt_heat * [0.0], heat_time_values)
        bc_heat_east = heatBC.NeumannBC("East", nt_heat * [0.0], heat_time_values)
        bc_heat_south = heatBC.NeumannBC("South", nt_heat * [0.0], heat_time_values)
        bc_heat_north = heatBC.NeumannBC("North", nt_heat * [0.0], heat_time_values)
        bc_heat_top = heatBC.DirichletBC("Top", nt_heat * [T_surface_K], heat_time_values)
        bc_heat_bottom = heatBC.NeumannBC("Bottom", nt_heat * [dTdz], heat_time_values)

        heat_bc_handler.add_boundary_condition(bc_heat_west)
        heat_bc_handler.add_boundary_condition(bc_heat_east)
        heat_bc_handler.add_boundary_condition(bc_heat_south)
        heat_bc_handler.add_boundary_condition(bc_heat_north)
        heat_bc_handler.add_boundary_condition(bc_heat_top)
        heat_bc_handler.add_boundary_condition(bc_heat_bottom)

        n_t_steps = len(t_pressure)
        T_gas_values = []
        for i, t in enumerate(t_pressure):
            period_s = tc_operation.t_final / N_CYCLES if SCHEDULE_MODE == "stretch" else 24.0 * ut.hour
            phase = 2.0 * math.pi * t / period_s
            if PRESSURE_SCENARIO == "industry":
                T_gas = T_cavern_init_K - T_GAS_AMPLITUDE_C * math.sin(phase)
            else:
                T_gas = T_cavern_init_K + T_GAS_AMPLITUDE_C * math.sin(phase)
            T_gas_values.append(T_gas)

        bc_heat_cavern = heatBC.RobinBC("Cavern", T_gas_values, H_CONVECTION, t_pressure)
        heat_bc_handler.add_boundary_condition(bc_heat_cavern)

        heat_eq.set_boundary_conditions(heat_bc_handler)

        output_heat_op = SparseSaveFields(heat_eq, interval=15)
        output_heat_op.set_output_folder(output_folder_operation)
        output_heat_op.add_output_field("T", "Temperature (K)")
        outputs_op.append(output_heat_op)

        if MPI.COMM_WORLD.rank == 0:
            print(f"[THERMAL] Enabled")

        sim_op = sf.Simulator_TM(mom_eq, heat_eq, tc_operation, outputs_op, False)
        sim_op.run()
    else:
        sim_op = sf.Simulator_M(mom_eq, tc_operation, outputs_op, False)
        sim_op.run()

    if MPI.COMM_WORLD.rank == 0:
        print("[OPERATION] Complete.")
        print("=" * 70)
        print("SIMULATION FINISHED")
        extra_d = DEBRINING_DAYS if USE_LEACHING else 0.0
        total_days = (LEACHING_DAYS if USE_LEACHING else EQUILIBRIUM_HOURS/24) + extra_d + OPERATION_DAYS
        print(f"  Total simulated time: {total_days:.1f} days")
        print(f"  Output folder: {output_folder}")
        print("=" * 70)


if __name__ == '__main__':
    main()
