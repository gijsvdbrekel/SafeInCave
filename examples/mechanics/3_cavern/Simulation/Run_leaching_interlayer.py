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
# ║  This script supports caverns with interlayer heterogeneity.                  ║
# ║  - For homogeneous salt: use CAVERN_TYPE = "nointerlayer"                     ║
# ║  - For heterogeneous with interlayers: use CAVERN_TYPE = "interlayer"         ║
# ║  - For irregular shapes with interlayers:                                     ║
# ║      * "bulbous_ledges_600" / "bulbous_ledges_1200" - 3 interlayers, bulging  ║
# ║      * "asymmetric_shelf_600" / "asymmetric_shelf_1200" - 2 interlayers, tilt ║
# ║      * "vertical_intrusion_600" / "vertical_intrusion_1200" - vertical effect ║
# ║                                                                               ║
# ║  Interlayers are purely elastic (no creep). You can choose the material       ║
# ║  for each interlayer independently:                                           ║
# ║  - Anhydrite: E = 61.5 GPa, nu = 0.32                                        ║
# ║  - Mudstone:  E = 19.33 GPa, nu = 0.223                                      ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── CAVERN SELECTION ───────────────────────────────────────────────────────────
# CAVERN_TYPE: Choose one of:
#   "nointerlayer"           - Homogeneous salt, same cavern shape but no interlayer volumes
#   "interlayer"             - Standard heterogeneous with 2 interlayers (5 material regions)
#   "bulbous_ledges_600"     - Irregular bulging shape with 3 interlayers (600,000 m³)
#   "bulbous_ledges_1200"    - Irregular bulging shape with 3 interlayers (1,200,000 m³)
#   "asymmetric_shelf_600"   - Asymmetric tilted shape with 2 interlayers (600,000 m³)
#   "asymmetric_shelf_1200"  - Asymmetric tilted shape with 2 interlayers (1,200,000 m³)
#   "vertical_intrusion_600" - Shape with vertical interlayer effect (600,000 m³)
#   "vertical_intrusion_1200"- Shape with vertical interlayer effect (1,200,000 m³)

CAVERN_TYPE = "interlayer"

# ── INTERLAYER MATERIAL SELECTION ──────────────────────────────────────────────
# Choose the material for each interlayer (used for all heterogeneous cavern types)
# Available materials:
#   "anhydrite" - E = 61.5 GPa, nu = 0.32  (stiffer, from Table 1)
#   "mudstone"  - E = 19.33 GPa, nu = 0.223 (softer, from Table 1)
#
# All interlayers are PURELY ELASTIC (no creep behavior).
# Note: INTERLAYER_3_MATERIAL is only used for "bulbous_ledges" caverns (3 interlayers)

INTERLAYER_1_MATERIAL = "anhydrite"  # Lower interlayer
INTERLAYER_2_MATERIAL = "mudstone"   # Middle/Upper interlayer
INTERLAYER_3_MATERIAL = "anhydrite"  # Top interlayer (bulbous_ledges only)

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
#   "sinus"     - Sinusoidal pressure variation around p_leach_end + amplitude
#   "linear"    - Piecewise linear pressure profile (ramp up, hold, ramp down)
#   "irregular" - Irregular/spline-smoothed pressure profile (realistic fluctuations)
#   "csv"       - Load pressure profile from CSV file (e.g., real operational data)
#
# All scenarios start from p_leach_end (= LEACHING_END_FRACTION * lithostatic pressure)

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
#   Note: Desai is only applied to SALT regions, not to interlayers
USE_DESAI = True

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                        END OF USER CONFIGURATION                              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


# ══════════════════════════════════════════════════════════════════════════════
# AUTOMATIC CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

VALID_CAVERN_TYPES = [
    "nointerlayer", "interlayer",
    "bulbous_ledges_600", "bulbous_ledges_1200",
    "asymmetric_shelf_600", "asymmetric_shelf_1200",
    "vertical_intrusion_600", "vertical_intrusion_1200",
]
VALID_SCENARIOS = ["sinus", "linear", "irregular", "csv"]
VALID_MODES_STANDARD = ["stretch", "repeat"]
VALID_MODES_CSV = ["stretch", "repeat", "direct"]
VALID_LEACHING_MODES = ["linear", "stepped"]
VALID_INTERLAYER_MATERIALS = ["anhydrite", "mudstone"]

# Material properties lookup table
INTERLAYER_PROPERTIES = {
    "anhydrite": {"E_GPa": 61.5, "nu": 0.32},
    "mudstone": {"E_GPa": 19.33, "nu": 0.223},
}

# Cavern parameters per type
# Each entry: {"z_max": top of cavern, "z_center": center z, "height": cavern height, "n_interlayers": number of interlayers}
CAVERN_PARAMS = {
    # Standard interlayer cavern (600k)
    "nointerlayer": {"z_max": 345.0, "z_center": 245.0, "height": 200.0, "n_interlayers": 0},
    "interlayer": {"z_max": 345.0, "z_center": 245.0, "height": 200.0, "n_interlayers": 2},
    # Bulbous ledges (3 interlayers)
    "bulbous_ledges_600": {"z_max": 340.0, "z_center": 240.0, "height": 200.0, "n_interlayers": 3},
    "bulbous_ledges_1200": {"z_max": 350.0, "z_center": 247.5, "height": 215.0, "n_interlayers": 3},
    # Asymmetric shelf (2 interlayers)
    "asymmetric_shelf_600": {"z_max": 345.0, "z_center": 240.0, "height": 210.0, "n_interlayers": 2},
    "asymmetric_shelf_1200": {"z_max": 355.0, "z_center": 247.5, "height": 225.0, "n_interlayers": 2},
    # Vertical intrusion effect (2 interlayers)
    "vertical_intrusion_600": {"z_max": 340.0, "z_center": 237.5, "height": 205.0, "n_interlayers": 2},
    "vertical_intrusion_1200": {"z_max": 355.0, "z_center": 247.5, "height": 225.0, "n_interlayers": 2},
}

# Domain dimensions
Z_SURFACE = 660.0  # meters (top of domain)

# Reference pressure (similar to 600k caverns)
P_REF_MPA = 17.5

# Grid folder mapping
GRID_FOLDERS = {
    "nointerlayer": "cavern_nointerlayer",
    "interlayer": "cavern_interlayer_600_3D",
    "bulbous_ledges_600": "cavern_bulbous_ledges_600_3D",
    "bulbous_ledges_1200": "cavern_bulbous_ledges_1200_3D",
    "asymmetric_shelf_600": "cavern_asymmetric_shelf_600_3D",
    "asymmetric_shelf_1200": "cavern_asymmetric_shelf_1200_3D",
    "vertical_intrusion_600": "cavern_vertical_intrusion_600_3D",
    "vertical_intrusion_1200": "cavern_vertical_intrusion_1200_3D",
}

# Mesh filename (standardized)
MESH_NAME = "geom"


def validate_configuration():
    """Validate user configuration and return derived values."""
    errors = []

    if CAVERN_TYPE not in VALID_CAVERN_TYPES:
        errors.append(f"CAVERN_TYPE '{CAVERN_TYPE}' invalid. Choose from: {VALID_CAVERN_TYPES}")
    if PRESSURE_SCENARIO not in VALID_SCENARIOS:
        errors.append(f"PRESSURE_SCENARIO '{PRESSURE_SCENARIO}' invalid. Choose from: {VALID_SCENARIOS}")
    if LEACHING_MODE not in VALID_LEACHING_MODES:
        errors.append(f"LEACHING_MODE '{LEACHING_MODE}' invalid. Choose from: {VALID_LEACHING_MODES}")

    # Get number of interlayers for this cavern type
    n_interlayers = CAVERN_PARAMS.get(CAVERN_TYPE, {}).get("n_interlayers", 0)

    # Validate interlayer materials (only relevant for heterogeneous cavern types)
    if n_interlayers >= 1:
        if INTERLAYER_1_MATERIAL not in VALID_INTERLAYER_MATERIALS:
            errors.append(f"INTERLAYER_1_MATERIAL '{INTERLAYER_1_MATERIAL}' invalid. Choose from: {VALID_INTERLAYER_MATERIALS}")
    if n_interlayers >= 2:
        if INTERLAYER_2_MATERIAL not in VALID_INTERLAYER_MATERIALS:
            errors.append(f"INTERLAYER_2_MATERIAL '{INTERLAYER_2_MATERIAL}' invalid. Choose from: {VALID_INTERLAYER_MATERIALS}")
    if n_interlayers >= 3:
        if INTERLAYER_3_MATERIAL not in VALID_INTERLAYER_MATERIALS:
            errors.append(f"INTERLAYER_3_MATERIAL '{INTERLAYER_3_MATERIAL}' invalid. Choose from: {VALID_INTERLAYER_MATERIALS}")

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

    grid_folder = GRID_FOLDERS[CAVERN_TYPE]
    cavern_params = CAVERN_PARAMS[CAVERN_TYPE]

    # Get interlayer properties
    il1_props = INTERLAYER_PROPERTIES.get(INTERLAYER_1_MATERIAL, INTERLAYER_PROPERTIES["anhydrite"])
    il2_props = INTERLAYER_PROPERTIES.get(INTERLAYER_2_MATERIAL, INTERLAYER_PROPERTIES["mudstone"])
    il3_props = INTERLAYER_PROPERTIES.get(INTERLAYER_3_MATERIAL, INTERLAYER_PROPERTIES["anhydrite"])

    return {
        "cavern_type": CAVERN_TYPE,
        "grid_folder": grid_folder,
        "z_max": cavern_params["z_max"],
        "z_center": cavern_params["z_center"],
        "cavern_height": cavern_params["height"],
        "n_interlayers": n_interlayers,
        "p_ref_mpa": P_REF_MPA,
        "interlayer_1_material": INTERLAYER_1_MATERIAL,
        "interlayer_2_material": INTERLAYER_2_MATERIAL,
        "interlayer_3_material": INTERLAYER_3_MATERIAL,
        "interlayer_1_E_GPa": il1_props["E_GPa"],
        "interlayer_1_nu": il1_props["nu"],
        "interlayer_2_E_GPa": il2_props["E_GPa"],
        "interlayer_2_nu": il2_props["nu"],
        "interlayer_3_E_GPa": il3_props["E_GPa"],
        "interlayer_3_nu": il3_props["nu"],
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


# ══════════════════════════════════════════════════════════════════════════════
# CSV PRESSURE PROFILE FUNCTIONS
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
    """Read CSV and return pressure array in MPa."""
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


def build_csv_pressure_schedule(tc, csv_file, *, days, mode, total_cycles=1,
                                rescale=False, rescale_min=None, rescale_max=None,
                                resample_at_dt=True):
    """Build pressure schedule from CSV file."""
    pressures_mpa = read_pressure_csv(csv_file)
    csv_hours = int(pressures_mpa.size)

    if MPI.COMM_WORLD.rank == 0:
        print(f"[CSV] Loaded '{os.path.basename(csv_file)}' with {csv_hours} hourly values")
        print(f"[CSV] Original range: [{pressures_mpa.min():.2f}, {pressures_mpa.max():.2f}] MPa")

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
    """Replace first part of schedule with a linear ramp."""
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


def build_linear_schedule_multi(tc, times_h, pressures_MPa, *, days, mode,
                                resample_at_dt=True, total_cycles=None):
    if mode not in ("repeat", "stretch"):
        raise ValueError("mode must be 'repeat' or 'stretch'")

    times_h = list(map(float, times_h))
    pressures_MPa = list(map(float, pressures_MPa))
    if len(times_h) != len(pressures_MPa) or len(times_h) < 2:
        raise ValueError("Provide at least two waypoints with matching lengths.")

    total_h = float(days) * DAY_H

    if mode == "stretch":
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
    else:  # repeat
        t_h = []
        p_h = []
        for d in range(int(days)):
            off = d * DAY_H
            for i, t in enumerate(times_h):
                if d > 0 and i == 0:
                    continue
                t_h.append(off + t)
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
# MATERIAL SETUP FOR INTERLAYER HETEROGENEITY
# ══════════════════════════════════════════════════════════════════════════════

def setup_material_homogeneous(mat, n_elems):
    """Setup homogeneous salt material (no interlayers)."""
    sec_per_year = 365.25 * 24 * 3600

    # Elastic (Spring)
    E0 = 20.425 * ut.GPa * to.ones(n_elems)
    nu0 = 0.25 * to.ones(n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")

    # Kelvin-Voigt (Viscoelastic)
    eta = 105e11 * to.ones(n_elems)
    E1 = 10 * ut.GPa * to.ones(n_elems)
    nu1 = 0.25 * to.ones(n_elems)
    kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")

    # Dislocation creep
    ndc = 4.6
    A_dc = (40.0 * (1e-6)**ndc / sec_per_year) * to.ones(n_elems)
    Q_dc = (6495.0 * 8.32) * to.ones(n_elems)
    n_dc = ndc * to.ones(n_elems)
    creep_0 = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation")

    # Pressure-solution creep
    A_ps = (14176.0 * 1e-9 / 1e6 / sec_per_year) * to.ones(n_elems)
    d_ps = 5.25e-3 * to.ones(n_elems)
    Q_ps = (3252.0 * 8.32) * to.ones(n_elems)
    creep_pressure = sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "creep_pressure")

    mat.add_to_elastic(spring_0)
    mat.add_to_non_elastic(kelvin)
    mat.add_to_non_elastic(creep_0)
    mat.add_to_non_elastic(creep_pressure)

    return mat


def setup_material_heterogeneous(mat, n_elems, grid, config):
    """
    Setup heterogeneous material with interlayers.

    Material properties:
    - Salt (all salt regions): Full creep model
    - Interlayer_1, _2, _3: Material specified by config (anhydrite or mudstone), PURELY ELASTIC
    """
    sec_per_year = 365.25 * 24 * 3600

    # Get interlayer properties from config
    il1_E = config["interlayer_1_E_GPa"]
    il1_nu = config["interlayer_1_nu"]
    il1_mat = config["interlayer_1_material"]
    il2_E = config["interlayer_2_E_GPa"]
    il2_nu = config["interlayer_2_nu"]
    il2_mat = config["interlayer_2_material"]
    il3_E = config["interlayer_3_E_GPa"]
    il3_nu = config["interlayer_3_nu"]
    il3_mat = config["interlayer_3_material"]
    n_interlayers = config["n_interlayers"]

    # Get region indices from grid
    region_indices = grid.region_indices

    # Identify salt and interlayer regions
    salt_cells = []
    interlayer_1_cells = []
    interlayer_2_cells = []
    interlayer_3_cells = []

    for region_name, indices in region_indices.items():
        region_lower = region_name.lower()
        if "interlayer_1" in region_lower or "interlayer1" in region_lower:
            interlayer_1_cells.extend(indices)
        elif "interlayer_2" in region_lower or "interlayer2" in region_lower:
            interlayer_2_cells.extend(indices)
        elif "interlayer_3" in region_lower or "interlayer3" in region_lower:
            interlayer_3_cells.extend(indices)
        elif "salt" in region_lower:
            salt_cells.extend(indices)

    salt_cells = np.array(salt_cells, dtype=int) if salt_cells else np.array([], dtype=int)
    interlayer_1_cells = np.array(interlayer_1_cells, dtype=int) if interlayer_1_cells else np.array([], dtype=int)
    interlayer_2_cells = np.array(interlayer_2_cells, dtype=int) if interlayer_2_cells else np.array([], dtype=int)
    interlayer_3_cells = np.array(interlayer_3_cells, dtype=int) if interlayer_3_cells else np.array([], dtype=int)

    # Combine all interlayer cells for operations that apply to all
    all_interlayer_cells = np.concatenate([interlayer_1_cells, interlayer_2_cells, interlayer_3_cells])

    total_assigned = len(salt_cells) + len(all_interlayer_cells)

    if MPI.COMM_WORLD.rank == 0:
        print(f"[MATERIAL] Region assignment:")
        print(f"  Salt cells: {len(salt_cells)}")
        print(f"  Interlayer_1 ({il1_mat}) cells: {len(interlayer_1_cells)}")
        print(f"  Interlayer_2 ({il2_mat}) cells: {len(interlayer_2_cells)}")
        if n_interlayers >= 3:
            print(f"  Interlayer_3 ({il3_mat}) cells: {len(interlayer_3_cells)}")
        print(f"  Total: {total_assigned} / {n_elems}")

    # ── ELASTIC PROPERTIES ────────────────────────────────────────────────────
    # Salt: E = 20.425 GPa, nu = 0.25
    # Interlayers: As specified in config

    E0 = to.zeros(n_elems, dtype=to.float64)
    nu0 = to.zeros(n_elems, dtype=to.float64)

    # Default to salt properties
    E0[:] = 20.425 * ut.GPa
    nu0[:] = 0.25

    # Assign interlayer properties from config
    if len(interlayer_1_cells) > 0:
        E0[interlayer_1_cells] = il1_E * ut.GPa
        nu0[interlayer_1_cells] = il1_nu
    if len(interlayer_2_cells) > 0:
        E0[interlayer_2_cells] = il2_E * ut.GPa
        nu0[interlayer_2_cells] = il2_nu
    if len(interlayer_3_cells) > 0:
        E0[interlayer_3_cells] = il3_E * ut.GPa
        nu0[interlayer_3_cells] = il3_nu

    spring_0 = sf.Spring(E0, nu0, "spring")
    mat.add_to_elastic(spring_0)

    # ── VISCOELASTIC (Kelvin-Voigt) ───────────────────────────────────────────
    # Only for salt regions; interlayers have eta -> infinity (or very high)
    # To disable viscoelastic behavior in interlayers, set eta to very high value

    eta = to.zeros(n_elems, dtype=to.float64)
    E1 = to.zeros(n_elems, dtype=to.float64)
    nu1 = to.zeros(n_elems, dtype=to.float64)

    # Salt: viscoelastic
    eta[:] = 105e11  # Salt viscosity
    E1[:] = 10 * ut.GPa
    nu1[:] = 0.25

    # Interlayers: very high viscosity (effectively rigid viscoelastic)
    # This effectively disables viscoelastic strain in interlayers
    if len(interlayer_1_cells) > 0:
        eta[interlayer_1_cells] = 1e30  # Very high -> no viscoelastic strain
        E1[interlayer_1_cells] = il1_E * ut.GPa
        nu1[interlayer_1_cells] = il1_nu
    if len(interlayer_2_cells) > 0:
        eta[interlayer_2_cells] = 1e30
        E1[interlayer_2_cells] = il2_E * ut.GPa
        nu1[interlayer_2_cells] = il2_nu
    if len(interlayer_3_cells) > 0:
        eta[interlayer_3_cells] = 1e30
        E1[interlayer_3_cells] = il3_E * ut.GPa
        nu1[interlayer_3_cells] = il3_nu

    kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")
    mat.add_to_non_elastic(kelvin)

    # ── DISLOCATION CREEP ─────────────────────────────────────────────────────
    # Only for salt; interlayers have A = 0 (no creep)

    ndc = 4.6
    A_dc = to.zeros(n_elems, dtype=to.float64)
    Q_dc = to.zeros(n_elems, dtype=to.float64)
    n_dc = to.zeros(n_elems, dtype=to.float64)

    # Salt: dislocation creep
    A_dc[:] = 40.0 * (1e-6)**ndc / sec_per_year
    Q_dc[:] = 6495.0 * 8.32
    n_dc[:] = ndc

    # Interlayers: A = 0 -> no dislocation creep
    if len(all_interlayer_cells) > 0:
        A_dc[all_interlayer_cells] = 0.0

    creep_0 = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation")
    mat.add_to_non_elastic(creep_0)

    # ── PRESSURE-SOLUTION CREEP ───────────────────────────────────────────────
    # Only for salt; interlayers have A = 0 (no creep)

    A_ps = to.zeros(n_elems, dtype=to.float64)
    d_ps = to.zeros(n_elems, dtype=to.float64)
    Q_ps = to.zeros(n_elems, dtype=to.float64)

    # Salt: pressure-solution creep
    A_ps[:] = 14176.0 * 1e-9 / 1e6 / sec_per_year
    d_ps[:] = 5.25e-3
    Q_ps[:] = 3252.0 * 8.32

    # Interlayers: A = 0 -> no pressure-solution creep
    if len(all_interlayer_cells) > 0:
        A_ps[all_interlayer_cells] = 0.0

    creep_pressure = sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "creep_pressure")
    mat.add_to_non_elastic(creep_pressure)

    if MPI.COMM_WORLD.rank == 0:
        print(f"[MATERIAL] Heterogeneous material setup complete:")
        print(f"  Salt: E=20.425 GPa, nu=0.25, full creep model")
        print(f"  Interlayer_1 ({il1_mat}): E={il1_E} GPa, nu={il1_nu}, purely elastic")
        print(f"  Interlayer_2 ({il2_mat}): E={il2_E} GPa, nu={il2_nu}, purely elastic")
        if n_interlayers >= 3:
            print(f"  Interlayer_3 ({il3_mat}): E={il3_E} GPa, nu={il3_nu}, purely elastic")

    return mat, salt_cells, interlayer_1_cells, interlayer_2_cells, interlayer_3_cells


def add_desai_heterogeneous(mat, mom_eq, salt_cells, interlayer_1_cells, interlayer_2_cells, interlayer_3_cells=None):
    """
    Add Desai viscoplastic model for operation phase.
    Only applied to salt regions; interlayers are purely elastic.
    """
    if interlayer_3_cells is None:
        interlayer_3_cells = np.array([], dtype=int)

    n_elems = mom_eq.n_elems

    # Combine all interlayer cells
    all_interlayer_cells = np.concatenate([interlayer_1_cells, interlayer_2_cells, interlayer_3_cells])

    # Initialize all Desai parameters
    mu_1 = to.zeros(n_elems, dtype=to.float64)
    N_1 = to.zeros(n_elems, dtype=to.float64)
    n = to.zeros(n_elems, dtype=to.float64)
    a_1 = to.zeros(n_elems, dtype=to.float64)
    eta_vp = to.zeros(n_elems, dtype=to.float64)
    beta_1 = to.zeros(n_elems, dtype=to.float64)
    beta = to.zeros(n_elems, dtype=to.float64)
    m = to.zeros(n_elems, dtype=to.float64)
    gamma = to.zeros(n_elems, dtype=to.float64)
    alpha_0 = to.zeros(n_elems, dtype=to.float64)
    sigma_t = to.zeros(n_elems, dtype=to.float64)

    # Salt regions: standard Desai parameters
    mu_1[:] = 5.3665857009859815e-11
    N_1[:] = 3.1
    n[:] = 3.0
    a_1[:] = 1.965018496922832e-05
    eta_vp[:] = 0.8275682807874163
    beta_1[:] = 0.0048
    beta[:] = 0.995
    m[:] = -0.5
    gamma[:] = 0.095
    alpha_0[:] = 0.0022
    sigma_t[:] = 5.0

    # Interlayers: Set mu_1 = 0 to disable viscoplasticity
    # (Setting mu_1 = 0 makes the viscoplastic strain rate zero)
    if len(all_interlayer_cells) > 0:
        mu_1[all_interlayer_cells] = 0.0
        eta_vp[all_interlayer_cells] = 1e30  # Very high viscosity

    desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta_vp, n, beta_1, beta, m, gamma, sigma_t, alpha_0, "desai")

    # Compute initial hardening from current stress state
    stress_to = ut.numpy2torch(mom_eq.sig.x.array.reshape((n_elems, 3, 3)))
    desai.compute_initial_hardening(stress_to, Fvp_0=0.0)

    mat.add_to_non_elastic(desai)

    if MPI.COMM_WORLD.rank == 0:
        print(f"[MATERIAL] Desai viscoplasticity added (salt only, interlayers excluded)")

    return mat


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    config = validate_configuration()

    # Material constants
    salt_density = 2200  # kg/m³
    g = -9.81  # m/s²

    # Compute lithostatic pressure
    p_lithostatic = compute_lithostatic_pressure(
        config["z_center"],
        config["p_ref_mpa"],
        salt_density,
        g
    )
    p_lithostatic_mpa = p_lithostatic / ut.MPa

    # Compute leaching end pressure
    p_leach_end_mpa = LEACHING_END_FRACTION * p_lithostatic_mpa
    p_leach_end = p_leach_end_mpa * ut.MPa

    # Derive operational pressures
    if PRESSURE_SCENARIO == "sinus":
        p_op_min_mpa = p_leach_end_mpa
        p_op_mean_mpa = p_leach_end_mpa + P_AMPLITUDE_MPA
        p_op_max_mpa = p_leach_end_mpa + 2 * P_AMPLITUDE_MPA
    elif PRESSURE_SCENARIO == "linear":
        p_op_min_mpa = p_leach_end_mpa
        p_op_max_mpa = p_leach_end_mpa + PRESSURE_SWING_MPA
    else:
        p_op_min_mpa = p_leach_end_mpa

    # Print configuration
    if MPI.COMM_WORLD.rank == 0:
        print("=" * 70)
        print("INTERLAYER CAVERN SIMULATION")
        print("=" * 70)
        print(f"  Cavern type:      {config['cavern_type']}")
        print(f"  Grid folder:      {config['grid_folder']}")
        print(f"  Cavern z_max:     {config['z_max']:.2f} m")
        print(f"  Cavern z_center:  {config['z_center']:.2f} m")
        print("-" * 70)
        print("  LEACHING PHASE:")
        print(f"    Mode:           {LEACHING_MODE}")
        print(f"    Duration:       {LEACHING_DAYS} days")
        print(f"    P_start:        {p_lithostatic_mpa:.2f} MPa (lithostatic)")
        print(f"    P_end:          {p_leach_end_mpa:.2f} MPa ({LEACHING_END_FRACTION*100:.0f}%)")
        print("-" * 70)
        print("  OPERATION PHASE:")
        print(f"    Scenario:       {PRESSURE_SCENARIO}")
        print(f"    Days:           {OPERATION_DAYS}")
        print(f"    Desai:          {'enabled' if USE_DESAI else 'disabled'}")
        if config["n_interlayers"] > 0:
            print("-" * 70)
            print("  MATERIAL HETEROGENEITY:")
            print("    Salt:         Full creep model (elastic + viscoelastic + creep + Desai)")
            print(f"    Interlayer_1: {config['interlayer_1_material'].capitalize()} - E={config['interlayer_1_E_GPa']} GPa, nu={config['interlayer_1_nu']}, PURELY ELASTIC")
            print(f"    Interlayer_2: {config['interlayer_2_material'].capitalize()} - E={config['interlayer_2_E_GPa']} GPa, nu={config['interlayer_2_nu']}, PURELY ELASTIC")
            if config["n_interlayers"] >= 3:
                print(f"    Interlayer_3: {config['interlayer_3_material'].capitalize()} - E={config['interlayer_3_E_GPa']} GPa, nu={config['interlayer_3_nu']}, PURELY ELASTIC")
        print("=" * 70)

    # Load grid
    grid_path = os.path.join("..", "..", "..", "..", "grids", config["grid_folder"])
    grid = sf.GridHandlerGMSH(MESH_NAME, grid_path)

    z_max = config["z_max"]
    p_ref = config["p_ref_mpa"] * ut.MPa

    # Output folder
    output_folder = os.path.join(
        "output",
        f"case_{config['cavern_type']}_{PRESSURE_SCENARIO}({N_CYCLES})_{OPERATION_DAYS}days"
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

    # Set material density (same for all regions for simplicity)
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    # Setup material based on cavern type
    n_interlayers = config["n_interlayers"]
    if n_interlayers == 0:
        mat = setup_material_homogeneous(mat, mom_eq.n_elems)
        salt_cells = np.arange(mom_eq.n_elems)
        interlayer_1_cells = np.array([], dtype=int)
        interlayer_2_cells = np.array([], dtype=int)
        interlayer_3_cells = np.array([], dtype=int)
    else:  # heterogeneous with interlayers
        mat, salt_cells, interlayer_1_cells, interlayer_2_cells, interlayer_3_cells = setup_material_heterogeneous(
            mat, mom_eq.n_elems, grid, config
        )

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

    t_leaching, p_leaching = build_leaching_pressure_schedule(
        tc_leaching,
        p_start_pa=p_lithostatic,
        p_end_pa=p_leach_end,
        mode=LEACHING_MODE,
        n_steps=STEPPED_N_STEPS
    )

    gas_density = 0.089  # Hydrogen

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

    leaching_save_interval = max(1, int(72 / LEACHING_DT_HOURS))
    output_mom_leaching = SparseSaveFields(mom_eq, interval=leaching_save_interval)
    output_mom_leaching.set_output_folder(output_folder_leaching)
    output_mom_leaching.add_output_field("u", "Displacement (m)")
    output_mom_leaching.add_output_field("eps_tot", "Total strain (-)")
    output_mom_leaching.add_output_field("sig", "Stress (Pa)")
    output_mom_leaching.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom_leaching.add_output_field("q_elems", "Von Mises stress (Pa)")
    outputs_leaching = [output_mom_leaching]

    # Save leaching schedule
    os.makedirs(output_folder, exist_ok=True)
    leaching_data = {
        "phase": "leaching",
        "cavern_type": config["cavern_type"],
        "leaching_mode": LEACHING_MODE,
        "leaching_days": LEACHING_DAYS,
        "p_lithostatic_mpa": p_lithostatic_mpa,
        "p_leach_end_mpa": p_leach_end_mpa,
        "t_hours": [float(t / ut.hour) for t in t_leaching],
        "p_MPa": [float(p / ut.MPa) for p in p_leaching],
    }
    with open(os.path.join(output_folder, "leaching_schedule.json"), 'w') as f:
        json.dump(leaching_data, f, indent=2)

    if MPI.COMM_WORLD.rank == 0:
        print(f"[LEACHING] Running for {LEACHING_DAYS} days...")

    sim_leaching = sf.Simulator_M(mom_eq, tc_leaching, outputs_leaching, True)
    sim_leaching.run()

    if MPI.COMM_WORLD.rank == 0:
        print("[LEACHING] Complete.")

    # ===== OPERATION PHASE =====
    if USE_DESAI:
        if n_interlayers > 0:
            mat = add_desai_heterogeneous(mat, mom_eq, salt_cells, interlayer_1_cells, interlayer_2_cells, interlayer_3_cells)
        else:
            # Homogeneous Desai
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
            rescale=False,
            resample_at_dt=True
        )

        if CSV_SHIFT_TO_LEACH_END:
            p_array = np.array(p_pressure)
            csv_min_pa = p_array.min()
            shift_pa = p_leach_end - csv_min_pa
            p_pressure = [p + shift_pa for p in p_pressure]
            if MPI.COMM_WORLD.rank == 0:
                new_min = min(p_pressure) / ut.MPa
                new_max = max(p_pressure) / ut.MPa
                print(f"[CSV] Shifted profile: [{new_min:.2f}, {new_max:.2f}] MPa")

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

    # Save pressure schedule
    pressure_data = {
        "phase": "operation",
        "cavern_type": config["cavern_type"],
        "scenario": PRESSURE_SCENARIO,
        "p_leach_end_mpa": p_leach_end_mpa,
        "mode": SCHEDULE_MODE,
        "n_cycles": N_CYCLES,
        "operation_days": OPERATION_DAYS,
        "use_desai": USE_DESAI,
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

    if MPI.COMM_WORLD.rank == 0:
        print(f"[OPERATION] Running for {OPERATION_DAYS} days...")

    sim_op = sf.Simulator_M(mom_eq, tc_operation, outputs_op, False)
    sim_op.run()

    if MPI.COMM_WORLD.rank == 0:
        print("[OPERATION] Complete.")
        print("=" * 70)
        print("SIMULATION FINISHED")
        print(f"  Output folder: {output_folder}")
        print("=" * 70)


if __name__ == '__main__':
    main()
