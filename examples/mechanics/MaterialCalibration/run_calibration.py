"""
Material parameter calibration: single-element creep simulation.

Reads the multi-stage triaxial creep test stress history from the ZWD CSV files,
then integrates the creep strain using either:
  (a) SafeInCave model  (Spring + Kelvin + DislocationCreep + Desai)
  (b) Munson-Dawson model (Spring + MunsonDawsonCreep)

No FEM mesh is needed — this is a 1-element, uniaxial-equivalent integration
under prescribed differential stress σ_diff(t) and confining pressure σ_3(t).

Output: JSON files with time vs. axial strain for each model configuration,
used by plot_calibration.py for comparison against lab data.
"""

import os
import json
import numpy as np

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          USER CONFIGURATION                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── SELECT WHICH LAB TESTS TO SIMULATE ────────────────────────────────────────
# Use test names like "TCC1", "TCC7", etc. Set to None to run all 12.
TESTS = ["TCC1", "TCC2", "TCC6", "TCC7", "TCC11", "TCC12"]

# ── SELECT WHICH MODEL(S) TO RUN ─────────────────────────────────────────────
RUN_SAFEINCAVE = True
RUN_MUNSONDAWSON = True

# ── TEMPERATURE (all tests ran at 100 °C) ────────────────────────────────────
T_CELSIUS = 100.0
T_KELVIN = T_CELSIUS + 273.15

# ── TIME INTEGRATION ─────────────────────────────────────────────────────────
DT_HOURS = 0.5       # integration time step (hours)

# ── GAS CONSTANT ─────────────────────────────────────────────────────────────
R_GAS = 8.314        # J/(mol·K) — standard value

# ══════════════════════════════════════════════════════════════════════════════
# CALIBRATION PARAMETERS — EDIT THESE TO CALIBRATE
# ══════════════════════════════════════════════════════════════════════════════

# ── SHARED: Dislocation creep (used by both models) ──────────────────────────
#   ε̇_ss = A · σ^n · exp(-Q/(R·T))
#   A is specified in [MPa^{-n} / yr] and converted to [Pa^{-n} / s]
Q_OVER_R = 6252.0                        # [K]  — CALIBRATED (range 6252–6495)
Q_DISLOC = Q_OVER_R * R_GAS              # [J/mol]

# -- Munson-Dawson dislocation creep (calibrated independently via search6) --
_sec_per_year = 365.25 * 24 * 3600
_A_MD_MPA_YR  = 188.1      # [MPa^{-n}/yr]  — search_weighted: TCC1, 7-param (beta_w free)
N_MD          = 4.4483     # [-]            — search_weighted
A_MD          = _A_MD_MPA_YR * (1e-6) ** N_MD / _sec_per_year   # [Pa^{-n}/s]

# -- SafeInCave dislocation creep (independent calibration via search_sic) ---
_A_SIC_MPA_YR = 712.68     # [MPa^{-n}/yr]  — search_weighted: TCC1, balanced weights
N_SIC         = 4.0557     # [-]            — search_weighted
A_SIC         = _A_SIC_MPA_YR * (1e-6) ** N_SIC / _sec_per_year # [Pa^{-n}/s]

# Backward-compat aliases (used by integrate_munsondawson)
N_DISLOC = N_MD
A_DISLOC = A_MD
_A_MPA_YR = _A_MD_MPA_YR

# ── SafeInCave: Kelvin element (transient creep) ─────────────────────────────
#   Kelvin dashpot viscosity η [Pa·s], spring E1 [Pa], Poisson ν1 [-]
ETA_KELVIN = 3.3504e12     # [Pa·s]   — CALIBRATED (search_weighted: TCC1, τ=1.1h)
E1_KELVIN  = 8.2841e8      # [Pa]     — CALIBRATED (search_weighted, eq@21MPa=2.53%)
NU1_KELVIN = 0.25          # [-]      — CALIBRATED

# ── SafeInCave: Desai viscoplastic ───────────────────────────────────────────
#   Only the parameters that need calibration are listed here.
#   The rest are fixed (same as in ScenarioTest.py).
MU1_DESAI     = 1e-15                     # [1/s]  — CALIBRATED (search3)
N1_DESAI      = 2.0                       # [-]    — CALIBRATED
A1_DESAI      = 1e-03                     # [-]    — CALIBRATED
ETA_DESAI     = 1.2                       # [-]    — CALIBRATED
ALPHA0_DESAI  = 0.005                     # [-]    — CALIBRATED

# Fixed Desai parameters
N_DESAI_POWER = 3.0
BETA1_DESAI   = 0.0048
BETA_DESAI    = 0.995
M_DESAI       = -0.5
GAMMA_DESAI   = 0.095
SIGMA_T_DESAI = 5.0       # [Pa]

# ── Munson-Dawson: Transient parameters ──────────────────────────────────────
K0_MD      = 0.0041      # [-]    — CALIBRATED (search_weighted: TCC1, 7-param)
C_MD       = 0.00902     # [1/K]  — fixed
M_MD       = 0.2000      # [-]    — CALIBRATED (search_weighted)
ALPHA_W_MD = -2.8951     # [-]    — CALIBRATED (search_weighted)

# Munson-Dawson transient parameters
BETA_W_MD  = -4.7745     # [-]    — CALIBRATED (search_weighted, was fixed at -7.738)
DELTA_MD   = 1.1312      # [-]    — CALIBRATED (search_weighted)

# ── Elastic parameters (for shear modulus in Munson-Dawson) ──────────────────
E_ELASTIC  = 20.425e9    # [Pa]
NU_ELASTIC = 0.25        # [-]
MU_SHEAR   = E_ELASTIC / (2.0 * (1.0 + NU_ELASTIC))  # ~8.17 GPa

# ── PATHS ────────────────────────────────────────────────────────────────────
DATA_DIR = os.path.join(os.path.dirname(__file__),
                        "ZWD_Creeptests_rawdata", "Creep tests")
OUT_DIR  = os.path.join(os.path.dirname(__file__), "output")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                      END OF USER CONFIGURATION                             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


HOUR = 3600.0
DAY = 24.0 * HOUR


# ══════════════════════════════════════════════════════════════════════════════
# STRESS ARTEFACT CORRECTION
# ══════════════════════════════════════════════════════════════════════════════

def correct_stress_artefacts(time_h, sigma_diff,
                              high_thresh=15.0,
                              dip_thresh=5.0,
                              max_dip_h=24.0):
    """
    Fill brief dips in differential stress that are measurement artefacts.

    In some TCC tests (TCC7, TCC11, TCC12) the raw data shows brief drops
    from a high-stress plateau (>15 MPa) to near-zero (<5 MPa) for <24 hours,
    followed by a return to the same high plateau.  The real test protocol
    only decreases stress; these dips are equipment artefacts and should be
    replaced by the preceding plateau value before model integration.

    Parameters
    ----------
    time_h : array  elapsed time [hours]
    sigma_diff : array  differential stress [MPa]
    high_thresh : float  minimum stress considered "high" [MPa]
    dip_thresh : float  maximum stress during a dip [MPa]
    max_dip_h : float  maximum duration of a valid artefact dip [hours]

    Returns
    -------
    array  corrected differential stress [MPa]
    """
    sig = sigma_diff.copy()
    n = len(time_h)
    i = 0
    while i < n - 1:
        if sig[i] >= high_thresh and sig[i + 1] < dip_thresh:
            dip_start = i + 1
            plateau_val = sig[i]
            j = dip_start
            while j < n and sig[j] < dip_thresh:
                j += 1
            if j < n and sig[j] >= high_thresh:
                dip_duration = time_h[j] - time_h[dip_start]
                if dip_duration <= max_dip_h:
                    sig[dip_start:j] = plateau_val
                    i = j
                    continue
        i += 1
    return sig


# ══════════════════════════════════════════════════════════════════════════════
# CSV READING
# ══════════════════════════════════════════════════════════════════════════════

def read_creep_csv(filepath):
    """
    Read a ZWD creep test CSV and return arrays of:
      time_h      — elapsed time [hours]
      strain_pct  — axial technical strain [%]
      sigma_diff  — differential stress σ1 - σ3 [MPa]
      sigma3      — confining pressure [MPa]
      temp_C      — cell temperature [°C]

    Skips header rows and empty trailing rows.
    """
    with open(filepath, "r", encoding="utf-8-sig") as f:
        lines = f.readlines()

    # Find the units row (starts with [TT.MM)
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith("[TT.MM"):
            data_start = i + 1
            break

    if data_start is None:
        raise ValueError(f"Could not find data start in {filepath}")

    time_h = []
    strain_pct = []
    sigma_diff = []
    sigma3 = []
    temp_C = []

    for line in lines[data_start:]:
        parts = line.strip().split(",")
        if len(parts) < 12 or parts[2].strip() == "":
            continue
        try:
            t = float(parts[2])    # time [h]
            e = float(parts[4])    # strain [%]
            sd = float(parts[9])   # sDiff-eff [MPa]
            s3 = float(parts[10])  # σ3 [MPa]
            tc = float(parts[11])  # cell temperature [°C]
        except (ValueError, IndexError):
            continue
        time_h.append(t)
        strain_pct.append(e)
        sigma_diff.append(sd)
        sigma3.append(s3)
        temp_C.append(tc)

    return (np.array(time_h), np.array(strain_pct),
            np.array(sigma_diff), np.array(sigma3), np.array(temp_C))


def build_stress_schedule(time_h_lab, sigma_diff_lab, dt_hours):
    """
    Build a fine time grid with piecewise-constant differential stress
    interpolated from the lab data.

    Returns:
        t_s     — time array [seconds], starting at the first deviatoric loading
        sigma   — differential stress [Pa] at each time step
    """
    # Find the first time point with significant differential stress (> 1 MPa)
    mask = sigma_diff_lab > 1.0
    if not np.any(mask):
        raise ValueError("No deviatoric loading found in data")

    idx_start = np.argmax(mask)
    t0_h = time_h_lab[idx_start]
    t_end_h = time_h_lab[-1]

    # Lab time/stress pairs (shifted so deviatoric loading starts at t=0)
    t_lab_s = (time_h_lab[idx_start:] - t0_h) * HOUR
    sigma_lab_Pa = sigma_diff_lab[idx_start:] * 1e6

    # Fine time grid
    n_steps = int(np.ceil((t_end_h - t0_h) / dt_hours))
    t_s = np.linspace(0.0, (t_end_h - t0_h) * HOUR, n_steps + 1)

    # Piecewise-constant interpolation (zero-order hold)
    sigma_Pa = np.interp(t_s, t_lab_s, sigma_lab_Pa)

    # Also return the lab reference arrays shifted to the same t=0
    t_lab_shifted_h = time_h_lab[idx_start:] - t0_h
    strain_lab_shifted = None  # caller handles this

    return t_s, sigma_Pa, idx_start


# ══════════════════════════════════════════════════════════════════════════════
# CREEP RATE FUNCTIONS (scalar, single-element)
# ══════════════════════════════════════════════════════════════════════════════

def disloc_rate(sigma_diff_Pa, A, n, Q, T):
    """
    Dislocation creep axial strain rate [1/s].
    For a uniaxial/triaxial test:  ε̇_axial = A · σ_diff^n · exp(-Q/(R·T))
    (deviatoric stress = σ_diff in triaxial test since q = σ1 - σ3)
    """
    return A * (abs(sigma_diff_Pa) ** n) * np.exp(-Q / (R_GAS * T))


def kelvin_rate(eps_kelvin, sigma_diff_Pa, eta, E1):
    """
    Kelvin element axial strain rate [1/s].
    ε̇_kelvin = (σ_dev - E1·ε_kelvin) / η

    In uniaxial: σ_dev = σ_diff (the differential stress drives the dashpot).
    """
    return (sigma_diff_Pa - E1 * eps_kelvin) / eta


def desai_rate(sigma_diff_Pa, sigma3_Pa, alpha, mu1, N1, a1, eta_vp,
               n_pow, beta1, beta, m, gamma, sigma_t):
    """
    Desai viscoplastic axial strain rate [1/s] (simplified scalar version).

    Evaluates the Desai yield function at the current stress state and returns
    the viscoplastic strain rate magnitude in the axial direction.

    NOTE: The Desai model works in MPa internally (matching MaterialProps.py).
    Stress inputs are in Pa and converted here.
    """
    MPA = 1e6

    # Compression-positive, in MPa (matching MaterialProps.py convention)
    s1 = (sigma3_Pa + sigma_diff_Pa) / MPA   # axial stress [MPa]
    s2 = sigma3_Pa / MPA                      # confining [MPa]
    s3_val = sigma3_Pa / MPA                  # confining [MPa]

    I1 = s1 + s2 + s3_val
    I2 = s1 * s2 + s2 * s3_val + s1 * s3_val
    I3 = s1 * s2 * s3_val
    J2 = (1.0 / 3.0) * I1**2 - I2
    J3 = (2.0 / 27.0) * I1**3 - (1.0 / 3.0) * I1 * I2 + I3

    if J2 < 1e-30:
        return 0.0

    # Lode parameter Sr
    sqJ2 = np.sqrt(J2)
    arg_sr = -(J3 * np.sqrt(27.0)) / (2.0 * sqJ2**3)
    arg_sr = np.clip(arg_sr, -1.0, 1.0)
    Sr = arg_sr  # simplified; full model uses arcsin

    # Shifted first invariant
    I1_star = I1 + sigma_t  # sigma_t is in MPa

    # Yield function (all in MPa, matching MaterialProps.py line 1214-1216)
    F1 = alpha * I1_star**n_pow - gamma * I1_star**2
    F2_base = np.exp(beta1 * I1_star) - beta * Sr
    # Handle negative base with fractional exponent
    if F2_base <= 0.0:
        F_vp = J2  # F2 term vanishes
    else:
        F_vp = J2 + F1 * F2_base**m

    if F_vp <= 0.0:
        return 0.0

    # Viscoplastic multiplier
    lam = mu1 * (F_vp) ** N1

    # Flow direction: ∂F/∂σ has a J2-derivative part ∂J2/∂σ = s_dev
    # In triaxial: s_dev_axial = s1 - I1/3
    p = I1 / 3.0
    dev1 = s1 - p

    # The strain rate in 1/s; dev1 is in MPa, but mu1 already absorbs the units
    # so the product lam * |dev1| gives 1/s
    eps_rate = lam * abs(dev1)
    return eps_rate


def desai_update_alpha(alpha_old, eps_vp_rate, dt, a1, alpha0, eta_vp):
    """Update the Desai hardening variable α."""
    zeta_increment = abs(eps_vp_rate) * dt
    # α evolves inversely with accumulated viscoplastic strain
    # Simplified: α decreases (hardens) as ζ grows
    # Full formula: α = a1 / [(a1/α0)^(1/η) + ζ]^η
    # We track ζ as a running sum
    return zeta_increment


def munsondawson_rate(sigma_diff_Pa, zeta, A, n, Q, T,
                      K0, c, m_md, alpha_w, beta_w, delta, mu):
    """
    Munson-Dawson creep axial strain rate [1/s] and updated zeta.

    Returns (eps_rate, F_value) where:
      eps_rate = F · ε̇_ss  (total creep rate including transient)
      F_value = transient multiplier
    """
    # Steady-state rate
    eps_ss = A * (abs(sigma_diff_Pa) ** n) * np.exp(-Q / (R_GAS * T))

    if eps_ss < 1e-50:
        return 0.0, 1.0

    # Transient threshold strain
    sigma_over_mu = abs(sigma_diff_Pa) / mu if mu > 0 else 0.0
    eps_t_star = K0 * np.exp(c * T) * (sigma_over_mu ** m_md)

    if eps_t_star < 1e-50:
        return eps_ss, 1.0

    # Work-hardening parameter Delta
    log_arg = max(sigma_over_mu, 1e-30)
    Delta = alpha_w + beta_w * np.log10(log_arg)

    # Transient function F
    ratio = zeta / eps_t_star
    if ratio <= 1.0:
        # Hardening branch
        F = np.exp(Delta * (1.0 - ratio) ** 2)
    else:
        # Recovery branch
        F = np.exp(-delta * (1.0 - ratio) ** 2)

    eps_rate = F * eps_ss
    return eps_rate, F


# ══════════════════════════════════════════════════════════════════════════════
# TIME INTEGRATION
# ══════════════════════════════════════════════════════════════════════════════

def _update_alpha(zeta, a1, alpha0, eta_vp):
    """Compute hardening variable α from accumulated VP strain ζ."""
    base = (a1 / alpha0) ** (1.0 / eta_vp) + zeta
    if base > 0:
        return max(a1 / (base ** eta_vp), 1e-30)
    return alpha0


# Maximum Desai strain increment per sub-step (prevents forward-Euler blow-up)
_DESAI_MAX_DEPS = 1e-4   # 0.01% strain per sub-step

# Maximum zeta change per sub-step for Munson-Dawson (prevents timestep dependency)
_MD_MAX_DZETA = 1e-4


def integrate_safeincave(t_s, sigma_Pa, sigma3_Pa_const):
    """
    Forward-Euler integration of the SafeInCave model:
      Spring (elastic, instant) + Kelvin (transient) + Dislocation (steady) + Desai (viscoplastic)

    The Desai component uses adaptive sub-stepping: if the predicted strain
    increment exceeds _DESAI_MAX_DEPS, the timestep is subdivided so that
    the hardening variable α updates gradually (avoiding the single-step
    blow-up that plagued the original implementation).

    Returns:
        eps_total  — total axial strain [-] at each time step
        eps_disloc — dislocation component
        eps_kelvin — Kelvin component
        eps_desai  — Desai component
    """
    n_steps = len(t_s)
    eps_elastic = np.zeros(n_steps)
    eps_disloc = np.zeros(n_steps)
    eps_kelvin_arr = np.zeros(n_steps)
    eps_desai_arr = np.zeros(n_steps)

    eps_k = 0.0       # current Kelvin strain
    eps_d = 0.0       # accumulated dislocation strain
    eps_vp = 0.0      # accumulated Desai viscoplastic strain
    alpha = ALPHA0_DESAI
    zeta_desai = 0.0  # accumulated VP strain for hardening

    for i in range(n_steps):
        sigma = sigma_Pa[i]
        dt = (t_s[i] - t_s[i - 1]) if i > 0 else 0.0

        # Elastic (instantaneous)
        eps_el = sigma / E_ELASTIC

        if dt > 0:
            # Dislocation creep
            rate_dc = disloc_rate(sigma, A_SIC, N_SIC, Q_DISLOC, T_KELVIN)
            eps_d += rate_dc * dt

            # Kelvin element
            rate_kv = kelvin_rate(eps_k, sigma, ETA_KELVIN, E1_KELVIN)
            eps_k += rate_kv * dt

            # Desai viscoplastic — adaptive sub-stepping
            if MU1_DESAI > 0:
                dt_remaining = dt
                while dt_remaining > 1e-10:
                    rate_desai = desai_rate(
                        sigma, sigma3_Pa_const, alpha,
                        MU1_DESAI, N1_DESAI, A1_DESAI, ETA_DESAI,
                        N_DESAI_POWER, BETA1_DESAI, BETA_DESAI, M_DESAI,
                        GAMMA_DESAI, SIGMA_T_DESAI
                    )
                    if rate_desai < 1e-30:
                        break

                    # Choose sub-step: limit strain increment to _DESAI_MAX_DEPS
                    dt_sub = min(dt_remaining, _DESAI_MAX_DEPS / rate_desai)
                    dt_sub = max(dt_sub, 1e-6)  # floor to avoid infinite loop

                    eps_vp += rate_desai * dt_sub
                    zeta_desai += rate_desai * dt_sub
                    alpha = _update_alpha(zeta_desai, A1_DESAI,
                                          ALPHA0_DESAI, ETA_DESAI)

                    dt_remaining -= dt_sub

        eps_elastic[i] = eps_el
        eps_disloc[i] = eps_d
        eps_kelvin_arr[i] = eps_k
        eps_desai_arr[i] = eps_vp

    eps_total = eps_elastic + eps_disloc + eps_kelvin_arr + eps_desai_arr
    return eps_total, eps_disloc, eps_kelvin_arr, eps_desai_arr


def integrate_munsondawson(t_s, sigma_Pa):
    """
    Forward-Euler integration of the Munson-Dawson model:
      Spring (elastic, instant) + MunsonDawson (steady + transient)

    Returns:
        eps_total    — total axial strain [-] at each time step
        eps_steady   — steady-state creep component
        eps_transient — transient component (total - steady - elastic)
    """
    n_steps = len(t_s)
    eps_elastic = np.zeros(n_steps)
    eps_creep = np.zeros(n_steps)
    eps_steady_arr = np.zeros(n_steps)

    eps_c = 0.0       # accumulated total creep strain
    eps_s = 0.0       # accumulated steady-state strain
    zeta = 0.0        # transient internal variable

    for i in range(n_steps):
        sigma = sigma_Pa[i]
        dt = (t_s[i] - t_s[i - 1]) if i > 0 else 0.0

        # Elastic
        eps_el = sigma / E_ELASTIC

        if dt > 0:
            # Steady-state rate (for decomposition)
            rate_ss = disloc_rate(sigma, A_MD, N_MD, Q_DISLOC, T_KELVIN)
            eps_s += rate_ss * dt

            # Adaptive sub-stepping for Munson-Dawson transient
            # (prevents timestep-dependent results when zeta changes rapidly)
            dt_remaining = dt
            while dt_remaining > 1e-10:
                rate_md, F = munsondawson_rate(
                    sigma, zeta, A_MD, N_MD, Q_DISLOC, T_KELVIN,
                    K0_MD, C_MD, M_MD, ALPHA_W_MD, BETA_W_MD,
                    DELTA_MD, MU_SHEAR
                )
                # Limit zeta change per sub-step
                zeta_rate = abs(F - 1.0) * rate_ss
                if zeta_rate > 1e-30:
                    dt_sub = min(dt_remaining, _MD_MAX_DZETA / zeta_rate)
                else:
                    dt_sub = dt_remaining
                dt_sub = max(dt_sub, 1e-6)
                dt_sub = min(dt_sub, dt_remaining)

                eps_c += rate_md * dt_sub
                zeta += (F - 1.0) * rate_ss * dt_sub
                zeta = max(zeta, 0.0)
                dt_remaining -= dt_sub

        eps_elastic[i] = eps_el
        eps_creep[i] = eps_c
        eps_steady_arr[i] = eps_s

    eps_total = eps_elastic + eps_creep
    eps_transient = eps_creep - eps_steady_arr
    return eps_total, eps_steady_arr, eps_transient


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Discover available test files
    available = {}
    for fname in sorted(os.listdir(DATA_DIR)):
        if fname.startswith("ZW_TCC") and fname.endswith(".csv"):
            test_id = fname.replace("ZW_", "").replace(".csv", "")
            available[test_id] = os.path.join(DATA_DIR, fname)

    tests_to_run = TESTS if TESTS else list(available.keys())

    print("=" * 70)
    print("CREEP TEST CALIBRATION")
    print("=" * 70)
    print(f"  Temperature:        {T_CELSIUS} °C ({T_KELVIN} K)")
    print(f"  Q/R (fixed):        {Q_OVER_R} K")
    print(f"  dt:                 {DT_HOURS} hours")
    print(f"  SafeInCave:         {'ON' if RUN_SAFEINCAVE else 'OFF'}")
    print(f"    A = {A_SIC:.4e}, n = {N_SIC}")
    print(f"    η_Kelvin = {ETA_KELVIN:.2e}, E1 = {E1_KELVIN:.2e}, ν1 = {NU1_KELVIN}")
    print(f"    μ1_Desai = {MU1_DESAI:.4e}, N1 = {N1_DESAI}, a1 = {A1_DESAI:.4e}")
    print(f"    η_Desai = {ETA_DESAI}, α0 = {ALPHA0_DESAI}")
    print(f"  Munson-Dawson:      {'ON' if RUN_MUNSONDAWSON else 'OFF'}")
    print(f"    A = {A_MD:.4e}, n = {N_MD}")
    print(f"    K0 = {K0_MD:.2e}, c = {C_MD}, m = {M_MD}, α_w = {ALPHA_W_MD}")
    print(f"  Tests:              {tests_to_run}")
    print("=" * 70)

    for test_id in tests_to_run:
        if test_id not in available:
            print(f"[WARN] {test_id} not found in {DATA_DIR}, skipping")
            continue

        print(f"\n[PROCESSING] {test_id}")
        filepath = available[test_id]

        # Read lab data
        time_h, strain_pct, sigma_diff, sigma3, temp_C = read_creep_csv(filepath)
        print(f"  Loaded {len(time_h)} data points, {time_h[-1]:.1f} hours")
        print(f"  σ_diff range: {sigma_diff.min():.1f} – {sigma_diff.max():.1f} MPa")
        print(f"  σ_3 range:    {sigma3.min():.1f} – {sigma3.max():.1f} MPa")

        # Apply artefact correction (fills brief dip-and-return events)
        sigma_diff_corr = correct_stress_artefacts(time_h, sigma_diff)

        # Build stress schedule on fine grid
        t_s, sigma_fine_Pa, idx_start = build_stress_schedule(
            time_h, sigma_diff_corr, DT_HOURS
        )
        sigma3_Pa = float(sigma3[idx_start]) * 1e6  # constant confining pressure

        # Shift lab data to same t=0
        t0_h = time_h[idx_start]
        lab_time_h = time_h[idx_start:] - t0_h
        lab_strain_pct = strain_pct[idx_start:]
        # Subtract the initial elastic strain so lab starts at 0
        lab_strain_pct = lab_strain_pct - lab_strain_pct[0]
        lab_sigma_diff = sigma_diff_corr[idx_start:]

        result = {
            "test_id": test_id,
            "temperature_C": T_CELSIUS,
            "sigma3_MPa": float(sigma3[idx_start]),
            "lab": {
                "time_hours": lab_time_h.tolist(),
                "strain_pct": lab_strain_pct.tolist(),
                "sigma_diff_MPa": lab_sigma_diff.tolist(),
            },
            "params": {
                "Q_over_R": Q_OVER_R,
                "A_mpa_yr": _A_MPA_YR,
                "n_disloc": N_DISLOC,
                "E_elastic": E_ELASTIC,
                "nu_elastic": NU_ELASTIC,
            }
        }

        if RUN_SAFEINCAVE:
            print("  Running SafeInCave model...")
            eps_tot, eps_dc, eps_kv, eps_desai = integrate_safeincave(
                t_s, sigma_fine_Pa, sigma3_Pa
            )  # uses corrected sigma schedule
            t_hours = t_s / HOUR
            result["safeincave"] = {
                "time_hours": t_hours.tolist(),
                "strain_total_pct": (eps_tot * 100).tolist(),
                "strain_disloc_pct": (eps_dc * 100).tolist(),
                "strain_kelvin_pct": (eps_kv * 100).tolist(),
                "strain_desai_pct": (eps_desai * 100).tolist(),
                "params": {
                    "A": A_SIC, "n": N_SIC,
                    "eta_kelvin": ETA_KELVIN, "E1_kelvin": E1_KELVIN,
                    "nu1_kelvin": NU1_KELVIN,
                    "mu1_desai": MU1_DESAI, "N1_desai": N1_DESAI,
                    "a1_desai": A1_DESAI, "eta_desai": ETA_DESAI,
                    "alpha0_desai": ALPHA0_DESAI,
                    "n_desai": N_DESAI_POWER, "beta1_desai": BETA1_DESAI,
                    "beta_desai": BETA_DESAI, "m_desai": M_DESAI,
                    "gamma_desai": GAMMA_DESAI, "sigma_t_desai": SIGMA_T_DESAI,
                },
            }
            print(f"    Final strain: {eps_tot[-1] * 100:.3f}%")

        if RUN_MUNSONDAWSON:
            print("  Running Munson-Dawson model...")
            eps_tot, eps_ss, eps_tr = integrate_munsondawson(t_s, sigma_fine_Pa)
            t_hours = t_s / HOUR
            result["munsondawson"] = {
                "time_hours": t_hours.tolist(),
                "strain_total_pct": (eps_tot * 100).tolist(),
                "strain_steady_pct": (eps_ss * 100).tolist(),
                "strain_transient_pct": (eps_tr * 100).tolist(),
                "params": {
                    "A": A_MD, "n": N_MD,
                    "K0": K0_MD, "c": C_MD, "m": M_MD,
                    "alpha_w": ALPHA_W_MD,
                    "beta_w": BETA_W_MD, "delta": DELTA_MD,
                },
            }
            print(f"    Final strain: {eps_tot[-1] * 100:.3f}%")

        # Save
        outpath = os.path.join(OUT_DIR, f"calibration_{test_id}.json")
        with open(outpath, "w") as f:
            json.dump(result, f, indent=2)
        print(f"  [SAVED] {outpath}")

    print("\n[DONE] Run plot_calibration.py to visualize results.")


if __name__ == "__main__":
    main()
