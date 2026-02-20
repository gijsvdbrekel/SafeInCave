"""
Material parameter calibration: single-element creep simulation for new lab data.

Reads the multi-stage cyclic triaxial creep test from data_processed.xlsx
(5 loading stages with unloading to isostatic between each), then integrates
the creep strain using either:
  (a) SafeInCave model  (Spring + Kelvin + DislocationCreep + Desai)
  (b) Munson-Dawson model (Spring + MunsonDawsonCreep)

Temperature: 21 C (294.15 K) — room temperature test from Herminio's dataset.
Confining pressure: ~12 MPa (read from data).

Output: JSON file with time vs. axial & radial strain for each model,
used by plot_newdata.py for comparison against lab data.
"""

import os
import json
import numpy as np

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          USER CONFIGURATION                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── SELECT WHICH MODEL(S) TO RUN ─────────────────────────────────────────────
RUN_SAFEINCAVE = True
RUN_MUNSONDAWSON = True

# ── TEMPERATURE (room temperature test) ──────────────────────────────────────
T_CELSIUS = 21.0
T_KELVIN = T_CELSIUS + 273.15  # 294.15 K

# ── TIME INTEGRATION ─────────────────────────────────────────────────────────
DT_HOURS = 0.25       # integration time step (hours) — finer for cyclic loading

# ── GAS CONSTANT ─────────────────────────────────────────────────────────────
R_GAS = 8.314        # J/(mol K)

# ══════════════════════════════════════════════════════════════════════════════
# CALIBRATION PARAMETERS — EDIT THESE TO CALIBRATE
# ══════════════════════════════════════════════════════════════════════════════

# ── SHARED: Dislocation creep ────────────────────────────────────────────────
#   eps_ss = A * sigma^n * exp(-Q/(R*T))
#   A is specified in [MPa^{-n} / yr] and converted to [Pa^{-n} / s]
Q_OVER_R = 6252                        # [K]  — literature range 6252-6495
Q_DISLOC = Q_OVER_R * R_GAS              # [J/mol]

_sec_per_year = 365.25 * 24 * 3600

# -- SafeInCave dislocation creep (calibrated independently) ------------------
_A_SIC_MPA_YR = 40.0       # [MPa^{-n}/yr]  — initial guess within [15, 60]
N_SIC         = 4.5         # [-]            — initial guess within [3, 6]
A_SIC         = _A_SIC_MPA_YR * (1e-6) ** N_SIC / _sec_per_year  # [Pa^{-n}/s]

# -- Munson-Dawson dislocation creep (calibrated independently) ---------------
_A_MD_MPA_YR  = 12285.55       # [MPa^{-n}/yr]  — initial guess within [15, 60]
N_MD          = 3.2468         # [-]            — initial guess within [3, 6]
A_MD          = _A_MD_MPA_YR * (1e-6) ** N_MD / _sec_per_year    # [Pa^{-n}/s]

# ── SafeInCave: Kelvin element (transient creep) ─────────────────────────────
ETA_KELVIN = 1e12       # [Pa s]   — initial guess
E1_KELVIN  = 5e8        # [Pa]     — initial guess
NU1_KELVIN = 0.25       # [-]

# ── SafeInCave: Desai viscoplastic ───────────────────────────────────────────
MU1_DESAI     = 1e-15      # [1/s]
N1_DESAI      = 2.0        # [-]
A1_DESAI      = 1e-03      # [-]
ETA_DESAI     = 1.2        # [-]
ALPHA0_DESAI  = 0.005      # [-]

# Fixed Desai parameters
N_DESAI_POWER = 3.0
BETA1_DESAI   = 0.0048
BETA_DESAI    = 0.995
M_DESAI       = -0.5
GAMMA_DESAI   = 0.095
SIGMA_T_DESAI = 5.0       # [Pa]

# ── Munson-Dawson: Transient parameters ──────────────────────────────────────
K0_MD      = 10       # [-]    — initial guess
C_MD       = 0.00902    # [1/K]  — fixed from literature
M_MD       = 1.4412        # [-]    — initial guess
ALPHA_W_MD = 2.2663       # [-]    — initial guess
BETA_W_MD  = -0.0001       # [-]    — initial guess
DELTA_MD   = 0.003281        # [-]    — initial guess

# ── Elastic parameters ──────────────────────────────────────────────────────
E_ELASTIC  = 20.425e9    # [Pa]
NU_ELASTIC = 0.25        # [-]
MU_SHEAR   = E_ELASTIC / (2.0 * (1.0 + NU_ELASTIC))  # ~8.17 GPa

# ── PATHS ────────────────────────────────────────────────────────────────────
XLSX_PATH = os.path.join(os.path.dirname(__file__),
                         "experimental_data", "data_processed.xlsx")
OUT_DIR   = os.path.join(os.path.dirname(__file__), "output")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                      END OF USER CONFIGURATION                             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


HOUR = 3600.0
DAY = 24.0 * HOUR

# Maximum strain increment per sub-step
_DESAI_MAX_DEPS = 1e-4
_MD_MAX_DZETA = 1e-4


# ══════════════════════════════════════════════════════════════════════════════
# EXCEL READING
# ══════════════════════════════════════════════════════════════════════════════

def read_excel_data(filepath):
    """
    Read the cyclic triaxial creep test from data_processed.xlsx.

    Columns: index, Time(hours), Sigma_1(MPa), Sigma_3(MPa),
             Epsilon_1(%), Epsilon_3(%)

    Returns:
        time_h      — elapsed time [hours]
        sigma1      — axial stress [MPa]
        sigma3      — confining pressure [MPa]
        epsilon1    — axial strain [%]
        epsilon3    — radial strain [%]
    """
    import openpyxl
    wb = openpyxl.load_workbook(filepath, read_only=True)
    ws = wb.active

    time_h, sigma1, sigma3, epsilon1, epsilon3 = [], [], [], [], []

    for i, row in enumerate(ws.iter_rows(min_row=2, values_only=True)):
        if row[1] is None:
            continue
        time_h.append(float(row[1]))
        sigma1.append(float(row[2]))
        sigma3.append(float(row[3]))
        epsilon1.append(float(row[4]))
        epsilon3.append(float(row[5]))

    wb.close()
    return (np.array(time_h), np.array(sigma1), np.array(sigma3),
            np.array(epsilon1), np.array(epsilon3))


def build_stress_schedule(time_h, sigma1, sigma3, dt_hours):
    """
    Build a fine time grid with stress interpolated from the lab data.

    For cyclic loading we need both sigma1 and sigma3 at each time step.
    We use piecewise-linear interpolation to capture the loading/unloading
    transitions accurately.

    Returns:
        t_s         — time array [seconds]
        sigma1_Pa   — axial stress [Pa] at each time step
        sigma3_Pa   — confining pressure [Pa] at each time step
        idx_start   — index in original data where deviatoric loading begins
    """
    # Find the first time point with significant differential stress (> 1 MPa)
    sigma_diff = sigma1 - sigma3
    mask = sigma_diff > 1.0
    if not np.any(mask):
        raise ValueError("No deviatoric loading found in data")

    idx_start = np.argmax(mask)
    t0_h = time_h[idx_start]
    t_end_h = time_h[-1]

    # Lab time/stress pairs (shifted so deviatoric loading starts at t=0)
    t_lab_s = (time_h[idx_start:] - t0_h) * HOUR

    # Fine time grid
    n_steps = int(np.ceil((t_end_h - t0_h) / dt_hours))
    t_s = np.linspace(0.0, (t_end_h - t0_h) * HOUR, n_steps + 1)

    # Piecewise-linear interpolation for both stress components
    sigma1_Pa = np.interp(t_s, t_lab_s, sigma1[idx_start:] * 1e6)
    sigma3_Pa = np.interp(t_s, t_lab_s, sigma3[idx_start:] * 1e6)

    return t_s, sigma1_Pa, sigma3_Pa, idx_start


# ══════════════════════════════════════════════════════════════════════════════
# CREEP RATE FUNCTIONS (scalar, single-element)
# ══════════════════════════════════════════════════════════════════════════════

def disloc_rate(sigma_diff_Pa, A, n, Q, T):
    """Dislocation creep axial strain rate [1/s]."""
    return A * (abs(sigma_diff_Pa) ** n) * np.exp(-Q / (R_GAS * T))


def kelvin_rate(eps_kelvin, sigma_diff_Pa, eta, E1):
    """Kelvin element axial strain rate [1/s]."""
    return (sigma_diff_Pa - E1 * eps_kelvin) / eta


def desai_rate(sigma_diff_Pa, sigma3_Pa, alpha, mu1, N1, a1, eta_vp,
               n_pow, beta1, beta, m, gamma, sigma_t):
    """Desai viscoplastic axial strain rate [1/s] (simplified scalar version)."""
    MPA = 1e6

    s1 = (sigma3_Pa + sigma_diff_Pa) / MPA
    s2 = sigma3_Pa / MPA
    s3_val = sigma3_Pa / MPA

    I1 = s1 + s2 + s3_val
    I2 = s1 * s2 + s2 * s3_val + s1 * s3_val
    I3 = s1 * s2 * s3_val
    J2 = (1.0 / 3.0) * I1**2 - I2
    J3 = (2.0 / 27.0) * I1**3 - (1.0 / 3.0) * I1 * I2 + I3

    if J2 < 1e-30:
        return 0.0

    sqJ2 = np.sqrt(J2)
    arg_sr = -(J3 * np.sqrt(27.0)) / (2.0 * sqJ2**3)
    arg_sr = np.clip(arg_sr, -1.0, 1.0)
    Sr = arg_sr

    I1_star = I1 + sigma_t

    F1 = alpha * I1_star**n_pow - gamma * I1_star**2
    F2_base = np.exp(beta1 * I1_star) - beta * Sr
    if F2_base <= 0.0:
        F_vp = J2
    else:
        F_vp = J2 + F1 * F2_base**m

    if F_vp <= 0.0:
        return 0.0

    lam = mu1 * (F_vp) ** N1

    p = I1 / 3.0
    dev1 = s1 - p

    eps_rate = lam * abs(dev1)
    return eps_rate


def _update_alpha(zeta, a1, alpha0, eta_vp):
    """Compute hardening variable alpha from accumulated VP strain zeta."""
    base = (a1 / alpha0) ** (1.0 / eta_vp) + zeta
    if base > 0:
        return max(a1 / (base ** eta_vp), 1e-30)
    return alpha0


def munsondawson_rate(sigma_diff_Pa, zeta, A, n, Q, T,
                      K0, c, m_md, alpha_w, beta_w, delta, mu):
    """
    Munson-Dawson creep axial strain rate [1/s] and transient multiplier.

    Returns (eps_rate, F_value).
    """
    eps_ss = A * (abs(sigma_diff_Pa) ** n) * np.exp(-Q / (R_GAS * T))

    if eps_ss < 1e-50:
        return 0.0, 1.0

    sigma_over_mu = abs(sigma_diff_Pa) / mu if mu > 0 else 0.0
    eps_t_star = K0 * np.exp(c * T) * (sigma_over_mu ** m_md)

    if eps_t_star < 1e-50:
        return eps_ss, 1.0

    log_arg = max(sigma_over_mu, 1e-30)
    Delta = alpha_w + beta_w * np.log10(log_arg)

    ratio = zeta / eps_t_star
    if ratio <= 1.0:
        F = np.exp(Delta * (1.0 - ratio) ** 2)
    else:
        F = np.exp(-delta * (1.0 - ratio) ** 2)

    eps_rate = F * eps_ss
    return eps_rate, F


# ══════════════════════════════════════════════════════════════════════════════
# TIME INTEGRATION
# ══════════════════════════════════════════════════════════════════════════════

def integrate_safeincave(t_s, sigma1_Pa, sigma3_Pa,
                         A=None, n=None, Q=None, T=None,
                         eta=None, E1=None,
                         mu1=None, N1=None, a1=None, eta_vp=None, alpha0=None):
    """
    Forward-Euler integration of the SafeInCave model with exact Kelvin update.

    Handles cyclic loading: both sigma1 and sigma3 vary with time.

    Returns:
        eps_total_ax  — total axial strain [-]
        eps_disloc    — dislocation component
        eps_kelvin    — Kelvin component
        eps_desai     — Desai component
        eps_total_rad — total radial strain [-] (elastic only for now)
    """
    # Default to module-level parameters
    if A is None: A = A_SIC
    if n is None: n = N_SIC
    if Q is None: Q = Q_DISLOC
    if T is None: T = T_KELVIN
    if eta is None: eta = ETA_KELVIN
    if E1 is None: E1 = E1_KELVIN
    if mu1 is None: mu1 = MU1_DESAI
    if N1 is None: N1 = N1_DESAI
    if a1 is None: a1 = A1_DESAI
    if eta_vp is None: eta_vp = ETA_DESAI
    if alpha0 is None: alpha0 = ALPHA0_DESAI

    n_steps = len(t_s)
    eps_elastic = np.zeros(n_steps)
    eps_elastic_rad = np.zeros(n_steps)
    eps_disloc = np.zeros(n_steps)
    eps_kelvin_arr = np.zeros(n_steps)
    eps_desai_arr = np.zeros(n_steps)

    eps_k = 0.0
    eps_d = 0.0
    eps_vp = 0.0
    alpha = alpha0
    zeta_desai = 0.0

    for i in range(n_steps):
        sigma_diff = sigma1_Pa[i] - sigma3_Pa[i]
        s3 = sigma3_Pa[i]
        dt = (t_s[i] - t_s[i - 1]) if i > 0 else 0.0

        # Elastic (instantaneous) — axial
        eps_el_ax = sigma_diff / E_ELASTIC
        # Elastic — radial: -nu * (sigma1 - sigma3) / E
        # More precisely for triaxial: eps_rad = (sigma3 - nu*(sigma1 + sigma3)) / E
        # but for differential: eps_rad_diff = -nu * sigma_diff / E
        eps_el_rad = -NU_ELASTIC * sigma_diff / E_ELASTIC

        if dt > 0:
            # Dislocation creep
            rate_dc = disloc_rate(sigma_diff, A, n, Q, T)
            eps_d += rate_dc * dt

            # Kelvin element — exact update (exponential integrator)
            eps_k_eq = sigma_diff / E1 if E1 > 0 else 0.0
            decay = np.exp(-E1 * dt / eta) if (eta > 0 and E1 > 0) else 0.0
            eps_k = eps_k_eq + (eps_k - eps_k_eq) * decay

            # Desai viscoplastic — adaptive sub-stepping
            if mu1 > 0:
                dt_remaining = dt
                while dt_remaining > 1e-10:
                    rate_desai = desai_rate(
                        sigma_diff, s3, alpha,
                        mu1, N1, a1, eta_vp,
                        N_DESAI_POWER, BETA1_DESAI, BETA_DESAI, M_DESAI,
                        GAMMA_DESAI, SIGMA_T_DESAI
                    )
                    if rate_desai < 1e-30:
                        break

                    dt_sub = min(dt_remaining, _DESAI_MAX_DEPS / rate_desai)
                    dt_sub = max(dt_sub, 1e-6)

                    eps_vp += rate_desai * dt_sub
                    zeta_desai += rate_desai * dt_sub
                    alpha = _update_alpha(zeta_desai, a1, alpha0, eta_vp)

                    dt_remaining -= dt_sub

        eps_elastic[i] = eps_el_ax
        eps_elastic_rad[i] = eps_el_rad
        eps_disloc[i] = eps_d
        eps_kelvin_arr[i] = eps_k
        eps_desai_arr[i] = eps_vp

    eps_total_ax = eps_elastic + eps_disloc + eps_kelvin_arr + eps_desai_arr
    # Radial: elastic + creep (assume volume-preserving creep: eps_rad_creep = -0.5 * eps_ax_creep)
    eps_creep_ax = eps_disloc + eps_kelvin_arr + eps_desai_arr
    eps_total_rad = eps_elastic_rad - 0.5 * eps_creep_ax

    return eps_total_ax, eps_disloc, eps_kelvin_arr, eps_desai_arr, eps_total_rad


def integrate_munsondawson(t_s, sigma1_Pa, sigma3_Pa,
                           A=None, n=None, Q=None, T=None,
                           K0=None, c=None, m_md=None,
                           alpha_w=None, beta_w=None, delta=None):
    """
    Forward-Euler integration of the Munson-Dawson model with adaptive sub-stepping.

    Handles cyclic loading: both sigma1 and sigma3 vary with time.

    Returns:
        eps_total_ax  — total axial strain [-]
        eps_steady    — steady-state creep component
        eps_transient — transient component
        eps_total_rad — total radial strain [-]
    """
    if A is None: A = A_MD
    if n is None: n = N_MD
    if Q is None: Q = Q_DISLOC
    if T is None: T = T_KELVIN
    if K0 is None: K0 = K0_MD
    if c is None: c = C_MD
    if m_md is None: m_md = M_MD
    if alpha_w is None: alpha_w = ALPHA_W_MD
    if beta_w is None: beta_w = BETA_W_MD
    if delta is None: delta = DELTA_MD

    n_steps = len(t_s)
    eps_elastic = np.zeros(n_steps)
    eps_elastic_rad = np.zeros(n_steps)
    eps_creep = np.zeros(n_steps)
    eps_steady_arr = np.zeros(n_steps)

    eps_c = 0.0
    eps_s = 0.0
    zeta = 0.0

    for i in range(n_steps):
        sigma_diff = sigma1_Pa[i] - sigma3_Pa[i]
        dt = (t_s[i] - t_s[i - 1]) if i > 0 else 0.0

        eps_el_ax = sigma_diff / E_ELASTIC
        eps_el_rad = -NU_ELASTIC * sigma_diff / E_ELASTIC

        if dt > 0:
            rate_ss = disloc_rate(sigma_diff, A, n, Q, T)
            eps_s += rate_ss * dt

            if rate_ss < 1e-60:
                eps_elastic[i] = eps_el_ax
                eps_elastic_rad[i] = eps_el_rad
                eps_creep[i] = eps_c
                eps_steady_arr[i] = eps_s
                continue

            dt_remaining = dt
            while dt_remaining > 1e-10:
                rate_md, F = munsondawson_rate(
                    sigma_diff, zeta, A, n, Q, T,
                    K0, c, m_md, alpha_w, beta_w, delta, MU_SHEAR
                )
                zeta_rate = abs(F - 1.0) * rate_ss
                if zeta_rate > 1e-30:
                    dt_sub = min(dt_remaining, _MD_MAX_DZETA / zeta_rate)
                else:
                    dt_sub = dt_remaining
                dt_sub = max(dt_sub, 10.0)
                dt_sub = min(dt_sub, dt_remaining)

                eps_c += rate_md * dt_sub
                zeta += (F - 1.0) * rate_ss * dt_sub
                zeta = max(zeta, 0.0)
                dt_remaining -= dt_sub

        eps_elastic[i] = eps_el_ax
        eps_elastic_rad[i] = eps_el_rad
        eps_creep[i] = eps_c
        eps_steady_arr[i] = eps_s

    eps_total_ax = eps_elastic + eps_creep
    eps_transient = eps_creep - eps_steady_arr
    # Radial strain
    eps_total_rad = eps_elastic_rad - 0.5 * eps_creep

    return eps_total_ax, eps_steady_arr, eps_transient, eps_total_rad


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=" * 70)
    print("CYCLIC CREEP TEST CALIBRATION — NEW LAB DATA")
    print("=" * 70)
    print(f"  Data:               {os.path.basename(XLSX_PATH)}")
    print(f"  Temperature:        {T_CELSIUS} C ({T_KELVIN} K)")
    print(f"  Q/R (fixed):        {Q_OVER_R} K")
    print(f"  dt:                 {DT_HOURS} hours")

    # Read lab data
    print("\n[READING] Excel data...")
    time_h, sigma1, sigma3, epsilon1, epsilon3 = read_excel_data(XLSX_PATH)
    sigma_diff = sigma1 - sigma3
    print(f"  Loaded {len(time_h)} data points, {time_h[-1]:.1f} hours")
    print(f"  sigma_diff range: {sigma_diff.min():.1f} - {sigma_diff.max():.1f} MPa")
    print(f"  sigma_3 range:    {sigma3.min():.1f} - {sigma3.max():.1f} MPa")
    print(f"  epsilon_1 range:  {epsilon1.min():.4f} - {epsilon1.max():.4f} %")
    print(f"  epsilon_3 range:  {epsilon3.min():.4f} - {epsilon3.max():.4f} %")

    # Build stress schedule on fine grid
    t_s, sigma1_fine, sigma3_fine, idx_start = build_stress_schedule(
        time_h, sigma1, sigma3, DT_HOURS
    )

    # Shift lab data to same t=0
    t0_h = time_h[idx_start]
    lab_time_h = time_h[idx_start:] - t0_h
    lab_eps1 = epsilon1[idx_start:] - epsilon1[idx_start]
    lab_eps3 = epsilon3[idx_start:] - epsilon3[idx_start]
    lab_sigma_diff = sigma_diff[idx_start:]
    lab_sigma1 = sigma1[idx_start:]
    lab_sigma3 = sigma3[idx_start:]

    result = {
        "test_id": "newdata_cyclic",
        "temperature_C": T_CELSIUS,
        "sigma3_MPa_mean": float(np.mean(sigma3[idx_start:])),
        "lab": {
            "time_hours": lab_time_h.tolist(),
            "strain_axial_pct": lab_eps1.tolist(),
            "strain_radial_pct": lab_eps3.tolist(),
            "sigma_diff_MPa": lab_sigma_diff.tolist(),
            "sigma1_MPa": lab_sigma1.tolist(),
            "sigma3_MPa": lab_sigma3.tolist(),
        },
        "params": {
            "Q_over_R": Q_OVER_R,
            "E_elastic": E_ELASTIC,
            "nu_elastic": NU_ELASTIC,
        }
    }

    if RUN_SAFEINCAVE:
        print(f"\n  SafeInCave:         ON")
        print(f"    A = {_A_SIC_MPA_YR:.2f} MPa^-n/yr, n = {N_SIC:.4f}")
        print(f"    eta_Kelvin = {ETA_KELVIN:.2e}, E1 = {E1_KELVIN:.2e}")
        print(f"    mu1_Desai = {MU1_DESAI:.2e}, N1 = {N1_DESAI}, a1 = {A1_DESAI:.2e}")
        print(f"    eta_Desai = {ETA_DESAI}, alpha0 = {ALPHA0_DESAI}")
        print("  Running SafeInCave model...")

        eps_ax, eps_dc, eps_kv, eps_desai, eps_rad = integrate_safeincave(
            t_s, sigma1_fine, sigma3_fine
        )
        t_hours = t_s / HOUR
        result["safeincave"] = {
            "time_hours": t_hours.tolist(),
            "strain_axial_pct": (eps_ax * 100).tolist(),
            "strain_radial_pct": (eps_rad * 100).tolist(),
            "strain_disloc_pct": (eps_dc * 100).tolist(),
            "strain_kelvin_pct": (eps_kv * 100).tolist(),
            "strain_desai_pct": (eps_desai * 100).tolist(),
            "params": {
                "A_mpa_yr": _A_SIC_MPA_YR, "n": N_SIC,
                "eta_kelvin": ETA_KELVIN, "E1_kelvin": E1_KELVIN,
                "nu1_kelvin": NU1_KELVIN,
                "mu1_desai": MU1_DESAI, "N1_desai": N1_DESAI,
                "a1_desai": A1_DESAI, "eta_desai": ETA_DESAI,
                "alpha0_desai": ALPHA0_DESAI,
            },
        }
        print(f"    Final axial strain:  {eps_ax[-1] * 100:.4f}%")
        print(f"    Final radial strain: {eps_rad[-1] * 100:.4f}%")

    if RUN_MUNSONDAWSON:
        print(f"\n  Munson-Dawson:      ON")
        print(f"    A = {_A_MD_MPA_YR:.2f} MPa^-n/yr, n = {N_MD:.4f}")
        print(f"    K0 = {K0_MD:.4f}, m = {M_MD}, alpha_w = {ALPHA_W_MD}")
        print(f"    beta_w = {BETA_W_MD}, delta = {DELTA_MD}")
        print("  Running Munson-Dawson model...")

        eps_ax, eps_ss, eps_tr, eps_rad = integrate_munsondawson(
            t_s, sigma1_fine, sigma3_fine
        )
        t_hours = t_s / HOUR
        result["munsondawson"] = {
            "time_hours": t_hours.tolist(),
            "strain_axial_pct": (eps_ax * 100).tolist(),
            "strain_radial_pct": (eps_rad * 100).tolist(),
            "strain_steady_pct": (eps_ss * 100).tolist(),
            "strain_transient_pct": (eps_tr * 100).tolist(),
            "params": {
                "A_mpa_yr": _A_MD_MPA_YR, "n": N_MD,
                "K0": K0_MD, "c": C_MD, "m": M_MD,
                "alpha_w": ALPHA_W_MD,
                "beta_w": BETA_W_MD, "delta": DELTA_MD,
            },
        }
        print(f"    Final axial strain:  {eps_ax[-1] * 100:.4f}%")
        print(f"    Final radial strain: {eps_rad[-1] * 100:.4f}%")

    # Save
    outpath = os.path.join(OUT_DIR, "calibration_newdata.json")
    with open(outpath, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\n  [SAVED] {outpath}")
    print("\n[DONE] Run plot_newdata.py to visualize results.")


if __name__ == "__main__":
    main()
