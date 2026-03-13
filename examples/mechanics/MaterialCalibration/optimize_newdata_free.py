"""
Material parameter optimizer with free Desai yield surface parameters.

Based on optimize_newdata.py but additionally calibrates the Desai yield surface
parameters (n_desai, beta1, beta, m_desai, gamma) that were previously fixed,
and the MD temperature constant c.

Elastic parameters remain fixed at E=20.425 GPa, nu=0.25.
Dislocation creep (A, n, Q/R) fitted from steady-state slopes without bounds.
The only constraint: A_SIC = 1.5 × A_MD (2/3 flow rule factor).

Strategy:
  Phase 1: Fit dislocation creep (A, n, Q/R) — literature-bounded.
  Phase 2: Optimize transient parameters via differential_evolution.
           SIC: Kelvin (eta, E1) + Desai (all 12 params free)
           MD:  All transient params including c (K0, c, m, alpha_w, beta_w, delta)

Usage:
    python optimize_newdata_free.py
    SKIP_PHASE1=1 python optimize_newdata_free.py  # skip Phase 1
"""

import os
import sys
import math
import time as _time
import numpy as np
import numba
from scipy.optimize import differential_evolution
from scipy.stats import linregress

sys.path.insert(0, os.path.dirname(__file__))
from calibrate_newdata import (
    read_excel_data, build_stress_schedule,
    R_GAS, E_ELASTIC, NU_ELASTIC, MU_SHEAR,
    HOUR,
    _DESAI_MAX_DEPS, _MD_MAX_DZETA,
    XLSX_PATH,
)


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          CONFIGURATION                                     ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

SKIP_PHASE1 = os.environ.get("SKIP_PHASE1", "0") == "1"

# Temperature
T_CELSIUS = 21.0
T_KELVIN = T_CELSIUS + 273.15

# Time integration
DT_HOURS = 0.5

# Resampling
N_PER_STAGE = 40

# Zone weights for cyclic loading
W_LOADING_TRANSIENT  = 5.0
W_LOADING_STEADY     = 1.5
W_UNLOADING_TRANSIENT = 3.0
W_UNLOADING_STEADY   = 0.5

T_LOADING_TRANS  = 24.0
T_UNLOADING_TRANS = 12.0

# Phase 1: dislocation creep — unconstrained wide ranges
QR_GRID = np.linspace(4000.0, 10000.0, 25)  # K
A_BOUNDS = (0.1, 1e6)         # MPa^-n/yr
N_BOUNDS = (1.0, 10.0)        # [-]
STEADY_STATE_TAIL_FRAC = 0.3

_SEC_PER_YR = 365.25 * 24 * 3600


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                      STAGE DETECTION FOR CYCLIC LOADING                    ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

def detect_cyclic_stages(time_h, sigma_diff, min_change_mpa=5.0, min_hold_h=20.0):
    LOADING_THRESH = 5.0
    stages = []
    n_init = min(100, len(sigma_diff))
    current_sigma = np.median(sigma_diff[:n_init])
    stage_start = time_h[0]

    for i in range(1, len(time_h)):
        if abs(sigma_diff[i] - current_sigma) > min_change_mpa:
            new_sigma = sigma_diff[i]
            hold_end = time_h[i] + min_hold_h
            held = True
            for j in range(i + 1, len(time_h)):
                if time_h[j] > hold_end:
                    break
                if abs(sigma_diff[j] - new_sigma) > min_change_mpa:
                    held = False
                    break
            else:
                if time_h[-1] - time_h[i] < 0.5 * min_hold_h:
                    held = False
            if held:
                stage_type = "loading" if current_sigma > LOADING_THRESH else "unloading"
                stages.append({'t_start': stage_start, 't_end': time_h[i],
                                'sigma_mpa': float(current_sigma), 'type': stage_type})
                current_sigma = new_sigma
                stage_start = time_h[i]

    stage_type = "loading" if current_sigma > LOADING_THRESH else "unloading"
    stages.append({'t_start': stage_start, 't_end': time_h[-1],
                    'sigma_mpa': float(current_sigma), 'type': stage_type})

    for stg in stages:
        t0, t1 = stg['t_start'], stg['t_end']
        t_mid = t0 + 0.5 * (t1 - t0)
        mask = (time_h >= t_mid) & (time_h <= t1)
        if np.any(mask):
            stg['sigma_mpa'] = float(np.median(sigma_diff[mask]))
        stg['type'] = "loading" if stg['sigma_mpa'] > LOADING_THRESH else "unloading"
    return stages


def resample_and_label_cyclic(time_h, strain_pct, sigma_diff, stages,
                               n_per_stage=N_PER_STAGE):
    t_all, e_all, s_all, w_all, sid_all = [], [], [], [], []
    for k, stg in enumerate(stages):
        t0, t1 = stg['t_start'], stg['t_end']
        mask = (time_h >= t0) & (time_h <= t1)
        t_stg, e_stg, s_stg = time_h[mask], strain_pct[mask], sigma_diff[mask]
        if len(t_stg) < 2:
            continue
        t_uni = np.linspace(t_stg[0], t_stg[-1], n_per_stage)
        e_uni = np.interp(t_uni, t_stg, e_stg)
        s_uni = np.interp(t_uni, t_stg, s_stg)
        offsets = t_uni - t_uni[0]
        weights = np.ones(n_per_stage)
        if stg['type'] == "loading":
            for i, dt_off in enumerate(offsets):
                weights[i] = W_LOADING_TRANSIENT if dt_off < T_LOADING_TRANS else W_LOADING_STEADY
        else:
            for i, dt_off in enumerate(offsets):
                weights[i] = W_UNLOADING_TRANSIENT if dt_off < T_UNLOADING_TRANS else W_UNLOADING_STEADY
        t_all.append(t_uni); e_all.append(e_uni); s_all.append(s_uni)
        w_all.append(weights); sid_all.append(np.full(n_per_stage, k, dtype=int))
    return (np.concatenate(t_all), np.concatenate(e_all),
            np.concatenate(s_all), np.concatenate(w_all), np.concatenate(sid_all))


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                    PHASE 1: DISLOCATION CREEP FROM SLOPES                  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

def extract_steady_state_rates(time_h, strain_pct, sigma_diff, stages,
                                tail_frac=STEADY_STATE_TAIL_FRAC):
    stress_mpa, rates_per_s = [], []
    for stg in stages:
        if stg['type'] != 'loading' or stg['sigma_mpa'] < 5.0:
            continue
        t0, t1 = stg['t_start'], stg['t_end']
        duration = t1 - t0
        if duration < 10.0:
            continue
        t_tail_start = t1 - tail_frac * duration
        t_tail_end = t1 - 2.0
        if t_tail_end <= t_tail_start:
            t_tail_end = t1
        mask = ((time_h >= t_tail_start) & (time_h <= t_tail_end) &
                (np.abs(sigma_diff - stg['sigma_mpa']) < 2.0))
        t_stg, e_stg = time_h[mask], strain_pct[mask]
        if len(t_stg) < 10:
            continue
        slope, _, _, _, _ = linregress(t_stg, e_stg)
        rate_per_s = slope / 100.0 / 3600.0
        if rate_per_s > 0:
            stress_mpa.append(stg['sigma_mpa'])
            rates_per_s.append(rate_per_s)
    return np.array(stress_mpa), np.array(rates_per_s)


def fit_dislocation_params(stress_mpa, rates_per_s, T, QR_values):
    sigma_Pa = stress_mpa * 1e6
    log_rate = np.log(rates_per_s)
    best_residual = np.inf
    best_params = None
    logA_grid = np.linspace(np.log10(A_BOUNDS[0]), np.log10(A_BOUNDS[1]), 80)
    n_grid = np.linspace(N_BOUNDS[0], N_BOUNDS[1], 50)
    for QR in QR_values:
        exp_QRT = np.exp(-QR / T)
        for n_val in n_grid:
            sigma_n = sigma_Pa ** n_val
            for logA in logA_grid:
                A_mpa = 10.0 ** logA
                A_Pa_s = A_mpa * (1e-6) ** n_val / _SEC_PER_YR
                predicted = A_Pa_s * sigma_n * exp_QRT
                if np.any(predicted <= 0):
                    continue
                residual = np.sum((log_rate - np.log(predicted)) ** 2)
                if residual < best_residual:
                    best_residual = residual
                    best_params = (A_mpa, n_val, QR)
    return best_params, best_residual


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                    PHASE 2: MODEL INTEGRATION (NUMBA JIT)                  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

_E_EL = E_ELASTIC
_MU_SH = MU_SHEAR
_R = R_GAS
_MAX_DEPS = _DESAI_MAX_DEPS
_MAX_DZETA = _MD_MAX_DZETA
_MAX_SUBSTEPS = 2000


@numba.njit(cache=True)
def integrate_sic_free(t_s, sigma_diff_Pa, A_pa_s, n, Q, T, eta, E1,
                       mu1, N1, a1, eta_vp, alpha0,
                       n_desai, beta1, beta, m_desai, gamma, sigma_t):
    """
    Numba-JIT SIC integration with ALL Desai parameters as inputs.
    Elastic params fixed at compile-time constants.
    """
    exp_QRT = math.exp(-Q / (_R * T))
    n_pts = len(t_s)
    eps_k = 0.0
    eps_d = 0.0
    eps_vp = 0.0
    alpha = alpha0
    zeta_desai = 0.0
    eps = np.zeros(n_pts)
    MPA = 1e6
    sigma3_approx = 12.0e6

    for i in range(n_pts):
        sigma = sigma_diff_Pa[i]
        dt = (t_s[i] - t_s[i - 1]) if i > 0 else 0.0
        e_el = sigma / _E_EL

        if dt > 0.0:
            rate_dc = (2.0 / 3.0) * A_pa_s * (abs(sigma) ** n) * exp_QRT
            eps_d += rate_dc * dt

            if E1 > 0.0 and eta > 0.0:
                eps_k_eq = sigma / E1
                decay = math.exp(-E1 * dt / eta)
                eps_k = eps_k_eq + (eps_k - eps_k_eq) * decay

            if mu1 > 0.0:
                dt_rem = dt
                n_sub = 0
                while dt_rem > 1e-10 and n_sub < _MAX_SUBSTEPS:
                    n_sub += 1
                    s1 = (sigma3_approx + sigma) / MPA
                    s2 = sigma3_approx / MPA
                    s3v = sigma3_approx / MPA
                    I1 = s1 + s2 + s3v
                    I2 = s1 * s2 + s2 * s3v + s1 * s3v
                    I3 = s1 * s2 * s3v
                    J2 = (1.0 / 3.0) * I1 * I1 - I2
                    J3 = (2.0 / 27.0) * I1 * I1 * I1 - (1.0 / 3.0) * I1 * I2 + I3
                    if J2 < 1e-30:
                        break
                    sqJ2 = math.sqrt(J2)
                    arg_sr = -(J3 * math.sqrt(27.0)) / (2.0 * sqJ2 * sqJ2 * sqJ2)
                    if arg_sr < -1.0:
                        arg_sr = -1.0
                    if arg_sr > 1.0:
                        arg_sr = 1.0
                    Sr = arg_sr
                    I1_star = I1 + sigma_t
                    F1 = alpha * I1_star ** n_desai - gamma * I1_star * I1_star
                    F2_base = math.exp(beta1 * I1_star) - beta * Sr
                    if F2_base <= 0.0:
                        F_vp = J2
                    else:
                        F_vp = J2 + F1 * F2_base ** m_desai
                    if F_vp <= 0.0:
                        break
                    lam = mu1 * (F_vp ** N1)
                    p = I1 / 3.0
                    dev1 = s1 - p
                    rate_desai = lam * abs(dev1)
                    if rate_desai < 1e-30:
                        break
                    dt_sub = min(dt_rem, _MAX_DEPS / rate_desai)
                    if dt_sub < 1e-6:
                        dt_sub = 1e-6
                    eps_vp += rate_desai * dt_sub
                    zeta_desai += rate_desai * dt_sub
                    base = (a1 / alpha0) ** (1.0 / eta_vp) + zeta_desai
                    if base > 0.0:
                        alpha = a1 / (base ** eta_vp)
                        if alpha < 1e-30:
                            alpha = 1e-30
                    else:
                        alpha = alpha0
                    dt_rem -= dt_sub

        eps[i] = (e_el + eps_d + eps_k + eps_vp) * 100.0

    return eps


@numba.njit(cache=True)
def integrate_sic_free_components(t_s, sigma_diff_Pa, A_pa_s, n, Q, T, eta, E1,
                                  mu1, N1, a1, eta_vp, alpha0,
                                  n_desai, beta1, beta, m_desai, gamma, sigma_t):
    """
    Same as integrate_sic_free but returns (eps_total, eps_disloc, eps_kelvin, eps_desai)
    all in % — for plotting component breakdown.
    """
    exp_QRT = math.exp(-Q / (_R * T))
    n_pts = len(t_s)
    eps_k = 0.0
    eps_d = 0.0
    eps_vp = 0.0
    alpha = alpha0
    zeta_desai = 0.0
    MPA = 1e6
    sigma3_approx = 12.0e6

    out_total = np.zeros(n_pts)
    out_dc = np.zeros(n_pts)
    out_kv = np.zeros(n_pts)
    out_de = np.zeros(n_pts)

    for i in range(n_pts):
        sigma = sigma_diff_Pa[i]
        dt = (t_s[i] - t_s[i - 1]) if i > 0 else 0.0
        e_el = sigma / _E_EL

        if dt > 0.0:
            rate_dc = (2.0 / 3.0) * A_pa_s * (abs(sigma) ** n) * exp_QRT
            eps_d += rate_dc * dt

            if E1 > 0.0 and eta > 0.0:
                eps_k_eq = sigma / E1
                decay = math.exp(-E1 * dt / eta)
                eps_k = eps_k_eq + (eps_k - eps_k_eq) * decay

            if mu1 > 0.0:
                dt_rem = dt
                n_sub = 0
                while dt_rem > 1e-10 and n_sub < _MAX_SUBSTEPS:
                    n_sub += 1
                    s1 = (sigma3_approx + sigma) / MPA
                    s2 = sigma3_approx / MPA
                    s3v = sigma3_approx / MPA
                    I1 = s1 + s2 + s3v
                    I2 = s1 * s2 + s2 * s3v + s1 * s3v
                    I3 = s1 * s2 * s3v
                    J2 = (1.0 / 3.0) * I1 * I1 - I2
                    J3 = (2.0 / 27.0) * I1 * I1 * I1 - (1.0 / 3.0) * I1 * I2 + I3
                    if J2 < 1e-30:
                        break
                    sqJ2 = math.sqrt(J2)
                    arg_sr = -(J3 * math.sqrt(27.0)) / (2.0 * sqJ2 * sqJ2 * sqJ2)
                    if arg_sr < -1.0:
                        arg_sr = -1.0
                    if arg_sr > 1.0:
                        arg_sr = 1.0
                    Sr = arg_sr
                    I1_star = I1 + sigma_t
                    F1 = alpha * I1_star ** n_desai - gamma * I1_star * I1_star
                    F2_base = math.exp(beta1 * I1_star) - beta * Sr
                    if F2_base <= 0.0:
                        F_vp = J2
                    else:
                        F_vp = J2 + F1 * F2_base ** m_desai
                    if F_vp <= 0.0:
                        break
                    lam = mu1 * (F_vp ** N1)
                    p = I1 / 3.0
                    dev1 = s1 - p
                    rate_desai = lam * abs(dev1)
                    if rate_desai < 1e-30:
                        break
                    dt_sub = min(dt_rem, _MAX_DEPS / rate_desai)
                    if dt_sub < 1e-6:
                        dt_sub = 1e-6
                    eps_vp += rate_desai * dt_sub
                    zeta_desai += rate_desai * dt_sub
                    base = (a1 / alpha0) ** (1.0 / eta_vp) + zeta_desai
                    if base > 0.0:
                        alpha = a1 / (base ** eta_vp)
                        if alpha < 1e-30:
                            alpha = 1e-30
                    else:
                        alpha = alpha0
                    dt_rem -= dt_sub

        out_total[i] = (e_el + eps_d + eps_k + eps_vp) * 100.0
        out_dc[i] = eps_d * 100.0
        out_kv[i] = eps_k * 100.0
        out_de[i] = eps_vp * 100.0

    return out_total, out_dc, out_kv, out_de


@numba.njit(cache=True)
def integrate_md_fast(t_s, sigma_diff_Pa, A_pa_s, n, Q, T,
                      K0, c, m_md, alpha_w, beta_w, delta):
    """
    Numba-JIT MD integration with c as free parameter.
    Elastic params fixed at compile-time constants.
    """
    exp_QRT = math.exp(-Q / (_R * T))
    exp_CT = math.exp(c * T)
    n_pts = len(t_s)
    eps_c = 0.0
    zeta = 0.0
    eps = np.zeros(n_pts)
    log10_e = 1.0 / math.log(10.0)

    for i in range(n_pts):
        sigma = sigma_diff_Pa[i]
        dt = (t_s[i] - t_s[i - 1]) if i > 0 else 0.0
        e_el = sigma / _E_EL

        if dt > 0.0:
            rate_ss = A_pa_s * (abs(sigma) ** n) * exp_QRT
            if rate_ss < 1e-60:
                eps[i] = (e_el + eps_c) * 100.0
                continue

            dt_rem = dt
            n_sub = 0
            while dt_rem > 1e-10 and n_sub < _MAX_SUBSTEPS:
                n_sub += 1
                som = abs(sigma) / _MU_SH if _MU_SH > 0.0 else 0.0
                if som < 1e-60:
                    som = 1e-60
                eps_t_star = K0 * exp_CT * (som ** m_md)
                if eps_t_star < 1e-50:
                    eps_c += rate_ss * dt_rem
                    break
                log10_som = math.log(max(som, 1e-30)) * log10_e
                Delta = alpha_w + beta_w * log10_som
                ratio = zeta / eps_t_star
                if ratio <= 1.0:
                    F = math.exp(Delta * (1.0 - ratio) * (1.0 - ratio))
                else:
                    F = math.exp(-delta * (1.0 - ratio) * (1.0 - ratio))

                zeta_rate = abs(F - 1.0) * rate_ss
                if zeta_rate > 1e-30:
                    dt_sub = min(dt_rem, _MAX_DZETA / zeta_rate)
                else:
                    dt_sub = dt_rem
                if dt_sub < 10.0:
                    dt_sub = 10.0
                if dt_sub > dt_rem:
                    dt_sub = dt_rem

                eps_c += F * rate_ss * dt_sub
                zeta += (F - 1.0) * rate_ss * dt_sub
                if zeta < 0.0:
                    zeta = 0.0
                dt_rem -= dt_sub

        eps[i] = (e_el + eps_c) * 100.0

    return eps


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                     PHASE 2 RUNNER (per model)                             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

def run_phase2(model, A_mpa_yr, N_disloc, QR_disloc,
               t_s, sigma_diff_fine, t_mod_h,
               t_resamp, e_resamp, weights, stage_ids, stages):
    Q_val = QR_disloc * R_GAS
    A_pa_s = A_mpa_yr * (1e-6) ** N_disloc / _SEC_PER_YR

    _call_count = [0]
    _best = [np.inf]
    _best_params_opt = [None]

    def objective(x):
        try:
            if model == "sic":
                (log_eta, log_E1, log_mu1, N1, log_a1, eta_vp, log_alpha0,
                 n_desai, log_beta1, beta, m_desai, log_gamma) = x
                eta = 10.0 ** log_eta
                E1 = 10.0 ** log_E1
                mu1 = 10.0 ** log_mu1
                a1 = 10.0 ** log_a1
                alpha0 = 10.0 ** log_alpha0
                beta1 = 10.0 ** log_beta1
                gamma = 10.0 ** log_gamma
                sigma_t = 5.0

                e_mod = integrate_sic_free(
                    t_s, sigma_diff_fine,
                    A_pa_s, N_disloc, Q_val, T_KELVIN,
                    eta, E1, mu1, N1, a1, eta_vp, alpha0,
                    n_desai, beta1, beta, m_desai, gamma, sigma_t
                )
            else:
                log_K0, c_md, m_md, alpha_w, beta_w, delta = x
                K0 = 10.0 ** log_K0

                e_mod = integrate_md_fast(
                    t_s, sigma_diff_fine,
                    A_pa_s, N_disloc, Q_val, T_KELVIN,
                    K0, c_md, m_md, alpha_w, beta_w, delta
                )
        except Exception:
            return 1e6

        e_mod_resamp = np.interp(t_resamp, t_mod_h, e_mod)
        eps_floor = 0.005
        abs_lab = np.maximum(np.abs(e_resamp), eps_floor)
        ape = np.abs(e_mod_resamp - e_resamp) / abs_lab
        wmape = np.sum(weights * ape) / np.sum(weights) * 100.0

        if not np.isfinite(wmape):
            return 1e6

        _call_count[0] += 1

        if wmape < _best[0]:
            _best[0] = wmape
            _best_params_opt[0] = x.copy()
            if model == "sic":
                print(f"  [{_call_count[0]:5d}] NEW BEST WMAPE={wmape:.2f}%  "
                      f"eta={10**log_eta:.2e} E1={10**log_E1:.2e} "
                      f"mu1={10**log_mu1:.2e} N1={N1:.2f} "
                      f"a1={10**log_a1:.2e} eta_vp={eta_vp:.2f} "
                      f"a0={10**log_alpha0:.3e} n_d={n_desai:.2f} "
                      f"b1={10**log_beta1:.4f} beta={beta:.3f} "
                      f"m_d={m_desai:.2f} gam={10**log_gamma:.4f}", flush=True)
            else:
                print(f"  [{_call_count[0]:5d}] NEW BEST WMAPE={wmape:.2f}%  "
                      f"K0={10**log_K0:.4f} c={c_md:.5f} m={m_md:.3f} "
                      f"aw={alpha_w:.2f} bw={beta_w:.2f} d={delta:.4f}",
                      flush=True)
        elif _call_count[0] % 500 == 0:
            print(f"  [{_call_count[0]:5d}] best so far WMAPE={_best[0]:.2f}%",
                  flush=True)

        return wmape

    # ── Search bounds ─────────────────────────────────────────────────────────

    if model == "sic":
        # 12 parameters: Kelvin (2) + Desai calibrated (5) + Desai yield surface (5)
        BOUNDS = [
            (8.0, 15.0),       # log10(eta [Pa·s])
            (7.0, 11.0),       # log10(E1 [Pa])
            (-15.0, -8.0),     # log10(mu1 [1/s])
            (1.5, 2.5),        # N1 [-]
            (-7.0, -1.0),      # log10(a1)
            (0.3, 2.0),        # eta_vp
            (-4.0, -1.0),      # log10(alpha0)
            (1.0, 6.0),        # n_desai (power)       — was fixed at 3.0
            (-4.0, 0.0),       # log10(beta1)          — was fixed at 0.0048
            (0.5, 1.0),        # beta [-]              — was fixed at 0.995
            (-2.0, 0.0),       # m_desai [-]           — was fixed at -0.5
            (-4.0, 0.0),       # log10(gamma)          — was fixed at 0.095
        ]
        _SEED_ROW = np.array([
            12.0, 9.5,                # log10(eta), log10(E1)
            -15.0, 2.0,               # log10(mu1), N1
            -3.0, 1.2,                # log10(a1), eta_vp
            -2.3,                      # log10(alpha0)
            3.0,                       # n_desai
            np.log10(0.0048),          # log10(beta1)
            0.995,                     # beta
            -0.5,                      # m_desai
            np.log10(0.095),           # log10(gamma)
        ])
        _n_params = 12
        _popsize = 15
        _maxiter = 200
    else:
        # 6 parameters: all MD transient params including c
        BOUNDS = [
            (-6.0, 8.0),      # log10(K0)
            (0.001, 0.05),    # c [1/K]               — was fixed at 0.00902
            (0.1, 6.0),       # m_md
            (-10.0, 200.0),   # alpha_w
            (-30.0, 60.0),    # beta_w
            (0.001, 300.0),   # delta
        ]
        _SEED_ROW = np.array([
            3.39,              # log10(K0)
            0.00902,           # c
            2.47,              # m
            93.3,              # alpha_w
            30.0,              # beta_w
            99.4,              # delta
        ])
        _n_params = 6
        _popsize = 20
        _maxiter = 300

    _n_pop = _popsize * _n_params
    rng = np.random.default_rng(42)
    _init_pop = rng.uniform(0, 1, size=(_n_pop, _n_params))
    for j, (lo, hi) in enumerate(BOUNDS):
        _init_pop[:, j] = lo + _init_pop[:, j] * (hi - lo)
    _init_pop[0] = _SEED_ROW

    print(f"  {_n_params} parameters, popsize={_popsize}, maxiter={_maxiter}")
    print(f"  Dislocation (fixed from Phase 1): A={A_mpa_yr:.2f} MPa^-n/yr, "
          f"n={N_disloc:.4f}, Q/R={QR_disloc:.0f} K")

    t0 = _time.time()
    result = differential_evolution(
        objective,
        bounds=BOUNDS,
        popsize=_popsize,
        maxiter=_maxiter,
        tol=1e-6,
        mutation=(0.5, 1.5),
        recombination=0.8,
        seed=42,
        init=_init_pop,
        polish=True,
        disp=False,
    )
    elapsed = _time.time() - t0

    return result.x, result.fun, elapsed, _call_count[0]


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          MAIN EXECUTION                                    ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

def main():
    print("=" * 70)
    print("OPTIMIZER — Desai yield surface + MD c free")
    print("  Elastic fixed, dislocation creep unconstrained")
    print("  Only constraint: A_SIC = 1.5 x A_MD")
    print("=" * 70)
    print("\nLoading new lab data ...")

    time_h_raw, sigma1_raw, sigma3_raw, epsilon1_raw, epsilon3_raw = read_excel_data(XLSX_PATH)
    sigma_diff_raw = sigma1_raw - sigma3_raw

    t_s, sigma1_fine, sigma3_fine, idx_start = build_stress_schedule(
        time_h_raw, sigma1_raw, sigma3_raw, DT_HOURS
    )
    sigma_diff_fine = sigma1_fine - sigma3_fine

    t0_h = time_h_raw[idx_start]
    t_lab = time_h_raw[idx_start:] - t0_h
    e_lab_ax = epsilon1_raw[idx_start:] - epsilon1_raw[idx_start]
    s_lab = sigma_diff_raw[idx_start:]
    t_mod_h = t_s / HOUR

    stages = detect_cyclic_stages(t_lab, s_lab)

    print(f"  {len(time_h_raw)} raw pts, {len(stages)} stages detected")
    for k, stg in enumerate(stages):
        dur = stg['t_end'] - stg['t_start']
        print(f"    Stage {k}: {stg['t_start']:.0f}-{stg['t_end']:.0f} h, "
              f"sigma={stg['sigma_mpa']:.1f} MPa, type={stg['type']}, dur={dur:.0f} h")

    # Warm up numba JIT
    print("  Warming up numba JIT ...", end="", flush=True)
    _t_warm = np.linspace(0, 3600.0, 10)
    _s_warm = np.full(10, 10e6)
    _ = integrate_sic_free(_t_warm, _s_warm, 1e-40, 4.0, 50000.0, 294.15,
                           1e12, 5e8, 1e-15, 2.0, 1e-3, 1.2, 0.005,
                           3.0, 0.0048, 0.995, -0.5, 0.095, 5.0)
    _ = integrate_sic_free_components(_t_warm, _s_warm, 1e-40, 4.0, 50000.0, 294.15,
                                      1e12, 5e8, 1e-15, 2.0, 1e-3, 1.2, 0.005,
                                      3.0, 0.0048, 0.995, -0.5, 0.095, 5.0)
    _ = integrate_md_fast(_t_warm, _s_warm, 1e-40, 4.0, 50000.0, 294.15,
                          1.0, 0.009, 2.0, 5.0, -2.0, 1.0)
    print(" done", flush=True)

    t_resamp, e_resamp, s_resamp, weights, stage_ids = \
        resample_and_label_cyclic(t_lab, e_lab_ax, s_lab, stages)
    print(f"  Resampled to {len(t_resamp)} points")

    # ══════════════════════════════════════════════════════════════════════════
    # PHASE 1: DISLOCATION CREEP — literature-bounded
    # ══════════════════════════════════════════════════════════════════════════

    if SKIP_PHASE1:
        print("\n[PHASE 1] SKIPPED — using preset dislocation parameters")
        A_SIC_MPA_YR = 60.0;  N_SIC = 4.5;  QR_SIC = 6252.0
        A_MD_MPA_YR  = 40.0;  N_MD  = 4.5;  QR_MD  = 6252.0
    else:
        print("\n" + "=" * 70)
        print("PHASE 1: Dislocation creep from steady-state slopes (UNCONSTRAINED)")
        print("=" * 70)
        print(f"  Ranges: A in [{A_BOUNDS[0]}, {A_BOUNDS[1]}] MPa^-n/yr, "
              f"n in [{N_BOUNDS[0]}, {N_BOUNDS[1]}], "
              f"Q/R in [{QR_GRID[0]:.0f}, {QR_GRID[-1]:.0f}] K")

        stress_mpa, rates = extract_steady_state_rates(t_lab, e_lab_ax, s_lab, stages)

        if len(stress_mpa) < 2:
            print("[ERROR] Need at least 2 loading stages for Phase 1 fit")
            sys.exit(1)

        print(f"\n  Extracted steady-state rates from {len(stress_mpa)} loading stages:")
        for sm, r in zip(stress_mpa, rates):
            print(f"    sigma_diff = {sm:.1f} MPa -> rate = {r:.3e} /s "
                  f"({r * 100 * 3600:.4e} %/h)")

        best_params, best_resid = fit_dislocation_params(
            stress_mpa, rates, T_KELVIN, QR_GRID
        )

        if best_params is not None:
            A_FIT, N_FIT, QR_FIT = best_params
            print(f"\n  Best fit:")
            print(f"    A   = {A_FIT:.2f} MPa^-n/yr")
            print(f"    n   = {N_FIT:.4f}")
            print(f"    Q/R = {QR_FIT:.0f} K")
            print(f"    log-space residual = {best_resid:.4f}")
        else:
            print("\n[WARNING] Could not fit — using defaults")
            A_FIT = 40.0; N_FIT = 4.5; QR_FIT = 6252.0

        A_SIC_MPA_YR = 1.5 * A_FIT;  N_SIC = N_FIT;  QR_SIC = QR_FIT
        A_MD_MPA_YR  = A_FIT;         N_MD  = N_FIT;  QR_MD  = QR_FIT

        print(f"\n  Scalar fit:      A={A_FIT:.2f}, n={N_FIT:.4f}, Q/R={QR_FIT:.0f}")
        print(f"  SIC dislocation: A={A_SIC_MPA_YR:.2f} (x1.5 for 2/3 flow rule)")
        print(f"  MD  dislocation: A={A_MD_MPA_YR:.2f}")

    # ══════════════════════════════════════════════════════════════════════════
    # PHASE 2a: SIC
    # ══════════════════════════════════════════════════════════════════════════

    print("\n" + "=" * 70)
    print("PHASE 2a: SIC transient optimization (12 params)")
    print("=" * 70)

    x_sic, wmape_sic, elapsed_sic, n_evals_sic = run_phase2(
        "sic", A_SIC_MPA_YR, N_SIC, QR_SIC,
        t_s, sigma_diff_fine, t_mod_h,
        t_resamp, e_resamp, weights, stage_ids, stages,
    )

    (log_eta, log_E1, log_mu1, N1, log_a1, eta_vp, log_alpha0,
     n_desai, log_beta1, beta, m_desai, log_gamma) = x_sic
    eta_sic = 10.0 ** log_eta
    E1_sic = 10.0 ** log_E1
    mu1_sic = 10.0 ** log_mu1
    a1_sic = 10.0 ** log_a1
    alpha0_sic = 10.0 ** log_alpha0
    beta1_sic = 10.0 ** log_beta1
    gamma_sic = 10.0 ** log_gamma

    print(f"\n{'='*80}")
    print(f"SIC OPTIMUM (elapsed {elapsed_sic/60:.1f} min, {n_evals_sic} evals)")
    print(f"{'='*80}")
    print(f"  A          = {A_SIC_MPA_YR:.2f} MPa^-n/yr")
    print(f"  n          = {N_SIC:.4f}")
    print(f"  Q/R        = {QR_SIC:.0f} K")
    tau_h = eta_sic / E1_sic / 3600
    kelvin_eq = 21e6 / E1_sic * 100
    print(f"  ETA_KELVIN = {eta_sic:.4e} Pa s")
    print(f"  E1_KELVIN  = {E1_sic:.4e} Pa")
    print(f"  Kelvin tau = {tau_h:.1f} h, eq@21MPa = {kelvin_eq:.2f}%")
    print(f"  MU1_DESAI  = {mu1_sic:.4e} 1/s")
    print(f"  N1_DESAI   = {N1:.4f}")
    print(f"  A1_DESAI   = {a1_sic:.4e}")
    print(f"  ETA_DESAI  = {eta_vp:.4f}")
    print(f"  ALPHA0     = {alpha0_sic:.4e}")
    print(f"  N_DESAI    = {n_desai:.4f}  (was fixed at 3.0)")
    print(f"  BETA1      = {beta1_sic:.6f}  (was fixed at 0.0048)")
    print(f"  BETA       = {beta:.4f}  (was fixed at 0.995)")
    print(f"  M_DESAI    = {m_desai:.4f}  (was fixed at -0.5)")
    print(f"  GAMMA      = {gamma_sic:.6f}  (was fixed at 0.095)")
    print(f"  WMAPE      = {wmape_sic:.2f}%")

    # ══════════════════════════════════════════════════════════════════════════
    # PHASE 2b: MD
    # ══════════════════════════════════════════════════════════════════════════

    print("\n" + "=" * 70)
    print("PHASE 2b: MD transient optimization (6 params)")
    print("=" * 70)

    x_md, wmape_md, elapsed_md, n_evals_md = run_phase2(
        "md", A_MD_MPA_YR, N_MD, QR_MD,
        t_s, sigma_diff_fine, t_mod_h,
        t_resamp, e_resamp, weights, stage_ids, stages,
    )

    log_K0, c_md, m_md, alpha_w, beta_w, delta = x_md
    K0_md = 10.0 ** log_K0

    print(f"\n{'='*80}")
    print(f"MD OPTIMUM (elapsed {elapsed_md/60:.1f} min, {n_evals_md} evals)")
    print(f"{'='*80}")
    print(f"  A          = {A_MD_MPA_YR:.2f} MPa^-n/yr")
    print(f"  n          = {N_MD:.4f}")
    print(f"  Q/R        = {QR_MD:.0f} K")
    print(f"  K0         = {K0_md:.6f}")
    print(f"  c          = {c_md:.6f}  (was fixed at 0.00902)")
    print(f"  m          = {m_md:.4f}")
    print(f"  alpha_w    = {alpha_w:.4f}")
    print(f"  beta_w     = {beta_w:.4f}")
    print(f"  delta      = {delta:.6f}")
    print(f"  WMAPE      = {wmape_md:.2f}%")

    # ══════════════════════════════════════════════════════════════════════════
    # PER-STAGE METRICS
    # ══════════════════════════════════════════════════════════════════════════

    Q_SIC_val = QR_SIC * R_GAS
    Q_MD_val = QR_MD * R_GAS
    A_SIC_PA_S = A_SIC_MPA_YR * (1e-6) ** N_SIC / _SEC_PER_YR
    A_MD_PA_S = A_MD_MPA_YR * (1e-6) ** N_MD / _SEC_PER_YR

    for label, model_name in [("SIC", "sic"), ("MD", "md")]:
        print(f"\n{'─'*80}")
        print(f"Per-stage metrics — {label}:")
        print(f"{'─'*80}")

        if model_name == "sic":
            e_mod = integrate_sic_free(
                t_s, sigma_diff_fine,
                A_SIC_PA_S, N_SIC, Q_SIC_val, T_KELVIN,
                eta_sic, E1_sic, mu1_sic, N1, a1_sic, eta_vp, alpha0_sic,
                n_desai, beta1_sic, beta, m_desai, gamma_sic, 5.0
            )
        else:
            e_mod = integrate_md_fast(
                t_s, sigma_diff_fine,
                A_MD_PA_S, N_MD, Q_MD_val, T_KELVIN,
                K0_md, c_md, m_md, alpha_w, beta_w, delta
            )

        e_mod_resamp = np.interp(t_resamp, t_mod_h, e_mod)
        residuals = e_mod_resamp - e_resamp

        for k, stg in enumerate(stages):
            mask_stg = stage_ids == k
            if not np.any(mask_stg):
                continue
            res_stg = residuals[mask_stg]
            e_lab_stg = e_resamp[mask_stg]
            rmse = np.sqrt(np.mean(res_stg ** 2))
            eps_floor = 0.005
            abs_lab = np.maximum(np.abs(e_lab_stg), eps_floor)
            mape = np.mean(np.abs(res_stg) / abs_lab) * 100.0
            print(f"  Stage {k} ({stg['type']:>9}, sigma={stg['sigma_mpa']:.1f} MPa): "
                  f"RMSE={rmse:.4f}%  MAPE={mape:.1f}%")

        e_mod_full = np.interp(t_lab, t_mod_h, e_mod)
        rmse_total = np.sqrt(np.mean((e_mod_full - e_lab_ax) ** 2))
        final_err = float(e_mod_full[-1] - e_lab_ax[-1])
        print(f"\n  Overall RMSE = {rmse_total:.4f}%")
        print(f"  Final strain error = {final_err:+.4f}%")

    # ══════════════════════════════════════════════════════════════════════════
    # COPY-PASTE PARAMETERS
    # ══════════════════════════════════════════════════════════════════════════

    print(f"\n{'─'*80}")
    print("Copy-paste parameters for simulation scripts:")
    print(f"{'─'*80}")

    print(f"\n  -- SIC --")
    print(f"  A_SIC         = {A_SIC_MPA_YR:.2f}")
    print(f"  N_SIC         = {N_SIC:.4f}")
    print(f"  Q_OVER_R      = {QR_SIC:.0f}")
    print(f"  ETA_KELVIN    = {eta_sic:.4e}")
    print(f"  E1_KELVIN     = {E1_sic:.4e}")
    print(f"  MU1_DESAI     = {mu1_sic:.4e}")
    print(f"  N1_DESAI      = {N1:.4f}")
    print(f"  A1_DESAI      = {a1_sic:.4e}")
    print(f"  ETA_DESAI     = {eta_vp:.4f}")
    print(f"  ALPHA0_DESAI  = {alpha0_sic:.4e}")
    print(f"  N_DESAI_POWER = {n_desai:.4f}")
    print(f"  BETA1_DESAI   = {beta1_sic:.6f}")
    print(f"  BETA_DESAI    = {beta:.4f}")
    print(f"  M_DESAI       = {m_desai:.4f}")
    print(f"  GAMMA_DESAI   = {gamma_sic:.6f}")

    print(f"\n  -- MD --")
    print(f"  A_MD          = {A_MD_MPA_YR:.2f}")
    print(f"  N_MD          = {N_MD:.4f}")
    print(f"  Q_OVER_R      = {QR_MD:.0f}")
    print(f"  K0_MD         = {K0_md:.6f}")
    print(f"  C_MD          = {c_md:.6f}")
    print(f"  M_MD          = {m_md:.4f}")
    print(f"  ALPHA_W_MD    = {alpha_w:.4f}")
    print(f"  BETA_W_MD     = {beta_w:.4f}")
    print(f"  DELTA_MD      = {delta:.6f}")

    # ══════════════════════════════════════════════════════════════════════════
    # AUTO-SAVE: full integration → JSON
    # ══════════════════════════════════════════════════════════════════════════

    print(f"\n{'='*70}")
    print("AUTO-SAVE: Running full integration with optimized parameters ...")
    print(f"{'='*70}")

    import calibrate_newdata as cn
    import json

    time_h_raw_full, sigma1_raw_full, sigma3_raw_full, eps1_raw, eps3_raw = \
        cn.read_excel_data(XLSX_PATH)
    t_s_full, sig1_full, sig3_full, idx0 = cn.build_stress_schedule(
        time_h_raw_full, sigma1_raw_full, sigma3_raw_full, cn.DT_HOURS
    )

    t0_h_full = time_h_raw_full[idx0]
    lab_time = time_h_raw_full[idx0:] - t0_h_full
    lab_eps1 = eps1_raw[idx0:] - eps1_raw[idx0]
    lab_eps3 = eps3_raw[idx0:] - eps3_raw[idx0]
    lab_sdiff = (sigma1_raw_full - sigma3_raw_full)[idx0:]
    lab_sig1 = sigma1_raw_full[idx0:]
    lab_sig3 = sigma3_raw_full[idx0:]

    t_hours_full = t_s_full / HOUR
    out_dir = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out_dir, exist_ok=True)

    # -- SIC JSON --
    # Must use integrate_sic_free (not cn.integrate_safeincave) because the
    # Desai yield surface params are calibrated here, not the hardcoded ones
    # in calibrate_newdata.py.
    A_sic_pa_s_full = A_SIC_MPA_YR * (1e-6) ** N_SIC / _SEC_PER_YR
    sigma_diff_full = sig1_full - sig3_full
    eps_ax_sic_pct, eps_dc_sic_pct, eps_kv_sic_pct, eps_de_sic_pct = \
        integrate_sic_free_components(
            t_s_full, sigma_diff_full,
            A_sic_pa_s_full, N_SIC, QR_SIC * R_GAS, T_KELVIN,
            eta_sic, E1_sic, mu1_sic, N1, a1_sic, eta_vp, alpha0_sic,
            n_desai, beta1_sic, beta, m_desai, gamma_sic, 5.0
        )
    print(f"  SIC final axial strain:  {eps_ax_sic_pct[-1]:.4f}%")

    result_sic = {
        "test_id": "newdata_cyclic_free",
        "temperature_C": T_CELSIUS,
        "sigma3_MPa_mean": float(np.mean(sigma3_raw_full[idx0:])),
        "lab": {
            "time_hours": lab_time.tolist(),
            "strain_axial_pct": lab_eps1.tolist(),
            "strain_radial_pct": lab_eps3.tolist(),
            "sigma_diff_MPa": lab_sdiff.tolist(),
            "sigma1_MPa": lab_sig1.tolist(),
            "sigma3_MPa": lab_sig3.tolist(),
        },
        "params": {"Q_over_R": QR_SIC, "E_elastic": E_ELASTIC, "nu_elastic": NU_ELASTIC},
        "safeincave": {
            "time_hours": t_hours_full.tolist(),
            "strain_axial_pct": eps_ax_sic_pct.tolist(),
            "strain_radial_pct": [0.0] * len(t_s_full),
            "strain_disloc_pct": eps_dc_sic_pct.tolist(),
            "strain_kelvin_pct": eps_kv_sic_pct.tolist(),
            "strain_desai_pct": eps_de_sic_pct.tolist(),
            "params": {
                "A_mpa_yr": A_SIC_MPA_YR, "n": N_SIC,
                "eta_kelvin": eta_sic, "E1_kelvin": E1_sic,
                "mu1_desai": mu1_sic, "N1_desai": N1,
                "a1_desai": a1_sic, "eta_desai": eta_vp,
                "alpha0_desai": alpha0_sic,
                "n_desai_power": n_desai,
                "beta1_desai": beta1_sic,
                "beta_desai": beta,
                "m_desai": m_desai,
                "gamma_desai": gamma_sic,
            },
        },
    }
    json_path_sic = os.path.join(out_dir, "calibration_newdata_free_sic.json")
    with open(json_path_sic, "w") as f:
        json.dump(result_sic, f, indent=2)
    print(f"  [SAVED] {json_path_sic}")

    # -- MD JSON --
    A_md_pa_s_full = A_MD_MPA_YR * (1e-6) ** N_MD / _SEC_PER_YR
    eps_ax_md, eps_ss_md, eps_tr_md, eps_rad_md = cn.integrate_munsondawson(
        t_s_full, sig1_full, sig3_full,
        A=A_md_pa_s_full, n=N_MD, Q=QR_MD * R_GAS, T=T_KELVIN,
        K0=K0_md, c=c_md, m_md=m_md,
        alpha_w=alpha_w, beta_w=beta_w, delta=delta,
    )
    print(f"  MD  final axial strain:  {eps_ax_md[-1] * 100:.4f}%")

    result_md = {
        "test_id": "newdata_cyclic_free",
        "temperature_C": T_CELSIUS,
        "sigma3_MPa_mean": float(np.mean(sigma3_raw_full[idx0:])),
        "lab": {
            "time_hours": lab_time.tolist(),
            "strain_axial_pct": lab_eps1.tolist(),
            "strain_radial_pct": lab_eps3.tolist(),
            "sigma_diff_MPa": lab_sdiff.tolist(),
            "sigma1_MPa": lab_sig1.tolist(),
            "sigma3_MPa": lab_sig3.tolist(),
        },
        "params": {"Q_over_R": QR_MD, "E_elastic": E_ELASTIC, "nu_elastic": NU_ELASTIC},
        "munsondawson": {
            "time_hours": t_hours_full.tolist(),
            "strain_axial_pct": (eps_ax_md * 100).tolist(),
            "strain_radial_pct": (eps_rad_md * 100).tolist(),
            "strain_steady_pct": (eps_ss_md * 100).tolist(),
            "strain_transient_pct": (eps_tr_md * 100).tolist(),
            "params": {
                "A_mpa_yr": A_MD_MPA_YR, "n": N_MD,
                "K0": K0_md, "c": c_md, "m": m_md,
                "alpha_w": alpha_w, "beta_w": beta_w, "delta": delta,
            },
        },
    }
    json_path_md = os.path.join(out_dir, "calibration_newdata_free_md.json")
    with open(json_path_md, "w") as f:
        json.dump(result_md, f, indent=2)
    print(f"  [SAVED] {json_path_md}")

    # ══════════════════════════════════════════════════════════════════════════
    # AUTO-PLOT
    # ══════════════════════════════════════════════════════════════════════════

    print("\nGenerating plots ...")
    fig_dir = os.path.join(os.path.dirname(__file__), "figures")
    os.makedirs(fig_dir, exist_ok=True)

    combined = dict(result_sic)
    combined["munsondawson"] = result_md["munsondawson"]

    from plot_newdata import plot_newdata
    plot_newdata(combined, fig_dir)
    plot_newdata(combined, fig_dir, thesis=True)

    print("\n[ALL DONE]")


if __name__ == "__main__":
    main()
