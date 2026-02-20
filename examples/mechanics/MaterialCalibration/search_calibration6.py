"""
Global MD parameter search using scipy.optimize.differential_evolution.

Key changes vs search_calibration5.py
--------------------------------------
1. Artefact correction: brief dips from >15 MPa to <5 MPa for <24 h (followed
   by return to >15 MPa) are filled with the preceding plateau value, matching
   the real monotonically-decreasing protocol.

2. A is parameterised via A_norm = rate at σ_ref=21 MPa so that n and A do
   not double-count the overall amplitude.

3. Wide bounds on all 6 MD parameters (n, A_norm, K0, alpha_w, delta, m_md).

4. Objective = mean RMSE over 6 tests, equally weighted.

Typical run time: 10-20 min (popsize=10, maxiter=120, 6 tests, DT_HOURS=1.0).
"""

import os, sys, time
import numpy as np
from scipy.optimize import differential_evolution

sys.path.insert(0, os.path.dirname(__file__))
from run_calibration import (
    read_creep_csv, build_stress_schedule,
    R_GAS, T_KELVIN, Q_DISLOC, E_ELASTIC,
    MU_SHEAR, C_MD, BETA_W_MD,
    DATA_DIR, HOUR,
    _MD_MAX_DZETA,
)

# ─── configuration ────────────────────────────────────────────────────────────
TESTS    = ["TCC1", "TCC2", "TCC6", "TCC7", "TCC11", "TCC12"]
WEIGHTS  = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   # equal weight per test
DT_HOURS = 1.0   # coarser for speed; verify best params with 0.5 afterwards

# Reference stress for A normalisation  [Pa]
_SIGMA_REF = 21.0e6
_EXP_QRT   = np.exp(-Q_DISLOC / (R_GAS * T_KELVIN))   # Arrhenius factor
_EXP_CT    = np.exp(C_MD * T_KELVIN)
_SEC_PER_YR = 365.25 * 24 * 3600

# ─── artefact correction ──────────────────────────────────────────────────────
def correct_stress_artefacts(time_h, sigma_diff,
                              high_thresh=15.0,
                              dip_thresh=5.0,
                              max_dip_h=24.0):
    """
    Fill brief dips from high stress (<dip_thresh MPa for <max_dip_h hours)
    that are surrounded by high-stress plateaus (>high_thresh MPa).
    These dips are measurement artefacts; the real protocol only ever
    decreases stress.
    """
    sig = sigma_diff.copy()
    n   = len(time_h)
    i   = 0
    while i < n - 1:
        # Step 1: find start of a dip (transition from high to low)
        if sig[i] >= high_thresh and sig[i + 1] < dip_thresh:
            dip_start = i + 1
            plateau_val = sig[i]
            # Step 2: find end of dip
            j = dip_start
            while j < n and sig[j] < dip_thresh:
                j += 1
            # Step 3: check whether stress recovers to high and dip is short
            if j < n and sig[j] >= high_thresh:
                dip_duration = time_h[j] - time_h[dip_start]
                if dip_duration <= max_dip_h:
                    sig[dip_start:j] = plateau_val   # fill artefact
                    i = j
                    continue
        i += 1
    return sig


# ─── MD integration ───────────────────────────────────────────────────────────
def integrate_md(t_s, sigma_Pa, n, A_norm, K0, alpha_w, delta, m_md):
    """
    Forward-Euler Munson-Dawson integration.

    Parameters
    ----------
    A_norm : float
        Creep rate at sigma_ref=21 MPa [1/s], i.e.
        eps_dot = A_norm * (sigma / sigma_ref)^n * exp(-Q/RT at T=100°C)
        (the exp factor is absorbed into A_norm for the reference condition)
    """
    # recover A_pa_s: A_norm = A_pa_s * sigma_ref^n * exp(-Q/RT)
    A_pa_s = A_norm / ((_SIGMA_REF ** n) * _EXP_QRT)

    n_pts = len(t_s)
    eps_c = 0.0
    zeta  = 0.0
    eps   = np.zeros(n_pts)

    for i in range(n_pts):
        sigma = sigma_Pa[i]
        dt    = (t_s[i] - t_s[i - 1]) if i > 0 else 0.0
        e_el  = sigma / E_ELASTIC

        if dt > 0:
            rate_ss = A_pa_s * (abs(sigma) ** n) * _EXP_QRT
            if rate_ss < 1e-60:
                eps[i] = (e_el + eps_c) * 100.0
                continue

            dt_rem = dt
            while dt_rem > 1e-10:
                som = abs(sigma) / MU_SHEAR if MU_SHEAR > 0 else 0.0
                eps_t_star = K0 * _EXP_CT * (max(som, 1e-60) ** m_md)
                if eps_t_star < 1e-50:
                    eps_c += rate_ss * dt_rem
                    break
                Delta = alpha_w + BETA_W_MD * np.log10(max(som, 1e-30))
                ratio = zeta / eps_t_star
                if ratio <= 1.0:
                    F = np.exp(Delta * (1.0 - ratio) ** 2)
                else:
                    F = np.exp(-delta * (1.0 - ratio) ** 2)

                zeta_rate = abs(F - 1.0) * rate_ss
                if zeta_rate > 1e-30:
                    dt_sub = min(dt_rem, _MD_MAX_DZETA / zeta_rate)
                else:
                    dt_sub = dt_rem
                dt_sub = max(dt_sub, 10.0)
                dt_sub = min(dt_sub, dt_rem)

                eps_c += F * rate_ss * dt_sub
                zeta  += (F - 1.0) * rate_ss * dt_sub
                zeta   = max(zeta, 0.0)
                dt_rem -= dt_sub

        eps[i] = (e_el + eps_c) * 100.0

    return eps


# ─── load lab data once ───────────────────────────────────────────────────────
print("Loading lab data …")
_lab = {}
for tid in TESTS:
    fp = os.path.join(DATA_DIR, f"ZW_{tid}.csv")
    time_h, strain_pct, sigma_diff, sigma3, _ = read_creep_csv(fp)

    # artefact correction on σ_diff
    sigma_diff_corr = correct_stress_artefacts(time_h, sigma_diff)

    t_s, sigma_Pa, idx = build_stress_schedule(time_h, sigma_diff_corr, DT_HOURS)
    t_lab = time_h[idx:] - time_h[idx]
    e_lab = strain_pct[idx:] - strain_pct[idx]
    _lab[tid] = dict(t_s=t_s, sigma_Pa=sigma_Pa,
                     t_mod_h=t_s / HOUR,
                     t_lab=t_lab, e_lab=e_lab)
    print(f"  {tid}: {len(time_h)} pts, corrected σ_diff max = "
          f"{sigma_diff_corr[idx:].max():.1f} MPa, "
          f"lab_final = {e_lab[-1]:.2f}%")


# ─── objective ────────────────────────────────────────────────────────────────
_call_count = [0]
_best       = [np.inf]
_best_params = [None]

def objective(x):
    n, log_A_norm, log_K0, alpha_w, delta, m_md = x
    A_norm = 10.0 ** log_A_norm
    K0     = 10.0 ** log_K0

    total = 0.0
    w_sum = sum(WEIGHTS)
    for tid, w in zip(TESTS, WEIGHTS):
        d = _lab[tid]
        try:
            e_mod = integrate_md(d['t_s'], d['sigma_Pa'],
                                 n, A_norm, K0, alpha_w, delta, m_md)
        except Exception:
            return 1e6
        e_int = np.interp(d['t_lab'], d['t_mod_h'], e_mod)
        rmse  = np.sqrt(np.mean((e_int - d['e_lab']) ** 2))
        total += w * rmse

    obj = total / w_sum

    _call_count[0] += 1
    if obj < _best[0]:
        _best[0] = obj
        _best_params[0] = x.copy()
        if _call_count[0] % 200 == 1:
            print(f"  [{_call_count[0]:5d}] NEW BEST {obj:.4f}%  "
                  f"n={n:.3f} log_A={log_A_norm:.2f} K0={K0:.3f} "
                  f"aw={alpha_w:.1f} d={delta:.3f} m={m_md:.3f}", flush=True)
    elif _call_count[0] % 500 == 0:
        print(f"  [{_call_count[0]:5d}] best so far {_best[0]:.4f}%", flush=True)

    return obj


# ─── search bounds ────────────────────────────────────────────────────────────
# [n, log10(A_norm), log10(K0), alpha_w, delta, m_md]
BOUNDS = [
    (2.5,  8.0),    # n
    (-9.0, -4.0),   # log10(A_norm)  [1/s at 21 MPa]
    (-2.5,  1.5),   # log10(K0)
    (-50.0, 5.0),   # alpha_w
    (0.01,  3.0),   # delta
    (0.2,   4.0),   # m_md
]

print(f"\nStarting differential_evolution …")
print(f"  6 parameters, popsize=12, maxiter=150, tol=1e-4")
print(f"  Estimated evaluations: 12*6*150 ≈ 10800  (~{10800*len(TESTS)*0.015/60:.0f} min)")

t0 = time.time()
result = differential_evolution(
    objective,
    bounds=BOUNDS,
    popsize=12,
    maxiter=150,
    tol=1e-4,
    mutation=(0.5, 1.5),
    recombination=0.8,
    seed=42,
    polish=True,      # final local polish with L-BFGS-B
    disp=False,
)
elapsed = time.time() - t0

# ─── report ───────────────────────────────────────────────────────────────────
x    = result.x
n, log_A, log_K0, alpha_w, delta, m_md = x
A_norm = 10.0 ** log_A
K0     = 10.0 ** log_K0

# Convert A_norm to A_mpa_yr for reporting
A_pa_s = A_norm / ((_SIGMA_REF ** n) * _EXP_QRT)
A_mpa_yr = A_pa_s * (1e6 ** n) * _SEC_PER_YR

print(f"\n{'='*80}")
print(f"GLOBAL OPTIMUM  (elapsed {elapsed/60:.1f} min, {_call_count[0]} evals)")
print(f"{'='*80}")
print(f"  n          = {n:.4f}")
print(f"  A_norm     = {A_norm:.3e} /s  (rate at 21 MPa, 100°C)")
print(f"  A_mpa_yr   = {A_mpa_yr:.4g} MPa^-n/yr")
print(f"  K0         = {K0:.4f}")
print(f"  m_md       = {m_md:.4f}")
print(f"  alpha_w    = {alpha_w:.4f}")
print(f"  delta      = {delta:.4f}")
print(f"  Combined RMSE = {result.fun:.4f}%")
print()

# Per-test RMSE with 0.5h timestep for final evaluation
print("Per-test RMSE (DT=0.5h, corrected stress):")
from run_calibration import A_MD, N_MD, K0_MD, ALPHA_W_MD, DELTA_MD, M_MD

for tid in TESTS:
    fp = os.path.join(DATA_DIR, f"ZW_{tid}.csv")
    time_h, strain_pct, sigma_diff, _, _ = read_creep_csv(fp)
    sigma_diff_corr = correct_stress_artefacts(time_h, sigma_diff)
    t_s, sigma_Pa, idx = build_stress_schedule(time_h, sigma_diff_corr, 0.5)
    t_lab = time_h[idx:] - time_h[idx]
    e_lab = strain_pct[idx:] - strain_pct[idx]
    t_mod_h = t_s / HOUR

    e_new = integrate_md(t_s, sigma_Pa, n, A_norm, K0, alpha_w, delta, m_md)
    r_new = float(np.sqrt(np.mean((np.interp(t_lab, t_mod_h, e_new) - e_lab)**2)))
    f_new = float(np.interp(t_lab[-1], t_mod_h, e_new) - e_lab[-1])

    # old (current run_calibration params with uncorrected stress)
    t_s2, sigma_Pa2, idx2 = build_stress_schedule(time_h, sigma_diff, 0.5)
    t_lab2 = time_h[idx2:] - time_h[idx2]
    e_lab2 = strain_pct[idx2:] - strain_pct[idx2]
    A_norm_old = A_MD * (_SIGMA_REF ** N_MD) * _EXP_QRT   # reverse of A_pa_s formula
    # Actually: A_pa_s_old * sigma_ref^n_old * EXP_QRT = A_norm_old
    A_pa_s_old = A_MD  # already in Pa^-n/s
    A_norm_old2 = A_pa_s_old * (_SIGMA_REF ** N_MD) * _EXP_QRT

    from run_calibration import integrate_munsondawson
    e_old_tot, _, _ = integrate_munsondawson(t_s2, sigma_Pa2)
    e_old = e_old_tot * 100
    r_old = float(np.sqrt(np.mean((np.interp(t_lab2, t_s2/HOUR, e_old) - e_lab2)**2)))
    f_old = float(np.interp(t_lab2[-1], t_s2/HOUR, e_old) - e_lab2[-1])

    print(f"  {tid}: NEW RMSE={r_new:.4f}% fin={f_new:+.3f}%  |  "
          f"OLD RMSE={r_old:.4f}% fin={f_old:+.3f}%")
