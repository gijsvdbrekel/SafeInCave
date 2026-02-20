"""
Global MD parameter search — extended bounds on K0 and delta.

Key changes vs search_calibration6.py
--------------------------------------
1. log10(K0) upper bound extended from 1.5 to 3.5 (K0 up to ~3162).
2. delta lower bound reduced from 0.01 to 0.001.
3. Runs with popsize=15, maxiter=200 for more thorough search.
4. Warm-started with search6 best params in the initial population.
5. Output written to search_calibration7.log for later review.

Context
-------
search6 found K0=31.62 (at upper bound) and delta=0.010 (at lower bound),
suggesting the true optimum lies outside those bounds.  This script
extends the search space to find the true global minimum.
"""

import os, sys, time
import numpy as np
from scipy.optimize import differential_evolution

sys.path.insert(0, os.path.dirname(__file__))
from run_calibration import (
    read_creep_csv, build_stress_schedule,
    correct_stress_artefacts,
    R_GAS, T_KELVIN, Q_DISLOC, E_ELASTIC,
    MU_SHEAR, C_MD, BETA_W_MD,
    DATA_DIR, HOUR,
    _MD_MAX_DZETA,
)

# ─── configuration ────────────────────────────────────────────────────────────
TESTS    = ["TCC1", "TCC2", "TCC6", "TCC7", "TCC11", "TCC12"]
WEIGHTS  = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
DT_HOURS = 1.0

_SIGMA_REF  = 21.0e6
_EXP_QRT    = np.exp(-Q_DISLOC / (R_GAS * T_KELVIN))
_EXP_CT     = np.exp(C_MD * T_KELVIN)
_SEC_PER_YR = 365.25 * 24 * 3600

# ─── MD integration (same as search6) ─────────────────────────────────────────
def integrate_md(t_s, sigma_Pa, n, A_norm, K0, alpha_w, delta, m_md):
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


# ─── load lab data ─────────────────────────────────────────────────────────────
print("Loading lab data …")
_lab = {}
for tid in TESTS:
    fp = os.path.join(DATA_DIR, f"ZW_{tid}.csv")
    time_h, strain_pct, sigma_diff, sigma3, _ = read_creep_csv(fp)
    sigma_diff_corr = correct_stress_artefacts(time_h, sigma_diff)
    t_s, sigma_Pa, idx = build_stress_schedule(time_h, sigma_diff_corr, DT_HOURS)
    t_lab = time_h[idx:] - time_h[idx]
    e_lab = strain_pct[idx:] - strain_pct[idx]
    _lab[tid] = dict(t_s=t_s, sigma_Pa=sigma_Pa,
                     t_mod_h=t_s / HOUR,
                     t_lab=t_lab, e_lab=e_lab)
    print(f"  {tid}: {len(time_h)} pts, max σ_diff = {sigma_diff_corr[idx:].max():.1f} MPa")


# ─── objective ─────────────────────────────────────────────────────────────────
_call_count  = [0]
_best        = [np.inf]
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
        print(f"  [{_call_count[0]:5d}] NEW BEST {obj:.4f}%  "
              f"n={n:.3f} log_A={log_A_norm:.2f} K0={K0:.3f} "
              f"aw={alpha_w:.1f} d={delta:.4f} m={m_md:.3f}", flush=True)
    elif _call_count[0] % 500 == 0:
        print(f"  [{_call_count[0]:5d}] best so far {_best[0]:.4f}%", flush=True)

    return obj


# ─── extended search bounds ────────────────────────────────────────────────────
# [n, log10(A_norm), log10(K0), alpha_w, delta, m_md]
BOUNDS = [
    (2.0,   8.0),    # n  — wider lower bound
    (-9.0, -4.0),    # log10(A_norm) [1/s at 21 MPa]
    (-2.5,  3.5),    # log10(K0)  — extended upper bound (K0 up to 3162)
    (-50.0, 5.0),    # alpha_w
    (0.001, 3.0),    # delta  — extended lower bound
    (0.2,   4.0),    # m_md
]

# Seed the initial population with the search6 best params
# n=2.7681, log_A=log10(7.731e-8)≈-7.11, log_K0=log10(31.62)=1.50,
# aw=-13.632, delta=0.010, m=1.4841
_SEED_ROW = np.array([2.7681, -7.11, 1.50, -13.632, 0.010, 1.4841])

popsize = 15
n_params = len(BOUNDS)
n_pop    = popsize * n_params    # total population size = 90

rng = np.random.default_rng(42)
# Latin-hypercube initial population
pop = rng.uniform(0, 1, size=(n_pop, n_params))
for j, (lo, hi) in enumerate(BOUNDS):
    pop[:, j] = lo + pop[:, j] * (hi - lo)
# Replace first member with seed (search6 best)
pop[0] = _SEED_ROW

print(f"\nStarting differential_evolution (extended bounds) …")
print(f"  Bounds: log_K0 up to 3.5, delta down to 0.001")
print(f"  6 parameters, popsize={popsize}, maxiter=200, tol=1e-5")
print(f"  Population seeded with search6 best params")

t0 = time.time()
result = differential_evolution(
    objective,
    bounds=BOUNDS,
    popsize=popsize,
    maxiter=200,
    tol=1e-5,
    mutation=(0.5, 1.5),
    recombination=0.8,
    seed=42,
    init=pop,
    polish=True,
    disp=False,
)
elapsed = time.time() - t0

# ─── report ───────────────────────────────────────────────────────────────────
x    = result.x
n, log_A, log_K0, alpha_w, delta, m_md = x
A_norm = 10.0 ** log_A
K0     = 10.0 ** log_K0

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
print(f"  delta      = {delta:.6f}")
print(f"  Combined RMSE = {result.fun:.4f}%")
print()

# Per-test RMSE at DT=0.5h
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

    print(f"  {tid}: RMSE={r_new:.4f}%  final_err={f_new:+.3f}%")

print(f"\nTo update run_calibration.py:")
print(f"  _A_MPA_YR  = {A_mpa_yr:.1f}")
print(f"  N_DISLOC   = {n:.4f}")
print(f"  K0_MD      = {K0:.4f}")
print(f"  M_MD       = {m_md:.4f}")
print(f"  ALPHA_W_MD = {alpha_w:.4f}")
print(f"  DELTA_MD   = {delta:.6f}")
