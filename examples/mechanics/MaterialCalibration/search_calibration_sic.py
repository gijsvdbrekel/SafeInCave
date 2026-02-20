"""
Global SafeInCave (SIC) parameter search using scipy.optimize.differential_evolution.

SIC model = elastic spring + Kelvin element (transient) + dislocation creep (steady).
Desai viscoplastic is kept disabled (MU1=0) for speed; it can be added later.

Key differences from the MD search:
- N_SIC and A_SIC are independent from the MD model parameters.
- The Kelvin element (E1, eta) provides the transient creep burst.
- The Kelvin ODE is integrated analytically per timestep for numerical stability.

Free parameters (4):
  n         : dislocation exponent           [3.0 – 7.0]
  log_A     : log10(A [MPa^-n / yr])         [0.0 – 5.5]
  log_eta   : log10(eta_Kelvin [Pa·s])       [10.0 – 15.0]
  log_E1    : log10(E1_Kelvin [Pa])          [7.5 – 10.3]

Objective: mean RMSE over 6 lab tests (artefact-corrected stress).

Typical run time: 10–20 min (popsize=12, maxiter=150).
"""

import os, sys, time
import numpy as np
from scipy.optimize import differential_evolution

sys.path.insert(0, os.path.dirname(__file__))
from run_calibration import (
    read_creep_csv, build_stress_schedule,
    correct_stress_artefacts,
    R_GAS, T_KELVIN, Q_DISLOC, E_ELASTIC,
    DATA_DIR, HOUR,
)

# ─── configuration ────────────────────────────────────────────────────────────
TESTS    = ["TCC1", "TCC2", "TCC6", "TCC7", "TCC11", "TCC12"]
WEIGHTS  = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
DT_HOURS = 1.0

_SEC_PER_YR = 365.25 * 24 * 3600
_EXP_QRT    = np.exp(-Q_DISLOC / (R_GAS * T_KELVIN))

# ─── SIC integration ──────────────────────────────────────────────────────────
def integrate_sic(t_s, sigma_Pa, n, A_pa_s, eta, E1):
    """
    Forward-Euler dislocation + exact-update Kelvin integration.

    Kelvin ODE: dε_k/dt = (σ - E1·ε_k) / η
    Exact solution for constant σ over dt:
        ε_k(t+dt) = σ/E1 + (ε_k(t) - σ/E1) * exp(-E1·dt/η)

    Returns axial strain [%] at each time step.
    """
    n_pts = len(t_s)
    eps_k = 0.0    # Kelvin strain
    eps_d = 0.0    # dislocation strain
    eps   = np.zeros(n_pts)

    for i in range(n_pts):
        sigma = sigma_Pa[i]
        dt    = (t_s[i] - t_s[i - 1]) if i > 0 else 0.0
        e_el  = sigma / E_ELASTIC

        if dt > 0:
            # Dislocation creep (Euler)
            rate_dc = A_pa_s * (abs(sigma) ** n) * _EXP_QRT
            eps_d  += rate_dc * dt

            # Kelvin element (exact analytical update)
            eps_k_eq  = sigma / E1 if E1 > 0 else 0.0
            decay     = np.exp(-E1 * dt / eta) if (eta > 0 and E1 > 0) else 0.0
            eps_k     = eps_k_eq + (eps_k - eps_k_eq) * decay

        eps[i] = (e_el + eps_d + eps_k) * 100.0

    return eps


# ─── load lab data once ───────────────────────────────────────────────────────
print("Loading lab data …")
_lab = {}
for tid in TESTS:
    fp = os.path.join(DATA_DIR, f"ZW_{tid}.csv")
    time_h, strain_pct, sigma_diff, sigma3, _ = read_creep_csv(fp)
    sigma_diff_corr = correct_stress_artefacts(time_h, sigma_diff)
    t_s, sigma_Pa, idx = build_stress_schedule(time_h, sigma_diff_corr, DT_HOURS)
    t_lab = time_h[idx:] - time_h[idx]
    e_lab = strain_pct[idx:] - strain_pct[idx]
    sigma3_Pa = float(sigma3[idx]) * 1e6
    _lab[tid] = dict(t_s=t_s, sigma_Pa=sigma_Pa,
                     t_mod_h=t_s / HOUR,
                     t_lab=t_lab, e_lab=e_lab,
                     sigma3_Pa=sigma3_Pa)
    print(f"  {tid}: {len(time_h)} pts, max σ_diff = {sigma_diff_corr[idx:].max():.1f} MPa, "
          f"σ3 = {sigma3_Pa/1e6:.1f} MPa")


# ─── objective ────────────────────────────────────────────────────────────────
_call_count  = [0]
_best        = [np.inf]
_best_params = [None]

def objective(x):
    n, log_A, log_eta, log_E1 = x
    A_mpa_yr = 10.0 ** log_A
    A_pa_s   = A_mpa_yr * (1e-6) ** n / _SEC_PER_YR
    eta      = 10.0 ** log_eta
    E1       = 10.0 ** log_E1

    total = 0.0
    w_sum = sum(WEIGHTS)
    for tid, w in zip(TESTS, WEIGHTS):
        d = _lab[tid]
        try:
            e_mod = integrate_sic(d['t_s'], d['sigma_Pa'], n, A_pa_s, eta, E1)
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
              f"n={n:.3f} A={A_mpa_yr:.2f} η={eta:.2e} E1={E1:.2e}", flush=True)
    elif _call_count[0] % 500 == 0:
        print(f"  [{_call_count[0]:5d}] best so far {_best[0]:.4f}%", flush=True)

    return obj


# ─── search bounds ────────────────────────────────────────────────────────────
# [n, log10(A [MPa^-n/yr]), log10(eta [Pa·s]), log10(E1 [Pa])]
BOUNDS = [
    (3.0,  7.0),    # n
    (0.0,  5.5),    # log10(A_mpa_yr)  [1 to 316000 MPa^-n/yr]
    (10.0, 15.0),   # log10(eta)       [1e10 to 1e15 Pa·s]
    (7.5,  10.3),   # log10(E1)        [30 MPa to 2 GPa]
]

print(f"\nStarting differential_evolution (SIC model) …")
print(f"  4 parameters: n, log(A), log(eta), log(E1)")
print(f"  popsize=12, maxiter=150, tol=1e-4")

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
    polish=True,
    disp=False,
)
elapsed = time.time() - t0

# ─── report ───────────────────────────────────────────────────────────────────
x                      = result.x
n, log_A, log_eta, log_E1 = x
A_mpa_yr = 10.0 ** log_A
eta      = 10.0 ** log_eta
E1       = 10.0 ** log_E1
A_pa_s   = A_mpa_yr * (1e-6) ** n / _SEC_PER_YR

# Kelvin equilibrium at reference stresses
kelvin_eq_21 = 21e6 / E1 * 100   # % at 21 MPa
kelvin_eq_28 = 28e6 / E1 * 100   # % at 28 MPa
tau_h        = eta / E1 / 3600   # time constant [hours]

print(f"\n{'='*80}")
print(f"SIC GLOBAL OPTIMUM  (elapsed {elapsed/60:.1f} min, {_call_count[0]} evals)")
print(f"{'='*80}")
print(f"  n_SIC      = {n:.4f}")
print(f"  A_SIC      = {A_mpa_yr:.2f} MPa^-n/yr")
print(f"  ETA_KELVIN = {eta:.3e} Pa·s")
print(f"  E1_KELVIN  = {E1:.3e} Pa")
print(f"  Kelvin equilibrium @ 21 MPa  = {kelvin_eq_21:.2f}%")
print(f"  Kelvin equilibrium @ 28 MPa  = {kelvin_eq_28:.2f}%")
print(f"  Kelvin time constant τ = {tau_h:.1f} h")
print(f"  Combined RMSE = {result.fun:.4f}%")
print()

# Per-test RMSE at DT=0.5h
print("Per-test RMSE (DT=0.5h, corrected stress):")
for tid in TESTS:
    fp = os.path.join(DATA_DIR, f"ZW_{tid}.csv")
    time_h, strain_pct, sigma_diff, _, _ = read_creep_csv(fp)
    sigma_diff_corr = correct_stress_artefacts(time_h, sigma_diff)
    t_s, sigma_Pa, idx = build_stress_schedule(time_h, sigma_diff_corr, 0.5)
    t_lab   = time_h[idx:] - time_h[idx]
    e_lab   = strain_pct[idx:] - strain_pct[idx]
    t_mod_h = t_s / HOUR

    e_new = integrate_sic(t_s, sigma_Pa, n, A_pa_s, eta, E1)
    r_new = float(np.sqrt(np.mean((np.interp(t_lab, t_mod_h, e_new) - e_lab)**2)))
    f_new = float(np.interp(t_lab[-1], t_mod_h, e_new) - e_lab[-1])
    print(f"  {tid}: RMSE={r_new:.4f}%  final_err={f_new:+.3f}%")

print(f"\nTo update run_calibration.py:")
print(f"  _A_SIC_MPA_YR = {A_mpa_yr:.2f}")
print(f"  N_SIC         = {n:.4f}")
print(f"  ETA_KELVIN    = {eta:.4e}")
print(f"  E1_KELVIN     = {E1:.4e}")
