"""
Weighted creep parameter optimizer — zone-based weighting for transition emphasis.

The standard flat-RMSE objective is dominated by the many data points in the
flat tail (low-stress stages) of TCC tests.  This script:

1. Detects stress stages from sigma_diff(t).
2. Resamples to N uniform points per stage (removes point-density bias).
3. Subtracts elastic strain to isolate viscous behaviour.
4. Labels each point as elastic / transition / steady-state zone.
5. Applies per-zone weights (heavy on transition where Kelvin/F-shape matters).
6. Normalizes residuals per-stage (prevents high-stress stage dominating).
7. Runs differential_evolution on the weighted objective.

Supports both SIC (4 params) and MD (6 params) models.

Usage:
    python search_calibration_weighted.py          # default: SIC model, TCC1
    MODEL="md" python search_calibration_weighted.py
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


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          CONFIGURATION                                     ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

MODEL = os.environ.get("MODEL", "sic").lower()  # "sic" or "md"
TESTS = os.environ.get("TESTS", "TCC1").split(",")

# Resampling
N_PER_STAGE = 30          # uniform points per detected stress stage

# Zone weights
W_ELASTIC    = 0.5        # first ~2 h after load change (E is fixed)
W_TRANSITION = 2.0        # where transient mechanism is visible
W_STEADY     = 2.0        # constrains A and n (steady-state rate)

# Per-stage normalization: if False, use absolute residuals (no division by strain range)
NORMALIZE_PER_STAGE = False

# Transition duration [hours from stage start]
T_ELASTIC    = 2.0        # elastic zone cutoff
T_TRANS_HIGH = 50.0       # transition zone end for high-stress stages (σ > 10 MPa)
T_TRANS_LOW  = 30.0       # transition zone end for low-stress stages

# Low-stress global multiplier
LOW_STRESS_THRESH = 5.0   # [MPa]
LOW_STRESS_MULT   = 0.2   # weight multiplier for stages with σ < threshold

# Integration time step
DT_HOURS = 1.0

# Precomputed constants
_SEC_PER_YR = 365.25 * 24 * 3600
_EXP_QRT    = np.exp(-Q_DISLOC / (R_GAS * T_KELVIN))
_SIGMA_REF  = 21.0e6
_EXP_CT     = np.exp(C_MD * T_KELVIN)


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                      STAGE DETECTION & RESAMPLING                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

def detect_stress_stages(time_h, sigma_diff, min_change_mpa=3.0, min_hold_h=5.0):
    """
    Detect sustained stress stages from sigma_diff(t).

    A new stage starts when sigma_diff changes by more than min_change_mpa
    and the new level is held for at least min_hold_h hours.

    Returns list of dicts: [{'t_start': h, 't_end': h, 'sigma_mpa': float}, ...]
    """
    stages = []
    current_sigma = sigma_diff[0]
    stage_start = time_h[0]

    for i in range(1, len(time_h)):
        if abs(sigma_diff[i] - current_sigma) > min_change_mpa:
            # Potential new stage — check if it holds
            new_sigma = sigma_diff[i]
            hold_end = time_h[i] + min_hold_h
            # Check that stress stays near new_sigma for min_hold_h
            held = True
            for j in range(i + 1, len(time_h)):
                if time_h[j] > hold_end:
                    break
                if abs(sigma_diff[j] - new_sigma) > min_change_mpa:
                    held = False
                    break
            else:
                # Reached end of data before hold_end — still counts if
                # remaining duration > 0.5 * min_hold_h
                if time_h[-1] - time_h[i] < 0.5 * min_hold_h:
                    held = False

            if held:
                # Close current stage
                stages.append({
                    't_start': stage_start,
                    't_end': time_h[i],
                    'sigma_mpa': float(current_sigma),
                })
                current_sigma = new_sigma
                stage_start = time_h[i]

    # Close final stage
    stages.append({
        't_start': stage_start,
        't_end': time_h[-1],
        'sigma_mpa': float(current_sigma),
    })

    return stages


def resample_and_label(time_h, strain_pct, sigma_diff, stages,
                       n_per_stage=N_PER_STAGE):
    """
    Resample lab data to uniform points per stage and assign zone labels/weights.

    Returns:
        t_resamp    — resampled time [hours]
        e_resamp    — resampled strain [%]
        sigma_resamp — resampled sigma_diff [MPa]
        weights     — per-point weights
        scales      — per-point normalization scale (strain range of stage)
        stage_ids   — per-point stage index
    """
    t_all, e_all, s_all, w_all, sc_all, sid_all = [], [], [], [], [], []

    for k, stg in enumerate(stages):
        t0, t1 = stg['t_start'], stg['t_end']
        sigma_mpa = stg['sigma_mpa']

        # Mask for this stage
        mask = (time_h >= t0) & (time_h <= t1)
        t_stg = time_h[mask]
        e_stg = strain_pct[mask]
        s_stg = sigma_diff[mask]

        if len(t_stg) < 2:
            continue

        # Uniform resampling within stage
        t_uni = np.linspace(t_stg[0], t_stg[-1], n_per_stage)
        e_uni = np.interp(t_uni, t_stg, e_stg)
        s_uni = np.interp(t_uni, t_stg, s_stg)

        # Stage strain range for normalization
        if NORMALIZE_PER_STAGE:
            strain_range = abs(e_stg[-1] - e_stg[0])
            scale = max(strain_range, 0.1)  # floor at 0.1% to avoid division by ~0
        else:
            scale = 1.0  # no normalization — use absolute residuals

        # Zone labelling based on time offset from stage start
        offsets = t_uni - t_uni[0]
        t_trans = T_TRANS_HIGH if sigma_mpa > 10.0 else T_TRANS_LOW

        weights = np.ones(n_per_stage)
        for i, dt in enumerate(offsets):
            if dt < T_ELASTIC:
                weights[i] = W_ELASTIC
            elif dt < t_trans:
                weights[i] = W_TRANSITION
            else:
                weights[i] = W_STEADY

        # Low-stress global multiplier
        if sigma_mpa < LOW_STRESS_THRESH:
            weights *= LOW_STRESS_MULT

        t_all.append(t_uni)
        e_all.append(e_uni)
        s_all.append(s_uni)
        w_all.append(weights)
        sc_all.append(np.full(n_per_stage, scale))
        sid_all.append(np.full(n_per_stage, k, dtype=int))

    return (np.concatenate(t_all), np.concatenate(e_all),
            np.concatenate(s_all), np.concatenate(w_all),
            np.concatenate(sc_all), np.concatenate(sid_all))


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          MODEL INTEGRATION                                 ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

def integrate_sic(t_s, sigma_Pa, n, A_pa_s, eta, E1):
    """
    Forward-Euler dislocation + exact-update Kelvin integration.
    Returns axial strain [%] at each time step.
    """
    n_pts = len(t_s)
    eps_k = 0.0
    eps_d = 0.0
    eps   = np.zeros(n_pts)

    for i in range(n_pts):
        sigma = sigma_Pa[i]
        dt    = (t_s[i] - t_s[i - 1]) if i > 0 else 0.0
        e_el  = sigma / E_ELASTIC

        if dt > 0:
            rate_dc = A_pa_s * (abs(sigma) ** n) * _EXP_QRT
            eps_d  += rate_dc * dt

            eps_k_eq  = sigma / E1 if E1 > 0 else 0.0
            decay     = np.exp(-E1 * dt / eta) if (eta > 0 and E1 > 0) else 0.0
            eps_k     = eps_k_eq + (eps_k - eps_k_eq) * decay

        eps[i] = (e_el + eps_d + eps_k) * 100.0

    return eps


def integrate_md(t_s, sigma_Pa, n, A_norm, K0, alpha_w, delta, m_md,
                 beta_w=None):
    """
    Munson-Dawson integration with adaptive sub-stepping.
    Returns axial strain [%] at each time step.
    """
    if beta_w is None:
        beta_w = BETA_W_MD
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
                Delta = alpha_w + beta_w * np.log10(max(som, 1e-30))
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


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          LAB DATA LOADING                                  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

print(f"Loading lab data for weighted optimizer ({MODEL.upper()} model) ...")
_lab = {}
for tid in TESTS:
    fp = os.path.join(DATA_DIR, f"ZW_{tid}.csv")
    time_h, strain_pct, sigma_diff, sigma3, _ = read_creep_csv(fp)
    sigma_diff_corr = correct_stress_artefacts(time_h, sigma_diff)

    # Build model time grid
    t_s, sigma_Pa, idx = build_stress_schedule(time_h, sigma_diff_corr, DT_HOURS)
    t_lab = time_h[idx:] - time_h[idx]
    e_lab = strain_pct[idx:] - strain_pct[idx]
    s_lab = sigma_diff_corr[idx:]

    # Detect stages
    stages = detect_stress_stages(t_lab, s_lab)

    # Resample and label
    t_resamp, e_resamp, s_resamp, weights, scales, stage_ids = \
        resample_and_label(t_lab, e_lab, s_lab, stages)

    # Subtract elastic strain from lab (to isolate creep)
    e_creep_lab = e_resamp - (s_resamp * 1e6 / E_ELASTIC) * 100.0

    _lab[tid] = dict(
        t_s=t_s, sigma_Pa=sigma_Pa,
        t_mod_h=t_s / HOUR,
        t_lab=t_lab, e_lab=e_lab,
        t_resamp=t_resamp, e_resamp=e_resamp,
        s_resamp=s_resamp,
        e_creep_lab=e_creep_lab,
        weights=weights, scales=scales,
        stage_ids=stage_ids,
        stages=stages,
    )

    print(f"  {tid}: {len(time_h)} raw pts -> {len(t_resamp)} resampled pts "
          f"({len(stages)} stages)")
    for k, stg in enumerate(stages):
        dur = stg['t_end'] - stg['t_start']
        print(f"    Stage {k}: {stg['t_start']:.0f}-{stg['t_end']:.0f} h, "
              f"sigma={stg['sigma_mpa']:.1f} MPa, dur={dur:.0f} h")


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          OBJECTIVE FUNCTION                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

_call_count  = [0]
_best        = [np.inf]
_best_params = [None]


def objective(x):
    if MODEL == "sic":
        n, log_A, log_eta, log_E1 = x
        A_mpa_yr = 10.0 ** log_A
        A_pa_s   = A_mpa_yr * (1e-6) ** n / _SEC_PER_YR
        eta      = 10.0 ** log_eta
        E1       = 10.0 ** log_E1
    else:
        n, log_A_norm, log_K0, alpha_w, delta, m_md, beta_w = x
        A_norm = 10.0 ** log_A_norm
        K0     = 10.0 ** log_K0

    total_weighted_sq = 0.0
    total_count = 0

    for tid in TESTS:
        d = _lab[tid]
        try:
            if MODEL == "sic":
                e_mod = integrate_sic(d['t_s'], d['sigma_Pa'], n, A_pa_s, eta, E1)
            else:
                e_mod = integrate_md(d['t_s'], d['sigma_Pa'],
                                     n, A_norm, K0, alpha_w, delta, m_md,
                                     beta_w=beta_w)
        except Exception:
            return 1e6

        # Interpolate model onto resampled lab times
        e_mod_resamp = np.interp(d['t_resamp'], d['t_mod_h'], e_mod)

        # Subtract elastic from model
        e_creep_mod = e_mod_resamp - (d['s_resamp'] * 1e6 / E_ELASTIC) * 100.0

        # Normalized, weighted residuals
        residuals = (e_creep_mod - d['e_creep_lab']) / d['scales']
        weighted_sq = (d['weights'] * residuals) ** 2

        total_weighted_sq += np.sum(weighted_sq)
        total_count += len(weighted_sq)

    obj = np.sqrt(total_weighted_sq / total_count)
    _call_count[0] += 1

    if obj < _best[0]:
        _best[0] = obj
        _best_params[0] = x.copy()
        if MODEL == "sic":
            print(f"  [{_call_count[0]:5d}] NEW BEST {obj:.6f}  "
                  f"n={n:.3f} A={10**log_A:.2f} eta={10**log_eta:.2e} "
                  f"E1={10**log_E1:.2e}", flush=True)
        else:
            print(f"  [{_call_count[0]:5d}] NEW BEST {obj:.6f}  "
                  f"n={n:.3f} log_A={log_A_norm:.2f} K0={10**log_K0:.3f} "
                  f"aw={alpha_w:.1f} d={delta:.4f} m={m_md:.3f} "
                  f"bw={beta_w:.2f}", flush=True)
    elif _call_count[0] % 500 == 0:
        print(f"  [{_call_count[0]:5d}] best so far {_best[0]:.6f}", flush=True)

    return obj


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          SEARCH BOUNDS & SEEDS                             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

if MODEL == "sic":
    BOUNDS = [
        (3.0,  7.0),    # n
        (0.0,  5.5),    # log10(A_mpa_yr)
        (10.0, 15.0),   # log10(eta [Pa·s])
        (7.5,  10.3),   # log10(E1 [Pa])
    ]
    # Seed with current best: n=4.034, A=736, eta=3.9e12, E1=6.81e8
    _SEED_ROW = np.array([4.034, np.log10(736.0), np.log10(3.9e12), np.log10(6.81e8)])
    _n_params = 4
    _popsize = 12
    _maxiter = 150
else:
    BOUNDS = [
        (2.0,   8.0),    # n
        (-9.0, -4.0),    # log10(A_norm)
        (-2.5,  3.5),    # log10(K0)
        (-50.0, 5.0),    # alpha_w
        (0.001, 3.0),    # delta
        (0.2,   4.0),    # m_md
        (-15.0, 0.0),    # beta_w (was fixed at -7.738)
    ]
    # Seed with previous weighted best + original beta_w
    _SEED_ROW = np.array([5.255, -6.49, np.log10(0.312), -10.176, 1.4508, 2.007, -7.738])
    _n_params = 7
    _popsize = 15
    _maxiter = 200

# Build initial population with seed
_n_pop = _popsize * _n_params
rng = np.random.default_rng(42)
_init_pop = rng.uniform(0, 1, size=(_n_pop, _n_params))
for j, (lo, hi) in enumerate(BOUNDS):
    _init_pop[:, j] = lo + _init_pop[:, j] * (hi - lo)
_init_pop[0] = _SEED_ROW


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          RUN OPTIMIZATION                                  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

print(f"\nStarting differential_evolution ({MODEL.upper()} model, weighted objective) ...")
print(f"  {_n_params} parameters, popsize={_popsize}, maxiter={_maxiter}")
print(f"  Weights: elastic={W_ELASTIC}, transition={W_TRANSITION}, steady={W_STEADY}")
print(f"  Low-stress mult={LOW_STRESS_MULT} (sigma < {LOW_STRESS_THRESH} MPa)")
print(f"  N_PER_STAGE={N_PER_STAGE}, T_ELASTIC={T_ELASTIC}h, "
      f"T_TRANS_HIGH={T_TRANS_HIGH}h, T_TRANS_LOW={T_TRANS_LOW}h")
print(f"  Seeded with current best params")

t0 = time.time()
result = differential_evolution(
    objective,
    bounds=BOUNDS,
    popsize=_popsize,
    maxiter=_maxiter,
    tol=1e-5,
    mutation=(0.5, 1.5),
    recombination=0.8,
    seed=42,
    init=_init_pop,
    polish=True,
    disp=False,
)
elapsed = time.time() - t0


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          REPORT & DIAGNOSTICS                              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

x = result.x

print(f"\n{'='*80}")
print(f"WEIGHTED OPTIMUM — {MODEL.upper()} MODEL  "
      f"(elapsed {elapsed/60:.1f} min, {_call_count[0]} evals)")
print(f"{'='*80}")

if MODEL == "sic":
    n, log_A, log_eta, log_E1 = x
    A_mpa_yr = 10.0 ** log_A
    eta      = 10.0 ** log_eta
    E1       = 10.0 ** log_E1
    A_pa_s   = A_mpa_yr * (1e-6) ** n / _SEC_PER_YR

    tau_h = eta / E1 / 3600
    kelvin_eq_21 = 21e6 / E1 * 100

    print(f"  n_SIC      = {n:.4f}")
    print(f"  A_SIC      = {A_mpa_yr:.2f} MPa^-n/yr")
    print(f"  ETA_KELVIN = {eta:.3e} Pa·s")
    print(f"  E1_KELVIN  = {E1:.3e} Pa")
    print(f"  Kelvin tau = {tau_h:.1f} h, eq@21MPa = {kelvin_eq_21:.2f}%")
else:
    n, log_A_norm, log_K0, alpha_w, delta, m_md, beta_w = x
    A_norm = 10.0 ** log_A_norm
    K0     = 10.0 ** log_K0
    A_pa_s = A_norm / ((_SIGMA_REF ** n) * _EXP_QRT)
    A_mpa_yr = A_pa_s * (1e6 ** n) * _SEC_PER_YR

    print(f"  n          = {n:.4f}")
    print(f"  A_norm     = {A_norm:.3e} /s")
    print(f"  A_mpa_yr   = {A_mpa_yr:.4g} MPa^-n/yr")
    print(f"  K0         = {K0:.4f}")
    print(f"  m_md       = {m_md:.4f}")
    print(f"  alpha_w    = {alpha_w:.4f}")
    print(f"  beta_w     = {beta_w:.4f}")
    print(f"  delta      = {delta:.6f}")

print(f"  Weighted objective = {result.fun:.6f}")


# ── Per-stage, per-zone RMSE breakdown ──────────────────────────────────────
print(f"\n{'─'*80}")
print("Per-stage, per-zone RMSE breakdown (creep strain, %):")
print(f"{'─'*80}")

for tid in TESTS:
    d = _lab[tid]
    if MODEL == "sic":
        e_mod = integrate_sic(d['t_s'], d['sigma_Pa'], n, A_pa_s, eta, E1)
    else:
        e_mod = integrate_md(d['t_s'], d['sigma_Pa'],
                             n, A_norm, K0, alpha_w, delta, m_md,
                             beta_w=beta_w)

    e_mod_resamp = np.interp(d['t_resamp'], d['t_mod_h'], e_mod)
    e_creep_mod = e_mod_resamp - (d['s_resamp'] * 1e6 / E_ELASTIC) * 100.0
    residuals = e_creep_mod - d['e_creep_lab']

    print(f"\n  {tid}:")
    for k, stg in enumerate(d['stages']):
        mask_stg = d['stage_ids'] == k
        if not np.any(mask_stg):
            continue
        res_stg = residuals[mask_stg]
        w_stg   = d['weights'][mask_stg]
        rmse_stg = np.sqrt(np.mean(res_stg ** 2))

        # Break down by zone (based on weight value)
        for zone_name, zone_w in [("elastic", W_ELASTIC),
                                   ("transition", W_TRANSITION),
                                   ("steady", W_STEADY)]:
            # Account for low-stress multiplier
            effective_w = zone_w
            if stg['sigma_mpa'] < LOW_STRESS_THRESH:
                effective_w *= LOW_STRESS_MULT
            mask_zone = mask_stg & np.isclose(d['weights'], effective_w, rtol=0.01)
            if not np.any(mask_zone):
                continue
            res_z = residuals[mask_zone]
            rmse_z = np.sqrt(np.mean(res_z ** 2))
            n_z = np.sum(mask_zone)
            print(f"    Stage {k} (sigma={stg['sigma_mpa']:.1f} MPa) "
                  f"{zone_name:>10}: RMSE={rmse_z:.4f}%  (n={n_z})")


# ── Traditional flat RMSE for comparison ────────────────────────────────────
print(f"\n{'─'*80}")
print("Traditional flat RMSE (for comparison with previous results):")
print(f"{'─'*80}")

for tid in TESTS:
    d = _lab[tid]
    # Re-run at DT=0.5h for fair comparison
    fp = os.path.join(DATA_DIR, f"ZW_{tid}.csv")
    time_h, strain_pct, sigma_diff, _, _ = read_creep_csv(fp)
    sigma_diff_corr = correct_stress_artefacts(time_h, sigma_diff)
    t_s_fine, sigma_Pa_fine, idx = build_stress_schedule(time_h, sigma_diff_corr, 0.5)
    t_lab = time_h[idx:] - time_h[idx]
    e_lab = strain_pct[idx:] - strain_pct[idx]
    t_mod_h = t_s_fine / HOUR

    if MODEL == "sic":
        e_mod = integrate_sic(t_s_fine, sigma_Pa_fine, n, A_pa_s, eta, E1)
    else:
        e_mod = integrate_md(t_s_fine, sigma_Pa_fine,
                             n, A_norm, K0, alpha_w, delta, m_md,
                             beta_w=beta_w)
    e_int = np.interp(t_lab, t_mod_h, e_mod)
    rmse = np.sqrt(np.mean((e_int - e_lab) ** 2))
    f_err = float(e_int[-1] - e_lab[-1])
    print(f"  {tid}: RMSE={rmse:.4f}%  final_err={f_err:+.3f}%")


# ── Copy-paste parameters ──────────────────────────────────────────────────
print(f"\n{'─'*80}")
print("To update run_calibration.py:")
print(f"{'─'*80}")

if MODEL == "sic":
    print(f"  _A_SIC_MPA_YR = {A_mpa_yr:.2f}")
    print(f"  N_SIC         = {n:.4f}")
    print(f"  ETA_KELVIN    = {eta:.4e}")
    print(f"  E1_KELVIN     = {E1:.4e}")
else:
    print(f"  _A_MPA_YR  = {A_mpa_yr:.1f}")
    print(f"  N_DISLOC   = {n:.4f}")
    print(f"  K0_MD      = {K0:.4f}")
    print(f"  M_MD       = {m_md:.4f}")
    print(f"  ALPHA_W_MD = {alpha_w:.4f}")
    print(f"  BETA_W_MD  = {beta_w:.4f}")
    print(f"  DELTA_MD   = {delta:.6f}")
