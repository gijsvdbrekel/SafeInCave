"""
Joint MD parameter search: TCC1 + TCC2 + TCC6.

Rationale
---------
Current calibration (n=5, A=40 MPa^-5/yr) fits TCC1 (RMSE=0.37%) and TCC2
(RMSE=0.92%) well.  TCC6 has the same σ3=18 MPa but Stage 1 lasts 203 h
(vs 142 h for TCC1), causing a 6.0% RMSE from over-predicted dislocation
strain.

Strategy: search (n, A_factor) jointly with MD transient params (K0, α_w,
δ).  Objective is a weighted combined RMSE across all three tests.  Lower n
reduces Stage-1 sensitivity; A_factor compensates to keep the overall level
correct.  Transient params control the initial-loading burst.

Grid sizes
----------
n         : 5 values   → dislocation exponent
A_factor  : 5 values   → multiplied by current A_MD from run_calibration.py
K0_factor : 3 values   → multiplied by K0_ref(m) so eps_t_star stays sane
alpha_w   : 3 values   → work-hardening amplitude
delta     : 3 values   → recovery exponent
m_md      : fixed at current M_MD (only affects eps_t_star scaling)

Total: 5×5×3×3×3 = 675 combos × 3 tests = 2025 integrations (~2–5 min).

Output: table sorted by combined RMSE, showing per-test RMSEs.
"""

import os
import sys
import numpy as np
from itertools import product

sys.path.insert(0, os.path.dirname(__file__))
from run_calibration import (
    read_creep_csv, build_stress_schedule,
    R_GAS, T_KELVIN, Q_DISLOC, E_ELASTIC,
    A_MD, N_MD, MU_SHEAR, C_MD, BETA_W_MD,
    M_MD, ALPHA_W_MD, DELTA_MD,
    DATA_DIR, HOUR,
    _MD_MAX_DZETA,
)

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

TESTS   = ["TCC1", "TCC2", "TCC6"]
WEIGHTS = {"TCC1": 1.0, "TCC2": 1.0, "TCC6": 1.0}   # equal weight
DT_HOURS = 0.5

# ── Primary grid (steady-state creep rate) ────────────────────────────────
N_GRID    = [4.5, 4.75, 5.0, 5.25, 5.5]
A_FACTORS = [0.60, 0.75, 1.0, 1.25, 1.50]   # multiplied by current A_MD

# ── Transient grid (m_md kept fixed at current value) ────────────────────
# K0 is scaled per m_md so eps_t_star(20 MPa) stays at _REF_EPS_T.
_REF_EPS_T = 0.028
_EXP_CT    = np.exp(C_MD * T_KELVIN)
_RATIO_20  = 20e6 / MU_SHEAR
K0_FACTORS = [0.5, 1.0, 2.0]
AW_GRID    = [-20.0, -17.0, -14.0]
DELTA_GRID = [0.15, 0.25, 0.40]
M_FIXED    = M_MD    # keep m fixed; override here if desired


# ══════════════════════════════════════════════════════════════════════════════
# INTEGRATION
# ══════════════════════════════════════════════════════════════════════════════

def integrate_md(t_s, sigma_Pa, n, A, K0, alpha_w, delta):
    """
    Forward-Euler Munson-Dawson integration with explicit (n, A, K0, α_w, δ).
    m_md is taken from M_FIXED (module level).  All other MD parameters
    (C_MD, BETA_W_MD, MU_SHEAR, Q_DISLOC, T_KELVIN, E_ELASTIC) are shared
    from run_calibration.py.

    Returns axial strain [%] at each timestep.
    """
    n_pts = len(t_s)
    eps_c  = 0.0   # accumulated total creep
    zeta   = 0.0   # MD transient internal variable
    eps = np.zeros(n_pts)

    for i in range(n_pts):
        sigma = sigma_Pa[i]
        dt    = (t_s[i] - t_s[i - 1]) if i > 0 else 0.0
        e_el  = sigma / E_ELASTIC

        if dt > 0:
            rate_ss = A * (abs(sigma) ** n) * np.exp(-Q_DISLOC / (R_GAS * T_KELVIN))
            if rate_ss < 1e-50:
                eps[i] = (e_el + eps_c) * 100.0
                continue

            dt_rem = dt
            while dt_rem > 1e-10:
                sigma_over_mu = abs(sigma) / MU_SHEAR if MU_SHEAR > 0 else 0.0
                eps_t_star    = K0 * _EXP_CT * (sigma_over_mu ** M_FIXED)
                if eps_t_star < 1e-50:
                    eps_c += rate_ss * dt_rem
                    break
                log_arg = max(sigma_over_mu, 1e-30)
                Delta   = alpha_w + BETA_W_MD * np.log10(log_arg)
                ratio   = zeta / eps_t_star
                if ratio <= 1.0:
                    F = np.exp(Delta * (1.0 - ratio) ** 2)
                else:
                    F = np.exp(-delta * (1.0 - ratio) ** 2)

                rate_md   = F * rate_ss
                zeta_rate = abs(F - 1.0) * rate_ss
                if zeta_rate > 1e-30:
                    dt_sub = min(dt_rem, _MD_MAX_DZETA / zeta_rate)
                else:
                    dt_sub = dt_rem
                dt_sub = max(dt_sub, 10.0)
                dt_sub = min(dt_sub, dt_rem)

                eps_c += rate_md * dt_sub
                zeta  += (F - 1.0) * rate_ss * dt_sub
                zeta   = max(zeta, 0.0)
                dt_rem -= dt_sub

        eps[i] = (e_el + eps_c) * 100.0

    return eps


# ══════════════════════════════════════════════════════════════════════════════
# DIAGNOSTICS
# ══════════════════════════════════════════════════════════════════════════════

def rmse(t_mod_h, e_mod, t_lab_h, e_lab):
    """RMSE of model vs lab at lab time points [%]."""
    e_int = np.interp(t_lab_h, t_mod_h, e_mod)
    return float(np.sqrt(np.mean((e_int - e_lab) ** 2)))


def final_err(t_mod_h, e_mod, t_lab_h, e_lab):
    """Signed error at the last lab time point [%]."""
    e_int = np.interp(t_lab_h, t_mod_h, e_mod)
    return float(e_int[-1] - e_lab[-1])


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    # ── Load lab data for each target test ───────────────────────────────────
    lab = {}
    for test_id in TESTS:
        filepath = os.path.join(DATA_DIR, f"ZW_{test_id}.csv")
        time_h, strain_pct, sigma_diff, sigma3, temp_C = read_creep_csv(filepath)
        t_s, sigma_Pa, idx = build_stress_schedule(time_h, sigma_diff, DT_HOURS)
        t0_h       = time_h[idx]
        t_lab_h    = time_h[idx:] - t0_h
        e_lab      = strain_pct[idx:] - strain_pct[idx]
        t_mod_h    = t_s / HOUR
        lab[test_id] = dict(t_s=t_s, sigma_Pa=sigma_Pa,
                            t_lab_h=t_lab_h, e_lab=e_lab,
                            t_mod_h=t_mod_h)
        print(f"  Loaded {test_id}: {len(time_h)} pts, "
              f"{time_h[-1]:.0f} h, σ3={sigma3[idx]:.0f} MPa")

    # Compute K0_ref for the fixed m value
    K0_ref = _REF_EPS_T / (_EXP_CT * (_RATIO_20 ** M_FIXED))
    W_total = sum(WEIGHTS[t] for t in TESTS)

    total_combos = (len(N_GRID) * len(A_FACTORS) *
                    len(K0_FACTORS) * len(AW_GRID) * len(DELTA_GRID))
    print(f"\n  m_md fixed at {M_FIXED:.3f}, K0_ref = {K0_ref:.2e}")
    print(f"  Total combinations: {total_combos}")
    print(f"  Objective: weighted combined RMSE (weights {dict(WEIGHTS)})\n")

    results = []
    done    = 0

    for n_val, a_fac, k0_fac, aw, dlt in product(
            N_GRID, A_FACTORS, K0_FACTORS, AW_GRID, DELTA_GRID):

        A_val = A_MD * a_fac            # effective A in Pa^{-n}/s
        # Re-scale A so that at the reference stress 20 MPa and reference n,
        # the rate is unchanged.  This keeps the effective creep speed in
        # the field comparable across different n values.
        # (optional: just scale by a_fac directly — done here without re-normalization)
        K0    = K0_ref * k0_fac

        per_test = {}
        combined = 0.0
        for test_id in TESTS:
            d = lab[test_id]
            e_mod = integrate_md(d["t_s"], d["sigma_Pa"],
                                 n_val, A_val, K0, aw, dlt)
            r = rmse(d["t_mod_h"], e_mod, d["t_lab_h"], d["e_lab"])
            f = final_err(d["t_mod_h"], e_mod, d["t_lab_h"], d["e_lab"])
            per_test[test_id] = (r, f)
            combined += WEIGHTS[test_id] * r

        combined /= W_total
        results.append((combined, n_val, a_fac, K0, aw, dlt, per_test))

        done += 1
        if done % 100 == 0:
            print(f"  [{done}/{total_combos}] best so far = "
                  f"{min(r[0] for r in results):.4f}%", flush=True)

    results.sort(key=lambda x: x[0])

    # ── Print results table ───────────────────────────────────────────────
    print("\n" + "=" * 130)
    print("TOP RESULTS — sorted by combined RMSE across TCC1 + TCC2 + TCC6")
    print(f"  Current baseline: n={N_MD:.2f}, A_fac=1.00, "
          f"K0={K0_ref:.2e}×1.0, α_w={ALPHA_W_MD:.1f}, δ={DELTA_MD:.3f}")
    print("=" * 130)
    hdr = (f"{'n':>5} {'A_fac':>6} {'K0':>10} {'α_w':>6} {'δ':>5} | "
           f"{'combined':>9} | "
           f"{'TCC1_rmse':>10} {'TCC1_fin':>9} | "
           f"{'TCC2_rmse':>10} {'TCC2_fin':>9} | "
           f"{'TCC6_rmse':>10} {'TCC6_fin':>9}")
    print(hdr)
    print("-" * 130)
    for row in results[:80]:
        comb, n_val, a_fac, K0, aw, dlt, per = row
        r1, f1 = per["TCC1"]
        r2, f2 = per["TCC2"]
        r6, f6 = per["TCC6"]
        print(f"{n_val:5.2f} {a_fac:6.2f} {K0:10.2e} {aw:6.1f} {dlt:5.2f} | "
              f"{comb:9.4f} | "
              f"{r1:10.4f} {f1:+9.3f} | "
              f"{r2:10.4f} {f2:+9.3f} | "
              f"{r6:10.4f} {f6:+9.3f}")

    # ── Baseline for comparison ───────────────────────────────────────────
    print("\n" + "-" * 130)
    print("BASELINE (current parameters):")
    A_base = A_MD
    K0_base = K0_ref   # k0_fac=1.0
    comb_base = 0.0
    for test_id in TESTS:
        d = lab[test_id]
        e_mod = integrate_md(d["t_s"], d["sigma_Pa"],
                             N_MD, A_base, K0_base, ALPHA_W_MD, DELTA_MD)
        r = rmse(d["t_mod_h"], e_mod, d["t_lab_h"], d["e_lab"])
        f = final_err(d["t_mod_h"], e_mod, d["t_lab_h"], d["e_lab"])
        comb_base += WEIGHTS[test_id] * r
        print(f"  {test_id}: RMSE={r:.4f}%  final_err={f:+.3f}%")
    comb_base /= W_total
    print(f"  Combined RMSE = {comb_base:.4f}%")

    # ── Best result detail ────────────────────────────────────────────────
    best = results[0]
    comb, n_val, a_fac, K0, aw, dlt, per = best
    print("\n" + "=" * 130)
    print("BEST PARAMETER SET:")
    print(f"  n       = {n_val:.4f}")
    # Recover effective A in MPa^{-n}/yr units for easy comparison
    sec_per_year = 365.25 * 24 * 3600
    A_eff_MPa_yr = A_MD * a_fac * (1e6 ** n_val) * sec_per_year
    print(f"  A_fac   = {a_fac:.4f}  (A_eff = {A_eff_MPa_yr:.1f} MPa^-n/yr)")
    print(f"  K0      = {K0:.4e}")
    print(f"  m_md    = {M_FIXED:.4f}  (fixed)")
    print(f"  α_w     = {aw:.4f}")
    print(f"  δ       = {dlt:.4f}")
    print(f"  Combined RMSE = {comb:.4f}%  (baseline: {comb_base:.4f}%)")
    print(f"  Per-test:")
    for test_id in TESTS:
        r, f = per[test_id]
        print(f"    {test_id}: RMSE={r:.4f}%  final_err={f:+.3f}%")
    print("=" * 130)


if __name__ == "__main__":
    main()
