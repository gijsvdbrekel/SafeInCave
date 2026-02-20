"""Search MD parameters including m and delta, with K0 scaled per m."""

import os, sys
import numpy as np
from itertools import product

sys.path.insert(0, os.path.dirname(__file__))
from run_calibration import (
    read_creep_csv, build_stress_schedule,
    disloc_rate,
    R_GAS, T_KELVIN, Q_DISLOC, E_ELASTIC,
    A_MD, N_MD, MU_SHEAR,
    C_MD, BETA_W_MD,
    DATA_DIR, HOUR,
    _MD_MAX_DZETA,
)

DT_HOURS = 0.5

# ── SEARCH GRIDS ──
# K0 must be scaled with m to keep eps_t_star in the right range.
# Reference: K0=6.5e4 at m=3 gives eps_t*(20MPa)=0.028.
# For each m, K0_ref = 0.028 / (exp(c*T) * (sigma_20/mu)^m)
_REF_EPS_T = 0.028
_EXP_CT = np.exp(C_MD * T_KELVIN)
_RATIO_20 = 20e6 / MU_SHEAR

M_GRID = [1.0, 1.25, 1.5, 2.0, 2.5, 3.0]
# For each m, search K0 around the reference value (factors 0.5, 0.75, 1.0, 1.5, 2.0)
K0_FACTORS = [0.5, 0.75, 1.0, 1.5, 2.0]
AW_GRID    = [-14, -10, -6, -2]
DELTA_GRID = [0.1, 0.2, 0.4, 0.58]


def integrate_md(t_s, sigma_Pa, K0, alpha_w, m_md, delta):
    n = len(t_s)
    ec = 0.0; zeta = 0.0
    eps = np.zeros(n)
    for i in range(n):
        sigma = sigma_Pa[i]
        dt = (t_s[i] - t_s[i-1]) if i > 0 else 0.0
        e_el = sigma / E_ELASTIC
        if dt > 0:
            rate_ss = disloc_rate(sigma, A_MD, N_MD, Q_DISLOC, T_KELVIN)
            if rate_ss < 1e-50:
                eps[i] = (e_el + ec) * 100
                continue
            dt_rem = dt
            while dt_rem > 1e-10:
                sigma_over_mu = abs(sigma) / MU_SHEAR if MU_SHEAR > 0 else 0.0
                eps_t_star = K0 * _EXP_CT * (sigma_over_mu ** m_md)
                if eps_t_star < 1e-50:
                    ec += rate_ss * dt_rem
                    break
                log_arg = max(sigma_over_mu, 1e-30)
                Delta = alpha_w + BETA_W_MD * np.log10(log_arg)
                ratio = zeta / eps_t_star
                if ratio <= 1.0:
                    F = np.exp(Delta * (1.0 - ratio) ** 2)
                else:
                    F = np.exp(-delta * (1.0 - ratio) ** 2)
                rate_md = F * rate_ss
                zeta_rate = abs(F - 1.0) * rate_ss
                dt_sub = min(dt_rem, _MD_MAX_DZETA / zeta_rate) if zeta_rate > 1e-30 else dt_rem
                dt_sub = max(dt_sub, 10.0); dt_sub = min(dt_sub, dt_rem)
                ec += rate_md * dt_sub
                zeta += (F - 1.0) * rate_ss * dt_sub
                zeta = max(zeta, 0.0)
                dt_rem -= dt_sub
        eps[i] = (e_el + ec) * 100
    return eps


def diag(t_mod_h, e_mod, t_lab_h, e_lab):
    e_int = np.interp(t_lab_h, t_mod_h, e_mod)
    rmse = np.sqrt(np.mean((e_int - e_lab)**2))
    m1 = t_lab_h < 142; m2 = (t_lab_h >= 142) & (t_lab_h < 670); m3 = t_lab_h >= 670
    rp1 = np.sqrt(np.mean((e_int[m1]-e_lab[m1])**2)) if m1.any() else 0
    rp2 = np.sqrt(np.mean((e_int[m2]-e_lab[m2])**2)) if m2.any() else 0
    rp3 = np.sqrt(np.mean((e_int[m3]-e_lab[m3])**2)) if m3.any() else 0
    i100 = np.argmin(np.abs(t_lab_h-100)); os100 = e_int[i100]-e_lab[i100]
    i660 = np.argmin(np.abs(t_lab_h-660)); i700 = np.argmin(np.abs(t_lab_h-700))
    drop_err = (e_int[i700]-e_int[i660]) - (e_lab[i700]-e_lab[i660])
    final = e_int[-1]-e_lab[-1]
    return rmse, rp1, rp2, rp3, os100, drop_err, final


def main():
    filepath = os.path.join(DATA_DIR, "ZW_TCC1.csv")
    time_h, strain_pct, sigma_diff, sigma3, temp_C = read_creep_csv(filepath)
    t_s, sigma_fine_Pa, idx_start = build_stress_schedule(time_h, sigma_diff, DT_HOURS)
    t0_h = time_h[idx_start]
    lab_t = time_h[idx_start:] - t0_h
    lab_e = strain_pct[idx_start:] - strain_pct[idx_start]
    t_mod_h = t_s / HOUR

    print("=" * 150)
    print("MD SEARCH — varying m, delta, K0 (scaled), alpha_w")
    print("=" * 150)

    # Print K0 reference values
    for m_md in M_GRID:
        K0_ref = _REF_EPS_T / (_EXP_CT * (_RATIO_20 ** m_md))
        print(f"  m={m_md:.2f}: K0_ref = {K0_ref:.2e}")
    print()

    results = []
    total = len(M_GRID)*len(K0_FACTORS)*len(AW_GRID)*len(DELTA_GRID)
    done = 0

    for m_md in M_GRID:
        K0_ref = _REF_EPS_T / (_EXP_CT * (_RATIO_20 ** m_md))
        for k_fac, aw, delta in product(K0_FACTORS, AW_GRID, DELTA_GRID):
            K0 = K0_ref * k_fac
            e_tot = integrate_md(t_s, sigma_fine_Pa, K0, aw, m_md, delta)
            rmse, rp1, rp2, rp3, os100, drp, fin = diag(t_mod_h, e_tot, lab_t, lab_e)
            results.append((rmse, K0, aw, m_md, delta, rp1, rp2, rp3, os100, drp, fin))
            done += 1
            if done % 50 == 0:
                print(f"  [{done}/{total}]", flush=True)

    results.sort(key=lambda x: x[0])

    hdr = f"{'K0':>10} {'aw':>6} {'m':>5} {'delta':>6} | {'RMSE':>6} {'P1':>6} {'P2':>6} {'P3':>6} {'OS100':>7} {'drpErr':>7} {'final':>7}"
    print(hdr); print("-"*150)
    for r in results[:80]:
        rmse, K0, aw, m_md, delta, rp1, rp2, rp3, os100, drp, fin = r
        print(f"{K0:10.2e} {aw:6.1f} {m_md:5.2f} {delta:6.2f} | "
              f"{rmse:6.3f} {rp1:6.3f} {rp2:6.3f} {rp3:6.3f} {os100:+7.3f} {drp:+7.3f} {fin:+7.3f}")


if __name__ == "__main__":
    main()
