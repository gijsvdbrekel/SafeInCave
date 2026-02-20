"""Refined search with wider Kelvin range + finer MD grid."""

import os, json, sys
import numpy as np
from itertools import product

sys.path.insert(0, os.path.dirname(__file__))
from run_calibration import (
    read_creep_csv, build_stress_schedule,
    disloc_rate, kelvin_rate, desai_rate, _update_alpha,
    munsondawson_rate,
    R_GAS, T_KELVIN, Q_DISLOC, E_ELASTIC, NU_ELASTIC,
    A_SIC, N_SIC, A_MD, N_MD, MU_SHEAR,
    N_DESAI_POWER, BETA1_DESAI, BETA_DESAI, M_DESAI,
    GAMMA_DESAI, SIGMA_T_DESAI,
    C_MD, M_MD, BETA_W_MD, DELTA_MD,
    DATA_DIR, HOUR,
    _DESAI_MAX_DEPS, _MD_MAX_DZETA,
)

DT_HOURS = 0.5
N1_DESAI = 2.0; A1_DESAI = 1e-3; ETA_DESAI = 1.2; NU1_KELVIN = 0.25

# ── REFINED SEARCH GRIDS ──
# Kelvin: go much higher on E1 to reduce amplitude, higher eta to slow it
ETA_KELVIN_GRID = [1e13, 2e13, 5e13, 1e14]
E1_KELVIN_GRID  = [1.0e9, 2.0e9, 3.0e9, 5.0e9]

# Desai: fine-tune
ALPHA0_GRID = [0.003, 0.005, 0.008, 0.01]
MU1_GRID    = [5e-16, 1e-15, 3e-15]

# MD: fine-tune around K0=7e4
K0_GRID      = [4e4, 6e4, 7e4, 8e4, 1e5]
ALPHA_W_GRID = [-14, -10, -6, -2, 2]


def integrate_sic(t_s, sigma_Pa, sigma3_Pa, eta_k, E1_k, alpha0, mu1):
    n = len(t_s)
    ek = ed = ev = 0.0
    alpha = alpha0; zeta = 0.0
    eps = np.zeros(n)
    eps_kv = np.zeros(n)
    
    for i in range(n):
        sigma = sigma_Pa[i]
        dt = (t_s[i] - t_s[i-1]) if i > 0 else 0.0
        e_el = sigma / E_ELASTIC
        if dt > 0:
            ed += disloc_rate(sigma, A_SIC, N_SIC, Q_DISLOC, T_KELVIN) * dt
            ek += kelvin_rate(ek, sigma, eta_k, E1_k) * dt
            if mu1 > 0:
                dt_rem = dt
                while dt_rem > 1e-10:
                    rate = desai_rate(sigma, sigma3_Pa, alpha, mu1, N1_DESAI,
                                     A1_DESAI, ETA_DESAI, N_DESAI_POWER,
                                     BETA1_DESAI, BETA_DESAI, M_DESAI,
                                     GAMMA_DESAI, SIGMA_T_DESAI)
                    if rate < 1e-30: break
                    dt_sub = min(dt_rem, _DESAI_MAX_DEPS / rate)
                    dt_sub = max(dt_sub, 1e-6)
                    ev += rate * dt_sub
                    zeta += rate * dt_sub
                    alpha = _update_alpha(zeta, A1_DESAI, alpha0, ETA_DESAI)
                    dt_rem -= dt_sub
        eps[i] = (e_el + ed + ek + ev) * 100
        eps_kv[i] = ek * 100
    return eps, eps_kv


def integrate_md(t_s, sigma_Pa, K0, alpha_w):
    n = len(t_s)
    ec = 0.0; zeta = 0.0
    eps = np.zeros(n)
    for i in range(n):
        sigma = sigma_Pa[i]
        dt = (t_s[i] - t_s[i-1]) if i > 0 else 0.0
        e_el = sigma / E_ELASTIC
        if dt > 0:
            rate_ss = disloc_rate(sigma, A_MD, N_MD, Q_DISLOC, T_KELVIN)
            dt_rem = dt
            while dt_rem > 1e-10:
                rate_md, F = munsondawson_rate(
                    sigma, zeta, A_MD, N_MD, Q_DISLOC, T_KELVIN,
                    K0, C_MD, M_MD, alpha_w, BETA_W_MD, DELTA_MD, MU_SHEAR)
                zeta_rate = abs(F - 1.0) * rate_ss
                dt_sub = min(dt_rem, _MD_MAX_DZETA / zeta_rate) if zeta_rate > 1e-30 else dt_rem
                dt_sub = max(dt_sub, 1e-6); dt_sub = min(dt_sub, dt_rem)
                ec += rate_md * dt_sub
                zeta += (F - 1.0) * rate_ss * dt_sub
                zeta = max(zeta, 0.0)
                dt_rem -= dt_sub
        eps[i] = (e_el + ec) * 100
    return eps


def diagnostics(t_mod_h, e_mod, t_lab_h, e_lab):
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
    sigma3_Pa = float(sigma3[idx_start]) * 1e6
    t0_h = time_h[idx_start]
    lab_t = time_h[idx_start:] - t0_h
    lab_e = strain_pct[idx_start:] - strain_pct[idx_start]
    t_mod_h = t_s / HOUR
    
    # ── SIC ──
    print("=" * 130)
    print("SIC REFINED SEARCH")
    print("=" * 130)
    hdr = f"{'eta_K':>10} {'E1_K':>10} {'a0':>6} {'mu1':>10} | {'RMSE':>6} {'P1':>6} {'P2':>6} {'P3':>6} {'OS100':>7} {'drpErr':>7} {'final':>7} {'kvMax':>6}"
    print(hdr); print("-"*130)
    
    results = []
    total = len(ETA_KELVIN_GRID)*len(E1_KELVIN_GRID)*len(ALPHA0_GRID)*len(MU1_GRID)
    done = 0
    
    for eta_k, E1_k, a0, mu1 in product(ETA_KELVIN_GRID, E1_KELVIN_GRID, ALPHA0_GRID, MU1_GRID):
        e_tot, e_kv = integrate_sic(t_s, sigma_fine_Pa, sigma3_Pa, eta_k, E1_k, a0, mu1)
        rmse, rp1, rp2, rp3, os100, drp, fin = diagnostics(t_mod_h, e_tot, lab_t, lab_e)
        kv_max = e_kv.max()
        results.append((rmse, eta_k, E1_k, a0, mu1, rp1, rp2, rp3, os100, drp, fin, kv_max))
        done += 1
        if done % 48 == 0:
            print(f"  [{done}/{total}]", flush=True)
    
    results.sort(key=lambda x: x[0])
    for rmse, eta_k, E1_k, a0, mu1, rp1, rp2, rp3, os100, drp, fin, kvmax in results[:40]:
        print(f"{eta_k:10.2e} {E1_k:10.2e} {a0:6.3f} {mu1:10.2e} | "
              f"{rmse:6.3f} {rp1:6.3f} {rp2:6.3f} {rp3:6.3f} {os100:+7.3f} {drp:+7.3f} {fin:+7.3f} {kvmax:6.2f}")
    
    # ── MD ──
    print("\n" + "="*100)
    print("MD REFINED SEARCH")
    print("="*100)
    hdr = f"{'K0':>10} {'aw':>6} | {'RMSE':>6} {'P1':>6} {'P2':>6} {'P3':>6} {'OS100':>7} {'drpErr':>7} {'final':>7}"
    print(hdr); print("-"*100)
    
    md_res = []
    for K0, aw in product(K0_GRID, ALPHA_W_GRID):
        e_tot = integrate_md(t_s, sigma_fine_Pa, K0, aw)
        rmse, rp1, rp2, rp3, os100, drp, fin = diagnostics(t_mod_h, e_tot, lab_t, lab_e)
        md_res.append((rmse, K0, aw, rp1, rp2, rp3, os100, drp, fin))
    
    md_res.sort(key=lambda x: x[0])
    for rmse, K0, aw, rp1, rp2, rp3, os100, drp, fin in md_res:
        print(f"{K0:10.2e} {aw:6.1f} | {rmse:6.3f} {rp1:6.3f} {rp2:6.3f} {rp3:6.3f} {os100:+7.3f} {drp:+7.3f} {fin:+7.3f}")


if __name__ == "__main__":
    main()
