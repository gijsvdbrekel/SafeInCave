"""
Systematic parameter search for TCC1 calibration.
Sweeps Kelvin (eta, E1) + Desai (alpha0, mu1) + MD (K0, alpha_w) parameters.
Reports full-curve RMSE and key diagnostics.
"""

import os, json, sys
import numpy as np
from itertools import product

# ── Import from run_calibration ──
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

# ── Fixed dislocation creep params (already calibrated) ──
# A=40, n=5.0, Q/R=6252 — these stay fixed

# ── SEARCH GRIDS ──
# Kelvin
ETA_KELVIN_GRID  = [5e12, 1e13, 2e13, 5e13]
E1_KELVIN_GRID   = [0.2e9, 0.35e9, 0.5e9, 0.8e9]
NU1_KELVIN = 0.25

# Desai
ALPHA0_GRID = [0.003, 0.005, 0.008, 0.012]
MU1_GRID    = [5e-16, 1e-15, 5e-15]
N1_DESAI    = 2.0
A1_DESAI    = 1e-3
ETA_DESAI   = 1.2

# MD
K0_GRID     = [5e4, 7e4, 1e5, 2e5]
ALPHA_W_GRID = [-14, -10, -5, 0]


def integrate_sic(t_s, sigma_Pa, sigma3_Pa, eta_k, E1_k, alpha0, mu1):
    n_steps = len(t_s)
    eps_el = np.zeros(n_steps)
    eps_dc = np.zeros(n_steps)
    eps_kv = np.zeros(n_steps)
    eps_vp = np.zeros(n_steps)
    
    ek = 0.0; ed = 0.0; ev = 0.0
    alpha = alpha0; zeta = 0.0
    
    for i in range(n_steps):
        sigma = sigma_Pa[i]
        dt = (t_s[i] - t_s[i-1]) if i > 0 else 0.0
        
        eps_el[i] = sigma / E_ELASTIC
        
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
        
        eps_dc[i] = ed; eps_kv[i] = ek; eps_vp[i] = ev
    
    return (eps_el + eps_dc + eps_kv + eps_vp) * 100, eps_dc * 100, eps_kv * 100, eps_vp * 100


def integrate_md(t_s, sigma_Pa, K0, alpha_w):
    n_steps = len(t_s)
    eps_el = np.zeros(n_steps)
    eps_c_arr = np.zeros(n_steps)
    
    ec = 0.0; zeta = 0.0
    
    for i in range(n_steps):
        sigma = sigma_Pa[i]
        dt = (t_s[i] - t_s[i-1]) if i > 0 else 0.0
        eps_el[i] = sigma / E_ELASTIC
        
        if dt > 0:
            rate_ss = disloc_rate(sigma, A_MD, N_MD, Q_DISLOC, T_KELVIN)
            dt_rem = dt
            while dt_rem > 1e-10:
                rate_md, F = munsondawson_rate(
                    sigma, zeta, A_MD, N_MD, Q_DISLOC, T_KELVIN,
                    K0, C_MD, M_MD, alpha_w, BETA_W_MD, DELTA_MD, MU_SHEAR)
                zeta_rate = abs(F - 1.0) * rate_ss
                if zeta_rate > 1e-30:
                    dt_sub = min(dt_rem, _MD_MAX_DZETA / zeta_rate)
                else:
                    dt_sub = dt_rem
                dt_sub = max(dt_sub, 1e-6)
                dt_sub = min(dt_sub, dt_rem)
                ec += rate_md * dt_sub
                zeta += (F - 1.0) * rate_ss * dt_sub
                zeta = max(zeta, 0.0)
                dt_rem -= dt_sub
        
        eps_c_arr[i] = ec
    
    return (eps_el + eps_c_arr) * 100


def compute_rmse(t_model_h, e_model_pct, t_lab_h, e_lab_pct):
    e_interp = np.interp(t_lab_h, t_model_h, e_model_pct)
    return np.sqrt(np.mean((e_interp - e_lab_pct)**2))


def compute_diagnostics(t_model_h, e_model_pct, t_lab_h, e_lab_pct, sig_lab):
    """Compute key diagnostic metrics beyond overall RMSE."""
    e_interp = np.interp(t_lab_h, t_model_h, e_model_pct)
    
    # Phase 1 RMSE: 20 MPa (0-142h)
    m1 = t_lab_h < 142
    rmse_p1 = np.sqrt(np.mean((e_interp[m1] - e_lab_pct[m1])**2)) if m1.any() else 0
    
    # Phase 2 RMSE: 10 MPa (142-670h)
    m2 = (t_lab_h >= 142) & (t_lab_h < 670)
    rmse_p2 = np.sqrt(np.mean((e_interp[m2] - e_lab_pct[m2])**2)) if m2.any() else 0
    
    # Phase 3 RMSE: 2-3 MPa (670h+)
    m3 = t_lab_h >= 670
    rmse_p3 = np.sqrt(np.mean((e_interp[m3] - e_lab_pct[m3])**2)) if m3.any() else 0
    
    # Overshoot at end of 20 MPa
    idx_100h = np.argmin(np.abs(t_lab_h - 100))
    overshoot_100h = e_interp[idx_100h] - e_lab_pct[idx_100h]
    
    # Drop at 10->2 MPa transition
    idx_660 = np.argmin(np.abs(t_lab_h - 660))
    idx_700 = np.argmin(np.abs(t_lab_h - 700))
    model_drop = e_interp[idx_700] - e_interp[idx_660]
    lab_drop = e_lab_pct[idx_700] - e_lab_pct[idx_660]
    drop_error = model_drop - lab_drop  # negative = model drops too much
    
    final_err = e_interp[-1] - e_lab_pct[-1]
    
    return {
        'rmse_total': compute_rmse(t_model_h, e_model_pct, t_lab_h, e_lab_pct),
        'rmse_p1': rmse_p1, 'rmse_p2': rmse_p2, 'rmse_p3': rmse_p3,
        'overshoot_100h': overshoot_100h,
        'drop_error': drop_error,
        'final_err': final_err,
    }


def main():
    # Load TCC1 lab data
    filepath = os.path.join(DATA_DIR, "ZW_TCC1.csv")
    time_h, strain_pct, sigma_diff, sigma3, temp_C = read_creep_csv(filepath)
    t_s, sigma_fine_Pa, idx_start = build_stress_schedule(time_h, sigma_diff, DT_HOURS)
    sigma3_Pa = float(sigma3[idx_start]) * 1e6
    
    t0_h = time_h[idx_start]
    lab_t_h = time_h[idx_start:] - t0_h
    lab_e_pct = strain_pct[idx_start:] - strain_pct[idx_start]
    lab_sig = sigma_diff[idx_start:]
    t_model_h = t_s / HOUR
    
    print("=" * 120)
    print("TCC1 PARAMETER SEARCH — SIC (Kelvin + Desai) + MD (K0, alpha_w)")
    print("=" * 120)
    
    # ── SIC SEARCH ──
    print("\n--- SafeInCave search ---")
    print(f"{'eta_K':>10} {'E1_K':>10} {'alpha0':>8} {'mu1':>10} | {'RMSE':>7} {'P1':>7} {'P2':>7} {'P3':>7} {'OS@100h':>8} {'drop_err':>9} {'final':>7}")
    print("-" * 120)
    
    sic_results = []
    
    for eta_k, E1_k, alpha0, mu1 in product(ETA_KELVIN_GRID, E1_KELVIN_GRID, ALPHA0_GRID, MU1_GRID):
        e_tot, _, _, _ = integrate_sic(t_s, sigma_fine_Pa, sigma3_Pa, eta_k, E1_k, alpha0, mu1)
        diag = compute_diagnostics(t_model_h, e_tot, lab_t_h, lab_e_pct, lab_sig)
        
        sic_results.append((eta_k, E1_k, alpha0, mu1, diag))
    
    # Sort by total RMSE
    sic_results.sort(key=lambda x: x[4]['rmse_total'])
    
    for eta_k, E1_k, alpha0, mu1, diag in sic_results[:30]:
        print(f"{eta_k:10.2e} {E1_k:10.2e} {alpha0:8.4f} {mu1:10.2e} | "
              f"{diag['rmse_total']:7.3f} {diag['rmse_p1']:7.3f} {diag['rmse_p2']:7.3f} "
              f"{diag['rmse_p3']:7.3f} {diag['overshoot_100h']:+8.3f} {diag['drop_error']:+9.3f} "
              f"{diag['final_err']:+7.3f}")
    
    # ── MD SEARCH ──
    print("\n--- Munson-Dawson search ---")
    print(f"{'K0':>10} {'alpha_w':>8} | {'RMSE':>7} {'P1':>7} {'P2':>7} {'P3':>7} {'OS@100h':>8} {'drop_err':>9} {'final':>7}")
    print("-" * 90)
    
    md_results = []
    
    for K0, aw in product(K0_GRID, ALPHA_W_GRID):
        e_tot = integrate_md(t_s, sigma_fine_Pa, K0, aw)
        diag = compute_diagnostics(t_model_h, e_tot, lab_t_h, lab_e_pct, lab_sig)
        md_results.append((K0, aw, diag))
    
    md_results.sort(key=lambda x: x[2]['rmse_total'])
    
    for K0, aw, diag in md_results[:20]:
        print(f"{K0:10.2e} {aw:8.1f} | "
              f"{diag['rmse_total']:7.3f} {diag['rmse_p1']:7.3f} {diag['rmse_p2']:7.3f} "
              f"{diag['rmse_p3']:7.3f} {diag['overshoot_100h']:+8.3f} {diag['drop_error']:+9.3f} "
              f"{diag['final_err']:+7.3f}")


if __name__ == "__main__":
    main()
