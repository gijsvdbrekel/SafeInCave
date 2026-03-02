"""
Quick comparison: MD calibration with constrained delta ranges.
Tests delta in [0.001, 2.0] and [0.001, 0.75] vs the unconstrained result (d=300).
"""
import os, sys, math, numpy as np, numba
from scipy.optimize import differential_evolution

sys.path.insert(0, os.path.dirname(__file__))
from calibrate_newdata import (
    read_excel_data, build_stress_schedule,
    R_GAS, E_ELASTIC, NU_ELASTIC, MU_SHEAR,
    C_MD, HOUR, _MD_MAX_DZETA, XLSX_PATH,
)
from optimize_newdata import fit_dislocation_params

_R = R_GAS
_E_EL = E_ELASTIC
_MU_SH = MU_SHEAR
_MAX_DZETA = _MD_MAX_DZETA
_MAX_SUBSTEPS = 5000
_SEC_PER_YR = 365.25 * 24 * 3600

T_CELSIUS = 21.0
T_KELVIN = T_CELSIUS + 273.15
DT_HOURS = 0.5
N_PER_STAGE = 40

W_LOADING_TRANSIENT  = 5.0
W_LOADING_STEADY     = 1.5
W_UNLOADING_TRANSIENT = 3.0
W_UNLOADING_STEADY   = 0.5
T_LOADING_TRANS  = 24.0
T_UNLOADING_TRANS = 12.0

QR_GRID = [6252.0, 6313.0, 6374.0, 6434.0, 6495.0]
A_BOUNDS = (15.0, 60.0)
N_BOUNDS = (3.0, 6.0)
STEADY_STATE_TAIL_FRAC = 0.3


@numba.njit(cache=True)
def integrate_md_fast(t_s, sigma_diff_Pa, A_pa_s, n, Q, T,
                      K0, c, m_md, alpha_w, beta_w, delta):
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


def main():
    # Read and prepare data
    time_h_raw, sigma1_raw, sigma3_raw, eps1_raw, eps3_raw = read_excel_data(XLSX_PATH)
    t_s, sig1_fine, sig3_fine, idx0 = build_stress_schedule(
        time_h_raw, sigma1_raw, sigma3_raw, DT_HOURS)

    sigma_diff_fine = (sig1_fine - sig3_fine) * 1e6
    t_mod_h = t_s / 3600.0
    t0_h = time_h_raw[idx0]
    t_lab = time_h_raw[idx0:] - t0_h
    e_lab_ax = eps1_raw[idx0:] - eps1_raw[idx0]

    # Detect stages
    sig_diff_lab = (sigma1_raw - sigma3_raw)[idx0:]
    stages = []
    diffs = np.abs(np.diff(sig_diff_lab))
    change_pts = np.where(diffs > 0.5)[0]
    boundaries = np.concatenate(([0], change_pts + 1, [len(sig_diff_lab)]))
    for j in range(len(boundaries) - 1):
        s, e = boundaries[j], boundaries[j + 1]
        if e - s < 3:
            continue
        smpa = np.median(sig_diff_lab[s:e])
        stages.append({
            'start': s, 'end': e, 'sigma_mpa': smpa,
            'type': 'loading' if smpa > 0.5 else 'unloading'
        })

    # Resample
    t_resamp, e_resamp, stage_ids, weights_arr = [], [], [], []
    for k, stg in enumerate(stages):
        s, e = stg['start'], stg['end']
        idx_stg = np.linspace(s, e - 1, N_PER_STAGE).astype(int)
        for ii in idx_stg:
            t_resamp.append(t_lab[ii])
            e_resamp.append(e_lab_ax[ii])
            stage_ids.append(k)
            dt_from_start = t_lab[ii] - t_lab[s]
            if stg['type'] == 'loading':
                w = W_LOADING_TRANSIENT if dt_from_start < T_LOADING_TRANS else W_LOADING_STEADY
            else:
                w = W_UNLOADING_TRANSIENT if dt_from_start < T_UNLOADING_TRANS else W_UNLOADING_STEADY
            weights_arr.append(w)

    t_resamp = np.array(t_resamp)
    e_resamp = np.array(e_resamp)
    weights = np.array(weights_arr)
    weights /= weights.sum()

    # Phase 1: fit dislocation
    stress_mpa = np.array([s['sigma_mpa'] for s in stages
                           if s['type'] == 'loading' and s['sigma_mpa'] > 1.0])
    rates = []
    for stg in stages:
        if stg['type'] != 'loading' or stg['sigma_mpa'] <= 1.0:
            continue
        s, e = stg['start'], stg['end']
        n_tail = max(2, int((e - s) * STEADY_STATE_TAIL_FRAC))
        t_tail = t_lab[e - n_tail:e] * 3600
        e_tail = e_lab_ax[e - n_tail:e] / 100
        sl = np.polyfit(t_tail, e_tail, 1)
        rates.append(max(sl[0], 1e-20))
    rates = np.array(rates)

    best_params, _ = fit_dislocation_params(stress_mpa, rates, T_KELVIN, QR_GRID)
    A_DISLOC_MPA_YR, N_DISLOC, QR_DISLOC = best_params
    Q_DISLOC_val = QR_DISLOC * R_GAS
    A_DISLOC_PA_S = A_DISLOC_MPA_YR * (1e-6) ** N_DISLOC / _SEC_PER_YR
    print(f"Phase 1: A={A_DISLOC_MPA_YR:.2f}, n={N_DISLOC:.4f}, Q/R={QR_DISLOC:.0f}")

    # Warm up numba
    _ = integrate_md_fast(t_s[:10], sigma_diff_fine[:10],
                          A_DISLOC_PA_S, N_DISLOC, Q_DISLOC_val, T_KELVIN,
                          1.0, C_MD, 1.5, -14.0, -8.0, 0.58)

    # Run for each delta constraint
    for label, d_bound in [
        ("delta [0.001, 300.0] (unconstrained)", (0.001, 300.0)),
        ("delta [0.001, 2.0]", (0.001, 2.0)),
        ("delta [0.001, 0.75]", (0.001, 0.75)),
    ]:
        print(f"\n{'='*70}")
        print(f"Optimizing: {label}")
        print(f"{'='*70}")

        BOUNDS = [
            (-6.0, 8.0),
            (0.1, 6.0),
            (-10.0, 200.0),
            (-30.0, 60.0),
            d_bound,
        ]

        _call_count = [0]
        _best = [np.inf]

        def make_objective(bounds_local):
            def objective(x):
                try:
                    log_K0, m_md, alpha_w, beta_w, delta = x
                    K0 = 10.0 ** log_K0
                    e_mod = integrate_md_fast(
                        t_s, sigma_diff_fine,
                        A_DISLOC_PA_S, N_DISLOC, Q_DISLOC_val, T_KELVIN,
                        K0, C_MD, m_md, alpha_w, beta_w, delta
                    )
                except Exception:
                    return 1e6

                e_mod_resamp = np.interp(t_resamp, t_mod_h, e_mod)
                eps_floor = 0.005
                abs_lab = np.maximum(np.abs(e_resamp), eps_floor)
                ape = np.abs(e_mod_resamp - e_resamp) / abs_lab
                wmape = np.sum(weights * ape) * 100.0

                _call_count[0] += 1
                if wmape < _best[0]:
                    _best[0] = wmape
                    log_K0v, m_mdv, alpha_wv, beta_wv, deltav = x
                    if _call_count[0] <= 5 or _call_count[0] % 200 == 0:
                        print(f"  [{_call_count[0]:5d}] WMAPE={wmape:.2f}% "
                              f"K0={10**log_K0v:.4f} m={m_mdv:.3f} "
                              f"aw={alpha_wv:.2f} bw={beta_wv:.2f} d={deltav:.4f}",
                              flush=True)
                elif _call_count[0] % 2000 == 0:
                    print(f"  [{_call_count[0]:5d}] best={_best[0]:.2f}%", flush=True)
                return wmape
            return objective

        rng = np.random.default_rng(42)
        _n_params = 5
        _popsize = 20
        _n_pop = _popsize * _n_params
        _init_pop = rng.uniform(0, 1, size=(_n_pop, _n_params))
        for j, (lo, hi) in enumerate(BOUNDS):
            _init_pop[:, j] = lo + _init_pop[:, j] * (hi - lo)
        _init_pop[0] = np.array([3.39, 2.47, 93.3, 30.0,
                                 min(d_bound[1], 0.58)])

        result = differential_evolution(
            make_objective(BOUNDS),
            bounds=BOUNDS,
            popsize=_popsize,
            maxiter=300,
            tol=1e-6,
            mutation=(0.5, 1.5),
            recombination=0.8,
            seed=42,
            init=_init_pop,
            polish=True,
            disp=False,
        )

        log_K0, m_md, alpha_w, beta_w, delta = result.x
        K0 = 10.0 ** log_K0
        print(f"\n  RESULT ({label}):")
        print(f"    K0      = {K0:.6f}")
        print(f"    m       = {m_md:.4f}")
        print(f"    alpha_w = {alpha_w:.4f}")
        print(f"    beta_w  = {beta_w:.4f}")
        print(f"    delta   = {delta:.6f}")
        print(f"    WMAPE   = {result.fun:.2f}%")


if __name__ == "__main__":
    main()
