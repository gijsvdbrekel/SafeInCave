"""
Plot calibration results for new cyclic triaxial creep data.

Reads JSON output from optimize_newdata.py and creates comparison figures:
  - Top panel:    Stress schedule (sigma_1, sigma_3, sigma_diff)
  - Bottom panel: Axial strain vs time — lab + model predictions with components
  - Right column: Parameter boxes
  - Per-stage RMSE and MAPE annotations

When run standalone, merges calibration_newdata_sic.json and
calibration_newdata_md.json (if both exist) into a single combined plot.
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          USER CONFIGURATION                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

DATA_DIR = os.path.join(os.path.dirname(__file__), "output")
FIG_DIR  = os.path.join(os.path.dirname(__file__), "figures")

SHOW_COMPONENTS = True
DPI = 180
SHOW = False

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                      END OF USER CONFIGURATION                             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


def _fmt(val):
    """Format a number compactly for parameter boxes."""
    if isinstance(val, float):
        if abs(val) < 1e-3 or abs(val) >= 1e4:
            return f"{val:.2e}"
        if val == int(val):
            return f"{val:.1f}"
        return f"{val:.4g}"
    return str(val)


def compute_stage_metrics(t_lab, e_lab, t_mod, e_mod, stages):
    """Compute per-stage RMSE and MAPE between lab and model."""
    e_mod_interp = np.interp(t_lab, t_mod, e_mod)
    metrics = []

    for stg in stages:
        t0, t1 = stg['t_start'], stg['t_end']
        mask = (t_lab >= t0) & (t_lab <= t1)
        if np.sum(mask) < 5:
            metrics.append(None)
            continue

        res = e_mod_interp[mask] - e_lab[mask]
        rmse = np.sqrt(np.mean(res ** 2))

        abs_lab = np.maximum(np.abs(e_lab[mask]), 0.005)
        mape = np.mean(np.abs(res) / abs_lab) * 100.0

        metrics.append({'rmse': rmse, 'mape': mape})

    return metrics


def detect_stages_from_stress(time_h, sigma_diff):
    """Simple stage detection for annotation on plots."""
    stages = []
    LOADING_THRESH = 5.0
    min_change = 5.0
    min_hold = 15.0

    n_init = min(100, len(sigma_diff))
    current_sigma = np.median(sigma_diff[:n_init])
    stage_start = time_h[0]

    for i in range(1, len(time_h)):
        if abs(sigma_diff[i] - current_sigma) > min_change:
            new_sigma = sigma_diff[i]
            hold_end = time_h[i] + min_hold
            held = True
            for j in range(i + 1, len(time_h)):
                if time_h[j] > hold_end:
                    break
                if abs(sigma_diff[j] - new_sigma) > min_change:
                    held = False
                    break
            else:
                if time_h[-1] - time_h[i] < 0.5 * min_hold:
                    held = False

            if held:
                stages.append({
                    't_start': stage_start,
                    't_end': time_h[i],
                    'sigma_mpa': float(current_sigma),
                })
                current_sigma = new_sigma
                stage_start = time_h[i]

    stages.append({
        't_start': stage_start,
        't_end': time_h[-1],
        'sigma_mpa': float(current_sigma),
    })

    # Use median of second half for sigma
    for stg in stages:
        t0, t1 = stg['t_start'], stg['t_end']
        t_mid = t0 + 0.5 * (t1 - t0)
        mask = (time_h >= t_mid) & (time_h <= t1)
        if np.any(mask):
            stg['sigma_mpa'] = float(np.median(sigma_diff[mask]))
        stg['type'] = "loading" if stg['sigma_mpa'] > LOADING_THRESH else "unloading"

    return stages


def plot_newdata(data, fig_dir):
    """Create a 2-panel figure (stress + axial strain) for the cyclic creep test."""
    test_id = data["test_id"]
    T_C = data["temperature_C"]
    s3_mean = data.get("sigma3_MPa_mean", 12.0)

    lab = data["lab"]
    t_lab = np.array(lab["time_hours"])
    e_ax_lab = np.array(lab["strain_axial_pct"])
    sig_diff_lab = np.array(lab["sigma_diff_MPa"])
    sig1_lab = np.array(lab.get("sigma1_MPa", sig_diff_lab + s3_mean))
    sig3_lab = np.array(lab.get("sigma3_MPa", np.full_like(sig_diff_lab, s3_mean)))

    has_sic = "safeincave" in data
    has_md = "munsondawson" in data

    # Detect stages for metrics
    stages = detect_stages_from_stress(t_lab, sig_diff_lab)

    # ── Figure layout: 2 rows (stress + axial strain) + parameter column ──
    fig = plt.figure(figsize=(20, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=[3, 1],
                          height_ratios=[1, 3],
                          hspace=0.08, wspace=0.2)

    ax_sig = fig.add_subplot(gs[0, 0])
    ax_ax = fig.add_subplot(gs[1, 0], sharex=ax_sig)
    ax_par = fig.add_subplot(gs[:, 1])
    ax_par.axis("off")

    # ── Top panel: stress schedule ───────────────────────────────────────────
    t_days = t_lab / 24
    ax_sig.plot(t_days, sig1_lab, "r-", linewidth=1.5, label="$\\sigma_1$")
    ax_sig.plot(t_days, sig3_lab, "b-", linewidth=1.5, label="$\\sigma_3$")
    ax_sig.plot(t_days, sig_diff_lab, "k--", linewidth=1, alpha=0.7,
                label="$\\sigma_{diff}$")
    ax_sig.set_ylabel("Stress (MPa)")
    ax_sig.set_title(f"Cyclic Triaxial Creep Test  |  T = {T_C} C  |  "
                     f"$\\sigma_3$ ~ {s3_mean:.0f} MPa")
    ax_sig.grid(True, alpha=0.3)
    ax_sig.legend(loc="upper right", fontsize=9)
    plt.setp(ax_sig.get_xticklabels(), visible=False)

    # ── Bottom panel: axial strain ───────────────────────────────────────────
    ax_ax.plot(t_days, e_ax_lab, "ko", markersize=1.5, alpha=0.5,
               label="Lab data", zorder=5)

    metrics_text = []

    if has_sic:
        sic = data["safeincave"]
        t_sic = np.array(sic["time_hours"])
        eps_sic = np.array(sic["strain_axial_pct"])
        ax_ax.plot(t_sic / 24, eps_sic, "-", color="#1f77b4", linewidth=2,
                   label="SafeInCave (total)")

        if SHOW_COMPONENTS:
            eps_dc = np.array(sic["strain_disloc_pct"])
            eps_kv = np.array(sic["strain_kelvin_pct"])
            eps_de = np.array(sic["strain_desai_pct"])
            ax_ax.plot(t_sic / 24, eps_dc, "--", color="#1f77b4", linewidth=1,
                       alpha=0.5, label="  dislocation")
            ax_ax.plot(t_sic / 24, eps_kv, "-.", color="#1f77b4", linewidth=1,
                       alpha=0.5, label="  Kelvin")
            if np.max(np.abs(eps_de)) > 1e-4:
                ax_ax.plot(t_sic / 24, eps_de, ":", color="#1f77b4",
                           linewidth=1, alpha=0.5, label="  Desai")

        # Per-stage metrics for SIC
        stg_metrics = compute_stage_metrics(
            t_lab, e_ax_lab, t_sic, eps_sic, stages
        )
        metrics_text.append("SIC:")
        for k, (stg, m) in enumerate(zip(stages, stg_metrics)):
            if m and stg['type'] == 'loading':
                metrics_text.append(
                    f"  S{k} ({stg['sigma_mpa']:.0f}MPa): "
                    f"RMSE={m['rmse']:.3f}% MAPE={m['mape']:.0f}%"
                )

    if has_md:
        md = data["munsondawson"]
        t_md = np.array(md["time_hours"])
        eps_md = np.array(md["strain_axial_pct"])
        ax_ax.plot(t_md / 24, eps_md, "-", color="#d62728", linewidth=2,
                   label="Munson-Dawson (total)")

        if SHOW_COMPONENTS:
            eps_ss = np.array(md["strain_steady_pct"])
            eps_tr = np.array(md["strain_transient_pct"])
            ax_ax.plot(t_md / 24, eps_ss, "--", color="#d62728", linewidth=1,
                       alpha=0.5, label="  steady-state")
            ax_ax.plot(t_md / 24, eps_tr, ":", color="#d62728", linewidth=1,
                       alpha=0.5, label="  transient")

        # Per-stage metrics for MD
        stg_metrics_md = compute_stage_metrics(
            t_lab, e_ax_lab, t_md, eps_md, stages
        )
        metrics_text.append("MD:")
        for k, (stg, m) in enumerate(zip(stages, stg_metrics_md)):
            if m and stg['type'] == 'loading':
                metrics_text.append(
                    f"  S{k} ({stg['sigma_mpa']:.0f}MPa): "
                    f"RMSE={m['rmse']:.3f}% MAPE={m['mape']:.0f}%"
                )

    ax_ax.set_xlabel("Time (days)")
    ax_ax.set_ylabel("Axial strain (%)")
    ax_ax.grid(True, alpha=0.3)
    ax_ax.legend(loc="upper left", fontsize=8, ncol=2)

    # Add metrics annotation
    if metrics_text:
        ax_ax.annotate("\n".join(metrics_text),
                       xy=(0.98, 0.02), xycoords="axes fraction",
                       ha="right", va="bottom", fontsize=7,
                       family="monospace",
                       bbox=dict(boxstyle="round,pad=0.3", fc="white",
                                 alpha=0.8, ec="0.7"))

    # ── Right column: parameter boxes ────────────────────────────────────────
    y = 0.98
    fs = 7.5
    lh = 0.018

    def _box(ax, x, y, title, lines, fc):
        header = f"  {title}\n"
        body = "\n".join(f"  {l}" for l in lines)
        text = header + body
        n_lines = 1 + len(lines)
        box_h = n_lines * lh + 0.02
        ax.annotate(text, xy=(x, y), xycoords="axes fraction",
                    va="top", ha="left", fontsize=fs, family="monospace",
                    bbox=dict(boxstyle="round,pad=0.4", fc=fc, ec="0.6",
                              alpha=0.85, linewidth=0.8))
        return y - box_h - 0.01

    # Elastic
    gp = data.get("params", {})
    el_lines = [
        f"E   = {gp.get('E_elastic', 0)/1e9:.3f} GPa",
        f"nu  = {gp.get('nu_elastic', 0):.2f}",
    ]
    y = _box(ax_par, 0.02, y, "ELASTIC (fixed)", el_lines, "#e8e8e8")

    # SafeInCave parameters
    if has_sic:
        sic = data["safeincave"]
        p = sic.get("params", {})
        sic_lines = [
            f"A   = {_fmt(p.get('A_mpa_yr', 0))} MPa^-n/yr",
            f"n   = {p.get('n', 0):.4f}",
            f"Q/R = {gp.get('Q_over_R', 6252):.0f} K",
            "-- Kelvin --",
            f"eta = {_fmt(p.get('eta_kelvin', 0))} Pa.s",
            f"E1  = {_fmt(p.get('E1_kelvin', 0))} Pa",
            "-- Desai --",
            f"mu1 = {_fmt(p.get('mu1_desai', 0))} 1/s",
            f"N1  = {p.get('N1_desai', 0):.1f}",
            f"a1  = {_fmt(p.get('a1_desai', 0))}",
            f"eta = {p.get('eta_desai', 0):.4g}",
            f"a0  = {p.get('alpha0_desai', 0):.4g}",
        ]
        y = _box(ax_par, 0.02, y, "SAFEINCAVE", sic_lines, "#dde8f7")

    # Munson-Dawson parameters
    if has_md:
        md = data["munsondawson"]
        p = md.get("params", {})
        md_lines = [
            f"A   = {_fmt(p.get('A_mpa_yr', 0))} MPa^-n/yr",
            f"n   = {p.get('n', 0):.4f}",
            f"Q/R = {gp.get('Q_over_R', 6252):.0f} K",
            "-- Transient --",
            f"K0  = {_fmt(p.get('K0', 0))}",
            f"c   = {p.get('c', 0):.5f} 1/K",
            f"m   = {p.get('m', 0):.4f}",
            f"a_w = {_fmt(p.get('alpha_w', 0))}",
            f"b_w = {_fmt(p.get('beta_w', 0))}",
            f"d   = {_fmt(p.get('delta', 0))}",
        ]
        y = _box(ax_par, 0.02, y, "MUNSON-DAWSON", md_lines, "#f7dddd")

    fig.subplots_adjust(left=0.06, right=0.97, top=0.95, bottom=0.07)

    outpath = os.path.join(fig_dir, f"calibration_{test_id}.png")
    fig.savefig(outpath, dpi=DPI)
    print(f"[SAVED] {outpath}")

    if SHOW:
        plt.show()
    plt.close(fig)


def load_and_merge():
    """Load SIC and MD JSON files and merge into a single dict."""
    sic_path = os.path.join(DATA_DIR, "calibration_newdata_sic.json")
    md_path = os.path.join(DATA_DIR, "calibration_newdata_md.json")

    data = None

    if os.path.exists(sic_path):
        with open(sic_path, "r") as f:
            data = json.load(f)
        print(f"  Loaded SIC: {sic_path}")

    if os.path.exists(md_path):
        with open(md_path, "r") as f:
            md_data = json.load(f)
        if data is None:
            data = md_data
        else:
            # Merge MD model into SIC data
            if "munsondawson" in md_data:
                data["munsondawson"] = md_data["munsondawson"]
        print(f"  Loaded MD:  {md_path}")

    if data is None:
        # Fallback: try old combined file
        old_path = os.path.join(DATA_DIR, "calibration_newdata.json")
        if os.path.exists(old_path):
            with open(old_path, "r") as f:
                data = json.load(f)
            print(f"  Loaded: {old_path}")

    return data


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    data = load_and_merge()
    if data is None:
        print("[ERROR] No calibration JSON files found in output/")
        print("        Run optimize_newdata.py first.")
        return

    models = []
    if "safeincave" in data:
        models.append("SIC")
    if "munsondawson" in data:
        models.append("MD")
    print(f"Plotting {data['test_id']} ({' + '.join(models)}) ...")

    plot_newdata(data, FIG_DIR)
    print("\n[DONE]")


if __name__ == "__main__":
    main()
