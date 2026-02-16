"""
Plot calibration results: lab data vs. SafeInCave vs. Munson-Dawson.

Reads JSON output from run_calibration.py and creates comparison figures.
Each test gets its own figure with:
  - Top panel:    Stress schedule (differential stress over time)
  - Bottom panel: Strain vs time — lab data + model predictions (with components)
"""

import os
import json
import glob
import numpy as np
import matplotlib.pyplot as plt

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                          USER CONFIGURATION                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── INPUT / OUTPUT ───────────────────────────────────────────────────────────
DATA_DIR = os.path.join(os.path.dirname(__file__), "output")
FIG_DIR  = os.path.join(os.path.dirname(__file__), "figures")

# ── PLOT OPTIONS ─────────────────────────────────────────────────────────────
SHOW_COMPONENTS = True      # show disloc/kelvin/desai/transient breakdown
DPI = 180
SHOW = False

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                      END OF USER CONFIGURATION                             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


def _fmt(val, fmt=":.2e"):
    """Format a number compactly for parameter boxes."""
    if isinstance(val, float):
        if abs(val) < 1e-3 or abs(val) >= 1e4:
            return f"{val:.2e}"
        if val == int(val):
            return f"{val:.1f}"
        return f"{val:.4g}"
    return str(val)


def _c(calibrated):
    """Return ' (calibr)' tag if True."""
    return " (calibr)" if calibrated else ""


def plot_single_test(data, fig_dir):
    """Create a 2-panel figure for one creep test, with full parameter boxes."""
    test_id = data["test_id"]
    T_C = data["temperature_C"]
    s3 = data["sigma3_MPa"]

    lab = data["lab"]
    t_lab = np.array(lab["time_hours"])
    eps_lab = np.array(lab["strain_pct"])
    sig_lab = np.array(lab["sigma_diff_MPa"])

    has_sic = "safeincave" in data
    has_md = "munsondawson" in data

    # Wider figure to make room for parameter boxes on the right
    fig = plt.figure(figsize=(18, 9))

    # Left column: stress + strain panels.  Right column: parameter boxes.
    gs = fig.add_gridspec(2, 2, width_ratios=[2.8, 1],
                          height_ratios=[1, 2.5], hspace=0.08, wspace=0.25)
    ax_sig = fig.add_subplot(gs[0, 0])
    ax_eps = fig.add_subplot(gs[1, 0], sharex=ax_sig)
    ax_par = fig.add_subplot(gs[:, 1])  # full-height param panel
    ax_par.axis("off")

    # ── Top panel: stress schedule ──────────────────────────────────────────
    ax_sig.plot(t_lab / 24, sig_lab, "k-", linewidth=2, label="Lab $\\sigma_{diff}$")
    ax_sig.set_ylabel("Differential stress (MPa)")
    ax_sig.set_title(f"{test_id}  |  T = {T_C} °C  |  $\\sigma_3$ = {s3} MPa")
    ax_sig.grid(True, alpha=0.3)
    ax_sig.legend(loc="upper right", fontsize=9)

    # ── Bottom panel: strain vs time ────────────────────────────────────────
    ax_eps.plot(t_lab / 24, eps_lab, "ko", markersize=3, label="Lab data", zorder=5)

    if has_sic:
        sic = data["safeincave"]
        t_sic = np.array(sic["time_hours"])
        eps_sic = np.array(sic["strain_total_pct"])
        ax_eps.plot(t_sic / 24, eps_sic, "-", color="#1f77b4", linewidth=2,
                    label="SafeInCave (total)")

        if SHOW_COMPONENTS:
            eps_dc = np.array(sic["strain_disloc_pct"])
            eps_kv = np.array(sic["strain_kelvin_pct"])
            eps_de = np.array(sic["strain_desai_pct"])
            ax_eps.plot(t_sic / 24, eps_dc, "--", color="#1f77b4", linewidth=1,
                        alpha=0.6, label="  dislocation")
            ax_eps.plot(t_sic / 24, eps_kv, "-.", color="#1f77b4", linewidth=1,
                        alpha=0.6, label="  Kelvin")
            ax_eps.plot(t_sic / 24, eps_de, ":", color="#1f77b4", linewidth=1,
                        alpha=0.6, label="  Desai")

    if has_md:
        md = data["munsondawson"]
        t_md = np.array(md["time_hours"])
        eps_md = np.array(md["strain_total_pct"])
        ax_eps.plot(t_md / 24, eps_md, "-", color="#d62728", linewidth=2,
                    label="Munson-Dawson (total)")

        if SHOW_COMPONENTS:
            eps_ss = np.array(md["strain_steady_pct"])
            eps_tr = np.array(md["strain_transient_pct"])
            ax_eps.plot(t_md / 24, eps_ss + (eps_md - eps_ss - eps_tr),
                        "--", color="#d62728", linewidth=1, alpha=0.6,
                        label="  steady-state")
            ax_eps.plot(t_md / 24, eps_tr, ":", color="#d62728", linewidth=1,
                        alpha=0.6, label="  transient")

    ax_eps.set_xlabel("Time (days)")
    ax_eps.set_ylabel("Axial strain (%)")
    ax_eps.grid(True, alpha=0.3)
    ax_eps.legend(loc="lower right", fontsize=8, ncol=2)

    # ── Right column: parameter boxes ───────────────────────────────────────
    gp = data.get("params", {})
    y = 0.98   # vertical cursor (top-down in axes fraction)
    fs = 7.5   # font size
    lh = 0.024 # line height in axes fraction

    def _box(ax, x, y, title, lines, fc):
        """Draw a parameter box and return the new y position below it."""
        header = f"  {title}\n"
        body = "\n".join(f"  {l}" for l in lines)
        text = header + body
        n_lines = 1 + len(lines)
        box_h = n_lines * lh + 0.02
        ax.annotate(text, xy=(x, y), xycoords="axes fraction",
                    va="top", ha="left", fontsize=fs, family="monospace",
                    bbox=dict(boxstyle="round,pad=0.4", fc=fc, ec="0.6",
                              alpha=0.85, linewidth=0.8))
        return y - box_h - 0.015

    # -- Elastic parameters --
    el_lines = [
        f"E   = {gp.get('E_elastic', 0)/1e9:.3f} GPa",
        f"nu  = {gp.get('nu_elastic', 0):.2f}",
    ]
    y = _box(ax_par, 0.02, y, "ELASTIC", el_lines, "#e8e8e8")

    # -- Shared dislocation creep --
    A_display = gp.get("A_mpa_yr", None)
    n_display = gp.get("n_disloc", None)
    if A_display is None and has_sic:
        A_display = sic["params"]["A"]
    if n_display is None and has_sic:
        n_display = sic["params"]["n"]
    shared_lines = [
        f"A   = {_fmt(A_display)} MPa^-n/yr{_c(True)}",
        f"n   = {n_display:.2f}{_c(True)}",
        f"Q/R = {gp.get('Q_over_R', 6445.0):.1f} K{_c(True)}",
    ]
    y = _box(ax_par, 0.02, y, "DISLOCATION CREEP (shared)", shared_lines, "#f0f0d8")

    # -- SafeInCave --
    if has_sic:
        p = sic["params"]
        sic_lines = [
            "── Kelvin ──",
            f"eta   = {_fmt(p['eta_kelvin'])} Pa.s{_c(True)}",
            f"E1    = {_fmt(p['E1_kelvin'])} Pa{_c(True)}",
            f"nu1   = {p.get('nu1_kelvin', 0.25):.2f}{_c(True)}",
            "── Desai ──",
            f"mu1   = {_fmt(p['mu1_desai'])} 1/s{_c(True)}",
            f"N1    = {p['N1_desai']:.1f}{_c(True)}",
            f"a1    = {_fmt(p['a1_desai'])}{_c(True)}",
            f"eta   = {p['eta_desai']:.4g}{_c(True)}",
            f"alpha0= {p['alpha0_desai']:.4g}{_c(True)}",
            f"n     = {p.get('n_desai', 3.0):.1f}",
            f"beta1 = {p.get('beta1_desai', 0.0048):.4f}",
            f"beta  = {p.get('beta_desai', 0.995):.3f}",
            f"m     = {p.get('m_desai', -0.5):.1f}",
            f"gamma = {p.get('gamma_desai', 0.095):.3f}",
            f"sig_t = {p.get('sigma_t_desai', 5.0):.1f} Pa",
        ]
        y = _box(ax_par, 0.02, y, "SAFEINCAVE", sic_lines, "#dde8f7")

    # -- Munson-Dawson --
    if has_md:
        p = md["params"]
        md_lines = [
            "── Transient ──",
            f"K0      = {_fmt(p['K0'])}{_c(True)}",
            f"c       = {p['c']:.5f} 1/K",
            f"m       = {p['m']:.1f}",
            f"alpha_w = {p['alpha_w']:.2f}{_c(True)}",
            f"beta_w  = {p.get('beta_w', -7.738):.3f}",
            f"delta   = {p.get('delta', 0.58):.2f}",
        ]
        y = _box(ax_par, 0.02, y, "MUNSON-DAWSON", md_lines, "#f7dddd")

    fig.subplots_adjust(left=0.06, right=0.97, top=0.94, bottom=0.08)

    outpath = os.path.join(fig_dir, f"calibration_{test_id}.png")
    fig.savefig(outpath, dpi=DPI)
    print(f"[SAVED] {outpath}")

    if SHOW:
        plt.show()
    plt.close(fig)


def plot_overview(all_data, fig_dir):
    """
    Overview figure: all tests in a single grid.
    Shows only total strain curves (no component breakdown).
    """
    n = len(all_data)
    if n == 0:
        return

    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows),
                              squeeze=False)

    for idx, data in enumerate(all_data):
        ax = axes[idx // ncols][idx % ncols]
        test_id = data["test_id"]
        lab = data["lab"]
        t_lab = np.array(lab["time_hours"]) / 24
        eps_lab = np.array(lab["strain_pct"])

        ax.plot(t_lab, eps_lab, "ko", markersize=2, label="Lab")

        if "safeincave" in data:
            sic = data["safeincave"]
            t = np.array(sic["time_hours"]) / 24
            ax.plot(t, sic["strain_total_pct"], "-", color="#1f77b4",
                    linewidth=1.5, label="SafeInCave")

        if "munsondawson" in data:
            md = data["munsondawson"]
            t = np.array(md["time_hours"]) / 24
            ax.plot(t, md["strain_total_pct"], "-", color="#d62728",
                    linewidth=1.5, label="Munson-Dawson")

        ax.set_title(test_id, fontsize=10)
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Strain (%)")
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(fontsize=8)

    # Hide unused axes
    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle("Creep Test Calibration Overview", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    outpath = os.path.join(fig_dir, "calibration_overview.png")
    fig.savefig(outpath, dpi=DPI)
    print(f"[SAVED] {outpath}")

    if SHOW:
        plt.show()
    plt.close(fig)


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    json_files = sorted(glob.glob(os.path.join(DATA_DIR, "calibration_TCC*.json")))
    if not json_files:
        print(f"[ERROR] No calibration_TCC*.json files found in {DATA_DIR}")
        print("        Run run_calibration.py first.")
        return

    print(f"Found {len(json_files)} result file(s)")

    all_data = []
    for jf in json_files:
        with open(jf, "r") as f:
            data = json.load(f)
        print(f"  Plotting {data['test_id']}...")
        plot_single_test(data, FIG_DIR)
        all_data.append(data)

    if len(all_data) > 1:
        plot_overview(all_data, FIG_DIR)

    print("\n[DONE]")


if __name__ == "__main__":
    main()
