import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

from case_index import (
    detect_layout_and_collect_cases,
    filter_cases,
    newest_case,
    make_case_tag,
    debug_print_inventory,
    path_field_xdmf,
    path_pressure_json,
    read_cell_scalar_timeseries,
    read_pressure_schedule,
    DAY,
    MPA,
)

# =============================================================================
# STYLING (consistent with other plot scripts)
# =============================================================================
CAVERN_COLORS = {
    "Asymmetric":    "#1f77b4",   # blue
    "Irregular":     "#ff7f0e",   # orange
    "IrregularFine": "#d62728",   # red
    "Multichamber":  "#2ca02c",   # green
    "Regular":       "#9467bd",   # purple
    "Teardrop":      "#8c564b",   # brown
    "Tilt":          "#e377c2",   # pink
}

SCENARIO_COLORS = {
    "disloc_old_only": "#1f77b4",   # blue
    "disloc_new_only": "#ff7f0e",   # orange
}

# =============================================================================
# DEFAULT CONFIG (can be overridden via command line)
# =============================================================================
DEFAULT_ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/3_cavern/output"

DEFAULT_SELECT_BASE = {
    "caverns": ["Regular"],
    "pressure": "sinus",
    "n_cycles": 8,
    "operation_days": 365,
    "case_contains": None,
}

SCEN_OLD = "disloc_old_only"
SCEN_NEW = "disloc_new_only"

DEFAULT_PICK_PEAK = "last"  # "first", "last", or "cycle_N" where N is cycle number
DEFAULT_DPI = 180
DEFAULT_SHOW = False


# =============================================================================
# HELPERS
# =============================================================================
def pick_peak_timestep_index(t_s, p_mpa, pick_peak: str = "last"):
    """
    Find timestep index at pressure peak.

    Parameters
    ----------
    t_s : np.ndarray
        Time in seconds
    p_mpa : np.ndarray
        Pressure in MPa
    pick_peak : str
        "first" for first peak, "last" for last peak,
        or "cycle_N" for Nth cycle peak (1-indexed)
    """
    p = np.asarray(p_mpa)
    peaks = np.where((p[1:-1] > p[:-2]) & (p[1:-1] > p[2:]))[0] + 1

    if peaks.size == 0:
        return int(np.argmax(p))

    if pick_peak == "first":
        return int(peaks[0])
    elif pick_peak == "last":
        return int(peaks[-1])
    elif pick_peak.startswith("cycle_"):
        try:
            cycle_num = int(pick_peak.split("_")[1]) - 1  # Convert to 0-indexed
            if 0 <= cycle_num < len(peaks):
                return int(peaks[cycle_num])
            else:
                print(f"[WARNING] Requested cycle {cycle_num + 1} but only {len(peaks)} peaks found. Using last.")
                return int(peaks[-1])
        except ValueError:
            print(f"[WARNING] Invalid pick_peak '{pick_peak}'. Using last.")
            return int(peaks[-1])
    else:
        return int(peaks[-1])


def bin_median_trend(x, y, nbins: int = 35):
    """Compute binned median trend for scatter data."""
    x = np.asarray(x)
    y = np.asarray(y)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < 100:
        return None, None
    edges = np.linspace(np.min(x), np.max(x), nbins + 1)
    xm, ym = [], []
    for i in range(nbins):
        mi = (x >= edges[i]) & (x < edges[i + 1])
        if np.count_nonzero(mi) < 30:
            continue
        xm.append(0.5 * (edges[i] + edges[i + 1]))
        ym.append(np.median(y[mi]))
    return np.asarray(xm), np.asarray(ym)


def pick_case_for_scenario(all_cases_meta, select_base: dict, scenario_name: str):
    """Find a matching case for the given scenario."""
    sel = dict(select_base)
    sel["scenario"] = scenario_name
    candidates = filter_cases(all_cases_meta, sel)

    needed_fields = ["q_elems", "eps_ne_rate_eq_disloc"]
    keep = []
    for m in candidates:
        cp = m["case_path"]
        ok = os.path.isfile(path_pressure_json(cp))
        ok = ok and all(os.path.isfile(path_field_xdmf(cp, f)) for f in needed_fields)
        if ok:
            keep.append(m)

    if not keep:
        print(f"\n[ERROR] No usable case for scenario='{scenario_name}' with select_base={select_base}")
        debug_print_inventory(all_cases_meta)
        raise RuntimeError("Selection matched nothing. Fix SELECT_BASE or case naming.")

    chosen = newest_case(keep)
    return chosen["case_path"], chosen


def plot_old_vs_new(case_old: str, case_new: str, meta_old: dict, pick_peak: str, use_hexbin: bool = False):
    """Generate comparison plots for old vs new dislocation model."""
    tag = make_case_tag(meta_old)

    t_sched_s, p_sched_mpa = read_pressure_schedule(case_old)
    i_peak = pick_peak_timestep_index(t_sched_s, p_sched_mpa, pick_peak)
    target_t = t_sched_s[i_peak]

    # Read timeseries data
    q_old_t, q_old_list = read_cell_scalar_timeseries(path_field_xdmf(case_old, "q_elems"))
    r_old_t, r_old_list = read_cell_scalar_timeseries(path_field_xdmf(case_old, "eps_ne_rate_eq_disloc"))

    q_new_t, q_new_list = read_cell_scalar_timeseries(path_field_xdmf(case_new, "q_elems"))
    r_new_t, r_new_list = read_cell_scalar_timeseries(path_field_xdmf(case_new, "eps_ne_rate_eq_disloc"))

    # Find closest timestep to target
    k_old = int(np.argmin(np.abs(q_old_t - target_t)))
    k_new = int(np.argmin(np.abs(q_new_t - target_t)))

    q_old = np.asarray(q_old_list[k_old], float) / MPA
    r_old = np.asarray(r_old_list[k_old], float)
    q_new = np.asarray(q_new_list[k_new], float) / MPA
    r_new = np.asarray(r_new_list[k_new], float)

    # Log transform
    eps = 1e-30
    x_old = np.log10(np.maximum(q_old, eps))
    y_old = np.log10(np.maximum(r_old, eps))
    x_new = np.log10(np.maximum(q_new, eps))
    y_new = np.log10(np.maximum(r_new, eps))

    # Main scatter/hexbin plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.set_title(f"Dislocation: log10(strain-rate) vs log10(von Mises)\npeak='{pick_peak}' @ t={target_t / DAY:.2f} days - {tag}")

    if use_hexbin:
        # Use hexbin for dense data visualization
        hb_old = ax.hexbin(x_old, y_old, gridsize=40, cmap='Blues', mincnt=1, alpha=0.6, label="OLD dislocation")
        hb_new = ax.hexbin(x_new, y_new, gridsize=40, cmap='Oranges', mincnt=1, alpha=0.6, label="NEW dislocation")
    else:
        ax.scatter(x_old, y_old, s=4, alpha=0.25, color=SCENARIO_COLORS["disloc_old_only"], label="OLD dislocation")
        ax.scatter(x_new, y_new, s=4, alpha=0.25, color=SCENARIO_COLORS["disloc_new_only"], label="NEW dislocation")

    # Add median trend lines
    xm_old, ym_old = bin_median_trend(x_old, y_old)
    if xm_old is not None:
        ax.plot(xm_old, ym_old, linewidth=2.5, color=SCENARIO_COLORS["disloc_old_only"], label="OLD median trend")
    xm_new, ym_new = bin_median_trend(x_new, y_new)
    if xm_new is not None:
        ax.plot(xm_new, ym_new, linewidth=2.5, color=SCENARIO_COLORS["disloc_new_only"], label="NEW median trend")

    ax.set_xlabel("log10(q_elems [MPa])")
    ax.set_ylabel("log10(epsdot_disloc_eq [1/s])")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    # Pressure schedule plot showing selected peak
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 3))
    ax2.set_title("Pressure schedule (MPa) - selected peak")
    ax2.plot(t_sched_s / DAY, p_sched_mpa, linewidth=1.5, color="#333333")
    ax2.scatter([t_sched_s[i_peak] / DAY], [p_sched_mpa[i_peak]], s=80, color="#d62728", zorder=5, label=f"Selected: {pick_peak}")
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Pressure (MPa)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="best")

    # Additional plot: rate ratio comparison over time (mean values)
    fig3, ax3 = plt.subplots(1, 1, figsize=(10, 4))
    ax3.set_title(f"Mean dislocation strain-rate over time - {tag}")

    # Compute mean rates at each timestep
    r_old_mean = np.array([np.nanmean(r) for r in r_old_list])
    r_new_mean = np.array([np.nanmean(r) for r in r_new_list])

    ax3.semilogy(r_old_t / DAY, np.maximum(r_old_mean, 1e-30), linewidth=1.8,
                 color=SCENARIO_COLORS["disloc_old_only"], label="OLD dislocation")
    ax3.semilogy(r_new_t / DAY, np.maximum(r_new_mean, 1e-30), linewidth=1.8,
                 color=SCENARIO_COLORS["disloc_new_only"], label="NEW dislocation")
    ax3.axvline(target_t / DAY, color="#d62728", linestyle="--", alpha=0.7, label=f"Selected peak")
    ax3.set_xlabel("Time (days)")
    ax3.set_ylabel("Mean eq strain-rate (1/s)")
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc="best")

    return fig, fig2, fig3, tag


def print_config_summary(args, select_base):
    """Print configuration summary at startup."""
    print("=" * 60)
    print("DISLOCATION EFFECT COMPARISON (OLD vs NEW)")
    print("=" * 60)
    print(f"  ROOT:           {args.root}")
    print(f"  Scenario A:     {SCEN_OLD}")
    print(f"  Scenario B:     {SCEN_NEW}")
    print(f"  Peak selection: {args.pick_peak}")
    print(f"  Caverns:        {select_base.get('caverns')}")
    print(f"  Pressure:       {select_base.get('pressure')}")
    print(f"  N cycles:       {select_base.get('n_cycles')}")
    print(f"  Operation days: {select_base.get('operation_days')}")
    print(f"  Use hexbin:     {args.hexbin}")
    print(f"  DPI:            {args.dpi}")
    print(f"  Show plots:     {args.show}")
    print("=" * 60)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare old vs new dislocation creep models",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--root", type=str, default=DEFAULT_ROOT,
                        help="Output root directory containing case folders")
    parser.add_argument("--caverns", nargs="+", default=DEFAULT_SELECT_BASE["caverns"],
                        help="Cavern types to select")
    parser.add_argument("--pressure", type=str, default=DEFAULT_SELECT_BASE["pressure"],
                        help="Pressure scenario")
    parser.add_argument("--n-cycles", type=int, default=DEFAULT_SELECT_BASE["n_cycles"],
                        help="Number of cycles")
    parser.add_argument("--operation-days", type=int, default=DEFAULT_SELECT_BASE["operation_days"],
                        help="Operation duration in days")
    parser.add_argument("--case-contains", type=str, default=None,
                        help="Substring filter on case name")
    parser.add_argument("--pick-peak", type=str, default=DEFAULT_PICK_PEAK,
                        help="Peak selection: 'first', 'last', or 'cycle_N' (e.g., 'cycle_3')")
    parser.add_argument("--hexbin", action="store_true",
                        help="Use hexbin instead of scatter for dense data visualization")
    parser.add_argument("--show", action="store_true", default=DEFAULT_SHOW,
                        help="Show plots interactively")
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI,
                        help="Output image DPI")
    return parser.parse_args()


def main():
    args = parse_args()

    # Build selection criteria from arguments
    select_base = {
        "caverns": args.caverns,
        "pressure": args.pressure,
        "n_cycles": args.n_cycles,
        "operation_days": args.operation_days,
        "case_contains": args.case_contains,
    }

    print_config_summary(args, select_base)

    out_dir = os.path.join(args.root, "_figures")
    os.makedirs(out_dir, exist_ok=True)

    all_cases = detect_layout_and_collect_cases(args.root)
    if not all_cases:
        raise RuntimeError(f"No cases found under ROOT={args.root}")

    case_old, meta_old = pick_case_for_scenario(all_cases, select_base, SCEN_OLD)
    case_new, _ = pick_case_for_scenario(all_cases, select_base, SCEN_NEW)

    print(f"\n[INFO] OLD case: {case_old}")
    print(f"[INFO] NEW case: {case_new}")

    fig, fig2, fig3, tag = plot_old_vs_new(case_old, case_new, meta_old, args.pick_peak, args.hexbin)

    p1 = os.path.join(out_dir, f"disloc_effect_scatter_{tag}_peak={args.pick_peak}.png")
    p2 = os.path.join(out_dir, f"disloc_effect_pressure_{tag}_peak={args.pick_peak}.png")
    p3 = os.path.join(out_dir, f"disloc_effect_timeseries_{tag}.png")

    fig.savefig(p1, dpi=args.dpi, bbox_inches="tight")
    fig2.savefig(p2, dpi=args.dpi, bbox_inches="tight")
    fig3.savefig(p3, dpi=args.dpi, bbox_inches="tight")

    print(f"\n[SAVED] {p1}")
    print(f"[SAVED] {p2}")
    print(f"[SAVED] {p3}")

    if args.show:
        plt.show()

    plt.close(fig)
    plt.close(fig2)
    plt.close(fig3)


if __name__ == "__main__":
    main()
