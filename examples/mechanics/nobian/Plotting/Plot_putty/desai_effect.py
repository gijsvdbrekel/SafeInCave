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
    read_cell_scalar_timeseries,
    read_ksp_jsonl,
    DAY,
    MPA,
)

# =============================================================================
# STYLING (consistent with other plot scripts)
# =============================================================================
CAVERN_COLORS = {
    "Asymmetric":          "#1f77b4",   # blue
    "Direct-circulation":  "#ff7f0e",   # orange
    "IrregularFine":       "#d62728",   # red
    "Regular":             "#9467bd",   # purple
    "Reversed-circulation":"#8c564b",   # brown
    "Tilt":                "#e377c2",   # pink
    "Fast-leached":        "#17becf",   # cyan
    "Tube-failure":        "#bcbd22",   # olive
}

SCENARIO_COLORS = {
    "full":             "#2ca02c",   # green
    "full_minus_desai": "#d62728",   # red
    "full_md":          "#1f77b4",   # blue
}

# =============================================================================
# DEFAULT CONFIG (can be overridden via command line)
# =============================================================================
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", "Simulation", "output"))

DEFAULT_SELECT_BASE = {
    "caverns": ["Regular"],
    "pressure": "sinus",
    "n_cycles": 8,
    "operation_days": 365,
    "case_contains": None,
}

SCEN_FULL = "full"
SCEN_NODESAI = "full_md"

DEFAULT_DPI = 180
DEFAULT_SHOW = True


# =============================================================================
# HELPERS
# =============================================================================
def series_mean(times, vals_list):
    """Compute mean of cell values at each timestep."""
    y = np.array([np.nanmean(v) for v in vals_list], dtype=float)
    return np.asarray(times, float), y


def pick_case_for_scenario(all_cases_meta, select_base: dict, scenario_name: str):
    """Find a matching case for the given scenario."""
    sel = dict(select_base)
    sel["scenario"] = scenario_name

    candidates = filter_cases(all_cases_meta, sel)

    # Only require q_elems - other fields are optional for flexible comparison
    needed_fields = ["q_elems"]
    keep = []
    for m in candidates:
        cp = m["case_path"]
        ok = os.path.isfile(os.path.join(cp, "ksp_operation.jsonl"))
        ok = ok and all(os.path.isfile(path_field_xdmf(cp, f)) for f in needed_fields)
        if ok:
            keep.append(m)

    if not keep:
        print(f"\n[ERROR] No usable case for scenario='{scenario_name}' with select_base={select_base}")
        debug_print_inventory(all_cases_meta)
        raise RuntimeError("Selection matched nothing. Fix ROOT/SELECT_BASE or case naming.")

    chosen = newest_case(keep)
    return chosen["case_path"], chosen


def plot_desai_effect(case_full: str, case_nodesai: str, meta_full: dict, dpi: int):
    """Generate comparison plots between two scenarios."""
    tag = make_case_tag(meta_full)
    cavern_label = meta_full.get("cavern_label", "Unknown")
    base_color = CAVERN_COLORS.get(cavern_label, "#333333")

    # Get scenario labels from the actual scenarios being compared
    label_A = SCEN_FULL.upper()
    label_B = SCEN_NODESAI.upper()
    color_A = SCENARIO_COLORS.get(SCEN_FULL, "#2ca02c")
    color_B = SCENARIO_COLORS.get(SCEN_NODESAI, "#d62728")

    # KSP iterations comparison
    tF, itsF, _, _ = read_ksp_jsonl(case_full)
    tN, itsN, _, _ = read_ksp_jsonl(case_nodesai)

    fig1, ax = plt.subplots(1, 1, figsize=(11, 4))
    ax.set_title(f"KSP iterations vs time (operation) - {tag}")
    ax.plot(tF / DAY, itsF, linewidth=1.8, color=color_A, label=label_A)
    ax.plot(tN / DAY, itsN, linewidth=1.8, color=color_B, label=label_B)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("KSP iterations")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    # Mean von Mises stress comparison
    tqF, qF_list = read_cell_scalar_timeseries(path_field_xdmf(case_full, "q_elems"))
    tqN, qN_list = read_cell_scalar_timeseries(path_field_xdmf(case_nodesai, "q_elems"))
    tqF, qF = series_mean(tqF, qF_list)
    tqN, qN = series_mean(tqN, qN_list)

    fig2, ax2 = plt.subplots(1, 1, figsize=(11, 4))
    ax2.set_title(f"Mean q_elems vs time - {tag}")
    ax2.plot(tqF / DAY, qF / MPA, linewidth=1.8, color=color_A, label=label_A)
    ax2.plot(tqN / DAY, qN / MPA, linewidth=1.8, color=color_B, label=label_B)
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Mean von Mises (MPa)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="best")

    # Strain rates comparison (log scale) - try total strain rate if available
    fig3, fig4 = None, None
    try:
        ttF, rtotF_list = read_cell_scalar_timeseries(path_field_xdmf(case_full, "eps_ne_rate_eq_total"))
        ttN, rtotN_list = read_cell_scalar_timeseries(path_field_xdmf(case_nodesai, "eps_ne_rate_eq_total"))
        ttF, rtotF = series_mean(ttF, rtotF_list)
        ttN, rtotN = series_mean(ttN, rtotN_list)

        fig3, ax3 = plt.subplots(1, 1, figsize=(11, 4))
        ax3.set_title(f"Mean eq strain-rates (log scale) - {tag}")
        ax3.semilogy(ttF / DAY, np.maximum(rtotF, 1e-30), linewidth=1.8, color=color_A, label=f"{label_A} total")
        ax3.semilogy(ttN / DAY, np.maximum(rtotN, 1e-30), linewidth=1.8, color=color_B, label=f"{label_B} total")

        # Try to add Desai component if available
        try:
            tdF, rdesF_list = read_cell_scalar_timeseries(path_field_xdmf(case_full, "eps_ne_rate_eq_desai"))
            tdF, rdesF = series_mean(tdF, rdesF_list)
            ax3.semilogy(tdF / DAY, np.maximum(rdesF, 1e-30), linewidth=1.8, color=color_A, linestyle="--", label=f"{label_A} Desai")
        except Exception:
            pass

        try:
            tdN, rdesN_list = read_cell_scalar_timeseries(path_field_xdmf(case_nodesai, "eps_ne_rate_eq_desai"))
            tdN, rdesN = series_mean(tdN, rdesN_list)
            ax3.semilogy(tdN / DAY, np.maximum(rdesN, 1e-30), linewidth=1.8, color=color_B, linestyle="--", label=f"{label_B} Desai")
        except Exception:
            pass

        ax3.set_xlabel("Time (days)")
        ax3.set_ylabel("Mean eq strain-rate (1/s)")
        ax3.grid(True, alpha=0.3)
        ax3.legend(loc="best")

        # Ratio plot only if Desai fields exist for both
        try:
            tdF, rdesF_list = read_cell_scalar_timeseries(path_field_xdmf(case_full, "eps_ne_rate_eq_desai"))
            tdN, rdesN_list = read_cell_scalar_timeseries(path_field_xdmf(case_nodesai, "eps_ne_rate_eq_desai"))
            tdF, rdesF = series_mean(tdF, rdesF_list)
            tdN, rdesN = series_mean(tdN, rdesN_list)

            fig4, ax4 = plt.subplots(1, 1, figsize=(11, 4))
            ax4.set_title(f"Desai contribution ratio (Desai / Total) - {tag}")
            ratio_F = np.where(rtotF > 1e-30, rdesF / rtotF, 0.0)
            ratio_N = np.where(rtotN > 1e-30, rdesN / rtotN, 0.0)
            ax4.plot(ttF / DAY, ratio_F * 100, linewidth=1.8, color=color_A, label=label_A)
            ax4.plot(ttN / DAY, ratio_N * 100, linewidth=1.8, color=color_B, label=label_B)
            ax4.set_xlabel("Time (days)")
            ax4.set_ylabel("Desai contribution (%)")
            ax4.set_ylim(0, 105)
            ax4.grid(True, alpha=0.3)
            ax4.legend(loc="best")
        except Exception:
            print("[INFO] Desai contribution ratio plot skipped (Desai fields not available)")

    except Exception as e:
        print(f"[INFO] Strain rate plots skipped: {e}")

    return fig1, fig2, fig3, fig4, tag


def print_config_summary(args, select_base):
    """Print configuration summary at startup."""
    print("=" * 60)
    print("DESAI EFFECT COMPARISON")
    print("=" * 60)
    print(f"  ROOT:           {args.root}")
    print(f"  Scenario A:     {SCEN_FULL}")
    print(f"  Scenario B:     {SCEN_NODESAI}")
    print(f"  Caverns:        {select_base.get('caverns')}")
    print(f"  Pressure:       {select_base.get('pressure')}")
    print(f"  N cycles:       {select_base.get('n_cycles')}")
    print(f"  Operation days: {select_base.get('operation_days')}")
    print(f"  DPI:            {args.dpi}")
    print(f"  Show plots:     {args.show}")
    print("=" * 60)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare simulations with and without Desai creep component",
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

    case_full, meta_full = pick_case_for_scenario(all_cases, select_base, SCEN_FULL)
    case_nodesai, _ = pick_case_for_scenario(all_cases, select_base, SCEN_NODESAI)

    print(f"\n[INFO] {SCEN_FULL} case:    {case_full}")
    print(f"[INFO] {SCEN_NODESAI} case:  {case_nodesai}")

    fig1, fig2, fig3, fig4, tag = plot_desai_effect(case_full, case_nodesai, meta_full, args.dpi)

    p1 = os.path.join(out_dir, f"scenario_compare_ksp_{tag}.png")
    p2 = os.path.join(out_dir, f"scenario_compare_qmean_{tag}.png")
    p3 = os.path.join(out_dir, f"scenario_compare_rates_{tag}.png")
    p4 = os.path.join(out_dir, f"scenario_compare_ratio_{tag}.png")

    fig1.savefig(p1, dpi=args.dpi, bbox_inches="tight")
    fig2.savefig(p2, dpi=args.dpi, bbox_inches="tight")
    print(f"\n[SAVED] {p1}")
    print(f"[SAVED] {p2}")

    if fig3 is not None:
        fig3.savefig(p3, dpi=args.dpi, bbox_inches="tight")
        print(f"[SAVED] {p3}")
    if fig4 is not None:
        fig4.savefig(p4, dpi=args.dpi, bbox_inches="tight")
        print(f"[SAVED] {p4}")

    if args.show:
        plt.show()

    plt.close(fig1)
    plt.close(fig2)
    if fig3 is not None:
        plt.close(fig3)
    if fig4 is not None:
        plt.close(fig4)


if __name__ == "__main__":
    main()
