import os
import numpy as np
import matplotlib.pyplot as plt
import meshio

import safeincave.PostProcessingTools as post

from case_index import (
    detect_layout_and_collect_cases,
    filter_cases,
)


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                           USER CONFIGURATION                                  ║
# ╠══════════════════════════════════════════════════════════════════════════════╣
# ║  Modify the settings below to configure which cases to plot.                  ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── OUTPUT FOLDER ──────────────────────────────────────────────────────────────
# ROOT: Path to the output folder containing simulation results
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/nobian/simulation/output"

# ── CASE SELECTION ─────────────────────────────────────────────────────────────
# SELECT: Filter which cases to plot (set to None to include all)
#
# Available filters:
#   "caverns"        - Cavern shapes to include: ["Regular", "Irregular", ...] or None for all
#   "pressure"       - Pressure scenario: "sinus", "linear", "irregular", "csv_profile", or None
#   "scenario"       - Material scenario(s): ["disloc_old_only", "full", ...] or None
#   "n_cycles"       - Number of cycles (int) or None
#   "operation_days" - Operation duration (int) or None
#   "case_contains"  - Substring match in case name or None

SELECT = {
    "caverns": None,                                      # e.g. ["Regular", "Irregular"]
    "pressure": "sinus",                                  # "sinus"/"irregular"/"linear"/"csv_profile"
    "scenario": ["full", "full_minus_desai"],   # from ScenarioTest.py
    "n_cycles": None,
    "operation_days": None,
    "case_contains": None,
}

# ── PLOT OPTIONS ───────────────────────────────────────────────────────────────
# ONE_CASE_PER_SERIES: If True, pick only one case per (cavern, scenario) to avoid duplicates
ONE_CASE_PER_SERIES = True

# PLOT_MODE: How to arrange the plots
#   "combined"  - All cases in one figure (overlay stress paths)
#   "separate"  - One figure per case
PLOT_MODE = "combined"

# ── DILATANCY BOUNDARIES ───────────────────────────────────────────────────────
# Choose which dilatancy boundaries to show on stress path plots
# Available: "ratigan_027", "ratigan_018", "spiers", "devries_comp", "devries_ext"
SHOW_DILATANCY = ["ratigan_027", "spiers", "devries_comp", "devries_ext"]

# ── OUTPUT SETTINGS ────────────────────────────────────────────────────────────
OUT_DIR = os.path.join(ROOT, "_figures")
SHOW = False       # Set True to display plots interactively (requires GUI)
DPI = 180

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                        END OF USER CONFIGURATION                              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


# ══════════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════
MPA = 1e6
HOUR = 3600.0
DAY = 24.0 * HOUR

# Ordering for consistent legend/colors
CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt", "IrregularFine"]
SCENARIO_ORDER = ["disloc_old_only", "disloc_new_only", "desai_only", "full_minus_desai", "full",
                  "md_only", "md_steady_only", "full_md", "interlayer", "nointerlayer", None]

# Color and linestyle definitions
CAVERN_COLORS = {
    "Regular": "#1f77b4",       # blue
    "Irregular": "#ff7f0e",     # orange
    "Tilted": "#2ca02c",        # green
    "Tilt": "#2ca02c",          # green (alias)
    "Teardrop": "#d62728",      # red
    "Asymmetric": "#9467bd",    # purple
    "Multichamber": "#8c564b",  # brown
    "IrregularFine": "#e377c2", # pink
}

SCENARIO_LINESTYLES = {
    "disloc_old_only": "-",
    "disloc_new_only": "--",
    "desai_only": "-.",
    "full_minus_desai": "-",   # solid
    "full": "-",               # solid
    # Munson-Dawson scenarios
    "md_only":          "-",
    "md_steady_only":   "--",
    "full_md":          "-.",
    # Interlayer scenarios
    "interlayer":       "-",
    "nointerlayer":     "--",
    None: "-",
}

SCENARIO_COLORS = {
    "disloc_old_only":  "#1f77b4",   # blue
    "disloc_new_only":  "#ff7f0e",   # orange
    "desai_only":       "#2ca02c",   # green
    "full_minus_desai": "#2ca02c",   # green
    "full":             "#d62728",   # red
    # Munson-Dawson scenarios
    "md_only":          "#17becf",   # cyan
    "md_steady_only":   "#bcbd22",   # olive
    "full_md":          "#e377c2",   # pink
    # Interlayer scenarios
    "interlayer":       "#7f7f7f",   # gray
    "nointerlayer":     "#8c564b",   # brown
    None:               "#333333",   # dark gray
}


# ------------------------
# Plot helpers
# ------------------------
def build_color_map_for_scenarios(scenarios):
    """
    Color is determined ONLY by scenario preset.
    """
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    scenarios = list(scenarios)
    return {sc: cycle[i % len(cycle)] for i, sc in enumerate(scenarios)}

def build_linestyle_map_for_caverns(cav_labels):
    """
    Linestyle is determined by cavern label.
    """
    styles = ["-", "--", "-.", ":"]
    cav_labels = list(cav_labels)
    return {lab: styles[i % len(styles)] for i, lab in enumerate(cav_labels)}

def read_pressure_schedule(case_folder: str):
    pjson = os.path.join(case_folder, "pressure_schedule.json")
    if not os.path.isfile(pjson):
        return None, None
    import json
    with open(pjson, "r") as f:
        data = json.load(f)

    if "t_hours" in data and "p_MPa" in data:
        t = np.asarray(data["t_hours"], float)
        p = np.asarray(data["p_MPa"], float)
        return t, p

    if "t_values_s" in data and "p_values_Pa" in data:
        t = np.asarray(data["t_values_s"], float) / HOUR
        p = np.asarray(data["p_values_Pa"], float) / MPA
        return t, p

    if "t_values" in data and "p_values" in data:
        t = np.asarray(data["t_values"], float) / HOUR
        p = np.asarray(data["p_values"], float) / MPA
        return t, p

    return None, None

def plot_dilatancy_boundaries(ax, show_boundaries=None, p_min=0.01, p_max=40.0, npts=500):
    """
    Plot selected dilatancy boundaries in p–q space (MPa).

    Available boundaries:
      - ratigan_027: Ratigan 1991 (D=0.27)
      - ratigan_018: Ratigan 1991 (D=0.18)
      - spiers:      Spiers 1988 (D=0.27, b=1.9 MPa)
      - devries_comp: De Vries 2005 - compression branch
      - devries_ext:  De Vries 2005 - extension branch

    Args:
        ax: matplotlib axis
        show_boundaries: list of boundary names to show, or None for all
        p_min, p_max: pressure range (MPa)
        npts: number of points for curves
    """
    if show_boundaries is None:
        show_boundaries = ["ratigan_027", "ratigan_018", "spiers", "devries_comp", "devries_ext"]

    p = np.linspace(p_min, p_max, npts)
    I1 = 3.0 * p  # First stress invariant

    def q_from_sqrtJ2(sqrtJ2):
        return np.sqrt(3.0) * sqrtJ2

    # Dilatancy boundary colors and styles
    boundary_styles = {
        "ratigan_027":  {"color": "#7570b3", "linestyle": "--", "linewidth": 1.3, "alpha": 0.85},
        "ratigan_018":  {"color": "#7570b3", "linestyle": ":",  "linewidth": 1.3, "alpha": 0.85},
        "spiers":       {"color": "#66a61e", "linestyle": "-.", "linewidth": 1.3, "alpha": 0.85},
        "devries_comp": {"color": "#e7298a", "linestyle": "-",  "linewidth": 1.5, "alpha": 0.90},
        "devries_ext":  {"color": "#e7298a", "linestyle": "--", "linewidth": 1.5, "alpha": 0.90},
    }

    # --- Ratigan (1991): sqrt(J2) = D * I1 ---
    if "ratigan_027" in show_boundaries:
        D = 0.27
        sqrtJ2 = D * I1
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), label="Ratigan 1991 (D=0.27)",
                **boundary_styles["ratigan_027"])

    if "ratigan_018" in show_boundaries:
        D = 0.18
        sqrtJ2 = D * I1
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), label="Ratigan 1991 (D=0.18)",
                **boundary_styles["ratigan_018"])

    # --- Spiers (1988): sqrt(J2) = D*I1 + b ---
    if "spiers" in show_boundaries:
        D_sp, b_sp = 0.27, 1.9  # MPa
        sqrtJ2 = D_sp * I1 + b_sp
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), label="Spiers 1988",
                **boundary_styles["spiers"])

    # --- De Vries (2005) ---
    def devries_q(I1_MPa, psi_rad, D1=0.683, D2=0.512, T0=1.50, m=0.75, sigma_ref=1.0):
        sgn = np.sign(I1_MPa)
        sgn[sgn == 0.0] = 1.0
        denom = (np.sqrt(3.0) * np.cos(psi_rad) - D2 * np.sin(psi_rad))
        sqrtJ2_ = D1 * ((I1_MPa / (sgn * sigma_ref)) ** m) / denom + T0
        return np.sqrt(3.0) * sqrtJ2_

    psi_comp = -np.pi / 6.0  # compression: ψ = -30°
    psi_ext = np.pi / 6.0    # extension: ψ = +30°

    if "devries_comp" in show_boundaries:
        ax.plot(p, devries_q(I1, psi_comp), label="De Vries 2005 (comp)",
                **boundary_styles["devries_comp"])

    if "devries_ext" in show_boundaries:
        ax.plot(p, devries_q(I1, psi_ext), label="De Vries 2005 (ext)",
                **boundary_styles["devries_ext"])


def get_wall_indices_from_msh(msh_path):
    msh = meshio.read(msh_path)
    wall_idx = None

    if hasattr(msh, "cells_dict") and "line" in msh.cells_dict:
        wall_idx = np.unique(np.asarray(msh.cells_dict["line"]).reshape(-1))

    if wall_idx is None and isinstance(getattr(msh, "cells", None), dict):
        if "line" in msh.cells:
            wall_idx = np.unique(np.asarray(msh.cells["line"]).reshape(-1))

    if wall_idx is None:
        for cb in msh.cells:
            if getattr(cb, "type", None) == "line":
                wall_idx = np.unique(np.asarray(cb.data).reshape(-1))
                break

    if wall_idx is None:
        raise ValueError(f"No 'line' cells found in {msh_path}")

    return msh.points, wall_idx

def path_u_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.xdmf")

def path_geom_msh(case_folder):
    return os.path.join(case_folder, "operation", "mesh", "geom.msh")

def path_p_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "p_elems", "p_elems.xdmf")

def path_q_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "q_elems", "q_elems.xdmf")

def case_has_required_files(case_path: str) -> bool:
    return all(os.path.isfile(p) for p in [
        path_u_xdmf(case_path),
        path_geom_msh(case_path),
        path_p_xdmf(case_path),
        path_q_xdmf(case_path),
    ])

def load_wall_points(case_folder):
    u_xdmf = path_u_xdmf(case_folder)
    msh_path = path_geom_msh(case_folder)

    points, _, _ = post.read_node_vector(u_xdmf)
    points_msh, wall_idx_msh = get_wall_indices_from_msh(msh_path)

    mapping = post.build_mapping(points_msh, points)
    wall_idx = np.array([mapping[i] for i in wall_idx_msh], dtype=int)

    wall_points = points[wall_idx]
    order = np.argsort(wall_points[:, 2])
    return wall_points[order]

def auto_generate_probes_from_wall_points(wall_points_sorted_z, n_bend_probes=2, min_gap_idx=5):
    pts = wall_points_sorted_z
    x = pts[:, 0]
    z = pts[:, 2]

    idx_bottom = int(np.argmin(z))
    idx_top    = int(np.argmax(z))
    z_mid      = 0.5 * (z[idx_bottom] + z[idx_top])
    idx_mid    = int(np.argmin(np.abs(z - z_mid)))

    base_idx = [idx_top, idx_mid, idx_bottom]

    dx  = np.gradient(x)
    dz_ = np.gradient(z)
    ddx = np.gradient(dx)
    ddz = np.gradient(dz_)

    denom = (dx*dx + dz_*dz_)**1.5
    denom[denom == 0.0] = np.inf
    curvature = np.abs(dx * ddz - dz_ * ddx) / denom
    curvature[np.isnan(curvature)] = 0.0

    idx_all = np.argsort(curvature)[::-1]

    def too_close(new_i, existing, gap=min_gap_idx):
        return any(abs(new_i - ei) < gap for ei in existing)

    bend_idx = []
    for i in idx_all:
        if len(bend_idx) >= n_bend_probes:
            break
        if i in base_idx:
            continue
        if too_close(i, base_idx + bend_idx):
            continue
        bend_idx.append(int(i))

    bend_idx = sorted(bend_idx, key=lambda k: z[k], reverse=True)

    probes = {
        "top": pts[idx_top],
        "bend1": pts[bend_idx[0]] if len(bend_idx) > 0 else pts[idx_mid],
        "bend2": pts[bend_idx[1]] if len(bend_idx) > 1 else pts[idx_mid],
        "mid": pts[idx_mid],
        "bottom": pts[idx_bottom],
    }
    return probes

def read_stress_paths(case_folder, probes_dict):
    p_path = path_p_xdmf(case_folder)
    q_path = path_q_xdmf(case_folder)

    points_p, time_list, p_elems = post.read_cell_scalar(p_path)
    _, time_list2, q_elems = post.read_cell_scalar(q_path)

    n = min(len(time_list), len(time_list2))
    time_list = time_list[:n]
    p_elems = p_elems[:n]
    q_elems = q_elems[:n]

    out = {}
    for key, probe_xyz in probes_dict.items():
        idx = post.find_closest_point(probe_xyz, points_p)
        p = -p_elems[:, idx] / MPA
        q =  q_elems[:, idx] / MPA
        out[key] = (p, q)

    return out

def pick_one_case_per_series(cases_meta):
    out = {}
    for m in cases_meta:
        key = (m.get("cavern_label"), m.get("scenario_preset"), m.get("pressure_scenario"))
        if key not in out:
            out[key] = m
    return list(out.values())

# =============================================================================
# MAIN
# =============================================================================
def get_case_color_and_style(cavern_label, scenario_preset):
    """Get color (by scenario) and linestyle (by scenario) for a case."""
    color = SCENARIO_COLORS.get(scenario_preset, CAVERN_COLORS.get(cavern_label, "#333333"))
    linestyle = SCENARIO_LINESTYLES.get(scenario_preset, "-")
    return color, linestyle


def plot_combined(cases_meta, stress_by_series):
    """Plot all cases combined in one figure."""
    probe_types = ["top", "mid", "bottom", "bend1", "bend2"]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    for i, ptype in enumerate(probe_types):
        ax = axes[i]

        # Plot dilatancy boundaries (only show legend on first plot)
        plot_dilatancy_boundaries(ax, show_boundaries=SHOW_DILATANCY)

        for (cav, sc), d in stress_by_series.items():
            p, q = d[ptype]
            color, linestyle = get_case_color_and_style(cav, sc)
            label = f"{cav} | {sc}" if sc is not None else f"{cav}"

            ax.plot(p, q, linewidth=2.0, color=color, linestyle=linestyle, label=label)
            ax.scatter(p[-1], q[-1], s=30, edgecolors="black", linewidths=0.6,
                       color=color, zorder=5)

        ax.set_title(f"p–q stress path: {ptype}")
        ax.set_xlabel("Mean stress p (MPa)")
        ax.set_ylabel("Von Mises q (MPa)")
        ax.grid(True, alpha=0.3)

    # Pressure schedule subplot: show first matched case
    axp = axes[5]
    tH, pMPa = read_pressure_schedule(cases_meta[0]["case_path"])
    if tH is None:
        axp.text(0.5, 0.5, "No pressure_schedule.json found.",
                 ha="center", va="center", transform=axp.transAxes)
        axp.axis("off")
    else:
        axp.plot(tH / 24.0, pMPa, linewidth=2.0, color="#1f77b4")
        axp.set_title("Pressure schedule")
        axp.set_xlabel("Time (days)")
        axp.set_ylabel("Pressure (MPa)")
        axp.grid(True, alpha=0.3)

    # Deduplicate legend (show unique case labels + dilatancy boundaries)
    handles, labels = axes[0].get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h
    fig.legend(uniq.values(), uniq.keys(), loc="upper center", bbox_to_anchor=(0.5, 0.96),
               ncol=min(5, len(uniq)), frameon=True, fontsize=9)

    fig.tight_layout(rect=[0, 0, 1, 0.92])

    outname = f"pq_paths_combined_pressure={SELECT.get('pressure')}_scenario={SELECT.get('scenario')}.png"
    outpath = os.path.join(OUT_DIR, outname.replace(" ", "").replace("'", "").replace("[", "").replace("]", ""))
    fig.savefig(outpath, dpi=DPI)
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close(fig)


def plot_separate(cases_meta, stress_by_series):
    """Plot each case in a separate figure."""
    probe_types = ["top", "mid", "bottom", "bend1", "bend2"]

    for (cav, sc), d in stress_by_series.items():
        fig, axes = plt.subplots(2, 3, figsize=(14, 8))
        axes = axes.flatten()

        case_label = f"{cav} | {sc}" if sc is not None else f"{cav}"
        color, linestyle = get_case_color_and_style(cav, sc)

        for i, ptype in enumerate(probe_types):
            ax = axes[i]

            # Plot dilatancy boundaries
            plot_dilatancy_boundaries(ax, show_boundaries=SHOW_DILATANCY)

            p, q = d[ptype]
            ax.plot(p, q, linewidth=2.5, color=color, linestyle=linestyle, label=case_label)
            ax.scatter(p[-1], q[-1], s=40, edgecolors="black", linewidths=0.8,
                       color=color, zorder=5)

            ax.set_title(f"p–q stress path: {ptype}")
            ax.set_xlabel("Mean stress p (MPa)")
            ax.set_ylabel("Von Mises q (MPa)")
            ax.grid(True, alpha=0.3)
            if i == 0:
                ax.legend(loc="upper left", fontsize=9, frameon=True)

        # Pressure schedule subplot
        axp = axes[5]
        # Find the case metadata for this series
        case_path = None
        for m in cases_meta:
            if m.get("cavern_label") == cav and m.get("scenario_preset") == sc:
                case_path = m["case_path"]
                break

        if case_path:
            tH, pMPa = read_pressure_schedule(case_path)
            if tH is not None:
                axp.plot(tH / 24.0, pMPa, linewidth=2.0, color=color)
                axp.set_title("Pressure schedule")
                axp.set_xlabel("Time (days)")
                axp.set_ylabel("Pressure (MPa)")
                axp.grid(True, alpha=0.3)
            else:
                axp.text(0.5, 0.5, "No pressure data", ha="center", va="center", transform=axp.transAxes)
                axp.axis("off")
        else:
            axp.axis("off")

        fig.suptitle(f"Stress State: {case_label}", fontsize=12, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])

        safe_name = f"{cav}_{sc}".replace(" ", "_").replace("/", "_")
        outname = f"pq_paths_{safe_name}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Print configuration summary
    print("=" * 70)
    print("PLOT STRESS STATE - Configuration")
    print("=" * 70)
    print(f"  ROOT:       {ROOT}")
    print(f"  Pressure:   {SELECT.get('pressure')}")
    print(f"  Scenario:   {SELECT.get('scenario')}")
    print(f"  Caverns:    {SELECT.get('caverns')}")
    print(f"  Plot mode:  {PLOT_MODE}")
    print(f"  Dilatancy:  {SHOW_DILATANCY}")
    print("=" * 70)

    all_cases = detect_layout_and_collect_cases(ROOT)
    # Keep only cases that actually have the fields we need for this plot
    all_cases = [m for m in all_cases if case_has_required_files(m["case_path"])]

    cases_meta = filter_cases(all_cases, SELECT)
    if not cases_meta:
        print("[DEBUG] Examples of what was found (first 15):")
        for m in all_cases[:15]:
            print(" -", m.get("cavern_label"), m.get("scenario_preset"), m.get("pressure_scenario"), m.get("case_name"))
        raise RuntimeError(f"No cases matched SELECT={SELECT}")

    if ONE_CASE_PER_SERIES:
        cases_meta = pick_one_case_per_series(cases_meta)

    print(f"[INFO] Found {len(cases_meta)} case(s) to plot:")
    for m in cases_meta:
        print(f"  - {m.get('cavern_label')} | {m.get('scenario_preset')} | {m.get('case_name')}")

    # Build stress paths
    stress_by_series = {}  # key -> dict(probe -> (p,q))
    for m in cases_meta:
        folder = m["case_path"]
        try:
            wall_points = load_wall_points(folder)
            probes = auto_generate_probes_from_wall_points(wall_points, n_bend_probes=2, min_gap_idx=5)
            series_key = (m.get("cavern_label"), m.get("scenario_preset"))
            stress_by_series[series_key] = read_stress_paths(folder, probes)
        except Exception as e:
            print(f"[WARN] Skipping {m.get('case_name')}: {e}")

    if not stress_by_series:
        raise RuntimeError("No valid stress data could be loaded.")

    # Plot based on mode
    if PLOT_MODE == "combined":
        plot_combined(cases_meta, stress_by_series)
    elif PLOT_MODE == "separate":
        plot_separate(cases_meta, stress_by_series)
    else:
        print(f"[WARN] Unknown PLOT_MODE '{PLOT_MODE}', defaulting to combined")
        plot_combined(cases_meta, stress_by_series)

    print("[DONE]")

if __name__ == "__main__":
    main()

