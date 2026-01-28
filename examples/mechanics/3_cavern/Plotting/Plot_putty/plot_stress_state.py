import os
import numpy as np
import matplotlib.pyplot as plt
import meshio

import safeincave.PostProcessingTools as post

from case_index import (
    detect_layout_and_collect_cases,
    filter_cases,
)

# =============================================================================
# USER SELECTION
# =============================================================================
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/3_cavern/output"

SELECT = {
    # can be None OR list OR string
    "caverns": None,                 # e.g. ["Regular"] or ["regular600"] or None
    "pressure": "sinus",             # "sinus"/"irregular"/"linear"/"csv_profile"/None
    "scenario": ["disloc_old_only", "disloc_new_only"],  # preset(s) from ScenarioTest, or None
    "n_cycles": None,
    "operation_days": None,
    "case_contains": None,
}

# If True: pick only one case per (cavern, scenario) to avoid duplicates
ONE_CASE_PER_SERIES = True

# Headless output
OUT_DIR = os.path.join(ROOT, "_figures")
SHOW = False
DPI = 180

MPA = 1e6
HOUR = 3600.0
DAY = 24.0 * HOUR

# Stable ordering for legend/grouping
CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt", "IrregularFine"]
SCENARIO_ORDER = ["disloc_old_only", "disloc_new_only", "desai_only", "full_minus_desai", "full", None]

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

def plot_dilatancy_boundary(ax,
                            D1=0.683, D2=0.512, m=0.75, T0=1.5,
                            sigma_ref=1.0,  # MPa
                            p_min=0.01, p_max=40.0, npts=400):
    p = np.linspace(p_min, p_max, npts)    # MPa
    I1 = 3.0 * p                           # MPa

    def q_from_I1(I1_MPa, psi_rad):
        sgn = np.sign(I1_MPa)
        sgn[sgn == 0.0] = 1.0
        denom = (np.sqrt(3.0) * np.cos(psi_rad) - D2 * np.sin(psi_rad))
        sqrtJ2 = (D1 * ((I1_MPa / (sgn * sigma_ref)) ** m) + T0) / denom
        return np.sqrt(3.0) * sqrtJ2

    psi_c = -np.pi/6.0
    psi_e =  np.pi/6.0

    q_c = q_from_I1(I1, psi_c)
    q_e = q_from_I1(I1, psi_e)

    ax.plot(p, q_c, "-", linewidth=1.2, alpha=0.9, label="RD – compression")
    ax.plot(p, q_e, "-", linewidth=1.2, alpha=0.9, label="RD – extension")

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
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

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

    # Order caverns for linestyle mapping
    cav_labels = []
    scenarios = []
    for m in cases_meta:
        if m.get("cavern_label") not in cav_labels:
            cav_labels.append(m.get("cavern_label"))
        if m.get("scenario_preset") not in scenarios:
            scenarios.append(m.get("scenario_preset"))

    cav_labels_sorted = [c for c in CAVERN_ORDER if c in cav_labels] + [c for c in cav_labels if c not in CAVERN_ORDER]
    scenarios_sorted = [s for s in SCENARIO_ORDER if s in scenarios] + [s for s in scenarios if s not in SCENARIO_ORDER]

    scenario_colors = build_color_map_for_scenarios(scenarios_sorted)
    cavern_styles = build_linestyle_map_for_caverns(cav_labels_sorted)

    # Build stress paths
    stress_by_series = {}  # key -> dict(probe -> (p,q))
    for m in cases_meta:
        folder = m["case_path"]
        wall_points = load_wall_points(folder)
        probes = auto_generate_probes_from_wall_points(wall_points, n_bend_probes=2, min_gap_idx=5)
        series_key = (m.get("cavern_label"), m.get("scenario_preset"))
        stress_by_series[series_key] = read_stress_paths(folder, probes)

    probe_types = ["top", "mid", "bottom", "bend1", "bend2"]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    for i, ptype in enumerate(probe_types):
        ax = axes[i]
        plot_dilatancy_boundary(ax)

        for (cav, sc), d in stress_by_series.items():
            p, q = d[ptype]
            col = scenario_colors.get(sc, "C0")
            ls = cavern_styles.get(cav, "-")
            label = f"{cav} | {sc}" if sc is not None else f"{cav} | (no scenario)"

            ax.plot(p, q, linewidth=2.0, color=col, linestyle=ls, label=label)
            ax.scatter(p[-1], q[-1], s=30, edgecolors="black", linewidths=0.6,
                       color=col, zorder=5)

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
        axp.plot(tH / 24.0, pMPa, linewidth=2.0)
        axp.set_title("Pressure schedule (from first selected case)")
        axp.set_xlabel("Time (days)")
        axp.set_ylabel("Pressure (MPa)")
        axp.grid(True, alpha=0.3)

    # Deduplicate legend
    handles, labels = axes[0].get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h
    fig.legend(uniq.values(), uniq.keys(), loc="upper center", bbox_to_anchor=(0.5, 0.94),
               ncol=min(6, len(uniq)), frameon=True)

    fig.tight_layout(rect=[0, 0, 1, 0.92])

    outname = f"pq_paths_pressure={SELECT.get('pressure')}_scenario={SELECT.get('scenario')}.png"
    outpath = os.path.join(OUT_DIR, outname.replace(" ", ""))
    fig.savefig(outpath, dpi=DPI)
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close(fig)

if __name__ == "__main__":
    main()

