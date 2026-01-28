#!/usr/bin/env python3
import os
import json
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meshio

import safeincave.PostProcessingTools as post
from case_index import detect_layout_and_collect_cases, filter_cases, one_case_per_cavern_label


MPA = 1e6
HOUR = 3600.0
DAY = 24.0 * HOUR

# For nicer legend order
CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt", "IrregularFine"]

# =============================================================================
# USER SELECTION
# =============================================================================
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/3_cavern/output"

SELECT = {
    "pressure": "sinus",          # "sinus"/"irregular"/"csv_profile"/"linear"/None
    "scenario": None,             # "full"/"desai_only"/... or None
    "caverns": None,              # e.g. ["Regular","regular600"] or None
    "n_cycles": None,
    "operation_days": None,
    "case_contains": None,

    # if multiple matching cases per cavern label:
    "one_case_per_cavern": True,
}

# Output settings
OUTDIR = os.path.join(ROOT, "_plots")
OUTNAME = "pq_stress_paths.png"
DPI = 220

# Probes
N_BEND_PROBES = 2
MIN_GAP_IDX = 5

# Optional: also save a shape+probe plot for the first selected cavern
SAVE_SHAPE_FIG = True
OUTNAME_SHAPE = "cavern_shape_probes.png"


# ------------------------
# Consistent colors per cavern label
# ------------------------
def build_color_map(labels):
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    return {lab: cycle[i % len(cycle)] for i, lab in enumerate(labels)}


# ------------------------
# Pressure schedule reader
# ------------------------
def read_pressure_schedule(case_folder: str):
    """
    Returns (t_hours, p_MPa).
    Supports:
      - new: t_hours + p_MPa
      - raw: t_values_s + p_values_Pa
      - old: t_values + p_values (s/Pa)
    """
    pjson = os.path.join(case_folder, "pressure_schedule.json")
    if not os.path.isfile(pjson):
        return None, None

    with open(pjson, "r") as f:
        data = json.load(f)

    if "t_hours" in data and "p_MPa" in data:
        t = np.asarray(data["t_hours"], dtype=float)
        p = np.asarray(data["p_MPa"], dtype=float)
        return t, p

    if "t_values_s" in data and "p_values_Pa" in data:
        t = np.asarray(data["t_values_s"], dtype=float) / HOUR
        p = np.asarray(data["p_values_Pa"], dtype=float) / MPA
        return t, p

    if "t_values" in data and "p_values" in data:
        t = np.asarray(data["t_values"], dtype=float) / HOUR
        p = np.asarray(data["p_values"], dtype=float) / MPA
        return t, p

    return None, None


# ------------------------
# RD dilatancy boundary
# ------------------------
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


# ------------------------
# Mesh wall extraction
# ------------------------
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


# ------------------------
# Probe generation
# ------------------------
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

    if len(bend_idx) < n_bend_probes:
        for i in range(len(pts)):
            if len(bend_idx) >= n_bend_probes:
                break
            if i in base_idx or i in bend_idx:
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


# ------------------------
# Paths (operation layout)
# ------------------------
def path_u_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.xdmf")

def path_geom_msh(case_folder):
    return os.path.join(case_folder, "operation", "mesh", "geom.msh")

def path_p_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "p_elems", "p_elems.xdmf")

def path_q_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "q_elems", "q_elems.xdmf")


def case_has_required_files(case_path: str) -> bool:
    req = [path_u_xdmf(case_path), path_geom_msh(case_path), path_p_xdmf(case_path), path_q_xdmf(case_path)]
    return all(os.path.isfile(p) for p in req)


# ------------------------
# Read wall points from u + geom
# ------------------------
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


def load_wall_points_and_u(case_folder):
    u_xdmf = path_u_xdmf(case_folder)
    msh_path = path_geom_msh(case_folder)

    points, time_list, u_field = post.read_node_vector(u_xdmf)
    points_msh, wall_idx_msh = get_wall_indices_from_msh(msh_path)

    mapping = post.build_mapping(points_msh, points)
    wall_idx = np.array([mapping[i] for i in wall_idx_msh], dtype=int)

    wall_points = points[wall_idx]
    wall_u = u_field[:, wall_idx, :]

    order = np.argsort(wall_points[:, 2])
    return wall_points[order], wall_u[:, order, :], time_list


def plot_cavern_shape_with_probes(ax, wall_points, wall_u, probes_dict, scale=1.0):
    wall_t0 = wall_points + scale * wall_u[0]
    wall_tf = wall_points + scale * wall_u[-1]

    ax.plot(wall_t0[:, 0], wall_t0[:, 2], "-", linewidth=2.0, label="Initial shape")
    ax.plot(wall_tf[:, 0], wall_tf[:, 2], "-", linewidth=2.0, label="Final shape")

    for key, pt in probes_dict.items():
        ax.scatter(pt[0], pt[2], s=70, edgecolors="black", linewidths=0.7, label=key)

    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.grid(True, alpha=0.3)
    ax.axis("equal")

    handles, labels = ax.get_legend_handles_labels()
    seen = set()
    h2, l2 = [], []
    for h, l in zip(handles, labels):
        if l in seen:
            continue
        seen.add(l)
        h2.append(h)
        l2.append(l)
    ax.legend(h2, l2, loc="best", fontsize=8, frameon=True)


# ------------------------
# Stress path reader for probes
# ------------------------
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


def main():
    # 1) Collect + filter cases (both layouts)
    all_cases = detect_layout_and_collect_cases(ROOT)
    cases_meta = filter_cases(all_cases, SELECT)

    # 2) ensure required files exist for this plot
    cases_meta = [c for c in cases_meta if case_has_required_files(c["case_path"])]

    if not cases_meta:
        print("[INFO] Found cases (unfiltered):")
        for m in all_cases[:40]:
            print(" -", m["case_name"],
                  "| pressure:", m.get("pressure_scenario"),
                  "| scenario:", m.get("scenario_preset"),
                  "| cavern:", m.get("cavern_label"))
        raise RuntimeError(f"No cases matched SELECT={SELECT} (or missing required operation files).")

    # Optionally keep one case per cavern label
    if SELECT.get("one_case_per_cavern", True):
        cases_meta = one_case_per_cavern_label(cases_meta)

    # 3) Labels + colors
    labels_present = []
    for m in cases_meta:
        lab = m["cavern_label"]
        if lab not in labels_present:
            labels_present.append(lab)

    labels_sorted = [l for l in CAVERN_ORDER if l in labels_present] + \
                    [l for l in labels_present if l not in CAVERN_ORDER]

    color_map = build_color_map(labels_sorted)

    # 4) Read stress paths per cavern label
    stress_by_label = {}
    probes_by_label = {}

    for m in cases_meta:
        lab = m["cavern_label"]
        folder = m["case_path"]

        wall_points = load_wall_points(folder)
        probes = auto_generate_probes_from_wall_points(wall_points, n_bend_probes=N_BEND_PROBES, min_gap_idx=MIN_GAP_IDX)

        stress_by_label[lab] = read_stress_paths(folder, probes)
        probes_by_label[lab] = probes

    # 5) Optional: shape+probes figure for first label
    if SAVE_SHAPE_FIG and labels_sorted:
        target_lab = labels_sorted[0]
        target_folder = None
        for m in cases_meta:
            if m["cavern_label"] == target_lab:
                target_folder = m["case_path"]
                break

        if target_folder is not None:
            wall_points, wall_u, _ = load_wall_points_and_u(target_folder)
            probes = probes_by_label[target_lab]

            figS, axS = plt.subplots(figsize=(6, 6))
            plot_cavern_shape_with_probes(axS, wall_points, wall_u, probes, scale=1.0)
            axS.set_title(f"{target_lab} cavern shape + probes")

            os.makedirs(OUTDIR, exist_ok=True)
            outS = os.path.join(OUTDIR, OUTNAME_SHAPE)
            figS.tight_layout()
            figS.savefig(outS, dpi=DPI, bbox_inches="tight")
            plt.close(figS)
            print("[OK] Saved:", outS)

    # 6) p–q plots per probe type
    probe_types = ["top", "mid", "bottom", "bend1", "bend2"]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    for i, ptype in enumerate(probe_types):
        ax = axes[i]
        plot_dilatancy_boundary(ax)

        for lab in labels_sorted:
            if lab not in stress_by_label:
                continue
            p, q = stress_by_label[lab][ptype]
            ax.plot(p, q, linewidth=2.0, color=color_map[lab])
            ax.scatter(p[-1], q[-1], s=28, edgecolors="black", linewidths=0.6,
                       color=color_map[lab], zorder=5)

        ax.set_title(f"p–q stress path: {ptype}")
        ax.set_xlabel("Mean stress p (MPa)")
        ax.set_ylabel("Von Mises q (MPa)")
        ax.grid(True, alpha=0.3)

    # 7) Pressure schedule plot from first selected case
    axp = axes[5]
    ref_case = cases_meta[0]["case_path"]
    tH, pMPa = read_pressure_schedule(ref_case)
    if tH is None:
        axp.text(0.5, 0.5, "No pressure_schedule.json found.", ha="center", va="center", transform=axp.transAxes)
        axp.axis("off")
    else:
        axp.plot(tH / 24.0, pMPa, linewidth=2.0)
        axp.set_title("Pressure schedule (selected case)")
        axp.set_xlabel("Time (days)")
        axp.set_ylabel("Pressure (MPa)")
        axp.grid(True, alpha=0.3)

    # 8) One legend for caverns
    handles = []
    labels = []
    for lab in labels_sorted:
        h, = axes[0].plot([], [], color=color_map[lab], linewidth=3)
        handles.append(h)
        labels.append(lab)

    fig.legend(
        handles, labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.94),
        ncol=min(6, len(labels)),
        frameon=True
    )

    ptxt = SELECT.get("pressure", None) or "ALL"
    stxt = SELECT.get("scenario", None) or "ANY"
    fig.suptitle(f"p–q stress paths | pressure={ptxt} | scenario={stxt}", y=0.98, fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.92])

    os.makedirs(OUTDIR, exist_ok=True)
    outpath = os.path.join(OUTDIR, OUTNAME)
    fig.savefig(outpath, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("[OK] Saved:", outpath)


if __name__ == "__main__":
    main()
