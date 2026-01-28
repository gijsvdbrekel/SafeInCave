#!/usr/bin/env python3
import os
import re
import json
import numpy as np

# --- Headless plotting (NO GUI required) ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import meshio
import safeincave.PostProcessingTools as post


# =============================================================================
# USER SELECTION (edit only this block)
# =============================================================================
ROOT = "/data/home/gbrekel/SafeInCave_new/examples/mechanics/3_cavern/output"

SELECT = {
    # Cavern type as used in your Run.py CAVERN_TYPE, e.g. "regular600", "tilted1200", ...
    # None = all
    "cavern_type": None,            # e.g. ["regular600", "irregular600"] or None

    # Pressure scheme from pressure_schedule.json: "sinus", "irregular", "linear", ...
    # None = all
    "pressure": "sinus",

    # Optional numeric filters
    "n_cycles": None,
    "operation_days": None,

    # Optional substring filter on the case folder name
    "case_name_contains": None,

    # If True -> make one figure per case (instead of one combined figure)
    "separate_per_case": False,
}

# Output behavior
SAVE_PNG = True
SAVE_PDF = True
DPI = 200

MPA = 1e6
HOUR = 3600.0
DAY = 24.0 * HOUR


# =============================================================================
# Helpers: reading schedule + paths in RUN output structure
# =============================================================================
def _safe_read_json(path: str):
    try:
        with open(path, "r") as f:
            return json.load(f)
    except Exception:
        return None


def read_pressure_schedule(case_folder: str):
    """
    Returns (t_hours, p_MPa) from case_folder/pressure_schedule.json.
    Supports the new keys (t_hours/p_MPa) and your newer (t_values_s/p_values_Pa).
    """
    pjson = os.path.join(case_folder, "pressure_schedule.json")
    if not os.path.isfile(pjson):
        return None, None, None

    data = _safe_read_json(pjson)
    if not isinstance(data, dict):
        return None, None, None

    # New format (preferred)
    if "t_hours" in data and "p_MPa" in data:
        t = np.asarray(data["t_hours"], dtype=float)
        p = np.asarray(data["p_MPa"], dtype=float)
        scenario = data.get("scenario", None)  # you store scenario="sinus"/...
        return t, p, scenario

    # Newer format in your Run.py
    if "t_values_s" in data and "p_values_Pa" in data:
        t = np.asarray(data["t_values_s"], dtype=float) / HOUR
        p = np.asarray(data["p_values_Pa"], dtype=float) / MPA
        scenario = data.get("scenario", None)
        return t, p, scenario

    # Older possible format
    if "t_values" in data and "p_values" in data:
        t = np.asarray(data["t_values"], dtype=float) / HOUR
        p = np.asarray(data["p_values"], dtype=float) / MPA
        scenario = data.get("scenario", None)
        return t, p, scenario

    raise KeyError(f"Unexpected pressure JSON keys: {list(data.keys())}")


# --- Run-output structure paths: output/case_.../operation/<field>/<field>.xdmf, mesh/geom.msh ---
def path_u_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.xdmf")

def path_geom_msh(case_folder):
    return os.path.join(case_folder, "operation", "mesh", "geom.msh")

def path_p_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "p_elems", "p_elems.xdmf")

def path_q_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "q_elems", "q_elems.xdmf")


def case_has_required_files(case_path: str) -> bool:
    required = [
        os.path.join(case_path, "pressure_schedule.json"),
        path_u_xdmf(case_path),
        path_geom_msh(case_path),
        path_p_xdmf(case_path),
        path_q_xdmf(case_path),
    ]
    return all(os.path.isfile(p) for p in required)


# =============================================================================
# Mesh wall extraction + probes (same logic as your original)
# =============================================================================
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
        ax.scatter(pt[0], pt[2], s=90, edgecolors="black", linewidths=0.8, label=key)

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


# =============================================================================
# Stress path extraction
# =============================================================================
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


# =============================================================================
# Indexing + metadata (RUN output structure)
# =============================================================================
def infer_cavern_type_from_case_name(case_name: str):
    """
    Your run folders look like:
      case_sinus(10)_365days_regular600
    So cavern_type is usually after the last underscore.
    """
    parts = case_name.split("_")
    if len(parts) < 2:
        return None
    return parts[-1]


def read_case_metadata(case_folder: str) -> dict:
    name = os.path.basename(case_folder)
    meta = {
        "case_folder": case_folder,
        "case_name": name,
        "cavern_type": infer_cavern_type_from_case_name(name),
        "pressure_scenario": None,
        "mode": None,
        "n_cycles": None,
        "operation_days": None,
    }

    data = _safe_read_json(os.path.join(case_folder, "pressure_schedule.json"))
    if isinstance(data, dict):
        # In your Run.py you store:
        # data["scenario"] = PRESSURE_SCENARIO
        meta["pressure_scenario"] = data.get("scenario", None)
        meta["mode"] = data.get("mode", None)
        meta["n_cycles"] = data.get("n_cycles", None)
        meta["operation_days"] = data.get("operation_days", None)

    # Fallback parse if json missing keys
    # n_cycles from "(N)"
    m = re.search(r"\((\d+)\)", name)
    if m and meta["n_cycles"] is None:
        meta["n_cycles"] = int(m.group(1))

    # operation_days from "_365days_"
    m = re.search(r"_(\d+)\s*days_", name)
    if m and meta["operation_days"] is None:
        meta["operation_days"] = int(m.group(1))

    # normalize
    if meta["pressure_scenario"] is not None:
        meta["pressure_scenario"] = str(meta["pressure_scenario"]).lower()
    if meta["mode"] is not None:
        meta["mode"] = str(meta["mode"]).lower()
    if meta["cavern_type"] is not None:
        meta["cavern_type"] = str(meta["cavern_type"]).lower()

    return meta


def index_all_cases(root: str) -> list:
    out = []
    if not os.path.isdir(root):
        raise RuntimeError(f"ROOT does not exist: {root}")

    for sub in sorted(os.listdir(root)):
        if not sub.lower().startswith("case_"):
            continue
        cpath = os.path.join(root, sub)
        if not os.path.isdir(cpath):
            continue
        if not case_has_required_files(cpath):
            continue
        out.append(read_case_metadata(cpath))
    return out


def filter_cases(cases: list, sel: dict) -> list:
    def ok(m: dict) -> bool:
        cavs = sel.get("cavern_type", None)
        if cavs is not None:
            cavs_l = [c.lower() for c in cavs]
            if (m.get("cavern_type") or "").lower() not in cavs_l:
                return False

        p = sel.get("pressure", None)
        if p is not None:
            if (m.get("pressure_scenario") or "").lower() != p.lower():
                return False

        nc = sel.get("n_cycles", None)
        if nc is not None and m.get("n_cycles", None) != nc:
            return False

        od = sel.get("operation_days", None)
        if od is not None and m.get("operation_days", None) != od:
            return False

        contains = sel.get("case_name_contains", None)
        if contains is not None and contains.lower() not in m["case_name"].lower():
            return False

        return True

    return [m for m in cases if ok(m)]


def build_color_map(labels):
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    return {lab: cycle[i % len(cycle)] for i, lab in enumerate(labels)}


# =============================================================================
# Plotting runners
# =============================================================================
def save_fig(fig, out_dir, base_name):
    os.makedirs(out_dir, exist_ok=True)
    saved = []
    if SAVE_PNG:
        p = os.path.join(out_dir, base_name + ".png")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        saved.append(p)
    if SAVE_PDF:
        p = os.path.join(out_dir, base_name + ".pdf")
        fig.savefig(p, bbox_inches="tight")
        saved.append(p)
    return saved


def plot_one_case(case_meta: dict):
    folder = case_meta["case_folder"]
    case_name = case_meta["case_name"]

    # probes + stress paths
    wall_points = load_wall_points(folder)
    probes = auto_generate_probes_from_wall_points(wall_points, n_bend_probes=2, min_gap_idx=5)
    stress = read_stress_paths(folder, probes)

    # optional shape plot
    wall_points2, wall_u, _ = load_wall_points_and_u(folder)

    probe_types = ["top", "mid", "bottom", "bend1", "bend2"]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    for i, ptype in enumerate(probe_types):
        ax = axes[i]
        plot_dilatancy_boundary(ax)
        p, q = stress[ptype]
        ax.plot(p, q, linewidth=2.0)
        ax.scatter(p[-1], q[-1], s=30, edgecolors="black", linewidths=0.6, zorder=5)
        ax.set_title(f"{case_name} | p–q: {ptype}")
        ax.set_xlabel("Mean stress p (MPa)")
        ax.set_ylabel("Von Mises q (MPa)")
        ax.grid(True, alpha=0.3)

    # pressure schedule (from this case)
    axp = axes[5]
    tH, pMPa, _ = read_pressure_schedule(folder)
    if tH is None:
        axp.text(0.5, 0.5, "No pressure_schedule.json", ha="center", va="center", transform=axp.transAxes)
        axp.axis("off")
    else:
        axp.plot(tH / 24.0, pMPa, linewidth=2.0)
        axp.set_title("Pressure schedule")
        axp.set_xlabel("Time (days)")
        axp.set_ylabel("Pressure (MPa)")
        axp.grid(True, alpha=0.3)

    fig.tight_layout()

    fig2, ax2 = plt.subplots(figsize=(6, 6))
    plot_cavern_shape_with_probes(ax2, wall_points2, wall_u, probes, scale=1.0)
    ax2.set_title(f"{case_name} | shape + probes")
    fig2.tight_layout()

    out_dir = os.path.join(folder, "plots")
    saved = []
    saved += save_fig(fig, out_dir, "pq_paths")
    saved += save_fig(fig2, out_dir, "shape_probes")

    plt.close(fig)
    plt.close(fig2)
    return saved


def plot_combined(cases_meta: list):
    # label by cavern_type (so you get multi-cavern comparison)
    labels = []
    for m in cases_meta:
        lab = m.get("cavern_type") or m["case_name"]
        if lab not in labels:
            labels.append(lab)
    color_map = build_color_map(labels)

    # store per label (if multiple cases same cavern_type match, last wins)
    stress_by_label = {}
    probes_by_label = {}
    folder_by_label = {}

    for m in cases_meta:
        lab = m.get("cavern_type") or m["case_name"]
        folder = m["case_folder"]
        wall_points = load_wall_points(folder)
        probes = auto_generate_probes_from_wall_points(wall_points, n_bend_probes=2, min_gap_idx=5)

        stress_by_label[lab] = read_stress_paths(folder, probes)
        probes_by_label[lab] = probes
        folder_by_label[lab] = folder

    probe_types = ["top", "mid", "bottom", "bend1", "bend2"]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    for i, ptype in enumerate(probe_types):
        ax = axes[i]
        plot_dilatancy_boundary(ax)

        for lab in labels:
            if lab not in stress_by_label:
                continue
            p, q = stress_by_label[lab][ptype]
            ax.plot(p, q, linewidth=2.0, color=color_map[lab])
            ax.scatter(p[-1], q[-1], s=30, edgecolors="black", linewidths=0.6,
                       color=color_map[lab], zorder=5)

        ax.set_title(f"p–q stress path: {ptype}")
        ax.set_xlabel("Mean stress p (MPa)")
        ax.set_ylabel("Von Mises q (MPa)")
        ax.grid(True, alpha=0.3)

    # pressure schedule from FIRST case
    axp = axes[5]
    ref_case = cases_meta[0]["case_folder"]
    tH, pMPa, _ = read_pressure_schedule(ref_case)
    if tH is None:
        axp.text(0.5, 0.5, "No pressure_schedule.json in selected case.", ha="center", va="center", transform=axp.transAxes)
        axp.axis("off")
    else:
        axp.plot(tH / 24.0, pMPa, linewidth=2.0)
        axp.set_title("Pressure schedule (from first selected case)")
        axp.set_xlabel("Time (days)")
        axp.set_ylabel("Pressure (MPa)")
        axp.grid(True, alpha=0.3)

    # legend once
    handles = []
    leg_labels = []
    for lab in labels:
        h, = axes[0].plot([], [], color=color_map[lab], linewidth=3)
        handles.append(h)
        leg_labels.append(lab)

    fig.legend(handles, leg_labels, loc="upper center", bbox_to_anchor=(0.5, 0.94),
               ncol=min(6, len(leg_labels)), frameon=True)

    ptxt = SELECT["pressure"] if SELECT["pressure"] is not None else "mixed"
    fig.suptitle(f"p–q stress paths | pressure={ptxt}", y=0.98, fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.92])

    # Save combined plot into ROOT/_plots
    out_dir = os.path.join(ROOT, "_plots")
    os.makedirs(out_dir, exist_ok=True)
    tag = f"pressure={ptxt}"
    if SELECT.get("n_cycles") is not None:
        tag += f"_cycles={SELECT['n_cycles']}"
    if SELECT.get("operation_days") is not None:
        tag += f"_days={SELECT['operation_days']}"

    saved = save_fig(fig, out_dir, f"pq_paths_{tag}")
    plt.close(fig)
    return saved


def main():
    all_cases = index_all_cases(ROOT)
    cases_meta = filter_cases(all_cases, SELECT)

    if not cases_meta:
        print("[INFO] Found cases (unfiltered):")
        for m in all_cases[:30]:
            print(" -", m["case_name"], "| cavern:", m.get("cavern_type"),
                  "| pressure:", m.get("pressure_scenario"),
                  "| cycles:", m.get("n_cycles"),
                  "| days:", m.get("operation_days"))
        raise RuntimeError("No cases matched your SELECT filters. Fix SELECT or check ROOT.")

    print(f"[OK] Selected {len(cases_meta)} case(s):")
    for m in cases_meta:
        print(" -", m["case_name"],
              "| cavern:", m.get("cavern_type"),
              "| pressure:", m.get("pressure_scenario"),
              "| cycles:", m.get("n_cycles"),
              "| days:", m.get("operation_days"))

    saved_all = []

    if SELECT.get("separate_per_case", False):
        for m in cases_meta:
            saved = plot_one_case(m)
            saved_all.extend(saved)
            print("[OK] Saved for", m["case_name"])
            for s in saved:
                print("   ", s)
    else:
        saved = plot_combined(cases_meta)
        saved_all.extend(saved)
        print("[OK] Saved combined figure(s):")
        for s in saved:
            print("  ", s)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
