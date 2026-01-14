import os
import json
import numpy as np
import matplotlib.pyplot as plt
import meshio

import safeincave.PostProcessingTools as post

MPA = 1e6
HOUR = 3600.0
DAY = 24.0 * HOUR

# ------------------------
# Naming + consistent colors
# ------------------------
CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt", "IrregularFine"]

def cavern_label_from_folder(folder_name: str) -> str:
    # e.g. "Asymmetric_irregular_600" -> "Asymmetric"
    return folder_name.split("_")[0]

def pressure_scheme_from_folder(folder_name: str) -> str:
    # e.g. "Asymmetric_irregular_600" -> "irregular"
    parts = folder_name.split("_")
    return parts[1].lower() if len(parts) > 1 else ""

def build_color_map(labels):
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    return {lab: cycle[i % len(cycle)] for i, lab in enumerate(labels)}

# ------------------------
# Pressure schedule reader (FIXED: actually uses prefix)
# ------------------------
def read_pressure_schedule_from_pressure_folder(root_folder: str, prefix: str):
    """
    Looks for a folder starting with `prefix` (case-insensitive) and reads the first file
    starting with 'pressure_schedule' inside it.
    """
    prefix_l = prefix.lower()

    for folder_name in sorted(os.listdir(root_folder)):
        if not folder_name.lower().startswith(prefix_l):
            continue

        folder_path = os.path.join(root_folder, folder_name)
        if not os.path.isdir(folder_path):
            continue

        candidate = None
        for fn in sorted(os.listdir(folder_path)):
            if fn.lower().startswith("pressure_schedule") and not fn.lower().endswith("zone.identifier"):
                candidate = os.path.join(folder_path, fn)
                break

        if candidate is None:
            print(f"[WARN] Found {folder_name} but no pressure_schedule file inside.")
            continue

        with open(candidate, "r") as f:
            data = json.load(f)

        t_days = np.asarray(data["t_values"], dtype=float) / DAY
        p_MPa = np.asarray(data["p_values"], dtype=float) / MPA

        print(f"Using pressure schedule from: {folder_name}/{os.path.basename(candidate)}")
        return t_days, p_MPa

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
        sqrtJ2 = D1 * ((I1_MPa / (sgn * sigma_ref)) ** m) / denom + T0
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

    # meshio <=4
    if hasattr(msh, "cells_dict") and "line" in msh.cells_dict:
        wall_idx = np.unique(np.asarray(msh.cells_dict["line"]).reshape(-1))

    # some versions: dict
    if wall_idx is None and isinstance(getattr(msh, "cells", None), dict):
        if "line" in msh.cells:
            wall_idx = np.unique(np.asarray(msh.cells["line"]).reshape(-1))

    # meshio >=5: list of CellBlock
    if wall_idx is None:
        for cb in msh.cells:
            if getattr(cb, "type", None) == "line":
                wall_idx = np.unique(np.asarray(cb.data).reshape(-1))
                break

    if wall_idx is None:
        raise ValueError(f"No 'line' cells found in {msh_path}")

    return msh.points, wall_idx

# ------------------------
# Probe generation (top/mid/bottom + bends)
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
# Read wall points from u + geom.msh
# ------------------------
def load_wall_points(case_folder):
    u_xdmf = os.path.join(case_folder, "u.xdmf")
    msh_path = os.path.join(case_folder, "geom.msh")

    points, _, _ = post.read_node_vector(u_xdmf)
    points_msh, wall_idx_msh = get_wall_indices_from_msh(msh_path)

    mapping = post.build_mapping(points_msh, points)
    wall_idx = np.array([mapping[i] for i in wall_idx_msh], dtype=int)

    wall_points = points[wall_idx]
    order = np.argsort(wall_points[:, 2])
    return wall_points[order]

def load_wall_points_and_u(case_folder):
    u_xdmf = os.path.join(case_folder, "u.xdmf")
    msh_path = os.path.join(case_folder, "geom.msh")

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

# ------------------------
# Stress path reader for probes
# ------------------------
def read_stress_paths(case_folder, probes_dict):
    p_path = os.path.join(case_folder, "p_elems.xdmf")
    q_path = os.path.join(case_folder, "q_elems.xdmf")

    points_p, time_list, p_elems = post.read_cell_scalar(p_path)
    _, time_list2, q_elems = post.read_cell_scalar(q_path)

    n = min(len(time_list), len(time_list2))
    time_list = time_list[:n]
    p_elems = p_elems[:n]
    q_elems = q_elems[:n]

    out = {}
    for key, probe_xyz in probes_dict.items():
        idx = post.find_closest_point(probe_xyz, points_p)
        p = -p_elems[:, idx] / (3.0 * MPA)   # your “/3” correction
        q =  q_elems[:, idx] / MPA
        out[key] = (p, q)

    return out

# ------------------------
# Main
# ------------------------
def main():
    ROOT = r"/home/gvandenbrekel/SafeInCave/OutputNobian"
    required = ["u.xdmf", "geom.msh", "p_elems.xdmf", "q_elems.xdmf"]

    TARGET_PRESSURE = "sinus"

    # collect cases
    cases = []
    for nm in sorted(os.listdir(ROOT)):
        fpath = os.path.join(ROOT, nm)
        if not os.path.isdir(fpath):
            continue
        if nm.lower().startswith("pressure_"):
            continue
        if pressure_scheme_from_folder(nm) != TARGET_PRESSURE:
            continue
        if all(os.path.isfile(os.path.join(fpath, r)) for r in required):
            cases.append((nm, fpath))

    if not cases:
        raise RuntimeError(f"No '{TARGET_PRESSURE}' case folders found with " + ", ".join(required))

    # consistent colors per cavern label
    labels_present = []
    for nm, _ in cases:
        lab = cavern_label_from_folder(nm)
        if lab not in labels_present:
            labels_present.append(lab)

    labels_sorted = [l for l in CAVERN_ORDER if l in labels_present] + \
                    [l for l in labels_present if l not in CAVERN_ORDER]
    color_map = build_color_map(labels_sorted)

    probe_types = ["top", "mid", "bottom", "bend1", "bend2"]

    # build stress paths
    stress_by_case = {}
    for nm, folder in cases:
        label = cavern_label_from_folder(nm)
        wall_points = load_wall_points(folder)
        probes = auto_generate_probes_from_wall_points(wall_points, n_bend_probes=2, min_gap_idx=5)
        stress_by_case[label] = read_stress_paths(folder, probes)

    # extra figure: shape + probes for one cavern
    TARGET = "IrregularFine"
    target_folder = None
    for nm, folder in cases:
        if cavern_label_from_folder(nm) == TARGET:
            target_folder = folder
            break

    if target_folder is not None:
        wall_points, wall_u, _ = load_wall_points_and_u(target_folder)
        probes = auto_generate_probes_from_wall_points(wall_points, n_bend_probes=2, min_gap_idx=5)
        fig2, ax2 = plt.subplots(figsize=(6, 6))
        plot_cavern_shape_with_probes(ax2, wall_points, wall_u, probes, scale=1.0)
        ax2.set_title(f"{TARGET} cavern shape + probes")
    else:
        print(f"[WARN] Could not find target cavern '{TARGET}' in cases. No shape plot made.")

    # p–q plots per probe type
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.flatten()

    for i, ptype in enumerate(probe_types):
        ax = axes[i]
        plot_dilatancy_boundary(ax)

        for lab in labels_sorted:
            if lab not in stress_by_case:
                continue
            p, q = stress_by_case[lab][ptype]
            ax.plot(p, q, linewidth=2.0, color=color_map[lab])
            ax.scatter(p[-1], q[-1], s=30, edgecolors="black", linewidths=0.6,
                       color=color_map[lab], zorder=5)

        ax.set_title(f"p–q stress path: {ptype}")
        ax.set_xlabel("Mean stress p (MPa)")
        ax.set_ylabel("Von Mises q (MPa)")
        ax.grid(True, alpha=0.3)

    # pressure schedule in bottom-right subplot (correct schedule!)
    axp = axes[5]
    tP, pP = read_pressure_schedule_from_pressure_folder(ROOT, f"Pressure_{TARGET_PRESSURE}")

    if tP is None:
        axp.text(0.5, 0.5, f"No pressure schedule found for Pressure_{TARGET_PRESSURE}*",
                 ha="center", va="center", transform=axp.transAxes)
        axp.axis("off")
    else:
        axp.plot(tP, pP, linewidth=2.0, color="darkred")
        axp.set_title(f"Pressure schedule ({TARGET_PRESSURE})")
        axp.set_xlabel("Time (days)")
        axp.set_ylabel("Pressure (MPa)")
        axp.grid(True, alpha=0.3)

    # legend once
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

    fig.suptitle(f"p–q stress paths per probe type ({TARGET_PRESSURE} cases)", y=0.98, fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.92])

    plt.show()

if __name__ == "__main__":
    main()

