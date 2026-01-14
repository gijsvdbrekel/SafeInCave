import os
import json
import numpy as np
import matplotlib.pyplot as plt
import meshio

import safeincave.PostProcessingTools as post

HOUR = 3600.0
DAY = 24.0 * HOUR
MPA = 1e6


# ---------- helpers: naming & consistent colors ----------

CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt"]

def cavern_label_from_folder(folder_name: str) -> str:
    # take text before first underscore, e.g. "Asymmetric_regular_600_convergence" -> "Asymmetric"
    return folder_name.split("_")[0]

def build_color_map(labels):
    # use Matplotlib default cycle, but pinned to labels so it stays consistent
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]

    cmap = {}
    for i, lab in enumerate(labels):
        cmap[lab] = cycle[i % len(cycle)]
    return cmap


# ---------- pressure reading ----------


def pressure_scheme_from_folder(folder_name: str) -> str:
    parts = folder_name.split("_")
    return parts[1].lower() if len(parts) > 1 else ""


def read_pressure_schedule_from_pressure_folder(root_folder, pressure_folder_prefix):
    """
    Finds a folder whose name starts with pressure_folder_prefix (case-insensitive),
    then reads a file that starts with 'pressure_schedule' (with or without .json).
    """
    prefix = pressure_folder_prefix.lower()

    for folder_name in sorted(os.listdir(root_folder)):
        if not folder_name.lower().startswith(prefix):
            continue

        folder_path = os.path.join(root_folder, folder_name)
        if not os.path.isdir(folder_path):
            continue

        # find any file like pressure_schedule or pressure_schedule.json
        candidate = None
        for fn in os.listdir(folder_path):
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



# ---------- mesh/wall extraction ----------

def get_wall_indices_from_msh(msh_path):
    msh = meshio.read(msh_path)
    wall_idx = None

    if hasattr(msh, "cells_dict"):
        if "line" in msh.cells_dict:
            wall_idx = np.unique(msh.cells_dict["line"].reshape(-1))

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


def trapezoidal_volume(z, r):
    vol = 0.0
    for i in range(1, len(z)):
        R = 0.5 * (r[i] + r[i - 1])
        A = np.pi * R**2
        dz = z[i] - z[i - 1]
        vol += A * dz
    return vol


def compute_convergence_percent(case_folder, scale=1.0):
    u_xdmf = os.path.join(case_folder, "u.xdmf")
    msh_path = os.path.join(case_folder, "geom.msh")

    if not os.path.isfile(u_xdmf):
        raise FileNotFoundError(f"Missing {u_xdmf}")
    if not os.path.isfile(msh_path):
        raise FileNotFoundError(f"Missing {msh_path}")

    points, time_list, u_field = post.read_node_vector(u_xdmf)
    points_msh, wall_idx_msh = get_wall_indices_from_msh(msh_path)

    mapping = post.build_mapping(points_msh, points)
    wall_idx = np.array([mapping[i] for i in wall_idx_msh], dtype=int)

    wall_points = points[wall_idx]
    wall_u = u_field[:, wall_idx, :]

    order = np.argsort(wall_points[:, 2])
    wall_points = wall_points[order]
    wall_u = wall_u[:, order, :]

    wall0 = wall_points + scale * wall_u[0]
    vol0 = trapezoidal_volume(wall0[:, 2], wall0[:, 0])

    conv = []
    for k in range(len(time_list)):
        wallk = wall_points + scale * wall_u[k]
        volk = trapezoidal_volume(wallk[:, 2], wallk[:, 0])
        conv.append(100.0 * abs(vol0 - volk) / abs(vol0))

    t_days = np.asarray(time_list, dtype=float) / DAY
    return t_days, np.asarray(conv, dtype=float)


# ---------- plotting ----------

def plot_set(fig_title, cases, color_map, pressure_t, pressure_p):
    fig = plt.figure(figsize=(14, 9))
    fig.suptitle(fig_title, fontsize=14)

    gs = fig.add_gridspec(2, 1, height_ratios=[2.2, 1.0], hspace=0.25)
    ax_conv = fig.add_subplot(gs[0, 0])
    ax_p = fig.add_subplot(gs[1, 0], sharex=ax_conv)

    # convergences
    for name, folder in cases:
        label = cavern_label_from_folder(name)
        col = color_map.get(label, None)
        try:
            t, c = compute_convergence_percent(folder, scale=1.0)
        except Exception as e:
            print(f"[SKIP] {name}: {e}")
            continue
        ax_conv.plot(t, c, linewidth=2.2, label=label, color=col)

    ax_conv.set_ylabel("Convergence (volume change, %)")
    ax_conv.set_xlabel("Time (days)")
    ax_conv.grid(True, alpha=0.3)

    # legend: unique labels only
    handles, labels = ax_conv.get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        uniq[l] = h
    ax_conv.legend(uniq.values(), uniq.keys(), loc="best", fontsize=9, frameon=True)

    # pressure
    if pressure_t is None:
        ax_p.text(0.5, 0.5, "pressure_schedule.json not found.", transform=ax_p.transAxes,
                  ha="center", va="center")
    else:
        ax_p.plot(pressure_t, pressure_p, linewidth=1.8)
        ax_p.set_ylabel("Cavern pressure (MPa)")
        ax_p.set_xlabel("Time (days)")
        ax_p.grid(True, alpha=0.3)

    return fig


def main():
    ROOT = r"/home/gvandenbrekel/SafeInCave/OutputNobian"

    # gather cases
    irregular_cases = []
    sinus_cases = []

    for name in sorted(os.listdir(ROOT)):
        fpath = os.path.join(ROOT, name)
        if not os.path.isdir(fpath):
            continue
        if not (os.path.isfile(os.path.join(fpath, "u.xdmf")) and os.path.isfile(os.path.join(fpath, "geom.msh"))):
            continue

        # classify by folder name
        scheme = pressure_scheme_from_folder(name)
        if scheme == "irregular":
            irregular_cases.append((name, fpath))
        elif scheme == "sinus":
            sinus_cases.append((name, fpath))


    if not irregular_cases and not sinus_cases:
        raise RuntimeError(f"No case folders with u.xdmf + geom.msh found in {ROOT}")

    # consistent colors by cavern label
    # use your known set if present; otherwise derive from cases
    labels_present = []
    for name, _ in irregular_cases + sinus_cases:
        lab = cavern_label_from_folder(name)
        if lab not in labels_present:
            labels_present.append(lab)

    # keep nicer order if possible
    labels_sorted = [l for l in CAVERN_ORDER if l in labels_present] + [l for l in labels_present if l not in CAVERN_ORDER]
    color_map = build_color_map(labels_sorted)

    # read both pressures
    tPi, pPi = read_pressure_schedule_from_pressure_folder(ROOT, "Pressure_irregular")
    tPr, pPr = read_pressure_schedule_from_pressure_folder(ROOT, "Pressure_sinus")

    # make two figures
    if irregular_cases:
        plot_set("Irregular pressure: convergence + pressure", irregular_cases, color_map, tPi, pPi)

    if sinus_cases:
        plot_set("Sinus pressure: convergence + pressure", sinus_cases, color_map, tPr, pPr)

    plt.show()


if __name__ == "__main__":
    main()
