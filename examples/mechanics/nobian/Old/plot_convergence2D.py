import os
import json
import numpy as np
import matplotlib.pyplot as plt

import meshio
import safeincave.PostProcessingTools as post

HOUR = 3600.0
MPA = 1e6

# ---------------------------
# Developer-style convergence
# ---------------------------

class WallProfileData:
    """
    Uses the existing SafeInCave 'wall profile' approach:
    - reads wall nodes from 1D line elements in geom.msh (Wall_profile)
    - maps to xdmf nodes
    - computes a volume-of-revolution proxy convergence (%)
    """
    def __init__(self, operation_folder, scale=1.0):
        # Read displacements from xdmf
        points, self.time_list, u_field = post.read_node_vector(
            os.path.join(operation_folder, "u", "u.xdmf")
        )

        # Read msh node coordinates from msh file
        reader_msh = meshio.read(os.path.join(operation_folder, "mesh", "geom.msh"))
        points_msh = reader_msh.points

        # Get wall indices from msh file (robust across meshio APIs)
        wall_idx = None

        # Old meshio API
        if hasattr(reader_msh, "cells_dict"):
            if "line" in reader_msh.cells_dict:
                wall_idx = np.unique(np.asarray(reader_msh.cells_dict["line"]).flatten())
            else:
                # fallback: any line* (line2/line3)
                for k, v in reader_msh.cells_dict.items():
                    if isinstance(k, str) and k.startswith("line"):
                        wall_idx = np.unique(np.asarray(v).flatten())
                        break

        # New meshio API (dict format)
        elif isinstance(reader_msh.cells, dict):
            if "line" in reader_msh.cells:
                wall_idx = np.unique(np.asarray(reader_msh.cells["line"]).flatten())
            else:
                for k, v in reader_msh.cells.items():
                    if isinstance(k, str) and k.startswith("line"):
                        wall_idx = np.unique(np.asarray(v).flatten())
                        break

        # New meshio API (list of CellBlock format)
        else:
            for cb in reader_msh.cells:
                ctype = getattr(cb, "type", None)
                if isinstance(ctype, str) and ctype.startswith("line"):
                    wall_idx = np.unique(np.asarray(cb.data).flatten())
                    break

        if wall_idx is None:
            avail = []
            try:
                if isinstance(reader_msh.cells, dict):
                    avail = list(reader_msh.cells.keys())
                else:
                    avail = [getattr(cb, "type", None) for cb in reader_msh.cells]
            except Exception:
                pass
            raise ValueError(f"No line-type cells found in geom.msh. Available cell types: {avail}")

        # Map msh node indices to xdmf node indices
        mapping = post.build_mapping(points_msh, points)
        wall_idx = np.array([mapping[i] for i in wall_idx], dtype=int)

        # Coordinates + displacement time series for wall nodes
        self.wall_points = points[wall_idx]
        self.wall_u = u_field[:, wall_idx]

        # Sort by z coordinate
        order = np.argsort(self.wall_points[:, 2])
        self.wall_points = self.wall_points[order]
        self.wall_u = self.wall_u[:, order]

        self.scale = float(scale)
        self._compute_volumes_percent()

    def get_wall_coordinates(self, time_step: int):
        return self.wall_points + self.scale * self.wall_u[time_step, :]

    def _compute_volumes_percent(self):
        wall0 = self.get_wall_coordinates(0)
        vol0 = self._trapezoidal_volume(wall0[:, 2], wall0[:, 0])

        vols = []
        for k in range(len(self.time_list)):
            wallk = self.get_wall_coordinates(k)
            volk = self._trapezoidal_volume(wallk[:, 2], wallk[:, 0])
            vols.append(100.0 * abs(vol0 - volk) / abs(vol0))
        self.volumes = np.asarray(vols, float)

    @staticmethod
    def _trapezoidal_volume(z, r):
        """Volume of revolution around r=0 axis via trapezoidal rule."""
        v = 0.0
        for i in range(1, len(z)):
            R = 0.5 * (r[i] + r[i - 1])
            A = np.pi * R**2
            dz = z[i] - z[i - 1]
            v += A * dz
        return v


# ---------------------------
# Naming helpers (like your stress-state script)
# ---------------------------

CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt", "IrregularFine"]

def cavern_label_from_group(group_folder: str) -> str:
    low = group_folder.lower()
    if low.startswith("asymmetric"):
        return "Asymmetric"
    if low.startswith("multichamber"):
        return "Multichamber"
    if low.startswith("teardrop"):
        return "Teardrop"
    if low.startswith("tilt") or low.startswith("tilted"):
        return "Tilt"
    if low.startswith("regular"):
        return "Regular"
    if low.startswith("irregularfine") or low.startswith("irregular_fine"):
        return "IrregularFine"
    if low.startswith("irregular"):
        return "Irregular"
    return group_folder.split("_")[0]

def scheme_from_case_folder(case_folder_name: str) -> str:
    low = case_folder_name.lower()
    if low.startswith("case_sinus"):
        return "sinus"
    if low.startswith("case_irregular"):
        return "irregular"
    return ""


def build_color_map(labels):
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    return {lab: cycle[i % len(cycle)] for i, lab in enumerate(labels)}


# ---------------------------
# Folder structure helpers
# ROOT/<GROUP>/<CASE>/operation/u/u.xdmf
# ROOT/<GROUP>/<CASE>/operation/mesh/geom.msh
# ROOT/<GROUP>/<CASE>/pressure_schedule.json  (optional for inclusion)
# ---------------------------

def path_u_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.xdmf")

def path_geom_msh(case_folder):
    return os.path.join(case_folder, "operation", "mesh", "geom.msh")

def path_pressure_json(case_folder):
    return os.path.join(case_folder, "pressure_schedule.json")


def collect_cases_nested(ROOT, target_scheme: str):
    """
    Returns list of dicts:
      {label, group, case_name, case_path, operation_path, has_pressure_json}
    Filters by scheme using CASE folder name:
      case_sinus(...) or case_irregular(...)
    """
    cases = []
    for group in sorted(os.listdir(ROOT)):
        group_path = os.path.join(ROOT, group)
        if not os.path.isdir(group_path):
            continue

        # skip Pressure_* folders
        if group.lower().startswith("pressure_"):
            continue

        label = cavern_label_from_group(group)

        for sub in sorted(os.listdir(group_path)):
            if not sub.lower().startswith("case_"):
                continue

            case_scheme = scheme_from_case_folder(sub)
            if case_scheme != target_scheme:
                continue

            case_path = os.path.join(group_path, sub)
            if not os.path.isdir(case_path):
                continue

            u_xdmf = path_u_xdmf(case_path)
            msh = path_geom_msh(case_path)

            if not (os.path.isfile(u_xdmf) and os.path.isfile(msh)):
                # don't exclude just because pressure json is missing
                print(f"[SKIP] {group}/{sub} missing u.xdmf or geom.msh")
                continue

            pjson = path_pressure_json(case_path)
            cases.append({
                "label": label,
                "group": group,
                "case_name": sub,
                "case_path": case_path,
                "operation_path": os.path.join(case_path, "operation"),
                "has_pressure_json": os.path.isfile(pjson),
                "pressure_json": pjson,
            })

    return cases


def read_pressure_from_case(case_path: str):
    pjson = path_pressure_json(case_path)
    if not os.path.isfile(pjson):
        return None, None
    with open(pjson, "r") as f:
        data = json.load(f)
    t_h = np.asarray(data["t_values"], float) / HOUR
    p_mpa = np.asarray(data["p_values"], float) / MPA
    return t_h, p_mpa


# ---------------------------
# Plotting
# ---------------------------

def plot_scheme(ROOT: str, scheme: str, *, scale=1.0, color_map=None):
    """
    Plot one figure for a given pressure scheme ("sinus" or "irregular"):
      - top: convergence (%) for all caverns in that scheme
      - bottom: pressure schedule (MPa), using the first available json and overlaying others only if different

    Pass a GLOBAL color_map (label -> color) to keep colors consistent across schemes.
    """
    cases = collect_cases_nested(ROOT, scheme)
    if not cases:
        raise RuntimeError(f"No '{scheme}' cases found under nested folders in {ROOT}.")

    # Fallback: if no global map is passed, build one from the labels in this scheme only
    if color_map is None:
        labels_present = []
        for c in cases:
            if c["label"] not in labels_present:
                labels_present.append(c["label"])

        labels_sorted = [l for l in CAVERN_ORDER if l in labels_present] + \
                        [l for l in labels_present if l not in CAVERN_ORDER]
        color_map = build_color_map(labels_sorted)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(13, 7), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12},
    )
    fig.suptitle(f"All caverns – pressure scheme: {scheme}", fontsize=12)

    # --- convergence ---
    plotted_any = False
    for c in cases:
        lab = c["label"]
        col = color_map.get(lab, None)

        try:
            wall = WallProfileData(c["operation_path"], scale=scale)
            t_h = np.asarray(wall.time_list, float) / HOUR
            conv = wall.volumes
        except Exception as e:
            print(f"[SKIP] convergence {c['group']}/{c['case_name']}: {e}")
            continue

        ax1.plot(t_h, conv, linewidth=2.0, alpha=0.95, color=col, label=lab)
        plotted_any = True

    if not plotted_any:
        ax1.text(
            0.5, 0.5,
            "No convergence curves could be plotted (all cases failed).",
            ha="center", va="center", transform=ax1.transAxes
        )

    ax1.set_ylabel("Convergence (%)")
    ax1.grid(True, alpha=0.3)

    # legend unique (so each cavern label appears once)
    handles, labs = ax1.get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labs):
        uniq[l] = h
    if uniq:
        ax1.legend(uniq.values(), uniq.keys(), loc="best", fontsize=9, frameon=True)

    # --- pressure schedule ---
    # Pick first case that has pressure_schedule.json
    ref_t, ref_p = None, None
    ref_case_name = None

    for c in cases:
        t, p = read_pressure_from_case(c["case_path"])
        if t is not None and p is not None:
            ref_t, ref_p = t, p
            ref_case_name = c["case_name"]
            break

    if ref_t is None:
        ax2.text(
            0.5, 0.5,
            "No pressure_schedule.json found in any case folder.",
            ha="center", va="center", transform=ax2.transAxes
        )
        ax2.set_ylabel("Pressure (MPa)")
        ax2.set_xlabel("Time (hours)")
        ax2.grid(True, alpha=0.3)
        return fig

    ax2.plot(ref_t, ref_p, linewidth=1.7, label=f"pressure ({ref_case_name})")

    # Overlay other pressure schedules only if they differ from the reference
    n_diff = 0
    for c in cases:
        t2, p2 = read_pressure_from_case(c["case_path"])
        if t2 is None or p2 is None:
            continue

        # skip identical reference (or numerically the same)
        if len(t2) == len(ref_t) and np.max(np.abs(t2 - ref_t)) <= 1e-9 and np.max(np.abs(p2 - ref_p)) <= 1e-6:
            continue

        n_diff += 1
        ax2.plot(t2, p2, linewidth=0.9, alpha=0.6, label=f"pressure ({c['case_name']})")

    ax2.set_ylabel("Pressure (MPa)")
    ax2.set_xlabel("Time (hours)")
    ax2.grid(True, alpha=0.3)
    if n_diff > 0:
        ax2.legend(loc="best", fontsize=8, frameon=True)

    return fig



def main():
    ROOT = "/home/gvandenbrekel/SafeInCave/OutputNobian"

    # --- collect all labels across *all* schemes first ---
    all_cases = []
    for scheme in ("sinus", "irregular"):
        all_cases.extend(collect_cases_nested(ROOT, scheme))

    all_labels = []
    for c in all_cases:
        if c["label"] not in all_labels:
            all_labels.append(c["label"])

    labels_sorted = [l for l in CAVERN_ORDER if l in all_labels] + \
                    [l for l in all_labels if l not in CAVERN_ORDER]

    global_cmap = build_color_map(labels_sorted)

    # --- now plot using the SAME color map for every scheme ---
    for scheme in ("sinus", "irregular"):
        try:
            plot_scheme(ROOT, scheme, scale=1.0, color_map=global_cmap)
        except Exception as e:
            print(f"[ERROR] {scheme}: {e}")

    plt.show()



if __name__ == "__main__":
    main()


