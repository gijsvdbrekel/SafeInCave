import os
import json
import numpy as np
import matplotlib.pyplot as plt
import meshio

import safeincave.PostProcessingTools as post


# =============================================================================
# USER SELECTION (edit only this block later)
# =============================================================================
ROOT = r"/home/gvandenbrekel/SafeInCave/OutputNobian"

SELECT = {
    # Cavern groups are folder names inside OutputNobian (Regular, Irregular, Tilt, ...)
    # None = all
    "caverns": ["Regular"],      # e.g. ["Regular"] or ["Irregular"] or None = all

    # Pressure scheme from pressure_schedule.json: "sinus", "irregular", "csv_profile", "linear", ...
    # None = all
    "pressure": None,

    # Scenario preset (ScenarioTest): "desai_only", "full", "full_minus_desai", "disloc_old_only", ...
    # None = all
    "scenario": None,

    # Optional numeric filters
    "n_cycles": None,
    "operation_days": None,

    # Optional substring filter on the case folder name
    "case_name_contains": None,

    # If True -> make one figure per cavern group
    "separate_per_cavern": True,
}



MPA = 1e6
HOUR = 3600.0
DAY = 24.0 * HOUR

# ------------------------
# Naming + consistent colors
# ------------------------
CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt", "IrregularFine"]

def cavern_label_from_group(group_folder: str) -> str:
    """
    group folder looks like: Asymmetric_sinus_600 or Tilt_irregular_600
    """
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



def read_pressure_schedule(case_folder: str):
    """
    Returns (t_hours, p_MPa).
    Supports both new JSON keys (t_hours/p_MPa) and old keys (t_values/p_values in s/Pa).
    """
    import os, json
    import numpy as np

    pjson = os.path.join(case_folder, "pressure_schedule.json")
    if not os.path.isfile(pjson):
        return None, None

    with open(pjson, "r") as f:
        data = json.load(f)

    # New format (preferred)
    if "t_hours" in data and "p_MPa" in data:
        t = np.asarray(data["t_hours"], dtype=float)
        p = np.asarray(data["p_MPa"], dtype=float)
        return t, p

    # Old format (seconds / Pa)
    if "t_values" in data and "p_values" in data:
        # Import your unit constants in the calling file, or define here:
        HOUR = 3600.0
        MPA = 1e6
        t = np.asarray(data["t_values"], dtype=float) / HOUR
        p = np.asarray(data["p_values"], dtype=float) / MPA
        return t, p
    
    if "t_values_s" in data and "p_values_Pa" in data:
        t = np.asarray(data["t_values_s"], dtype=float) / HOUR
        p = np.asarray(data["p_values_Pa"], dtype=float) / MPA
        return t, p


    raise KeyError(f"Pressure JSON has unexpected keys: {list(data.keys())}")


def build_color_map(labels):
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    return {lab: cycle[i % len(cycle)] for i, lab in enumerate(labels)}

# ------------------------
# Pressure schedule reader (Pressure_sinus / Pressure_irregular folders)
# ------------------------
def read_pressure_schedule_from_pressure_folder(root_folder: str, prefix: str):
    """
    Finds folder starting with prefix (case-insensitive), reads pressure_schedule*.json inside it.
    Example:
      prefix="Pressure_sinus" -> OutputNobian/Pressure_sinus/pressure_schedule.json
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

        t_hours, p_MPa = read_pressure_schedule(folder_path)  # or the folder containing the json
        if t_hours is None:
            return None, None
        t_days = t_hours / 24.0

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
# Paths for your folder-per-field layout inside operation/
# ------------------------
def path_u_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.xdmf")

def path_geom_msh(case_folder):
    return os.path.join(case_folder, "operation", "mesh", "geom.msh")

def path_p_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "p_elems", "p_elems.xdmf")

def path_q_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "q_elems", "q_elems.xdmf")

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
        p = -p_elems[:, idx] / MPA   # Here was the old /3 corrections
        q =  q_elems[:, idx] / MPA
        out[key] = (p, q)

    return out


import re

def _safe_read_json(path: str):
    try:
        with open(path, "r") as f:
            return json.load(f)
    except Exception:
        return None


def read_case_metadata(case_folder: str) -> dict:
    """
    Reads pressure_schedule.json and returns metadata.

    Key point:
      - cavern_group comes from folder name: OutputNobian/<cavern_group>/<case_...>
      - pressure_scenario comes ONLY from pressure_schedule.json (if possible)

    This avoids "Irregular" cavern vs "irregular" pressure collisions.
    """
    meta = {
        "case_folder": case_folder,
        "case_name": os.path.basename(case_folder),
        "cavern_group": os.path.basename(os.path.dirname(case_folder)),
        "pressure_scenario": None,
        "scenario": None,
        "n_cycles": None,
        "operation_days": None,
        "mode": None,
    }

    pjson = os.path.join(case_folder, "pressure_schedule.json")
    data = _safe_read_json(pjson)

    if isinstance(data, dict):
        # --- Scenario preset (ScenarioTest) ---
        # In ScenarioTest.py this is stored as "scenario": "desai_only"/"full"/...
        # In some older scripts "scenario" was used for pressure scheme -> handle that below.
        meta["scenario"] = data.get("scenario", None)

        # --- Pressure scheme ---
        # Preferred key in your newer ScenarioTest: "pressure_scenario"
        if "pressure_scenario" in data:
            meta["pressure_scenario"] = data.get("pressure_scenario")
        else:
            # Some scripts store pressure scheme under "scenario"
            sc = data.get("scenario", None)
            if sc in ("sinus", "irregular", "csv_profile", "linear"):
                meta["pressure_scenario"] = sc
                meta["scenario"] = None  # don't treat pressure scheme as scenario preset

        meta["mode"] = data.get("mode", None)
        meta["n_cycles"] = data.get("n_cycles", None)
        meta["operation_days"] = data.get("operation_days", None)
        return meta

    # ---- Fallback only if json missing (less reliable) ----
    # Still: DO NOT use cavern_group to infer pressure.
    name = meta["case_name"].lower()

    # pressure scheme guessed from CASE folder name only
    if "sinus" in name:
        meta["pressure_scenario"] = "sinus"
    elif "irregular" in name:
        meta["pressure_scenario"] = "irregular"
    elif "csv" in name:
        meta["pressure_scenario"] = "csv_profile"
    elif "linear" in name:
        meta["pressure_scenario"] = "linear"

    # scenario preset guessed from case name
    for s in ["desai_only", "disloc_old_only", "disloc_new_only", "full_minus_desai", "full"]:
        if f"_{s}_" in name or name.startswith(f"case_{s}_"):
            meta["scenario"] = s

    # numbers guessed from case name
    m = re.search(r"(\d+)\s*cyc", name)
    if m:
        meta["n_cycles"] = int(m.group(1))
    m = re.search(r"(\d+)\s*days", name)
    if m:
        meta["operation_days"] = int(m.group(1))
        
        
    # --- normalize strings (avoid filter mismatches) ---
    if meta["pressure_scenario"] is not None:
        meta["pressure_scenario"] = str(meta["pressure_scenario"]).lower()

    if meta["scenario"] is not None:
        meta["scenario"] = str(meta["scenario"]).lower()

    if meta["mode"] is not None:
        meta["mode"] = str(meta["mode"]).lower()


    return meta


def case_has_required_files(case_path: str) -> bool:
    required = [
        path_u_xdmf(case_path),
        path_geom_msh(case_path),
        path_p_xdmf(case_path),
        path_q_xdmf(case_path),
    ]
    return all(os.path.isfile(p) for p in required)


def index_all_cases(root: str) -> list[dict]:
    """
    Walks: ROOT/<cavern_group>/<case_*> and returns list of metadata dicts.
    Skips Pressure_* folders.
    """
    out = []
    for group in sorted(os.listdir(root)):
        gpath = os.path.join(root, group)
        if not os.path.isdir(gpath):
            continue
        if group.lower().startswith("pressure_"):
            continue

        for sub in sorted(os.listdir(gpath)):
            if not sub.lower().startswith("case_"):
                continue
            cpath = os.path.join(gpath, sub)
            if not os.path.isdir(cpath):
                continue
            if not case_has_required_files(cpath):
                continue

            out.append(read_case_metadata(cpath))

    return out


def filter_cases(cases: list[dict], sel: dict) -> list[dict]:
    def ok(m: dict) -> bool:
        cavs = sel.get("caverns", None)
        if cavs is not None:
            # accept either exact group folder name OR the pretty label
            group = m.get("cavern_group", "")
            label = cavern_label_from_group(group)
            if (group not in cavs) and (label not in cavs):
                return False


        p = sel.get("pressure", None)
        if p is not None:
            if (m.get("pressure_scenario") or "").lower() != p.lower():
                return False

        sc = sel.get("scenario", None)
        if sc is not None:
            if (m.get("scenario") or "").lower() != sc.lower():
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




def main():
    # 1) Index + filter cases
    all_cases = index_all_cases(ROOT)
    cases_meta = filter_cases(all_cases, SELECT)

    if not cases_meta:
        # Laat meteen zien wat er wél gevonden is (handig debuggen)
        print("[INFO] Found cases (unfiltered):")
        for m in all_cases[:20]:
            print(" -", m["cavern_group"], m["case_name"], "| pressure:", m.get("pressure_scenario"), "| scenario:", m.get("scenario"))
        raise RuntimeError("No cases matched your SELECT filters. Check SELECT block at the top.")

    print(f"[OK] Selected {len(cases_meta)} case(s):")
    for m in cases_meta:
        print(" -", m["cavern_group"], "/", m["case_name"],
              "| pressure:", m.get("pressure_scenario"),
              "| scenario:", m.get("scenario"),
              "| n_cycles:", m.get("n_cycles"),
              "| days:", m.get("operation_days"))

    # 2) Decide labels (= cavern groups)
    labels_present = []
    for m in cases_meta:
        if m["cavern_group"] not in labels_present:
            labels_present.append(m["cavern_group"])

    # Keep your preferred ordering if you want
    labels_sorted = [l for l in CAVERN_ORDER if l in labels_present] + \
                    [l for l in labels_present if l not in CAVERN_ORDER]

    color_map = build_color_map(labels_sorted)

    # 3) Build stress paths per label (if multiple cases per label match, last one wins)
    # If you expect multiple cases per same cavern, you should tighten SELECT or store list per label.
    stress_by_label = {}
    probes_by_label = {}

    for m in cases_meta:
        lab = m["cavern_group"]
        folder = m["case_folder"]

        wall_points = load_wall_points(folder)
        probes = auto_generate_probes_from_wall_points(wall_points, n_bend_probes=2, min_gap_idx=5)

        stress_by_label[lab] = read_stress_paths(folder, probes)
        probes_by_label[lab] = probes

    # 4) Optional: shape plot for one cavern
    TARGET_LABEL_FOR_SHAPE = labels_sorted[0] if labels_sorted else None
    if TARGET_LABEL_FOR_SHAPE is not None:
        # find that case folder
        target_folder = None
        for m in cases_meta:
            if m["cavern_group"] == TARGET_LABEL_FOR_SHAPE:
                target_folder = m["case_folder"]
                break

        if target_folder is not None:
            wall_points, wall_u, _ = load_wall_points_and_u(target_folder)
            probes = probes_by_label[TARGET_LABEL_FOR_SHAPE]
            fig2, ax2 = plt.subplots(figsize=(6, 6))
            plot_cavern_shape_with_probes(ax2, wall_points, wall_u, probes, scale=1.0)
            ax2.set_title(f"{TARGET_LABEL_FOR_SHAPE} cavern shape + probes")

    # 5) p–q plots per probe type
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
            ax.scatter(p[-1], q[-1], s=30, edgecolors="black", linewidths=0.6,
                       color=color_map[lab], zorder=5)

        ax.set_title(f"p–q stress path: {ptype}")
        ax.set_xlabel("Mean stress p (MPa)")
        ax.set_ylabel("Von Mises q (MPa)")
        ax.grid(True, alpha=0.3)

    # 6) Pressure schedule plot: take it from the first selected case (most robust)
    axp = axes[5]
    ref_case = cases_meta[0]["case_folder"]
    tH, pMPa = read_pressure_schedule(ref_case)
    if tH is None:
        axp.text(0.5, 0.5, "No pressure_schedule.json found in selected case.",
                 ha="center", va="center", transform=axp.transAxes)
        axp.axis("off")
    else:
        axp.plot(tH / 24.0, pMPa, linewidth=2.0)
        axp.set_title("Pressure schedule (from selected case)")
        axp.set_xlabel("Time (days)")
        axp.set_ylabel("Pressure (MPa)")
        axp.grid(True, alpha=0.3)

    # 7) Legend once
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

    ptxt = SELECT["pressure"] if SELECT["pressure"] is not None else "mixed"
    stxt = SELECT["scenario"] if SELECT["scenario"] is not None else "mixed"
    fig.suptitle(f"p–q stress paths | pressure={ptxt} | scenario={stxt}", y=0.98, fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.92])

    plt.show()

if __name__ == "__main__":
    main()




