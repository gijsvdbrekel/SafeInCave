#!/usr/bin/env python3
import os
import json
import numpy as np

# ---- headless plotting (NO GUI needed) ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import meshio
import safeincave.PostProcessingTools as post


# -------------------------
# Global constants
# -------------------------
GEOMETRY_TYPE = "irregular"  # only affects default bend probes count
hour = 3600.0
day = 24.0 * hour
MPa = 1e6

# Which caverns to show (order) (pretty labels)
CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt"]

# =============================================================================
# USER SELECTION (edit only this block)
# =============================================================================
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/3_cavern/output"

SELECT = {
    "pressure": "sinus",        # "sinus", "irregular", "csv_profile", "linear", or None (=all)

    # Flat output layout: cavern type is inferred from case folder name suffix, e.g. *_regular600
    # Use None (=all) or list like ["regular600","irregular600","tilted1200"]
    "caverns": ["regular600"],

    "case_contains": None,      # e.g. "365days" or "8cyc" or "desai_only" or None
}

# Saving (no interactive slider/buttons)
SAVE_DIR_NAME = "_dashboard_outputs"   # saved under ROOT
DPI = 200
SAVE_PNG = True
SAVE_PDF = False


# -------------------------
# UI / layout (kept, but NO widgets)
# -------------------------
def apply_grey_theme(fig, axes, transparent=True, grid_color="0.92", back_color="0.85"):
    fig.patch.set_facecolor("#212121ff")
    if transparent:
        fig.patch.set_alpha(0.0)
    for ax in axes:
        if ax is None:
            continue
        ax.grid(True, color=grid_color)
        ax.set_axisbelow(True)
        ax.spines["bottom"].set_color("black")
        ax.spines["top"].set_color("black")
        ax.spines["right"].set_color("black")
        ax.spines["left"].set_color("black")
        ax.tick_params(axis="x", colors="black", which="both")
        ax.tick_params(axis="y", colors="black", which="both")
        ax.yaxis.label.set_color("black")
        ax.xaxis.label.set_color("black")
        ax.set_facecolor(back_color)


def create_layout():
    fig = plt.figure(figsize=(16, 9))
    fig.subplots_adjust(top=0.975, bottom=0.120, left=0.060, right=0.986, hspace=0.44, wspace=0.64)
    gs = GridSpec(18, 19, figure=fig)
    ax_logo = fig.add_subplot(gs[0:2, 0:4])
    ax_info_1 = fig.add_subplot(gs[0:2, 5:9])
    ax_info_2 = fig.add_subplot(gs[0:2, 10:14])
    ax_info_3 = fig.add_subplot(gs[0:2, 15:])
    ax0 = fig.add_subplot(gs[3:12, 0:4])
    ax00 = fig.add_subplot(gs[3:7, 5:9])
    ax01 = fig.add_subplot(gs[3:7, 10:14])
    ax02 = fig.add_subplot(gs[3:7, 15:])
    ax10 = fig.add_subplot(gs[8:12, 5:9])
    ax11 = fig.add_subplot(gs[8:12, 10:14])
    ax12 = fig.add_subplot(gs[8:12, 15:])
    ax30 = fig.add_subplot(gs[14:, 0:5])
    ax31 = fig.add_subplot(gs[14:, 6:12])
    ax32 = fig.add_subplot(gs[14:, 13:19])

    apply_grey_theme(fig, [ax0, ax00, ax01, ax02, ax10, ax11, ax12, ax30, ax31, ax32], transparent=True)
    return fig, ax_logo, ax_info_1, ax_info_2, ax_info_3, ax0, ax00, ax01, ax02, ax10, ax11, ax12, ax30, ax31, ax32


def plot_logo(ax):
    try:
        img = plt.imread(os.path.join("..", "..", "..", "assets", "logo_2.png"))
        ax.imshow(img)
        ax.text(910, 295, "Version 2.0.0")
    except Exception:
        pass
    ax.axis("off")


def show_metadata(ax_mesh, ax_model, ax_solver, operation_folder):
    try:
        with open(os.path.join(operation_folder, "log.txt")) as f:
            log = f.readlines()

        dh = 0.3
        h = 0.8

        for i, line in enumerate(log):
            if "| Mesh info:" in line:
                line_split = log[i + 4].split("|")
                n_elems = int(line_split[1])
                n_nodes = int(line_split[2])
                location = line_split[3].strip(" ")
                ax_mesh.text(0, h - 0 * dh, "Mesh info:", size=12, fontname="serif")
                ax_mesh.text(0, h - 1 * dh, f"- Location: {location}", size=10, fontname="serif")
                ax_mesh.text(0, h - 2 * dh, f"- Number of elems: {n_elems}", size=10, fontname="serif")
                ax_mesh.text(0, h - 3 * dh, f"- Number of nodes: {n_nodes}", size=10, fontname="serif")
                ax_mesh.axis("off")

            elif "| Solver info:" in line:
                line_split = log[i + 4].split("|")
                ksp_solver = line_split[1].strip(" ")
                ksp_pc = line_split[2].strip(" ")
                tol = float(line_split[3])
                max_ite = int(line_split[4])
                cpu_time = log[-1].strip(" ").strip("Total time: ")[:8] if len(log) > 0 else "N/A"
                ax_solver.text(0, h - 0 * dh, "Simulation Info: ", size=12, fontname="serif")
                ax_solver.text(0, h - 1 * dh, f"- Solver: {ksp_solver}, PC: {ksp_pc}", size=10, fontname="serif")
                ax_solver.text(0, h - 2 * dh, f"- tol: {tol}, max_ite: {max_ite}", size=10, fontname="serif")
                ax_solver.text(0, h - 3 * dh, f"- CPU Time: {cpu_time}", size=10, fontname="serif")
                ax_solver.axis("off")

            elif "| Constitutive model:" in line:
                elastic_elems = log[i + 4].split("|")[2].strip(" ") if i + 4 < len(log) else "N/A"
                nonelastic_elems = log[i + 5].split("|")[2].strip(" ") if i + 5 < len(log) else "N/A"
                thermoelastic_elems = log[i + 6].split("|")[2].strip(" ") if i + 6 < len(log) else "N/A"
                ax_model.text(0, h - 0 * dh, "Constitutive model: ", size=12, fontname="serif")
                ax_model.text(0, h - 1 * dh, f"- Elastic: {elastic_elems}", size=10, fontname="serif")
                ax_model.text(0, h - 2 * dh, f"- Non-elastic: {nonelastic_elems}", size=10, fontname="serif")
                ax_model.text(0, h - 3 * dh, f"- Thermoelastic: {thermoelastic_elems}", size=10, fontname="serif")
                ax_model.axis("off")
    except Exception:
        pass


# -------------------------
# Flat output readers (ROOT/case_*/...)
# -------------------------
def infer_cavern_type_from_case_name(case_name: str) -> str:
    parts = case_name.split("_")
    return parts[-1].lower() if len(parts) > 1 else case_name.lower()


def nice_cavern_label_from_cavern_type(cavern_type: str) -> str:
    low = cavern_type.lower()
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
    if low.startswith("irregular"):
        return "Irregular"
    return cavern_type


def path_u_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.xdmf")

def path_geom_msh(case_folder):
    return os.path.join(case_folder, "operation", "mesh", "geom.msh")

def path_p_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "p_elems", "p_elems.xdmf")

def path_q_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "q_elems", "q_elems.xdmf")

def path_pressure_json(case_folder):
    return os.path.join(case_folder, "pressure_schedule.json")


def read_case_pressure_scenario(case_path: str) -> str | None:
    pjson = path_pressure_json(case_path)
    if not os.path.isfile(pjson):
        return None
    try:
        with open(pjson, "r") as f:
            data = json.load(f)
    except Exception:
        return None

    if isinstance(data, dict) and "pressure_scenario" in data:
        v = data.get("pressure_scenario")
        return str(v).lower() if v is not None else None

    if isinstance(data, dict) and "scenario" in data:
        v = str(data.get("scenario")).lower()
        if v in ("sinus", "irregular", "csv_profile", "linear"):
            return v

    return None


def collect_cases_flat(ROOT: str, target_pressure: str | None):
    """
    ROOT/case_*/operation/...
    Filter by pressure scheme from pressure_schedule.json
    """
    cases = []
    for sub in sorted(os.listdir(ROOT)):
        if not sub.lower().startswith("case_"):
            continue

        case_path = os.path.join(ROOT, sub)
        if not os.path.isdir(case_path):
            continue

        needed = [
            path_u_xdmf(case_path),
            path_geom_msh(case_path),
            path_p_xdmf(case_path),
            path_q_xdmf(case_path),
        ]
        if not all(os.path.isfile(p) for p in needed):
            continue

        pres = read_case_pressure_scenario(case_path)
        if target_pressure is not None and (pres or "").lower() != target_pressure.lower():
            continue

        cavern_type = infer_cavern_type_from_case_name(sub)
        label = nice_cavern_label_from_cavern_type(cavern_type)

        cases.append({
            "label": label,
            "cavern_type": cavern_type,        # e.g. regular600
            "case_name": sub,
            "case_path": case_path,
            "operation_folder": os.path.join(case_path, "operation"),
            "pressure_scenario": pres,
        })
    return cases


def read_pressure_schedule(case_path: str):
    pjson = path_pressure_json(case_path)
    if not os.path.isfile(pjson):
        return None, None

    with open(pjson, "r") as f:
        data = json.load(f)

    if "t_hours" in data and "p_MPa" in data:
        t_hours = np.asarray(data["t_hours"], dtype=float)
        p_mpa = np.asarray(data["p_MPa"], dtype=float)
        return t_hours, p_mpa

    if "t_values_s" in data and "p_values_Pa" in data:
        t_hours = np.asarray(data["t_values_s"], dtype=float) / hour
        p_mpa = np.asarray(data["p_values_Pa"], dtype=float) / MPa
        return t_hours, p_mpa

    if "t_values" in data and "p_values" in data:
        t_hours = np.asarray(data["t_values"], dtype=float) / hour
        p_mpa = np.asarray(data["p_values"], dtype=float) / MPa
        return t_hours, p_mpa

    return None, None


# -------------------------
# WallProfile + probes
# -------------------------
class WallProfileData:
    def __init__(self, operation_folder, scale=1.0):
        points, self.time_list, u_field = post.read_node_vector(os.path.join(operation_folder, "u", "u.xdmf"))

        reader_msh = meshio.read(os.path.join(operation_folder, "mesh", "geom.msh"))
        points_msh = reader_msh.points

        wall_idx = None
        if hasattr(reader_msh, "cells_dict"):
            wall_idx = np.unique(reader_msh.cells_dict["line"].flatten())
        elif isinstance(reader_msh.cells, dict):
            if "line" in reader_msh.cells:
                wall_idx = np.unique(reader_msh.cells["line"].flatten())
        else:
            for cell_block in reader_msh.cells:
                if hasattr(cell_block, "type") and cell_block.type == "line":
                    wall_idx = np.unique(cell_block.data.flatten())
                    break
        if wall_idx is None:
            raise ValueError("No 'line' cells found in mesh (wall profile)")

        mapping = post.build_mapping(points_msh, points)
        wall_idx = np.asarray([mapping[int(i)] for i in wall_idx], dtype=int)

        self.wall_points = points[wall_idx]
        self.wall_u = u_field[:, wall_idx]

        sorted_idx = np.argsort(self.wall_points[:, 2])
        self.wall_points = self.wall_points[sorted_idx]
        self.wall_u = self.wall_u[:, sorted_idx]

        self.scale = scale

    def get_wall_coordinates(self, time_step: int):
        return self.wall_points + self.scale * self.wall_u[time_step, :]


def auto_generate_probes(wall_profile: WallProfileData, n_bend_probes: int):
    pts = wall_profile.wall_points
    x = pts[:, 0]
    z = pts[:, 2]

    idx_bottom = int(np.argmin(z))
    idx_top = int(np.argmax(z))
    z_mid = 0.5 * (z[idx_bottom] + z[idx_top])
    idx_mid = int(np.argmin(np.abs(z - z_mid)))

    probes_idx = [idx_top, idx_mid, idx_bottom]

    dx = np.gradient(x)
    dz_ = np.gradient(z)
    ddx = np.gradient(dx)
    ddz = np.gradient(dz_)

    denom = (dx * dx + dz_ * dz_) ** 1.5
    denom[denom == 0.0] = np.inf
    curvature = np.abs(dx * ddz - dz_ * ddx) / denom
    curvature[np.isnan(curvature)] = 0.0

    idx_all = np.argsort(curvature)[::-1]

    def too_close(new_i, existing_idx, min_gap=5):
        return any(abs(int(new_i) - int(ei)) < min_gap for ei in existing_idx)

    for i in idx_all:
        if len(probes_idx) >= 3 + int(n_bend_probes):
            break
        if int(i) in probes_idx:
            continue
        if too_close(int(i), probes_idx):
            continue
        probes_idx.append(int(i))

    probes_idx = sorted(probes_idx, key=lambda k: z[k], reverse=True)
    probes = pts[probes_idx]
    return probes


def plot_cavern_shape(ax, wall: WallProfileData, t_step: int):
    wall_t0 = wall.get_wall_coordinates(time_step=0)
    wall_tf = wall.get_wall_coordinates(time_step=t_step)

    ax.plot(wall_t0[:, 0], wall_t0[:, 2], "-", color="black", label="Initial shape")
    ax.plot(wall_tf[:, 0], wall_tf[:, 2], "-", color="steelblue", label="Final shape")
    ax.set_xlabel("x (m)", size=12, fontname="serif")
    ax.set_ylabel("z (m)", size=12, fontname="serif")
    ax.legend(loc=1, shadow=True, fancybox=True, prop={"size": 8})
    ax.axis("equal")


def plot_pressure_schedule(ax, case_path: str, t_step: int, n_steps: int):
    t_hours, p_mpa = read_pressure_schedule(case_path)
    if t_hours is None:
        ax.text(0.5, 0.5, "pressure_schedule.json not found", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return

    # Map t_step -> pressure index safely
    idx = min(int(t_step), len(t_hours) - 1)
    ax.plot(t_hours, p_mpa, "-", color="darkred", linewidth=0.8, label="Pressure schedule")
    stride = max(1, len(t_hours) // 200)
    ax.plot(t_hours[::stride], p_mpa[::stride], "o", color="darkred", markersize=2)
    ax.scatter(t_hours[idx], p_mpa[idx], c="white", edgecolors="black", zorder=10000)

    ax.set_xlabel("Time (hours)", size=12, fontname="serif")
    ax.set_ylabel("Cavern pressure (MPa)", size=12, fontname="serif")
    ax.legend(loc="best", shadow=True, fancybox=True, prop={"size": 8})
    ax.grid(True, alpha=0.3)


# -------------------------
# Stress paths
# -------------------------
class StressPath:
    def __init__(self, operation_folder, probes):
        self.probes = probes
        self.points, self.time_list, self.p_elems = post.read_cell_scalar(os.path.join(operation_folder, "p_elems", "p_elems.xdmf"))
        self.points, self.time_list, self.q_elems = post.read_cell_scalar(os.path.join(operation_folder, "q_elems", "q_elems.xdmf"))

        self.p_probes = np.zeros((probes.shape[0], self.time_list.size))
        self.q_probes = np.zeros((probes.shape[0], self.time_list.size))
        for i, probe in enumerate(probes):
            idx = post.find_closest_point(probe, self.points)
            self.p_probes[i, :] = -self.p_elems[:, idx] / MPa
            self.q_probes[i, :] = self.q_elems[:, idx] / MPa


def plot_dilatancy_boundaries_all(ax, p_min=0.01, p_max=40.0, npts=500, show_legend=True):
    p = np.linspace(p_min, p_max, npts)
    I1 = 3.0 * p

    def q_from_sqrtJ2(sqrtJ2):
        return np.sqrt(3.0) * sqrtJ2

    for D, lab, ls in [
        (0.27, "Ratigan 1991 (D=0.27)", "--"),
        (0.18, "Ratigan 1991 (D=0.18)", ":"),
    ]:
        sqrtJ2 = D * I1
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), linestyle=ls, linewidth=1.3, alpha=0.9, label=lab)

    D_sp, b_sp = 0.27, 1.9
    sqrtJ2 = D_sp * I1 + b_sp
    ax.plot(p, q_from_sqrtJ2(sqrtJ2), linestyle="-.", linewidth=1.3, alpha=0.9, label="Spiers 1988")

    def devries_q(I1_MPa, psi_rad, *, D1, D2, T0, m, sigma_ref=1.0):
        sgn = np.sign(I1_MPa)
        sgn[sgn == 0.0] = 1.0
        denom = (np.sqrt(3.0) * np.cos(psi_rad) - D2 * np.sin(psi_rad))
        sqrtJ2_ = (D1 * ((I1_MPa / (sgn * sigma_ref)) ** m) + T0) / denom
        return np.sqrt(3.0) * sqrtJ2_

    psi_c = -np.pi / 6.0
    psi_e = np.pi / 6.0

    D1, D2, T0, m = 0.683, 0.512, 1.50, 0.75
    ax.plot(p, devries_q(I1, psi_c, D1=D1, D2=D2, T0=T0, m=m), "-", linewidth=1.5, alpha=0.95, label="De Vries 2005 – comp")
    ax.plot(p, devries_q(I1, psi_e, D1=D1, D2=D2, T0=T0, m=m), "-", linewidth=1.5, alpha=0.95, label="De Vries 2005 – ext")

    if show_legend:
        ax.legend(loc="best", fontsize=8, frameon=True, fancybox=True)


COLORS = ["deepskyblue", "tomato", "orange", "steelblue", "purple", "magenta"]


def plot_stress_paths(ax, stress_path: StressPath, i: int, t_step: int, show_legend: bool):
    ax.plot(stress_path.p_probes[i], stress_path.q_probes[i], "-", color=COLORS[i % len(COLORS)], linewidth=2.0)
    ax.scatter(
        stress_path.p_probes[i, t_step],
        stress_path.q_probes[i, t_step],
        c="white",
        edgecolors="black",
        zorder=10000,
    )
    plot_dilatancy_boundaries_all(ax, p_min=0.01, p_max=40.0, npts=500, show_legend=show_legend)
    ax.set_xlabel("Mean stress (MPa)", size=10, fontname="serif")
    ax.set_ylabel("Von Mises stress (MPa)", size=10, fontname="serif")


def plot_probes(ax, probes):
    for i, probe in enumerate(probes):
        ax.scatter(probe[0], probe[2], c=COLORS[i % len(COLORS)], edgecolors="black", zorder=10000, s=100)


# =============================================================================
# Dashboard (static) - render to PNG/PDF at a chosen time step
# =============================================================================
def render_dashboard_for_case(case_dict, *, t_step: int = -1, n_bend_probes=None, offset_r=0.15):
    """
    Builds ONE dashboard figure for ONE case at ONE time step and SAVES it.
    No interactive widgets (Slider/Button removed).
    """
    case_path = case_dict["case_path"]
    operation_folder = case_dict["operation_folder"]
    label = case_dict["label"]
    case_name = case_dict["case_name"]

    fig, ax_logo, ax_info_1, ax_info_2, ax_info_3, ax0, ax00, ax01, ax02, ax10, ax11, ax12, ax30, ax31, ax32 = create_layout()
    fig.suptitle(f"{label} — {case_name}", fontsize=12)

    # Turn off unused axes (layout placeholders)
    ax31.axis("off")
    ax32.axis("off")

    # Load wall for shape + probes
    wall_profile = WallProfileData(operation_folder, scale=1.0)

    # Decide bend probes
    if n_bend_probes is None:
        n_bend_probes = 3 if GEOMETRY_TYPE == "irregular" else 1

    probes = auto_generate_probes(wall_profile, n_bend_probes=n_bend_probes)

    # Push interior probes a bit into salt (except top & bottom)
    n_probes = probes.shape[0]
    for k in range(n_probes):
        if 0 < k < n_probes - 1:
            probes[k, 0] += float(offset_r)

    stress_path = StressPath(operation_folder, probes)
    n_steps = int(stress_path.time_list.size)

    # Resolve t_step
    if t_step < 0:
        t_step_use = n_steps - 1
    else:
        t_step_use = min(int(t_step), n_steps - 1)

    plot_logo(ax_logo)
    show_metadata(ax_info_1, ax_info_2, ax_info_3, operation_folder)

    plot_cavern_shape(ax0, wall_profile, t_step=t_step_use)
    plot_probes(ax0, probes)

    plot_pressure_schedule(ax30, case_path, t_step=t_step_use, n_steps=n_steps)

    stress_axes = [ax00, ax01, ax02, ax10, ax11, ax12]
    n_axes = len(stress_axes)
    n_use = min(n_probes, n_axes)

    for i in range(n_use):
        plot_stress_paths(stress_axes[i], stress_path, i=i, t_step=t_step_use, show_legend=(i == 0))
    for j in range(n_use, n_axes):
        stress_axes[j].axis("off")

    apply_grey_theme(fig, [ax0, ax00, ax01, ax02, ax10, ax11, ax12, ax30], transparent=True)

    # Save
    out_dir = os.path.join(ROOT, "_dashboard_outputs")
    os.makedirs(out_dir, exist_ok=True)

    base = f"dashboard_{case_name}_t{t_step_use:05d}"
    saved = []

    if SAVE_PNG:
        p = os.path.join(out_dir, base + ".png")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        saved.append(p)

    if SAVE_PDF:
        p = os.path.join(out_dir, base + ".pdf")
        fig.savefig(p, bbox_inches="tight")
        saved.append(p)

    plt.close(fig)
    return saved


# =============================================================================
# Main: select cases and render one static dashboard per cavern label
# =============================================================================
def main():
    # 1) Collect matching cases (pressure from JSON)
    target_pressure = SELECT.get("pressure", None)
    cases = collect_cases_flat(ROOT, target_pressure)

    if not cases:
        raise RuntimeError(f"No cases found for pressure='{target_pressure}' under {ROOT}")

    # 2) Optional extra filters (caverns + substring)
    cav_filter = SELECT.get("caverns", None)
    contains = SELECT.get("case_contains", None)

    filtered = []
    for c in cases:
        if cav_filter is not None:
            cav_filter_l = [x.lower() for x in cav_filter]
            if c["cavern_type"].lower() not in cav_filter_l and c["label"].lower() not in cav_filter_l:
                continue
        if contains is not None and contains.lower() not in c["case_name"].lower():
            continue
        filtered.append(c)

    if not filtered:
        raise RuntimeError(f"No cases left after filters SELECT={SELECT}")

    # 3) Pick ONE case per cavern label (first match)
    by_label = {}
    for c in filtered:
        if c["label"] not in by_label:
            by_label[c["label"]] = c

    # 4) Render dashboards in preferred order
    saved_all = []
    for lab in CAVERN_ORDER:
        if lab in by_label:
            saved = render_dashboard_for_case(by_label[lab], t_step=-1)
            saved_all.extend(saved)
        else:
            # not an error; maybe your SELECT filtered it out
            pass

    # 5) Any extras not in CAVERN_ORDER
    extras = [lab for lab in by_label.keys() if lab not in CAVERN_ORDER]
    for lab in sorted(extras):
        saved = render_dashboard_for_case(by_label[lab], t_step=-1)
        saved_all.extend(saved)

    print("\n[OK] Saved dashboards:")
    for p in saved_all:
        print(" -", p)


if __name__ == "__main__":
    main()
