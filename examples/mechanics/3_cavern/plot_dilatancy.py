import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider, Button
import meshio

import safeincave.PostProcessingTools as post

# -------------------------
# Global constants
# -------------------------
GEOMETRY_TYPE = "irregular"  # "regular" / "irregular" / "tilted" (only affects default bend probes)
hour = 3600.0
day = 24.0 * hour
MPa = 1e6

# Which caverns to show (order)
CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt"]


# -------------------------
# UI / layout (unchanged)
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

    # We will turn off ax31/ax32 later (no convergence/subsidence), but keep them for layout.
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
    # Uses operation/log.txt
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
# Data readers (adapted to OutputNobian structure)
# -------------------------
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


def collect_cases_nested(ROOT: str, target_scheme: str):
    """
    ROOT/<GROUP>/<CASE>/operation/...
    Select based on CASE folder name (case_sinus... / case_irregular...)
    """
    cases = []
    for group in sorted(os.listdir(ROOT)):
        group_path = os.path.join(ROOT, group)
        if not os.path.isdir(group_path):
            continue
        if group.lower().startswith("pressure_"):
            continue

        label = cavern_label_from_group(group)

        for sub in sorted(os.listdir(group_path)):
            if not sub.lower().startswith("case_"):
                continue
            if scheme_from_case_folder(sub) != target_scheme:
                continue

            case_path = os.path.join(group_path, sub)
            if not os.path.isdir(case_path):
                continue

            # Required for stress dashboard
            needed = [
                path_u_xdmf(case_path),
                path_geom_msh(case_path),
                path_p_xdmf(case_path),
                path_q_xdmf(case_path),
            ]
            if not all(os.path.isfile(p) for p in needed):
                print(f"[SKIP] {group}/{sub} missing required output files")
                continue

            operation_folder = os.path.join(case_path, "operation")
            cases.append({
                "label": label,
                "group": group,
                "case_name": sub,
                "case_path": case_path,
                "operation_folder": operation_folder,
            })

    return cases


def read_pressure_schedule(case_path: str):
    pjson = path_pressure_json(case_path)
    if not os.path.isfile(pjson):
        return None, None
    with open(pjson, "r") as f:
        data = json.load(f)
    t_hours = np.asarray(data["t_values"], float) / hour
    p_mpa = np.asarray(data["p_values"], float) / MPa
    return t_hours, p_mpa


# -------------------------
# WallProfile + probes (unchanged logic)
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


def plot_pressure_schedule(ax, case_path: str):
    t_hours, p_mpa = read_pressure_schedule(case_path)
    if t_hours is None:
        ax.text(0.5, 0.5, "pressure_schedule.json not found", ha="center", va="center", transform=ax.transAxes)
        return None, None, None

    ax.plot(t_hours, p_mpa, "-", color="darkred", linewidth=0.8, label="Pressure schedule")
    ax.plot(t_hours[::50], p_mpa[::50], "o", color="darkred", markersize=2)
    marker = ax.scatter(t_hours[0], p_mpa[0], c="white", edgecolors="black", zorder=10000)

    ax.set_xlabel("Time (hours)", size=12, fontname="serif")
    ax.set_ylabel("Cavern pressure (MPa)", size=12, fontname="serif")
    ax.legend(loc="best", shadow=True, fancybox=True, prop={"size": 8})
    ax.grid(True, alpha=0.3)
    return t_hours, p_mpa, marker


# -------------------------
# Stress paths (unchanged)
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
        sqrtJ2_ = D1 * ((I1_MPa / (sgn * sigma_ref)) ** m) / denom + T0
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


# -------------------------
# NEW: one dashboard per cavern label
# -------------------------
def show_dashboard_for_case(case_dict, *, n_bend_probes=None, offset_r=0.15):
    """
    Builds ONE dashboard figure for ONE case (one cavern geometry).
    No convergence/subsidence. Keeps stress plots + shape + pressure.
    """
    case_path = case_dict["case_path"]
    operation_folder = case_dict["operation_folder"]
    label = case_dict["label"]
    case_name = case_dict["case_name"]

    fig, ax_logo, ax_info_1, ax_info_2, ax_info_3, ax0, ax00, ax01, ax02, ax10, ax11, ax12, ax30, ax31, ax32 = create_layout()
    fig.suptitle(f"{label} — {case_name}", fontsize=12)

    # Turn off the unused axes (we don't want convergence/subsidence)
    ax31.axis("off")
    ax32.axis("off")

    # Load wall for shape + probes
    wall_profile = WallProfileData(operation_folder, scale=1.0)

    # Decide bend probes
    if n_bend_probes is None:
        if GEOMETRY_TYPE == "irregular":
            n_bend_probes = 3
        else:
            n_bend_probes = 1

    probes = auto_generate_probes(wall_profile, n_bend_probes=n_bend_probes)

    # Push interior probes a bit into salt (except top & bottom)
    n_probes = probes.shape[0]
    for k in range(n_probes):
        if 0 < k < n_probes - 1:
            probes[k, 0] += float(offset_r)

    stress_path = StressPath(operation_folder, probes)
    n_steps = stress_path.time_list.size

    plot_logo(ax_logo)
    show_metadata(ax_info_1, ax_info_2, ax_info_3, operation_folder)

    plot_cavern_shape(ax0, wall_profile, t_step=-1)
    plot_probes(ax0, probes)

    t_press, p_press, press_marker = plot_pressure_schedule(ax30, case_path)

    stress_axes = [ax00, ax01, ax02, ax10, ax11, ax12]
    n_axes = len(stress_axes)
    n_use = min(n_probes, n_axes)

    for i in range(n_use):
        plot_stress_paths(stress_axes[i], stress_path, i=i, t_step=0, show_legend=(i == 0))
    for j in range(n_use, n_axes):
        stress_axes[j].axis("off")

    def update_plot(val):
        t_step = int(val)

        # update pressure marker
        if t_press is not None and p_press is not None and press_marker is not None:
            idx = min(t_step, len(t_press) - 1)
            press_marker.set_offsets([[t_press[idx], p_press[idx]]])

        # update stress plots
        for i in range(n_use):
            ax = stress_axes[i]
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            ax.cla()
            plot_stress_paths(ax, stress_path, i=i, t_step=t_step, show_legend=(i == 0))
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)

        apply_grey_theme(fig, [ax0, ax00, ax01, ax02, ax10, ax11, ax12, ax30], transparent=True)

    # Slider
    axtime = fig.add_axes([0.09, 0.02, 0.75, 0.03])
    time_slider = Slider(ax=axtime, label="Time step", valmin=0, valmax=n_steps - 1, valinit=0, valstep=1)
    time_slider.on_changed(update_plot)

    # Play / Pause
    playing = False
    timer = fig.canvas.new_timer(interval=12)
    PLAY_STEP = 50

    def advance_frame():
        nonlocal playing
        if not playing:
            return
        k = int(time_slider.val)
        if k < n_steps - 1:
            time_slider.set_val(min(k + PLAY_STEP, n_steps - 1))
            timer.start()
        else:
            playing = False

    timer.add_callback(advance_frame)

    def on_play(event):
        nonlocal playing
        if not playing:
            playing = True
            timer.start()

    def on_pause(event):
        nonlocal playing
        playing = False

    ax_play = fig.add_axes([0.88, 0.02, 0.05, 0.03])
    ax_pause = fig.add_axes([0.93, 0.02, 0.05, 0.03])
    btn_play = Button(ax_play, "Play")
    btn_pause = Button(ax_pause, "Pause")
    btn_play.on_clicked(on_play)
    btn_pause.on_clicked(on_pause)

    plt.show()


# -------------------------
# Main: choose scheme and show 6 figures
# -------------------------
def main():
    ROOT = r"/home/gvandenbrekel/SafeInCave/OutputNobian"

    # Choose your scheme here:
    PRESSURE_SCHEME = "sinus"   # "sinus" or "irregular"

    cases = collect_cases_nested(ROOT, PRESSURE_SCHEME)
    if not cases:
        raise RuntimeError(f"No cases found for scheme='{PRESSURE_SCHEME}' under {ROOT}")

    # Choose ONE case per cavern label (if multiple exist, just take the first)
    by_label = {}
    for c in cases:
        if c["label"] not in by_label:
            by_label[c["label"]] = c

    # Open up to 6 figures (in consistent order)
    for lab in CAVERN_ORDER:
        if lab in by_label:
            show_dashboard_for_case(by_label[lab])
        else:
            print(f"[WARN] No case found for cavern label '{lab}' in scheme '{PRESSURE_SCHEME}'")

    # If there are labels not in CAVERN_ORDER, still show them (optional)
    extras = [lab for lab in by_label.keys() if lab not in CAVERN_ORDER]
    for lab in sorted(extras):
        show_dashboard_for_case(by_label[lab])


if __name__ == "__main__":
    main()
