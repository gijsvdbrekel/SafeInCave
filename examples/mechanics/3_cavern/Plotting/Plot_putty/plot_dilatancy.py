#!/usr/bin/env python3
import os
import json
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import meshio

import safeincave.PostProcessingTools as post
from case_index import detect_layout_and_collect_cases, filter_cases, one_case_per_cavern_label


# -------------------------
# Global constants
# -------------------------
hour = 3600.0
day = 24.0 * hour
MPa = 1e6

CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt", "IrregularFine"]

# =============================================================================
# USER SELECTION
# =============================================================================
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/3_cavern/output"

SELECT = {
    "pressure": "sinus",
    "scenario": None,        # <-- NOW POSSIBLE
    "caverns": ["Irregular"],  # or None
    "n_cycles": None,
    "operation_days": None,
    "case_contains": None,
    "one_case_per_cavern": True,

    # NEW: which timestep to render in the dashboard
    # "last" or an int (0..Nt-1)
    "t_step": "last",
}

OUTDIR = os.path.join(ROOT, "_plots_dashboard")
DPI = 220


# -------------------------
# layout helpers
# -------------------------
def apply_theme(fig, axes):
    fig.patch.set_facecolor("white")
    for ax in axes:
        if ax is None:
            continue
        ax.grid(True, alpha=0.25)
        ax.set_axisbelow(True)

def create_layout():
    fig = plt.figure(figsize=(16, 9))
    fig.subplots_adjust(top=0.95, bottom=0.08, left=0.06, right=0.98, hspace=0.44, wspace=0.64)
    gs = GridSpec(18, 19, figure=fig)

    ax_info = fig.add_subplot(gs[0:2, 0:19])

    ax_shape = fig.add_subplot(gs[3:12, 0:4])

    ax00 = fig.add_subplot(gs[3:7, 5:9])
    ax01 = fig.add_subplot(gs[3:7, 10:14])
    ax02 = fig.add_subplot(gs[3:7, 15:])

    ax10 = fig.add_subplot(gs[8:12, 5:9])
    ax11 = fig.add_subplot(gs[8:12, 10:14])
    ax12 = fig.add_subplot(gs[8:12, 15:])

    ax_pressure = fig.add_subplot(gs[14:, 0:7])
    ax_empty1 = fig.add_subplot(gs[14:, 8:13])
    ax_empty2 = fig.add_subplot(gs[14:, 14:19])
    ax_empty1.axis("off")
    ax_empty2.axis("off")

    apply_theme(fig, [ax_info, ax_shape, ax00, ax01, ax02, ax10, ax11, ax12, ax_pressure])
    return fig, ax_info, ax_shape, [ax00, ax01, ax02, ax10, ax11, ax12], ax_pressure


def path_pressure_json(case_folder):
    return os.path.join(case_folder, "pressure_schedule.json")

def read_pressure_schedule(case_path: str):
    pjson = path_pressure_json(case_path)
    if not os.path.isfile(pjson):
        return None, None
    with open(pjson, "r") as f:
        data = json.load(f)

    if "t_hours" in data and "p_MPa" in data:
        return np.asarray(data["t_hours"], float), np.asarray(data["p_MPa"], float)
    if "t_values_s" in data and "p_values_Pa" in data:
        return np.asarray(data["t_values_s"], float)/hour, np.asarray(data["p_values_Pa"], float)/MPa
    if "t_values" in data and "p_values" in data:
        return np.asarray(data["t_values"], float)/hour, np.asarray(data["p_values"], float)/MPa
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
        if hasattr(reader_msh, "cells_dict") and "line" in reader_msh.cells_dict:
            wall_idx = np.unique(reader_msh.cells_dict["line"].flatten())
        else:
            for cell_block in reader_msh.cells:
                if getattr(cell_block, "type", None) == "line":
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


def auto_generate_probes(wall_profile: WallProfileData, n_bend_probes: int = 2, min_gap=5):
    pts = wall_profile.wall_points
    x = pts[:, 0]
    z = pts[:, 2]

    idx_bottom = int(np.argmin(z))
    idx_top = int(np.argmax(z))
    z_mid = 0.5 * (z[idx_bottom] + z[idx_top])
    idx_mid = int(np.argmin(np.abs(z - z_mid)))

    base = [idx_top, idx_mid, idx_bottom]

    dx = np.gradient(x)
    dz = np.gradient(z)
    ddx = np.gradient(dx)
    ddz = np.gradient(dz)

    denom = (dx*dx + dz*dz)**1.5
    denom[denom == 0.0] = np.inf
    curvature = np.abs(dx*ddz - dz*ddx) / denom
    curvature[np.isnan(curvature)] = 0.0

    idx_all = np.argsort(curvature)[::-1]

    def too_close(i, existing):
        return any(abs(int(i) - int(e)) < min_gap for e in existing)

    bend = []
    for i in idx_all:
        if len(bend) >= n_bend_probes:
            break
        if int(i) in base:
            continue
        if too_close(i, base + bend):
            continue
        bend.append(int(i))

    pick = sorted(base + bend, key=lambda k: z[k], reverse=True)
    return pts[pick]


# -------------------------
# Stress paths
# -------------------------
class StressPath:
    def __init__(self, operation_folder, probes):
        self.probes = probes
        self.points, self.time_list, self.p_elems = post.read_cell_scalar(os.path.join(operation_folder, "p_elems", "p_elems.xdmf"))
        _, _, self.q_elems = post.read_cell_scalar(os.path.join(operation_folder, "q_elems", "q_elems.xdmf"))

        self.p_probes = np.zeros((probes.shape[0], self.time_list.size))
        self.q_probes = np.zeros((probes.shape[0], self.time_list.size))
        for i, probe in enumerate(probes):
            idx = post.find_closest_point(probe, self.points)
            self.p_probes[i, :] = -self.p_elems[:, idx] / MPa
            self.q_probes[i, :] = self.q_elems[:, idx] / MPa


def plot_shape(ax, wall: WallProfileData, t_step: int, probes):
    w0 = wall.get_wall_coordinates(0)
    wt = wall.get_wall_coordinates(t_step)
    ax.plot(w0[:, 0], w0[:, 2], "-", linewidth=2, label="Initial")
    ax.plot(wt[:, 0], wt[:, 2], "-", linewidth=2, label=f"t_step={t_step}")
    for i, pr in enumerate(probes):
        ax.scatter(pr[0], pr[2], s=60, edgecolors="black")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.axis("equal")
    ax.legend(fontsize=8, frameon=True)


def plot_pressure(ax, case_path: str, t_step: int, n_steps: int):
    tH, pMPa = read_pressure_schedule(case_path)
    if tH is None:
        ax.text(0.5, 0.5, "No pressure_schedule.json", ha="center", va="center", transform=ax.transAxes)
        return
    ax.plot(tH, pMPa, linewidth=2)
    idx = min(int(t_step), len(tH)-1)
    ax.scatter([tH[idx]], [pMPa[idx]], s=50, edgecolors="black", zorder=10)
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Pressure (MPa)")


def plot_stress(ax, sp: StressPath, i_probe: int, t_step: int):
    ax.plot(sp.p_probes[i_probe], sp.q_probes[i_probe], linewidth=2)
    ax.scatter([sp.p_probes[i_probe, t_step]], [sp.q_probes[i_probe, t_step]], s=40, edgecolors="black", zorder=10)
    ax.set_xlabel("p (MPa)")
    ax.set_ylabel("q (MPa)")
    ax.grid(True, alpha=0.25)


def render_dashboard(case_meta: dict):
    case_path = case_meta["case_path"]
    op = os.path.join(case_path, "operation")

    fig, ax_info, ax_shape, stress_axes, ax_pressure = create_layout()

    wall = WallProfileData(op, scale=1.0)
    Nt = len(wall.time_list)

    tsel = SELECT.get("t_step", "last")
    if isinstance(tsel, str) and tsel.lower() == "last":
        t_step = Nt - 1
    else:
        t_step = int(tsel)
        t_step = max(0, min(t_step, Nt - 1))

    probes = auto_generate_probes(wall, n_bend_probes=2)
    sp = StressPath(op, probes)

    ax_info.axis("off")
    ax_info.text(
        0.01, 0.65,
        f"Case: {case_meta['case_name']} | cavern={case_meta.get('cavern_label')} | "
        f"pressure={case_meta.get('pressure_scenario')} | scenario={case_meta.get('scenario_preset')} | "
        f"t_step={t_step}/{Nt-1}",
        fontsize=11
    )

    plot_shape(ax_shape, wall, t_step, probes)
    plot_pressure(ax_pressure, case_path, t_step, Nt)

    n_use = min(len(stress_axes), probes.shape[0])
    for i in range(n_use):
        plot_stress(stress_axes[i], sp, i, t_step)
        stress_axes[i].set_title(f"Probe {i}", fontsize=9)
    for j in range(n_use, len(stress_axes)):
        stress_axes[j].axis("off")

    os.makedirs(OUTDIR, exist_ok=True)
    outname = f"dashboard_{case_meta.get('cavern_label','Unknown')}_{case_meta['case_name']}.png"
    outpath = os.path.join(OUTDIR, outname)
    fig.savefig(outpath, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("[OK] Saved:", outpath)


def main():
    all_cases = detect_layout_and_collect_cases(ROOT)
    selected = filter_cases(all_cases, SELECT)
    if not selected:
        print("[INFO] Found cases (unfiltered):")
        for m in all_cases[:40]:
            print(" -", m["case_name"], "| pressure:", m.get("pressure_scenario"),
                  "| scenario:", m.get("scenario_preset"), "| cavern:", m.get("cavern_label"))
        raise RuntimeError(f"No cases match SELECT={SELECT}")

    if SELECT.get("one_case_per_cavern", True):
        selected = one_case_per_cavern_label(selected)

    # order by preferred cavern order
    by_label = {c["cavern_label"]: c for c in selected}
    for lab in CAVERN_ORDER:
        if lab in by_label:
            render_dashboard(by_label[lab])

    # extras
    extras = [k for k in by_label.keys() if k not in CAVERN_ORDER]
    for lab in sorted(extras):
        render_dashboard(by_label[lab])


if __name__ == "__main__":
    main()

