import os
import numpy as np
import matplotlib.pyplot as plt
import meshio


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                           USER CONFIGURATION                                  ║
# ╠══════════════════════════════════════════════════════════════════════════════╣
# ║  Plot 2x3 overview of cavern wall profiles with probe points.                 ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GRIDS_ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", "..", "..", "..", "grids"))

CAVERNS = [
    ("Regular",             "cavern_regular_1200_3D"),
    ("Reversed-circulation","cavern_reversedcirculation_1200_3D"),
    ("Direct-circulation",  "cavern_directcirculation_1200_3D"),
    ("Tube-failure",        "cavern_tubefailure_1200_3D"),
    ("Tilt",                "cavern_tilted_1200_3D"),
    ("Fast-leached",        "cavern_fastleached_1200_3D"),
]

OUT_DIR = os.path.join(GRIDS_ROOT, "..", "examples", "mechanics", "nobian", "Simulation", "output", "_figures")
SHOW = True
DPI = 180

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                        END OF USER CONFIGURATION                              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

PROBE_COLORS = {
    "top":          "#e41a1c",
    "quarter":      "#377eb8",
    "mid":          "#4daf4a",
    "threequarter": "#984ea3",
    "bottom":       "#ff7f00",
}


# ══════════════════════════════════════════════════════════════════════════════
# GEOMETRY HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def get_wall_indices_from_msh(msh_path):
    """Extract wall node indices from mesh file."""
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


def load_wall_points_from_msh(msh_path):
    """Load cavern wall points directly from mesh file, sorted by z-coordinate."""
    points, wall_idx = get_wall_indices_from_msh(msh_path)
    wall_points = points[wall_idx]
    order = np.argsort(wall_points[:, 2])
    return wall_points[order]


def auto_generate_probes_from_wall_points(wall_points_sorted_z):
    """Generate probe locations at evenly-spaced heights along the wall."""
    pts = wall_points_sorted_z
    z = pts[:, 2]
    z_min, z_max = z.min(), z.max()

    fractions = {"bottom": 0.0, "threequarter": 0.25, "mid": 0.5, "quarter": 0.75, "top": 1.0}

    probes = {}
    for name, frac in fractions.items():
        z_target = z_min + frac * (z_max - z_min)
        idx = int(np.argmin(np.abs(z - z_target)))
        probes[name] = pts[idx]

    return probes


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    probe_order = ["top", "quarter", "mid", "threequarter", "bottom"]

    # First pass: load all wall profiles to determine shared z-limits
    wall_data = []
    for label, folder_name in CAVERNS:
        msh_path = os.path.join(GRIDS_ROOT, folder_name, "geom.msh")
        if os.path.isfile(msh_path):
            wp = load_wall_points_from_msh(msh_path)
            probes = auto_generate_probes_from_wall_points(wp)
            wall_data.append((label, wp, probes))
        else:
            wall_data.append((label, None, None))

    # Compute shared z-limits across all caverns for consistent vertical axes
    all_z_min = min(wp[:, 2].min() for _, wp, _ in wall_data if wp is not None)
    all_z_max = max(wp[:, 2].max() for _, wp, _ in wall_data if wp is not None)
    z_pad = 0.03 * (all_z_max - all_z_min)
    z_range = (all_z_max + z_pad) - (all_z_min - z_pad)

    # Find the max x-range across all caverns to set uniform subplot width
    max_x_range = 0
    for _, wp, _ in wall_data:
        if wp is not None:
            xr = wp[:, 0].max() - wp[:, 0].min()
            if xr > max_x_range:
                max_x_range = xr
    x_pad = 0.06 * max_x_range
    x_total = max_x_range + 2 * x_pad

    # Figure size: 3 columns, 2 rows, matching the data aspect ratio
    col_width = 3.0  # inches per column
    row_height = col_width * (z_range / x_total)
    fig_w = col_width * 3 + 1.5  # extra for labels/margins
    fig_h = row_height * 2 + 2.0  # extra for legend/titles/labels

    fig, axes = plt.subplots(2, 3, figsize=(fig_w, fig_h))
    axes = axes.flatten()

    for i, (label, wp, probes) in enumerate(wall_data):
        ax = axes[i]

        if wp is None:
            ax.text(0.5, 0.5, f"Missing", ha="center", va="center",
                    transform=ax.transAxes, fontsize=9)
            ax.set_title(label, fontsize=10, fontweight='bold')
            continue

        # Plot wall profile
        ax.plot(wp[:, 0], wp[:, 2], 'k-', linewidth=1.5)

        # Plot probe points
        for pname in probe_order:
            pt = probes[pname]
            color = PROBE_COLORS[pname]
            ax.scatter(pt[0], pt[2], c=color, s=50, marker='o',
                       edgecolors='black', linewidths=0.5, zorder=5, label=pname)

        # Tight x-limits centered on profile, uniform width
        x_center = 0.5 * (wp[:, 0].min() + wp[:, 0].max())
        ax.set_xlim(x_center - x_total / 2, x_center + x_total / 2)

        # Shared z-limits
        ax.set_ylim(all_z_min - z_pad, all_z_max + z_pad)

        ax.set_title(label, fontsize=10, fontweight='bold')
        ax.set_xlabel('x (m)', fontsize=8)
        ax.set_ylabel('z (m)', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

    # Shared legend on top of figure
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=5, fontsize=9,
           frameon=True, bbox_to_anchor=(0.5, 0.965))

    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.04, top=0.93,
                    wspace=0.25, hspace=0.25)

    outpath = os.path.join(OUT_DIR, "cavern_profiles_overview.png")
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight')
    print(f"[SAVED] {outpath}")

    if SHOW:
        plt.show()
    plt.close(fig)

    print("[DONE]")


if __name__ == "__main__":
    main()
