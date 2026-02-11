import os
import numpy as np
import matplotlib.pyplot as plt
import meshio


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                           USER CONFIGURATION                                  ║
# ╠══════════════════════════════════════════════════════════════════════════════╣
# ║  Plot 2x3 overview of cavern wall profiles with probe points.                 ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

GRIDS_ROOT = r"/data/home/gbrekel/SafeInCave_new/grids"

CAVERNS = [
    ("Regular",             "cavern_regular_1200_3D"),
    ("Reversed-circulation","cavern_reversedcirculation_1200_3D"),
    ("Direct-circulation",  "cavern_directcirculation_1200_3D"),
    ("Tube-failure",        "cavern_tubefailure_1200_3D"),
    ("Tilt",                "cavern_tilted_1200_3D"),
    ("Fast-leached",        "cavern_fastleached_1200_3D"),
]

OUT_DIR = os.path.join(GRIDS_ROOT, "..", "examples", "mechanics", "nobian", "Simulation", "output", "_figures")
SHOW = False
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

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    axes = axes.flatten()

    probe_order = ["top", "quarter", "mid", "threequarter", "bottom"]

    for i, (label, folder_name) in enumerate(CAVERNS):
        ax = axes[i]
        msh_path = os.path.join(GRIDS_ROOT, folder_name, "geom.msh")

        if not os.path.isfile(msh_path):
            ax.text(0.5, 0.5, f"Missing:\n{folder_name}", ha="center", va="center",
                    transform=ax.transAxes, fontsize=9)
            ax.set_title(label)
            continue

        wall_points = load_wall_points_from_msh(msh_path)
        probes = auto_generate_probes_from_wall_points(wall_points)

        # Plot wall profile
        ax.plot(wall_points[:, 0], wall_points[:, 2], 'k-', linewidth=1.5)

        # Plot probe points
        for pname in probe_order:
            pt = probes[pname]
            color = PROBE_COLORS[pname]
            ax.scatter(pt[0], pt[2], c=color, s=60, marker='o',
                       edgecolors='black', linewidths=0.5, zorder=5, label=pname)
            ax.annotate(pname, (pt[0], pt[2]),
                        textcoords="offset points", xytext=(8, 0),
                        fontsize=7, color=color, fontweight='bold')

        ax.set_title(label, fontsize=11, fontweight='bold')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('z (m)')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)

    # Shared legend from first subplot
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=5, fontsize=9,
               frameon=True, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle('Cavern Wall Profiles with Probe Locations', fontsize=13, fontweight='bold')
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])

    outpath = os.path.join(OUT_DIR, "cavern_profiles_overview.png")
    fig.savefig(outpath, dpi=DPI, bbox_inches='tight')
    print(f"[SAVED] {outpath}")

    if SHOW:
        plt.show()
    plt.close(fig)

    print("[DONE]")


if __name__ == "__main__":
    main()
