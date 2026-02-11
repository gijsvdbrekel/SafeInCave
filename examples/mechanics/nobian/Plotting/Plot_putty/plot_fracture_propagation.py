import os
import numpy as np
import matplotlib.pyplot as plt
import meshio

import safeincave.PostProcessingTools as post

from case_index import (
    detect_layout_and_collect_cases,
    filter_cases,
)


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                           USER CONFIGURATION                                  ║
# ╠══════════════════════════════════════════════════════════════════════════════╣
# ║  Visualize Factor of Safety (FOS) propagation into the rock mass from the     ║
# ║  cavern wall. Points are sampled at increasing radial distances from the      ║
# ║  wall to check if damage (FOS < 1) is spreading inward (fracture propagation).║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── OUTPUT FOLDER ──────────────────────────────────────────────────────────────
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/nobian/Simulation/output"

# ── CASE SELECTION ─────────────────────────────────────────────────────────────
SELECT = {
    "caverns": None,
    "pressure": None,
    "scenario": None,
    "n_cycles": None,
    "operation_days": None,
    "case_contains": None,
}

# ── PLOT OPTIONS ───────────────────────────────────────────────────────────────
ONE_CASE_PER_SERIES = True

# PLOT_MODE: "combined" or "separate"
PLOT_MODE = "separate"

# ── RADIAL SAMPLING ────────────────────────────────────────────────────────────
# RADIAL_DISTANCES: Distances from cavern wall into the rock (meters)
#   The first value (0.0) is on the wall, subsequent values go into the rock
RADIAL_DISTANCES = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]

# ── FOS THRESHOLD ─────────────────────────────────────────────────────────────
# FOS (Factor of Safety) = q_boundary / q
#   FOS < 1.0 means dilating (unsafe)
#   FOS > 1.0 means safe
# Uses De Vries 2005 criterion with proper Lode angle from stress tensor
FOS_THRESHOLD = 1.0

# ── OUTPUT SETTINGS ────────────────────────────────────────────────────────────
OUT_DIR = os.path.join(ROOT, "_figures")
SHOW = False
DPI = 180

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                        END OF USER CONFIGURATION                              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝


# ══════════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════
MPA = 1e6
HOUR = 3600.0
DAY = 24.0 * HOUR

# --- Color coding by cavern shape (consistent across all plot scripts) ---
CAVERN_COLORS = {
    "Asymmetric":          "#1f77b4",   # blue
    "Direct-circulation":  "#ff7f0e",   # orange
    "IrregularFine":       "#d62728",   # red
    "Regular":             "#9467bd",   # purple
    "Reversed-circulation":"#8c564b",   # brown
    "Tilt":                "#e377c2",   # pink
    "Fast-leached":        "#17becf",   # cyan
    "Tube-failure":        "#bcbd22",   # olive
}

SCENARIO_COLORS = {
    "disloc_old_only":        "#1f77b4",
    "disloc_new_only":        "#ff7f0e",
    "desai_only":             "#2ca02c",
    "full_minus_desai":       "#d62728",
    "full":                   "#9467bd",
    "full_minus_ps": "#d62728",   # red
    "md_only":                "#1f77b4",   # blue
    "md_steady_only":         "#bcbd22",
    "full_md":                "#e377c2",
    "interlayer":             "#7f7f7f",
    "nointerlayer":           "#8c564b",
    None:                     "#333333",
}

SCENARIO_LINESTYLES = {
    "disloc_old_only":        "-",
    "disloc_new_only":        "--",
    "desai_only":             "-.",
    "full_minus_desai":       ":",
    "full":                   "-",
    "full_minus_ps": "-",
    "md_only":                "-",
    "md_steady_only":         "--",
    "full_md":                "-.",
    "interlayer":             "-",
    "nointerlayer":           "--",
    None:                     "-",
}

# Probe location colors
PROBE_COLORS = {
    "top": "#e41a1c",
    "quarter": "#377eb8",
    "mid": "#4daf4a",
    "threequarter": "#984ea3",
    "bottom": "#ff7f00",
}


# ══════════════════════════════════════════════════════════════════════════════
# FOS COMPUTATION FUNCTIONS (from stress tensor with Lode angle)
# ══════════════════════════════════════════════════════════════════════════════

def compute_psi_from_sigma(sig_voigt):
    """
    Compute only the Lode angle psi from stress tensor.

    The Lode angle requires the third invariant J3 which cannot be recovered
    from p and q alone, so it must be computed directly from the stress tensor.

    sig_voigt: array of shape (..., 6) in Voigt notation [s11, s22, s33, s12, s13, s23]
               (compression positive convention)

    Returns:
        psi: Lode angle (radians, range [-pi/6, pi/6])
    """
    s11 = sig_voigt[..., 0]
    s22 = sig_voigt[..., 1]
    s33 = sig_voigt[..., 2]
    s12 = sig_voigt[..., 3]
    s13 = sig_voigt[..., 4]
    s23 = sig_voigt[..., 5]

    # Mean stress (compression positive)
    p = (s11 + s22 + s33) / 3.0

    # Deviatoric stress components
    dev11 = s11 - p
    dev22 = s22 - p
    dev33 = s33 - p

    # J2 = 0.5 * s_ij * s_ij
    J2 = 0.5 * (dev11**2 + dev22**2 + dev33**2) + s12**2 + s13**2 + s23**2
    J2 = np.maximum(J2, 1e-30)  # Avoid numerical issues

    # J3 = det(deviatoric stress tensor)
    J3 = (dev11 * (dev22 * dev33 - s23**2)
          - s12 * (s12 * dev33 - s23 * s13)
          + s13 * (s12 * s23 - dev22 * s13))

    # Lode angle from J2 and J3
    sqrtJ2 = np.sqrt(J2)
    arg = (3.0 * np.sqrt(3.0) / 2.0) * J3 / (sqrtJ2**3)
    arg = np.clip(arg, -1.0, 1.0)
    psi = (1.0 / 3.0) * np.arcsin(arg)

    return psi


def q_dil_devries(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    """
    De Vries (2005) dilatancy boundary with Lode angle dependence.

    Parameters:
        p_MPa: mean stress in MPa (compression positive)
        psi: Lode angle in radians
        D1, D2, m, T0, sigma_ref: De Vries parameters

    Returns:
        q_boundary in MPa
    """
    I1 = 3.0 * p_MPa

    # Denominator with Lode angle dependence
    denom = np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi)
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)

    # Numerator
    num = D1 * (np.abs(I1) / sigma_ref) ** m + T0

    sqrtJ2_dil = num / denom
    return np.sqrt(3.0) * sqrtJ2_dil


def compute_dilatancy_boundary_q(p_mpa, psi=None):
    """
    Compute the dilatancy boundary q value for given mean stress p and Lode angle.
    Uses De Vries 2005 criterion with proper Lode angle.

    If psi is None, uses compression branch (psi = -pi/6).
    """
    if psi is None:
        psi = -np.pi / 6.0  # Compression branch
    return q_dil_devries(p_mpa, psi)


def compute_FOS(p_mpa, q_mpa, psi=None, q_tol_MPa=1e-3):
    """
    Compute Factor of Safety = q_boundary / q.

    FOS < 1.0 means dilating (unsafe)
    FOS > 1.0 means safe

    Parameters:
        p_mpa: mean stress in MPa
        q_mpa: von Mises stress in MPa
        psi: Lode angle in radians (if None, uses compression branch)
        q_tol_MPa: minimum q value to avoid division by zero

    Returns:
        FOS array
    """
    q_boundary = compute_dilatancy_boundary_q(p_mpa, psi)

    # Handle near-zero q values
    q_safe = np.where(q_mpa < q_tol_MPa, q_tol_MPa, q_mpa)

    FOS = q_boundary / q_safe

    # For very small q, set FOS to large value (safe)
    FOS = np.where(q_mpa < q_tol_MPa, 100.0, FOS)

    return FOS


# ══════════════════════════════════════════════════════════════════════════════
# FILE PATH HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def path_u_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.xdmf")

def path_geom_msh(case_folder):
    return os.path.join(case_folder, "operation", "mesh", "geom.msh")

def path_sig_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "sig", "sig.xdmf")

def path_p_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "p_elems", "p_elems.xdmf")

def path_q_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "q_elems", "q_elems.xdmf")

def case_has_required_files(case_path: str) -> bool:
    return all(os.path.isfile(p) for p in [
        path_u_xdmf(case_path),
        path_geom_msh(case_path),
        path_sig_xdmf(case_path),
        path_p_xdmf(case_path),
        path_q_xdmf(case_path),
    ])


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


def load_wall_points(case_folder):
    """Load cavern wall points, sorted by z-coordinate."""
    u_xdmf = path_u_xdmf(case_folder)
    msh_path = path_geom_msh(case_folder)

    points, _, _ = post.read_node_vector(u_xdmf)
    points_msh, wall_idx_msh = get_wall_indices_from_msh(msh_path)

    mapping = post.build_mapping(points_msh, points)
    wall_idx = np.array([mapping[i] for i in wall_idx_msh], dtype=int)

    wall_points = points[wall_idx]
    order = np.argsort(wall_points[:, 2])
    return wall_points[order]


def compute_outward_normal_2d(wall_points):
    """
    Compute outward normal direction at each wall point.
    Assumes wall is in x-z plane (y=0 symmetry) and cavern is to the left.
    Returns unit vectors pointing INTO the rock (radially outward from cavern).
    """
    n_pts = len(wall_points)
    normals = np.zeros((n_pts, 3))

    for i in range(n_pts):
        # Get neighboring points for tangent calculation
        if i == 0:
            tangent = wall_points[1] - wall_points[0]
        elif i == n_pts - 1:
            tangent = wall_points[-1] - wall_points[-2]
        else:
            tangent = wall_points[i+1] - wall_points[i-1]

        # Normalize tangent (in x-z plane)
        tangent[1] = 0  # Ensure y=0
        t_len = np.sqrt(tangent[0]**2 + tangent[2]**2)
        if t_len > 1e-10:
            tangent = tangent / t_len

        # Normal is perpendicular to tangent in x-z plane
        # Rotate 90° clockwise to point outward (into the rock, away from cavern)
        # tangent = (tx, 0, tz) -> normal = (tz, 0, -tx) for outward
        # But we want to point INTO the rock (positive x direction generally)
        normal = np.array([tangent[2], 0.0, -tangent[0]])

        # Ensure normal points outward (positive x direction for typical cavern)
        # The cavern center is typically at x~0, so outward is +x
        if normal[0] < 0:
            normal = -normal

        normals[i] = normal

    return normals


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
        probes[name] = {"wall_idx": idx, "wall_point": pts[idx]}

    return probes


def generate_radial_sample_points(wall_points, probes, normals, radial_distances):
    """
    Generate sample points at radial distances from the wall for each probe.

    Returns dict: probe_name -> list of (distance, xyz_point)
    """
    sample_points = {}

    for probe_name, probe_info in probes.items():
        wall_idx = probe_info["wall_idx"]
        wall_pt = probe_info["wall_point"]
        normal = normals[wall_idx]

        points_at_distances = []
        for dist in radial_distances:
            pt = wall_pt + dist * normal
            points_at_distances.append((dist, pt))

        sample_points[probe_name] = points_at_distances

    return sample_points


# ══════════════════════════════════════════════════════════════════════════════
# STRESS DATA READING
# ══════════════════════════════════════════════════════════════════════════════

def tensor33_to_voigt6(sig33):
    """
    Convert (...,3,3) tensor to (...,6) Voigt: [s11,s22,s33,s12,s13,s23]
    """
    s11 = sig33[..., 0, 0]
    s22 = sig33[..., 1, 1]
    s33 = sig33[..., 2, 2]
    s12 = sig33[..., 0, 1]
    s13 = sig33[..., 0, 2]
    s23 = sig33[..., 1, 2]
    return np.stack([s11, s22, s33, s12, s13, s23], axis=-1)


def read_stress_at_points(case_folder, sample_points):
    """
    Read stress state at all sample points over time.

    Uses STORED p_elems and q_elems (smoothed values written during simulation)
    for consistency with the solver output. The Lode angle psi is computed
    directly from the stress tensor since it requires J3 which cannot be
    recovered from p and q alone.

    Returns:
        time_days: array of time values in days
        stress_data: dict[probe_name] -> dict[distance] -> (p_array, q_array, psi_array)
    """
    # Read stored p_elems and q_elems (smoothed values from simulation)
    p_path = path_p_xdmf(case_folder)
    q_path = path_q_xdmf(case_folder)

    centroids_p, time_list_p, p_elems = post.read_cell_scalar(p_path)  # (nT, nC)
    centroids_q, time_list_q, q_elems = post.read_cell_scalar(q_path)  # (nT, nC)

    # Read stress tensor for Lode angle computation
    sig_path = path_sig_xdmf(case_folder)
    centroids_sig, time_list_sig, sig33 = post.read_cell_tensor(sig_path)  # (nT, nC, 3, 3)

    # Synchronize time steps (use the minimum common length)
    n_times = min(len(time_list_p), len(time_list_q), len(time_list_sig))
    time_list = time_list_p[:n_times]
    p_elems = p_elems[:n_times]
    q_elems = q_elems[:n_times]
    sig33 = sig33[:n_times]

    time_days = np.array(time_list) / DAY

    # Convert p to compression-positive convention and to MPa:
    # SafeInCave stores tension-positive => negate
    p_elems_MPa = -p_elems / MPA
    q_elems_MPa = q_elems / MPA  # q is always positive (magnitude)

    # Convert stress tensor to compression-positive for psi calculation
    sig33_MPa = -sig33 / MPA

    stress_data = {}
    for probe_name, dist_points in sample_points.items():
        stress_data[probe_name] = {}
        for dist, xyz in dist_points:
            # Find closest cell for each field (should be same cell, but be robust)
            idx_p = post.find_closest_point(xyz, centroids_p)
            idx_sig = post.find_closest_point(xyz, centroids_sig)

            # Get stored p and q values (smoothed)
            p = p_elems_MPa[:, idx_p]
            q = q_elems_MPa[:, idx_p]

            # Compute psi from stress tensor (requires J3)
            sig_point33 = sig33_MPa[:, idx_sig, :, :]
            sig_point6 = tensor33_to_voigt6(sig_point33)
            psi = compute_psi_from_sigma(sig_point6)

            stress_data[probe_name][dist] = (p, q, psi)

    return time_days, stress_data



def compute_FOS_data(time_days, stress_data):
    """
    Compute Factor of Safety for all points over time.

    Returns:
        fos_data: dict[probe_name] -> dict[distance] -> FOS_array
    """
    fos_data = {}
    for probe_name, dist_data in stress_data.items():
        fos_data[probe_name] = {}
        for dist, (p, q, psi) in dist_data.items():
            fos = compute_FOS(p, q, psi)
            fos_data[probe_name][dist] = fos

    return fos_data


# ══════════════════════════════════════════════════════════════════════════════
# PLOTTING
# ══════════════════════════════════════════════════════════════════════════════

def plot_cavern_with_probes(ax, wall_points, sample_points):
    """Plot cavern wall profile with probe locations and radial sample lines."""
    # Plot wall profile
    ax.plot(wall_points[:, 0], wall_points[:, 2], 'k-', linewidth=2, label='Cavern wall')

    # Plot each probe and its radial sample points
    for probe_name, dist_points in sample_points.items():
        color = PROBE_COLORS.get(probe_name, 'gray')

        # Extract coordinates
        xs = [pt[0] for _, pt in dist_points]
        zs = [pt[2] for _, pt in dist_points]

        # Plot line connecting radial points
        ax.plot(xs, zs, '-', color=color, linewidth=1.5, alpha=0.7)

        # Plot points
        for i, (dist, pt) in enumerate(dist_points):
            marker = 'o' if dist == 0 else 's'
            size = 80 if dist == 0 else 40
            ax.scatter(pt[0], pt[2], c=color, s=size, marker=marker,
                      edgecolors='black', linewidths=0.5, zorder=5)

        # Label at wall point
        wall_pt = dist_points[0][1]
        ax.annotate(probe_name, (wall_pt[0], wall_pt[2]),
                   textcoords="offset points", xytext=(10, 0),
                   fontsize=9, color=color, fontweight='bold')

    ax.set_xlabel('Radial distance (m)')
    ax.set_ylabel('Height z (m)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_title('Cavern profile with sample points')


def plot_FOS_time_evolution(ax, time_days, fos_data, probe_name, radial_distances):
    """Plot Factor of Safety evolution over time for one probe location."""
    color = PROBE_COLORS.get(probe_name, 'gray')

    for dist in radial_distances:
        fos = fos_data[probe_name][dist]
        alpha = 1.0 if dist == 0 else 0.4 + 0.6 * (1 - dist / max(radial_distances))
        linestyle = '-' if dist == 0 else '--'
        linewidth = 2.0 if dist == 0 else 1.2

        label = f'{dist:.1f}m' if dist > 0 else 'wall'
        ax.plot(time_days, fos, linestyle=linestyle, linewidth=linewidth,
               color=color, alpha=alpha, label=label)

    # FOS threshold line (FOS < 1 = dilating/unsafe)
    ax.axhline(y=1.0, color='red', linestyle=':', linewidth=1.5, alpha=0.8, label='FOS = 1')

    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Factor of Safety (q_boundary/q)')
    ax.set_title(f'{probe_name}')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=7, ncol=2)


def plot_propagation_depth(ax, time_days, fos_data, radial_distances, threshold=1.0):
    """
    Plot the maximum depth at which FOS < threshold (dilating) over time.
    This shows if/how far damage is propagating into the rock.
    """
    for probe_name in fos_data.keys():
        color = PROBE_COLORS.get(probe_name, 'gray')

        max_depth = np.zeros(len(time_days))
        for t_idx in range(len(time_days)):
            depth = 0.0
            for dist in radial_distances:
                fos = fos_data[probe_name][dist][t_idx]
                if fos < threshold:  # FOS < 1 means dilating
                    depth = dist
            max_depth[t_idx] = depth

        ax.plot(time_days, max_depth, linewidth=2, color=color, label=probe_name)

    ax.set_xlabel('Time (days)')
    ax.set_ylabel('FOS<1 penetration depth (m)')
    ax.set_title('Connected dilating zone size')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=9)


def plot_dilating_points_count(ax, time_days, fos_data, radial_distances, threshold=1.0):
    """
    Plot the number of points with FOS < threshold (dilating) over time.
    Grouped by probe location.
    """
    total_points = len(radial_distances)

    for probe_name in fos_data.keys():
        color = PROBE_COLORS.get(probe_name, 'gray')

        count = np.zeros(len(time_days))
        for t_idx in range(len(time_days)):
            n_dilating = 0
            for dist in radial_distances:
                if fos_data[probe_name][dist][t_idx] < threshold:  # FOS < 1 = dilating
                    n_dilating += 1
            count[t_idx] = n_dilating

        ax.plot(time_days, count, linewidth=2, color=color, label=probe_name)

    ax.set_xlabel('Time (days)')
    ax.set_ylabel(f'Number of dilating points (of {total_points})')
    ax.set_title('Connected dilating zone size (FOS<1)')
    ax.set_ylim(0, total_points + 0.5)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=9)


def get_case_color_and_style(cavern_label, scenario_preset):
    """Get color and linestyle for a case. Uses cavern color when no scenario."""
    if scenario_preset is not None:
        color = SCENARIO_COLORS.get(scenario_preset, "#333333")
    else:
        color = CAVERN_COLORS.get(cavern_label, "#333333")
    linestyle = SCENARIO_LINESTYLES.get(scenario_preset, "-")
    return color, linestyle


def pick_one_case_per_series(cases_meta):
    out = {}
    for m in cases_meta:
        key = (m.get("cavern_label"), m.get("scenario_preset"), m.get("pressure_scenario"))
        if key not in out:
            out[key] = m
    return list(out.values())


def plot_single_case(case_meta, wall_points, sample_points, time_days, stress_data, fos_data):
    """Create visualization for a single case: cavern profile + dilating zone depth."""
    cav = case_meta.get("cavern_label")
    sc = case_meta.get("scenario_preset")
    case_label = f"{cav} | {sc}" if sc else cav

    fig, (ax_cavern, ax_depth) = plt.subplots(1, 2, figsize=(16, 8))

    # Left: cavern profile with probe sample points
    plot_cavern_with_probes(ax_cavern, wall_points, sample_points)

    # Right: connected dilating zone size over time
    plot_propagation_depth(ax_depth, time_days, fos_data, RADIAL_DISTANCES, FOS_THRESHOLD)

    fig.suptitle(f'Dilatancy Zone Analysis: {case_label}\n'
                 f'(De Vries criterion, radial distances: {RADIAL_DISTANCES} m)',
                 fontsize=12, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    return fig


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=" * 70)
    print("DILATANCY ZONE ANALYSIS (Factor of Safety)")
    print("=" * 70)
    print(f"  ROOT:              {ROOT}")
    print(f"  Pressure:          {SELECT.get('pressure')}")
    print(f"  Scenario:          {SELECT.get('scenario')}")
    print(f"  Radial distances:  {RADIAL_DISTANCES} m")
    print(f"  FOS threshold:     {FOS_THRESHOLD} (FOS < threshold = dilating)")
    print(f"  Criterion:         De Vries 2005 with Lode angle")
    print(f"  p, q source:       Stored p_elems/q_elems (smoothed)")
    print(f"  psi source:        Computed from stress tensor (J3)")
    print(f"  Plot mode:         {PLOT_MODE}")
    print("=" * 70)

    all_cases = detect_layout_and_collect_cases(ROOT)
    all_cases = [m for m in all_cases if case_has_required_files(m["case_path"])]

    cases_meta = filter_cases(all_cases, SELECT)
    if not cases_meta:
        print("[DEBUG] Examples of found cases (first 15):")
        for m in all_cases[:15]:
            print(" -", m.get("cavern_label"), m.get("scenario_preset"),
                  m.get("pressure_scenario"), m.get("case_name"))
        raise RuntimeError(f"No cases matched SELECT={SELECT}")

    if ONE_CASE_PER_SERIES:
        cases_meta = pick_one_case_per_series(cases_meta)

    print(f"[INFO] Found {len(cases_meta)} case(s) to analyze:")
    for m in cases_meta:
        print(f"  - {m.get('cavern_label')} | {m.get('scenario_preset')} | {m.get('case_name')}")

    for case_meta in cases_meta:
        folder = case_meta["case_path"]
        cav = case_meta.get("cavern_label")
        sc = case_meta.get("scenario_preset")

        print(f"\n[PROCESSING] {cav} | {sc}")

        try:
            # Load geometry
            wall_points = load_wall_points(folder)
            normals = compute_outward_normal_2d(wall_points)
            probes = auto_generate_probes_from_wall_points(wall_points)

            # Generate radial sample points
            sample_points = generate_radial_sample_points(
                wall_points, probes, normals, RADIAL_DISTANCES
            )

            # Read stress data (p, q from stored fields; psi from stress tensor)
            time_days, stress_data = read_stress_at_points(folder, sample_points)

            # Compute Factor of Safety
            fos_data = compute_FOS_data(time_days, stress_data)

            # Create plot
            fig = plot_single_case(case_meta, wall_points, sample_points,
                                   time_days, stress_data, fos_data)

            # Save
            safe_name = f"{cav}_{sc}".replace(" ", "_").replace("/", "_").replace("None", "none")
            outname = f"fracture_propagation_FOS_{safe_name}.png"
            outpath = os.path.join(OUT_DIR, outname)
            fig.savefig(outpath, dpi=DPI)
            print(f"[SAVED] {outpath}")

            if SHOW:
                plt.show()
            plt.close(fig)

        except Exception as e:
            print(f"[ERROR] Failed to process: {e}")
            import traceback
            traceback.print_exc()

    print("\n[DONE]")


if __name__ == "__main__":
    main()
