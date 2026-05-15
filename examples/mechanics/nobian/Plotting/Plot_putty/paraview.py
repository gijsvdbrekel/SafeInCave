#!/usr/bin/env python3
"""
Post-processing: extract Paraview visualizations from simulation output.

Generates multiple XDMF files for Paraview visualization:

  1. cavern_surface.xdmf  — cavern wall colored by stress, FOS, displacement, creep
  2. cross_section.xdmf   — vertical slice with principal stress arrows
  3. convergence.xdmf     — displacement vectors on cross-section (wall movement)

Usage:
    python paraview.py
"""

import os
import re
import sys
import argparse
import numpy as np
import meshio
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "..", "..", ".."))
sys.path.insert(0, PROJECT_ROOT)

from case_index import detect_layout_and_collect_cases, filter_cases

# =============================================================================
# USER CONFIGURATION — edit these settings before running
# =============================================================================

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", "Simulation", "output"))

# ── CASE SELECTION ─────────────────────────────────────────────────────────────
# SELECT: Filter which cases to process (same format as plot_results.py)
#   Set to None to include all cases.
#
# Available filters:
#   "caverns"        - Cavern shapes to include (e.g. "regular1200", "tilted1200")
#   "pressure"       - Pressure scenario: "industry", "transport", "power_generation", "csv", or None
#   "scenario"       - Material scenario(s): "A_SIC", "A_MD", "B_SIC", "B_MD" or None
#   "n_cycles"       - Number of cycles (int) or None
#   "operation_days" - Operation duration (int) or None
#   "case_contains"  - Substring match in case name or None

SELECT = {
    "caverns": ["spike_none", "spike_lower_il4x", "spike_upper_il4x"],
    "pressure": None,
    "scenario": None,
    "n_cycles": None,
    "operation_days": 1095,
    "case_contains": None,
}

# Phase: which phase to read from ("operation" or leaching phase name)
PHASE = "operation"

# Timesteps to process: list of indices, or None for all timesteps
#   Examples:  [-1]        → last timestep only
#              [0, -1]     → first and last
#              None        → all (can be slow)
#   Special keywords (can be mixed with indices):
#              "p_min"     → timestep at minimum cavern pressure
#              "p_max"     → timestep at maximum cavern pressure
#              "p_min_top_<N>" → N deepest pressure dips, e.g. "p_min_top_3"
TIMESTEPS = [0, "p_min", "p_max", "p_min_top_3", -1]

# Cross-section settings
CROSS_SECTION_Y = 225.0       # y-coordinate for vertical slice (m)
CROSS_SECTION_THICKNESS = 5.0  # half-thickness of slice (m)

# Mohr-Coulomb (anhydrite interlayer) parameters used for f_MC visualisation.
# Failure when f_MC > 0, with f_MC = q - alpha_F * p - k_F.
MC_C_MPA = 4.0
MC_PHI_DEG = 35.0

# Minimum separation (in days) between dips reported by p_min_top_N — prevents
# multiple "dips" inside the same cycle from being reported as separate events.
P_MIN_DIP_SEP_DAYS = 5.0

# Interlayer physical-group cell tags in run_interlayer.py meshes
INTERLAYER_TAGS = (32, 34)

# =============================================================================

MPA = 1e6


# =============================================================================
# PRESSURE-BASED TIMESTEP SELECTION
# =============================================================================

def find_lowest_pressure_dips(t_s, p_mpa, n=3, min_sep_days=P_MIN_DIP_SEP_DAYS,
                              t_min_s=0.0):
    """
    Return the times (s) of the n deepest pressure dips, with a minimum
    separation in days so multiple dips in the same cycle don't collapse.

    Uses a greedy approach: pick the global minimum, then mask out a window
    of ±min_sep_days around it, repeat until n dips are found or no
    candidates remain.
    """
    p = np.asarray(p_mpa, dtype=float).copy()
    t = np.asarray(t_s, dtype=float)
    p[t < t_min_s] = np.inf
    sep_s = float(min_sep_days) * 86400.0
    picks = []
    for _ in range(int(n)):
        if not np.isfinite(p).any():
            break
        idx = int(np.argmin(p))
        picks.append(idx)
        # Mask a window around this dip
        lo, hi = t[idx] - sep_s, t[idx] + sep_s
        p[(t >= lo) & (t <= hi)] = np.inf
    return picks


def resolve_timesteps(timestep_sel, n_steps, xdmf_times, case_folder):
    """
    Resolve timestep selection to a list of integer indices.

    Handles integer indices (including negative) and special keywords:
      "p_min"          → timestep closest to minimum cavern pressure
      "p_max"          → timestep closest to maximum cavern pressure
      "p_min_top_<N>"  → the N deepest pressure dips (separated in time)
    """
    if timestep_sel is None:
        return list(range(n_steps))

    # Read pressure schedule if any keyword is used
    needs_pressure = any(isinstance(s, str) for s in timestep_sel)
    p_min_step = None
    p_max_step = None
    top_dip_steps = {}  # N -> list of XDMF step indices
    t_p_s = p_mpa = None

    if needs_pressure:
        try:
            from case_index import read_pressure_schedule
            t_p_s, p_mpa = read_pressure_schedule(case_folder)

            # For p_min: find the lowest dip DURING cycling (skip the initial
            # ramp-up period by requiring t > first XDMF step).  This avoids
            # p_min collapsing to step 0 when the schedule starts at low pressure.
            t_first_xdmf = xdmf_times[1] if len(xdmf_times) > 1 else 0.0
            active_mask = t_p_s > t_first_xdmf
            if active_mask.any():
                idx_pmin = np.where(active_mask)[0][np.argmin(p_mpa[active_mask])]
            else:
                idx_pmin = np.argmin(p_mpa)
            t_pmin = t_p_s[idx_pmin]

            # For p_max: global maximum is fine
            t_pmax = t_p_s[np.argmax(p_mpa)]

            # Find closest XDMF timestep to each
            p_min_step = int(np.argmin(np.abs(xdmf_times - t_pmin)))
            p_max_step = int(np.argmin(np.abs(xdmf_times - t_pmax)))

            print(f"[INFO] Pressure min: {p_mpa[idx_pmin]:.2f} MPa at t={t_pmin/86400:.1f} days → step {p_min_step}")
            print(f"[INFO] Pressure max: {p_mpa.max():.2f} MPa at t={t_pmax/86400:.1f} days → step {p_max_step}")
        except Exception as e:
            print(f"[WARN] Could not read pressure schedule: {e}")
            print(f"[WARN] Skipping p_min/p_max keywords")

    steps = []
    top_pat = re.compile(r"^p_min_top_(\d+)$")
    for s in timestep_sel:
        if s == "p_min":
            if p_min_step is not None:
                steps.append(p_min_step)
        elif s == "p_max":
            if p_max_step is not None:
                steps.append(p_max_step)
        elif isinstance(s, str) and top_pat.match(s):
            n_top = int(top_pat.match(s).group(1))
            if t_p_s is None or p_mpa is None:
                continue
            t_first_xdmf = xdmf_times[1] if len(xdmf_times) > 1 else 0.0
            dip_idx = find_lowest_pressure_dips(t_p_s, p_mpa, n=n_top,
                                                t_min_s=t_first_xdmf)
            for di in dip_idx:
                step = int(np.argmin(np.abs(xdmf_times - t_p_s[di])))
                steps.append(step)
                print(f"[INFO] Dip: {p_mpa[di]:.2f} MPa at t={t_p_s[di]/86400:.1f} days → step {step}")
        else:
            steps.append(int(s) % n_steps)

    # Remove duplicates while preserving order
    seen = set()
    unique = []
    for s in steps:
        if s not in seen:
            seen.add(s)
            unique.append(s)

    return unique


# =============================================================================
# INTERLAYER MASK LOADER
# =============================================================================

def load_interlayer_mask(case_folder, n_cells_xdmf, xdmf_centroids):
    """
    Return a boolean mask (length n_cells_xdmf, XDMF cell order) marking
    interlayer cells, or None if no interlayer cells / msh missing.
    """
    msh_path = os.path.join(case_folder, "operation", "mesh", "geom.msh")
    if not os.path.isfile(msh_path):
        return None
    try:
        import dolfinx as dfx
        from mpi4py import MPI
        from scipy.spatial import cKDTree
    except Exception as e:
        print(f"[WARN] Cannot load interlayer mask (missing dolfinx/mpi4py): {e}")
        return None

    try:
        mesh, cell_tags, _ = dfx.io.gmshio.read_from_msh(msh_path,
                                                         MPI.COMM_WORLD, 0)
    except Exception as e:
        print(f"[WARN] Could not read {msh_path}: {e}")
        return None

    il_dfx_cells = np.where(np.isin(cell_tags.values,
                                    list(INTERLAYER_TAGS)))[0]
    if len(il_dfx_cells) == 0:
        return None

    tdim = mesh.topology.dim
    c2v = mesh.topology.connectivity(tdim, 0)
    X = mesh.geometry.x
    n_cells_dfx = mesh.topology.index_map(tdim).size_local
    dfx_centroids = np.zeros((n_cells_dfx, 3))
    for i in range(n_cells_dfx):
        verts = X[c2v.links(i)]
        dfx_centroids[i] = verts.mean(axis=0)

    tree = cKDTree(xdmf_centroids)
    _, xdmf_idx = tree.query(dfx_centroids[il_dfx_cells])
    il_mask = np.zeros(n_cells_xdmf, dtype=bool)
    il_mask[xdmf_idx] = True
    print(f"[INFO] Interlayer cells: {il_mask.sum()} of {n_cells_xdmf} "
          f"(tags {list(INTERLAYER_TAGS)})")
    return il_mask


def compute_mc_fields(sig_3x3, il_mask, c_MPa=MC_C_MPA, phi_deg=MC_PHI_DEG):
    """
    Compute Mohr-Coulomb yield function f_MC and binary `failed` flag.
    Returns (f_MC, failed) — both length n_cells. Salt cells get NaN/0.
    """
    n_cells = sig_3x3.shape[0]
    f_MC = np.full(n_cells, np.nan)
    failed = np.zeros(n_cells, dtype=np.float64)
    if il_mask is None or not il_mask.any():
        return f_MC, failed

    sig_il = sig_3x3[il_mask]
    I1 = sig_il[:, 0, 0] + sig_il[:, 1, 1] + sig_il[:, 2, 2]
    I2 = (sig_il[:, 0, 0] * sig_il[:, 1, 1]
          + sig_il[:, 1, 1] * sig_il[:, 2, 2]
          + sig_il[:, 0, 0] * sig_il[:, 2, 2]
          - sig_il[:, 0, 1] ** 2
          - sig_il[:, 0, 2] ** 2
          - sig_il[:, 1, 2] ** 2)
    J2 = np.maximum((1.0 / 3.0) * I1 ** 2 - I2, 0.0)
    p_MPa = -I1 / (3.0 * MPA)        # compression-positive
    q_MPa = np.sqrt(3.0 * J2) / MPA

    phi = np.radians(phi_deg)
    sin_phi, cos_phi = np.sin(phi), np.cos(phi)
    alpha_F = 6.0 * sin_phi / (3.0 - sin_phi)
    k_F = 6.0 * c_MPa * cos_phi / (3.0 - sin_phi)
    f_il = q_MPa - alpha_F * p_MPa - k_F

    f_MC[il_mask] = f_il
    failed[il_mask] = (f_il > 0).astype(np.float64)
    return f_MC, failed


# =============================================================================
# BOUNDARY FACE DETECTION
# =============================================================================

def find_boundary_faces_with_tets(cells_tetra):
    """Find boundary faces and which tet they belong to."""
    face_to_tets = defaultdict(list)
    for tet_idx, tet in enumerate(cells_tetra):
        faces = [
            tuple(sorted([tet[0], tet[1], tet[2]])),
            tuple(sorted([tet[0], tet[1], tet[3]])),
            tuple(sorted([tet[0], tet[2], tet[3]])),
            tuple(sorted([tet[1], tet[2], tet[3]])),
        ]
        for f in faces:
            face_to_tets[f].append(tet_idx)
    return [(f, tets[0]) for f, tets in face_to_tets.items() if len(tets) == 1]


def extract_cavern_faces(points, cells_tetra, tol=1.0):
    """Extract cavern wall triangles (boundary faces not on the domain box)."""
    boundary = find_boundary_faces_with_tets(cells_tetra)
    x_min, x_max = points[:, 0].min(), points[:, 0].max()
    y_min, y_max = points[:, 1].min(), points[:, 1].max()
    z_min, z_max = points[:, 2].min(), points[:, 2].max()

    cavern_tris, parent_tets = [], []
    for face_nodes, tet_idx in boundary:
        verts = points[list(face_nodes)]
        on_box = False
        for dim, (lo, hi) in enumerate([(x_min, x_max), (y_min, y_max), (z_min, z_max)]):
            if np.all(np.abs(verts[:, dim] - lo) < tol) or np.all(np.abs(verts[:, dim] - hi) < tol):
                on_box = True
                break
        if not on_box:
            cavern_tris.append(face_nodes)
            parent_tets.append(tet_idx)

    return cavern_tris, parent_tets


# =============================================================================
# STRESS AND FOS COMPUTATIONS
# =============================================================================

def orient_eigenvectors(evecs, centroids, axis_xy=(225.0, 225.0)):
    """
    Fix eigenvector sign ambiguity: ensure consistent orientation.

    For each eigenvector, compute the radial direction from the cavern axis
    to the cell centroid. If the eigenvector has a negative dot product with
    the radial direction, flip it. This makes all arrows point consistently
    (outward from the cavern axis for radial-like vectors).

    evecs: (n_cells, 3, 3) — columns are eigenvectors
    centroids: (n_cells, 3)
    """
    # Radial direction from cavern axis to cell (in xy-plane)
    radial = centroids[:, :2] - np.array(axis_xy)
    radial_3d = np.zeros((centroids.shape[0], 3))
    radial_3d[:, :2] = radial
    # Normalize (avoid division by zero for cells on the axis)
    r_norm = np.linalg.norm(radial_3d, axis=1, keepdims=True)
    r_norm = np.maximum(r_norm, 1e-10)
    radial_3d = radial_3d / r_norm

    for col in range(3):
        v = evecs[:, :, col]  # (n_cells, 3)
        dot = np.sum(v * radial_3d, axis=1)  # (n_cells,)
        flip = dot < 0
        evecs[flip, :, col] *= -1

    return evecs


def compute_stress_quantities(sig_3x3, centroids):
    """Compute derived stress quantities from full 3x3 tensor (Pa)."""
    sig = 0.5 * (sig_3x3 + np.swapaxes(sig_3x3, -1, -2))
    eigenvalues, eigenvectors = np.linalg.eigh(sig)
    evals_mpa = eigenvalues / MPA

    # Fix eigenvector sign ambiguity
    eigenvectors = orient_eigenvectors(eigenvectors.copy(), centroids)

    sig_3, sig_2, sig_1 = evals_mpa[:, 0], evals_mpa[:, 1], evals_mpa[:, 2]

    p_mpa = -(sig_1 + sig_2 + sig_3) / 3.0
    q_mpa = np.sqrt(0.5 * ((sig_1 - sig_2)**2 + (sig_2 - sig_3)**2 + (sig_3 - sig_1)**2))
    p_safe = np.where(np.abs(p_mpa) > 1e-6, p_mpa, 1e-6)
    stress_ratio = q_mpa / np.abs(p_safe)

    sig_1_vec = eigenvectors[:, :, 2] * np.abs(sig_1)[:, None]
    sig_3_vec = eigenvectors[:, :, 0] * np.abs(sig_3)[:, None]

    return {
        "sig_1_MPa": sig_1, "sig_2_MPa": sig_2, "sig_3_MPa": sig_3,
        "p_MPa": p_mpa, "q_MPa": q_mpa, "stress_ratio": stress_ratio,
        "sig_1_vec_MPa": sig_1_vec, "sig_3_vec_MPa": sig_3_vec,
    }


def compute_lode_angle(sig_3x3):
    """Compute Lode angle psi from stress tensor (Pa). Shape: (n_cells,)."""
    sig = 0.5 * (sig_3x3 + np.swapaxes(sig_3x3, -1, -2))
    sig_eff = -sig  # compression positive
    vals = np.linalg.eigvalsh(sig_eff)
    s1, s2, s3 = vals[:, 2], vals[:, 1], vals[:, 0]
    I1 = s1 + s2 + s3
    mean = I1 / 3.0
    s1d, s2d, s3d = s1 - mean, s2 - mean, s3 - mean
    J2 = (1.0 / 6.0) * ((s1d - s2d)**2 + (s2d - s3d)**2 + (s3d - s1d)**2)
    J3 = s1d * s2d * s3d
    J2_safe = np.maximum(J2, 1e-30)
    x = (3.0 * np.sqrt(3.0) / 2.0) * (J3 / (J2_safe**1.5))
    x = np.clip(x, -1.0, 1.0)
    theta = (1.0 / 3.0) * np.arccos(x)
    return theta - np.pi / 6.0


def q_dil_rd(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    """Compute dilatancy boundary q_dil (MPa) from p and Lode angle."""
    I1 = 3.0 * np.asarray(p_MPa, float)
    psi = np.asarray(psi, float)
    denom = np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi)
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)
    num = D1 * (np.abs(I1) / sigma_ref) ** m + T0
    return np.sqrt(3.0) * num / denom


def compute_fos(sig_3x3, q_tol_MPa=1e-3):
    """Compute Factor of Safety from stress tensor. FOS = q_dil / q."""
    sig = 0.5 * (sig_3x3 + np.swapaxes(sig_3x3, -1, -2))
    evals = np.linalg.eigvalsh(sig) / MPA
    sig_1, sig_2, sig_3 = evals[:, 2], evals[:, 1], evals[:, 0]

    p_mpa = -(sig_1 + sig_2 + sig_3) / 3.0
    q_mpa = np.sqrt(0.5 * ((sig_1 - sig_2)**2 + (sig_2 - sig_3)**2 + (sig_3 - sig_1)**2))

    psi = compute_lode_angle(sig_3x3)
    q_dil = q_dil_rd(p_mpa, psi)

    fos = np.full_like(q_mpa, 100.0)  # cap at 100 for safe elements
    mask = q_mpa >= q_tol_MPa
    fos[mask] = q_dil[mask] / q_mpa[mask]
    fos = np.clip(fos, 0.0, 100.0)
    return fos


def compute_eps_vp_equivalent(eps_vp_3x3):
    """Compute equivalent viscoplastic strain from tensor (scalar magnitude)."""
    eps = 0.5 * (eps_vp_3x3 + np.swapaxes(eps_vp_3x3, -1, -2))
    # Equivalent strain = sqrt(2/3 * eps_ij * eps_ij)
    return np.sqrt(2.0 / 3.0 * np.sum(eps * eps, axis=(-2, -1)))


# =============================================================================
# DATA READERS
# =============================================================================

def read_tensor_timestep(xdmf_path, step_idx):
    """Read a single timestep of a tensor field, return (time, n_cells, 3, 3)."""
    reader = meshio.xdmf.TimeSeriesReader(xdmf_path)
    pts, cells = reader.read_points_cells()
    n_cells = cells["tetra"].shape[0]
    n_steps = reader.num_steps
    k = step_idx % n_steps

    for i in range(n_steps):
        t, pd, cd = reader.read_data(i)
        if i == k:
            fname = list(cd["tetra"].keys())[0]
            return float(t), cd["tetra"][fname].reshape((n_cells, 3, 3))

    raise RuntimeError(f"Step {k} not found in {xdmf_path}")


def read_vector_timestep(xdmf_path, step_idx):
    """Read a single timestep of a vector field (point data), return (time, n_nodes, 3)."""
    reader = meshio.xdmf.TimeSeriesReader(xdmf_path)
    pts, cells = reader.read_points_cells()
    n_steps = reader.num_steps
    k = step_idx % n_steps

    for i in range(n_steps):
        t, pd, cd = reader.read_data(i)
        if i == k:
            fname = list(pd.keys())[0]
            return float(t), pd[fname]

    raise RuntimeError(f"Step {k} not found in {xdmf_path}")


def read_scalar_cell_timestep(xdmf_path, step_idx):
    """Read a single timestep of a scalar cell field."""
    reader = meshio.xdmf.TimeSeriesReader(xdmf_path)
    pts, cells = reader.read_points_cells()
    n_steps = reader.num_steps
    k = step_idx % n_steps

    for i in range(n_steps):
        t, pd, cd = reader.read_data(i)
        if i == k:
            fname = list(cd["tetra"].keys())[0]
            return float(t), cd["tetra"][fname].flatten()

    raise RuntimeError(f"Step {k} not found in {xdmf_path}")


# =============================================================================
# XDMF WRITERS
# =============================================================================

def write_surface_xdmf(output_path, points, tri_array, times, data_list):
    """Write cavern surface mesh with mapped scalar fields."""
    with meshio.xdmf.TimeSeriesWriter(output_path) as writer:
        writer.write_points_cells(points, {"triangle": tri_array})
        for k, t in enumerate(times):
            cell_data = {"triangle": {
                name: arr for name, arr in data_list[k].items()
                if arr.ndim == 1
            }}
            writer.write_data(t, cell_data=cell_data)
    print(f"[SAVED] {output_path}")


def write_section_xdmf(output_path, points, cells_tetra, cell_mask,
                        times, data_list):
    """Write cross-section cells with all fields."""
    with meshio.xdmf.TimeSeriesWriter(output_path) as writer:
        writer.write_points_cells(points, {"tetra": cells_tetra[cell_mask]})
        for k, t in enumerate(times):
            cell_data = {"tetra": data_list[k]}
            writer.write_data(t, cell_data=cell_data)
    print(f"[SAVED] {output_path}")


def write_convergence_xdmf(output_path, points, cells_tetra, cell_mask,
                            times, disp_list, cell_data_list):
    """Write cross-section with displacement vectors (point data) + cell scalars."""
    sub_cells = cells_tetra[cell_mask]
    with meshio.xdmf.TimeSeriesWriter(output_path) as writer:
        writer.write_points_cells(points, {"tetra": sub_cells})
        for k, t in enumerate(times):
            writer.write_data(t, point_data=disp_list[k],
                              cell_data={"tetra": cell_data_list[k]})
    print(f"[SAVED] {output_path}")


# =============================================================================
# MAIN
# =============================================================================

def process_case(case_folder, phase, timestep_sel):
    """Process a single case folder: extract Paraview XDMF files."""

    # ---- Verify paths ----
    sig_path = os.path.join(case_folder, phase, "sig", "sig.xdmf")
    u_path = os.path.join(case_folder, phase, "u", "u.xdmf")
    eps_vp_path = os.path.join(case_folder, phase, "eps_vp", "eps_vp.xdmf")

    for p, name in [(sig_path, "sig"), (u_path, "u")]:
        if not os.path.isfile(p):
            print(f"[ERROR] Missing {name}: {p}")
            sys.exit(1)

    has_eps_vp = os.path.isfile(eps_vp_path)
    if not has_eps_vp:
        print(f"[WARN] No eps_vp found — skipping viscoplastic strain")

    # ---- Read mesh and scan timestep times from sig ----
    print(f"[INFO] Reading mesh from: {sig_path}")
    reader = meshio.xdmf.TimeSeriesReader(sig_path)
    points, cells_dict = reader.read_points_cells()
    cells_tetra = cells_dict["tetra"]
    n_cells = cells_tetra.shape[0]
    n_steps = reader.num_steps

    # Quick scan: read all timestep times (needed for p_min/p_max matching)
    xdmf_times = np.zeros(n_steps)
    for k in range(n_steps):
        t, _, _ = reader.read_data(k)
        xdmf_times[k] = float(t)

    # Re-open reader (it's consumed after reading all steps)
    reader = meshio.xdmf.TimeSeriesReader(sig_path)
    reader.read_points_cells()

    # Resolve timestep selection (handles integers, negatives, "p_min", "p_max")
    steps = resolve_timesteps(timestep_sel, n_steps, xdmf_times, case_folder)

    print(f"[INFO] Mesh: {points.shape[0]} nodes, {n_cells} cells, {n_steps} timesteps")
    print(f"[INFO] Processing {len(steps)} timestep(s): {steps}")

    # ---- Extract cavern surface ----
    print("[INFO] Extracting cavern surface...")
    cavern_tris, parent_tets = extract_cavern_faces(points, cells_tetra)
    tri_array = np.array(cavern_tris, dtype=np.int64)
    parent_idx = np.array(parent_tets, dtype=np.int64)
    print(f"[INFO] Cavern surface: {len(cavern_tris)} triangles")

    # Cavern surface node set (for mapping nodal data to surface)
    cavern_node_set = sorted({n for tri in cavern_tris for n in tri})

    # ---- Distance to cavern wall ----
    print("[INFO] Computing distance to cavern wall...")
    centroids = np.mean(points[cells_tetra], axis=1)
    from scipy.spatial import KDTree
    cavern_pts = points[cavern_node_set]
    dist_to_wall = KDTree(cavern_pts).query(centroids)[0]

    # ---- Interlayer mask (None for non-interlayer cases) ----
    il_mask = load_interlayer_mask(case_folder, n_cells, centroids)
    has_interlayer = il_mask is not None and il_mask.any()

    # ---- Cross-section mask ----
    section_mask = np.abs(centroids[:, 1] - CROSS_SECTION_Y) < CROSS_SECTION_THICKNESS
    n_section = np.sum(section_mask)
    print(f"[INFO] Cross-section (y={CROSS_SECTION_Y}±{CROSS_SECTION_THICKNESS}m): {n_section} cells")

    # ---- Z-coordinate for surface triangles (for top/mid/bottom context) ----
    tri_z = np.mean(points[tri_array, 2], axis=1)

    # ---- Read reference displacement (t=0) for incremental convergence ----
    _, u_ref = read_vector_timestep(u_path, 0)
    print(f"[INFO] Reference displacement (step 0): |u_ref| max = {np.linalg.norm(u_ref, axis=1).max()*1000:.1f} mm")

    # ---- Process timesteps ----
    times = []
    surface_data_list = []
    section_data_list = []
    convergence_point_list = []
    convergence_cell_list = []
    volume_data_list = []

    # Per-cell max-over-time accumulators (for headline damage-map XDMF written
    # after the loop). Tracks the worst state each cell ever reached.
    n_surf = len(parent_idx)
    max_acc_surf = {
        "fos_min":   np.full(n_surf, np.inf),
        "u_mag_max": np.zeros(n_surf),
    }
    max_acc_vol = {
        "fos_min":   np.full(n_cells, np.inf),
        "u_mag_max": np.zeros(n_cells),
    }
    if has_eps_vp:
        max_acc_surf["eps_vp_eq_max"] = np.zeros(n_surf)
        max_acc_vol["eps_vp_eq_max"]  = np.zeros(n_cells)
    if has_interlayer:
        max_acc_vol["f_MC_max"]    = np.full(n_cells, -np.inf)
        max_acc_vol["ever_failed"] = np.zeros(n_cells, dtype=np.float64)

    for step_k in steps:
        # -- Read sig tensor --
        t_sig, sig_3x3 = read_tensor_timestep(sig_path, step_k)
        t_days = t_sig / 86400.0
        times.append(t_days)  # store in days so Paraview shows days in toolbar
        print(f"\n  Step {step_k} | t = {t_days:.1f} days")

        # -- Stress quantities --
        quantities = compute_stress_quantities(sig_3x3, centroids)
        print(f"    p: [{quantities['p_MPa'].min():.2f}, {quantities['p_MPa'].max():.2f}] MPa")
        print(f"    q: [{quantities['q_MPa'].min():.2f}, {quantities['q_MPa'].max():.2f}] MPa")

        # -- FOS --
        fos = compute_fos(sig_3x3)
        n_dilating = np.sum(fos < 1.0)
        print(f"    FOS: [{fos.min():.3f}, {fos[fos < 100].max():.3f}] "
              f"| {n_dilating} cells with FOS < 1")

        # -- Displacement --
        _, u_data = read_vector_timestep(u_path, step_k)
        u_mag = np.linalg.norm(u_data, axis=1)  # (n_nodes,)
        u_mag_cell = np.mean(u_mag[cells_tetra], axis=1)
        u_incr_check = u_data - u_ref
        print(f"    |u|: [{u_mag.min()*1000:.2f}, {u_mag.max()*1000:.2f}] mm "
              f"(incr: [{np.linalg.norm(u_incr_check, axis=1).max()*1000:.2f}] mm)")

        # -- Viscoplastic strain --
        if has_eps_vp:
            _, eps_vp_3x3 = read_tensor_timestep(eps_vp_path, step_k)
            eps_vp_eq = compute_eps_vp_equivalent(eps_vp_3x3)
            print(f"    eps_vp_eq: [{eps_vp_eq.min():.2e}, {eps_vp_eq.max():.2e}]")

        # ---- Build full volume data (all cells) ----
        vol = {}
        for name, arr in quantities.items():
            vol[name] = arr  # scalars and vectors
        vol["FOS"] = fos
        vol["dist_to_wall"] = dist_to_wall
        vol["u_mag_mm"] = u_mag_cell * 1000.0
        if has_eps_vp:
            vol["eps_vp_eq"] = eps_vp_eq

        # Mohr-Coulomb fields (interlayer only; salt = NaN/0). Lets you
        # Threshold on `is_interlayer == 1` in ParaView and color by
        # `f_MC` or `failed`, while the surrounding salt mesh remains
        # available as a translucent gray context.
        if has_interlayer:
            f_MC, failed = compute_mc_fields(sig_3x3, il_mask)
            n_failed = int(failed.sum())
            n_il = int(il_mask.sum())
            print(f"    f_MC (interlayer): max={np.nanmax(f_MC):.2f} MPa | "
                  f"{n_failed}/{n_il} cells failed (f>0)")
            vol["is_interlayer"] = il_mask.astype(np.float64)
            vol["f_MC"] = f_MC
            vol["failed"] = failed

        volume_data_list.append(vol)

        # ---- Build surface data ----
        surf = {}
        for name, arr in quantities.items():
            if arr.ndim == 1:
                surf[name] = arr[parent_idx]
        surf["FOS"] = fos[parent_idx]
        surf["u_mag_mm"] = u_mag_cell[parent_idx] * 1000.0
        surf["z_coord"] = tri_z
        if has_eps_vp:
            surf["eps_vp_eq"] = eps_vp_eq[parent_idx]

        surface_data_list.append(surf)

        # ---- Build cross-section data ----
        sect = {}
        for name, arr in quantities.items():
            sect[name] = arr[section_mask] if arr.ndim == 1 else arr[section_mask]
        sect["FOS"] = fos[section_mask]
        sect["dist_to_wall"] = dist_to_wall[section_mask]
        sect["u_mag_mm"] = u_mag_cell[section_mask] * 1000.0
        if has_eps_vp:
            sect["eps_vp_eq"] = eps_vp_eq[section_mask]

        section_data_list.append(sect)

        # ---- Update per-cell max-over-time accumulators ----
        np.minimum(max_acc_vol["fos_min"],  fos,                       out=max_acc_vol["fos_min"])
        np.maximum(max_acc_vol["u_mag_max"], u_mag_cell * 1000.0,      out=max_acc_vol["u_mag_max"])
        np.minimum(max_acc_surf["fos_min"], fos[parent_idx],           out=max_acc_surf["fos_min"])
        np.maximum(max_acc_surf["u_mag_max"], u_mag_cell[parent_idx] * 1000.0,
                   out=max_acc_surf["u_mag_max"])
        if has_eps_vp:
            np.maximum(max_acc_vol["eps_vp_eq_max"],  eps_vp_eq,             out=max_acc_vol["eps_vp_eq_max"])
            np.maximum(max_acc_surf["eps_vp_eq_max"], eps_vp_eq[parent_idx], out=max_acc_surf["eps_vp_eq_max"])
        if has_interlayer:
            np.maximum(max_acc_vol["f_MC_max"], np.nan_to_num(f_MC, nan=-np.inf),
                       out=max_acc_vol["f_MC_max"])
            max_acc_vol["ever_failed"] = np.maximum(max_acc_vol["ever_failed"], failed.astype(np.float64))

        # ---- Build convergence data (incremental from t=0) ----
        u_incr = u_data - u_ref  # subtract reference to get convergence signal
        u_incr_mag = np.linalg.norm(u_incr, axis=1)
        conv_pt = {
            "displacement_mm": u_incr * 1000.0,
            "u_magnitude_mm": u_incr_mag * 1000.0,
        }
        conv_cell = {
            "dist_to_wall": dist_to_wall[section_mask],
        }
        convergence_point_list.append(conv_pt)
        convergence_cell_list.append(conv_cell)

    # ---- Write outputs ----
    print("\n" + "=" * 60)

    # 1. Full volume (for threshold + clip cavern shape visualization)
    volume_out = os.path.join(case_folder, "paraview_volume.xdmf")
    write_section_xdmf(volume_out, points, cells_tetra,
                       np.ones(n_cells, dtype=bool), times, volume_data_list)

    # 2. Cavern surface (clean surface colored by scalars)
    surface_out = os.path.join(case_folder, "cavern_surface.xdmf")
    write_surface_xdmf(surface_out, points, tri_array, times, surface_data_list)

    # 2b. Cavern surface — max-over-time aggregation (single timestep).
    # Direct headline "damage map": min FOS, max u, max eps_vp per surface cell.
    # Replace inf in fos_min with a finite sentinel so ParaView colormaps cleanly.
    surf_max_data = {k: v.copy() for k, v in max_acc_surf.items()}
    surf_max_data["fos_min"][~np.isfinite(surf_max_data["fos_min"])] = 1e3
    surface_max_out = os.path.join(case_folder, "cavern_surface_max.xdmf")
    write_surface_xdmf(surface_max_out, points, tri_array,
                       [times[-1] if times else 0.0], [surf_max_data])

    # 2c. Full-volume max-over-time (for cross-sections / threshold filters).
    vol_max_data = {k: v.copy() for k, v in max_acc_vol.items()}
    vol_max_data["fos_min"][~np.isfinite(vol_max_data["fos_min"])] = 1e3
    if has_interlayer:
        vol_max_data["f_MC_max"][~np.isfinite(vol_max_data["f_MC_max"])] = 0.0
    volume_max_out = os.path.join(case_folder, "paraview_volume_max.xdmf")
    write_section_xdmf(volume_max_out, points, cells_tetra,
                       np.ones(n_cells, dtype=bool),
                       [times[-1] if times else 0.0], [vol_max_data])

    # 3. Cross-section (vertical slice with stress arrows)
    section_out = os.path.join(case_folder, "cross_section.xdmf")
    write_section_xdmf(section_out, points, cells_tetra, section_mask,
                       times, section_data_list)

    # 4. Convergence (displacement vectors on cross-section)
    convergence_out = os.path.join(case_folder, "convergence.xdmf")
    write_convergence_xdmf(convergence_out, points, cells_tetra, section_mask,
                           times, convergence_point_list, convergence_cell_list)

    # ---- Instructions ----
    print(f"\n{'=' * 70}")
    print(f"  PARAVIEW VISUALIZATION GUIDE")
    print(f"{'=' * 70}")
    print(f"  Output folder: {case_folder}")
    print()
    print("─" * 70)
    print("  A. CAVERN SHAPE VISUALIZATION (paraview_volume.xdmf)")
    print("─" * 70)
    print("""
  This shows the 3D cavern shape by thresholding + clipping the full mesh.
  Great for showing different cavern geometries in your thesis.

  Step 1: Open paraview_volume.xdmf → click Apply
  Step 2: Threshold to show only near-cavern cells
          • Select paraview_volume.xdmf in Pipeline Browser
          • Filters menu → Alphabetical → Threshold
          • Properties panel (bottom-left):
              Scalars:  dist_to_wall
              Minimum:  0
              Maximum:  15       (adjust: smaller = thinner shell)
          • Click Apply
  Step 3: Clip to see inside the cavern
          • Select the Threshold in Pipeline Browser
          • Filters menu → Alphabetical → Clip
          • Properties panel:
              Clip Type: Plane
              Normal:    [0, 1, 0]   (clips along y-axis)
              Origin:    [225, 225, 330]
          • Click Apply
  Step 4: Color by a field
          • In the toolbar, change the color dropdown from
            "Solid Color" to: FOS, q_MPa, p_MPa, sig_1_MPa, etc.
""")
    print("─" * 70)
    print("  B. CAVERN SURFACE COLORING (cavern_surface.xdmf)")
    print("─" * 70)
    print("""
  Clean cavern wall surface — ideal for side-by-side comparison of shapes.

  Step 1: Open cavern_surface.xdmf → click Apply
  Step 2: In the toolbar color dropdown, select:
            • FOS           — Factor of Safety (< 1 = dilating zone)
            • q_MPa         — deviatoric stress
            • p_MPa         — mean stress (positive = compression)
            • stress_ratio  — q/|p| (higher = closer to dilatancy)
            • u_mag_mm      — displacement magnitude (mm)""")
    if has_eps_vp:
        print("            • eps_vp_eq     — accumulated viscoplastic creep strain")
    print("""
  FOS yellow/grey highlight (dilating zones in yellow):
          • Color by FOS
          • Click the "Edit Color Map" button (icon with colored bar)
          • In the Color Map Editor:
              - Set "Data Range": min=0, max=3
              - Delete all existing control points
              - Add point at value=0,     color=yellow (R=1, G=1, B=0)
              - Add point at value=1,     color=yellow
              - Add point at value=1.01,  color=light grey (R=0.8, G=0.8, B=0.8)
              - Add point at value=3,     color=light grey
          • Result: FOS < 1 = yellow, FOS > 1 = grey
""")
    print("─" * 70)
    print("  C. STRESS DIRECTION ARROWS (cross_section.xdmf)")
    print("─" * 70)
    print("""
  Vertical slice through the cavern center with principal stress arrows.
  Shows WHY certain locations have high deviatoric stress.

  Step 1: Open cross_section.xdmf → click Apply
  Step 2: Threshold near cavern wall
          • Filters → Threshold
          • Scalars: dist_to_wall,  Maximum: 10
          • Click Apply
  Step 3: Cell Data to Point Data
          • Select the Threshold in Pipeline Browser
          • Filters → Alphabetical → Cell Data to Point Data
          • Click Apply
  Step 4: Glyph (arrows)
          • Select CellDatatoPointData in Pipeline Browser
          • Filters → Alphabetical → Glyph
          • Properties:
              Glyph Type:        Arrow
              Orientation Array: sig_3_vec_MPa  (most compressive)
                            or:  sig_1_vec_MPa  (least compressive)
              Scale Array:       same vector field
              Scale Mode:        vector
              Scale Factor:      start with 0.5, adjust until arrows
                                 are visible but not overlapping
              Glyph Mode:        All Points
          • Click Apply
          • Color the arrows: toolbar dropdown → sig_3_MPa or sig_1_MPa

  Tip: sig_3 = most compressive principal stress (tangential to wall)
       sig_1 = least compressive (often radial to wall)
       Where sig_1 and sig_3 differ a lot → high q → dilatancy risk
""")
    print("─" * 70)
    print("  D. CONVERGENCE / DISPLACEMENT (convergence.xdmf)")
    print("─" * 70)
    print("""
  Shows how the cavern wall moves inward (convergence) on the cross-section.
  Compare different cavern shapes to see different deformation patterns.

  Step 1: Open convergence.xdmf → click Apply
  Step 2: Glyph (arrows showing wall movement)
          • Filters → Glyph
          • Properties:
              Orientation Array: displacement_mm
              Scale Array:       displacement_mm
              Scale Mode:        vector
              Scale Factor:      start with 1.0, adjust to taste
              Glyph Mode:        All Points  (or Every Nth Point if too dense)
          • Click Apply
          • Color by: u_magnitude_mm (displacement size)

  Tip: Displacement is INCREMENTAL (relative to start of operation).
       Near the cavern wall, arrows point inward = convergence.
       At t=0 arrows are zero; they grow as creep progresses.
       Larger arrows = more displacement = more creep.
""")
    print("─" * 70)
    print("  E. MOHR-COULOMB FAILURE IN INTERLAYERS (paraview_volume.xdmf)")
    print("─" * 70)
    print(f"""
  Cell fields `is_interlayer`, `f_MC`, `failed` are written for runs that
  contain interlayer physical groups (tags {list(INTERLAYER_TAGS)}) — anhydrite
  c={MC_C_MPA} MPa, phi={MC_PHI_DEG} deg.  Salt cells get f_MC = NaN.

  Use TIMESTEPS = ["p_min_top_3"] to capture the three deepest pressure dips.

  Step 1: Open paraview_volume.xdmf → Apply
  Step 2: Context (translucent salt block + cavern void)
          • Select the dataset → Filters → Extract Surface
          • Representation: Surface; Color: solid gray; Opacity: 0.15
  Step 3: Failure overlay (interlayer-only, colored by f_MC or failed)
          • Re-select paraview_volume.xdmf → Filters → Threshold
              Scalars:  is_interlayer
              Minimum:  1
              Maximum:  1
          • Apply, then color by `f_MC` (sequential, centered at 0)
            or `failed` (binary 0/1).
          • Optional: another Threshold on `failed >= 1` to show ONLY
            failed cells.
  Step 4: Use the time slider to flip between the three pressure dips.
""")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Paraview visualizations from simulation output."
    )
    parser.add_argument("case_folder", nargs="?", default=None,
                        help="Single case folder (overrides SELECT filter)")
    parser.add_argument("--phase", default=None)
    parser.add_argument("--timesteps", nargs="*", type=int, default=None)
    args = parser.parse_args()

    phase = args.phase or PHASE
    timestep_sel = args.timesteps if args.timesteps is not None else TIMESTEPS

    # If a specific case folder is given on the command line, use only that
    if args.case_folder:
        case_raw = args.case_folder
        if not os.path.isabs(case_raw):
            case_folder = os.path.normpath(os.path.join(SCRIPT_DIR, case_raw))
        else:
            case_folder = os.path.abspath(case_raw)
        print(f"[INFO] Processing single case: {case_folder}")
        process_case(case_folder, phase, timestep_sel)
        return

    # Otherwise, use SELECT filter to discover cases (same as plot_results.py)
    all_cases = detect_layout_and_collect_cases(ROOT)
    cases = filter_cases(all_cases, SELECT)

    if not cases:
        print(f"[ERROR] No cases matched SELECT filter in {ROOT}")
        print(f"[INFO]  SELECT = {SELECT}")
        print(f"[INFO]  Total cases found: {len(all_cases)}")
        for m in all_cases[:10]:
            print(f"         {m.get('case_name')}")
        return

    print(f"[INFO] Found {len(cases)} case(s) matching SELECT filter:")
    for m in cases:
        print(f"  - {m['case_name']}")
    print()

    for i, meta in enumerate(cases):
        case_folder = meta["case_path"]
        print(f"\n{'='*70}")
        print(f"  CASE {i+1}/{len(cases)}: {meta['case_name']}")
        print(f"{'='*70}")
        try:
            process_case(case_folder, phase, timestep_sel)
        except Exception as e:
            print(f"[ERROR] Failed: {e}")
            continue


if __name__ == "__main__":
    main()
