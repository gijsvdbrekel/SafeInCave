import os
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
import dolfinx as do

import safeincave.PostProcessingTools as post

MPA  = 1e6
HOUR = 3600.0
DAY  = 24.0 * HOUR

# ---------- consistent naming + colors ----------
CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt"]

def cavern_label_from_folder(folder_name: str) -> str:
    return folder_name.split("_")[0]

def build_color_map(labels):
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    return {lab: cycle[i % len(cycle)] for i, lab in enumerate(labels)}

def pressure_scheme_from_folder(folder_name: str) -> str:
    parts = folder_name.split("_")
    return parts[1].lower() if len(parts) > 1 else ""

# ---------- RD model ----------
def q_dil_rd(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    """
    RD dilation boundary q_dil(p, psi) in MPa.

    Uses p = mean stress (compression-positive), psi = Lode angle.
    Implements the robust form with |I1| inside the power to avoid sign/power issues.
    """
    p = np.asarray(p_MPa, dtype=float)
    psi = np.asarray(psi, dtype=float)

    I1 = 3.0 * p  # MPa
    absI1 = np.abs(I1)

    denom = (np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi))
    # avoid blow-up; preserve sign if denom is tiny
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)

    sqrtJ2_dil = D1 * ((absI1 / sigma_ref) ** m) / denom + T0  # MPa
    q_dil = np.sqrt(3.0) * sqrtJ2_dil                          # MPa
    return q_dil

def compute_lode_angle_from_sigma(sig_Pa):
    """
    Returns psi in [-pi/6, +pi/6] per cell.
    sig_Pa: (ncells, 3, 3)
    """
    sig = 0.5 * (sig_Pa + np.swapaxes(sig_Pa, -1, -2))
    vals = np.linalg.eigvalsh(sig)
    s1, s2, s3 = vals[:, 2], vals[:, 1], vals[:, 0]

    I1 = s1 + s2 + s3
    mean = I1 / 3.0
    s1d = s1 - mean
    s2d = s2 - mean
    s3d = s3 - mean

    J2 = (1.0/6.0) * ((s1d - s2d)**2 + (s2d - s3d)**2 + (s3d - s1d)**2)
    J3 = s1d * s2d * s3d

    tiny = 1e-30
    J2_safe = np.maximum(J2, tiny)

    x = (3.0*np.sqrt(3.0)/2.0) * (J3 / (J2_safe**1.5))
    x = np.clip(x, -1.0, 1.0)

    theta = (1.0/3.0) * np.arccos(x)
    psi = theta - np.pi/6.0
    return psi

def compute_p_q_from_sigma(sig_Pa, compression_positive=True):
    """
    Compute p (mean stress, MPa, compression-positive) and q (von Mises, MPa)
    directly from the full stress tensor, per cell.

    sig_Pa: (ncells, 3, 3) in Pa
    """
    sig = 0.5 * (sig_Pa + np.swapaxes(sig_Pa, -1, -2))
    if compression_positive:
        sig = -sig  # align with SafeInCave's constitutive convention (see Desai: stress_vec = -stress)

    # p = tr(sig)/3
    I1 = np.trace(sig, axis1=1, axis2=2)          # Pa
    p = (I1 / 3.0) / MPA                           # MPa

    # deviatoric + J2 -> q = sqrt(3*J2)
    mean = (I1 / 3.0)[:, None, None]
    dev = sig - mean * np.eye(3)[None, :, :]
    J2 = 0.5 * np.sum(dev * dev, axis=(1, 2))      # Pa^2
    q = np.sqrt(3.0 * np.maximum(J2, 0.0)) / MPA   # MPa

    return p, q

# ---------- robust path helpers (flat vs nested output) ----------
def pick_existing(*paths):
    for p in paths:
        if os.path.isfile(p):
            return p
    return None

def load_sig(case_folder):
    """
    Supports both:
      case_folder/sig.xdmf
      case_folder/sig/sig.xdmf
    """
    sig_path = pick_existing(
        os.path.join(case_folder, "sig.xdmf"),
        os.path.join(case_folder, "sig", "sig.xdmf"),
    )
    if sig_path is None:
        raise FileNotFoundError(f"Missing sig in {case_folder}")

    pts_s, t, sig_vals = post.read_cell_tensor(sig_path)
    return np.asarray(t, float), np.asarray(sig_vals, float), sig_path

def compute_FOS(time_list, sig_vals, *, compression_positive=True, q_tol_MPa=1e-3):
    """
    Compute FoS per time step, per cell:
        FoS = q_dil(p, psi) / q

    IMPORTANT:
    - Computes p, q, psi consistently from the SAME stress tensor (sig_vals).
    - Uses compression-positive convention by default to match SafeInCave (Desai flips sign).
    - Masks cells with q < q_tol_MPa to avoid meaningless huge ratios in near-hydrostatic zones.
      Those cells get FoS = +inf so they DO NOT control your global minimum.
    """
    nt = sig_vals.shape[0]
    nc = sig_vals.shape[1]
    FOS = np.full((nt, nc), np.inf, dtype=float)

    for it in range(nt):
        sig_t = sig_vals[it]  # (nc, 3, 3) Pa

        # compute consistently from sig
        p_MPa, q_MPa = compute_p_q_from_sigma(sig_t, compression_positive=compression_positive)
        psi = compute_lode_angle_from_sigma((-sig_t) if compression_positive else sig_t)

        q_dil = q_dil_rd(p_MPa, psi)  # MPa

        mask = q_MPa >= float(q_tol_MPa)
        FOS[it, mask] = q_dil[mask] / q_MPa[mask]

    # clean-up
    FOS[~np.isfinite(FOS)] = np.inf
    FOS = np.clip(FOS, 0.0, 1e6)
    return FOS

def write_FOS_paraview(case_folder, time_list, FOS, mesh_source_xdmf, out_folder, field_name="FOS"):
    """
    Writes DG0 cell field so ParaView shows it under Cell Data with correct name.

    NOTE (strict): This assumes the cell ordering in meshio's TimeSeriesReader output matches
    dolfinx's XDMF mesh cell ordering. If you see scrambled spatial patterns, this is why.
    """
    os.makedirs(out_folder, exist_ok=True)
    out_path = os.path.join(out_folder, f"{os.path.basename(case_folder)}_{field_name}.xdmf")

    with do.io.XDMFFile(MPI.COMM_SELF, mesh_source_xdmf, "r") as xdmf:
        mesh = xdmf.read_mesh()
        mesh.name = "mesh"

    V = do.fem.functionspace(mesh, ("DG", 0))
    f = do.fem.Function(V)
    f.name = field_name

    ncells = FOS.shape[1]
    if f.x.array.size != ncells:
        raise RuntimeError(
            f"DG0 dofs ({f.x.array.size}) != ncells in FOS ({ncells}). "
            "You are almost certainly reading/writing with different meshes or different cell ordering."
        )

    with do.io.XDMFFile(MPI.COMM_SELF, out_path, "w") as xdmf_out:
        xdmf_out.write_mesh(mesh)
        for it, t in enumerate(time_list):
            f.x.array[:] = FOS[it]
            xdmf_out.write_function(f, float(t))

    print(f"[OK] Wrote ParaView FoS: {out_path}")
    return out_path

def main():
    if MPI.COMM_WORLD.rank != 0:
        return

    ROOT = r"/home/gvandenbrekel/SafeInCave/OutputNobian"
    TARGET_PRESSURE = "irregular"

    # SafeInCave constitutive convention (Desai): compression-positive uses -sigma
    COMPRESSION_POSITIVE = True

    # Avoid FoS blowing up where q ~ 0 (near hydrostatic); adjust if needed
    Q_TOL_MPA = 1e-3  # 0.001 MPa = 1 kPa

    # collect cases
    cases = []
    for nm in sorted(os.listdir(ROOT)):
        fpath = os.path.join(ROOT, nm)
        if not os.path.isdir(fpath):
            continue
        if nm.lower().startswith("pressure_"):
            continue
        if pressure_scheme_from_folder(nm) != TARGET_PRESSURE:
            continue

        try:
            load_sig(fpath)
            cases.append((nm, fpath))
        except Exception:
            continue

    if not cases:
        raise RuntimeError(f"No cases found for pressure scheme '{TARGET_PRESSURE}' with sig.")

    # consistent colors
    labels_present = []
    for nm, _ in cases:
        lab = cavern_label_from_folder(nm)
        if lab not in labels_present:
            labels_present.append(lab)

    labels_sorted = [l for l in CAVERN_ORDER if l in labels_present] + \
                    [l for l in labels_present if l not in CAVERN_ORDER]
    color_map = build_color_map(labels_sorted)

    OUT = os.path.join(ROOT, "_FOS_outputs", TARGET_PRESSURE)
    os.makedirs(OUT, exist_ok=True)

    plt.figure(figsize=(12, 6))

    for nm, folder in cases:
        lab = cavern_label_from_folder(nm)
        col = color_map.get(lab, None)

        print(f"\n=== {nm} ===")
        time_list, sig_vals, mesh_source = load_sig(folder)

        FOS = compute_FOS(
            time_list, sig_vals,
            compression_positive=COMPRESSION_POSITIVE,
            q_tol_MPa=Q_TOL_MPA
        )

        write_FOS_paraview(folder, time_list, FOS, mesh_source, OUT, field_name="FOS")

        # global minimum FoS across all cells, ignoring q<tol cells (they are inf)
        fos_min = np.min(FOS, axis=1)
        fos_min[~np.isfinite(fos_min)] = np.nan  # if everything was masked, show nan

        t_days = time_list / DAY
        plt.plot(t_days, fos_min, linewidth=2.2, label=lab, color=col)

    plt.axhline(1.0, linewidth=1.5)
    plt.xlabel("Time (days)")
    plt.ylabel("FoS = q_dil(p, ψ) / q   (min over cells; cells with q<tol ignored)")
    plt.grid(True, alpha=0.3)
    plt.title(f"Factor of Safety over time (pressure scheme = {TARGET_PRESSURE})")

    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h
    ax.legend(uniq.values(), uniq.keys(), fontsize=9, frameon=True, loc="best")

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
