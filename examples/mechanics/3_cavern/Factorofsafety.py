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

def cavern_label_from_group(group_folder: str) -> str:
    """
    group folder looks like: Asymmetric_sinus_600 or Tilt_irregular_600
    We want the cavern type label robustly.
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
    if low.startswith("irregular"):
        return "Irregular"
    return group_folder.split("_")[0]

def build_color_map(labels):
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    return {lab: cycle[i % len(cycle)] for i, lab in enumerate(labels)}

def scheme_from_case_folder(case_folder_name: str) -> str:
    low = case_folder_name.lower()
    if low.startswith("case_sinus"):
        return "sinus"
    if low.startswith("case_irregular"):
        return "irregular"
    return ""

# ---------- RD model ----------
def q_dil_rd(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    """
    RD dilation boundary q_dil(p, psi) in MPa.
    p_MPa: mean stress p (compression-positive) in MPa
    psi  : Lode angle in [-pi/6, +pi/6]
    """
    p = np.asarray(p_MPa, dtype=float)
    psi = np.asarray(psi, dtype=float)

    I1 = 3.0 * p                   # MPa
    absI1 = np.abs(I1)

    denom = (np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi))
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)

    sqrtJ2_dil = D1 * ((absI1 / sigma_ref) ** m) / denom + T0
    q_dil = np.sqrt(3.0) * sqrtJ2_dil
    return q_dil

# ---------- invariants from sigma (ONE consistent pathway) ----------
def compute_p_q_psi_from_sigma(sig_Pa, compression_positive=True):
    """
    sig_Pa: (ncells, 3, 3) in Pa.

    If compression_positive=True, use sig_eff = -sig.
    Returns:
      p in MPa, q in MPa, psi in [-pi/6, +pi/6]
    """
    sig = 0.5 * (sig_Pa + np.swapaxes(sig_Pa, -1, -2))  # enforce symmetry
    sig_eff = -sig if compression_positive else sig

    # principal stresses
    vals = np.linalg.eigvalsh(sig_eff)
    s1, s2, s3 = vals[:, 2], vals[:, 1], vals[:, 0]     # s1 >= s2 >= s3

    I1 = s1 + s2 + s3                                   # Pa
    p_MPa = (I1 / 3.0) / MPA                            # MPa

    mean = I1 / 3.0
    s1d = s1 - mean
    s2d = s2 - mean
    s3d = s3 - mean

    J2 = (1.0/6.0) * ((s1d - s2d)**2 + (s2d - s3d)**2 + (s3d - s1d)**2)   # Pa^2
    J3 = s1d * s2d * s3d                                                  # Pa^3

    q_MPa = np.sqrt(3.0 * np.maximum(J2, 0.0)) / MPA

    tiny = 1e-30
    J2_safe = np.maximum(J2, tiny)
    x = (3.0*np.sqrt(3.0)/2.0) * (J3 / (J2_safe**1.5))
    x = np.clip(x, -1.0, 1.0)
    theta = (1.0/3.0) * np.arccos(x)
    psi = theta - np.pi/6.0

    return p_MPa, q_MPa, psi

# ---------- robust path helpers ----------
def pick_existing(*paths):
    for p in paths:
        if os.path.isfile(p):
            return p
    return None

def load_sig(case_folder):
    """
    OutputNobian layout:
      case_folder/operation/sig/sig.xdmf
      case_folder/operation/sig.xdmf          (fallback)
    """
    sig_path = pick_existing(
        os.path.join(case_folder, "operation", "sig", "sig.xdmf"),
        os.path.join(case_folder, "operation", "sig.xdmf"),
    )
    if sig_path is None:
        raise FileNotFoundError(f"Missing sig in {case_folder} (checked operation/sig/sig.xdmf and operation/sig.xdmf)")

    _, t, sig_vals = post.read_cell_tensor(sig_path)
    return np.asarray(t, float), np.asarray(sig_vals, float), sig_path

def compute_FOS(time_list, sig_vals, *, compression_positive=True, q_tol_MPa=1e-3):
    """
    FoS per time step, per cell:
        FoS = q_dil(p, psi) / q

    Strict:
      - p,q,psi from sig only
      - q<tol -> FoS=+inf (ignored in min)
    """
    nt = sig_vals.shape[0]
    nc = sig_vals.shape[1]
    FOS = np.full((nt, nc), np.inf, dtype=float)

    for it in range(nt):
        p_MPa, q_MPa, psi = compute_p_q_psi_from_sigma(
            sig_vals[it],
            compression_positive=compression_positive
        )

        q_dil = q_dil_rd(p_MPa, psi)  # MPa

        mask = q_MPa >= float(q_tol_MPa)
        FOS[it, mask] = q_dil[mask] / q_MPa[mask]

    FOS[~np.isfinite(FOS)] = np.inf
    FOS = np.clip(FOS, 0.0, 1e6)
    return FOS

def write_FOS_paraview(case_folder, time_list, FOS, mesh_source_xdmf, out_folder, field_name="FOS"):
    """
    Writes DG0 cell field to XDMF for ParaView.
    Uses mesh from mesh_source_xdmf (sig.xdmf usually contains it).
    """
    os.makedirs(out_folder, exist_ok=True)

    # make output filename informative: <group>__<case>_FOS.xdmf
    group = os.path.basename(os.path.dirname(case_folder))
    case  = os.path.basename(case_folder)
    out_path = os.path.join(out_folder, f"{group}__{case}_{field_name}.xdmf")

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
            "Likely wrong mesh source or mismatched cell ordering."
        )

    with do.io.XDMFFile(MPI.COMM_SELF, out_path, "w") as xdmf_out:
        xdmf_out.write_mesh(mesh)
        for it, t in enumerate(time_list):
            f.x.array[:] = FOS[it]
            xdmf_out.write_function(f, float(t))

    print(f"[OK] Wrote ParaView FoS: {out_path}")
    return out_path

# ---------- nested case collector for OutputNobian ----------
def collect_cases_nested(ROOT: str, target_scheme: str):
    """
    ROOT/<GROUP>/<CASE>/...
    Filter by scheme using CASE folder name.
    Returns list of dicts: {label, group, case_name, case_path}
    """
    cases = []

    for group in sorted(os.listdir(ROOT)):
        group_path = os.path.join(ROOT, group)
        if not os.path.isdir(group_path):
            continue
        if group.lower().startswith("pressure_"):
            continue
        if group.lower().startswith("_fos_outputs"):
            continue

        lab = cavern_label_from_group(group)

        for sub in sorted(os.listdir(group_path)):
            if not sub.lower().startswith("case_"):
                continue

            sc = scheme_from_case_folder(sub)
            if sc != target_scheme:
                continue

            case_path = os.path.join(group_path, sub)
            if not os.path.isdir(case_path):
                continue

            # must have sig
            try:
                load_sig(case_path)
            except Exception as e:
                # keep it silent unless you want verbose:
                # print(f"[SKIP] {group}/{sub}: {e}")
                continue

            cases.append({
                "label": lab,
                "group": group,
                "case_name": sub,
                "case_path": case_path,
            })

    return cases

def main():
    if MPI.COMM_WORLD.rank != 0:
        return

    ROOT = r"/home/gvandenbrekel/SafeInCave/OutputNobian"
    TARGET_PRESSURE = "irregular"   # "sinus" or "irregular"

    COMPRESSION_POSITIVE = True
    Q_TOL_MPA = 1e-3  # 1 kPa

    cases = collect_cases_nested(ROOT, TARGET_PRESSURE)
    if not cases:
        raise RuntimeError(f"No cases found for pressure scheme '{TARGET_PRESSURE}' with sig under {ROOT}")

    # consistent colors based on labels present
    labels_present = []
    for c in cases:
        if c["label"] not in labels_present:
            labels_present.append(c["label"])

    labels_sorted = [l for l in CAVERN_ORDER if l in labels_present] + \
                    [l for l in labels_present if l not in CAVERN_ORDER]
    color_map = build_color_map(labels_sorted)

    OUT = os.path.join(ROOT, "_FOS_outputs", TARGET_PRESSURE)
    os.makedirs(OUT, exist_ok=True)

    plt.figure(figsize=(12, 6))

    for c in cases:
        lab = c["label"]
        col = color_map.get(lab, None)

        print(f"\n=== {c['group']}/{c['case_name']} ===")
        time_list, sig_vals, mesh_source = load_sig(c["case_path"])

        FOS = compute_FOS(
            time_list, sig_vals,
            compression_positive=COMPRESSION_POSITIVE,
            q_tol_MPa=Q_TOL_MPA
        )

        write_FOS_paraview(c["case_path"], time_list, FOS, mesh_source, OUT, field_name="FOS")

        fos_min = np.min(FOS, axis=1)
        fos_min[~np.isfinite(fos_min)] = np.nan  # if everything masked

        t_days = time_list / DAY
        plt.plot(t_days, fos_min, linewidth=2.2, label=lab, color=col)

    plt.axhline(1.0, linewidth=1.5)
    plt.xlabel("Time (days)")
    plt.ylabel("FoS = q_dil(p, ψ) / q   (min over cells; q<tol ignored)")
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
