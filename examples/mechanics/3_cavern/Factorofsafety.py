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
    RD dilation boundary q_dil(p,psi) in MPa.
    """
    p = np.asarray(p_MPa, dtype=float)
    psi = np.asarray(psi, dtype=float)

    I1 = 3.0 * p
    sgn = np.sign(I1)
    sgn[sgn == 0.0] = 1.0

    denom = (np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi))
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)

    sqrtJ2 = D1 * ((I1 / (sgn * sigma_ref)) ** m) / denom + T0
    q = np.sqrt(3.0) * sqrtJ2
    return q

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

# ---------- robust path helpers (flat vs nested output) ----------
def pick_existing(*paths):
    for p in paths:
        if os.path.isfile(p):
            return p
    return None

def load_p_q_sig(case_folder):
    """
    Supports both:
      case_folder/p_elems.xdmf
      case_folder/p_elems/p_elems.xdmf
    same for q_elems and sig.
    """
    p_path = pick_existing(
        os.path.join(case_folder, "p_elems.xdmf"),
        os.path.join(case_folder, "p_elems", "p_elems.xdmf"),
    )
    q_path = pick_existing(
        os.path.join(case_folder, "q_elems.xdmf"),
        os.path.join(case_folder, "q_elems", "q_elems.xdmf"),
    )
    sig_path = pick_existing(
        os.path.join(case_folder, "sig.xdmf"),
        os.path.join(case_folder, "sig", "sig.xdmf"),
    )

    if p_path is None or q_path is None or sig_path is None:
        raise FileNotFoundError(f"Missing p_elems/q_elems/sig in {case_folder}")

    pts_p, t1, p_vals = post.read_cell_scalar(p_path)
    pts_q, t2, q_vals = post.read_cell_scalar(q_path)
    pts_s, t3, sig_vals = post.read_cell_tensor(sig_path)

    if not (np.allclose(t1, t2) and np.allclose(t1, t3)):
        # clip to common length as safe fallback
        nt = min(len(t1), len(t2), len(t3))
        t1 = t1[:nt]
        p_vals = p_vals[:nt]
        q_vals = q_vals[:nt]
        sig_vals = sig_vals[:nt]

    return np.asarray(t1, float), np.asarray(p_vals, float), np.asarray(q_vals, float), np.asarray(sig_vals, float), p_path

def compute_FOS(time_list, p_vals, q_vals, sig_vals, divide_p_by_3=False):
    """
    p_vals, q_vals from PostProcessingTools are in Pa or MPa depending on your pipeline.
    In your setup: p_vals and q_vals are in Pa -> convert to MPa by /MPA.
    Convention: p positive in compression.
    """
    # convert to MPa and sign convention
    # p stored negative in compression in your plots -> use minus
    p_MPa = -p_vals / MPA
    q_MPa =  q_vals / MPA

    if divide_p_by_3:
        p_MPa = p_MPa / 3.0

    nt, nc = p_MPa.shape
    FOS = np.zeros((nt, nc), dtype=float)

    for it in range(nt):
        psi = compute_lode_angle_from_sigma(sig_vals[it])          # (nc,)
        q_dil = q_dil_rd(p_MPa[it], psi)                           # (nc,)
        q_safe = np.where(q_MPa[it] <= 0.0, 1e-12, q_MPa[it])
        FOS[it] = q_dil / q_safe

    FOS[~np.isfinite(FOS)] = 0.0
    FOS = np.clip(FOS, 0.0, 1e6)
    return FOS

def write_FOS_paraview(case_folder, time_list, FOS, mesh_source_xdmf, out_folder, field_name="FOS"):
    """
    Writes DG0 cell field so ParaView shows it under Cell Data with correct name.
    """
    os.makedirs(out_folder, exist_ok=True)
    out_path = os.path.join(out_folder, f"{os.path.basename(case_folder)}_{field_name}.xdmf")

    # read mesh SERIAL
    with do.io.XDMFFile(MPI.COMM_SELF, mesh_source_xdmf, "r") as xdmf:
        mesh = xdmf.read_mesh()
        mesh.name = "mesh"

    V = do.fem.functionspace(mesh, ("DG", 0))
    f = do.fem.Function(V)
    f.name = field_name

    # sanity: DG0 dofs should match ncells
    ncells = FOS.shape[1]
    if f.x.array.size != ncells:
        raise RuntimeError(
            f"DG0 dofs ({f.x.array.size}) != ncells in FOS ({ncells}). "
            "This usually means you're reading the wrong mesh source."
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

    # if your no-GUI run laptop has the old p bug:
    DIVIDE_P_BY_3 = True

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

        # minimal check: must have p/q/sig somewhere
        try:
            load_p_q_sig(fpath)
            cases.append((nm, fpath))
        except Exception:
            continue

    if not cases:
        raise RuntimeError(f"No cases found for pressure scheme '{TARGET_PRESSURE}' with p/q/sig.")

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

    # plot: ONLY one line per cavern = global min over all cells
    plt.figure(figsize=(12, 6))

    for nm, folder in cases:
        lab = cavern_label_from_folder(nm)
        col = color_map.get(lab, None)

        print(f"\n=== {nm} ===")
        time_list, p_vals, q_vals, sig_vals, mesh_source = load_p_q_sig(folder)
        FOS = compute_FOS(time_list, p_vals, q_vals, sig_vals, divide_p_by_3=DIVIDE_P_BY_3)

        # write for ParaView (this WILL show the array as "FOS")
        write_FOS_paraview(folder, time_list, FOS, mesh_source, OUT, field_name="FOS")

        # summary: global minimum FoS across all volume cells
        fos_min = FOS.min(axis=1)
        t_days = time_list / DAY
        plt.plot(t_days, fos_min, linewidth=2.2, label=lab, color=col)

    plt.axhline(1.0, linewidth=1.5)
    plt.xlabel("Time (days)")
    plt.ylabel("FoS = q_dil(p, ψ) / q  (plotted: min over all cells)")
    plt.grid(True, alpha=0.3)
    plt.title(f"Factor of Safety over time (pressure scheme = {TARGET_PRESSURE})")

    # legend with only cavern names (unique)
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
