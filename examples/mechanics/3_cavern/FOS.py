import os
import numpy as np
from mpi4py import MPI
import dolfinx as do
import safeincave.PostProcessingTools as post

# Units
MPa = 1e6



# Units
MPa = 1e6

def q_dil_rd(p_MPa, psi,
             D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    """
    Compute q_dil(p, psi) in MPa for the Rider–Dieterich dilatancy boundary.

    p_MPa : mean stress in MPa, compression-positive (zoals in je plots)
    psi   : Lode angle in radians (zelfde definitie als in je slide)
    """
    p = np.asarray(p_MPa, dtype=float)
    psi = np.asarray(psi, dtype=float)

    # I1 in MPa (let op: jouw conventie p>0 in compressie)
    I1 = 3.0 * p  # MPa

    # sign-truc zoals in de literatuur
    sgn = np.sign(I1)
    sgn[sgn == 0.0] = 1.0

    denom = (np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi))
    sqrtJ2 = D1 * ((I1 / (sgn * sigma_ref)) ** m) / denom + T0  # MPa

    q = np.sqrt(3.0) * sqrtJ2  # MPa
    return q



def compute_lode_angle_from_sigma(sig_Pa):
    """
    sig_Pa : (ncells, 3, 3) stress tensor in Pa
    returns psi : (ncells,) Lode angle in radians
    """
    # Eigenwaarden (principal total stresses)
    # omdat sigma symmetrisch is: eigvalsh
    vals = np.linalg.eigvalsh(sig_Pa)   # shape (ncells, 3)
    s1, s2, s3 = vals[:, 2], vals[:, 1], vals[:, 0]  # geordend 

    # deviatorische principalen
    I1 = s1 + s2 + s3
    mean = I1 / 3.0
    s1d = s1 - mean
    s2d = s2 - mean
    s3d = s3 - mean

    # invarianten J2, J3 van deviatorische spanningen
    J2 = (1.0 / 6.0) * ((s1d - s2d) ** 2 + (s2d - s3d) ** 2 + (s3d - s1d) ** 2)
    J3 = s1d * s2d * s3d

    # vermijd deling door nul
    tiny = 1e-30
    J2_safe = np.maximum(J2, tiny)

    arg = - (np.sqrt(27.0) / 2.0) * J3 / (J2_safe ** 1.5)
    arg = np.clip(arg, -1.0, 1.0)

    psi = (1.0 / 3.0) * np.arccos(arg) # Later checken of je hier niet Sinus moet nemen!!!
    return psi




# ---------------------------------------------------------------------
# Load p and q like your StressPath class does
# ---------------------------------------------------------------------
def load_p_q(operation_folder):
    """
    Returns:
        time_list : (nt,) array of times [s]
        p_MPa     : (nt, ncells) mean stress, compression-positive [MPa]
        q_MPa     : (nt, ncells) von Mises stress [MPa]
    """
    p_path = os.path.join(operation_folder, "p_elems", "p_elems.xdmf")
    q_path = os.path.join(operation_folder, "q_elems", "q_elems.xdmf")

    pts_p, time_list, p_vals = post.read_cell_scalar(p_path)
    pts_q, time_list2, q_vals = post.read_cell_scalar(q_path)

    if not np.allclose(time_list, time_list2):
        raise RuntimeError("Time lists of p_elems and q_elems do not match.")

    # Same convention as your StressPath: compression -> positive p
    p_MPa = -(p_vals / 3) / MPa   # (nt, ncells)
    q_MPa =  q_vals / MPa   # (nt, ncells)

    return time_list, p_MPa, q_MPa, p_path





def load_p_q_sig(operation_folder):
    p_path  = os.path.join(operation_folder, "p_elems", "p_elems.xdmf")
    q_path  = os.path.join(operation_folder, "q_elems", "q_elems.xdmf")
    sig_path = os.path.join(operation_folder, "sig", "sig.xdmf")

    pts_p,  time_list,  p_vals   = post.read_cell_scalar(p_path)
    pts_q,  time_list2, q_vals   = post.read_cell_scalar(q_path)
    pts_s,  time_list3, sig_vals = post.read_cell_tensor(sig_path)

    if not (np.allclose(time_list, time_list2) and np.allclose(time_list, time_list3)):
        raise RuntimeError("Time lists of p_elems, q_elems and sig do not match.")

    # jouw conventie: compressie-positief p
    p_MPa = -(p_vals / 3) / MPa       # (nt, ncells)
    q_MPa =  q_vals / MPa       # (nt, ncells)

    # sig_vals in Pa, shape (nt, ncells, 3, 3)
    return time_list, p_MPa, q_MPa, sig_vals, p_path


def compute_FOS_with_lode(p_MPa, q_MPa, sig_Pa_series):
    """
    p_MPa         : (nt, ncells)
    q_MPa         : (nt, ncells)
    sig_Pa_series : (nt, ncells, 3, 3)

    returns FOS : (nt, ncells)
    """
    nt, nc = p_MPa.shape
    FOS = np.zeros_like(p_MPa)

    for it in range(nt):
        # Lode angle per element op deze tijdstap
        sig_t = sig_Pa_series[it, :, :, :]    # (ncells, 3, 3)
        psi_t = compute_lode_angle_from_sigma(sig_t)  # (ncells,)

        p_t = p_MPa[it, :]   # (ncells,)
        q_t = q_MPa[it, :]   # (ncells,)

        q_dil_t = q_dil_rd(p_t, psi_t)  # MPa

        q_safe = np.where(q_t <= 0.0, 1e-12, q_t)
        FOS[it, :] = q_dil_t / q_safe

    return FOS



# ---------------------------------------------------------------------
# Write FOS as DG0 cell field to XDMF using mesh from p_elems.xdmf
# ---------------------------------------------------------------------
def write_FOS(operation_folder, time_list, FOS, mesh_source_xdmf):
    print("Writing FOS field...\n")

    # Mesh uit p_elems.xdmf halen (die heeft p_elems.h5 ernaast)
    with do.io.XDMFFile(MPI.COMM_WORLD, mesh_source_xdmf, "r") as xdmf:
        mesh = xdmf.read_mesh()
        mesh.name = "mesh"

    # ⬇️ HIER zit de bug – gebruik functionspace i.p.v. FunctionSpace
    V = do.fem.functionspace(mesh, ("DG", 0))
    fos_fun = do.fem.Function(V)
    fos_fun.name = "FOS"

    fos_folder = os.path.join(operation_folder, "FOS")
    os.makedirs(fos_folder, exist_ok=True)
    fos_path = os.path.join(fos_folder, "FOS.xdmf")

    with do.io.XDMFFile(MPI.COMM_WORLD, fos_path, "w") as xdmf_out:
        xdmf_out.write_mesh(mesh)

        nt = len(time_list)
        for it, t in enumerate(time_list):
            # FOS heeft vorm (nt, ncells); DG0 heeft 1 dof per cel
            fos_fun.x.array[:] = FOS[it, :]
            xdmf_out.write_function(fos_fun, float(t))

    if MPI.COMM_WORLD.rank == 0:
        print(f"FOS written to: {fos_path}")




def main():
    base_folder = os.path.join("output", "case_sinus_70days_irregular", "operation")

    if MPI.COMM_WORLD.rank == 0:
        print("Loading stresses...")
    time_list, p_MPa, q_MPa, sig_vals, p_source = load_p_q_sig(base_folder)

    if MPI.COMM_WORLD.rank == 0:
        print("Computing FOS with Lode angle...")
    FOS = compute_FOS_with_lode(p_MPa, q_MPa, sig_vals)

    if MPI.COMM_WORLD.rank == 0:
        print("Writing FOS field...")
    write_FOS(base_folder, time_list, FOS, p_source)


if __name__ == "__main__":
    main()
