import os
import numpy as np
from mpi4py import MPI
import dolfinx as do
import safeincave.PostProcessingTools as post

# Units
MPa = 1e6



# Units
MPa = 1e6

def q_dil_rd(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0): # tensile strenght is changed from 1.5 to 5 to stay consistent with the desai model
    p = np.asarray(p_MPa, dtype=float)
    psi = np.asarray(psi, dtype=float)

    I1 = 3.0 * p
    sgn = np.sign(I1)
    sgn[sgn == 0.0] = 1.0

    denom = (np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi))
    # Prevent blow-up if denom gets very small
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)

    sqrtJ2 = D1 * ((I1 / (sgn * sigma_ref)) ** m) / denom + T0
    q = np.sqrt(3.0) * sqrtJ2
    return q




def compute_lode_angle_from_sigma(sig_Pa):
    # Enforce symmetry per cell (robust against tiny numerical asymmetry)
    sig = 0.5 * (sig_Pa + np.swapaxes(sig_Pa, -1, -2))

    vals = np.linalg.eigvalsh(sig)
    s1, s2, s3 = vals[:, 2], vals[:, 1], vals[:, 0]

    I1 = s1 + s2 + s3
    mean = I1 / 3.0
    s1d = s1 - mean
    s2d = s2 - mean
    s3d = s3 - mean

    J2 = (1.0 / 6.0) * ((s1d - s2d) ** 2 + (s2d - s3d) ** 2 + (s3d - s1d) ** 2)
    J3 = s1d * s2d * s3d

    tiny = 1e-30
    J2_safe = np.maximum(J2, tiny)

    x = (3.0 * np.sqrt(3.0) / 2.0) * (J3 / (J2_safe ** 1.5))
    x = np.clip(x, -1.0, 1.0)

    theta = (1.0 / 3.0) * np.arccos(x)   # <-- FIXED
    psi = theta - np.pi / 6.0            # map to [-pi/6, +pi/6]
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
    p_MPa = -p_vals  / MPa   # (nt, ncells)
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
    p_MPa = -p_vals / MPa       # (nt, ncells)
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



def write_FOS_serial(operation_folder, time_list, FOS, mesh_source_xdmf):
    """
    Serial writer: reads mesh and writes XDMF with COMM_SELF.
    Works even if the simulation was run in parallel.
    """
    # Read mesh in SERIAL communicator
    with do.io.XDMFFile(MPI.COMM_SELF, mesh_source_xdmf, "r") as xdmf:
        mesh = xdmf.read_mesh()
        mesh.name = "mesh"

    V = do.fem.functionspace(mesh, ("DG", 0))
    fos_fun = do.fem.Function(V)
    fos_fun.name = "FOS"

    fos_folder = os.path.join(operation_folder, "FOS")
    os.makedirs(fos_folder, exist_ok=True)
    fos_path = os.path.join(fos_folder, "FOS.xdmf")

    # Write in SERIAL communicator
    with do.io.XDMFFile(MPI.COMM_SELF, fos_path, "w") as xdmf_out:
        xdmf_out.write_mesh(mesh)
        for it, t in enumerate(time_list):
            fos_fun.x.array[:] = FOS[it, :]  # now DG0 dofs == global ncells in serial mesh
            xdmf_out.write_function(fos_fun, float(t))

    print(f"FOS written to: {fos_path}")





def main():
    # Only rank 0 postprocesses (MPI-safe, simplest)
    if MPI.COMM_WORLD.rank != 0:
        return

    base_folder = os.path.join("output", "case_sinus_70days_irregular", "operation")

    print("Loading stresses...")
    time_list, p_MPa, q_MPa, sig_vals, p_source = load_p_q_sig(base_folder)

    print("Computing FOS with Lode angle...")
    FOS = compute_FOS_with_lode(p_MPa, q_MPa, sig_vals)

    print("Writing FOS field...")
    write_FOS_serial(base_folder, time_list, FOS, p_source)



if __name__ == "__main__":
    main()
