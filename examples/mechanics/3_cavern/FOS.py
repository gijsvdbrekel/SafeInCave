import os
import numpy as np
from mpi4py import MPI
import dolfinx as do
import safeincave.PostProcessingTools as post

# Units
MPa = 1e6

# ---------------------------------------------------------------------
# RD dilatancy boundary  (same parameters as your plotting function)
# ---------------------------------------------------------------------
def q_dil_rd(p_MPa,
             D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    """
    Compute q_dil(p) in MPa for the Rider–Dieterich dilatancy boundary,
    compression branch (psi = -30°), using same form as plot_dilatancy_boundary.
    p_MPa must be compression-positive (as in your plots).
    """
    p = np.asarray(p_MPa, dtype=float)
    # mean stress -> I1
    I1 = 3.0 * p  # MPa

    psi_c = -np.pi / 6.0  # triaxial compression branch

    # sign trick as in your plotting code
    sgn = np.sign(I1)
    sgn[sgn == 0.0] = 1.0

    denom = (np.sqrt(3.0) * np.cos(psi_c) - D2 * np.sin(psi_c))
    sqrtJ2 = D1 * ((I1 / (sgn * sigma_ref)) ** m) / denom + T0  # MPa

    q = np.sqrt(3.0) * sqrtJ2  # MPa
    return q


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
    p_MPa = -p_vals / MPa   # (nt, ncells)
    q_MPa =  q_vals / MPa   # (nt, ncells)

    return time_list, p_MPa, q_MPa, p_path


# ---------------------------------------------------------------------
# Compute FOS = q_dil(p) / q
# ---------------------------------------------------------------------
def compute_FOS(p_MPa, q_MPa):
    """
    p_MPa, q_MPa: (nt, ncells)
    returns FOS with same shape.
    """
    nt, nc = p_MPa.shape
    FOS = np.zeros_like(p_MPa)

    for it in range(nt):
        p_t = p_MPa[it, :]
        q_t = q_MPa[it, :]

        q_dil = q_dil_rd(p_t)  # MPa
        q_safe = np.where(q_t <= 0.0, 1e-12, q_t)  # avoid 0 division

        FOS[it, :] = q_dil / q_safe

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




# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------
def main():
    # Must match your Test.py output folder
    base_folder = os.path.join("output", "case_sin_1day_asymmetric", "operation")

    if MPI.COMM_WORLD.rank == 0:
        print("Loading stresses...")
    time_list, p_MPa, q_MPa, p_source = load_p_q(base_folder)

    if MPI.COMM_WORLD.rank == 0:
        print("Computing Factor of Safety...")
    FOS = compute_FOS(p_MPa, q_MPa)

    if MPI.COMM_WORLD.rank == 0:
        print("Writing FOS field...")
    write_FOS(base_folder, time_list, FOS, p_source)


if __name__ == "__main__":
    main()
