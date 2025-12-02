import os
import numpy as np
import dolfinx as do
from mpi4py import MPI
from dolfinx.io import XDMFFile
import ufl
import safeincave.PostProcessingTools as post

# Unit conversions
MPa = 1e6      # Pa
hour = 3600.0  # s


# ============================================================
#  RD DILATANCY MODEL (same as in your plotting script)
# ============================================================

def q_dilatancy(p_MPa, *, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    """
    Compute q_dil(p) for the Moss Bluff RD dilatancy surface.
    p_MPa must be >= 0 for compression.
    Returns q_dil in MPa.
    Uses the compression branch (psi = -30°).
    """
    p = np.maximum(p_MPa, 0.0)  # negative or small p → clipped for safety
    I1 = 3.0 * p                # MPa

    psi_c = -np.pi/6.0  # triaxial compression

    # denominator
    denom = (np.sqrt(3.0)*np.cos(psi_c) - D2*np.sin(psi_c))

    sgn = np.sign(I1)
    sgn[sgn == 0.0] = 1.0

    sqrtJ2 = D1 * ( (I1/(sgn*sigma_ref))**m ) / denom + T0       # MPa
    q = np.sqrt(3.0) * sqrtJ2                                    # MPa

    return q


# ============================================================
#   MAIN FOS CALCULATION
# ============================================================

def compute_FOS(operation_folder):
    """
    Reads p_elems and q_elems from XDMF, computes FOS for each cell and time step.
    Returns:
        times     (numpy array of floats)
        FOS_all   (array [num_times, num_cells])
    """

    # ------------- Read DG0 cell pressure and deviatoric stress -------------
    p_path = os.path.join(operation_folder, "p_elems", "p_elems.xdmf")
    q_path = os.path.join(operation_folder, "q_elems", "q_elems.xdmf")

    points_p, time_list, p_elems = post.read_cell_scalar(p_path)
    points_q, time_list2, q_elems = post.read_cell_scalar(q_path)

    assert np.allclose(time_list, time_list2), "Time lists in p and q fields do not match."

    num_times = len(time_list)
    num_cells = p_elems.shape[1]

    # ------------------ Allocate output ------------------------
    FOS = np.zeros((num_times, num_cells))

    # ------------------ Compute FOS ----------------------------
    for t in range(num_times):

        # Convert to MPa and flip sign for p (compression positive)
        p_MPa = -p_elems[t] / MPa
        q_MPa =  q_elems[t] / MPa

        # Dilatancy boundary q_dil(p)
        q_dil = q_dilatancy(p_MPa)       # MPa

        # Avoid division by zero
        q_MPa_safe = np.maximum(q_MPa, 1e-12)

        # Factor of Safety
        FOS[t,:] = q_dil / q_MPa_safe

    return time_list, FOS



# ============================================================
#   WRITE FOS TO XDMF FOR PARAVIEW
# ============================================================

def write_FOS_to_xdmf(operation_folder, time_list, FOS):
    """
    Writes FOS as a DG0 cell field into:
        operation/FOS/FOS.xdmf
    """

    mesh_path = os.path.join(operation_folder, "mesh", "geom.xdmf")
    mesh = do.io.XDMFFile(MPI.COMM_WORLD, mesh_path, "r").read_mesh(name="Grid")

    V = do.fem.FunctionSpace(mesh, ("DG", 0))
    FOS_field = do.fem.Function(V)
    FOS_field.name = "FOS"

    fos_folder = os.path.join(operation_folder, "FOS")
    os.makedirs(fos_folder, exist_ok=True)
    out_path = os.path.join(fos_folder, "FOS.xdmf")

    with XDMFFile(mesh.comm, out_path, "w") as xdmf:
        xdmf.write_mesh(mesh)

        for i, t in enumerate(time_list):
            if i % 10 == 0 and mesh.comm.rank == 0:
                print(f"Writing timestep {i+1}/{len(time_list)} ...")

            FOS_field.x.array[:] = FOS[i,:]
            xdmf.write_function(FOS_field, t)


# ============================================================
#   MAIN ENTRY POINT
# ============================================================

def main():
    # Folder that contains p_elems, q_elems, mesh, etc.
    output_folder = os.path.join("output", "case_sin_1day_asymmetric", "operation")

    if MPI.COMM_WORLD.rank == 0:
        print("Computing Factor of Safety...\n")

    time_list, FOS = compute_FOS(output_folder)

    if MPI.COMM_WORLD.rank == 0:
        print("Writing FOS to XDMF...\n")

    write_FOS_to_xdmf(output_folder, time_list, FOS)

    if MPI.COMM_WORLD.rank == 0:
        print("Done. FOS.xdmf is ready for ParaView.")


if __name__ == "__main__":
    main()
