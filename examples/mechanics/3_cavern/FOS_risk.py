import os
import numpy as np
from mpi4py import MPI
import dolfinx as do
import safeincave.PostProcessingTools as post

# Units
MPa = 1e6

# ---------------------------------------------------------------------
# RD dilatancy boundary
# ---------------------------------------------------------------------
def q_dil_rd(p_MPa,
             D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    """
    Compute q_dil(p) in MPa for the Rider–Dieterich dilatancy boundary,
    compression branch (psi = -30°).
    p_MPa must be compression-positive.
    """
    p = np.asarray(p_MPa, dtype=float)
    I1 = 3.0 * p  # MPa

    psi_c = -np.pi / 6.0  # triaxial compression branch

    # zorg dat sign array is (werkt ook voor scalars)
    sgn = np.sign(I1)
    sgn[sgn == 0.0] = 1.0

    denom = (np.sqrt(3.0) * np.cos(psi_c) - D2 * np.sin(psi_c))
    sqrtJ2 = D1 * ((I1 / (sgn * sigma_ref)) ** m) / denom + T0  # MPa

    q = np.sqrt(3.0) * sqrtJ2  # MPa
    return q

# ---------------------------------------------------------------------
# Load p and q
# ---------------------------------------------------------------------
def load_p_q(operation_folder):
    """
    Returns:
        time_list : (nt,) array of times
        p_MPa     : (nt, ncells_local) mean stress, compression-positive [MPa]
        q_MPa     : (nt, ncells_local) von Mises stress [MPa]
        p_path    : path to p_elems.xdmf (for mesh)
    """
    p_path = os.path.join(operation_folder, "p_elems", "p_elems.xdmf")
    q_path = os.path.join(operation_folder, "q_elems", "q_elems.xdmf")

    pts_p, time_list, p_vals = post.read_cell_scalar(p_path)
    pts_q, time_list2, q_vals = post.read_cell_scalar(q_path)

    if not np.allclose(time_list, time_list2):
        raise RuntimeError("Time lists of p_elems and q_elems do not match.")

    # zelfde conventie als stress path: compressie → positieve p
    p_MPa = -p_vals / MPa   # (nt, ncells_local)
    q_MPa =  q_vals / MPa   # (nt, ncells_local)

    return time_list, p_MPa, q_MPa, p_path

# ---------------------------------------------------------------------
# Compute FOS + RiskLevel
# ---------------------------------------------------------------------
def compute_FOS_and_risk(p_MPa, q_MPa,
                         fos_crit=1.0,
                         near_factor=1.2):
    """
    p_MPa, q_MPa: (nt, ncells_local)
    fos_crit:     FOS <= fos_crit  => tertiary creep
    near_factor:  fos_crit < FOS <= fos_crit*near_factor => near-tertiary

    Returns:
        FOS       : (nt, ncells_local)
        RiskLevel : (nt, ncells_local), int
                    0 = safe
                    1 = near-tertiary
                    2 = tertiary
    """
    nt, nc = p_MPa.shape
    FOS = np.zeros_like(p_MPa)
    RiskLevel = np.zeros_like(p_MPa, dtype=np.int32)

    for it in range(nt):
        p_t = p_MPa[it, :]
        q_t = q_MPa[it, :]

        # optioneel: p naar klein positief minimum clippen
        p_t_clip = np.maximum(p_t, 1e-3)
        q_dil = q_dil_rd(p_t_clip)  # MPa

        q_safe = np.where(q_t <= 0.0, 1e-12, q_t)  # avoid 0 division

        FOS_t = q_dil / q_safe
        FOS[it, :] = FOS_t

        # Risk levels
        tertiary_mask = FOS_t <= fos_crit
        near_mask = (FOS_t > fos_crit) & (FOS_t <= fos_crit * near_factor)
        # safe: default 0
        RiskLevel[it, tertiary_mask] = 2
        RiskLevel[it, near_mask] = 1

    return FOS, RiskLevel

# ---------------------------------------------------------------------
# Write FOS + RiskLevel as DG0 cell fields
# ---------------------------------------------------------------------
def write_FOS_and_risk(operation_folder, time_list,
                       FOS, RiskLevel, mesh_source_xdmf):
    print("Writing FOS and RiskLevel fields...\n")

    # Mesh uit p_elems.xdmf
    with do.io.XDMFFile(MPI.COMM_WORLD, mesh_source_xdmf, "r") as xdmf:
        mesh = xdmf.read_mesh(name="mesh")
        mesh.name = "mesh"

    V = do.fem.functionspace(mesh, ("DG", 0))
    fos_fun = do.fem.Function(V)
    fos_fun.name = "FOS"
    risk_fun = do.fem.Function(V)
    risk_fun.name = "RiskLevel"

    # sanity check voor serieel/parallel
    local_size = fos_fun.x.array.size
    assert FOS.shape[1] == local_size, (
        f"Shape mismatch: FOS has {FOS.shape[1]} cols, "
        f"but DG0 has {local_size} local dofs on this rank."
    )

    fos_folder = os.path.join(operation_folder, "FOS")
    os.makedirs(fos_folder, exist_ok=True)
    fos_path = os.path.join(fos_folder, "FOS_and_Risk.xdmf")

    with do.io.XDMFFile(MPI.COMM_WORLD, fos_path, "w") as xdmf_out:
        xdmf_out.write_mesh(mesh)

        nt = len(time_list)
        for it, t in enumerate(time_list):
            fos_fun.x.array[:] = FOS[it, :]
            risk_fun.x.array[:] = RiskLevel[it, :]

            # schrijf beide velden voor dezelfde tijdstap
            xdmf_out.write_function(fos_fun, float(t))
            xdmf_out.write_function(risk_fun, float(t))

    if MPI.COMM_WORLD.rank == 0:
        print(f"FOS and RiskLevel written to: {fos_path}")

# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------
def main():
    base_folder = os.path.join("output", "case_sinus_70days_regular", "operation")

    if MPI.COMM_WORLD.rank == 0:
        print("Loading stresses...")
    time_list, p_MPa, q_MPa, p_source = load_p_q(base_folder)

    if MPI.COMM_WORLD.rank == 0:
        print("Computing FOS and Risk levels...")
    FOS, RiskLevel = compute_FOS_and_risk(
        p_MPa, q_MPa,
        fos_crit=1.0,     # FOS <= 1.0 -> tertiary
        near_factor=1.2   # 1.0 < FOS <= 1.2 -> near-tertiary
    )

    if MPI.COMM_WORLD.rank == 0:
        print("Writing fields for ParaView...")
    write_FOS_and_risk(base_folder, time_list, FOS, RiskLevel, p_source)


if __name__ == "__main__":
    main()
