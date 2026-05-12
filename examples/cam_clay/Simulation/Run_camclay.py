"""
Cyclic loading of a salt cavern in a clay-rich formation using the
Modified Cam-Clay viscoplastic element (`sf.ModifiedCamClayViscoplastic`).

This is a minimal driver, modelled on
`examples/mechanics/nobian/Simulation/Run.py` but trimmed to the
essentials needed to exercise the new constitutive element:
    1. equilibrium phase at lithostatic pressure (short)
    2. operation phase with a sinusoidal cavern pressure for N cycles

Default material parameters are taken from the notebook
`examples/cam_clay/Calibration_Clay_Rich.ipynb` (Cell 4, `table_guess`)
for the Brine 100 bar case — a clay-rich formation.
"""
import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import dolfinx as do
import os
import math
import torch as to


# ── USER CONFIGURATION ────────────────────────────────────────────────────────
# Grid (relative to this script)
GRID_FOLDER = "cavern_regular_1200_3D"

# Cavern/domain geometry (matches the salt grids)
Z_MAX = 660.0                # top of domain (m)
Z_CAVERN_CENTER = 450.0      # cavern center depth from surface convention (m)
SALT_DENSITY = 2200.0        # kg/m^3 (overburden + side rock)
GAS_DENSITY = 0.089          # kg/m^3 (light gas — H2)

# Pressure schedule
P_REF_MPA = 25.0             # reference at z = Z_MAX, lithostatic top (MPa)
P_MEAN_MPA = 12.0            # mean operating pressure (MPa)
P_AMPL_MPA = 4.0             # amplitude (MPa)  → P swings 8..16 MPa
PERIOD_HOURS = 24.0          # one cycle = 24 h
N_CYCLES = 8                 # operation cycles
EQUILIBRIUM_HOURS = 2.0      # short equilibrium hold before cycling
DT_HOURS = 0.5               # time step

# Cam-Clay parameters (notebook Brine 100 bar baseline)
PHI0   = 0.32      # initial porosity
LAM    = 0.00285   # compression index
KAP    = 0.00055   # swelling index
M_CSL  = 0.702     # CSL slope
THETA  = 0.0045    # softening parameter
PC0_PA = 50.0e6    # initial preconsolidation pressure (Pa)
                   # Notebook used 13 MPa for low-confining triaxial (sigma3=10 MPa);
                   # cavern lithostatic at z~450 m is ~35 MPa, so pc0 must be ≥ that
                   # for the in-situ state to start inside the yield surface.
ETA_V  = 1.0e11    # Perzyna viscosity
N_RATE = 2.0       # Perzyna exponent

# Elastic skeleton (clay-rich — softer than salt)
E0_GPA = 2.0
NU0    = 0.30


def main():
    # ── Domain & lithostatic reference ───────────────────────────────────────
    g = -9.81
    g_vec = [0.0, 0.0, g]
    p_ref = P_REF_MPA * ut.MPa
    side_burden = p_ref
    over_burden = p_ref

    # Output folder
    output_folder = os.path.join(
        "output",
        f"camclay_sinus({N_CYCLES})_{int(N_CYCLES * PERIOD_HOURS / 24)}days"
    )

    if MPI.COMM_WORLD.rank == 0:
        print("=" * 70)
        print("MODIFIED CAM-CLAY VISCOPLASTIC — FEM DRIVER")
        print("=" * 70)
        print(f"  Grid:           {GRID_FOLDER}")
        print(f"  P_ref:          {P_REF_MPA:.2f} MPa")
        print(f"  P swing:        {P_MEAN_MPA - P_AMPL_MPA:.2f} .. "
              f"{P_MEAN_MPA + P_AMPL_MPA:.2f} MPa")
        print(f"  Period:         {PERIOD_HOURS:.1f} h × {N_CYCLES} cycles")
        print(f"  MCC params:     phi0={PHI0}, lam={LAM}, kap={KAP}, M={M_CSL}, "
              f"theta={THETA}, pc0={PC0_PA/1e6:.1f} MPa, eta_v={ETA_V:.1e}, "
              f"n={N_RATE}")
        print(f"  Output:         {output_folder}")
        print("=" * 70)

    # ── Grid ──────────────────────────────────────────────────────────────────
    grid_path = os.path.join("..", "..", "..", "grids", GRID_FOLDER)
    grid = sf.GridHandlerGMSH("geom", grid_path)

    mom_eq = sf.LinearMomentum(grid, theta=0.5)
    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("cg")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-10, max_it=100)
    mom_eq.set_solver(mom_solver)

    # ── Material: elastic spring + Modified Cam-Clay viscoplastic ───────────
    mat = sf.Material(mom_eq.n_elems)
    rho = SALT_DENSITY * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    E0 = E0_GPA * ut.GPa * to.ones(mom_eq.n_elems, dtype=to.float64)
    nu0 = NU0 * to.ones(mom_eq.n_elems, dtype=to.float64)
    spring = sf.Spring(E0, nu0, "spring")
    mat.add_to_elastic(spring)

    # Void ratio from porosity: e = phi / (1 - phi)
    e0_val = PHI0 / (1.0 - PHI0)

    mcc = sf.ModifiedCamClayViscoplastic(
        M=     M_CSL  * to.ones(mom_eq.n_elems, dtype=to.float64),
        lam=   LAM    * to.ones(mom_eq.n_elems, dtype=to.float64),
        kap=   KAP    * to.ones(mom_eq.n_elems, dtype=to.float64),
        theta= THETA  * to.ones(mom_eq.n_elems, dtype=to.float64),
        pc0=   PC0_PA * to.ones(mom_eq.n_elems, dtype=to.float64),
        e0=    e0_val * to.ones(mom_eq.n_elems, dtype=to.float64),
        eta_v= ETA_V  * to.ones(mom_eq.n_elems, dtype=to.float64),
        n_rate=N_RATE * to.ones(mom_eq.n_elems, dtype=to.float64),
        name="cam_clay",
    )
    mat.add_to_non_elastic(mcc)
    mom_eq.set_material(mat)
    mom_eq.build_body_force(g_vec)
    mom_eq.expect_vp_state = True

    # ── Equilibrium phase ───────────────────────────────────────────────────
    tc_init = sf.TimeController(
        dt=DT_HOURS, initial_time=0.0,
        final_time=EQUILIBRIUM_HOURS, time_unit="hour",
    )
    p_gas_init = P_MEAN_MPA * ut.MPa
    t_init = [0.0, tc_init.t_final]
    p_init = [p_gas_init, p_gas_init]

    bc_init = momBC.BcHandler(mom_eq)
    for bc in _build_boundary_conditions(tc_init.t_final, side_burden,
                                         over_burden, t_init, p_init,
                                         g_vec):
        bc_init.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_init)

    out_init = sf.SaveFields(mom_eq)
    out_init.set_output_folder(os.path.join(output_folder, "equilibrium"))
    for f, lbl in [("u", "Displacement (m)"),
                   ("eps_tot", "Total strain (-)"),
                   ("sig", "Stress (Pa)"),
                   ("p_elems", "Mean stress (Pa)"),
                   ("q_elems", "Von Mises stress (Pa)")]:
        out_init.add_output_field(f, lbl)
    os.makedirs(output_folder, exist_ok=True)

    sim_init = sf.Simulator_M(mom_eq, tc_init, [out_init], True)
    sim_init.run()

    if MPI.COMM_WORLD.rank == 0:
        print("[EQUILIBRIUM] complete — starting cyclic operation phase")

    # ── Operation phase (sinusoidal cavern pressure) ────────────────────────
    total_hours = N_CYCLES * PERIOD_HOURS
    tc_op = sf.TimeController(
        dt=DT_HOURS, initial_time=tc_init.t_final,
        final_time=tc_init.t_final + total_hours, time_unit="hour",
    )

    # Build sinus schedule sampled at dt
    n_steps = int(round(total_hours / DT_HOURS)) + 1
    t_op = [tc_init.t_final + i * DT_HOURS for i in range(n_steps)]
    p_mean = P_MEAN_MPA * ut.MPa
    p_amp = P_AMPL_MPA * ut.MPa
    omega = 2.0 * math.pi / PERIOD_HOURS
    p_op = [p_mean - p_amp * math.cos(omega * (t - tc_init.t_final))
            for t in t_op]

    bc_op = momBC.BcHandler(mom_eq)
    for bc in _build_boundary_conditions(tc_op.t_final, side_burden,
                                         over_burden, t_op, p_op, g_vec):
        bc_op.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_op)

    out_op = sf.SaveFields(mom_eq)
    out_op.set_output_folder(os.path.join(output_folder, "operation"))
    for f, lbl in [("u", "Displacement (m)"),
                   ("eps_tot", "Total strain (-)"),
                   ("sig", "Stress (Pa)"),
                   ("p_elems", "Mean stress (Pa)"),
                   ("q_elems", "Von Mises stress (Pa)")]:
        out_op.add_output_field(f, lbl)

    sim_op = sf.Simulator_M(mom_eq, tc_op, [out_op], False)
    sim_op.run()

    if MPI.COMM_WORLD.rank == 0:
        print("[OPERATION] complete.")
        print(f"  Output written to {output_folder}/")


def _build_boundary_conditions(t_final, side_burden, over_burden,
                               t_press, p_press, g_vec):
    """
    Standard symmetry-quadrant BCs for the cavern grid:
        West   — u_x = 0   (Dirichlet, dof 0)
        South  — u_y = 0   (Dirichlet, dof 1)
        Bottom — u_z = 0   (Dirichlet, dof 2)
        East   — sigma_xx = side_burden (Neumann, dof 0 axis, height-varying)
        North  — sigma_yy = side_burden (Neumann, dof 1 axis, height-varying)
        Top    — sigma_zz = over_burden (Neumann, dof 2 axis, constant)
        Cavern — sigma_n  = p_gas(t)
    """
    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, t_final])

    bc_east = momBC.NeumannBC("East", 2, SALT_DENSITY, Z_MAX,
                              [side_burden, side_burden],
                              [0.0, t_final], g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, SALT_DENSITY, Z_MAX,
                               [side_burden, side_burden],
                               [0.0, t_final], g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                             [over_burden, over_burden],
                             [0.0, t_final], g=g_vec[2])
    bc_cavern = momBC.NeumannBC("Cavern", 2, GAS_DENSITY, Z_MAX,
                                p_press, t_press, g=g_vec[2])

    return [bc_west, bc_south, bc_bottom, bc_east, bc_north, bc_top, bc_cavern]


if __name__ == "__main__":
    main()
