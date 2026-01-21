import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import dolfinx as do
import os
import torch as to
import json
import math


# ============================================================
#  SCENARIO SWITCHES (zet precies ÉÉN preset True)
# ============================================================
RUN_DESAI_ONLY = True            # spring + desai (alleen viscoplastisch)
RUN_DISLOC_OLD_ONLY = False       # spring + dislocation (oude params)
RUN_DISLOC_NEW_ONLY = False       # spring + dislocation (nieuwe params)
RUN_FULL = False                   # spring + kelvin + disloc + pressure-solution + desai
RUN_FULL_MINUS_DESAI = False      # spring + kelvin + disloc + pressure-solution (zonder desai)

# Equilibrium: meestal wil je puur elastisch initialiseren (aanrader)
EQUILIBRIUM_ELASTIC_ONLY = True

# ============================================================
#  FIXED RUN SETTINGS (zoals jij vroeg)
# ============================================================
CAVERN_TYPE = "regular600"
OPERATION_DAYS = 365
N_CYCLES = 10
dt_hours = 2
PRESSURE_SCENARIO = "sinus"   # fixed
SCHEDULE_MODE = "stretch"     # fixed (N_CYCLES over totaalduur)
p_gas_MPa_equilibrium = 15.0  # keep equilibrium consistent


# ============================================================
#  Helpers
# ============================================================
DAY_H = 24.0  # hours per day


def build_sinus_pressure_schedule(tc, *, p_mean, p_ampl, period_hours, phase_hours=0.0,
                                  clamp_min=None, clamp_max=None):
    period = period_hours * ut.hour
    phase = phase_hours * ut.hour

    n_steps = int(math.floor(tc.t_final / tc.dt))
    t_vals = [k * tc.dt for k in range(n_steps + 1)]
    if abs(t_vals[-1] - tc.t_final) > 1e-12:
        t_vals.append(tc.t_final)

    two_pi_over_T = (2.0 * math.pi / period) if period > 0.0 else 0.0

    p_vals = []
    for t in t_vals:
        p = p_mean if period <= 0.0 else p_mean + p_ampl * math.sin(two_pi_over_T * (t - phase))
        if clamp_min is not None:
            p = max(p, clamp_min)
        if clamp_max is not None:
            p = min(p, clamp_max)
        p_vals.append(p)

    return t_vals, p_vals


def build_sinus_schedule_multi(tc, *, p_mean, p_ampl, days, mode,
                               daily_period_hours=24.0,
                               total_cycles=1,
                               clamp_min=None, clamp_max=None):
    total_hours = days * DAY_H

    if mode == "repeat":
        T_hours = daily_period_hours
    elif mode == "stretch":
        total_cycles = max(1, int(total_cycles))
        T_hours = total_hours / float(total_cycles)
    else:
        raise ValueError("mode must be 'repeat' or 'stretch'")

    return build_sinus_pressure_schedule(
        tc, p_mean=p_mean, p_ampl=p_ampl,
        period_hours=T_hours, phase_hours=0.0,
        clamp_min=clamp_min, clamp_max=clamp_max
    )


def _scenario_name():
    if RUN_DESAI_ONLY:
        return "desai_only"
    if RUN_DISLOC_OLD_ONLY:
        return "disloc_old_only"
    if RUN_DISLOC_NEW_ONLY:
        return "disloc_new_only"
    if RUN_FULL:
        return "full"
    if RUN_FULL_MINUS_DESAI:
        return "full_minus_desai"
    return "unknown"


def _validate_scenario():
    flags = [
        RUN_DESAI_ONLY,
        RUN_DISLOC_OLD_ONLY,
        RUN_DISLOC_NEW_ONLY,
        RUN_FULL,
        RUN_FULL_MINUS_DESAI,
    ]
    if sum(bool(x) for x in flags) != 1:
        raise ValueError(
            "Zet precies één scenario True: "
            "RUN_DESAI_ONLY, RUN_DISLOC_OLD_ONLY, RUN_DISLOC_NEW_ONLY, RUN_FULL, RUN_FULL_MINUS_DESAI"
        )


class LinearMomentumMod(sf.LinearMomentum):
    def __init__(self, grid, theta):
        super().__init__(grid, theta)
        self.expect_vp_state = False

    def initialize(self) -> None:
        self.C.x.array[:] = to.flatten(self.mat.C)
        self.Fvp = do.fem.Function(self.DG0_1)
        self.alpha = do.fem.Function(self.DG0_1)
        self.eps_vp = do.fem.Function(self.DG0_3x3)

    def run_after_solve(self):
        if not hasattr(self, "eps_vp"):
            return

        elems = getattr(self.mat, "elems_ne", None)
        if not elems:
            return
        st = elems[-1]

        # eps_vp (werkt als je inelastic update eps_ne_k opslaat)
        if hasattr(st, "eps_ne_k"):
            self.eps_vp.x.array[:] = to.flatten(st.eps_ne_k)

        # Fvp/alpha alleen zinvol wanneer Desai actief is
        if self.expect_vp_state:
            if hasattr(st, "Fvp") and hasattr(st, "alpha"):
                self.Fvp.x.array[:] = st.Fvp
                self.alpha.x.array[:] = st.alpha


class SparseSaveFields(sf.SaveFields):
    """SaveFields dat alleen elke `interval`-ste stap schrijft (t=0 altijd)."""
    def __init__(self, mom_eq, interval: int):
        super().__init__(mom_eq)
        self.interval = max(1, int(interval))
        self._counter = 0

    def save_fields(self, t):
        if t == 0:
            return super().save_fields(t)
        self._counter += 1
        if self._counter % self.interval == 0:
            return super().save_fields(t)


def make_solver(grid):
    ksp = PETSc.KSP().create(grid.mesh.comm)
    ksp.setType("cg")
    ksp.getPC().setType("asm")
    ksp.setTolerances(rtol=1e-12, max_it=100)
    return ksp


def build_base_material_parts(mom_eq):
    """Definieer density + alle element-objects; kies later wat je toevoegt."""
    salt_density = 2200
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=to.float64)

    # Elastic spring
    E0 = 20.425 * ut.GPa * to.ones(mom_eq.n_elems)
    nu0 = 0.25 * to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")

    # Kelvin-Voigt
    eta = 105e11 * to.ones(mom_eq.n_elems)
    E1 = 10 * ut.GPa * to.ones(mom_eq.n_elems)
    nu1 = 0.25 * to.ones(mom_eq.n_elems)
    kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")

    sec_per_year = 365.25 * 24 * 3600

    # Dislocation creep (oude params)
    ndc = 4.6
    A_dc = (40.0 * (1e-6) ** ndc / sec_per_year) * to.ones(mom_eq.n_elems)
    Q_dc = (6495.0 * 8.32) * to.ones(mom_eq.n_elems)
    n_dc = (ndc) * to.ones(mom_eq.n_elems)
    creep_disloc_old = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation_old")

    # Dislocation creep (nieuwe params)
    A_new = 1.9e-20 * to.ones(mom_eq.n_elems)
    Q_new = 51600 * to.ones(mom_eq.n_elems)
    n_new = 3.0 * to.ones(mom_eq.n_elems)
    creep_disloc_new = sf.DislocationCreep(A_new, Q_new, n_new, "creep_dislocation_new")

    # Pressure-solution creep
    A_ps = (14176.0 * 1e-9 / 1e6 / sec_per_year) * to.ones(mom_eq.n_elems)
    d_ps = 5.25e-3 * to.ones(mom_eq.n_elems)
    Q_ps = (3252.0 * 8.32) * to.ones(mom_eq.n_elems)
    creep_pressure = sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "creep_pressure")

    # Desai viscoplastic
    mu_1 = 5.3665857009859815e-11 * to.ones(mom_eq.n_elems)
    N_1 = 3.1 * to.ones(mom_eq.n_elems)
    n = 3.0 * to.ones(mom_eq.n_elems)
    a_1 = 1.965018496922832e-05 * to.ones(mom_eq.n_elems)
    eta_vp = 0.8275682807874163 * to.ones(mom_eq.n_elems)
    beta_1 = 0.0048 * to.ones(mom_eq.n_elems)
    beta = 0.995 * to.ones(mom_eq.n_elems)
    m = -0.5 * to.ones(mom_eq.n_elems)
    gamma = 0.095 * to.ones(mom_eq.n_elems)
    alpha_0 = 0.0022 * to.ones(mom_eq.n_elems)
    sigma_t = 5.0 * to.ones(mom_eq.n_elems)

    desai = sf.ViscoplasticDesai(
        mu_1, N_1, a_1, eta_vp, n, beta_1, beta, m, gamma, sigma_t, alpha_0, "desai"
    )

    return {
        "salt_density": salt_density,
        "rho": rho,
        "spring": spring_0,
        "kelvin": kelvin,
        "disloc_old": creep_disloc_old,
        "disloc_new": creep_disloc_new,
        "pressure_solution": creep_pressure,
        "desai": desai,
    }


def make_material(mom_eq, parts, *, use_kelvin, use_disloc_old, use_disloc_new, use_pressure_solution, use_desai):
    mat = sf.Material(mom_eq.n_elems)
    mat.set_density(parts["rho"])
    mat.add_to_elastic(parts["spring"])

    if use_kelvin:
        mat.add_to_non_elastic(parts["kelvin"])

    # Zorg dat je niet per ongeluk beide dislocation creeps tegelijk toevoegt
    if use_disloc_old and use_disloc_new:
        raise ValueError("Kies óf disloc_old óf disloc_new, niet allebei.")

    if use_disloc_old:
        mat.add_to_non_elastic(parts["disloc_old"])
    if use_disloc_new:
        mat.add_to_non_elastic(parts["disloc_new"])

    if use_pressure_solution:
        mat.add_to_non_elastic(parts["pressure_solution"])

    if use_desai:
        mat.add_to_non_elastic(parts["desai"])

    return mat


def main():
    _validate_scenario()

    grid_path = os.path.join("..", "..", "..", "grids", "cavern_regular_600_3D")
    grid = sf.GridHandlerGMSH("geom", grid_path)

    Z_MAX_BY_CAVERN = {
        "regular600": 315.26,
        "tilted600": 345.67,
        "teardrop600": 353.15,
        "asymmetric600": 338.89,
        "irregular600": 319.86,
        "multichamber600": 334.14,
        "regular1200": 393.21,
        "tilted1200":  430.78,
        "teardrop1200":  445.06,
        "asymmetric1200":  422.76,
        "irregular1200":  402.21,
        "multichamber1200":  420.82,
    }
    z_max = Z_MAX_BY_CAVERN[CAVERN_TYPE]

    # Overburden & side burden
    if CAVERN_TYPE.endswith("600"):
        p_ref = 18.2 * ut.MPa
    elif CAVERN_TYPE.endswith("1200"):
        p_ref = 20.1 * ut.MPa
    else:
        raise ValueError(f"Cannot infer cavern set (600/1200) from CAVERN_TYPE='{CAVERN_TYPE}'")

    side_burden = p_ref
    over_burden = p_ref

    # Momentum equation + solver
    mom_eq = LinearMomentumMod(grid, theta=0.5)
    mom_eq.set_solver(make_solver(grid))

    # Build all material parts once
    parts = build_base_material_parts(mom_eq)
    salt_density = parts["salt_density"]

   # Body forces (PAS NA set_material!)
    g = -9.81
    g_vec = [0.0, 0.0, g]

    # Temperature
    T0_field = 298 * to.ones(mom_eq.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)


    # ============================================================
    # EQUILIBRIUM STAGE
    # ============================================================
    tc_equilibrium = sf.TimeController(dt=0.5, initial_time=0.0, final_time=10, time_unit="hour")

    # Equilibrium material: meestal puur elastisch
    if EQUILIBRIUM_ELASTIC_ONLY:
        mat_eq = make_material(
            mom_eq, parts,
            use_kelvin=False, use_disloc_old=False, use_disloc_new=False,
            use_pressure_solution=False, use_desai=False
        )
        mom_eq.expect_vp_state = False
    else:
        # als je tóch inelastic in equilibrium wil (niet aangeraden voor isolatie-tests)
        mat_eq = make_material(
            mom_eq, parts,
            use_kelvin=True, use_disloc_old=True, use_disloc_new=False,
            use_pressure_solution=True, use_desai=False
        )
        mom_eq.expect_vp_state = False

    mom_eq.set_material(mat_eq)
    mom_eq.build_body_force(g_vec)   # <-- HIER


    # Equilibrium BCs
    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_equilibrium.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_equilibrium.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_equilibrium.t_final])

    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0,
                              [side_burden, side_burden],
                              [0.0, tc_equilibrium.t_final],
                              g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0,
                               [side_burden, side_burden],
                               [0.0, tc_equilibrium.t_final],
                               g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                             [over_burden, over_burden],
                             [0.0, tc_equilibrium.t_final],
                             g=g_vec[2])

    gas_density = 0.089
    p_gas = p_gas_MPa_equilibrium * ut.MPa

    bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, z_max,
                                [p_gas, p_gas],
                                [0.0, tc_equilibrium.t_final],
                                g=g_vec[2])

    bc_equilibrium = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
        bc_equilibrium.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_equilibrium)

    scenario = _scenario_name()
    output_folder = os.path.join(
        "output",
        f"case_{scenario}_{PRESSURE_SCENARIO}_{N_CYCLES}cyc_{OPERATION_DAYS}days_{CAVERN_TYPE}"
    )

    output_folder_equilibrium = os.path.join(output_folder, "equilibrium")
    if MPI.COMM_WORLD.rank == 0:
        print("Equilibrium out:", output_folder_equilibrium)

    out_eq = sf.SaveFields(mom_eq)
    out_eq.set_output_folder(output_folder_equilibrium)
    out_eq.add_output_field("u", "Displacement (m)")
    out_eq.add_output_field("eps_tot", "Total strain (-)")
    out_eq.add_output_field("sig", "Stress (Pa)")
    out_eq.add_output_field("p_elems", "Mean stress (Pa)")
    out_eq.add_output_field("q_elems", "Von Mises stress (Pa)")

    sim_eq = sf.Simulator_M(mom_eq, tc_equilibrium, [out_eq], True)
    sim_eq.run()

    # ============================================================
    # OPERATION STAGE (scenario material + sinus schedule)
    # ============================================================
    tc_operation = sf.TimeController(
        dt=dt_hours, initial_time=0.0,
        final_time=OPERATION_DAYS * 24.0,
        time_unit="hour"
    )

    # Pressure schedule (fixed sinus, stretch, N_CYCLES)
    p_mean = p_gas_MPa_equilibrium * ut.MPa
    p_ampl = 6.5 * ut.MPa

    t_pressure, p_pressure = build_sinus_schedule_multi(
        tc_operation,
        p_mean=p_mean, p_ampl=p_ampl,
        days=OPERATION_DAYS, mode=SCHEDULE_MODE,
        daily_period_hours=24.0,
        total_cycles=N_CYCLES,
        clamp_min=None, clamp_max=None
    )

    # Choose scenario toggles
    use_kelvin = False
    use_disloc_old = False
    use_disloc_new = False
    use_pressure_solution = False
    use_desai = False

    if RUN_DESAI_ONLY:
        use_desai = True

    elif RUN_DISLOC_OLD_ONLY:
        use_disloc_old = True

    elif RUN_DISLOC_NEW_ONLY:
        use_disloc_new = True

    elif RUN_FULL:
        use_kelvin = True
        use_disloc_old = True
        use_pressure_solution = True
        use_desai = True

    elif RUN_FULL_MINUS_DESAI:
        use_kelvin = True
        use_disloc_old = True
        use_pressure_solution = True
        use_desai = False

    # Build operation material
    mat_op = make_material(
        mom_eq, parts,
        use_kelvin=use_kelvin,
        use_disloc_old=use_disloc_old,
        use_disloc_new=use_disloc_new,
        use_pressure_solution=use_pressure_solution,
        use_desai=use_desai
    )
    mom_eq.set_material(mat_op)
    mom_eq.expect_vp_state = bool(use_desai)

    # If Desai is active: initialize hardening from equilibrium stress state
    if use_desai:
        stress_to = ut.numpy2torch(mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3)))
        parts["desai"].compute_initial_hardening(stress_to, Fvp_0=0.0)

    # Operation BCs (same as equilibrium but cavern pressure is time series)
    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_operation.t_final])

    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0,
                              [side_burden, side_burden],
                              [0.0, tc_operation.t_final],
                              g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0,
                               [side_burden, side_burden],
                               [0.0, tc_operation.t_final],
                               g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                             [over_burden, over_burden],
                             [0.0, tc_operation.t_final],
                             g=g_vec[2])

    bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, z_max,
                                p_pressure,
                                t_pressure,
                                g=g_vec[2])

    bc_operation = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
        bc_operation.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_operation)

    # Save pressure schedule once
    os.makedirs(output_folder, exist_ok=True)
    pressure_data = {
        "scenario": scenario,
        "pressure_scenario": PRESSURE_SCENARIO,
        "mode": SCHEDULE_MODE,
        "operation_days": OPERATION_DAYS,
        "n_cycles": N_CYCLES,
        "dt_hours": dt_hours,
        "t_values": [float(t) for t in t_pressure],
        "p_values": [float(p) for p in p_pressure],
        "active_elements": {
            "kelvin": bool(use_kelvin),
            "dislocation_old": bool(use_disloc_old),
            "dislocation_new": bool(use_disloc_new),
            "pressure_solution": bool(use_pressure_solution),
            "desai": bool(use_desai),
        }
    }
    with open(os.path.join(output_folder, "pressure_schedule.json"), "w") as f:
        json.dump(pressure_data, f, indent=2)

    # Outputs
    output_folder_operation = os.path.join(output_folder, "operation")
    if MPI.COMM_WORLD.rank == 0:
        print("Operation out:", output_folder_operation)

    out_op = SparseSaveFields(mom_eq, interval=15)
    out_op.set_output_folder(output_folder_operation)
    out_op.add_output_field("u", "Displacement (m)")
    out_op.add_output_field("eps_tot", "Total strain (-)")
    out_op.add_output_field("p_elems", "Mean stress (Pa)")
    out_op.add_output_field("q_elems", "Von Mises stress (Pa)")
    out_op.add_output_field("sig", "Stress (Pa)")

    # Alleen zinvol als Desai aan staat; anders blijven ze leeg/0
    out_op.add_output_field("eps_vp", "Viscoplastic strain (-)")
    out_op.add_output_field("alpha", "Hardening parameter (-)")
    out_op.add_output_field("Fvp", "Yield function (-)")

    sim_op = sf.Simulator_M(mom_eq, tc_operation, [out_op], False)
    sim_op.run()


if __name__ == "__main__":
    main()
