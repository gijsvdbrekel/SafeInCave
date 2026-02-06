import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC

from petsc4py import PETSc
from mpi4py import MPI
import dolfinx as do
import torch as to

import os
import json
import math


# ============================================================
#  SCENARIO SWITCHES (zet precies ÉÉN preset True)
# ============================================================
RUN_FULL                      = False  # Kelvin + disloc(CCC) + pressure-solution + Desai
RUN_MD_ONLY                   = False   # Elastic + Munson–Dawson (steady + transient)
RUN_MD_STEADY_ONLY            = False  # Elastic + Munson–Dawson (steady-only) => transient uit
RUN_FULL_MD                   = False  # Kelvin + pressure-solution + Munson–Dawson (steady + transient)
RUN_FULL_PRESSURE_SOLUTION    = True  # Elastic + Kelvin + Dislocation(CCC) + Desai  (no pressure-solution, no Munson–Dawson)

# Equilibrium stage: puur elastisch is het meest consistent voor vergelijking
EQUILIBRIUM_ELASTIC_ONLY = False


# ============================================================
#  FIXED RUN SETTINGS
# ============================================================
CAVERN_TYPE = "regular600"
OPERATION_DAYS = 365
N_CYCLES = 8
dt_hours = 2.0

PRESSURE_SCENARIO = "sinus"
SCHEDULE_MODE = "stretch"
p_gas_MPa_equilibrium = 15.0


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
    if RUN_FULL:
        return "full"
    if RUN_MD_ONLY:
        return "md_only"
    if RUN_MD_STEADY_ONLY:
        return "md_steady_only"
    if RUN_FULL_MD:
        return "full_md"
    if RUN_FULL_PRESSURE_SOLUTION:
        # name kept for folder compatibility; contents are Kelvin+Disloc+Desai (no PS)
        return "full_pressure_solution"
    return "unknown"


def _validate_scenario():
    flags = [RUN_FULL, RUN_MD_ONLY, RUN_MD_STEADY_ONLY, RUN_FULL_MD, RUN_FULL_PRESSURE_SOLUTION]
    if sum(bool(x) for x in flags) != 1:
        raise ValueError(
            "Zet precies één scenario True: "
            "RUN_FULL, RUN_MD_ONLY, RUN_MD_STEADY_ONLY, RUN_FULL_MD, RUN_FULL_PRESSURE_SOLUTION"
        )


def _equivalent_rate_from_tensor(eps_rate_3x3: to.Tensor) -> to.Tensor:
    """
    Equivalent strain-rate scalar per element.
    eps_rate_3x3: (N,3,3)
    Return: (N,) in 1/s
    """
    # Deviatoric part
    tr = eps_rate_3x3[:, 0, 0] + eps_rate_3x3[:, 1, 1] + eps_rate_3x3[:, 2, 2]
    mean = tr / 3.0
    dev = eps_rate_3x3.clone()
    dev[:, 0, 0] -= mean
    dev[:, 1, 1] -= mean
    dev[:, 2, 2] -= mean

    # J2-equivalent: sqrt(2/3 * dev:dev)
    dd = (dev * dev).sum(dim=(1, 2))
    eq = to.sqrt((2.0 / 3.0) * dd + 1e-30)
    return eq


def _to_tensor_epsrate(er, n_elems: int) -> to.Tensor | None:
    """
    Convert eps_ne_rate-like object to torch tensor (N,3,3) float64.
    Return None if shape/type not usable.
    """
    if isinstance(er, to.Tensor):
        t = er.to(dtype=to.float64)
    else:
        try:
            t = to.as_tensor(er, dtype=to.float64)
        except Exception:
            return None

    if t.ndim != 3 or t.shape[0] != n_elems or t.shape[1:] != (3, 3):
        return None
    return t


class LinearMomentumMod(sf.LinearMomentum):
    """
    Uitbreiding om extra DG0 velden weg te schrijven:
      - eps_vp (tensor) uit inelastic update (waar beschikbaar)
      - Desai state: Fvp, alpha (waar beschikbaar)
      - Equivalent strain-rate scalars (1/s):
          eps_ne_rate_eq_total
          eps_ne_rate_eq_disloc
          eps_ne_rate_eq_desai
    """
    def __init__(self, grid, theta):
        super().__init__(grid, theta)
        self.expect_vp_state = False

    def initialize(self) -> None:
        self.C.x.array[:] = to.flatten(self.mat.C)

        # Existing fields
        self.Fvp = do.fem.Function(self.DG0_1)
        self.alpha = do.fem.Function(self.DG0_1)
        self.eps_vp = do.fem.Function(self.DG0_3x3)

        # New strain-rate fields (DG0 scalar)
        self.eps_ne_rate_eq_total = do.fem.Function(self.DG0_1)
        self.eps_ne_rate_eq_disloc = do.fem.Function(self.DG0_1)
        self.eps_ne_rate_eq_desai = do.fem.Function(self.DG0_1)

        # Initialize to 0
        self.eps_ne_rate_eq_total.x.array[:] = 0.0
        self.eps_ne_rate_eq_disloc.x.array[:] = 0.0
        self.eps_ne_rate_eq_desai.x.array[:] = 0.0

    def run_after_solve(self):
        # --- eps_vp / Desai state via laatste non-elastic element ---
        elems = getattr(self.mat, "elems_ne", None)
        if elems:
            st_last = elems[-1]

            if hasattr(st_last, "eps_ne_k"):
                self.eps_vp.x.array[:] = to.flatten(st_last.eps_ne_k)

            if self.expect_vp_state:
                if hasattr(st_last, "Fvp") and hasattr(st_last, "alpha"):
                    self.Fvp.x.array[:] = st_last.Fvp
                    self.alpha.x.array[:] = st_last.alpha

        # --- strain-rate outputs ---
        n = self.n_elems

        # tensor sums (N,3,3) -> then take equivalent
        er_total = to.zeros((n, 3, 3), dtype=to.float64)
        er_disloc = to.zeros((n, 3, 3), dtype=to.float64)
        er_desai = to.zeros((n, 3, 3), dtype=to.float64)

        if elems:
            for el in elems:
                if not hasattr(el, "eps_ne_rate"):
                    continue

                er = _to_tensor_epsrate(el.eps_ne_rate, n)
                if er is None:
                    continue

                er_total += er

                # Match on element name for categorization
                name = str(getattr(el, "name", "")).lower()
                if "dislocation" in name or "munson" in name:
                    er_disloc += er
                if "desai" in name:
                    er_desai += er

        rate_total = _equivalent_rate_from_tensor(er_total)
        rate_disloc = _equivalent_rate_from_tensor(er_disloc)
        rate_desai = _equivalent_rate_from_tensor(er_desai)

        self.eps_ne_rate_eq_total.x.array[:] = rate_total.detach().cpu().numpy()
        self.eps_ne_rate_eq_disloc.x.array[:] = rate_disloc.detach().cpu().numpy()
        self.eps_ne_rate_eq_desai.x.array[:] = rate_desai.detach().cpu().numpy()


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


class KSPConvergenceLogger:
    """
    Log PETSc KSP stats per time step to a JSONL file (1 json per line).
    Call .record(t) after each solve.
    """
    def __init__(self, ksp: PETSc.KSP, filepath: str):
        self.ksp = ksp
        self.filepath = filepath
        self._fh = None
        if MPI.COMM_WORLD.rank == 0:
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
            self._fh = open(filepath, "w", buffering=1)

    def record(self, t: float):
        if MPI.COMM_WORLD.rank != 0 or self._fh is None:
            return
        try:
            its = int(self.ksp.getIterationNumber())
        except Exception:
            its = None
        try:
            rnorm = float(self.ksp.getResidualNorm())
        except Exception:
            rnorm = None
        try:
            reason = int(self.ksp.getConvergedReason())
        except Exception:
            reason = None

        rec = {"t": float(t), "ksp_its": its, "ksp_rnorm": rnorm, "ksp_reason": reason}
        self._fh.write(json.dumps(rec) + "\n")

    def close(self):
        if MPI.COMM_WORLD.rank == 0 and self._fh is not None:
            self._fh.close()
            self._fh = None


class SimulatorWithKSPLog(sf.Simulator_M):
    """
    Minimal wrapper: after each time step, log KSP stats.
    """
    def __init__(self, mom_eq, tc, outputs, is_equilibrium, ksp_logger: KSPConvergenceLogger):
        super().__init__(mom_eq, tc, outputs, is_equilibrium)
        self._ksp_logger = ksp_logger

    def run(self):
        for out in self.outputs:
            if hasattr(out, "save_fields"):
                orig = out.save_fields

                def wrapped_save_fields(t, _orig=orig, _logger=self._ksp_logger):
                    _logger.record(t)
                    return _orig(t)

                out.save_fields = wrapped_save_fields

        return super().run()


def make_solver(grid):
    """
    GMRES + ASM preconditioner (robust baseline).
    Als je later AMG wil (hypre/gamg), kunnen we dat toevoegen.
    """
    ksp = PETSc.KSP().create(grid.mesh.comm)
    ksp.setType("gmres")
    pc = ksp.getPC()
    pc.setType("asm")
    ksp.setTolerances(rtol=1e-12, max_it=200)
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

    # Dislocation creep (CCC oude params)
    ndc = 4.6
    A_dc = (40.0 * (1e-6) ** ndc / sec_per_year) * to.ones(mom_eq.n_elems)
    Q_dc = (6495.0 * 8.32) * to.ones(mom_eq.n_elems)
    n_dc = (ndc) * to.ones(mom_eq.n_elems)
    creep_disloc_old = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation_old")

    # Pressure-solution creep (CCC)
    A_ps = (14176.0 * 1e-9 / 1e6 / sec_per_year) * to.ones(mom_eq.n_elems)
    d_ps = 5.25e-3 * to.ones(mom_eq.n_elems)
    Q_ps = (3252.0 * 8.32) * to.ones(mom_eq.n_elems)
    creep_pressure = sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "creep_pressure")

    # ---- Munson–Dawson (steady + transient) ----
    # shear modulus mu = E / (2(1+nu))
    mu_md = E0 / (2.0 * (1.0 + nu0))

    # steady-state (CCC values)
    n_md = 4.99 * to.ones(mom_eq.n_elems)

    # A = 18.31 (1/MPa^n / yr)  ->  1/Pa^n / s
    A_md_val = 18.31
    A_md = (A_md_val * (1e-6) ** n_md / sec_per_year) * to.ones(mom_eq.n_elems)

    # Q/R = 6356 K  ->  Q = (Q/R)*R in J/mol
    R = 8.32
    Q_md = (6356.0 * R) * to.ones(mom_eq.n_elems)

    # transient params (CCC Zuidwending)
    K0_md     = 7e-7    * to.ones(mom_eq.n_elems)
    c_md      = 0.00902 * to.ones(mom_eq.n_elems)
    m_md      = 3.0     * to.ones(mom_eq.n_elems)
    alpha_w   = -13.2   * to.ones(mom_eq.n_elems)
    beta_w    = -7.738  * to.ones(mom_eq.n_elems)
    delta_md  = 0.58    * to.ones(mom_eq.n_elems)

    md = sf.MunsonDawsonCreep(
        A=A_md, Q=Q_md, n=n_md,
        K0=K0_md, c=c_md, m=m_md,
        alpha_w=alpha_w, beta_w=beta_w, delta=delta_md,
        mu=mu_md,
        name="munson_dawson"
    )

    # ---- Desai Viscoplastic ----
    mu_1 = 5.3665857009859815e-11 * to.ones(mom_eq.n_elems)
    N_1 = 3.1 * to.ones(mom_eq.n_elems)
    n_desai = 3.0 * to.ones(mom_eq.n_elems)
    a_1 = 1.965018496922832e-05 * to.ones(mom_eq.n_elems)
    eta_vp = 0.8275682807874163 * to.ones(mom_eq.n_elems)
    beta_1 = 0.0048 * to.ones(mom_eq.n_elems)
    beta_desai = 0.995 * to.ones(mom_eq.n_elems)
    m_desai = -0.5 * to.ones(mom_eq.n_elems)
    gamma_desai = 0.095 * to.ones(mom_eq.n_elems)
    alpha_0 = 0.0022 * to.ones(mom_eq.n_elems)
    sigma_t = 5.0 * to.ones(mom_eq.n_elems)
    desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta_vp, n_desai, beta_1, beta_desai,
                                  m_desai, gamma_desai, sigma_t, alpha_0, "desai")

    return {
        "salt_density": salt_density,
        "rho": rho,
        "spring": spring_0,
        "kelvin": kelvin,
        "disloc_old": creep_disloc_old,
        "pressure_solution": creep_pressure,  # kept available for RUN_FULL / RUN_FULL_MD
        "munson_dawson": md,
        "desai": desai,

        # ook handig om deze te bewaren voor "steady-only" variant:
        "md_K0_transient": K0_md,
    }


def make_material(mom_eq, parts, *,
                  use_kelvin: bool,
                  use_disloc_old: bool,
                  use_pressure_solution: bool,
                  use_munson_dawson: bool,
                  md_transient_on: bool,
                  use_desai: bool = False):
    mat = sf.Material(mom_eq.n_elems)
    mat.set_density(parts["rho"])
    mat.add_to_elastic(parts["spring"])

    if use_kelvin:
        mat.add_to_non_elastic(parts["kelvin"])

    if use_pressure_solution:
        mat.add_to_non_elastic(parts["pressure_solution"])

    if use_munson_dawson:
        # transient off? => zet K0 = 0 (aanname: K0 is transient amplitude)
        if not md_transient_on:
            parts["munson_dawson"].K0 = 0.0 * parts["md_K0_transient"]
        else:
            parts["munson_dawson"].K0 = 1.0 * parts["md_K0_transient"]

        mat.add_to_non_elastic(parts["munson_dawson"])

        # Munson–Dawson bevat steady-state al -> dislocation creep mag dan NIET aan
        if use_disloc_old:
            raise ValueError("Gebruik Munson–Dawson NIET samen met dislocation creep (double-count steady-state).")
    else:
        if use_disloc_old:
            mat.add_to_non_elastic(parts["disloc_old"])

    if use_desai:
        mat.add_to_non_elastic(parts["desai"])

    return mat


def main():
    _validate_scenario()

    grid_path = os.path.join("..", "..", "..", "..", "grids", "cavern_regular_600_3D")
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
        p_ref = 17.5 * ut.MPa
    elif CAVERN_TYPE.endswith("1200"):
        p_ref = 19.3 * ut.MPa
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

    if EQUILIBRIUM_ELASTIC_ONLY:
        mat_eq = make_material(
            mom_eq, parts,
            use_kelvin=False,
            use_disloc_old=False,
            use_pressure_solution=False,
            use_munson_dawson=False,
            md_transient_on=True,
            use_desai=False
        )
        mom_eq.expect_vp_state = False
    else:
        # Use the same inelastic elements as the operation stage for equilibrium
        eq_use_kelvin = False
        eq_use_disloc_old = False
        eq_use_pressure_solution = False
        eq_use_munson_dawson = False
        eq_md_transient_on = True
        eq_use_desai = False

        if RUN_FULL:
            eq_use_kelvin = True
            eq_use_disloc_old = True
            eq_use_pressure_solution = True
            eq_use_desai = False  # Desai typically not in equilibrium
        elif RUN_MD_ONLY:
            eq_use_munson_dawson = True
            eq_md_transient_on = True
        elif RUN_MD_STEADY_ONLY:
            eq_use_munson_dawson = True
            eq_md_transient_on = False
        elif RUN_FULL_MD:
            eq_use_kelvin = True
            eq_use_pressure_solution = True
            eq_use_munson_dawson = True
            eq_md_transient_on = True
        elif RUN_FULL_PRESSURE_SOLUTION:
            # Kelvin + Dislocation only (NO pressure-solution), no Desai in equilibrium
            eq_use_kelvin = True
            eq_use_disloc_old = True
            eq_use_pressure_solution = False
            eq_use_desai = False

        mat_eq = make_material(
            mom_eq, parts,
            use_kelvin=eq_use_kelvin,
            use_disloc_old=eq_use_disloc_old,
            use_pressure_solution=eq_use_pressure_solution,
            use_munson_dawson=eq_use_munson_dawson,
            md_transient_on=eq_md_transient_on,
            use_desai=eq_use_desai
        )
        mom_eq.expect_vp_state = eq_use_desai

    mom_eq.set_material(mat_eq)
    mom_eq.build_body_force(g_vec)

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

    # Pressure schedule
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

    # Scenario toggles
    use_kelvin = False
    use_disloc_old = False
    use_pressure_solution = False
    use_munson_dawson = False
    md_transient_on = True
    use_desai = False

    if RUN_FULL:
        use_kelvin = True
        use_disloc_old = True
        use_pressure_solution = True
        use_munson_dawson = False
        md_transient_on = True  # irrelevant
        use_desai = True

    elif RUN_MD_ONLY:
        use_kelvin = False
        use_disloc_old = False
        use_pressure_solution = False
        use_munson_dawson = True
        md_transient_on = True
        use_desai = False

    elif RUN_MD_STEADY_ONLY:
        use_kelvin = False
        use_disloc_old = False
        use_pressure_solution = False
        use_munson_dawson = True
        md_transient_on = False   # K0 -> 0
        use_desai = False

    elif RUN_FULL_MD:
        use_kelvin = True
        use_disloc_old = False
        use_pressure_solution = True
        use_munson_dawson = True
        md_transient_on = True
        use_desai = False

    elif RUN_FULL_PRESSURE_SOLUTION:
        # Kelvin + Dislocation + Desai (NO pressure-solution)
        use_kelvin = True
        use_disloc_old = True
        use_pressure_solution = False
        use_munson_dawson = False
        md_transient_on = True  # irrelevant
        use_desai = True

    # Build operation material
    mat_op = make_material(
        mom_eq, parts,
        use_kelvin=use_kelvin,
        use_disloc_old=use_disloc_old,
        use_pressure_solution=use_pressure_solution,
        use_munson_dawson=use_munson_dawson,
        md_transient_on=md_transient_on,
        use_desai=use_desai
    )

    mom_eq.set_material(mat_op)
    mom_eq.expect_vp_state = use_desai

    # Initialize Desai hardening from current stress state (after equilibrium)
    if use_desai:
        stress_to = ut.numpy2torch(mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3)))
        parts["desai"].compute_initial_hardening(stress_to, Fvp_0=0.0)

    # Operation BCs
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

    # Save pressure schedule once (new format)
    os.makedirs(output_folder, exist_ok=True)
    pressure_data = {
        "scenario": scenario,
        "pressure_scenario": PRESSURE_SCENARIO,
        "mode": SCHEDULE_MODE,
        "operation_days": OPERATION_DAYS,
        "n_cycles": N_CYCLES,
        "dt_hours": float(dt_hours),
        "units": {"t_raw": "s", "p_raw": "Pa", "t": "hour", "p": "MPa"},
        "t_values_s": [float(t) for t in t_pressure],
        "p_values_Pa": [float(p) for p in p_pressure],
        "t_hours": [float(t / ut.hour) for t in t_pressure],
        "p_MPa": [float(p / ut.MPa) for p in p_pressure],
        "active_elements": {
            "kelvin": bool(use_kelvin),
            "dislocation_old": bool(use_disloc_old),
            "pressure_solution": bool(use_pressure_solution),
            "munson_dawson": bool(use_munson_dawson),
            "md_transient_on": bool(md_transient_on) if use_munson_dawson else None,
            "desai": bool(use_desai),
        }
    }
    with open(os.path.join(output_folder, "pressure_schedule.json"), "w") as f:
        json.dump(pressure_data, f, indent=2)

    # Outputs
    output_folder_operation = os.path.join(output_folder, "operation")
    if MPI.COMM_WORLD.rank == 0:
        print("Operation out:", output_folder_operation)

    # KSP logger for operation stage
    ksp_log_op = KSPConvergenceLogger(
        mom_eq.solver,
        os.path.join(output_folder, "ksp_operation.jsonl")
    )

    # interval: bij lange runs niet alles wegschrijven
    n_steps = int(math.ceil(tc_operation.t_final / tc_operation.dt))
    target_frames = 300
    interval = max(1, int(round(n_steps / target_frames)))

    out_op = SparseSaveFields(mom_eq, interval=interval)
    out_op.set_output_folder(output_folder_operation)
    out_op.add_output_field("u", "Displacement (m)")
    out_op.add_output_field("eps_tot", "Total strain (-)")
    out_op.add_output_field("p_elems", "Mean stress (Pa)")
    out_op.add_output_field("q_elems", "Von Mises stress (Pa)")
    out_op.add_output_field("sig", "Stress (Pa)")

    # Desai fields (remain empty/0 if Desai is off)
    out_op.add_output_field("eps_vp", "Viscoplastic strain (-)")
    out_op.add_output_field("alpha", "Hardening parameter (-)")
    out_op.add_output_field("Fvp", "Yield function (-)")

    # Strain-rate outputs (1/s) - for plotting scripts compatibility
    out_op.add_output_field("eps_ne_rate_eq_total", "Non-elastic eq strain rate total (1/s)")
    out_op.add_output_field("eps_ne_rate_eq_disloc", "Dislocation eq strain rate (1/s)")
    out_op.add_output_field("eps_ne_rate_eq_desai", "Desai eq strain rate (1/s)")

    sim_op = SimulatorWithKSPLog(mom_eq, tc_operation, [out_op], False, ksp_log_op)
    sim_op.run()
    ksp_log_op.close()


if __name__ == "__main__":
    main()

