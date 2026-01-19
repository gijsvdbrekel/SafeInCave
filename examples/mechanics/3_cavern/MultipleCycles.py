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
import numpy as np


class LinearMomentumMod(sf.LinearMomentum):
    def __init__(self, grid, theta):
        super().__init__(grid, theta)

    def initialize(self) -> None:
        self.C.x.array[:] = to.flatten(self.mat.C)
        self.Fvp = do.fem.Function(self.DG0_1)
        self.alpha = do.fem.Function(self.DG0_1)
        self.eps_vp = do.fem.Function(self.DG0_3x3)

    def run_after_solve(self):
        try:
            self.eps_vp.x.array[:] = to.flatten(self.mat.elems_ne[-1].eps_ne_k)
            self.Fvp.x.array[:] = self.mat.elems_ne[-1].Fvp
            self.alpha.x.array[:] = self.mat.elems_ne[-1].alpha
        except Exception:
            pass


class SparseSaveFields(sf.SaveFields):
    """
    SaveFields that only writes every `interval`-th call after t=0.
    t = 0 is always saved.
    """
    def __init__(self, mom_eq, interval: int):
        super().__init__(mom_eq)
        self.interval = max(1, int(interval))
        self._counter = 0

    def save_fields(self, t):
        # Always save the initial state at t=0
        if t == 0:
            return super().save_fields(t)

        # Count calls and only forward every `interval`-th one
        self._counter += 1
        if self._counter % self.interval == 0:
            return super().save_fields(t)
        # otherwise: skip this step (no write)


def build_sinus_pressure_schedule(tc, *, p_mean, p_ampl, period_hours, phase_hours=0.0,
                                  clamp_min=None, clamp_max=None):
    """Sinus schedule sampled at simulation time steps."""
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


def _cardinal_segment(p0, p1, p2, p3, u, tension):
    """Cardinal (Catmull–Rom with tension) for scalar values, u in [0,1]."""
    m1 = (1 - tension) * 0.5 * (p2 - p0)
    m2 = (1 - tension) * 0.5 * (p3 - p1)
    u2 = u * u
    u3 = u2 * u
    h00 = 2*u3 - 3*u2 + 1
    h10 = u3 - 2*u2 + u
    h01 = -2*u3 + 3*u2
    h11 = u3 - u2
    return h00*p1 + h10*m1 + h01*p2 + h11*m2


def _cardinal_interp(ts, ps, t, tension):
    """Evaluate cardinal spline at time t (scalar)."""
    if t <= ts[0]:
        return ps[0]
    if t >= ts[-1]:
        return ps[-1]

    i = np.searchsorted(ts, t) - 1
    t0 = ts[i]
    t1 = ts[i+1]

    i_m1 = max(i-1, 0)
    i_p2 = min(i+2, len(ts)-1)
    u = (t - t0) / (t1 - t0)
    return _cardinal_segment(ps[i_m1], ps[i], ps[i+1], ps[i_p2], u, tension)


def build_irregular_pressure_schedule(tc,
                                      times_hours, pressures_MPa,
                                      *, smooth=0.3,
                                      clamp_min=None, clamp_max=None,
                                      resample_at_dt=True):
    if len(times_hours) != len(pressures_MPa) or len(times_hours) < 2:
        raise ValueError("Provide at least two waypoints with matching lengths.")

    knots_t = np.array([h * ut.hour for h in times_hours], dtype=float)
    knots_p = np.array([p * ut.MPa for p in pressures_MPa], dtype=float)

    if knots_t[0] > 0.0:
        knots_t = np.insert(knots_t, 0, 0.0)
        knots_p = np.insert(knots_p, 0, knots_p[0])
    if knots_t[-1] < tc.t_final:
        knots_t = np.append(knots_t, tc.t_final)
        knots_p = np.append(knots_p, knots_p[-1])

    if resample_at_dt:
        n_steps = int(math.floor(tc.t_final / tc.dt))
        t_vals = [k * tc.dt for k in range(n_steps + 1)]
        if abs(t_vals[-1] - tc.t_final) > 1e-12:
            t_vals.append(tc.t_final)
    else:
        t_vals = knots_t.tolist()

    p_vals = []
    if smooth is None:
        for t in t_vals:
            p = np.interp(t, knots_t, knots_p)
            if clamp_min is not None:
                p = max(p, clamp_min)
            if clamp_max is not None:
                p = min(p, clamp_max)
            p_vals.append(p)
    else:
        tau = float(np.clip(smooth, 0.0, 1.0))
        for t in t_vals:
            p = _cardinal_interp(knots_t, knots_p, t, tau)
            if clamp_min is not None:
                p = max(p, clamp_min)
            if clamp_max is not None:
                p = min(p, clamp_max)
            p_vals.append(p)

    return t_vals, p_vals


DAY_H = 24.0  # hours per day


def _sample_at_dt(tc, t_end=None):
    t_end = tc.t_final if t_end is None else t_end
    n_steps = int(math.floor(t_end / tc.dt))
    t_vals = [k * tc.dt for k in range(n_steps + 1)]
    if abs(t_vals[-1] - t_end) > 1e-12:
        t_vals.append(t_end)
    return t_vals


def _stretch_hours(times_h, factor_days):
    return [t * float(factor_days) for t in times_h]


def _repeat_hours(times_h, days):
    times_h = list(map(float, times_h))
    out = []
    for d in range(int(days)):
        off = d * DAY_H
        for i, t in enumerate(times_h):
            if d > 0 and i == 0 and abs(t - 0.0) < 1e-12 and abs(times_h[-1] - DAY_H) < 1e-12:
                # skip the 0h knot because previous day already had 24h knot
                continue
            out.append(off + t)
    return out



def build_linear_schedule_multi(tc, times_h, pressures_MPa, *, days, mode,
                                resample_at_dt=True, total_cycles=None):
    """
    Extend your existing linear multi-day schedule builder.
    - mode="repeat": repeat daily pattern each day (total_cycles ignored)
    - mode="stretch":
        * if total_cycles is None -> 1 cycle over total duration (old behavior)
        * if total_cycles >= 1 -> repeat the base cycle total_cycles times over total duration
    """
    if mode not in ("repeat", "stretch"):
        raise ValueError("mode must be 'repeat' or 'stretch'")

    times_h = list(map(float, times_h))
    pressures_MPa = list(map(float, pressures_MPa))
    if len(times_h) != len(pressures_MPa) or len(times_h) < 2:
        raise ValueError("Provide at least two waypoints with matching lengths.")

    total_h = float(days) * DAY_H

    if mode == "repeat":
        t_h = _repeat_hours(times_h, days)
        p_h = []
        for d in range(int(days)):
            start = 0 if d == 0 else 1
            p_h.extend(pressures_MPa[start:])

    else:  # mode == "stretch"
        if total_cycles is None:
            total_cycles = 1
        total_cycles = max(1, int(total_cycles))

        base_start = times_h[0]
        base_end = times_h[-1]
        base_duration = base_end - base_start
        if base_duration <= 0.0:
            raise ValueError("times_h must span a positive duration.")

        cycle_duration = total_h / float(total_cycles)
        scale = cycle_duration / base_duration

        t_h = []
        p_h = []
        for k in range(total_cycles):
            offset = k * cycle_duration
            for i, t in enumerate(times_h):
                if k > 0 and i == 0:
                    continue  # avoid duplicate knot at cycle boundaries
                t_scaled = offset + (t - base_start) * scale
                t_h.append(t_scaled)
                p_h.append(pressures_MPa[i])

    # Ensure coverage [0, total_h]
    if t_h[0] > 0.0:
        t_h.insert(0, 0.0)
        p_h.insert(0, p_h[0])
    if t_h[-1] < total_h:
        t_h.append(total_h)
        p_h.append(p_h[-1])

    knots_t = np.array([h * ut.hour for h in t_h], dtype=float)
    knots_p = np.array([p * ut.MPa for p in p_h], dtype=float)

    if resample_at_dt:
        t_vals = _sample_at_dt(tc)
        p_vals = np.interp(t_vals, knots_t, knots_p).tolist()
    else:
        t_vals = knots_t.tolist()
        p_vals = knots_p.tolist()

    return t_vals, p_vals


def build_irregular_schedule_multi(tc, *, base_waypoints_h, base_pressures_MPa,
                                   days, mode, smooth=0.25,
                                   clamp_min=0.0, clamp_max=None,
                                   resample_at_dt=True,
                                   total_cycles=None):
    """
    Irregular pressure schedule over meerdere dagen.

    - mode="repeat":
        herhaalt het 0–24h patroon elke dag (zoals nu al).
        'total_cycles' wordt genegeerd.

    - mode="stretch" & total_cycles is None:
        base_waypoints_h wordt opgeschaald zodat je
        precies 1 cycle over de totale simulatie (days * 24h) krijgt
        (huidige gedrag).

    - mode="stretch" & total_cycles >= 1:
        je krijgt 'total_cycles' cycles over de totale simulatie
        (days * 24h). De basisvorm blijft hetzelfde, maar wordt in
        de tijd geschaald en achter elkaar geplakt.
    """
    times_h = np.asarray(base_waypoints_h, dtype=float)
    pressures = np.asarray(base_pressures_MPa, dtype=float)

    if len(times_h) != len(pressures):
        raise ValueError("base_waypoints_h and base_pressures_MPa must have same length")

    total_hours = days * DAY_H

    if mode == "repeat":
        # Huidig gedrag: gewoon per dag herhalen
        times_h_multi = _repeat_hours(times_h, days)
        pressures_multi = []
        for d in range(int(days)):
            start = 0 if d == 0 else 1
            pressures_multi.extend(pressures[start:])

    elif mode == "stretch":
        # total_cycles: hoe vaak de basis-cyclus over de hele simulatie herhaald wordt
        if total_cycles is None:
            total_cycles = 1
        total_cycles = max(1, int(total_cycles))

        base_start = times_h[0]
        base_end = times_h[-1]
        base_duration = base_end - base_start
        if base_duration <= 0.0:
            raise ValueError("base_waypoints_h must span a positive duration.")

        # duur van één cycle zodat total_cycles * cycle_duration = total_hours
        cycle_duration = total_hours / float(total_cycles)
        scale = cycle_duration / base_duration

        times_h_multi = []
        pressures_multi = []

        for k in range(total_cycles):
            offset = k * cycle_duration
            for i, t in enumerate(times_h):
                # voorkom dubbele tijdstap aan de grenzen
                if k > 0 and i == 0:
                    continue
                t_scaled = offset + (t - base_start) * scale
                times_h_multi.append(t_scaled)
                pressures_multi.append(pressures[i])

    else:
        raise ValueError("mode must be 'repeat' or 'stretch'")

    # Zorg dat we precies de totale simulatie afdekken
    # (clip/extend eventueel een klein beetje)
    times_h_multi = np.asarray(times_h_multi, dtype=float)
    pressures_multi = np.asarray(pressures_multi, dtype=float)

    # Bouw uiteindelijk de schedule met jouw bestaande helper
    return build_irregular_pressure_schedule(
        tc,
        times_hours=times_h_multi.tolist(),
        pressures_MPa=pressures_multi.tolist(),
        smooth=smooth,
        clamp_min=clamp_min, clamp_max=clamp_max,
        resample_at_dt=resample_at_dt
    )



def build_sinus_schedule_multi(tc, *, p_mean, p_ampl, days, mode,
                               daily_period_hours=24.0,
                               total_cycles=1,
                               clamp_min=0.0, clamp_max=None):
    """
    Sinusdrukschema over meerdere dagen.

    mode = "repeat":
        - Zelfde betekenis als nu: periode = daily_period_hours
        - total_cycles wordt genegeerd (backwards compatible)

    mode = "stretch":
        - Er komen `total_cycles` volledige sinussen over de totale periode
          van `days * 24` uur.
        - Voor total_cycles=1 krijg je exact je oude gedrag.
    """
    total_hours = days * DAY_H

    if mode == "repeat":
        # Zelfde als voorheen, je bepaalt zelf de periode per cycle
        T_hours = daily_period_hours

    elif mode == "stretch":
        # Verdeel de totaalduur over een gegeven aantal cycli
        total_cycles = max(1, int(total_cycles))
        T_hours = total_hours / float(total_cycles)

    else:
        raise ValueError("mode must be 'repeat' or 'stretch'")

    return build_sinus_pressure_schedule(
        tc, p_mean=p_mean, p_ampl=p_ampl,
        period_hours=T_hours, phase_hours=0.0,
        clamp_min=clamp_min, clamp_max=clamp_max
    )



def main():
    # Read grid
    grid_path = os.path.join("..", "..", "..", "grids", "cavern_irregular_600_3D")
    grid = sf.GridHandlerGMSH("geom", grid_path)


     # --- Cavern-specific z_max ---
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

    CAVERN_TYPE = "irregular600"
    z_max = Z_MAX_BY_CAVERN[CAVERN_TYPE]

    OPERATION_DAYS = 100
    SCHEDULE_MODE = "stretch"
    N_CYCLES = 1
    dt_hours = 0.1

    PRESSURE_SCENARIO = "linear"

    
    output_folder = os.path.join(
        "output",
        f"case_{PRESSURE_SCENARIO}({N_CYCLES})_{OPERATION_DAYS}days_{CAVERN_TYPE}"
    )



    # --- Overburden & side burden depend on cavern set (600 vs 1200) ---
    if CAVERN_TYPE.endswith("600"):
        p_ref = 18.2 * ut.MPa
    elif CAVERN_TYPE.endswith("1200"):
        p_ref = 20.1 * ut.MPa
    else:
        raise ValueError(f"Cannot infer cavern set (600/1200) from CAVERN_TYPE='{CAVERN_TYPE}'")

    side_burden = p_ref
    over_burden = p_ref

    # Define momentum equation
    mom_eq = LinearMomentumMod(grid, theta=0.5)

    # Define solver
    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("cg")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-12, max_it=100)
    mom_eq.set_solver(mom_solver)

    # Define material properties
    mat = sf.Material(mom_eq.n_elems)

    # Set material density
    salt_density = 2200
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    # Elastic
    E0 = 20.425 * ut.GPa * to.ones(mom_eq.n_elems)
    nu0 = 0.25 * to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")

    # Kelvin-Voigt
    eta = 105e11 * to.ones(mom_eq.n_elems)
    E1 = 10 * ut.GPa * to.ones(mom_eq.n_elems)
    nu1 = 0.25 * to.ones(mom_eq.n_elems)
    kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")

    sec_per_year = 365.25 * 24 * 3600

    # Dislocation creep
    ndc = 4.6
    A_dc = (40.0 * (1e-6)**ndc / sec_per_year) * to.ones(mom_eq.n_elems)
    Q_dc = (6495.0 * 8.32) * to.ones(mom_eq.n_elems)
    n_dc = ndc * to.ones(mom_eq.n_elems)
    creep_0 = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation")

    # Pressure-solution creep
    A_ps = (14176.0 * 1e-9 / 1e6 / sec_per_year) * to.ones(mom_eq.n_elems)
    d_ps = 5.25e-3 * to.ones(mom_eq.n_elems)
    Q_ps = (3252.0 * 8.32) * to.ones(mom_eq.n_elems)
    creep_pressure = sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "creep_pressure")

    # Constitutive model
    mat.add_to_elastic(spring_0)
    mat.add_to_non_elastic(kelvin)
    mat.add_to_non_elastic(creep_0)
    mat.add_to_non_elastic(creep_pressure)

    mom_eq.set_material(mat)

    # Body forces
    g = -9.81
    g_vec = [0.0, 0.0, g]
    mom_eq.build_body_force(g_vec)

    # Initial temperature
    T0_field = 298 * to.ones(mom_eq.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # ===== EQUILIBRIUM STAGE =====
    tc_equilibrium = sf.TimeController(dt=0.5, initial_time=0.0, final_time=10, time_unit="hour")

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

    gas_density = 7.6  # Hydrogen density in kg/m3 at specific T and P conditions
    p_gas = 14.7 * ut.MPa
    bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, z_max,
                                [p_gas, p_gas],
                                [0.0, tc_equilibrium.t_final],
                                g=g_vec[2])

    bc_equilibrium = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
        bc_equilibrium.add_boundary_condition(bc)

    mom_eq.set_boundary_conditions(bc_equilibrium)

    output_folder_equilibrium = os.path.join(output_folder, "equilibrium")
    if MPI.COMM_WORLD.rank == 0:
        print(output_folder_equilibrium)

    output_mom = sf.SaveFields(mom_eq)
    output_mom.set_output_folder(output_folder_equilibrium)
    output_mom.add_output_field("u", "Displacement (m)")
    output_mom.add_output_field("eps_tot", "Total strain (-)")
    output_mom.add_output_field("sig", "Stress (Pa)")
    output_mom.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
    outputs = [output_mom]

    sim = sf.Simulator_M(mom_eq, tc_equilibrium, outputs, True)
    sim.run()

    # ===== OPERATION STAGE =====
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
    desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta_vp, n, beta_1, beta, m, gamma, sigma_t, alpha_0, "desai")

    stress_to = ut.numpy2torch(mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3)))
    desai.compute_initial_hardening(stress_to, Fvp_0=0.0)

    mat.add_to_non_elastic(desai)
    mom_eq.set_material(mat)


    tc_operation = sf.TimeController(dt=dt_hours, initial_time=0.0,
                                     final_time=OPERATION_DAYS*24.0,
                                     time_unit="hour")


    if PRESSURE_SCENARIO == "linear":
        p_min = 8.0
        p_max = 21.4 

        base_times_h = [0.0, 2.0, 14.0, 16.0, 24.0]
        base_pressures_MPa = [p_max, p_min, p_min, p_max, p_max]

        t_pressure, p_pressure = build_linear_schedule_multi(
            tc_operation,
            base_times_h, base_pressures_MPa,
            days=OPERATION_DAYS,
            mode=SCHEDULE_MODE,
            resample_at_dt=True,
            total_cycles=N_CYCLES,   # <-- now it actually does something
        )

    elif PRESSURE_SCENARIO == "sinus":
        p_mean = 14.7 * ut.MPa
        p_ampl = 6.7 * ut.MPa
        t_pressure, p_pressure = build_sinus_schedule_multi(
            tc_operation,
            p_mean=p_mean, p_ampl=p_ampl,
            days=OPERATION_DAYS, mode=SCHEDULE_MODE,
            daily_period_hours=24.0,   # wordt genegeerd bij "stretch"
            total_cycles=N_CYCLES,     # <— nieuw
            clamp_min=0.0, clamp_max=None
        )

    elif PRESSURE_SCENARIO == "irregular":
        base_waypoints_h = [0, 1.0, 2.0, 3.2, 4.0, 5.0, 6.4, 7.1, 9.0, 11.5,
                            13.0, 16.0, 18.0, 21.0, 24.0]
        base_pressures_MPa = [10.0, 12.0, 8.5, 11.8, 7.6, 10.2, 8.8, 11.4,
                              9.3, 10.7, 8.9, 11.6, 9.5, 10.2, 11.0]
        t_pressure, p_pressure = build_irregular_schedule_multi(
            tc_operation,
            base_waypoints_h=base_waypoints_h,
            base_pressures_MPa=base_pressures_MPa,
            days=OPERATION_DAYS, mode=SCHEDULE_MODE,
            smooth=0.25, clamp_min=0.0, clamp_max=None,
            resample_at_dt=True,
            total_cycles=N_CYCLES,   # <— voeg dit toe als je ook hier 2 cycles wilt
        )


    else:
        raise ValueError(f"Unknown PRESSURE_SCENARIO: {PRESSURE_SCENARIO}")
    

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

    output_folder_operation = os.path.join(output_folder, "operation")
    if MPI.COMM_WORLD.rank == 0:
        print(output_folder_operation)

    pressure_data = {
        "scenario": PRESSURE_SCENARIO,
        "mode": SCHEDULE_MODE,
        "operation_days": OPERATION_DAYS,
        "t_values": [float(t) for t in t_pressure],
        "p_values": [float(p) for p in p_pressure],
    }

    os.makedirs(output_folder, exist_ok=True)
    with open(os.path.join(output_folder, "pressure_schedule.json"), 'w') as f:
        json.dump(pressure_data, f, indent=2)

    # Operation outputs – use sparse saver (every 15th step)
    output_mom_op = SparseSaveFields(mom_eq, interval=15)
    output_mom_op.set_output_folder(output_folder_operation)
    output_mom_op.add_output_field("u", "Displacement (m)")
    output_mom_op.add_output_field("eps_tot", "Total strain (-)")
    output_mom_op.add_output_field("eps_vp", "Viscoplastic strain (-)")
    output_mom_op.add_output_field("alpha", "Hardening parameter (-)")
    output_mom_op.add_output_field("Fvp", "Yield function (-)")
    output_mom_op.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom_op.add_output_field("p_nodes", "Mean stress (Pa) (nodes)")
    output_mom_op.add_output_field("q_elems", "Von Mises stress (Pa)")
    output_mom_op.add_output_field("sig", "Stress (Pa)")
    outputs_op = [output_mom_op]

    sim_op = sf.Simulator_M(mom_eq, tc_operation, outputs_op, False)
    sim_op.run()


if __name__ == '__main__':
    main()
