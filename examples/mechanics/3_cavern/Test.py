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
        except:
            pass




def build_sinus_pressure_schedule(tc, *, p_mean, p_ampl, period_hours, phase_hours=0.0,
                                  clamp_min=None, clamp_max=None):
    """Sinus schedule sampled at simulation time steps."""
    period = period_hours * ut.hour
    phase  = phase_hours  * ut.hour

    n_steps = int(math.floor(tc.t_final / tc.dt))
    t_vals = [k * tc.dt for k in range(n_steps + 1)]
    if abs(t_vals[-1] - tc.t_final) > 1e-12:
        t_vals.append(tc.t_final)

    two_pi_over_T = (2.0 * math.pi / period) if period > 0.0 else 0.0

    p_vals = []
    for t in t_vals:
        p = p_mean if period <= 0.0 else p_mean + p_ampl * math.sin(two_pi_over_T * (t - phase))
        if clamp_min is not None: p = max(p, clamp_min)
        if clamp_max is not None: p = min(p, clamp_max)
        p_vals.append(p)

    return t_vals, p_vals


def _cardinal_segment(p0, p1, p2, p3, u, tension):
    """Cardinal (Catmull–Rom with tension) for scalar values, u in [0,1]."""
    # tangents with tension τ (0 = Catmull-Rom; 1 = straight lines)
    m1 = (1 - tension) * 0.5 * (p2 - p0)
    m2 = (1 - tension) * 0.5 * (p3 - p1)
    u2 = u * u
    u3 = u2 * u
    h00 =  2*u3 - 3*u2 + 1
    h10 =      u3 - 2*u2 + u
    h01 = -2*u3 + 3*u2
    h11 =      u3 -   u2
    return h00*p1 + h10*m1 + h01*p2 + h11*m2


def _cardinal_interp(ts, ps, t, tension):
    """Evaluate cardinal spline at time t (scalar)."""
    # clamp outside range
    if t <= ts[0]: return ps[0]
    if t >= ts[-1]: return ps[-1]
    # find interval i with ts[i] <= t < ts[i+1]
    i = np.searchsorted(ts, t) - 1
    t0 = ts[i]; t1 = ts[i+1]
    # neighbors for tangents (edge-replicate)
    i_m1 = max(i-1, 0)
    i_p2 = min(i+2, len(ts)-1)
    u = (t - t0) / (t1 - t0)
    return _cardinal_segment(ps[i_m1], ps[i], ps[i+1], ps[i_p2], u, tension)


def build_irregular_pressure_schedule(tc,
                                      times_hours, pressures_MPa,
                                      *, smooth=0.3,  # 0=Catmull-Rom; ~0.2–0.5 is gentle smoothing
                                      clamp_min=None, clamp_max=None,
                                      resample_at_dt=True):
    """
    Irregular (non-periodic) schedule from user waypoints, optionally smoothed by a Cardinal spline.
    If smooth is None, use piecewise-linear interpolation.
    """
    if len(times_hours) != len(pressures_MPa) or len(times_hours) < 2:
        raise ValueError("Provide at least two waypoints with matching lengths.")

    # convert to solver units
    knots_t = np.array([h * ut.hour for h in times_hours], dtype=float)
    knots_p = np.array([p * ut.MPa  for p in pressures_MPa], dtype=float)

    # ensure coverage over [0, t_final] by padding if needed
    if knots_t[0] > 0.0:
        knots_t   = np.insert(knots_t, 0, 0.0)
        knots_p   = np.insert(knots_p, 0, knots_p[0])
    if knots_t[-1] < tc.t_final:
        knots_t   = np.append(knots_t, tc.t_final)
        knots_p   = np.append(knots_p, knots_p[-1])

    # sample times
    if resample_at_dt:
        n_steps = int(math.floor(tc.t_final / tc.dt))
        t_vals = [k * tc.dt for k in range(n_steps + 1)]
        if abs(t_vals[-1] - tc.t_final) > 1e-12:
            t_vals.append(tc.t_final)
    else:
        # just use knots (BC will interpolate between)
        t_vals = knots_t.tolist()

    # build values
    p_vals = []
    if smooth is None:
        # pure piecewise-linear
        for t in t_vals:
            p = np.interp(t, knots_t, knots_p)
            if clamp_min is not None: p = max(p, clamp_min)
            if clamp_max is not None: p = min(p, clamp_max)
            p_vals.append(p)
    else:
        # Cardinal spline with tension=smooth
        τ = float(np.clip(smooth, 0.0, 1.0))
        for t in t_vals:
            p = _cardinal_interp(knots_t, knots_p, t, τ)
            if clamp_min is not None: p = max(p, clamp_min)
            if clamp_max is not None: p = min(p, clamp_max)
            p_vals.append(p)

    return t_vals, p_vals




DAY_H = 24.0  # hours per day

def _sample_at_dt(tc, t_end=None):
    """Times aligned with solver steps, inclusive of t_final."""
    t_end = tc.t_final if t_end is None else t_end
    n_steps = int(math.floor(t_end / tc.dt))
    t_vals = [k * tc.dt for k in range(n_steps + 1)]
    if abs(t_vals[-1] - t_end) > 1e-12:
        t_vals.append(t_end)
    return t_vals

def _stretch_hours(times_h, factor_days):
    return [t * float(factor_days) for t in times_h]

def _repeat_hours(times_h, days):
    out = []
    for d in range(int(days)):
        off = d * DAY_H
        # avoid duplicate at day boundaries: skip first item except on day 0
        start = 0 if d == 0 else 1
        out.extend([off + t for t in times_h[start:]])
    return out

def build_linear_schedule_multi(tc, times_h, pressures_MPa, *, days, mode, resample_at_dt=True):
    """Linear (piecewise) schedule over multiple days."""
    if mode not in ("repeat", "stretch"):
        raise ValueError("mode must be 'repeat' or 'stretch'")
    if mode == "repeat":
        t_h = _repeat_hours(times_h, days)
    else:  # stretch
        t_h = _stretch_hours(times_h, days)
    # ensure covers [0, total_hours]
    total_h = days * DAY_H
    if t_h[0] > 0.0:
        t_h.insert(0, 0.0)
        pressures_MPa = [pressures_MPa[0]] + pressures_MPa
    if t_h[-1] < total_h:
        t_h.append(total_h)
        pressures_MPa = pressures_MPa + [pressures_MPa[-1]]
    # to seconds
    knots_t = np.array([h * ut.hour for h in t_h], dtype=float)
    knots_p = np.array([p * ut.MPa  for p in pressures_MPa], dtype=float)
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
                                   resample_at_dt=True):
    """Irregular schedule over multiple days using your cardinal-spline/piecewise helper."""
    times_h = base_waypoints_h
    if mode == "repeat":
        times_h_multi = _repeat_hours(times_h, days)
        pressures_multi = []
        for d in range(int(days)):
            start = 0 if d == 0 else 1
            pressures_multi.extend(base_pressures_MPa[start:])
    else:  # stretch
        times_h_multi = _stretch_hours(times_h, days)
        pressures_multi = list(base_pressures_MPa)

    # Use your existing build_irregular_pressure_schedule
    return build_irregular_pressure_schedule(
        tc,
        times_hours=times_h_multi,
        pressures_MPa=pressures_multi,
        smooth=smooth,
        clamp_min=clamp_min, clamp_max=clamp_max,
        resample_at_dt=resample_at_dt
    )

def build_sinus_schedule_multi(tc, *, p_mean, p_ampl, days, mode,
                               daily_period_hours=24.0, clamp_min=0.0, clamp_max=None):
    """
    Sinus that adapts to multi-day:
      - mode='repeat': period = 24 h (daily cycles)
      - mode='stretch': period = total duration (one cycle over all days)
    """
    if mode == "repeat":
        T_hours = daily_period_hours
    elif mode == "stretch":
        T_hours = days * DAY_H
    else:
        raise ValueError("mode must be 'repeat' or 'stretch'")

    return build_sinus_pressure_schedule(
        tc, p_mean=p_mean, p_ampl=p_ampl,
        period_hours=T_hours, phase_hours=0.0,
        clamp_min=clamp_min, clamp_max=clamp_max
    )




def main():
    # Read grid
    grid_path = os.path.join("..", "..", "..", "grids", "cavern_irregular")
    grid = sf.GridHandlerGMSH("geom", grid_path)

    # Define output folder
    output_folder = os.path.join("output", "case_sinus_50days")

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
    rho = salt_density*to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    # Constitutive model
    #E0 = 102*ut.GPa*to.ones(mom_eq.n_elems)
    #nu0 = 0.3*to.ones(mom_eq.n_elems)
    #spring_0 = sf.Spring(E0, nu0, "spring")

    E0 = 20.425 * ut.GPa * to.ones(mom_eq.n_elems)
    nu0 = 0.25 * to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")


    # Create Kelvin-Voigt viscoelastic element
    eta = 105e11*to.ones(mom_eq.n_elems)
    E1 = 10*ut.GPa*to.ones(mom_eq.n_elems)
    nu1 = 0.25*to.ones(mom_eq.n_elems)
    kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")

    # Create dislocation creep
    #A = 1.9e-20*to.ones(mom_eq.n_elems)
    #Q = 51600*to.ones(mom_eq.n_elems)
    #n = 3.0*to.ones(mom_eq.n_elems)
    #creep_0 = sf.DislocationCreep(A, Q, n, "creep_dislocation")

    sec_per_year = 365.25 * 24 * 3600  # 31_557_600 s

# -------------------------
    # Dislocation creep (SI)
    # Literature: A = 40 / (MPa^n * yr), n = 4.6, Q/R = 6495 K
    # Convert A to 1 / (Pa^n * s): multiply by (1e-6)^n and divide by sec_per_year
    ndc = 4.6
    A_dc = (40.0 * (1e-6)**ndc / sec_per_year) * to.ones(mom_eq.n_elems)   # 1/(Pa^n·s)
    Q_dc = (6495.0 * 8.32) * to.ones(mom_eq.n_elems)                       # J/mol
    n_dc = ndc * to.ones(mom_eq.n_elems)
    creep_0 = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation")

    # -------------------------
    # Pressure-solution creep (SI)
    # Literature: A = 14176 K·mm^3/(MPa·yr), d = 5.25 mm, Q/R = 3252 K, n = 1 (built-in)
    # Convert A to K·m^3/(Pa·s): mm^3 -> 1e-9 m^3, MPa -> 1e6 Pa, yr -> sec_per_year
    A_ps = (14176.0 * 1e-9 / 1e6 / sec_per_year) * to.ones(mom_eq.n_elems) # K·m^3/(Pa·s)
    d_ps = (5.25e-3) * to.ones(mom_eq.n_elems)                              # m
    Q_ps = (3252.0 * 8.32) * to.ones(mom_eq.n_elems)                        # J/mol
    creep_pressure = sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "creep_pressure")



    # Create constitutive model
    mat.add_to_elastic(spring_0)
    mat.add_to_non_elastic(kelvin)
    mat.add_to_non_elastic(creep_0)
    mat.add_to_non_elastic(creep_pressure)

    # Set constitutive model
    mom_eq.set_material(mat)

    # Set body forces
    g = -9.81
    g_vec = [0.0, 0.0, g]
    mom_eq.build_body_force(g_vec)

    # Set initial temperature field
    T0_field = 298*to.ones(mom_eq.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # Time settings for equilibrium stage
    tc_equilibrium = sf.TimeController(dt=0.5, initial_time=0.0, final_time=10, time_unit="hour")

    # Boundary conditions for equilibrium
    bc_west = momBC.DirichletBC(boundary_name = "West", 
                     		component = 0,
                            values = [0.0, 0.0],
                            time_values = [0.0, tc_equilibrium.t_final])

    bc_bottom = momBC.DirichletBC(boundary_name = "Bottom", 
                     	  component = 2,
                     	  values = [0.0, 0.0],
                     	  time_values = [0.0, tc_equilibrium.t_final])

    bc_south = momBC.DirichletBC(boundary_name = "South", 
                     	  component = 1,
                     	  values = [0.0, 0.0],
                     	  time_values = [0.0, tc_equilibrium.t_final])

    side_burden = 10.0*ut.MPa
    bc_east = momBC.NeumannBC(boundary_name = "East",
                        direction = 2,
                        density = salt_density,
                        ref_pos = 660.0,
                        values =      [side_burden, side_burden],
                        time_values = [0.0,            tc_equilibrium.t_final],
                        g = g_vec[2])

    bc_north = momBC.NeumannBC(boundary_name = "North",
                        direction = 2,
                        density = salt_density,
                        ref_pos = 660.0,
                        values =      [side_burden, side_burden],
                        time_values = [0.0,            tc_equilibrium.t_final],
                        g = g_vec[2])

    over_burden = 10.0*ut.MPa
    bc_top = momBC.NeumannBC(boundary_name = "Top",
                        direction = 2,
                        density = 0.0,
                        ref_pos = 0.0,
                        values =      [over_burden, over_burden],
                        time_values = [0.0,            tc_equilibrium.t_final],
                        g = g_vec[2])

    gas_density = 0.082
    p_gas = 10.0*ut.MPa
    bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
                        direction = 2,
                        density = gas_density,
                        ref_pos = 430.0,
                        values =      [p_gas, p_gas],
                        time_values = [0.0,            tc_equilibrium.t_final],
                        g = g_vec[2])

    bc_equilibrium = momBC.BcHandler(mom_eq)
    bc_equilibrium.add_boundary_condition(bc_west)
    bc_equilibrium.add_boundary_condition(bc_bottom)
    bc_equilibrium.add_boundary_condition(bc_south)
    bc_equilibrium.add_boundary_condition(bc_east)
    bc_equilibrium.add_boundary_condition(bc_north)
    bc_equilibrium.add_boundary_condition(bc_top)
    bc_equilibrium.add_boundary_condition(bc_cavern)

    # Set boundary conditions
    mom_eq.set_boundary_conditions(bc_equilibrium)

    # Equilibrium output folder
    output_folder_equilibrium = os.path.join(output_folder, "equilibrium")

    # Print output folder
    if MPI.COMM_WORLD.rank == 0:
        print(output_folder_equilibrium)

    # Create output handlers for equilibrium
    output_mom = sf.SaveFields(mom_eq)
    output_mom.set_output_folder(output_folder_equilibrium)
    output_mom.add_output_field("u", "Displacement (m)")
    output_mom.add_output_field("eps_tot", "Total strain (-)")
    output_mom.add_output_field("sig", "Stress (Pa)")
    output_mom.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
    outputs = [output_mom]

    # Define simulator for equilibrium
    sim = sf.Simulator_M(mom_eq, tc_equilibrium, outputs, True)
    sim.run()

    # ==================== OPERATION STAGE ====================

    # Create Desai's viscoplastic model
    mu_1 = 5.3665857009859815e-11*to.ones(mom_eq.n_elems)
    N_1 = 3.1*to.ones(mom_eq.n_elems)
    n = 3.0*to.ones(mom_eq.n_elems)
    a_1 = 1.965018496922832e-05*to.ones(mom_eq.n_elems)
    eta = 0.8275682807874163*to.ones(mom_eq.n_elems)
    beta_1 = 0.0048*to.ones(mom_eq.n_elems)
    beta = 0.995*to.ones(mom_eq.n_elems)
    m = -0.5*to.ones(mom_eq.n_elems)
    gamma = 0.095*to.ones(mom_eq.n_elems)
    alpha_0 = 0.0022*to.ones(mom_eq.n_elems)
    sigma_t = 1.5*to.ones(mom_eq.n_elems)
    desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0, "desai")


    stress_to = ut.numpy2torch(mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3)))
    desai.compute_initial_hardening(stress_to, Fvp_0=0.0)

    # Add Desai to the existing material and REBIND so internals are rebuilt
    mat.add_to_non_elastic(desai)
    mom_eq.set_material(mat) 

        # ==================== OPERATION STAGE ====================
    OPERATION_DAYS  = 50                  # <— set how many days you want
    SCHEDULE_MODE   = "stretch"           # "repeat" or "stretch"
    dt_hours        = 1                # time resolution

    tc_operation = sf.TimeController(dt=dt_hours, initial_time=0.0,
                                     final_time=OPERATION_DAYS*24.0,
                                     time_unit="hour")

    PRESSURE_SCENARIO = "sinus"      # "linear" | "sinus" | "irregular"

    if PRESSURE_SCENARIO == "linear":
        # base 1-day pattern (hours, MPa) — your original
        base_times_h       = [0.0, 2.0, 14.0, 16.0, 24.0]
        base_pressures_MPa = [10.0, 7.0,  7.0, 10.0, 10.0]
        t_pressure, p_pressure = build_linear_schedule_multi(
            tc_operation,
            base_times_h, base_pressures_MPa,
            days=OPERATION_DAYS, mode=SCHEDULE_MODE,
            resample_at_dt=True
        )

    elif PRESSURE_SCENARIO == "sinus":
        p_mean = 10.0 * ut.MPa
        p_ampl =  3.0 * ut.MPa
        t_pressure, p_pressure = build_sinus_schedule_multi(
            tc_operation,
            p_mean=p_mean, p_ampl=p_ampl,
            days=OPERATION_DAYS, mode=SCHEDULE_MODE,
            daily_period_hours=24.0,       # daily cycle for "repeat"
            clamp_min=0.0, clamp_max=None
        )

    elif PRESSURE_SCENARIO == "irregular":
        # Waypoints that define the *shape for one day* (hours, MPa).
        base_waypoints_h    = [0,  1.0,  2.0,  3.2,  4.0,  5.0,  6.4,  7.1,  9.0, 11.5, 13.0, 16.0, 18.0, 21.0, 24.0]
        base_pressures_MPa  = [9.0,12.0, 8.5, 11.8, 7.6, 10.2, 8.8, 11.4, 9.3, 10.7, 8.9, 11.6, 9.5, 10.2, 11.0]
        # Choose smoothing: None (piecewise linear), 0.0 (Catmull–Rom), 0.1–0.4 gentle smoothing
        t_pressure, p_pressure = build_irregular_schedule_multi(
            tc_operation,
            base_waypoints_h=base_waypoints_h,
            base_pressures_MPa=base_pressures_MPa,
            days=OPERATION_DAYS, mode=SCHEDULE_MODE,
            smooth=0.25, clamp_min=0.0, clamp_max=None,
            resample_at_dt=True
        )

    else:
        raise ValueError(f"Unknown PRESSURE_SCENARIO: {PRESSURE_SCENARIO}")



    # Boundary conditions for operation
    bc_west = momBC.DirichletBC(boundary_name = "West", 
                     		component = 0,
                            values = [0.0, 0.0],
                            time_values = [0.0, tc_operation.t_final])

    bc_bottom = momBC.DirichletBC(boundary_name = "Bottom", 
                     	  component = 2,
                     	  values = [0.0, 0.0],
                     	  time_values = [0.0, tc_operation.t_final])

    bc_south = momBC.DirichletBC(boundary_name = "South", 
                     	  component = 1,
                     	  values = [0.0, 0.0],
                     	  time_values = [0.0, tc_operation.t_final])

    bc_east = momBC.NeumannBC(boundary_name = "East",
                        direction = 2,
                        density = salt_density,
                        ref_pos = 660.0,
                        values =      [side_burden, side_burden],
                        time_values = [0.0,            tc_operation.t_final],
                        g = g_vec[2])

    bc_north = momBC.NeumannBC(boundary_name = "North",
                        direction = 2,
                        density = salt_density,
                        ref_pos = 660.0,
                        values =      [side_burden, side_burden],
                        time_values = [0.0,            tc_operation.t_final],
                        g = g_vec[2])

    bc_top = momBC.NeumannBC(boundary_name = "Top",
                        direction = 2,
                        density = 0.0,
                        ref_pos = 0.0,
                        values =      [over_burden, over_burden],
                        time_values = [0.0,            tc_operation.t_final],
                        g = g_vec[2])

    bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
                        direction = 2,
                        density = gas_density,
                        ref_pos = 430.0,
                        values = p_pressure,
                        time_values = t_pressure,
                        g = g_vec[2])

    bc_operation = momBC.BcHandler(mom_eq)
    bc_operation.add_boundary_condition(bc_west)
    bc_operation.add_boundary_condition(bc_bottom)
    bc_operation.add_boundary_condition(bc_south)
    bc_operation.add_boundary_condition(bc_east)
    bc_operation.add_boundary_condition(bc_north)
    bc_operation.add_boundary_condition(bc_top)
    bc_operation.add_boundary_condition(bc_cavern)

    # Set boundary conditions
    mom_eq.set_boundary_conditions(bc_operation)

    # Define output folder for operation
    output_folder_operation = os.path.join(output_folder, "operation")

    # Print output folder
    if MPI.COMM_WORLD.rank == 0:
        print(output_folder_operation)

    # Save pressure schedule to JSON for plotting
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

    # Create output handlers for operation
    output_mom = sf.SaveFields(mom_eq)
    output_mom.set_output_folder(output_folder_operation)
    output_mom.add_output_field("u", "Displacement (m)")
    output_mom.add_output_field("eps_tot", "Total strain (-)")
    output_mom.add_output_field("eps_vp", "Viscoplastic strain (-)")
    output_mom.add_output_field("alpha", "Hardening parameter (-)")
    output_mom.add_output_field("Fvp", "Yield function (-)")
    output_mom.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
    outputs = [output_mom]

    # Define simulator for operation
    sim = sf.Simulator_M(mom_eq, tc_operation, outputs, False)
    sim.run()


if __name__ == '__main__':
    main()