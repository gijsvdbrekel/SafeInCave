import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import sys
import numpy as np
import torch as to

import safeincave as sf
import safeincave.Utils as ut
from safeincave.Utils import GPa, MPa, day, hour, create_field_elems, create_field_nodes
import safeincave.HeatBC as heatBC
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI


# ------------ Scenario helpers (from v1) ----------------
def tile_profile(times_s, pressures_pa, weeks):
    """Repeat a 1-week profile (0..7 days) for 'weeks' weeks. Returns (t, p) strictly increasing."""
    W = 7 * day
    t_out, p_out = [], []
    for w in range(weeks):
        shift = w * W
        for t, p in zip(times_s, pressures_pa):
            t_out.append(float(t + shift))
            p_out.append(float(p))
    # sort + dedup
    order = np.argsort(t_out)
    t_sorted = [t_out[i] for i in order]
    p_sorted = [p_out[i] for i in order]
    t_final, p_final = [], []
    last_t = None
    for t, p in zip(t_sorted, p_sorted):
        if last_t is None or t > last_t:
            t_final.append(t); p_final.append(p); last_t = t
    return t_final, p_final


def build_week_profile_fraction(knots_days, fractions, p_roof):
    """Single-week (0..7 d) piecewise-linear profile from day knots & fractions of p_roof."""
    assert len(knots_days) == len(fractions), "knots and fractions length mismatch"
    times_s = [float(d * day) for d in knots_days]
    pressures_pa = [float(f * p_roof) for f in fractions]
    return times_s, pressures_pa


def max_daily_swing_bar(times_s, p_pa):
    """Max |Δp| in any 24h window, in bar."""
    times = np.array(times_s, dtype=float)
    press = np.array(p_pa, dtype=float)
    max_drop = 0.0
    for t0 in times:
        mask = (times >= t0) & (times <= t0 + 24 * hour)
        window = press[mask]
        if len(window) > 1:
            drop = np.max(window) - np.min(window)
            if drop > max_drop:
                max_drop = drop
    return max_drop / 1e5  # Pa → bar


# ------------ Geometry helpers (same as your v2) ---------------
def get_geometry_parameters(path_to_grid):
    with open(os.path.join(path_to_grid, "geom.geo"), "r") as f:
        data = f.readlines()
    ovb_thickness = float(data[10][len("ovb_thickness = "):-2])
    salt_thickness = float(data[11][len("salt_thickness = "):-2])
    hanging_wall = float(data[12][len("hanging_wall = "):-2])
    return ovb_thickness, salt_thickness, hanging_wall


def main():
    # -------------------- Grid --------------------
    grid_path = os.path.join("..", "..", "..", "grids", "cavern_overburden")
    grid = sf.GridHandlerGMSH("geom", grid_path)

    # -------------------- Output root --------------------
    output_folder = os.path.join("output", "case_1_chemical_scenario")

    # -------------------- Regions --------------------
    ind_salt = grid.region_indices["Salt"]
    ind_ovb = grid.region_indices["Overburden"]

    # -------------------- Momentum equation & solver --------------------
    mom_eq = sf.LinearMomentum(grid, theta=0.0)  # fixed-point for nonelastic

    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("bcgs")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-12, max_it=100)
    mom_eq.set_solver(mom_solver)

    # -------------------- Material model --------------------
    mat = sf.Material(mom_eq.n_elems)

    gas_density = 0.082  # kg/m^3 (used only in hydrostatic add-on)
    salt_density = 2200.0
    ovb_density = 2800.0

    rho = to.zeros(mom_eq.n_elems, dtype=to.float64)
    rho[ind_salt] = salt_density
    rho[ind_ovb] = ovb_density
    mat.set_density(rho)

    # Elastic
    E0 = to.zeros(mom_eq.n_elems)
    E0[ind_salt] = 102 * GPa
    E0[ind_ovb] = 180 * GPa
    nu0 = 0.3 * to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")

    # Kelvin–Voigt
    E1 = 10 * GPa * to.ones(mom_eq.n_elems)
    nu1 = 0.32 * to.ones(mom_eq.n_elems)
    eta = to.zeros(mom_eq.n_elems)
    eta[ind_salt] = 105e11
    eta[ind_ovb] = 105e21
    kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")

    # Dislocation creep (salt)
    A_ds = to.zeros(mom_eq.n_elems)
    A_ds[ind_salt] = 1.9e-20
    Q_ds = 51600 * to.ones(mom_eq.n_elems)
    n_ds = 3.0 * to.ones(mom_eq.n_elems)
    creep_ds = sf.DislocationCreep(A_ds, Q_ds, n_ds, "ds_creep")

    # Pressure-solution creep (salt)
    A_ps = to.zeros(mom_eq.n_elems)
    A_ps[ind_salt] = 1.29e-19
    Q_ps = 13184 * to.ones(mom_eq.n_elems)
    d_ps = 0.01 * to.ones(mom_eq.n_elems)
    creep_ps = sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "ps_creep")

    # Thermoelastic
    alpha = to.zeros(mom_eq.n_elems)
    alpha[ind_salt] = 44e-6
    alpha[ind_ovb] = 0.0
    thermo = sf.Thermoelastic(alpha, "thermo")

    mat.add_to_elastic(spring_0)
    mat.add_to_thermoelastic(thermo)
    mat.add_to_non_elastic(kelvin)
    mat.add_to_non_elastic(creep_ds)
    mat.add_to_non_elastic(creep_ps)

    mom_eq.set_material(mat)

    # -------------------- Body force (downward) --------------------
    g_down = -9.81
    mom_eq.build_body_force([0.0, 0.0, g_down])

    # -------------------- Initial temperature field --------------------
    km = 1000.0
    dTdZ = 27.0 / km
    T_top = 273.0 + 20.0
    T_field_fun = lambda x, y, z: T_top + dTdZ * (660.0 - z)
    T0_field = create_field_elems(grid, T_field_fun)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # -------------------- Equilibrium (pseudo-time) --------------------
    tc_eq = sf.TimeControllerParabolic(n_time_steps=20, initial_time=0.0, final_time=10, time_unit="day")

    # Kinematic clamps
    bc_west_salt = momBC.DirichletBC("West_salt", component=0, values=[0.0, 0.0],
                                     time_values=[0.0, tc_eq.t_final])
    bc_west_ovb = momBC.DirichletBC("West_ovb", component=0, values=[0.0, 0.0],
                                    time_values=[0.0, tc_eq.t_final])
    bc_east_salt = momBC.DirichletBC("East_salt", component=0, values=[0.0, 0.0],
                                     time_values=[0.0, tc_eq.t_final])
    bc_east_ovb = momBC.DirichletBC("East_ovb", component=0, values=[0.0, 0.0],
                                    time_values=[0.0, tc_eq.t_final])
    bc_south_salt = momBC.DirichletBC("South_salt", component=1, values=[0.0, 0.0],
                                      time_values=[0.0, tc_eq.t_final])
    bc_south_ovb = momBC.DirichletBC("South_ovb", component=1, values=[0.0, 0.0],
                                     time_values=[0.0, tc_eq.t_final])
    bc_north_salt = momBC.DirichletBC("North_salt", component=1, values=[0.0, 0.0],
                                      time_values=[0.0, tc_eq.t_final])
    bc_north_ovb = momBC.DirichletBC("North_ovb", component=1, values=[0.0, 0.0],
                                     time_values=[0.0, tc_eq.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", component=2, values=[0.0, 0.0],
                                  time_values=[0.0, tc_eq.t_final])

    # Geometry → reference pressures
    z_surface = 0.0
    ovb_thickness, salt_thickness, hanging_wall = get_geometry_parameters(grid_path)
    cavern_roof = ovb_thickness + hanging_wall
    g_mag = 9.81  # positive for magnitudes
    p_roof = salt_density * g_mag * hanging_wall + ovb_density * g_mag * ovb_thickness  # Pa

    # Top traction-free + cavern constant 0.8 p_roof during equilibrium
    bc_top_eq = momBC.NeumannBC("Top", direction=2, density=0.0, ref_pos=z_surface,
                                values=[0.0, 0.0], time_values=[0.0, 10 * day], g=g_down)
    bc_cav_eq = momBC.NeumannBC("Cavern", direction=2, density=gas_density, ref_pos=cavern_roof,
                                values=[0.8 * p_roof, 0.8 * p_roof], time_values=[0.0, tc_eq.t_final], g=g_down)

    bc_equilibrium = momBC.BcHandler(mom_eq)
    for bc in [bc_west_salt, bc_west_ovb, bc_east_salt, bc_east_ovb,
               bc_south_salt, bc_south_ovb, bc_north_salt, bc_north_ovb,
               bc_bottom, bc_top_eq, bc_cav_eq]:
        bc_equilibrium.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_equilibrium)

    # Outputs eq
    out_eq = sf.SaveFields(mom_eq)
    out_eq.set_output_folder(os.path.join(output_folder, "equilibrium"))
    out_eq.add_output_field("u", "Displacement (m)")
    out_eq.add_output_field("eps_tot", "Total strain (-)")
    out_eq.add_output_field("sig", "Stress (Pa)")
    out_eq.add_output_field("p_elems", "Mean stress (Pa)")
    out_eq.add_output_field("q_elems", "Von Mises stress (Pa)")
    sim_eq = sf.Simulator_M(mom_eq, tc_eq, [out_eq], True)
    sim_eq.run()

    # -------------------- Build cavern pressure scenario (OPERATION) --------------------
    # options: "baseline", "chemical", "power", "weekly_weekday_draw", "weekly_twice"
    SCENARIO = "chemical"

    if SCENARIO == "baseline":
        t_values = [0 * day, 2 * day, 6 * day, 8 * day, 10 * day]
        p_values = [0.8 * p_roof, 0.2 * p_roof, 0.2 * p_roof, 0.8 * p_roof, 0.8 * p_roof]

    elif SCENARIO == "chemical":
        blocks = 10  # 10×48h ≈ 20 days (scenario defines its own length)
        t_block_h = [0, 12, 24, 36, 48]
        f_block = [0.90, 0.895, 0.890, 0.895, 0.90]
        raw_t, raw_p = [], []
        for b in range(blocks):
            base = b * 48 * hour
            for h, f in zip(t_block_h, f_block):
                raw_t.append(base + h * hour)
                raw_p.append(f * p_roof)
        order = np.argsort(raw_t)
        t_values, p_values = [], []
        last = None
        for i in order:
            t = float(raw_t[i]); p = float(raw_p[i])
            if last is None or t > last:
                t_values.append(t); p_values.append(p); last = t

    elif SCENARIO == "power":
        N_days = 30
        hours_markers = [0, 2, 6, 10, 14, 18, 24]
        f_day = [0.95, 0.85, 0.78, 0.82, 0.75, 0.83, 0.95]
        raw_t, raw_p = [], []
        for d in range(N_days):
            base = d * 24 * hour
            for h, f in zip(hours_markers, f_day):
                raw_t.append(base + h * hour)
                raw_p.append(f * p_roof)
        order = np.argsort(raw_t)
        t_values, p_values = [], []
        last = None
        for i in order:
            t = float(raw_t[i]); p = float(raw_p[i])
            if last is None or t > last:
                t_values.append(t); p_values.append(p); last = t

    elif SCENARIO == "weekly_weekday_draw":
        knots_days = [0, 5, 7]  # Mon00 -> Fri24 -> Mon00
        fractions = [0.90, 0.82, 0.90]
        week_t, week_p = build_week_profile_fraction(knots_days, fractions, p_roof)
        weeks = 34  # ~34 weeks ≈ 238 d
        t_values, p_values = tile_profile(week_t, week_p, weeks)

    elif SCENARIO == "weekly_twice":
        knots_days = [0, 3, 5, 7]
        fractions = [0.92, 0.84, 0.80, 0.92]
        week_t, week_p = build_week_profile_fraction(knots_days, fractions, p_roof)
        weeks = 34
        t_values, p_values = tile_profile(week_t, week_p, weeks)

    else:
        raise ValueError("Unknown SCENARIO")

    # Safety check: ≤ 10 bar in any 24h
    swing_bar = max_daily_swing_bar(t_values, p_values)
    print(f"Max Δp in any 24h: {swing_bar:.2f} bar")
    if swing_bar > 10.0:
        raise ValueError(f"Pressure swing too large: {swing_bar:.2f} bar > 10.0 bar")

    # -------------------- Decide operation duration from scenario --------------------
    if len(t_values) == 0:
        raise ValueError("Scenario produced empty pressure profile")
    max_t_s = float(max(t_values))  # seconds
    op_final_days = max_t_s / day

    # create time controller for operation using scenario duration
    tc_op = sf.TimeController(dt=0.5, initial_time=0.0, final_time=op_final_days, time_unit="day")

    # -------------------- Heat equation for operation (now using tc_op) --------------------
    heat_eq = sf.HeatDiffusion(grid)

    solver_heat = PETSc.KSP().create(grid.mesh.comm)
    solver_heat.setType("cg")
    solver_heat.getPC().setType("asm")
    solver_heat.setTolerances(rtol=1e-12, max_it=100)
    heat_eq.set_solver(solver_heat)

    cp = 850 * to.ones(heat_eq.n_elems, dtype=to.float64)
    k = 7 * to.ones(heat_eq.n_elems, dtype=to.float64)
    mat.set_specific_heat_capacity(cp)
    mat.set_thermal_conductivity(k)
    heat_eq.set_material(mat)

    T0_field_nodes = create_field_nodes(grid, T_field_fun)
    heat_eq.set_initial_T(T0_field_nodes)

    # Heat BCs
    nt_h = 2
    bcH = heatBC.BcHandler(heat_eq)
    bcH.add_boundary_condition(heatBC.DirichletBC("Top", nt_h * [T_top], [tc_op.t_initial, tc_op.t_final]))
    bcH.add_boundary_condition(heatBC.NeumannBC("Bottom", nt_h * [dTdZ], [tc_op.t_initial, tc_op.t_final]))
    for name in ["East_salt", "East_ovb", "West_salt", "West_ovb", "South_salt", "South_ovb", "North_salt",
                 "North_ovb"]:
        bcH.add_boundary_condition(heatBC.NeumannBC(name, nt_h * [0.0], [tc_op.t_initial, tc_op.t_final]))
    T_gas = T_top
    h_conv = 5.0
    bcH.add_boundary_condition(heatBC.RobinBC("Cavern", nt_h * [T_gas], h_conv, [tc_op.t_initial, tc_op.t_final]))
    heat_eq.set_boundary_conditions(bcH)

    # -------------------- Operation BCs (v2) --------------------
    # Kinematic (same as eq, but spanning full op duration)
    bc_west_salt = momBC.DirichletBC("West_salt", component=0, values=[0.0, 0.0],
                                     time_values=[0.0, tc_op.t_final])
    bc_west_ovb = momBC.DirichletBC("West_ovb", component=0, values=[0.0, 0.0],
                                    time_values=[0.0, tc_op.t_final])
    bc_east_salt = momBC.DirichletBC("East_salt", component=0, values=[0.0, 0.0],
                                     time_values=[0.0, tc_op.t_final])
    bc_east_ovb = momBC.DirichletBC("East_ovb", component=0, values=[0.0, 0.0],
                                    time_values=[0.0, tc_op.t_final])
    bc_south_salt = momBC.DirichletBC("South_salt", component=1, values=[0.0, 0.0],
                                      time_values=[0.0, tc_op.t_final])
    bc_south_ovb = momBC.DirichletBC("South_ovb", component=1, values=[0.0, 0.0],
                                     time_values=[0.0, tc_op.t_final])
    bc_north_salt = momBC.DirichletBC("North_salt", component=1, values=[0.0, 0.0],
                                      time_values=[0.0, tc_op.t_final])
    bc_north_ovb = momBC.DirichletBC("North_ovb", component=1, values=[0.0, 0.0],
                                     time_values=[0.0, tc_op.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", component=2, values=[0.0, 0.0],
                                  time_values=[0.0, tc_op.t_final])
    bc_top_op = momBC.NeumannBC("Top", direction=2, density=0.0, ref_pos=z_surface,
                                values=[0.0, 0.0], time_values=[0.0, tc_op.t_final], g=g_down)

    # Cavern pressure scenario (important: pass g_down to keep hydrostatic compression with depth)
    bc_cavern_op = momBC.NeumannBC("Cavern", direction=2, density=gas_density, ref_pos=cavern_roof,
                                   values=p_values, time_values=t_values, g=g_down)

    bc_operation = momBC.BcHandler(mom_eq)
    for bc in [bc_west_salt, bc_west_ovb, bc_east_salt, bc_east_ovb,
               bc_south_salt, bc_south_ovb, bc_north_salt, bc_north_ovb,
               bc_bottom, bc_top_op, bc_cavern_op]:
        bc_operation.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_operation)

    # -------------------- Outputs (operation) --------------------
    # Add after setting up mom_eq and heat_eq, before running simulation:

# -------------------- Outputs (operation) --------------------
    out_m = sf.SaveFields(mom_eq)
    out_m.set_output_folder(os.path.join(output_folder, "operation"))
    out_m.add_output_field("u", "Displacement (m)")
    out_m.add_output_field("p_elems", "Mean stress (Pa)")
    out_m.add_output_field("q_elems", "Von Mises stress (Pa)")

    out_T = sf.SaveFields(heat_eq)
    out_T.set_output_folder(os.path.join(output_folder, "operation"))
    out_T.add_output_field("T", "Temperature (K)")

    
    out_m.save_mesh()




    # -------------------- Run coupled TM --------------------
    sim = sf.Simulator_TM(mom_eq, heat_eq, tc_op, [out_m, out_T], compute_elastic_response=False)
    sim.run()


if __name__ == '__main__':
    main()