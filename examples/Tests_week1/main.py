import os
from pathlib import Path
import json
import numpy as np
import torch as to
from mpi4py import MPI
from petsc4py import PETSc
import dolfinx as do

import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC


# ---------- Units ----------
HOUR = ut.hour          # seconds
DAY  = 24 * HOUR
MPa  = ut.MPa


# ---------- Helpers for weekly / cyclic profiles ----------
def tile_profile(times_s, pressures_pa, weeks: int):
    """Repeat a base profile (start at t=0) for 'weeks' times a 7-day window."""
    W = 7 * DAY
    t_out, p_out = [], []
    for w in range(weeks):
        shift = w * W
        for t, p in zip(times_s, pressures_pa):
            t_out.append(float(t + shift))
            p_out.append(float(p))
    order = np.argsort(t_out)
    t_sorted = [t_out[i] for i in order]
    p_sorted = [p_out[i] for i in order]
    t_final, p_final, last = [], [], None
    for t, p in zip(t_sorted, p_sorted):
        if last is None or t > last:
            t_final.append(t)
            p_final.append(p)
            last = t
    return t_final, p_final


def build_week_profile_fraction(knots_days, fractions, p_ref_pa):
    """
    One-week (0..7 days) piecewise-linear profile from day-knot fractions of p_ref.
    Returns (times_s, pressures_pa) for a single week.
    """
    assert len(knots_days) == len(fractions), "knots and fractions length mismatch"
    times_s = [float(d * DAY) for d in knots_days]
    pressures_pa = [float(f * p_ref_pa) for f in fractions]
    return times_s, pressures_pa


def safety_check_24h(time_list_s, p_list_pa, limit_bar=10.0):
    """Abort if any 24h window shows Δp > limit_bar."""
    times = np.array(time_list_s, dtype=float)
    press = np.array(p_list_pa, dtype=float)
    max_drop = 0.0
    for t0 in times:
        mask = (times >= t0) & (times <= t0 + 24 * HOUR)
        window = press[mask]
        if window.size > 1:
            drop = float(np.max(window) - np.min(window))
            max_drop = max(max_drop, drop)
    if max_drop > limit_bar * 1e5:
        raise ValueError(f"Pressure swing too large: {max_drop/1e5:.2f} bar > {limit_bar:.1f} bar")
    return max_drop / 1e5  # bar


def find_root_with_grids(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "grids").is_dir():
            return p
    raise FileNotFoundError("Could not find a 'grids' directory upwards from: " + str(start))


class LinearMomentumMod(sf.LinearMomentum):
    def __init__(self, grid, theta):
        super().__init__(grid, theta)

    def initialize(self) -> None:
        self.C.x.array[:] = to.flatten(self.mat.C)
        self.Fvp = do.fem.Function(self.DG0_1)
        self.alpha = do.fem.Function(self.DG0_1)
        self.eps_vp = do.fem.Function(self.DG0_3x3)

    def run_after_solve(self):
        # export Desai variables if present
        desai = next((e for e in self.mat.elems_ne if getattr(e, "name", "") == "desai"), None)
        if desai is not None:
            self.eps_vp.x.array[:] = to.flatten(desai.eps_ne_k)
            self.Fvp.x.array[:] = desai.Fvp
            self.alpha.x.array[:] = desai.alpha


def main():
    # ---------- Paths ----------
    HERE = Path(__file__).resolve().parent
    ROOT = find_root_with_grids(HERE)
    grid_path = (ROOT / "grids" / "cavern_irregular").resolve()
    output_folder = (HERE / "output" / "case_0").resolve()
    os.makedirs(output_folder, exist_ok=True)

    grid = sf.GridHandlerGMSH("geom", str(grid_path))

    # ---------- Momentum equation ----------
    mom_eq = LinearMomentumMod(grid, theta=0.5)

    # ---------- Solver ----------
    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("cg")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-12, max_it=100)
    mom_eq.set_solver(mom_solver)

    # ---------- Material ----------
    dtype = to.float64
    mat = sf.Material(mom_eq.n_elems)

    salt_density = 2000.0
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=dtype)
    mat.set_density(rho)

    E0 = (102 * ut.GPa) * to.ones(mom_eq.n_elems, dtype=dtype)
    nu0 = 0.3 * to.ones(mom_eq.n_elems, dtype=dtype)
    spring_0 = sf.Spring(E0, nu0, "spring")

    eta = (105e11) * to.ones(mom_eq.n_elems, dtype=dtype)
    E1 = (10 * ut.GPa) * to.ones(mom_eq.n_elems, dtype=dtype)
    nu1 = 0.32 * to.ones(mom_eq.n_elems, dtype=dtype)
    kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")

    A = (1.9e-20) * to.ones(mom_eq.n_elems, dtype=dtype)
    Q = 51600.0 * to.ones(mom_eq.n_elems, dtype=dtype)
    n = 3.0 * to.ones(mom_eq.n_elems, dtype=dtype)
    creep_0 = sf.DislocationCreep(A, Q, n, "creep")

    mat.add_to_elastic(spring_0)
    mat.add_to_non_elastic(kelvin)
    mat.add_to_non_elastic(creep_0)

    mom_eq.set_material(mat)

    # ---------- Body force ----------
    g = -9.81
    g_vec = [0.0, 0.0, g]
    mom_eq.build_body_force(g_vec)

    # ---------- Initial T ----------
    T0_field = 298.0 * to.ones(mom_eq.n_elems, dtype=dtype)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # ===================== EQUILIBRIUM STAGE =====================
    tc_equilibrium = sf.TimeController(dt=0.5, initial_time=0.0, final_time=10, time_unit="hour")

    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_equilibrium.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_equilibrium.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_equilibrium.t_final])

    side_burden = 10.0 * MPa
    over_burden = 10.0 * MPa
    gas_density = 0.082
    p_gas_eq = 10.0 * MPa

    bc_east = momBC.NeumannBC("East", 2, density=salt_density, ref_pos=660.0,
                              values=[side_burden, side_burden],
                              time_values=[0.0, tc_equilibrium.t_final], g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, density=salt_density, ref_pos=660.0,
                               values=[side_burden, side_burden],
                               time_values=[0.0, tc_equilibrium.t_final], g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, density=0.0, ref_pos=0.0,
                             values=[over_burden, over_burden],
                             time_values=[0.0, tc_equilibrium.t_final], g=g_vec[2])
    bc_cavern_eq = momBC.NeumannBC("Cavern", 2, density=gas_density, ref_pos=430.0,
                                   values=[p_gas_eq, p_gas_eq],
                                   time_values=[0.0, tc_equilibrium.t_final], g=g_vec[2])

    bc_equilibrium = momBC.BcHandler(mom_eq)
    for bc in (bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern_eq):
        bc_equilibrium.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_equilibrium)

    output_folder_equil = os.path.join(output_folder, "equilibrium")
    if MPI.COMM_WORLD.rank == 0:
        print(output_folder_equil)

    out_eq = sf.SaveFields(mom_eq)
    out_eq.set_output_folder(output_folder_equil)
    out_eq.add_output_field("u", "Displacement (m)")
    out_eq.add_output_field("eps_tot", "Total strain (-)")
    out_eq.add_output_field("sig", "Stress (Pa)")
    out_eq.add_output_field("p_elems", "Mean stress (Pa)")
    out_eq.add_output_field("q_elems", "Von Mises stress (Pa)")

    sf.Simulator_M(mom_eq, tc_equilibrium, [out_eq], True).run()

    # ===================== OPERATION STAGE (CYCLIC SCENARIOS) =====================

    p_ref = 10.0 * MPa
    scenario = "chemical"  # "baseline" | "chemical" | "power" | "weekly_weekday_draw" | "weekly_twice"

    if scenario == "baseline":
        time_list = [0*DAY, 2*DAY, 6*DAY, 8*DAY, 10*DAY]
        p_vals    = [0.8*p_ref, 0.2*p_ref, 0.2*p_ref, 0.8*p_ref, 0.8*p_ref]

    elif scenario == "chemical":
        blocks = 2
        t_block = [0, 12, 24, 36, 48]  # hours
        f_block = [0.90, 0.895, 0.890, 0.895, 0.90]
        raw_t, raw_p = [], []
        for b in range(blocks):
            base = b * 48 * HOUR
            for h, f in zip(t_block, f_block):
                raw_t.append(base + h * HOUR)
                raw_p.append(f * p_ref)
        order = np.argsort(raw_t)
        time_list, p_vals, last = [], [], None
        for i in order:
            t = float(raw_t[i]); p = float(raw_p[i])
            if last is None or t > last:
                time_list.append(t); p_vals.append(p); last = t

    elif scenario == "power":
        N_days = 3
        hours_marker = [0, 2, 6, 10, 14, 18, 24]
        f_day = [0.95, 0.85, 0.78, 0.82, 0.75, 0.83, 0.95]
        raw_t, raw_p = [], []
        for d in range(N_days):
            base = d * 24 * HOUR
            for h, f in zip(hours_marker, f_day):
                raw_t.append(base + h * HOUR)
                raw_p.append(f * p_ref)
        order = np.argsort(raw_t)
        time_list, p_vals, last = [], [], None
        for i in order:
            t = float(raw_t[i]); p = float(raw_p[i])
            if last is None or t > last:
                time_list.append(t); p_vals.append(p); last = t

    elif scenario == "weekly_weekday_draw":
        knots_days = [0, 5, 7]
        fractions  = [0.90, 0.82, 0.90]
        week_t, week_p = build_week_profile_fraction(knots_days, fractions, p_ref)
        weeks = 6
        time_list, p_vals = tile_profile(week_t, week_p, weeks)

    elif scenario == "weekly_twice":
        knots_days = [0, 3, 5, 7]
        fractions  = [0.92, 0.84, 0.80, 0.92]
        week_t, week_p = build_week_profile_fraction(knots_days, fractions, p_ref)
        weeks = 4
        time_list, p_vals = tile_profile(week_t, week_p, weeks)

    else:
        raise ValueError(f"Unknown scenario: {scenario}")

    max_drop_bar = safety_check_24h(time_list, p_vals, limit_bar=10.0)
    if MPI.COMM_WORLD.rank == 0:
        print(f"Max Δp in any 24h (bar): {max_drop_bar:.2f}")

    tc_operation = sf.TimeController(
        dt=0.1,
        initial_time=0.0,
        final_time=float(time_list[-1] / HOUR),
        time_unit="hour",
    )

    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_operation.t_final])

    bc_east = momBC.NeumannBC("East", 2, density=salt_density, ref_pos=660.0,
                              values=[side_burden, side_burden],
                              time_values=[0.0, tc_operation.t_final], g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, density=salt_density, ref_pos=660.0,
                               values=[side_burden, side_burden],
                               time_values=[0.0, tc_operation.t_final], g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, density=0.0, ref_pos=0.0,
                             values=[over_burden, over_burden],
                             time_values=[0.0, tc_operation.t_final], g=g_vec[2])

    bc_cavern = momBC.NeumannBC("Cavern", 2, density=gas_density, ref_pos=430.0,
                                values=list(map(float, p_vals)),
                                time_values=list(map(float, time_list)),
                                g=g_vec[2])

    bc_operation = momBC.BcHandler(mom_eq)
    for bc in (bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern):
        bc_operation.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_operation)

    output_folder_operation = os.path.join(output_folder, "operation")
    if MPI.COMM_WORLD.rank == 0:
        print(output_folder_operation)

    out_op = sf.SaveFields(mom_eq)
    out_op.set_output_folder(output_folder_operation)
    out_op.add_output_field("u", "Displacement (m)")
    out_op.add_output_field("eps_tot", "Total strain (-)")
    out_op.add_output_field("eps_vp", "Viscoplastic strain (-)")
    out_op.add_output_field("alpha", "Hardening parameter (-)")
    out_op.add_output_field("Fvp", "Yield function (-)")
    out_op.add_output_field("p_elems", "Mean stress (Pa)")
    out_op.add_output_field("q_elems", "Von Mises stress (Pa)")

    # Write schedule for plotter
    if MPI.COMM_WORLD.rank == 0:
        sched_path = os.path.join(output_folder_operation, "pressure_schedule.json")
        sched = {"time_list": list(map(float, time_list)),
                 "pressure": list(map(float, p_vals))}
        os.makedirs(output_folder_operation, exist_ok=True)
        with open(sched_path, "w") as f:
            json.dump(sched, f)
        print(f"[schedule] wrote {sched_path} with {len(time_list)} points")

    sf.Simulator_M(mom_eq, tc_operation, [out_op], False).run()


if __name__ == '__main__':
    main()
