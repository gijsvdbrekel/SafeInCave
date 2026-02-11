import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import dolfinx as do
import os
import sys
import torch as to
import json
import math
import numpy as np


# ══════════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

Z_SURFACE = 660.0
DAY_H = 24.0

# Cavern: regular 1200k
CAVERN_KEY = "regular1200"
GRID_FOLDER = "cavern_regular_1200_3D"
Z_MAX = 393.21
CAVERN_HEIGHT = 200.0
P_REF_MPA = 19.3

# Leaching
LEACHING_MODE = "stepped"
LEACHING_DAYS = 180
LEACHING_DT_HOURS = 12
STEPPED_N_STEPS = 6
LEACHING_END_FRACTION = 0.40

# Debrining
DEBRINING_DAYS = 20

# Fade-in
RAMP_UP_HOURS = 336  # 2 weeks

# Operation
OPERATION_DAYS = 365
DT_HOURS = 2

# Material flags
USE_DESAI = True
USE_THERMAL = False

# ── SWING SELECTION ──────────────────────────────────────────────────────────
# Set the swing value in bar for this run. Change this before each run.
# Valid values: 10, 12, 15, 16, 18, 20, 30
SWING_BAR = 10

# All swing values for reference (used only when running all sequentially via --all)
SWING_VALUES_MPA = [1.0, 1.2, 1.5, 1.6, 1.8, 2.0, 3.0]


# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS (copied from Run.py for standalone operation)
# ══════════════════════════════════════════════════════════════════════════════

def compute_lithostatic_pressure(z_center, p_ref_mpa, rho_salt, g):
    """Compute lithostatic pressure at cavern center."""
    depth = Z_SURFACE - z_center
    p_lithostatic = p_ref_mpa * ut.MPa + rho_salt * abs(g) * depth
    return p_lithostatic


def build_leaching_pressure_schedule(tc, *, p_start_pa, p_end_pa, mode, n_steps=6):
    """Build pressure schedule for leaching phase."""
    n_time_steps = int(math.floor(tc.t_final / tc.dt))
    t_vals = [k * tc.dt for k in range(n_time_steps + 1)]
    if abs(t_vals[-1] - tc.t_final) > 1e-12:
        t_vals.append(tc.t_final)

    if mode == "linear":
        p_vals = []
        for t in t_vals:
            frac = t / tc.t_final if tc.t_final > 0 else 1.0
            p = p_start_pa + frac * (p_end_pa - p_start_pa)
            p_vals.append(p)

    elif mode == "stepped":
        step_duration = tc.t_final / n_steps
        p_step_values = np.linspace(p_start_pa, p_end_pa, n_steps + 1)

        p_vals = []
        for t in t_vals:
            step_idx = min(int(t / step_duration), n_steps - 1)
            if t < tc.t_final:
                p = p_step_values[step_idx]
            else:
                p = p_end_pa
            p_vals.append(p)

    else:
        raise ValueError(f"Unknown leaching mode: {mode}")

    return t_vals, p_vals


def prepend_debrining(t_pressure, p_pressure, *, p_leach_end_pa, debrining_days):
    """Prepend debrining plateau to operational schedule."""
    debrining_s = debrining_days * 24.0 * 3600.0

    if debrining_s <= 0.0:
        return t_pressure, p_pressure

    t_pre = [0.0, debrining_s]
    p_pre = [p_leach_end_pa, p_leach_end_pa]

    t_shifted = [t + debrining_s for t in t_pressure[1:]]
    p_shifted = list(p_pressure[1:])

    return t_pre + t_shifted, p_pre + p_shifted


def apply_fade_in(t_pressure, p_pressure, *, p_start_pa, fade_in_hours):
    """Apply smooth cosine fade-in to the beginning of a pressure schedule."""
    if fade_in_hours <= 0.0:
        return

    fade_in_s = fade_in_hours * 3600.0

    for i in range(len(t_pressure)):
        t = t_pressure[i]
        if t >= fade_in_s:
            break
        alpha = 0.5 * (1.0 - math.cos(math.pi * t / fade_in_s))
        p_pressure[i] = (1.0 - alpha) * p_start_pa + alpha * p_pressure[i]


def _sample_at_dt(tc, t_end=None):
    t_end = tc.t_final if t_end is None else t_end
    n_steps = int(math.floor(t_end / tc.dt))
    t_vals = [k * tc.dt for k in range(n_steps + 1)]
    if abs(t_vals[-1] - t_end) > 1e-12:
        t_vals.append(t_end)
    return t_vals


def _repeat_hours(times_h, days):
    times_h = list(map(float, times_h))
    out = []
    for d in range(int(days)):
        off = d * DAY_H
        for i, t in enumerate(times_h):
            if d > 0 and i == 0 and abs(t - 0.0) < 1e-12 and abs(times_h[-1] - DAY_H) < 1e-12:
                continue
            out.append(off + t)
    return out


def build_linear_schedule_multi(tc, times_h, pressures_MPa, *, days, mode,
                                resample_at_dt=True, total_cycles=None):
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
                    continue
                t_scaled = offset + (t - base_start) * scale
                t_h.append(t_scaled)
                p_h.append(pressures_MPa[i])

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


# ══════════════════════════════════════════════════════════════════════════════
# MOMENTUM EQUATION CLASSES
# ══════════════════════════════════════════════════════════════════════════════

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

        if hasattr(st, "eps_ne_k"):
            self.eps_vp.x.array[:] = to.flatten(st.eps_ne_k)

        if self.expect_vp_state:
            if not (hasattr(st, "Fvp") and hasattr(st, "alpha")):
                if MPI.COMM_WORLD.rank == 0:
                    print("[WARN] Expected Fvp/alpha but missing.")
                return
            self.Fvp.x.array[:] = st.Fvp
            self.alpha.x.array[:] = st.alpha
        else:
            if hasattr(st, "Fvp") and hasattr(st, "alpha"):
                self.Fvp.x.array[:] = st.Fvp
                self.alpha.x.array[:] = st.alpha


class SparseSaveFields(sf.SaveFields):
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


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def run_single_swing(swing_mpa):
    """Run a complete leaching + operation simulation for one swing amplitude.

    Pressure profile: symmetric triangle with 48h period.
      - 24 hours ramp UP from p_min to p_min + swing_mpa
      - 24 hours ramp DOWN from p_min + swing_mpa back to p_min
    For swing_bar >= 30 the time step is reduced to 1 hour.
    """
    swing_bar = int(round(swing_mpa * 10))
    dt_hours = 1.0 if swing_bar >= 30 else DT_HOURS

    salt_density = 2200
    g = -9.81
    g_vec = [0.0, 0.0, g]

    z_center = Z_MAX - CAVERN_HEIGHT / 2.0

    p_lithostatic = compute_lithostatic_pressure(z_center, P_REF_MPA, salt_density, g)
    p_lithostatic_mpa = p_lithostatic / ut.MPa
    p_leach_end_mpa = LEACHING_END_FRACTION * p_lithostatic_mpa
    p_leach_end = p_leach_end_mpa * ut.MPa

    output_folder = os.path.join("output", f"case_swing_{swing_bar}bar_regular1200")

    if MPI.COMM_WORLD.rank == 0:
        print("=" * 70)
        print(f"PRESSURE SWING STUDY — {swing_bar} bar/day ({swing_mpa:.1f} MPa/day)")
        print("=" * 70)
        print(f"  Cavern:           regular (1200k m³)")
        print(f"  Cavern z_max:     {Z_MAX:.2f} m")
        print(f"  Cavern z_center:  {z_center:.2f} m")
        print(f"  P_ref (at z=660): {P_REF_MPA:.2f} MPa")
        print(f"  P_lithostatic:    {p_lithostatic_mpa:.2f} MPa (at cavern center)")
        print(f"  P_leach_end:      {p_leach_end_mpa:.2f} MPa ({LEACHING_END_FRACTION*100:.0f}%)")
        print(f"  Swing:            {swing_mpa:.1f} MPa = {swing_bar} bar")
        print(f"  Profile:          48h triangle (24h up, 24h down)")
        print(f"  P_min:            {p_leach_end_mpa:.2f} MPa")
        print(f"  P_max:            {p_leach_end_mpa + swing_mpa:.2f} MPa")
        print(f"  dt:               {dt_hours:.1f} hours")
        print(f"  Output:           {output_folder}")
        print("=" * 70)

    os.makedirs(output_folder, exist_ok=True)

    # ── Load grid ────────────────────────────────────────────────────────────
    grid_path = os.path.join("..", "..", "..", "..", "grids", GRID_FOLDER)
    grid = sf.GridHandlerGMSH("geom", grid_path)

    p_ref = P_REF_MPA * ut.MPa
    side_burden = p_ref
    over_burden = p_ref
    gas_density = 0.089

    # ── Momentum equation ────────────────────────────────────────────────────
    mom_eq = LinearMomentumMod(grid, theta=0.5)

    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("cg")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-10, max_it=100)
    mom_eq.set_solver(mom_solver)

    # ── Material ─────────────────────────────────────────────────────────────
    mat = sf.Material(mom_eq.n_elems)
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    E0 = 20.425 * ut.GPa * to.ones(mom_eq.n_elems)
    nu0 = 0.25 * to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")

    eta = 105e11 * to.ones(mom_eq.n_elems)
    E1 = 10 * ut.GPa * to.ones(mom_eq.n_elems)
    nu1 = 0.25 * to.ones(mom_eq.n_elems)
    kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")

    sec_per_year = 365.25 * 24 * 3600
    ndc = 4.6
    A_dc = (40.0 * (1e-6)**ndc / sec_per_year) * to.ones(mom_eq.n_elems)
    Q_dc = (6495.0 * 8.32) * to.ones(mom_eq.n_elems)
    n_dc = ndc * to.ones(mom_eq.n_elems)
    creep_0 = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation")

    A_ps = (14176.0 * 1e-9 / 1e6 / sec_per_year) * to.ones(mom_eq.n_elems)
    d_ps = 5.25e-3 * to.ones(mom_eq.n_elems)
    Q_ps = (3252.0 * 8.32) * to.ones(mom_eq.n_elems)
    creep_pressure = sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "creep_pressure")

    mat.add_to_elastic(spring_0)
    mat.add_to_non_elastic(kelvin)
    mat.add_to_non_elastic(creep_0)
    mat.add_to_non_elastic(creep_pressure)

    mom_eq.set_material(mat)
    mom_eq.build_body_force(g_vec)

    T0_field = 298 * to.ones(mom_eq.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # ═════ LEACHING PHASE ═════
    tc_init = sf.TimeController(
        dt=LEACHING_DT_HOURS,
        initial_time=0.0,
        final_time=LEACHING_DAYS * 24.0,
        time_unit="hour"
    )

    t_init, p_init = build_leaching_pressure_schedule(
        tc_init,
        p_start_pa=p_lithostatic,
        p_end_pa=p_leach_end,
        mode=LEACHING_MODE,
        n_steps=STEPPED_N_STEPS
    )

    save_interval_init = max(1, int(72 / LEACHING_DT_HOURS))

    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_init.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_init.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_init.t_final])

    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0,
                              [side_burden, side_burden],
                              [0.0, tc_init.t_final], g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0,
                               [side_burden, side_burden],
                               [0.0, tc_init.t_final], g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                             [over_burden, over_burden],
                             [0.0, tc_init.t_final], g=g_vec[2])
    bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, Z_MAX,
                                p_init, t_init, g=g_vec[2])

    bc_init = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
        bc_init.add_boundary_condition(bc)

    mom_eq.set_boundary_conditions(bc_init)

    output_folder_init = os.path.join(output_folder, "leaching")
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[LEACHING] Output: {output_folder_init}")

    output_mom_init = SparseSaveFields(mom_eq, interval=save_interval_init)
    output_mom_init.set_output_folder(output_folder_init)
    output_mom_init.add_output_field("u", "Displacement (m)")
    output_mom_init.add_output_field("eps_tot", "Total strain (-)")
    output_mom_init.add_output_field("sig", "Stress (Pa)")
    output_mom_init.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom_init.add_output_field("q_elems", "Von Mises stress (Pa)")

    sim_init = sf.Simulator_M(mom_eq, tc_init, [output_mom_init], True)
    sim_init.run()

    if MPI.COMM_WORLD.rank == 0:
        print("[LEACHING] Complete.")

    # ═════ ADD DESAI ═════
    if USE_DESAI:
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
        mom_eq.expect_vp_state = True
    else:
        mom_eq.expect_vp_state = False

    # ═════ BUILD TRIANGULAR WAVE PRESSURE ═════
    # 48h-period symmetric triangle: 24h ramp UP, 24h ramp DOWN
    p_min_mpa = p_leach_end_mpa
    p_max_mpa = p_min_mpa + swing_mpa

    tc_cycling = sf.TimeController(
        dt=dt_hours,
        initial_time=0.0,
        final_time=OPERATION_DAYS * 24.0,
        time_unit="hour"
    )

    # Build knots for the triangle wave
    n_full_cycles = OPERATION_DAYS // 2
    remainder_days = OPERATION_DAYS % 2

    knots_t_h = []
    knots_p_mpa = []
    for c in range(n_full_cycles):
        t_off = c * 48.0
        if c == 0:
            knots_t_h.append(t_off + 0.0)
            knots_p_mpa.append(p_min_mpa)
        knots_t_h.append(t_off + 24.0)
        knots_p_mpa.append(p_max_mpa)
        knots_t_h.append(t_off + 48.0)
        knots_p_mpa.append(p_min_mpa)

    if remainder_days > 0:
        t_off = n_full_cycles * 48.0
        if n_full_cycles == 0:
            knots_t_h.append(t_off)
            knots_p_mpa.append(p_min_mpa)
        knots_t_h.append(t_off + 24.0)
        knots_p_mpa.append(p_max_mpa)

    knots_t_s = np.array(knots_t_h) * ut.hour
    knots_p_pa = np.array(knots_p_mpa) * ut.MPa

    t_pressure = _sample_at_dt(tc_cycling)
    p_pressure = np.interp(t_pressure, knots_t_s, knots_p_pa).tolist()

    # Apply smooth fade-in
    if RAMP_UP_HOURS > 0:
        apply_fade_in(t_pressure, p_pressure,
                      p_start_pa=p_leach_end,
                      fade_in_hours=RAMP_UP_HOURS)
        if MPI.COMM_WORLD.rank == 0:
            print(f"[TRANSITION] Fade-in: {RAMP_UP_HOURS:.1f} hours from {p_leach_end_mpa:.2f} MPa")

    # Prepend debrining
    extra_hours = 0.0
    if DEBRINING_DAYS > 0:
        extra_hours = DEBRINING_DAYS * 24.0
        t_pressure, p_pressure = prepend_debrining(
            t_pressure, p_pressure,
            p_leach_end_pa=p_leach_end,
            debrining_days=DEBRINING_DAYS
        )
        if MPI.COMM_WORLD.rank == 0:
            print(f"[TRANSITION] Debrining: {DEBRINING_DAYS} days at {p_leach_end_mpa:.2f} MPa")

    # ═════ OPERATION PHASE ═════
    total_operation_hours = OPERATION_DAYS * 24.0 + extra_hours
    tc_operation = sf.TimeController(
        dt=dt_hours,
        initial_time=0.0,
        final_time=total_operation_hours,
        time_unit="hour"
    )

    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_operation.t_final])

    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0,
                              [side_burden, side_burden],
                              [0.0, tc_operation.t_final], g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0,
                               [side_burden, side_burden],
                               [0.0, tc_operation.t_final], g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                             [over_burden, over_burden],
                             [0.0, tc_operation.t_final], g=g_vec[2])
    bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, Z_MAX,
                                p_pressure, t_pressure, g=g_vec[2])

    bc_operation = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
        bc_operation.add_boundary_condition(bc)

    mom_eq.set_boundary_conditions(bc_operation)

    output_folder_operation = os.path.join(output_folder, "operation")
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[OPERATION] Output: {output_folder_operation}")

    # Save pressure schedule JSON
    pressure_data = {
        "cavern_key": CAVERN_KEY,
        "cavern_label": "regular (1200k m³)",
        "use_leaching": True,
        "leaching_mode": LEACHING_MODE,
        "leaching_days": LEACHING_DAYS,
        "leaching_end_fraction": LEACHING_END_FRACTION,
        "debrining_days": DEBRINING_DAYS,
        "ramp_up_hours": RAMP_UP_HOURS,
        "scenario": "triangle_48h",
        "mode": "triangle_48h",
        "swing_mpa": swing_mpa,
        "swing_bar": swing_bar,
        "operation_days": OPERATION_DAYS,
        "dt_hours": dt_hours,
        "use_desai": USE_DESAI,
        "p_lithostatic_mpa": p_lithostatic_mpa,
        "p_leach_end_mpa": p_leach_end_mpa,
        "p_min_mpa": p_min_mpa,
        "p_max_mpa": p_min_mpa + swing_mpa,
        "units": {"t_raw": "s", "p_raw": "Pa", "t": "hour", "p": "MPa"},
        "t_values_s": [float(t) for t in t_pressure],
        "p_values_Pa": [float(p) for p in p_pressure],
        "t_hours": [float(t / ut.hour) for t in t_pressure],
        "p_MPa": [float(p / ut.MPa) for p in p_pressure],
    }

    with open(os.path.join(output_folder, "pressure_schedule.json"), 'w') as f:
        json.dump(pressure_data, f, indent=2)

    output_mom_op = SparseSaveFields(mom_eq, interval=15)
    output_mom_op.set_output_folder(output_folder_operation)
    output_mom_op.add_output_field("u", "Displacement (m)")
    output_mom_op.add_output_field("eps_tot", "Total strain (-)")
    output_mom_op.add_output_field("eps_vp", "Viscoplastic strain (-)")
    output_mom_op.add_output_field("alpha", "Hardening parameter (-)")
    output_mom_op.add_output_field("Fvp", "Yield function (-)")
    output_mom_op.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom_op.add_output_field("q_elems", "Von Mises stress (Pa)")
    output_mom_op.add_output_field("sig", "Stress (Pa)")

    sim_op = sf.Simulator_M(mom_eq, tc_operation, [output_mom_op], False)
    sim_op.run()

    if MPI.COMM_WORLD.rank == 0:
        extra_d = DEBRINING_DAYS
        total_days = LEACHING_DAYS + extra_d + OPERATION_DAYS
        print("[OPERATION] Complete.")
        print("=" * 70)
        print(f"SWING {swing_bar} bar FINISHED")
        print(f"  Total simulated time: {total_days:.1f} days")
        print(f"  Output folder: {output_folder}")
        print("=" * 70)


def main():
    run_all = "--all" in sys.argv

    if run_all:
        # Sequential mode: run all swing values one after another
        valid_bar = [int(round(s * 10)) for s in SWING_VALUES_MPA]
        if MPI.COMM_WORLD.rank == 0:
            print("#" * 70)
            print("# PRESSURE SWING PARAMETER STUDY — ALL CASES")
            print(f"# Swing values: {valid_bar} bar/day")
            print(f"# Cavern: regular 1200k m³")
            print("#" * 70)

        for swing_mpa in SWING_VALUES_MPA:
            run_single_swing(swing_mpa)
    else:
        # Default: run the single SWING_BAR value set at the top of this file
        swing_mpa = SWING_BAR / 10.0
        if MPI.COMM_WORLD.rank == 0:
            print("#" * 70)
            print(f"# PRESSURE SWING — {SWING_BAR} bar/day")
            print(f"# Cavern: regular 1200k m³")
            print("#" * 70)

        run_single_swing(swing_mpa)

    if MPI.COMM_WORLD.rank == 0:
        print("\n" + "#" * 70)
        print("# COMPLETE")
        print("#" * 70)


if __name__ == '__main__':
    main()
