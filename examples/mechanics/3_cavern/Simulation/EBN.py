"""
EBN_CSV.py
==========
Volledig simulatiescript met CSV-drukprofiel, ZONDER pandas (dus geen "module pandas not found"),
met:
- GMRES + ASM (ipv CG)
- continuity fix: p_operation[0] = p_eq
- startup ramp (bijv. 24 uur) om NaNs bij start operation te voorkomen
- duidelijke prints die NA de fix/ramp komen

PAS AAN BOVENAAN:
- CSV_FILE_PATH
- OPERATION_DAYS, dt_hours
- SCHEDULE_MODE: "direct" | "stretch" | "repeat"
- RAMP_HOURS
- (optioneel) USE_DESAI tijdelijk False zetten om te isoleren
"""

import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC

from petsc4py import PETSc
from mpi4py import MPI
import dolfinx as do

import os
import json
import math
import torch as to
import numpy as np
import csv


# =========================
# USER SETTINGS
# =========================
CAVERN_TYPE = "regular600"

PRESSURE_SCENARIO = "csv_profile"   # fixed here
SCHEDULE_MODE = "direct"            # "direct" | "stretch" | "repeat"
OPERATION_DAYS = 365
dt_hours = 2.0

# Start-up smoothing to avoid NaNs at operation start
RAMP_HOURS = 24.0                   # try 6, 12, 24, 48
CLAMP_MIN_MPA = 6.0
CLAMP_MAX_MPA = 20.0

# Equilibrium pressure (MPa) must match your equilibrium BC cavern pressure
P_EQ_MPA = 15.0

# If you suspect Desai triggers NaNs, set False to test
USE_DESAI = True

# CSV file path (relative or absolute)
CSV_FILE_PATH = "drukprofiel_zoutcaverne_2035_8760u.csv"

# Output
OUTPUT_ROOT = "output"


# =========================
# Helpers / Units
# =========================
DAY_H = 24.0


def make_solver(comm):
    ksp = PETSc.KSP().create(comm)
    ksp.setType("gmres")
    pc = ksp.getPC()
    pc.setType("asm")
    ksp.setTolerances(rtol=1e-10, max_it=200)
    return ksp


class LinearMomentumMod(sf.LinearMomentum):
    """
    Minimal extension to keep Desai state fields (if present).
    """
    def __init__(self, grid, theta):
        super().__init__(grid, theta)
        self.expect_vp_state = False

    def initialize(self) -> None:
        self.C.x.array[:] = to.flatten(self.mat.C)
        self.Fvp = do.fem.Function(self.DG0_1)
        self.alpha = do.fem.Function(self.DG0_1)
        self.eps_vp = do.fem.Function(self.DG0_3x3)

    def run_after_solve(self):
        elems = getattr(self.mat, "elems_ne", None)
        if not elems:
            return
        st = elems[-1]

        if hasattr(st, "eps_ne_k"):
            self.eps_vp.x.array[:] = to.flatten(st.eps_ne_k)

        if self.expect_vp_state and hasattr(st, "Fvp") and hasattr(st, "alpha"):
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


# =========================
# CSV parsing (no pandas)
# =========================
def _parse_float_auto(s: str) -> float:
    """
    Parses floats with either '.' or ',' decimal. Strips spaces.
    """
    s = s.strip()
    if not s:
        return np.nan
    # If there is a comma and no dot, treat comma as decimal separator
    if "," in s and "." not in s:
        s = s.replace(",", ".")
    # remove thousand separators like "1.234,56" (rare). Best effort:
    # if both exist, assume dot is thousand and comma decimal.
    if "," in s and "." in s:
        s = s.replace(".", "").replace(",", ".")
    try:
        return float(s)
    except Exception:
        return np.nan


def read_pressure_csv(csv_file: str):
    """
    Reads CSV and returns pressure array in MPa.
    Supported columns (case-insensitive):
      - Druk_MPa
      - Druk_bar  (converted to MPa via /10)
    If header not recognized, tries first numeric column.
    """
    if not os.path.isfile(csv_file):
        raise FileNotFoundError(f"CSV not found: {csv_file}")

    # sniff delimiter
    with open(csv_file, "r", newline="", encoding="utf-8") as f:
        sample = f.read(2048)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=";, \t,")
            delim = dialect.delimiter
        except Exception:
            delim = ";"

    rows = []
    with open(csv_file, "r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=delim)
        for r in reader:
            if len(r) == 0:
                continue
            rows.append(r)

    if len(rows) < 2:
        raise ValueError("CSV has too few rows")

    header = [c.strip() for c in rows[0]]
    data = rows[1:]

    header_low = [h.lower() for h in header]

    idx_mpa = None
    idx_bar = None

    for i, h in enumerate(header_low):
        if h == "druk_mpa" or h.endswith("druk_mpa"):
            idx_mpa = i
        if h == "druk_bar" or h.endswith("druk_bar"):
            idx_bar = i

    pressures_mpa = []

    if idx_mpa is not None:
        for r in data:
            if idx_mpa >= len(r):
                continue
            pressures_mpa.append(_parse_float_auto(r[idx_mpa]))
        pressures_mpa = np.asarray(pressures_mpa, dtype=float)

    elif idx_bar is not None:
        for r in data:
            if idx_bar >= len(r):
                continue
            v_bar = _parse_float_auto(r[idx_bar])
            pressures_mpa.append(v_bar / 10.0)
        pressures_mpa = np.asarray(pressures_mpa, dtype=float)

    else:
        # fallback: find first column with many numeric values
        ncols = len(header)
        best_i = None
        best_count = -1
        for i in range(ncols):
            vals = [_parse_float_auto(r[i]) for r in data if i < len(r)]
            arr = np.asarray(vals, dtype=float)
            count = np.isfinite(arr).sum()
            if count > best_count:
                best_count = count
                best_i = i
        if best_i is None or best_count < 2:
            raise ValueError("Could not find a numeric pressure column in CSV")
        for r in data:
            if best_i >= len(r):
                continue
            pressures_mpa.append(_parse_float_auto(r[best_i]))
        pressures_mpa = np.asarray(pressures_mpa, dtype=float)

    pressures_mpa = pressures_mpa[np.isfinite(pressures_mpa)]
    if pressures_mpa.size < 2:
        raise ValueError("Parsed pressure series has <2 numeric values")

    return pressures_mpa


# =========================
# Build schedule from CSV
# =========================
def build_csv_pressure_schedule(tc, csv_file, *,
                                days, mode="direct",
                                total_cycles=1,
                                clamp_min_mpa=None, clamp_max_mpa=None,
                                resample_at_dt=True):
    """
    Returns:
      t_vals: list of seconds (float)
      p_vals: list of Pa (float)
    """
    druk_values_MPa = read_pressure_csv(csv_file)
    csv_hours = int(druk_values_MPa.size)

    if MPI.COMM_WORLD.rank == 0:
        print(f"[CSV] Loaded '{os.path.basename(csv_file)}' with {csv_hours} hourly values")
        print(f"[CSV] Pressure range: {druk_values_MPa.min():.2f} - {druk_values_MPa.max():.2f} MPa")

    total_hours = float(days) * 24.0

    if mode == "direct":
        # sample hourly, periodic over csv_hours
        sim_hours = np.arange(0.0, total_hours + 1e-12, 1.0)
        idx = (sim_hours % csv_hours).astype(int)
        times_hours = sim_hours
        pressures_MPa = druk_values_MPa[idx]

    elif mode == "stretch":
        total_cycles = max(1, int(total_cycles))
        cycle_duration_h = total_hours / float(total_cycles)
        scale = cycle_duration_h / float(csv_hours)

        times_list = []
        pres_list = []
        for k in range(total_cycles):
            off = k * cycle_duration_h
            for i in range(csv_hours):
                if k > 0 and i == 0:
                    continue
                times_list.append(off + i * scale)
                pres_list.append(druk_values_MPa[i])
        times_hours = np.asarray(times_list, dtype=float)
        pressures_MPa = np.asarray(pres_list, dtype=float)

    elif mode == "repeat":
        # repeat year-pattern until reach total_hours
        n_rep = int(np.ceil(total_hours / float(csv_hours)))
        times_list = []
        pres_list = []
        for r in range(n_rep):
            off = r * csv_hours
            for i in range(csv_hours):
                if r > 0 and i == 0:
                    continue
                t = off + i
                if t > total_hours:
                    break
                times_list.append(t)
                pres_list.append(druk_values_MPa[i])
        times_hours = np.asarray(times_list, dtype=float)
        pressures_MPa = np.asarray(pres_list, dtype=float)

    else:
        raise ValueError("mode must be 'direct', 'stretch', or 'repeat'")

    # clamp in MPa
    if clamp_min_mpa is not None:
        pressures_MPa = np.maximum(pressures_MPa, float(clamp_min_mpa))
    if clamp_max_mpa is not None:
        pressures_MPa = np.minimum(pressures_MPa, float(clamp_max_mpa))

    times_s = times_hours * 3600.0

    # enforce start/end
    if times_s[0] > 0.0:
        times_s = np.insert(times_s, 0, 0.0)
        pressures_MPa = np.insert(pressures_MPa, 0, pressures_MPa[0])
    if times_s[-1] < tc.t_final:
        times_s = np.append(times_s, tc.t_final)
        pressures_MPa = np.append(pressures_MPa, pressures_MPa[-1])

    if resample_at_dt:
        n_steps = int(math.floor(tc.t_final / tc.dt))
        t_vals = [k * tc.dt for k in range(n_steps + 1)]
        if abs(t_vals[-1] - tc.t_final) > 1e-12:
            t_vals.append(tc.t_final)
        p_vals_MPa = np.interp(t_vals, times_s, pressures_MPa)
    else:
        t_vals = times_s.tolist()
        p_vals_MPa = pressures_MPa.tolist()

    p_vals = [float(p) * ut.MPa for p in p_vals_MPa]

    if MPI.COMM_WORLD.rank == 0:
        print(f"[CSV] Schedule built: {len(t_vals)} points, mode={mode}")

    return t_vals, p_vals


def apply_startup_ramp(t_pressure, p_pressure, *, p_start_pa, ramp_hours, dt_hours):
    """
    Replace first part of schedule with a linear ramp from p_start_pa to the existing schedule.
    Operates in-place on p_pressure list.
    """
    if ramp_hours is None or ramp_hours <= 0.0:
        p_pressure[0] = p_start_pa
        return

    ramp_steps = max(1, int(round(float(ramp_hours) / float(dt_hours))))
    ramp_steps = min(ramp_steps, len(p_pressure) - 1)

    # target after ramp window
    p_target = p_pressure[ramp_steps]

    p_pressure[0] = p_start_pa
    for k in range(1, ramp_steps + 1):
        a = k / float(ramp_steps)
        p_pressure[k] = (1.0 - a) * p_start_pa + a * p_target


# =========================
# Main
# =========================
def main():
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

    # burdens
    if CAVERN_TYPE.endswith("600"):
        p_ref = 18.2 * ut.MPa
    elif CAVERN_TYPE.endswith("1200"):
        p_ref = 20.1 * ut.MPa
    else:
        raise ValueError(f"Cannot infer cavern set from {CAVERN_TYPE}")

    side_burden = p_ref
    over_burden = p_ref

    # momentum
    mom_eq = LinearMomentumMod(grid, theta=0.5)
    mom_eq.set_solver(make_solver(grid.mesh.comm))

    # material (elastic + kelvin + disloc + pressure solution)
    mat = sf.Material(mom_eq.n_elems)

    salt_density = 2200
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    # spring
    E0 = 20.425 * ut.GPa * to.ones(mom_eq.n_elems)
    nu0 = 0.25 * to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")

    # kelvin
    eta = 105e11 * to.ones(mom_eq.n_elems)
    E1 = 10 * ut.GPa * to.ones(mom_eq.n_elems)
    nu1 = 0.25 * to.ones(mom_eq.n_elems)
    kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")

    sec_per_year = 365.25 * 24 * 3600

    # dislocation (old example)
    ndc = 4.6
    A_dc = (40.0 * (1e-6) ** ndc / sec_per_year) * to.ones(mom_eq.n_elems)
    Q_dc = (6495.0 * 8.32) * to.ones(mom_eq.n_elems)
    n_dc = ndc * to.ones(mom_eq.n_elems)
    creep_disloc = sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation")

    # pressure solution
    A_ps = (14176.0 * 1e-9 / 1e6 / sec_per_year) * to.ones(mom_eq.n_elems)
    d_ps = 5.25e-3 * to.ones(mom_eq.n_elems)
    Q_ps = (3252.0 * 8.32) * to.ones(mom_eq.n_elems)
    creep_pressure = sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "creep_pressure")

    mat.add_to_elastic(spring_0)
    mat.add_to_non_elastic(kelvin)
    mat.add_to_non_elastic(creep_disloc)
    mat.add_to_non_elastic(creep_pressure)

    mom_eq.set_material(mat)

    # body force + temperature
    g_vec = [0.0, 0.0, -9.81]
    mom_eq.build_body_force(g_vec)

    T0_field = 298 * to.ones(mom_eq.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # =========================
    # EQUILIBRIUM
    # =========================
    tc_equilibrium = sf.TimeController(dt=0.5, initial_time=0.0, final_time=10, time_unit="hour")

    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_equilibrium.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_equilibrium.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_equilibrium.t_final])

    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0, [side_burden, side_burden],
                              [0.0, tc_equilibrium.t_final], g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0, [side_burden, side_burden],
                               [0.0, tc_equilibrium.t_final], g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0, [over_burden, over_burden],
                             [0.0, tc_equilibrium.t_final], g=g_vec[2])

    gas_density = 0.089
    p_gas = float(P_EQ_MPA) * ut.MPa

    bc_cavern_eq = momBC.NeumannBC("Cavern", 2, gas_density, z_max, [p_gas, p_gas],
                                   [0.0, tc_equilibrium.t_final], g=g_vec[2])

    bc_equilibrium = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern_eq]:
        bc_equilibrium.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_equilibrium)

    output_folder = os.path.join(
        OUTPUT_ROOT,
        f"case_{PRESSURE_SCENARIO}_{SCHEDULE_MODE}_{OPERATION_DAYS}days_{CAVERN_TYPE}"
    )

    out_eq_dir = os.path.join(output_folder, "equilibrium")
    if MPI.COMM_WORLD.rank == 0:
        print("[EQ OUT]", out_eq_dir)

    out_eq = sf.SaveFields(mom_eq)
    out_eq.set_output_folder(out_eq_dir)
    out_eq.add_output_field("u", "Displacement (m)")
    out_eq.add_output_field("eps_tot", "Total strain (-)")
    out_eq.add_output_field("sig", "Stress (Pa)")
    out_eq.add_output_field("p_elems", "Mean stress (Pa)")
    out_eq.add_output_field("q_elems", "Von Mises stress (Pa)")

    sim_eq = sf.Simulator_M(mom_eq, tc_equilibrium, [out_eq], True)
    sim_eq.run()

    # =========================
    # OPERATION (optioneel Desai toevoegen)
    # =========================
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

        desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta_vp, n, beta_1, beta,
                                     m, gamma, sigma_t, alpha_0, "desai")

        stress_to = ut.numpy2torch(mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3)))
        desai.compute_initial_hardening(stress_to, Fvp_0=0.0)

        mat.add_to_non_elastic(desai)
        mom_eq.set_material(mat)
        mom_eq.expect_vp_state = True
    else:
        mom_eq.expect_vp_state = False

    tc_operation = sf.TimeController(
        dt=dt_hours,
        initial_time=0.0,
        final_time=OPERATION_DAYS * 24.0,
        time_unit="hour"
    )

    # --- build CSV schedule ---
    t_pressure, p_pressure = build_csv_pressure_schedule(
        tc_operation,
        csv_file=CSV_FILE_PATH,
        days=OPERATION_DAYS,
        mode=SCHEDULE_MODE,
        total_cycles=1,
        clamp_min_mpa=CLAMP_MIN_MPA,
        clamp_max_mpa=CLAMP_MAX_MPA,
        resample_at_dt=True
    )

    # --- continuity fix (must happen BEFORE BC creation) ---
    p_pressure[0] = float(P_EQ_MPA) * ut.MPa

    # --- startup ramp (high impact for NaN-at-start problems) ---
    apply_startup_ramp(
        t_pressure, p_pressure,
        p_start_pa=float(P_EQ_MPA) * ut.MPa,
        ramp_hours=RAMP_HOURS,
        dt_hours=dt_hours
    )

    if MPI.COMM_WORLD.rank == 0:
        print("[CHECK AFTER FIX/RAMP] p_eq (MPa):", float(P_EQ_MPA))
        print("[CHECK AFTER FIX/RAMP] p_op[0] (MPa):", float(p_pressure[0] / ut.MPa))
        if len(p_pressure) > 1:
            print("[CHECK AFTER FIX/RAMP] p_op[1] (MPa):", float(p_pressure[1] / ut.MPa))
        print("[CHECK AFTER FIX/RAMP] t_op[0] (h):", float(t_pressure[0] / ut.hour))
        if len(t_pressure) > 1:
            print("[CHECK AFTER FIX/RAMP] t_op[1] (h):", float(t_pressure[1] / ut.hour))

    # Operation BCs
    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_operation.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_operation.t_final])

    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0, [side_burden, side_burden],
                              [0.0, tc_operation.t_final], g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0, [side_burden, side_burden],
                               [0.0, tc_operation.t_final], g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0, [over_burden, over_burden],
                             [0.0, tc_operation.t_final], g=g_vec[2])

    bc_cavern_op = momBC.NeumannBC("Cavern", 2, gas_density, z_max, p_pressure, t_pressure, g=g_vec[2])

    bc_operation = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern_op]:
        bc_operation.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_operation)

    # Save schedule json
    os.makedirs(output_folder, exist_ok=True)
    pressure_data = {
        "scenario": PRESSURE_SCENARIO,
        "mode": SCHEDULE_MODE,
        "operation_days": OPERATION_DAYS,
        "csv_file": os.path.basename(CSV_FILE_PATH),
        "dt_hours": dt_hours,
        "ramp_hours": RAMP_HOURS,
        "clamp_min_mpa": CLAMP_MIN_MPA,
        "clamp_max_mpa": CLAMP_MAX_MPA,
        "p_eq_mpa": P_EQ_MPA,
        "units": {"t_raw": "s", "p_raw": "Pa", "t": "hour", "p": "MPa"},
        "t_values_s": [float(t) for t in t_pressure],
        "p_values_Pa": [float(p) for p in p_pressure],
        "t_hours": [float(t / ut.hour) for t in t_pressure],
        "p_MPa": [float(p / ut.MPa) for p in p_pressure],
        "use_desai": bool(USE_DESAI),
        "solver": {"ksp": "gmres", "pc": "asm", "rtol": 1e-10, "max_it": 200},
    }
    with open(os.path.join(output_folder, "pressure_schedule.json"), "w") as f:
        json.dump(pressure_data, f, indent=2)

    # Outputs
    out_op_dir = os.path.join(output_folder, "operation")
    if MPI.COMM_WORLD.rank == 0:
        print("[OP OUT]", out_op_dir)

    out_op = SparseSaveFields(mom_eq, interval=15)
    out_op.set_output_folder(out_op_dir)
    out_op.add_output_field("u", "Displacement (m)")
    out_op.add_output_field("eps_tot", "Total strain (-)")
    out_op.add_output_field("p_elems", "Mean stress (Pa)")
    out_op.add_output_field("q_elems", "Von Mises stress (Pa)")
    out_op.add_output_field("sig", "Stress (Pa)")

    # desai fields (will be default/empty if USE_DESAI=False)
    out_op.add_output_field("eps_vp", "Viscoplastic strain (-)")
    out_op.add_output_field("alpha", "Hardening parameter (-)")
    out_op.add_output_field("Fvp", "Yield function (-)")

    sim_op = sf.Simulator_M(mom_eq, tc_operation, [out_op], False)
    sim_op.run()


if __name__ == "__main__":
    main()
