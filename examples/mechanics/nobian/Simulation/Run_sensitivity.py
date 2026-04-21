"""
Run_sensitivity.py — single-case runner for sensitivity analyses.

Based on Run.py but driven by a dict of named cases. Each case defines
overrides onto a _BASELINE config. Exactly ONE case is executed per
invocation so that multiple PuTTY screens can launch different cases in
parallel. Writes a profiling.json alongside each output folder.

Profiling is non-invasive (monkey-patched method wrappers on mom_eq /
mom_eq.solver). Baseline simulation reproducibility is preserved.

Usage
-----
    python Run_sensitivity.py                    # runs CASE_TO_RUN from top of file
    python Run_sensitivity.py --case md_A_baseline
    python Run_sensitivity.py --list             # show available case names
    python Run_sensitivity.py --force            # rerun even if output exists
"""

import os
import sys
import json
import time
import copy
import math
import traceback
import numpy as np
import torch as to
from mpi4py import MPI
from petsc4py import PETSc
import dolfinx as do

import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
import safeincave.HeatBC as heatBC

from Run import (
    compute_lithostatic_pressure,
    build_leaching_pressure_schedule,
    prepend_debrining,
    apply_fade_in,
    build_csv_pressure_schedule,
    apply_startup_ramp,
    TimeControllerFromList,
    build_time_list_by_dp_limit,
    build_linear_schedule_multi,
    build_sinus_schedule_multi,
    build_power_generation_schedule,
    LinearMomentumMod,
    SparseSaveFields,
    Z_MAX_BY_CAVERN,
    CAVERN_HEIGHT_BY_TYPE,
    P_REF_BY_SIZE,
    P_REF_BY_CAVERN,
    GRID_FOLDERS,
    Z_SURFACE,
)


# ══════════════════════════════════════════════════════════════════════════════
# USER CONFIGURATION — edit SELECTED_CASES and FORCE_RERUN as needed.
# ══════════════════════════════════════════════════════════════════════════════

_BASELINE = {
    # Initialization
    "USE_LEACHING": True,
    "LEACHING_MODE": "linear",
    "LEACHING_DAYS": 91,
    "LEACHING_DT_HOURS": 12,
    "STEPPED_N_STEPS": 6,
    "LEACHING_END_FRACTION": 0.30,
    "DEBRINING_DAYS": 30,
    "RAMP_UP_HOURS": 336,
    "EQUILIBRIUM_HOURS": 10.0,
    "EQUILIBRIUM_DT_HOURS": 0.5,

    # Cavern
    "CAVERN_TYPE": "regular",
    "CAVERN_SIZE": 1200,

    # Pressure scenario
    "PRESSURE_SCENARIO": "industry",
    "P_MEAN_MPA": 15.0,
    "P_AMPLITUDE_MPA": 5.0,
    "P_HIGH_OFFSET_MPA": 5.0,
    "P_LOW_OFFSET_MPA": 1.5,
    "N_EVENTS": 8,
    "P_BASE_OFFSET_MPA": 10.0,
    "RECOVERY_TAU_HOURS": 200.0,
    "P_MIN_MPA": 8.5,
    "P_EQUILIBRIUM_MPA": 11.0,
    "CSV_FILE_PATH": "drukprofiel_zoutcaverne_2035_8760u.csv",
    "CSV_SHIFT_TO_LEACH_END": False,
    "RESCALE_PRESSURE": True,
    "RESCALE_MIN_MPA": 6.0,
    "RESCALE_MAX_MPA": 20.0,
    "RAMP_HOURS": 24.0,
    "AMPLITUDE_SCALE": 1.0,

    # Schedule
    "SCHEDULE_MODE": "stretch",
    "OPERATION_DAYS": 365,
    "N_CYCLES": 22,
    "dt_hours": 2.0,

    # Variable dt
    "USE_VARIABLE_DT": True,
    "DT_FINE_HOURS": 0.1,
    "DT_COARSE_HOURS": 2.0,
    "MAX_DP_MPA": 0.2,

    # Material
    "MATERIAL_SCENARIO": "B",
    "USE_MUNSON_DAWSON": True,
    "MATERIAL_OVERRIDES": {},

    # Thermal (disabled — thermal dimension not included in thesis scope)
    "USE_THERMAL": False,
}


SENSITIVITY_CASES = {
    # ── Dimension 1: MD creep mechanism toggles (thesis core) ──────────────
    "md_B_baseline":       {},
    "md_B_no_transient":   {"MATERIAL_OVERRIDES": {"md_alpha_w": 0.0, "md_beta_w": 0.0}},
    "md_B_no_reverse":     {"MATERIAL_OVERRIDES": {"md_delta": 0.0}},
    "md_B_steady_state":   {"MATERIAL_OVERRIDES": {"md_alpha_w": 0.0, "md_beta_w": 0.0, "md_delta": 0.0}},

    # ── Dimension 2: Constitutive model ─────────────────────────────────────
    "sic_B_baseline":      {"USE_MUNSON_DAWSON": False},
    "md_A_baseline":       {"MATERIAL_SCENARIO": "A"},

    # ── Dimension 3: Startup ────────────────────────────────────────────────
    "no_leaching":         {"USE_LEACHING": False},

    # ── Dimension 4: Geometry ───────────────────────────────────────────────
    "tilted_cavern":       {"CAVERN_TYPE": "tilted", "CAVERN_SIZE": 1200},
    "volume_600":          {"CAVERN_SIZE": 600},

    # ── Dimension 5: Secondary creep contributions ──────────────────────────
    "no_pressure_solution": {"MATERIAL_OVERRIDES": {"disable_pressure_solution": True}},
    "no_kelvin":            {"USE_MUNSON_DAWSON": False, "MATERIAL_OVERRIDES": {"disable_kelvin": True}},
    "no_desai":             {"USE_MUNSON_DAWSON": False, "MATERIAL_OVERRIDES": {"disable_desai": True}},

    # ── Dimension 6: Numerical dt ───────────────────────────────────────────
    "dt_half":             {"dt_hours": 1.0},
    "dt_double":           {"dt_hours": 4.0},

    # ── Dimension 7: Pressure amplitude ─────────────────────────────────────
    "amplitude_half":      {"AMPLITUDE_SCALE": 0.5},
    "amplitude_double":    {"AMPLITUDE_SCALE": 2.0},

    # ── Smoke-test (minutes, not hours) — use to validate infrastructure ────
    "smoke": {
        "USE_LEACHING": False,
        "OPERATION_DAYS": 2,
        "dt_hours": 6.0,
        "N_CYCLES": 2,
        "RAMP_UP_HOURS": 12.0,
        "EQUILIBRIUM_HOURS": 6.0,
        "EQUILIBRIUM_DT_HOURS": 2.0,
    },
}


# ── WHICH CASE TO RUN ─────────────────────────────────────────────────────────
# Exactly ONE case is executed per script invocation. Open multiple PuTTY
# screens to run cases in parallel — each screen runs its own case.
#
# Switch the case by either:
#   (a) editing CASE_TO_RUN below, OR
#   (b) passing --case <name> on the command line (overrides the variable):
#         python Run_sensitivity.py --case md_B_no_transient
#
# Other CLI helpers:
#   python Run_sensitivity.py --list       # print all available case names
#   python Run_sensitivity.py --force      # rerun even if output exists

CASE_TO_RUN = "md_B_baseline"

# If False, skip cases whose output already has both log.txt ("Total time" line)
# and profiling.json — lets interrupted sweeps resume cheaply.
FORCE_RERUN = False

# ══════════════════════════════════════════════════════════════════════════════
# END OF USER CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════


VALID_SHAPES = ["regular", "tilted", "directcirculation", "asymmetric",
                "reversedcirculation", "fastleached", "tubefailure"]
VALID_SPECIAL_CAVERNS = ["A5"]
VALID_SCENARIOS = ["industry", "transport", "power_generation", "csv"]


def resolve_case(name):
    if name not in SENSITIVITY_CASES:
        raise KeyError(f"Unknown case '{name}'. Available: {sorted(SENSITIVITY_CASES)}")
    return {**_BASELINE, **SENSITIVITY_CASES[name]}


def resolve_cavern_metadata(cfg):
    ctype = cfg["CAVERN_TYPE"]
    csize = cfg["CAVERN_SIZE"]
    is_special = ctype in VALID_SPECIAL_CAVERNS
    if is_special:
        cavern_key = ctype
        p_ref_mpa = P_REF_BY_CAVERN.get(cavern_key, 17.5)
        cavern_label = f"{ctype} (~965k m³)"
    else:
        if ctype not in VALID_SHAPES:
            raise ValueError(f"CAVERN_TYPE '{ctype}' invalid")
        cavern_key = f"{ctype}{csize}"
        p_ref_mpa = P_REF_BY_SIZE[csize]
        cavern_label = f"{ctype} ({csize}k m³)"
    z_max = Z_MAX_BY_CAVERN[cavern_key]
    cavern_height = CAVERN_HEIGHT_BY_TYPE[cavern_key]
    z_center = z_max - cavern_height / 2.0
    grid_folder = GRID_FOLDERS[cavern_key]
    return {
        "cavern_key": cavern_key,
        "cavern_label": cavern_label,
        "grid_folder": grid_folder,
        "z_max": z_max,
        "z_center": z_center,
        "cavern_height": cavern_height,
        "p_ref_mpa": p_ref_mpa,
        "is_special": is_special,
    }


def build_output_folder(name, cfg, meta):
    scen_tag = cfg["MATERIAL_SCENARIO"]
    model_tag = "MD" if cfg["USE_MUNSON_DAWSON"] else "SIC"
    if cfg["USE_LEACHING"]:
        prefix = (f"case_leaching_{cfg['LEACHING_MODE']}_{cfg['PRESSURE_SCENARIO']}"
                  f"({cfg['N_CYCLES']})_{cfg['OPERATION_DAYS']}days_S{scen_tag}_{model_tag}_"
                  f"{meta['cavern_key']}")
    else:
        prefix = (f"case_{cfg['PRESSURE_SCENARIO']}"
                  f"({cfg['N_CYCLES']})_{cfg['OPERATION_DAYS']}days_S{scen_tag}_{model_tag}_"
                  f"{meta['cavern_key']}")
    return os.path.join("output", f"{prefix}_sens_{name}")


def is_case_complete(folder):
    log_path = os.path.join(folder, "operation", "log.txt")
    prof_path = os.path.join(folder, "profiling.json")
    if not (os.path.isfile(log_path) and os.path.isfile(prof_path)):
        return False
    try:
        with open(log_path, "r") as f:
            tail = f.readlines()[-5:]
        return any("Total time" in line for line in tail)
    except Exception:
        return False


# ══════════════════════════════════════════════════════════════════════════════
# PROFILING
# ══════════════════════════════════════════════════════════════════════════════

class ProfilingRecorder:
    """
    Non-invasive CPU profiler for a single simulation case.

    Wraps bound methods on the LinearMomentum instance (and its KSP solver).
    Collects per-LinearMomentum.solve() timings: total / compute_CT / ksp.solve.
    Assembly time is derived: assemble = total - ct - ksp - other_overhead.

    After the simulation, call finalize(log_path, wall_time_s) to compute
    summary stats and write profiling.json.
    """

    def __init__(self):
        self.per_iter = []
        self._scratch = {"ct": 0.0, "ksp": 0.0}
        self._originals = []
        self._installed = False

    def install(self, mom_eq):
        if self._installed:
            raise RuntimeError("ProfilingRecorder already installed")
        if not hasattr(mom_eq, "solver"):
            raise RuntimeError("mom_eq has no solver attribute — set_solver first")

        scratch = self._scratch
        orig_compute_CT = mom_eq.compute_CT
        self._originals.append((mom_eq, "compute_CT", orig_compute_CT))

        def wrapped_compute_CT(*a, **kw):
            t0 = time.perf_counter()
            try:
                return orig_compute_CT(*a, **kw)
            finally:
                scratch["ct"] += time.perf_counter() - t0
        mom_eq.compute_CT = wrapped_compute_CT

        orig_ksp_solve = mom_eq.solver.solve
        self._originals.append((mom_eq.solver, "solve", orig_ksp_solve))

        def wrapped_ksp_solve(*a, **kw):
            t0 = time.perf_counter()
            try:
                return orig_ksp_solve(*a, **kw)
            finally:
                scratch["ksp"] += time.perf_counter() - t0
        mom_eq.solver.solve = wrapped_ksp_solve

        orig_solve = mom_eq.solve
        self._originals.append((mom_eq, "solve", orig_solve))
        per_iter = self.per_iter

        def wrapped_solve(stress_k, t, dt):
            scratch["ct"] = 0.0
            scratch["ksp"] = 0.0
            t0 = time.perf_counter()
            try:
                orig_solve(stress_k, t, dt)
            finally:
                total = time.perf_counter() - t0
                per_iter.append({
                    "t_s": float(t),
                    "dt_s": float(dt),
                    "total_s": total,
                    "ct_s": scratch["ct"],
                    "ksp_s": scratch["ksp"],
                    "assemble_s": max(0.0, total - scratch["ct"] - scratch["ksp"]),
                })
        mom_eq.solve = wrapped_solve

        self._installed = True

    def uninstall(self):
        for obj, name, orig in reversed(self._originals):
            setattr(obj, name, orig)
        self._originals.clear()
        self._installed = False

    @staticmethod
    def _parse_log_newton_iters(log_path):
        """Extract Newton-iter count per time-step from Simulator_M log.txt."""
        newton_per_step = []
        if not os.path.isfile(log_path):
            return newton_per_step
        with open(log_path, "r") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("|") is False:
                    continue
                parts = [p.strip() for p in s.strip("|").split("|")]
                if len(parts) < 4:
                    continue
                try:
                    int(parts[0])
                except ValueError:
                    continue
                try:
                    newton_per_step.append(int(parts[3]))
                except (ValueError, IndexError):
                    pass
        return newton_per_step

    def finalize(self, case_name, output_folder, cfg, wall_time_s, op_log_path=None):
        per_iter = self.per_iter
        n_iter_total = len(per_iter)

        ksp_total = sum(x["ksp_s"] for x in per_iter)
        ct_total = sum(x["ct_s"] for x in per_iter)
        assemble_total = sum(x["assemble_s"] for x in per_iter)
        iter_total_s = sum(x["total_s"] for x in per_iter)
        other_total = max(0.0, wall_time_s - iter_total_s)

        ksp_list = [x["ksp_s"] for x in per_iter]
        assemble_list = [x["assemble_s"] for x in per_iter]
        ct_list = [x["ct_s"] for x in per_iter]

        def _stats(xs):
            if not xs:
                return {"mean": 0.0, "p50": 0.0, "p90": 0.0, "max": 0.0}
            arr = np.asarray(xs, dtype=float)
            return {
                "mean": float(arr.mean()),
                "p50": float(np.percentile(arr, 50)),
                "p90": float(np.percentile(arr, 90)),
                "max": float(arr.max()),
            }

        newton_per_step = self._parse_log_newton_iters(op_log_path) if op_log_path else []
        n_steps = len(newton_per_step)
        n_newton_total = int(sum(newton_per_step)) if newton_per_step else n_iter_total
        newton_stats = _stats(newton_per_step) if newton_per_step else {"mean": 0.0, "p50": 0.0, "p90": 0.0, "max": 0.0}

        total_accounted = ksp_total + ct_total + assemble_total + other_total
        frac = lambda x: (x / wall_time_s) if wall_time_s > 0 else 0.0

        payload = {
            "schema_version": 1,
            "case_name": case_name,
            "output_folder": output_folder,
            "config": {k: (v if not isinstance(v, dict) else dict(v)) for k, v in cfg.items()},
            "wall_time_s": wall_time_s,
            "n_linear_solves": n_iter_total,
            "n_steps": n_steps,
            "n_newton_total": n_newton_total,
            "newton_per_step": newton_per_step,
            "per_iter": per_iter,
            "summary": {
                "ksp_s_total": ksp_total,
                "assemble_s_total": assemble_total,
                "ct_s_total": ct_total,
                "other_s_total": other_total,
                "ksp_stats_s": _stats(ksp_list),
                "assemble_stats_s": _stats(assemble_list),
                "ct_stats_s": _stats(ct_list),
                "newton_per_step_stats": newton_stats,
                "fraction_ksp": frac(ksp_total),
                "fraction_assemble": frac(assemble_total),
                "fraction_ct": frac(ct_total),
                "fraction_other": frac(other_total),
                "total_accounted_s": total_accounted,
            },
        }

        out_path = os.path.join(output_folder, "profiling.json")
        if MPI.COMM_WORLD.rank == 0:
            with open(out_path, "w") as f:
                json.dump(payload, f, indent=2)
            print(f"[PROFILING] Wrote {out_path}")
        return payload


# ══════════════════════════════════════════════════════════════════════════════
# MATERIAL ASSEMBLY (with overrides)
# ══════════════════════════════════════════════════════════════════════════════

def build_spring_and_creep(cfg, mom_eq, E0, nu0):
    """Scenario-branched creep assembly with override hooks.

    Returns the list of non-elastic elements to add in the LEACHING phase.
    (Desai is added later in the OPERATION phase only, see build_desai.)

    Overrides (from cfg["MATERIAL_OVERRIDES"]):
        md_alpha_w, md_beta_w, md_delta   — MD parameter scalars (scalar replaces per-elem tensor)
        disable_pressure_solution          — skip PressureSolutionCreep
        disable_kelvin                     — skip Viscoelastic (SIC only)
        disable_desai                      — consumed by build_desai
    """
    overrides = cfg.get("MATERIAL_OVERRIDES", {})
    scenario = cfg["MATERIAL_SCENARIO"]
    use_md = cfg["USE_MUNSON_DAWSON"]
    n = mom_eq.n_elems
    sec_per_year = 365.25 * 24 * 3600

    elems_to_add = []

    if scenario == "A":
        if not use_md:
            if not overrides.get("disable_kelvin", False):
                eta = 2.5e5 * ut.GPa * to.ones(n)
                E1 = 42.0 * ut.GPa * to.ones(n)
                nu1 = 0.32 * to.ones(n)
                elems_to_add.append(sf.Viscoelastic(eta, E1, nu1, "kelvin"))
            ndc = 4.6
            A_dc = (40.0 * (1e-6) ** ndc / sec_per_year) * to.ones(n, dtype=to.float64)
            Q_dc = (6495.0 * 8.32) * to.ones(n)
            n_dc = ndc * to.ones(n)
            elems_to_add.append(sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation"))
        else:
            nmd = 4.99
            A_md = (18.31 * (1e-6) ** nmd / sec_per_year) * to.ones(n, dtype=to.float64)
            Q_md = (6356.0 * 8.32) * to.ones(n)
            n_md = nmd * to.ones(n)
            K0_md = 7.0e-7 * to.ones(n)
            c_md = 9.02e-3 * to.ones(n)
            m_md = 3.0 * to.ones(n)
            aw_md = overrides.get("md_alpha_w", -13.2) * to.ones(n)
            bw_md = overrides.get("md_beta_w", -7.738) * to.ones(n)
            d_md = overrides.get("md_delta", 0.58) * to.ones(n)
            mu_md = E0 / (2.0 * (1.0 + nu0))
            elems_to_add.append(sf.MunsonDawsonCreep(
                A=A_md, Q=Q_md, n=n_md, K0=K0_md, c=c_md, m=m_md,
                alpha_w=aw_md, beta_w=bw_md, delta=d_md, mu=mu_md,
                name="munson_dawson"))
    elif scenario == "B":
        if not use_md:
            if not overrides.get("disable_kelvin", False):
                eta = 1.0e15 * to.ones(n)
                E1 = 1.0e11 * to.ones(n)
                nu1 = 0.25 * to.ones(n)
                elems_to_add.append(sf.Viscoelastic(eta, E1, nu1, "kelvin"))
            ndc = 5.6897
            A_dc = (25.92 * (1e-6) ** ndc / sec_per_year) * to.ones(n, dtype=to.float64)
            Q_dc = (6495.0 * 8.32) * to.ones(n)
            n_dc = ndc * to.ones(n)
            elems_to_add.append(sf.DislocationCreep(A_dc, Q_dc, n_dc, "creep_dislocation"))
        else:
            nmd = 5.6897
            A_md = (17.28 * (1e-6) ** nmd / sec_per_year) * to.ones(n, dtype=to.float64)
            Q_md = (6495.0 * 8.32) * to.ones(n)
            n_md = nmd * to.ones(n)
            K0_md = 2253.87 * to.ones(n)
            c_md = 9.02e-3 * to.ones(n)
            m_md = 2.466 * to.ones(n)
            aw_md = overrides.get("md_alpha_w", 179.70) * to.ones(n)
            bw_md = overrides.get("md_beta_w", 60.00) * to.ones(n)
            d_md = overrides.get("md_delta", 299.95) * to.ones(n)
            mu_md = E0 / (2.0 * (1.0 + nu0))
            elems_to_add.append(sf.MunsonDawsonCreep(
                A=A_md, Q=Q_md, n=n_md, K0=K0_md, c=c_md, m=m_md,
                alpha_w=aw_md, beta_w=bw_md, delta=d_md, mu=mu_md,
                name="munson_dawson"))
    else:
        raise ValueError(f"MATERIAL_SCENARIO '{scenario}' not supported in sensitivity runner "
                         f"(only 'A' and 'B')")

    # PressureSolutionCreep — added unless disabled
    if not overrides.get("disable_pressure_solution", False):
        A_ps = (14176.0 * 1e-9 / 1e6 / sec_per_year) * to.ones(n)
        d_ps = 5.25e-3 * to.ones(n)
        Q_ps = (3252.0 * 8.32) * to.ones(n)
        elems_to_add.append(sf.PressureSolutionCreep(A_ps, d_ps, Q_ps, "creep_pressure"))

    return elems_to_add


def build_desai(cfg, mom_eq):
    """Desai viscoplastic element for operation phase (SIC model only).

    Returns the element, or None if USE_MUNSON_DAWSON or disable_desai.
    """
    overrides = cfg.get("MATERIAL_OVERRIDES", {})
    if cfg["USE_MUNSON_DAWSON"]:
        return None
    if overrides.get("disable_desai", False):
        return None

    scenario = cfg["MATERIAL_SCENARIO"]
    n = mom_eq.n_elems
    if scenario == "A":
        mu_1 = 6.89e-12 * to.ones(n)
        N_1 = 3.0 * to.ones(n)
        a_1 = 1.80e-5 * to.ones(n)
        eta_vp = 0.82 * to.ones(n)
        alpha_0 = 2.0e-3 * to.ones(n)
    elif scenario == "B":
        mu_1 = 1.016e-15 * to.ones(n)
        N_1 = 4.2515 * to.ones(n)
        a_1 = 1.101e-6 * to.ones(n)
        eta_vp = 1.7902 * to.ones(n)
        alpha_0 = 1.781e-3 * to.ones(n)
    else:
        raise ValueError(f"MATERIAL_SCENARIO '{scenario}' not supported")
    n_desai = 3.0 * to.ones(n)
    beta_1 = 0.0048 * to.ones(n)
    beta = 0.995 * to.ones(n)
    m_desai = -0.5 * to.ones(n)
    gamma = 0.095 * to.ones(n)
    sigma_t = 5.0 * to.ones(n)
    desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta_vp, n_desai, beta_1, beta,
                                 m_desai, gamma, sigma_t, alpha_0, "desai")
    return desai


# ══════════════════════════════════════════════════════════════════════════════
# PRESSURE SCHEDULE (with AMPLITUDE_SCALE)
# ══════════════════════════════════════════════════════════════════════════════

def build_op_pressure_schedule(cfg, tc_cycling, p_leach_end_mpa, p_leach_end, p_gas_mpa):
    """Build the operation-phase pressure schedule, applying AMPLITUDE_SCALE.

    AMPLITUDE_SCALE multiplies the pressure-swing amplitude around the scenario's
    reference pressure (leach_end for leaching cases, P_MEAN/P_EQUILIBRIUM for
    equilibrium cases).
    """
    amp_scale = float(cfg["AMPLITUDE_SCALE"])

    if cfg["PRESSURE_SCENARIO"] == "industry":
        if cfg["USE_LEACHING"]:
            p_mean = (p_leach_end_mpa + cfg["P_AMPLITUDE_MPA"]) * ut.MPa
        else:
            p_mean = cfg["P_MEAN_MPA"] * ut.MPa
        p_ampl = cfg["P_AMPLITUDE_MPA"] * amp_scale * ut.MPa

        return build_sinus_schedule_multi(
            tc_cycling, p_mean=p_mean, p_ampl=p_ampl,
            days=cfg["OPERATION_DAYS"], mode=cfg["SCHEDULE_MODE"],
            daily_period_hours=24.0, total_cycles=cfg["N_CYCLES"],
            clamp_min=None, clamp_max=None,
        )

    elif cfg["PRESSURE_SCENARIO"] == "transport":
        if cfg["USE_LEACHING"]:
            p_ref = p_leach_end_mpa
        else:
            p_ref = cfg["P_MEAN_MPA"]
        # scale offsets around p_ref by amp_scale
        p_high = p_ref + cfg["P_HIGH_OFFSET_MPA"] * amp_scale
        p_low = p_ref + cfg["P_LOW_OFFSET_MPA"] * amp_scale
        base_times_h = [0.0, 8.0, 12.0, 28.0, 32.0, 48.0]
        base_pressures_MPa = [p_high, p_high, p_low, p_low, p_high, p_high]
        return build_linear_schedule_multi(
            tc_cycling, base_times_h, base_pressures_MPa,
            days=cfg["OPERATION_DAYS"], mode=cfg["SCHEDULE_MODE"],
            resample_at_dt=True, total_cycles=cfg["N_CYCLES"],
        )

    elif cfg["PRESSURE_SCENARIO"] == "power_generation":
        if cfg["USE_LEACHING"]:
            p_base_pa = (p_leach_end_mpa + cfg["P_BASE_OFFSET_MPA"] * amp_scale) * ut.MPa
        else:
            p_base_pa = (cfg["P_EQUILIBRIUM_MPA"] + cfg["P_BASE_OFFSET_MPA"] * amp_scale) * ut.MPa
        return build_power_generation_schedule(
            tc_cycling, p_base_pa=p_base_pa,
            n_events=cfg["N_EVENTS"],
            operation_days=cfg["OPERATION_DAYS"],
            recovery_tau_hours=cfg["RECOVERY_TAU_HOURS"],
            p_min_pa=cfg["P_MIN_MPA"] * ut.MPa,
        )

    elif cfg["PRESSURE_SCENARIO"] == "csv":
        t_pressure, p_pressure = build_csv_pressure_schedule(
            tc_cycling, csv_file=cfg["CSV_FILE_PATH"],
            days=cfg["OPERATION_DAYS"], mode=cfg["SCHEDULE_MODE"],
            total_cycles=cfg["N_CYCLES"],
            rescale=cfg["RESCALE_PRESSURE"],
            rescale_min=cfg["RESCALE_MIN_MPA"],
            rescale_max=cfg["RESCALE_MAX_MPA"],
            resample_at_dt=True,
        )
        if not cfg["RESCALE_PRESSURE"] and cfg["USE_LEACHING"] and cfg["CSV_SHIFT_TO_LEACH_END"]:
            p_array = np.array(p_pressure)
            csv_min_pa = p_array.min()
            shift_pa = p_leach_end - csv_min_pa
            p_pressure = [p + shift_pa for p in p_pressure]
        if cfg["RAMP_UP_HOURS"] <= 0:
            apply_startup_ramp(
                t_pressure, p_pressure,
                p_start_pa=p_leach_end if cfg["USE_LEACHING"] else p_gas_mpa * ut.MPa,
                ramp_hours=cfg["RAMP_HOURS"], dt_hours=cfg["dt_hours"],
            )
        return t_pressure, p_pressure
    else:
        raise ValueError(f"Unknown PRESSURE_SCENARIO: {cfg['PRESSURE_SCENARIO']}")


# ══════════════════════════════════════════════════════════════════════════════
# CASE ORCHESTRATOR — mirrors Run.py main()
# ══════════════════════════════════════════════════════════════════════════════

def run_single_case(name):
    """Execute one sensitivity case. Returns profiling payload on success."""
    cfg = resolve_case(name)
    meta = resolve_cavern_metadata(cfg)
    output_folder = build_output_folder(name, cfg, meta)

    is_rank0 = MPI.COMM_WORLD.rank == 0

    if not FORCE_RERUN and is_case_complete(output_folder):
        if is_rank0:
            print(f"[SKIP] Case '{name}' already complete — {output_folder}")
        return None

    if is_rank0:
        print("\n" + "═" * 78)
        print(f"[CASE] {name}")
        print(f"[CASE] Output folder: {output_folder}")
        print("═" * 78)

    salt_density = 2200
    g = -9.81
    p_lithostatic = compute_lithostatic_pressure(meta["z_center"], meta["p_ref_mpa"], salt_density, g)
    p_lithostatic_mpa = p_lithostatic / ut.MPa

    # Derived pressures
    if cfg["USE_LEACHING"]:
        p_leach_end_mpa = cfg["LEACHING_END_FRACTION"] * p_lithostatic_mpa
        p_leach_end = p_leach_end_mpa * ut.MPa
        p_gas_mpa = None
        p_gas = None
    else:
        if cfg["PRESSURE_SCENARIO"] in ("industry", "transport"):
            p_gas_mpa = cfg["P_MEAN_MPA"]
        else:
            p_gas_mpa = cfg["P_EQUILIBRIUM_MPA"]
        p_gas = p_gas_mpa * ut.MPa
        p_leach_end_mpa = None
        p_leach_end = None

    # Load grid
    grid_path = os.path.join("..", "..", "..", "..", "grids", meta["grid_folder"])
    grid = sf.GridHandlerGMSH("geom", grid_path)
    z_max = meta["z_max"]
    p_ref = meta["p_ref_mpa"] * ut.MPa
    side_burden = p_ref
    over_burden = p_ref

    # Momentum equation
    mom_eq = LinearMomentumMod(grid, theta=0.5)
    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("cg")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-10, max_it=100)
    mom_eq.set_solver(mom_solver)

    # Material
    mat = sf.Material(mom_eq.n_elems)
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    E0 = 20.425 * ut.GPa * to.ones(mom_eq.n_elems)
    nu0 = 0.25 * to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")
    mat.add_to_elastic(spring_0)

    for elem in build_spring_and_creep(cfg, mom_eq, E0, nu0):
        mat.add_to_non_elastic(elem)

    mom_eq.set_material(mat)
    g_vec = [0.0, 0.0, g]
    mom_eq.build_body_force(g_vec)

    T0_field = 298 * to.ones(mom_eq.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)
    gas_density = 0.089

    # ─── INITIALIZATION PHASE (leaching or equilibrium) ──────────────────────
    if cfg["USE_LEACHING"]:
        tc_init = sf.TimeController(
            dt=cfg["LEACHING_DT_HOURS"], initial_time=0.0,
            final_time=cfg["LEACHING_DAYS"] * 24.0, time_unit="hour",
        )
        t_init, p_init = build_leaching_pressure_schedule(
            tc_init, p_start_pa=p_lithostatic, p_end_pa=p_leach_end,
            mode=cfg["LEACHING_MODE"], n_steps=cfg["STEPPED_N_STEPS"],
        )
        init_phase_name = "leaching"
        save_interval_init = max(1, int(72 / cfg["LEACHING_DT_HOURS"]))
    else:
        tc_init = sf.TimeController(
            dt=cfg["EQUILIBRIUM_DT_HOURS"], initial_time=0.0,
            final_time=cfg["EQUILIBRIUM_HOURS"], time_unit="hour",
        )
        t_init = [0.0, tc_init.t_final]
        p_init = [p_gas, p_gas]
        init_phase_name = "equilibrium"
        save_interval_init = 1

    bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_init.t_final])
    bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_init.t_final])
    bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_init.t_final])
    bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0,
                               [side_burden, side_burden], [0.0, tc_init.t_final], g=g_vec[2])
    bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0,
                                [side_burden, side_burden], [0.0, tc_init.t_final], g=g_vec[2])
    bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                              [over_burden, over_burden], [0.0, tc_init.t_final], g=g_vec[2])
    bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, z_max,
                                 p_init, t_init, g=g_vec[2])
    bc_init = momBC.BcHandler(mom_eq)
    for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
        bc_init.add_boundary_condition(bc)
    mom_eq.set_boundary_conditions(bc_init)

    output_folder_init = os.path.join(output_folder, init_phase_name)
    if is_rank0:
        print(f"[{init_phase_name.upper()}] Output: {output_folder_init}")

    output_mom_init = SparseSaveFields(mom_eq, interval=save_interval_init)
    output_mom_init.set_output_folder(output_folder_init)
    output_mom_init.add_output_field("u", "Displacement (m)")
    output_mom_init.add_output_field("eps_tot", "Total strain (-)")
    output_mom_init.add_output_field("sig", "Stress (Pa)")
    output_mom_init.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom_init.add_output_field("q_elems", "Von Mises stress (Pa)")
    os.makedirs(output_folder, exist_ok=True)

    # Install profiler BEFORE running (covers both leaching and operation)
    profiler = ProfilingRecorder()
    profiler.install(mom_eq)

    t_wall_start = time.perf_counter()
    try:
        sim_init = sf.Simulator_M(mom_eq, tc_init, [output_mom_init], True)
        sim_init.run()
        if is_rank0:
            print(f"[{init_phase_name.upper()}] Complete.")

        # ─── OPERATION PHASE ──────────────────────────────────────────────────
        desai = build_desai(cfg, mom_eq)
        if desai is not None:
            stress_to = ut.numpy2torch(mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3)))
            desai.compute_initial_hardening(stress_to, Fvp_0=0.0)
            mat.add_to_non_elastic(desai)
            mom_eq.set_material(mat)
            mom_eq.expect_vp_state = True
        else:
            mom_eq.expect_vp_state = False

        tc_cycling = sf.TimeController(
            dt=cfg["dt_hours"], initial_time=0.0,
            final_time=cfg["OPERATION_DAYS"] * 24.0, time_unit="hour",
        )
        t_pressure, p_pressure = build_op_pressure_schedule(
            cfg, tc_cycling, p_leach_end_mpa, p_leach_end, p_gas_mpa,
        )

        # Fade-in + debrining
        if cfg["RAMP_UP_HOURS"] > 0:
            p_fade_start = p_leach_end if cfg["USE_LEACHING"] else p_gas
            apply_fade_in(t_pressure, p_pressure, p_start_pa=p_fade_start,
                          fade_in_hours=cfg["RAMP_UP_HOURS"])
        extra_hours = 0.0
        if cfg["USE_LEACHING"] and cfg["DEBRINING_DAYS"] > 0:
            extra_hours = cfg["DEBRINING_DAYS"] * 24.0
            t_pressure, p_pressure = prepend_debrining(
                t_pressure, p_pressure, p_leach_end_pa=p_leach_end,
                debrining_days=cfg["DEBRINING_DAYS"],
            )

        total_operation_hours = cfg["OPERATION_DAYS"] * 24.0 + extra_hours

        if cfg["PRESSURE_SCENARIO"] in ("power_generation", "csv") and cfg["USE_VARIABLE_DT"]:
            t_arr = np.array(t_pressure, dtype=float)
            p_arr = np.array(p_pressure, dtype=float)
            p_of_t = lambda t, _ta=t_arr, _pa=p_arr: float(np.interp(t, _ta, _pa))
            time_list = build_time_list_by_dp_limit(
                total_operation_hours * ut.hour, p_of_t,
                dt_min_s=cfg["DT_FINE_HOURS"] * ut.hour,
                dt_max_s=cfg["DT_COARSE_HOURS"] * ut.hour,
                dp_max_pa=cfg["MAX_DP_MPA"] * ut.MPa,
            )
            tc_operation = TimeControllerFromList(time_list)
        else:
            tc_operation = sf.TimeController(
                dt=cfg["dt_hours"], initial_time=0.0,
                final_time=total_operation_hours, time_unit="hour",
            )

        # Operation BCs
        bc_west = momBC.DirichletBC("West", 0, [0.0, 0.0], [0.0, tc_operation.t_final])
        bc_bottom = momBC.DirichletBC("Bottom", 2, [0.0, 0.0], [0.0, tc_operation.t_final])
        bc_south = momBC.DirichletBC("South", 1, [0.0, 0.0], [0.0, tc_operation.t_final])
        bc_east = momBC.NeumannBC("East", 2, salt_density, 660.0,
                                   [side_burden, side_burden], [0.0, tc_operation.t_final], g=g_vec[2])
        bc_north = momBC.NeumannBC("North", 2, salt_density, 660.0,
                                    [side_burden, side_burden], [0.0, tc_operation.t_final], g=g_vec[2])
        bc_top = momBC.NeumannBC("Top", 2, 0.0, 0.0,
                                  [over_burden, over_burden], [0.0, tc_operation.t_final], g=g_vec[2])
        bc_cavern = momBC.NeumannBC("Cavern", 2, gas_density, z_max,
                                     p_pressure, t_pressure, g=g_vec[2])
        bc_operation = momBC.BcHandler(mom_eq)
        for bc in [bc_west, bc_bottom, bc_south, bc_east, bc_north, bc_top, bc_cavern]:
            bc_operation.add_boundary_condition(bc)
        mom_eq.set_boundary_conditions(bc_operation)

        output_folder_operation = os.path.join(output_folder, "operation")
        if is_rank0:
            print(f"[OPERATION] Output: {output_folder_operation}")

        # Save pressure schedule
        pressure_data = {
            "cavern_key": meta["cavern_key"],
            "cavern_label": meta["cavern_label"],
            "use_leaching": cfg["USE_LEACHING"],
            "debrining_days": cfg["DEBRINING_DAYS"] if cfg["USE_LEACHING"] else 0,
            "ramp_up_hours": cfg["RAMP_UP_HOURS"] if cfg["USE_LEACHING"] else 0,
            "pressure_scenario": cfg["PRESSURE_SCENARIO"],
            "mode": cfg["SCHEDULE_MODE"],
            "n_cycles": cfg["N_CYCLES"],
            "operation_days": cfg["OPERATION_DAYS"],
            "dt_hours": cfg["dt_hours"],
            "material_scenario": cfg["MATERIAL_SCENARIO"],
            "model": "munson_dawson" if cfg["USE_MUNSON_DAWSON"] else "safeincave",
            "amplitude_scale": cfg["AMPLITUDE_SCALE"],
            "p_lithostatic_mpa": p_lithostatic_mpa,
            "units": {"t_raw": "s", "p_raw": "Pa", "t": "hour", "p": "MPa"},
            "t_values_s": [float(t) for t in t_pressure],
            "p_values_Pa": [float(p) for p in p_pressure],
            "t_hours": [float(t / ut.hour) for t in t_pressure],
            "p_MPa": [float(p / ut.MPa) for p in p_pressure],
            "sensitivity_case": name,
            "material_overrides": cfg.get("MATERIAL_OVERRIDES", {}),
        }
        if is_rank0:
            with open(os.path.join(output_folder, "pressure_schedule.json"), "w") as f:
                json.dump(pressure_data, f, indent=2)

        _t_p = np.array(t_pressure, dtype=float)
        _p_p = np.array(p_pressure, dtype=float)
        p_interp_save = lambda t, _tp=_t_p, _pp=_p_p: float(np.interp(t, _tp, _pp))

        output_mom_op = SparseSaveFields(mom_eq, interval=15, p_interp=p_interp_save)
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

        if is_rank0:
            print("[OPERATION] Complete.")
    finally:
        profiler.uninstall()

    wall_time_s = time.perf_counter() - t_wall_start

    op_log_path = os.path.join(output_folder, "operation", "log.txt")
    payload = profiler.finalize(name, output_folder, cfg, wall_time_s, op_log_path=op_log_path)

    if is_rank0:
        summary = payload["summary"]
        print(f"[CASE] {name} complete: wall={wall_time_s:.1f}s  "
              f"ksp={summary['fraction_ksp']:.1%}  "
              f"assemble={summary['fraction_assemble']:.1%}  "
              f"ct={summary['fraction_ct']:.1%}  "
              f"newton_iters={payload['n_newton_total']}")

    return payload


# ══════════════════════════════════════════════════════════════════════════════
# SINGLE-CASE MAIN
# ══════════════════════════════════════════════════════════════════════════════

def _parse_cli():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", default=None,
                        help="Case name to run. Overrides CASE_TO_RUN at top of file.")
    parser.add_argument("--list", action="store_true",
                        help="Print all available case names and exit.")
    parser.add_argument("--force", action="store_true",
                        help="Rerun even if the output folder appears complete.")
    return parser.parse_args()


def main():
    global FORCE_RERUN
    args = _parse_cli()
    is_rank0 = MPI.COMM_WORLD.rank == 0

    if args.list:
        if is_rank0:
            print("Available sensitivity cases:")
            for name in SENSITIVITY_CASES:
                print(f"  - {name}")
        return

    if args.force:
        FORCE_RERUN = True

    case_name = args.case if args.case is not None else CASE_TO_RUN

    if case_name not in SENSITIVITY_CASES:
        if is_rank0:
            print(f"[ERROR] Unknown case '{case_name}'.")
            print("Available cases:")
            for name in SENSITIVITY_CASES:
                print(f"  - {name}")
        sys.exit(1)

    if is_rank0:
        print("\n" + "=" * 78)
        print("SENSITIVITY RUNNER — SINGLE CASE")
        print(f"  Case:        {case_name}")
        print(f"  FORCE_RERUN: {FORCE_RERUN}")
        print("=" * 78)

    try:
        payload = run_single_case(case_name)
    except Exception as exc:
        tb = traceback.format_exc()
        if is_rank0:
            print("\n" + "!" * 78)
            print(f"[ERROR] Case '{case_name}' FAILED: {exc}")
            print(tb)
            print("!" * 78)
            # Best-effort: write error.txt to the case folder
            try:
                cfg = resolve_case(case_name)
                meta = resolve_cavern_metadata(cfg)
                folder = build_output_folder(case_name, cfg, meta)
                os.makedirs(folder, exist_ok=True)
                with open(os.path.join(folder, "error.txt"), "w") as f:
                    f.write(f"Case: {case_name}\n{tb}")
            except Exception:
                pass
        sys.exit(2)

    if is_rank0:
        print("\n" + "=" * 78)
        if payload is None:
            print(f"[skipped]  {case_name}  (already complete; use --force to rerun)")
        else:
            s = payload["summary"]
            print(f"[ok]       {case_name:28s}  wall={payload['wall_time_s']:7.1f}s  "
                  f"ksp={s['fraction_ksp']:5.1%}  "
                  f"newton={payload['n_newton_total']:5d}")
        print("=" * 78)


if __name__ == "__main__":
    main()
