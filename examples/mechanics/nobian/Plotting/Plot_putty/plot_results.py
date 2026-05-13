"""
plot_results.py — Combined plotting script for salt cavern simulations.

Merges the functionality of:
  1. plot_convergence.py      — volume convergence + pressure schedule
  2. plot_stress_state.py     — p-q stress paths at 5 probes + dilatancy boundaries
  3. plot_FOS.py              — full-field FOS: global min + roof/mid/floor stats + XDMF
  4. plot_fracture_propagation.py — FOS at radial distances from wall + dilating zone depth

Usage:
  python plot_results.py
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
import meshio
import xml.etree.ElementTree as ET

import dolfinx as dfx
from mpi4py import MPI

import safeincave.PostProcessingTools as post
from case_index import detect_layout_and_collect_cases, filter_cases

# =============================================================================
# POSTER-FRIENDLY GLOBAL STYLE
# =============================================================================
plt.rcParams.update({
    'font.size':        18,
    'axes.titlesize':   22,
    'axes.labelsize':   20,
    'xtick.labelsize':  18,
    'ytick.labelsize':  18,
    'legend.fontsize':  18,
    'figure.titlesize': 24,
    'lines.linewidth':  2.0,
})

# =============================================================================
# USER CONFIGURATION
# =============================================================================

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", "Simulation", "output"))

# ── CASE SELECTION ─────────────────────────────────────────────────────────────
# SELECT: Filter which cases to plot (set to None to include all)
#
# Available filters:
#   "caverns"        - Cavern shapes to include: (e.g. "regular1200", "tilted1200","directcirculation600", "reversedcirculation1200")                                                   "fastleached1200", "tubefailure1200").
#   "pressure"       - Pressure scenario: "industry", "transport", "power_generation", "csv", or None
#   "scenario"       - Material scenario(s): "MD_A", "MD_B", "TUD2023_B" or None
#   "n_cycles"       - Number of cycles (int) or None
#   "operation_days" - Operation duration (int) or None
#   "case_contains"  - Substring match in case name or None

SELECT = {
    "caverns": ["spike_upper", "spike_upper_il2x", "spike_upper_il4x"],
    "pressure": ["industry", "transport", "power_generation"],
    "scenario": ["MD_B"],
    "n_cycles": None,
    "operation_days": 365,
    "case_contains": None,
}

# ── PLOT MODE ─────────────────────────────────────────────────────────────────
# PLOT_MODE controls how figures are grouped:
#
#   "compare_shapes"    - One figure per plot type, all cavern shapes overlaid.
#                         Use when running all shapes with the same settings
#                         (same pressure, scenario, cycles).
#                         Color = cavern shape, linestyle = scenario.
#
#   "compare_scenarios" - One figure PER CAVERN SHAPE, with different
#                         scenarios / pressures / cycles overlaid on it.
#                         Use when comparing what-if settings for the same cavern.
#                         Color = scenario, linestyle = scenario.
#
#   "compare_pressures" - One figure PER CAVERN SHAPE, with different
#                         pressure schemes overlaid on it.
#                         Use when comparing pressure schedules for the same
#                         cavern and scenario (e.g. industry vs power_generation).
#                         Color = pressure scheme, linestyle = pressure scheme.
#
#   "compare_sizes"     - One figure PER BASE SHAPE (e.g. "Regular"), with
#                         600,000 m³ and 1,200,000 m³ overlaid on the same axes.
#                         Color = scenario (TUD2023_A, MD_B, etc.).
#                         Linestyle = solid for 1,200,000 m³, dashed for 600,000 m³.
#                         Label = "TUD2023_B (1,200,000 m³)" etc.

PLOT_MODE = "compare_pressures"    # "compare_shapes", "compare_scenarios", "compare_pressures", or "compare_sizes"

FIGURES = {
    "convergence": False,           # Figure 1: volume convergence
    "stress_state": False,          # Figure 2: p-q stress paths
    "fos": False,                   # Figure 3: FOS over time
    "fracture_propagation": False, # Figure 4: dilatancy zone analysis
    "fos_summary": True,          # Figure 5: global min FOS + 4 pressure profiles
    "mc_failure": False,           # Figure 6: Mohr-Coulomb failure (interlayer cases)
}

# Stress state options
# Available: "ratigan_027", "ratigan_018", "spiers", "devries_comp", "devries_ext"
SHOW_DILATANCY = ["ratigan_027", "spiers", "devries_comp", "devries_ext"]

# FOS plot tuning (rolling quantile bands)
FOS_MAX_POINTS = 1200       # downsample for cleaner plots
FOS_BAND_WINDOW = 7        # rolling window size (odd)
FOS_SHOW_BAND = True        # show P10-P90 band + median
PROBE_ORDER = ["top", "quarter", "mid", "threequarter", "bottom"]

# Fracture propagation options
RADIAL_DISTANCES = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

# FOS (Factor of Safety) = q_boundary / q
#   FOS < 1.0 means dilating (unsafe)
#   FOS > 1.0 means safe

FOS_THRESHOLD = 1.0

# FOS options
WRITE_FOS_XDMF = False
XDMF_SUBDIR = os.path.join("operation", "_fos_outputs")
XDMF_FILENAME = "fos.xdmf"

OUT_DIR = os.path.join(ROOT, "_figures")
SHOW = False
DPI = 300
CAVERN_PHYS_TAG = 29

# Amplification factor for the initial-vs-final wall profile subplot in
# plot_convergence_per_cavern. Raise to exaggerate sub-percent radial changes
# so they are visible at print size; factor is annotated on the plot.
WALL_PROFILE_AMPLIFICATION = 100.0

# =============================================================================
# END OF USER CONFIGURATION
# =============================================================================


# =============================================================================
# 1. CONSTANTS & STYLE DICTS
# =============================================================================

MPA = 1e6
HOUR = 3600.0
DAY = 24.0 * HOUR

# Cavern colours picked to be visually distinct from the grayscale dilatancy
# palette below. Tilt moved off pink (which now belongs to De Vries 2005).
CAVERN_COLORS = {
    "Asymmetric":              "#ff7f0e",
    "Direct-circulation":      "#2ca02c",
    "IrregularFine":           "#9467bd",
    "Regular":                 "#1f77b4",
    "Reversed-circulation":    "#8c564b",
    "Tilt":                    "#17becf",   # cyan (was pink — clashed with De Vries)
    "Fast-leached":            "#d62728",
    "String-failure":          "#bcbd22",
    "Homogeneous":             "#1f77b4",
    # Reference (il4x, mesh-converged) keeps the primary red / green.
    # Coarser-interlayer variants use categorically distinct hues so that
    # an above-vs-above (or below-vs-below) mesh-comparison plot is readable.
    "Heterogeneous_above":                  "#d62728",  # red
    "Heterogeneous_above (2× coarser)":     "#ff7f0e",  # orange
    "Heterogeneous_above (4× coarser)":     "#e377c2",  # pink/magenta
    "Heterogeneous_below":                  "#2ca02c",  # green
    "Heterogeneous_below (2× coarser)":     "#17becf",  # cyan
    "Heterogeneous_below (4× coarser)":     "#9467bd",  # purple
    # Size-specific variants (solid = 1,200,000 m³, lighter = 600,000 m³)
    "Asymmetric (1,200,000 m³)":          "#ff7f0e",
    "Direct-circulation (1,200,000 m³)":  "#2ca02c",
    "IrregularFine (1,200,000 m³)":       "#9467bd",
    "Regular (1,200,000 m³)":             "#1f77b4",
    "Reversed-circulation (1,200,000 m³)":"#8c564b",
    "Tilt (1,200,000 m³)":                "#17becf",
    "Fast-leached (1,200,000 m³)":        "#d62728",
    "String-failure (1,200,000 m³)":      "#bcbd22",
    "Asymmetric (600,000 m³)":            "#ffbb78",
    "Direct-circulation (600,000 m³)":    "#98df8a",
    "IrregularFine (600,000 m³)":         "#c5b0d5",
    "Regular (600,000 m³)":               "#aec7e8",
    "Reversed-circulation (600,000 m³)":  "#c49c94",
    "Tilt (600,000 m³)":                  "#9edae5",
    "Fast-leached (600,000 m³)":          "#ff9896",
    "String-failure (600,000 m³)":        "#dbdb8d",
}

SCENARIO_COLORS = {
    # Current: <model>_<scenario>. TUD2023 = SafeInCave constitutive model.
    "TUD2023_A":            "#1f77b4",   # SafeInCave, Scenario A (blue)
    "MD_A":                 "#ff7f0e",   # Munson-Dawson, Scenario A (orange)
    "TUD2023_B":            "#2ca02c",   # SafeInCave, Scenario B (green)
    "MD_B":                 "#d62728",   # Munson-Dawson, Scenario B (red)
    "TUD2023_B_freecalibr": "#9467bd",   # SafeInCave, B_freecalibr (purple)
    "MD_B_freecalibr":      "#17becf",   # Munson-Dawson, B_freecalibr (cyan)
    # Legacy spelling (scenario-first) kept so old folders still match
    "A_SIC":            "#1f77b4",
    "A_MD":             "#ff7f0e",
    "B_SIC":            "#2ca02c",
    "B_MD":             "#d62728",
    "B_freecalibr_SIC": "#9467bd",
    "B_freecalibr_MD":  "#17becf",
    # Legacy names (kept for backward compatibility)
    "disloc_old_only":  "#1f77b4",
    "disloc_new_only":  "#ff7f0e",
    "desai_only":       "#2ca02c",
    "full_minus_desai": "#2ca02c",
    "full":             "#d62728",
    "full_minus_ps":    "#d62728",
    "md_only":          "#1f77b4",
    "md_steady_only":   "#bcbd22",
    "full_md":          "#e377c2",
    "interlayer":       "#7f7f7f",
    "nointerlayer":     "#8c564b",
    None:               "#333333",
}

SCENARIO_LINESTYLES = {
    # Current: <model>_<scenario>
    "TUD2023_A":            "-",
    "MD_A":                 "--",
    "TUD2023_B":            "-",
    "MD_B":                 "--",
    "TUD2023_B_freecalibr": "-",
    "MD_B_freecalibr":      "--",
    # Legacy spelling
    "A_SIC":            "-",
    "A_MD":             "--",
    "B_SIC":            "-",
    "B_MD":             "--",
    "B_freecalibr_SIC": "-",
    "B_freecalibr_MD":  "--",
    # Legacy
    "disloc_old_only":  "-",
    "disloc_new_only":  "--",
    "desai_only":       "-.",
    "full_minus_desai": "-",
    "full":             "-",
    "full_minus_ps":    "-",
    "md_only":          "-",
    "md_steady_only":   "--",
    "full_md":          "-",
    "interlayer":       "-",
    "nointerlayer":     "--",
    None:               "-",
}

PRESSURE_COLORS = {
    "industry":         "#1f77b4",   # blue
    "transport":        "#ff7f0e",   # orange
    "power_generation": "#2ca02c",   # green
    "csv":              "#d62728",   # red
    # Legacy
    "sinus":            "#9467bd",   # purple
    "irregular":        "#8c564b",   # brown
    "linear":           "#e377c2",   # pink
    None:               "#333333",
}

PRESSURE_LINESTYLES = {
    "industry":         "-",
    "transport":        "--",
    "power_generation": "-.",
    "csv":              ":",
    # Legacy
    "sinus":            "-",
    "irregular":        "--",
    "linear":           "-.",
    None:               "-",
}

PRESSURE_LABELS = {
    "industry":         "Industry",
    "transport":        "Transport",
    "power_generation": "Power generation",
    "csv":              "CSV profile",
    "sinus":            "Sinusoidal",
    "irregular":        "Irregular",
    "linear":           "Linear",
}

PROBE_COLORS = {
    "top":          "#e41a1c",
    "quarter":      "#377eb8",
    "mid":          "#4daf4a",
    "threequarter": "#984ea3",
    "bottom":       "#ff7f00",
}

PROBE_LINESTYLES = {
    "top":          "-",
    "quarter":      "-",
    "mid":          ":",
    "threequarter": "-.",
    "bottom":       (0, (5, 1)),
}

PROBE_LINEWIDTHS = {
    "top":          2.2,
    "quarter":      2.0,
    "mid":          2.6,
    "threequarter": 2.0,
    "bottom":       2.2,
}

CAVERN_ORDER = [
    "Asymmetric", "Direct-circulation", "Regular", "Reversed-circulation",
    "Tilt", "Fast-leached", "String-failure", "IrregularFine",
    "Homogeneous",
    "Heterogeneous_above",
    "Heterogeneous_above (2× coarser)",
    "Heterogeneous_above (4× coarser)",
    "Heterogeneous_below",
    "Heterogeneous_below (2× coarser)",
    "Heterogeneous_below (4× coarser)",
    "Asymmetric (1,200,000 m³)", "Direct-circulation (1,200,000 m³)", "Regular (1,200,000 m³)",
    "Reversed-circulation (1,200,000 m³)", "Tilt (1,200,000 m³)", "Fast-leached (1,200,000 m³)",
    "String-failure (1,200,000 m³)", "IrregularFine (1,200,000 m³)",
    "Asymmetric (600,000 m³)", "Direct-circulation (600,000 m³)", "Regular (600,000 m³)",
    "Reversed-circulation (600,000 m³)", "Tilt (600,000 m³)", "Fast-leached (600,000 m³)",
    "String-failure (600,000 m³)", "IrregularFine (600,000 m³)",
]
SCENARIO_ORDER = [
    # Current (model_scenario)
    "TUD2023_A", "MD_A", "TUD2023_B", "MD_B", "TUD2023_B_freecalibr", "MD_B_freecalibr",
    # Legacy (scenario_model)
    "A_SIC", "A_MD", "B_SIC", "B_MD", "B_freecalibr_SIC", "B_freecalibr_MD",
    # Legacy presets
    "disloc_old_only", "disloc_new_only", "desai_only", "full_minus_desai",
    "full", "full_minus_ps", "md_only", "md_steady_only", "full_md",
    "interlayer", "nointerlayer", None,
]


# =============================================================================
# 2. PATH HELPERS
# =============================================================================

def path_u_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.xdmf")

def path_u_h5(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.h5")

def path_geom_msh(case_folder):
    return os.path.join(case_folder, "operation", "mesh", "geom.msh")

def path_p_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "p_elems", "p_elems.xdmf")

def path_q_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "q_elems", "q_elems.xdmf")

def path_sig_xdmf(case_folder):
    candidates = [
        os.path.join(case_folder, "operation", "sig", "sig.xdmf"),
        os.path.join(case_folder, "operation", "sig.xdmf"),
    ]
    for p in candidates:
        if os.path.isfile(p):
            return p
    return candidates[0]  # default path even if missing

def path_pressure_json(case_folder):
    return os.path.join(case_folder, "pressure_schedule.json")

def is_valid_hdf5(h5_path: str) -> bool:
    try:
        with open(h5_path, "rb") as f:
            return f.read(8) == b"\x89HDF\r\n\x1a\n"
    except Exception:
        return False

def read_pressure_schedule(case_folder: str):
    pjson = path_pressure_json(case_folder)
    if not os.path.isfile(pjson):
        return None, None
    with open(pjson, "r") as f:
        data = json.load(f)
    if "t_hours" in data and "p_MPa" in data:
        t = np.asarray(data["t_hours"], float)
        p = np.asarray(data["p_MPa"], float)
        return t, p
    if "t_values_s" in data and "p_values_Pa" in data:
        t = np.asarray(data["t_values_s"], float) / HOUR
        p = np.asarray(data["p_values_Pa"], float) / MPA
        return t, p
    if "t_values" in data and "p_values" in data:
        t = np.asarray(data["t_values"], float) / HOUR
        p = np.asarray(data["p_values"], float) / MPA
        return t, p
    return None, None


# =============================================================================
# 3. CASE HELPERS
# =============================================================================

SIZE_LINESTYLES = {"1200k": "-", "600k": "--"}


def _extract_model_only(scenario_preset):
    """Return just the constitutive-model name from a preset key.

    'TUD2023_B' -> 'TUD2023', 'MD_B_freecalibr' -> 'MD'. Legacy forms where the
    parameter scenario comes first ('A_SIC', 'B_MD') are mapped to the new
    model name ('TUD2023', 'MD') so figures stay consistent across runs.
    """
    if not scenario_preset:
        return ""
    s = str(scenario_preset)
    head = s.split("_", 1)[0].upper()
    if head == "SIC":
        return "TUD2023"
    if head in ("A", "B"):
        # legacy "A_SIC" / "B_MD": the model is the second token
        parts = s.split("_")
        if len(parts) >= 2:
            tail = parts[-1].upper()
            if tail == "SIC":
                return "TUD2023"
            if tail == "MD":
                return "MD"
        return head  # fallback
    return head


# ──────────────────────────────────────────────────────────────────────────
# Depth conversion: mesh z-coordinate (0 at bottom, 660 at top) → actual
# depth-from-surface (m), based on the overburden above each cavern size.
# 0.6 Mm³ caverns: 200 m sand + 625 m salt above model → model top at 825 m.
# 1.2 Mm³ caverns: 200 m sand + 708 m salt above model → model top at 908 m.
# ──────────────────────────────────────────────────────────────────────────
DOMAIN_HEIGHT_M = 660.0
OVERBURDEN_TOP_DEPTH_M = {"1200k": 908.0, "600k": 825.0}


def z_mesh_to_depth_m(z_mesh, cavern_label):
    """Convert mesh z (0 = bottom, 660 = top) to depth from surface (m)."""
    tag = _extract_size_tag(cavern_label) or "1200k"
    top_depth = OVERBURDEN_TOP_DEPTH_M.get(tag, 908.0)
    return top_depth + (DOMAIN_HEIGHT_M - np.asarray(z_mesh, dtype=float))


def add_panel_label(ax, letter, fontsize=20):
    """Add a bold (A)/(B)/... label to the top-left of a subplot."""
    ax.text(0.02, 0.98, f"({letter})", transform=ax.transAxes,
            fontsize=fontsize, fontweight='bold', va='top', ha='left',
            bbox=dict(facecolor='white', alpha=0.85, edgecolor='none', pad=3),
            zorder=100)


def _extract_size_tag(cavern_label):
    """Extract '600k' or '1200k' from a cavern label.

    Handles both the current 'Regular (1,200,000 m³)' format and the legacy
    'Regular (1200k)' short form.
    """
    label = cavern_label or ""
    if "1,200,000" in label or "(1200k)" in label:
        return "1200k"
    if "600,000" in label or "(600k)" in label:
        return "600k"
    return None


def _strip_size(cavern_label):
    """Strip volume suffix: 'Regular (1,200,000 m³)' -> 'Regular' (also handles legacy '1200k')."""
    import re
    s = cavern_label or ""
    s = re.sub(r"\s*\([\d,]+\s*m³\)\s*$", "", s)
    s = re.sub(r"\s*\(\d+k\)\s*$", "", s)
    return s.strip()


def get_case_color_and_style(cavern_label, scenario_preset, pressure_scenario=None, mode=None):
    """Return (color, linestyle) based on the active PLOT_MODE."""
    if mode is None:
        mode = PLOT_MODE.lower().replace(" ", "_")
    if mode == "compare_pressures":
        color = PRESSURE_COLORS.get(pressure_scenario, "#333333")
        linestyle = PRESSURE_LINESTYLES.get(pressure_scenario, "-")
    elif mode == "compare_sizes":
        color = SCENARIO_COLORS.get(scenario_preset, "#333333")
        size_tag = _extract_size_tag(cavern_label)
        linestyle = SIZE_LINESTYLES.get(size_tag, "-")
    elif mode == "compare_scenarios":
        if scenario_preset is not None:
            color = SCENARIO_COLORS.get(scenario_preset, "#333333")
        else:
            color = CAVERN_COLORS.get(cavern_label, "#333333")
        linestyle = SCENARIO_LINESTYLES.get(scenario_preset, "-")
    else:  # compare_shapes
        color = CAVERN_COLORS.get(cavern_label, "#333333")
        linestyle = SCENARIO_LINESTYLES.get(scenario_preset, "-")
    return color, linestyle


def get_case_label(meta, mode=None):
    """Return a human-readable label for a case based on the active PLOT_MODE.

    Rules (thesis convention):
      - Never include the parameter scenario (A / B / B_freecalibr) — that goes
        in the figure caption.
      - compare_scenarios mode: show the constitutive model only (TUD2023 / MD).
      - All other modes: show the cavern shape / volume only.
      - Pressure scenario still shown in compare_pressures (distinguishes the
        lines when the whole point of the figure IS comparing pressures).
    """
    if mode is None:
        mode = PLOT_MODE.lower().replace(" ", "_")
    cav = meta.get("cavern_label", "")
    sc = meta.get("scenario_preset")
    ps = meta.get("pressure_scenario")
    if mode == "compare_pressures":
        return PRESSURE_LABELS.get(ps, ps or "unknown")
    elif mode == "compare_sizes":
        return cav
    elif mode == "compare_scenarios":
        return _extract_model_only(sc) or meta.get("case_name", "")
    else:  # compare_shapes
        return _strip_size(cav)


def get_series_key(meta):
    """Return a unique series key for a case (used in stress state dicts)."""
    return (meta.get("cavern_label"), meta.get("scenario_preset"), meta.get("pressure_scenario"))


def pick_one_case_per_series(cases_meta):
    """In compare_shapes mode, keep one case per (cavern, scenario, pressure) combo."""
    out = {}
    for m in cases_meta:
        key = (m.get("cavern_label"), m.get("scenario_preset"), m.get("pressure_scenario"))
        if key not in out:
            out[key] = m
    return list(out.values())


def group_cases_by_cavern(cases):
    """Group cases by cavern_label. Returns dict: cavern_label -> [cases]."""
    groups = {}
    for c in cases:
        cav = c.get("cavern_label", "unknown")
        groups.setdefault(cav, []).append(c)
    return groups


def group_cases_by_base_shape(cases):
    """Group cases by base shape (without volume suffix).
    Returns dict: base_shape -> [cases], e.g. 'Regular' -> [600k cases, 1200k cases]."""
    groups = {}
    for c in cases:
        base = _strip_size(c.get("cavern_label", "unknown"))
        groups.setdefault(base, []).append(c)
    return groups


def case_has_convergence_files(case_path):
    return (
        os.path.isfile(path_u_xdmf(case_path)) and
        os.path.isfile(path_u_h5(case_path)) and
        os.path.isfile(path_geom_msh(case_path)) and
        is_valid_hdf5(path_u_h5(case_path))
    )

def case_has_stress_files(case_path):
    return all(os.path.isfile(p) for p in [
        path_u_xdmf(case_path),
        path_geom_msh(case_path),
        path_p_xdmf(case_path),
        path_q_xdmf(case_path),
    ])

def case_has_fos_files(case_path):
    sig_ok = os.path.isfile(path_sig_xdmf(case_path))
    return sig_ok and all(os.path.isfile(p) for p in [
        path_geom_msh(case_path),
        path_p_xdmf(case_path),
        path_q_xdmf(case_path),
    ])

def case_has_fracture_files(case_path):
    return all(os.path.isfile(p) for p in [
        path_u_xdmf(case_path),
        path_geom_msh(case_path),
        path_sig_xdmf(case_path),
        path_p_xdmf(case_path),
        path_q_xdmf(case_path),
    ])


# =============================================================================
# 4. GEOMETRY HELPERS
# =============================================================================

def get_wall_indices_from_msh(msh_path):
    msh = meshio.read(msh_path)
    wall_idx = None

    if hasattr(msh, "cells_dict") and "line" in msh.cells_dict:
        wall_idx = np.unique(np.asarray(msh.cells_dict["line"]).reshape(-1))

    if wall_idx is None and isinstance(getattr(msh, "cells", None), dict):
        if "line" in msh.cells:
            wall_idx = np.unique(np.asarray(msh.cells["line"]).reshape(-1))

    if wall_idx is None:
        for cb in msh.cells:
            if getattr(cb, "type", None) == "line":
                wall_idx = np.unique(np.asarray(cb.data).reshape(-1))
                break

    if wall_idx is None:
        raise ValueError(f"No 'line' cells found in {msh_path}")

    return msh.points, wall_idx


def load_wall_points(case_folder):
    u_xdmf = path_u_xdmf(case_folder)
    msh_path = path_geom_msh(case_folder)

    points, _, _ = post.read_node_vector(u_xdmf)
    points_msh, wall_idx_msh = get_wall_indices_from_msh(msh_path)

    mapping = post.build_mapping(points_msh, points)
    wall_idx = np.array([mapping[i] for i in wall_idx_msh], dtype=int)

    wall_points = points[wall_idx]
    order = np.argsort(wall_points[:, 2])
    return wall_points[order]


def load_wall_points_initial_final(case_folder):
    """Return (initial, final) wall profile points, each sorted by z.

    Initial = mesh wall nodes at operation t=0 (post-leaching); final = same nodes
    after applying the displacement field at the last saved timestep.
    """
    u_xdmf = path_u_xdmf(case_folder)
    msh_path = path_geom_msh(case_folder)

    points, _, u_field = post.read_node_vector(u_xdmf)
    points_msh, wall_idx_msh = get_wall_indices_from_msh(msh_path)

    mapping = post.build_mapping(points_msh, points)
    wall_idx = np.array([mapping[i] for i in wall_idx_msh], dtype=int)

    base = points[wall_idx]
    u0 = u_field[0, wall_idx, :]
    uT = u_field[-1, wall_idx, :]

    initial = base + u0
    final = base + uT

    order = np.argsort(initial[:, 2])
    return initial[order], final[order]


def auto_generate_probes_from_wall_points(wall_points_sorted_z):
    pts = wall_points_sorted_z
    z = pts[:, 2]
    z_min, z_max = z.min(), z.max()

    fractions = {"bottom": 0.0, "threequarter": 0.25, "mid": 0.5, "quarter": 0.75, "top": 1.0}

    probes = {}
    for name, frac in fractions.items():
        z_target = z_min + frac * (z_max - z_min)
        idx = int(np.argmin(np.abs(z - z_target)))
        probes[name] = pts[idx]

    return probes


def auto_generate_probes_with_index(wall_points_sorted_z):
    """Like auto_generate_probes but also returns wall_idx (for fracture propagation)."""
    pts = wall_points_sorted_z
    z = pts[:, 2]
    z_min, z_max = z.min(), z.max()

    fractions = {"bottom": 0.0, "threequarter": 0.25, "mid": 0.5, "quarter": 0.75, "top": 1.0}

    probes = {}
    for name, frac in fractions.items():
        z_target = z_min + frac * (z_max - z_min)
        idx = int(np.argmin(np.abs(z - z_target)))
        probes[name] = {"wall_idx": idx, "wall_point": pts[idx]}

    return probes


def compute_outward_normal_2d(wall_points):
    n_pts = len(wall_points)
    normals = np.zeros((n_pts, 3))

    for i in range(n_pts):
        if i == 0:
            tangent = wall_points[1] - wall_points[0]
        elif i == n_pts - 1:
            tangent = wall_points[-1] - wall_points[-2]
        else:
            tangent = wall_points[i+1] - wall_points[i-1]

        tangent[1] = 0
        t_len = np.sqrt(tangent[0]**2 + tangent[2]**2)
        if t_len > 1e-10:
            tangent = tangent / t_len

        normal = np.array([tangent[2], 0.0, -tangent[0]])

        if normal[0] < 0:
            normal = -normal

        normals[i] = normal

    return normals


# =============================================================================
# 5. FIGURE 1 FUNCTIONS — Convergence
# =============================================================================

def extract_cavern_facets_from_msh_dolfinx(geom_msh_path: str, cavern_tag: int = 29):
    mesh, cell_tags, facet_tags = dfx.io.gmshio.read_from_msh(geom_msh_path, MPI.COMM_WORLD, 0)
    if facet_tags is None:
        raise RuntimeError("No facet tags were read from the .msh")
    fdim = mesh.topology.dim - 1
    cavern_facets = facet_tags.find(cavern_tag)
    if cavern_facets.size == 0:
        raise RuntimeError(f"No facets found with physical tag {cavern_tag}")
    mesh.topology.create_connectivity(fdim, 0)
    f2v = mesh.topology.connectivity(fdim, 0)
    facets = []
    for f in cavern_facets:
        verts = f2v.links(int(f))
        if len(verts) == 3:
            facets.append(verts)
    facets = np.asarray(facets, dtype=int)
    X = mesh.geometry.x.copy()
    return X, facets


def _area_vectors(points_xyz: np.ndarray, tri: np.ndarray):
    p0 = points_xyz[tri[:, 0]]
    p1 = points_xyz[tri[:, 1]]
    p2 = points_xyz[tri[:, 2]]
    return 0.5 * np.cross(p1 - p0, p2 - p0)


def orient_area_vectors_outward(points: np.ndarray, tris: np.ndarray, area_vecs: np.ndarray):
    center = np.mean(points[np.unique(tris.reshape(-1))], axis=0)
    tri_centroids = (points[tris[:, 0]] + points[tris[:, 1]] + points[tris[:, 2]]) / 3.0
    v = tri_centroids - center
    s = np.einsum("ij,ij->i", area_vecs, v)
    flip = s < 0.0
    area_vecs2 = area_vecs.copy()
    area_vecs2[flip] *= -1.0
    return area_vecs2


def compute_convergence_3d_percent(case_path: str, cavern_phys_tag: int = 29):
    geom_msh = path_geom_msh(case_path)
    u_xdmf = path_u_xdmf(case_path)

    pts_msh, tris_msh = extract_cavern_facets_from_msh_dolfinx(geom_msh, cavern_tag=cavern_phys_tag)

    area_vecs = _area_vectors(pts_msh, tris_msh)
    area_vecs = orient_area_vectors_outward(pts_msh, tris_msh, area_vecs)

    tri_centroids = (pts_msh[tris_msh[:, 0]] + pts_msh[tris_msh[:, 1]] + pts_msh[tris_msh[:, 2]]) / 3.0
    V0 = (1.0 / 3.0) * np.sum(np.einsum("ij,ij->i", tri_centroids, area_vecs))
    V0 = float(abs(V0))
    if not np.isfinite(V0) or V0 <= 0.0:
        raise RuntimeError(f"Invalid V0={V0}")

    points_xdmf, time_list, u_field = post.read_node_vector(u_xdmf)

    mapping = post.build_mapping(pts_msh, points_xdmf)
    tri_xdmf = np.vectorize(mapping.__getitem__)(tris_msh).astype(int)
    v0, v1, v2 = tri_xdmf[:, 0], tri_xdmf[:, 1], tri_xdmf[:, 2]

    u_ref = u_field[0]
    Nt = u_field.shape[0]
    dV = np.zeros(Nt, dtype=float)
    for k in range(Nt):
        u0 = u_field[k, v0, :] - u_ref[v0, :]
        u1 = u_field[k, v1, :] - u_ref[v1, :]
        u2 = u_field[k, v2, :] - u_ref[v2, :]
        uc = (u0 + u1 + u2) / 3.0
        dV[k] = np.sum(np.einsum("ij,ij->i", uc, area_vecs))

    conv_pct = -100.0 * (dV / V0)
    t_days = np.asarray(time_list, float) / DAY
    conv_pct[0] = 0.0
    return t_days, conv_pct


# =============================================================================
# 6. FIGURE 2 FUNCTIONS — Stress state
# =============================================================================

def read_stress_paths(case_folder, probes_dict):
    p_path = path_p_xdmf(case_folder)
    q_path = path_q_xdmf(case_folder)

    points_p, time_list, p_elems = post.read_cell_scalar(p_path)
    _, time_list2, q_elems = post.read_cell_scalar(q_path)

    n = min(len(time_list), len(time_list2))
    time_list = time_list[:n]
    p_elems = p_elems[:n]
    q_elems = q_elems[:n]

    out = {}
    for key, probe_xyz in probes_dict.items():
        idx = post.find_closest_point(probe_xyz, points_p)
        p = -p_elems[:, idx] / MPA
        q =  q_elems[:, idx] / MPA
        out[key] = (p, q)

    return out


def plot_dilatancy_boundaries(ax, show_boundaries=None, p_min=0.01, p_max=40.0, npts=500):
    if show_boundaries is None:
        show_boundaries = ["ratigan_027", "ratigan_018", "spiers", "devries_comp", "devries_ext"]

    p = np.linspace(p_min, p_max, npts)
    I1 = 3.0 * p

    def q_from_sqrtJ2(sqrtJ2):
        return np.sqrt(3.0) * sqrtJ2

    # Dilatancy boundaries use a dedicated grayscale palette so they never
    # collide with cavern or scenario colours in the legend. Line-style alone
    # distinguishes individual boundary criteria.
    boundary_styles = {
        "ratigan_027":  {"color": "#000000", "linestyle": "--", "linewidth": 1.6, "alpha": 0.95},
        "ratigan_018":  {"color": "#000000", "linestyle": ":",  "linewidth": 1.6, "alpha": 0.95},
        "spiers":       {"color": "#555555", "linestyle": "-.", "linewidth": 1.6, "alpha": 0.95},
        "devries_comp": {"color": "#888888", "linestyle": "-",  "linewidth": 1.8, "alpha": 1.00},
        "devries_ext":  {"color": "#888888", "linestyle": "--", "linewidth": 1.8, "alpha": 1.00},
        "mc_anhydrite": {"color": "#333333", "linestyle": "-",  "linewidth": 1.8, "alpha": 0.95},
        "mc_mudstone":  {"color": "#333333", "linestyle": "--", "linewidth": 1.6, "alpha": 0.90},
    }

    if "ratigan_027" in show_boundaries:
        D = 0.27
        sqrtJ2 = D * I1
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), label="Ratigan 1991 (D=0.27)",
                **boundary_styles["ratigan_027"])

    if "ratigan_018" in show_boundaries:
        D = 0.18
        sqrtJ2 = D * I1
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), label="Ratigan 1991 (D=0.18)",
                **boundary_styles["ratigan_018"])

    if "spiers" in show_boundaries:
        D_sp, b_sp = 0.27, 1.9
        sqrtJ2 = D_sp * I1 + b_sp
        ax.plot(p, q_from_sqrtJ2(sqrtJ2), label="Spiers 1988",
                **boundary_styles["spiers"])

    def devries_q(I1_MPa, psi_rad, D1=0.683, D2=0.512, T0=1.50, m=0.75, sigma_ref=1.0):
        sgn = np.sign(I1_MPa)
        sgn[sgn == 0.0] = 1.0
        denom = (np.sqrt(3.0) * np.cos(psi_rad) - D2 * np.sin(psi_rad))
        sqrtJ2_ = D1 * ((I1_MPa / (sgn * sigma_ref)) ** m) / denom + T0
        return np.sqrt(3.0) * sqrtJ2_

    psi_comp = -np.pi / 6.0
    psi_ext = np.pi / 6.0

    if "devries_comp" in show_boundaries:
        ax.plot(p, devries_q(I1, psi_comp), label="De Vries 2005 (comp)",
                **boundary_styles["devries_comp"])

    if "devries_ext" in show_boundaries:
        ax.plot(p, devries_q(I1, psi_ext), label="De Vries 2005 (ext)",
                **boundary_styles["devries_ext"])

    if "mc_anhydrite" in show_boundaries:
        ax.plot(p, q_dil_mc(p, c_MPa=4.0, phi_deg=35.0),
                label=r"Mohr-Coulomb anhydrite ($c$=4 MPa, $\varphi$=35°)",
                **boundary_styles["mc_anhydrite"])

    if "mc_mudstone" in show_boundaries:
        ax.plot(p, q_dil_mc(p, c_MPa=2.0, phi_deg=25.0),
                label=r"Mohr-Coulomb mudstone ($c$=2 MPa, $\varphi$=25°)",
                **boundary_styles["mc_mudstone"])


# =============================================================================
# 6b. FOS PLOT HELPERS (rolling bands, downsampling)
# =============================================================================

def _downsample_xy(x, y, max_points=1200):
    """Thin x, y arrays to at most max_points evenly-spaced samples."""
    n = len(x)
    if n <= max_points:
        return np.asarray(x), np.asarray(y)
    idx = np.linspace(0, n - 1, max_points).astype(int)
    idx[0] = 0
    idx[-1] = n - 1
    return np.asarray(x)[idx], np.asarray(y)[idx]


def _rolling_quantile_band(y, window=31, qlo=0.0, qhi=1.0):
    """Rolling quantile band + median for oscillatory series.

    Default qlo=0.0, qhi=1.0 makes the shaded band span the rolling minimum to
    rolling maximum (100% of values per window), so cycle troughs and peaks are
    fully visible — essential for FoS, where the trough is the safety-relevant
    instant.
    """
    y = np.asarray(y)
    n = len(y)
    lo = np.empty(n)
    hi = np.empty(n)
    med = np.empty(n)
    half = window // 2
    for i in range(n):
        a = max(0, i - half)
        b = min(n, i + half + 1)
        w = y[a:b]
        lo[i] = np.quantile(w, qlo)
        hi[i] = np.quantile(w, qhi)
        med[i] = np.quantile(w, 0.5)
    return lo, med, hi


def _compute_fos_probes_for_case(c):
    """Compute FOS time series at 5 wall probes + global field minimum.

    Returns (t_days, dict[probe_name -> FOS_array]).
    The dict includes a special key "global_min" with the per-timestep
    minimum FOS across ALL cells in the mesh.
    """
    folder = c["case_path"]
    p_path = path_p_xdmf(folder)
    q_path = path_q_xdmf(folder)
    sig_path = path_sig_xdmf(folder)

    centroids_p, time_list_p, p_elems = post.read_cell_scalar(p_path)
    _, time_list_q, q_elems = post.read_cell_scalar(q_path)
    centroids_sig, time_list_sig, sig33 = post.read_cell_tensor(sig_path)

    n_times = min(len(time_list_p), len(time_list_q), len(time_list_sig))
    time_days = np.array(time_list_p[:n_times]) / DAY
    p_elems = p_elems[:n_times]
    q_elems = q_elems[:n_times]
    sig33 = sig33[:n_times]

    p_MPa = -p_elems / MPA
    q_MPa = q_elems / MPA
    sig33_MPa = -sig33 / MPA

    # --- Global field FOS (all cells) ---
    FOS_field = compute_FOS_field(p_elems, q_elems, sig33, compression_positive=True, q_tol_MPa=1e-3)
    fos_global_min = np.nanmin(np.where(np.isfinite(FOS_field), FOS_field, np.nan), axis=1)

    # --- Probe-based FOS ---
    wall_points = load_wall_points(folder)
    probes = auto_generate_probes_from_wall_points(wall_points)

    fos_by_probe = {}
    for probe_name, probe_xyz in probes.items():
        idx_p = post.find_closest_point(probe_xyz, centroids_p)
        idx_sig = post.find_closest_point(probe_xyz, centroids_sig)

        p_t = p_MPa[:, idx_p]
        q_t = q_MPa[:, idx_p]
        sig_point6 = tensor33_to_voigt6(sig33_MPa[:, idx_sig, :, :])
        psi_t = _psi_from_voigt6(sig_point6)

        q_boundary = q_dil_devries(p_t, psi_t)
        q_safe = np.where(q_t < 1e-3, 1e-3, q_t)
        fos = np.where(q_t < 1e-3, 100.0, q_boundary / q_safe)
        fos_by_probe[probe_name] = fos

    fos_by_probe["global_min"] = fos_global_min

    return time_days, fos_by_probe


# =============================================================================
# 7. FIGURE 3 FUNCTIONS — FOS field
# =============================================================================

def load_p_q(case_folder):
    p_path = path_p_xdmf(case_folder)
    q_path = path_q_xdmf(case_folder)
    if not os.path.isfile(p_path):
        raise FileNotFoundError(f"Missing p_elems for {case_folder}")
    if not os.path.isfile(q_path):
        raise FileNotFoundError(f"Missing q_elems for {case_folder}")
    _, t_p, p_elems = post.read_cell_scalar(p_path)
    _, t_q, q_elems = post.read_cell_scalar(q_path)
    n = min(len(t_p), len(t_q))
    return np.asarray(t_p[:n], float), np.asarray(p_elems[:n], float), np.asarray(q_elems[:n], float)


def load_sig(case_folder):
    sig_path = path_sig_xdmf(case_folder)
    if not os.path.isfile(sig_path):
        raise FileNotFoundError(f"Missing sig for {case_folder}")
    _, t, sig_vals = post.read_cell_tensor(sig_path)
    return np.asarray(t, float), np.asarray(sig_vals, float), sig_path


def _psi_from_tensor33(sig_Pa, compression_positive=True):
    """Compute Lode angle from (nt, nc, 3, 3) stress tensor via eigenvalues."""
    sig = 0.5 * (sig_Pa + np.swapaxes(sig_Pa, -1, -2))
    sig_eff = -sig if compression_positive else sig
    vals = np.linalg.eigvalsh(sig_eff)
    s1, s2, s3 = vals[:, 2], vals[:, 1], vals[:, 0]
    I1 = s1 + s2 + s3
    mean = I1 / 3.0
    s1d, s2d, s3d = s1 - mean, s2 - mean, s3 - mean
    J2 = (1.0/6.0) * ((s1d - s2d)**2 + (s2d - s3d)**2 + (s3d - s1d)**2)
    J3 = s1d * s2d * s3d
    J2_safe = np.maximum(J2, 1e-30)
    x = (3.0*np.sqrt(3.0)/2.0) * (J3 / (J2_safe**1.5))
    x = np.clip(x, -1.0, 1.0)
    theta = (1.0/3.0) * np.arccos(x)
    psi = theta - np.pi/6.0
    return psi


def q_dil_rd(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    p = np.asarray(p_MPa, dtype=float)
    psi = np.asarray(psi, dtype=float)
    I1 = 3.0 * p
    denom = (np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi))
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)
    num = D1 * (np.abs(I1) / sigma_ref) ** m + T0
    sqrtJ2_dil = num / denom
    return np.sqrt(3.0) * sqrtJ2_dil


def compute_FOS_field(p_elems, q_elems, sig_vals, *, compression_positive=True, q_tol_MPa=1e-3):
    """Compute FOS for full field using stored p/q and psi from tensor."""
    nt, nc = p_elems.shape[0], p_elems.shape[1]
    FOS = np.full((nt, nc), np.inf, dtype=float)

    for it in range(nt):
        p_MPa = -p_elems[it] / MPA
        q_MPa = q_elems[it] / MPA

        psi = _psi_from_tensor33(sig_vals[it], compression_positive=compression_positive)

        q_dil = q_dil_rd(p_MPa, psi)

        mask = q_MPa >= float(q_tol_MPa)
        FOS[it, mask] = q_dil[mask] / q_MPa[mask]

    FOS[~np.isfinite(FOS)] = np.inf
    return np.clip(FOS, 0.0, 1e6)


def extract_cavern_wall_cells_slices(geom_msh_path, cavern_phys_tag=29, top_fraction=0.2, bottom_fraction=0.2):
    mesh, cell_tags, facet_tags = dfx.io.gmshio.read_from_msh(geom_msh_path, MPI.COMM_SELF, 0)
    if facet_tags is None:
        raise RuntimeError("No facet tags present in geom.msh")
    dim = mesh.topology.dim
    fdim = dim - 1
    facets = facet_tags.find(cavern_phys_tag)
    if facets.size == 0:
        raise RuntimeError(f"No facets with tag {cavern_phys_tag}")
    mesh.topology.create_connectivity(fdim, 0)
    f2v = mesh.topology.connectivity(fdim, 0)
    X = mesh.geometry.x
    zc = np.empty(facets.size, dtype=float)
    for i, f in enumerate(facets):
        vs = f2v.links(int(f))
        zc[i] = X[vs, 2].mean()
    z_min, z_max = float(zc.min()), float(zc.max())
    H = z_max - z_min
    z_top_thr = z_max - float(top_fraction) * H
    z_bot_thr = z_min + float(bottom_fraction) * H
    roof_mask  = zc >= z_top_thr
    floor_mask = zc <= z_bot_thr
    mid_mask   = ~(roof_mask | floor_mask)
    mesh.topology.create_connectivity(fdim, dim)
    f2c = mesh.topology.connectivity(fdim, dim)
    def cells_from_mask(mask):
        subset = facets[mask]
        cells = set()
        for f in subset:
            for c in f2c.links(int(f)):
                cells.add(int(c))
        return np.array(sorted(cells), dtype=int)
    return {"roof": cells_from_mask(roof_mask), "mid": cells_from_mask(mid_mask), "floor": cells_from_mask(floor_mask)}


def fos_series_stats(FOS, cell_idx):
    Nt = FOS.shape[0]
    out = {"min": np.full(Nt, np.nan), "mean": np.full(Nt, np.nan), "p05": np.full(Nt, np.nan)}
    if cell_idx.size == 0:
        return out
    vals = FOS[:, cell_idx].copy()
    vals[~np.isfinite(vals)] = np.nan
    out["min"]  = np.nanmin(vals, axis=1)
    out["mean"] = np.nanmean(vals, axis=1)
    out["p05"]  = np.nanpercentile(vals, 5.0, axis=1)
    return out


# XDMF mesh helpers

def _infer_xdmf_mesh_grid_name(xdmf_path: str) -> str:
    fallbacks = ["Grid", "mesh", "Mesh", "domain", "Domain"]
    try:
        tree = ET.parse(xdmf_path)
        root = tree.getroot()
        for grid in root.iter():
            if grid.tag.lower().endswith("grid"):
                name = grid.attrib.get("Name") or grid.attrib.get("name")
                if name:
                    return name
    except Exception:
        pass
    return fallbacks[0]


def _read_mesh_from_sig_xdmf(sig_xdmf_path: str):
    grid_name = _infer_xdmf_mesh_grid_name(sig_xdmf_path)
    candidates = [grid_name, "Grid", "mesh", "Mesh", "domain", "Domain"]
    tried = []
    with dfx.io.XDMFFile(MPI.COMM_SELF, sig_xdmf_path, "r") as xf:
        last_err = None
        for name in candidates:
            if name in tried:
                continue
            tried.append(name)
            try:
                mesh = xf.read_mesh(name=name)
                return mesh, name
            except Exception as e:
                last_err = e
                continue
    raise RuntimeError(
        f"Could not read mesh from {sig_xdmf_path}. Tried Grid names: {tried}. "
        f"Last error: {last_err}"
    )


def write_fos_xdmf_from_sig(case_folder: str, sig_xdmf_path: str, time_list: np.ndarray, FOS: np.ndarray):
    out_dir = os.path.join(case_folder, XDMF_SUBDIR)
    os.makedirs(out_dir, exist_ok=True)
    out_xdmf = os.path.join(out_dir, XDMF_FILENAME)

    mesh, used_name = _read_mesh_from_sig_xdmf(sig_xdmf_path)

    V0 = dfx.fem.functionspace(mesh, ("DG", 0))
    fos_fun = dfx.fem.Function(V0)
    fos_fun.name = "FoS"

    ncells = mesh.topology.index_map(mesh.topology.dim).size_local
    if FOS.shape[1] != ncells:
        raise RuntimeError(
            f"[FoS XDMF] Cell count mismatch for case={case_folder}\n"
            f"  FOS has ncells={FOS.shape[1]}\n"
            f"  mesh (from sig.xdmf grid '{used_name}') has ncells={ncells}\n"
            "This means you are NOT aligned with the mesh used to write sig.xdmf."
        )

    FOS_write = np.asarray(FOS, float).copy()
    FOS_write[~np.isfinite(FOS_write)] = np.nan

    with dfx.io.XDMFFile(MPI.COMM_SELF, out_xdmf, "w") as xw:
        xw.write_mesh(mesh)
        for it, t in enumerate(time_list):
            fos_fun.x.array[:] = FOS_write[it, :]
            xw.write_function(fos_fun, t=float(t))

    print("[SAVED XDMF]", out_xdmf, "(+ .h5)")


# =============================================================================
# 8. FIGURE 4 FUNCTIONS — Fracture propagation
# =============================================================================

def _psi_from_voigt6(sig_voigt):
    """Compute Lode angle from (..., 6) Voigt notation [s11,s22,s33,s12,s13,s23]."""
    s11 = sig_voigt[..., 0]
    s22 = sig_voigt[..., 1]
    s33 = sig_voigt[..., 2]
    s12 = sig_voigt[..., 3]
    s13 = sig_voigt[..., 4]
    s23 = sig_voigt[..., 5]

    p = (s11 + s22 + s33) / 3.0
    dev11 = s11 - p
    dev22 = s22 - p
    dev33 = s33 - p

    J2 = 0.5 * (dev11**2 + dev22**2 + dev33**2) + s12**2 + s13**2 + s23**2
    J2 = np.maximum(J2, 1e-30)

    J3 = (dev11 * (dev22 * dev33 - s23**2)
          - s12 * (s12 * dev33 - s23 * s13)
          + s13 * (s12 * s23 - dev22 * s13))

    sqrtJ2 = np.sqrt(J2)
    arg = (3.0 * np.sqrt(3.0) / 2.0) * J3 / (sqrtJ2**3)
    arg = np.clip(arg, -1.0, 1.0)
    psi = -(1.0 / 3.0) * np.arcsin(arg)
    return psi


def q_dil_devries(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    I1 = 3.0 * p_MPa
    denom = np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi)
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)
    num = D1 * (np.abs(I1) / sigma_ref) ** m + T0
    sqrtJ2_dil = num / denom
    return np.sqrt(3.0) * sqrtJ2_dil


def q_dil_mc(p_MPa, c_MPa=4.0, phi_deg=35.0):
    """Mohr-Coulomb failure boundary in p-q space (compression matching).

    Returns the von Mises stress q at the Drucker-Prager yield surface
    (matched to MC in triaxial compression) for given mean stress p.

    Parameters
    ----------
    p_MPa : array_like
        Mean stress in MPa (compression-positive).
    c_MPa : float
        Cohesion in MPa. Default: 4.0 (anhydrite, Cała et al. 2018).
    phi_deg : float
        Friction angle in degrees. Default: 35° (anhydrite).

    Returns
    -------
    q_MPa : ndarray
        Von Mises stress at the MC boundary (MPa).
    """
    phi = np.radians(phi_deg)
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    # Drucker-Prager in p-q space (compression matching):
    # q = 6*sin(phi)/(3-sin(phi)) * p + 6*c*cos(phi)/(3-sin(phi))
    q = 6.0 * sin_phi / (3.0 - sin_phi) * np.asarray(p_MPa) + 6.0 * c_MPa * cos_phi / (3.0 - sin_phi)
    return np.maximum(q, 0.0)


def compute_FOS_mc_field(p_elems, q_elems, c_MPa=4.0, phi_deg=35.0, q_tol_MPa=1e-3):
    """Compute Mohr-Coulomb safety factor for full field.

    FOS_MC = q_boundary_MC / q.  FOS < 1 means the MC yield criterion is exceeded.

    Parameters
    ----------
    p_elems, q_elems : ndarray, shape (nt, nc)
        Mean stress and von Mises stress in Pa (SafeInCave convention: p < 0 in compression).
    c_MPa : float
        Cohesion in MPa.
    phi_deg : float
        Friction angle in degrees.
    q_tol_MPa : float
        Minimum q threshold to avoid division by zero.

    Returns
    -------
    FOS : ndarray, shape (nt, nc)
    """
    nt, nc = p_elems.shape
    FOS = np.full((nt, nc), np.inf, dtype=float)
    for it in range(nt):
        p_MPa = -p_elems[it] / MPA   # compression-positive
        q_MPa = q_elems[it] / MPA
        q_boundary = q_dil_mc(p_MPa, c_MPa, phi_deg)
        mask = q_MPa >= q_tol_MPa
        FOS[it, mask] = q_boundary[mask] / q_MPa[mask]
    FOS[~np.isfinite(FOS)] = np.inf
    return np.clip(FOS, 0.0, 1e6)


def tensor33_to_voigt6(sig33):
    s11 = sig33[..., 0, 0]
    s22 = sig33[..., 1, 1]
    s33 = sig33[..., 2, 2]
    s12 = sig33[..., 0, 1]
    s13 = sig33[..., 0, 2]
    s23 = sig33[..., 1, 2]
    return np.stack([s11, s22, s33, s12, s13, s23], axis=-1)


def compute_FOS_point(p_mpa, q_mpa, psi=None, q_tol_MPa=1e-3):
    """Compute FOS at individual points (for fracture propagation)."""
    if psi is None:
        psi = -np.pi / 6.0
    q_boundary = q_dil_devries(p_mpa, psi)
    q_safe = np.where(q_mpa < q_tol_MPa, q_tol_MPa, q_mpa)
    FOS = q_boundary / q_safe
    FOS = np.where(q_mpa < q_tol_MPa, 100.0, FOS)
    return FOS


def generate_radial_sample_points(wall_points, probes, normals, radial_distances):
    sample_points = {}
    for probe_name, probe_info in probes.items():
        wall_idx = probe_info["wall_idx"]
        wall_pt = probe_info["wall_point"]
        normal = normals[wall_idx]
        points_at_distances = []
        for dist in radial_distances:
            pt = wall_pt + dist * normal
            points_at_distances.append((dist, pt))
        sample_points[probe_name] = points_at_distances
    return sample_points


def compute_global_dilation_penetration(case_folder, wall_points,
                                         region_names=None, region_fractions=None):
    """Volumetric scan of dilatancy zone, broken down by 5 wall regions.

    For every salt cell in the mesh, computes its perpendicular distance to the
    cavern wall (nearest wall node) and assigns it to one of 5 regions based on
    the z-fraction of that nearest wall node (top / quarter / mid / threequarter
    / bottom). At each timestep, the function returns, per region, the maximum
    distance-to-wall at which any cell has FOS < 1 — i.e. the actual depth the
    dilatancy zone has penetrated into the salt within that region.

    Returns
    -------
    time_days : ndarray (Nt,)
    region_depths : dict[name] -> ndarray (Nt,)
    """
    from scipy.spatial import cKDTree

    if region_names is None:
        region_names = ["bottom", "threequarter", "mid", "quarter", "top"]
    if region_fractions is None:
        region_fractions = [0.0, 0.25, 0.5, 0.75, 1.0]

    p_path = path_p_xdmf(case_folder)
    q_path = path_q_xdmf(case_folder)
    sig_path = path_sig_xdmf(case_folder)

    centroids_p, time_list_p, p_elems = post.read_cell_scalar(p_path)
    _, time_list_q, q_elems = post.read_cell_scalar(q_path)
    centroids_sig, time_list_sig, sig33 = post.read_cell_tensor(sig_path)

    n_times = min(len(time_list_p), len(time_list_q), len(time_list_sig))
    p_elems = p_elems[:n_times]
    q_elems = q_elems[:n_times]
    sig33 = sig33[:n_times]
    time_days = np.array(time_list_p[:n_times]) / DAY

    p_MPa = -p_elems / MPA
    q_MPa = q_elems / MPA
    sig33_MPa = -sig33 / MPA

    # Map sigma cell ordering to p/q ordering via centroid match (usually identical)
    if centroids_sig.shape[0] != centroids_p.shape[0] or \
       not np.allclose(centroids_sig, centroids_p, atol=1e-6):
        tree_c = cKDTree(centroids_sig)
        _, sig_to_p = tree_c.query(centroids_p)
        sig33_MPa = sig33_MPa[:, sig_to_p, :, :]

    sig_voigt = tensor33_to_voigt6(sig33_MPa)            # (Nt, nc, 6)
    psi = _psi_from_voigt6(sig_voigt)                    # (Nt, nc)
    q_boundary = q_dil_devries(p_MPa, psi)               # (Nt, nc)
    q_safe = np.where(q_MPa < 1e-3, 1e-3, q_MPa)
    fos = q_boundary / q_safe
    fos = np.where(q_MPa < 1e-3, 100.0, fos)
    fos_mask = fos < FOS_THRESHOLD                       # (Nt, nc)

    # Distance from each cell centroid to the nearest cavern-wall node
    tree = cKDTree(wall_points)
    dist_to_wall, nearest_idx = tree.query(centroids_p)  # (nc,), (nc,)

    # Region by z-fraction of the nearest wall node
    z_wall = wall_points[:, 2]
    z_min, z_max = float(z_wall.min()), float(z_wall.max())
    span = max(z_max - z_min, 1e-12)
    cell_z_frac = (wall_points[nearest_idx, 2] - z_min) / span

    region_idx = np.argmin(
        np.abs(cell_z_frac[:, None] - np.asarray(region_fractions)[None, :]),
        axis=1
    )

    # Contiguity parameters. The salt fine_size near the cavern wall is
    # ~4.5 m, so we use a tight gap threshold and only consider cells within
    # MAX_RADIAL_M of the wall. Cells beyond that distance are far-field and
    # not part of the dilatancy zone even if their FOS<1 (typically a sign
    # of outer-boundary numerical artifacts).
    GAP_THRESHOLD_M = 5.0
    NEAR_WALL_THRESH_M = 5.0
    MAX_RADIAL_M = 30.0

    Nt, nc = fos_mask.shape
    region_depths = {name: np.zeros(Nt) for name in region_names}
    for r, name in enumerate(region_names):
        sel = (region_idx == r) & (dist_to_wall <= MAX_RADIAL_M)
        if not sel.any():
            continue
        d_sel = dist_to_wall[sel]
        fos_sel = fos_mask[:, sel]                       # (Nt, n_in_region)
        for t in range(Nt):
            f = fos_sel[t]
            if not f.any():
                continue
            d_failing = np.sort(d_sel[f])
            # Dilation must reach close to the wall to count as connected
            if d_failing[0] > NEAR_WALL_THRESH_M:
                continue
            depth = float(d_failing[0])
            prev = depth
            for i in range(1, d_failing.size):
                if d_failing[i] - prev <= GAP_THRESHOLD_M:
                    depth = float(d_failing[i])
                    prev = depth
                else:
                    break
            region_depths[name][t] = depth

    return time_days, region_depths


def read_stress_at_points(case_folder, sample_points):
    p_path = path_p_xdmf(case_folder)
    q_path = path_q_xdmf(case_folder)

    centroids_p, time_list_p, p_elems = post.read_cell_scalar(p_path)
    centroids_q, time_list_q, q_elems = post.read_cell_scalar(q_path)

    sig_path = path_sig_xdmf(case_folder)
    centroids_sig, time_list_sig, sig33 = post.read_cell_tensor(sig_path)

    n_times = min(len(time_list_p), len(time_list_q), len(time_list_sig))
    time_list = time_list_p[:n_times]
    p_elems = p_elems[:n_times]
    q_elems = q_elems[:n_times]
    sig33 = sig33[:n_times]

    time_days = np.array(time_list) / DAY

    p_elems_MPa = -p_elems / MPA
    q_elems_MPa = q_elems / MPA
    sig33_MPa = -sig33 / MPA

    stress_data = {}
    for probe_name, dist_points in sample_points.items():
        stress_data[probe_name] = {}
        for dist, xyz in dist_points:
            idx_p = post.find_closest_point(xyz, centroids_p)
            idx_sig = post.find_closest_point(xyz, centroids_sig)

            p = p_elems_MPa[:, idx_p]
            q = q_elems_MPa[:, idx_p]

            sig_point33 = sig33_MPa[:, idx_sig, :, :]
            sig_point6 = tensor33_to_voigt6(sig_point33)
            psi = _psi_from_voigt6(sig_point6)

            stress_data[probe_name][dist] = (p, q, psi)

    return time_days, stress_data


def compute_FOS_data(time_days, stress_data):
    fos_data = {}
    for probe_name, dist_data in stress_data.items():
        fos_data[probe_name] = {}
        for dist, (p, q, psi) in dist_data.items():
            fos = compute_FOS_point(p, q, psi)
            fos_data[probe_name][dist] = fos
    return fos_data


# =============================================================================
# 9. FIGURE 1 PLOTTING — Convergence
# =============================================================================

def plot_convergence_combined(cases):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 9), sharex=True,
                                   gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12})

    single_scenario = len({c.get("scenario_preset") for c in cases}) == 1

    for c in cases:
        cav = c.get("cavern_label")
        sc = c.get("scenario_preset")
        ps = c.get("pressure_scenario")
        col, ls = get_case_color_and_style(cav, sc, ps)
        if single_scenario:
            ls = "-"
        label = get_case_label(c)
        try:
            t_days, conv = compute_convergence_3d_percent(c["case_path"], cavern_phys_tag=CAVERN_PHYS_TAG)
        except Exception as e:
            print(f"[SKIP] {cav}/{sc}: {e}")
            continue
        ax1.plot(t_days, conv, linewidth=2.0, color=col, linestyle=ls, alpha=0.95, label=label)

    ax1.set_ylabel("Convergence (ΔV/V₀) (%)")
    ax1.grid(True, alpha=0.3)

    h, l = ax1.get_legend_handles_labels()
    uniq = {}
    for hh, ll in zip(h, l):
        if ll not in uniq:
            uniq[ll] = hh
    if uniq:
        # Place the legend ABOVE the figure (not on the axes) so it never gets
        # clipped when there are many cases (e.g. Homogeneous vs Heterogeneous).
        fig.legend(uniq.values(), uniq.keys(), loc="upper center",
                   bbox_to_anchor=(0.5, 0.99),
                   ncol=min(4, len(uniq)), fontsize=14, frameon=True)

    # Plot pressure schedule(s) — overlay all distinct schedules in compare_pressures
    _plotted_pressures = set()
    for c in cases:
        ps = c.get("pressure_scenario")
        if ps in _plotted_pressures:
            continue
        _plotted_pressures.add(ps)
        tH, pMPa = read_pressure_schedule(c["case_path"])
        if tH is not None:
            ax2.plot(tH / 24.0, pMPa, linewidth=1.7, color='black',
                     label=PRESSURE_LABELS.get(ps, ps))
    if not _plotted_pressures:
        ax2.text(0.5, 0.5, "No pressure_schedule.json found.", ha="center", va="center", transform=ax2.transAxes)
    if len(_plotted_pressures) > 1:
        ax2.legend(fontsize=18, frameon=True, loc="upper right")
    ax2.set_ylabel("Pressure (MPa)")
    ax2.set_xlabel("Time (days)")
    ax2.grid(True, alpha=0.3)

    add_panel_label(ax1, "A")
    add_panel_label(ax2, "B")

    fig.tight_layout(rect=[0, 0, 1, 0.90])

    outname = f"convergence_combined_pressure={SELECT.get('pressure')}_scenario={SELECT.get('scenario')}.png"
    outpath = os.path.join(OUT_DIR, outname.replace(" ", ""))
    fig.savefig(outpath, dpi=DPI)
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close(fig)


def plot_convergence_separate(cases):
    for c in cases:
        cav = c.get("cavern_label")
        sc = c.get("scenario_preset")
        ps = c.get("pressure_scenario")
        col, ls = get_case_color_and_style(cav, sc, ps)
        label = get_case_label(c)

        try:
            t_days, conv = compute_convergence_3d_percent(c["case_path"], cavern_phys_tag=CAVERN_PHYS_TAG)
        except Exception as e:
            print(f"[SKIP] {cav}/{sc}: {e}")
            continue

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 10), sharex=True,
                                       gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12})

        ax1.plot(t_days, conv, linewidth=2.0, color=col, alpha=0.95, label=label)
        ax1.set_ylabel("Convergence (ΔV/V0) (%)")
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc="upper left", fontsize=18, frameon=True)

        ax1.text(0.985, 0.03, cav, transform=ax1.transAxes,
                 fontsize=18, fontweight='bold', ha='right', va='bottom',
                 bbox=dict(facecolor='white', alpha=0.85, edgecolor='gray', boxstyle='round,pad=0.4'))

        tH, pMPa = read_pressure_schedule(c["case_path"])
        if tH is None:
            ax2.text(0.5, 0.5, "No pressure_schedule.json found.", ha="center", va="center", transform=ax2.transAxes)
        else:
            ax2.plot(tH / 24.0, pMPa, linewidth=1.7, color='black')
        ax2.set_ylabel("Pressure (MPa)")
        ax2.set_xlabel("Time (days)")
        ax2.grid(True, alpha=0.3)

        add_panel_label(ax1, "A")
        add_panel_label(ax2, "B")

        safe_name = c.get("case_name", "unknown").replace(" ", "_")
        outname = f"convergence_{safe_name}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


def plot_convergence_per_cavern(cases, group_fn=None):
    """One convergence figure per cavern, scenarios/pressures overlaid.

    Layout: left column = convergence (top) + pressure (bottom); right column
    (spanning both rows) = wall profile, initial vs final. The wall panel uses
    non-equal aspect so sub-percent radial displacements remain visible.
    """
    groups = (group_fn or group_cases_by_cavern)(cases)
    for cav_label, cav_cases in groups.items():
        fig = plt.figure(figsize=(22, 10))
        gs = fig.add_gridspec(2, 2, width_ratios=[2.4, 1.0],
                              height_ratios=[2.2, 1.0], hspace=0.12, wspace=0.18)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
        ax3 = fig.add_subplot(gs[:, 1])

        for c in cav_cases:
            cav_full = c.get("cavern_label", cav_label)
            sc = c.get("scenario_preset")
            ps = c.get("pressure_scenario")
            col, ls = get_case_color_and_style(cav_full, sc, ps)
            label = get_case_label(c)
            try:
                t_days, conv = compute_convergence_3d_percent(c["case_path"], cavern_phys_tag=CAVERN_PHYS_TAG)
            except Exception as e:
                print(f"[SKIP] {cav_label}/{sc}: {e}")
                continue
            ax1.plot(t_days, conv, linewidth=2.0, color=col, linestyle=ls, alpha=0.95, label=label)

        ax1.set_ylabel("Convergence (ΔV/V0) (%)")
        ax1.grid(True, alpha=0.3)
        h, l = ax1.get_legend_handles_labels()
        uniq = {}
        for hh, ll in zip(h, l):
            if ll not in uniq:
                uniq[ll] = hh
        if uniq:
            ax1.legend(uniq.values(), uniq.keys(), loc="upper center",
                       bbox_to_anchor=(0.5, 1.0), ncol=min(3, len(uniq)),
                       fontsize=18, frameon=True)

        ax1.text(0.985, 0.03, cav_label, transform=ax1.transAxes,
                 fontsize=18, fontweight='bold', ha='right', va='bottom',
                 bbox=dict(facecolor='white', alpha=0.85, edgecolor='gray', boxstyle='round,pad=0.4'))

        # Plot pressure schedule(s) — overlay all distinct schedules.
        # Use a neutral colour (black) so it never clashes with a scenario line above.
        _plotted_pressures = set()
        for c in cav_cases:
            ps = c.get("pressure_scenario")
            if ps in _plotted_pressures:
                continue
            _plotted_pressures.add(ps)
            tH, pMPa = read_pressure_schedule(c["case_path"])
            if tH is not None:
                ax2.plot(tH / 24.0, pMPa, linewidth=1.7, color='black',
                         label=PRESSURE_LABELS.get(ps, ps))
        if not _plotted_pressures:
            ax2.text(0.5, 0.5, "No pressure_schedule.json found.", ha="center", va="center", transform=ax2.transAxes)
        if len(_plotted_pressures) > 1:
            ax2.legend(fontsize=18, frameon=True, loc="upper right")
        ax2.set_ylabel("Pressure (MPa)")
        ax2.set_xlabel("Time (days)")
        ax2.grid(True, alpha=0.3)

        # Wall profile: initial vs final shape. Displacement is amplified by
        # WALL_PROFILE_AMPLIFICATION so that <1 % radial changes are visible;
        # the factor is annotated on the subplot so readers aren't misled.
        # Y-axis is actual depth from the surface (m), growing downward.
        any_drawn = False
        initial_drawn = False
        for c in cav_cases:
            cav_full = c.get("cavern_label", cav_label)
            sc = c.get("scenario_preset")
            ps = c.get("pressure_scenario")
            col, ls = get_case_color_and_style(cav_full, sc, ps)
            label = get_case_label(c)
            try:
                wp_init, wp_final = load_wall_points_initial_final(c["case_path"])
            except Exception as e:
                print(f"[WALL-PROFILE SKIP] {cav_label}/{sc}/{ps}: {e}")
                continue
            u_delta = wp_final - wp_init
            wp_final_amp = wp_init + WALL_PROFILE_AMPLIFICATION * u_delta
            depth_init = z_mesh_to_depth_m(wp_init[:, 2], cav_full)
            depth_final = z_mesh_to_depth_m(wp_final_amp[:, 2], cav_full)
            if not initial_drawn:
                ax3.plot(wp_init[:, 0], depth_init, color='black', linewidth=2.0,
                         label='Initial')
                initial_drawn = True
            ax3.plot(wp_final_amp[:, 0], depth_final, color=col, linestyle=ls,
                     linewidth=1.8, alpha=0.9, label=f'Final — {label}')
            any_drawn = True

        if any_drawn:
            ax3.set_xlabel("x (m)")
            ax3.set_ylabel("Depth (m)")
            ax3.grid(True, alpha=0.3)
            ax3.invert_yaxis()  # depth grows downward
            ax3.legend(fontsize=13, frameon=True, loc='best')
            ax3.text(0.98, 0.02,
                     f"Displacement ×{WALL_PROFILE_AMPLIFICATION:g}",
                     transform=ax3.transAxes, fontsize=12, style='italic',
                     color='#444444', va='bottom', ha='right',
                     bbox=dict(facecolor='white', alpha=0.85, edgecolor='gray', boxstyle='round,pad=0.3'))
        else:
            ax3.text(0.5, 0.5, "No wall-profile data", ha="center", va="center",
                     transform=ax3.transAxes)
            ax3.axis("off")

        add_panel_label(ax1, "A")
        add_panel_label(ax2, "B")
        add_panel_label(ax3, "C")

        safe_cav = cav_label.replace(" ", "_")
        outname = f"convergence_{safe_cav}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


def plot_initial_vs_final_shapes(cases):
    """2x3 grid of cavern profiles: one cavern per subplot, initial wall in
    black and final (amplified) wall in the cavern's signature colour. The
    cavern order is fixed: top row = regular, direct-circulation,
    reversed-circulation; bottom row = tube-failure, tilted, fast-leached.
    Y-axis is depth from surface (m), with the axis inverted so the cavern is
    drawn the right way up.
    """
    SHAPE_ORDER = [
        ("regular",            "Regular"),
        ("directcirculation",  "Direct-circulation"),
        ("reversedcirculation","Reversed-circulation"),
        ("stringfailure",      "String-failure"),
        ("tilt",               "Tilt"),
        ("fastleached",        "Fast-leached"),
    ]

    def _norm(label):
        return _strip_size(label).lower().replace("-", "").replace(" ", "")

    by_base = {}
    for c in cases:
        base = _norm(c.get("cavern_label", ""))
        if base and base not in by_base:
            by_base[base] = c

    nrows, ncols = 2, 3
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.5 * ncols, 5.6 * nrows),
                             squeeze=False)

    for i, (base_key, display) in enumerate(SHAPE_ORDER):
        r, k = divmod(i, ncols)
        ax = axes[r][k]

        c = by_base.get(base_key)
        if c is None:
            ax.text(0.5, 0.5, f"No data\n({display})", ha="center",
                    va="center", transform=ax.transAxes, fontsize=12,
                    color="#666666")
            ax.set_title(display, fontsize=14)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        cav_full = c.get("cavern_label")
        sc = c.get("scenario_preset")
        ps = c.get("pressure_scenario")
        col, _ls = get_case_color_and_style(cav_full, sc, ps)

        try:
            wp_init, wp_final = load_wall_points_initial_final(c["case_path"])
        except Exception as e:
            print(f"[WALL-PROFILE SKIP] {cav_full}: {e}")
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(display, fontsize=14)
            continue

        u_delta = wp_final - wp_init
        wp_final_amp = wp_init + WALL_PROFILE_AMPLIFICATION * u_delta

        depth_init = z_mesh_to_depth_m(wp_init[:, 2], cav_full)
        depth_final = z_mesh_to_depth_m(wp_final_amp[:, 2], cav_full)

        ax.plot(wp_init[:, 0], depth_init, color='black',
                linewidth=1.8, label='Initial')
        ax.plot(wp_final_amp[:, 0], depth_final, color=col,
                linewidth=1.8, label='Final')

        ax.set_title(display, fontsize=14)
        ax.set_xlabel("x (m)")
        if k == 0:
            ax.set_ylabel("Depth (m)")
        ax.grid(True, alpha=0.3)
        ax.invert_yaxis()
        ax.legend(loc='best', fontsize=10, frameon=True)

    fig.text(0.99, 0.005,
             f"Final wall displacement ×{WALL_PROFILE_AMPLIFICATION:g}",
             fontsize=10, style='italic', color='#444444',
             ha='right', va='bottom',
             bbox=dict(facecolor='white', alpha=0.85, edgecolor='gray',
                       boxstyle='round,pad=0.3'))

    fig.tight_layout(rect=[0, 0.02, 1, 1])

    outname = (f"shapes_initial_final_pressure={SELECT.get('pressure')}"
               f"_scenario={SELECT.get('scenario')}.png")
    outpath = os.path.join(OUT_DIR, outname.replace(" ", ""))
    fig.savefig(outpath, dpi=DPI)
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close(fig)


# =============================================================================
# 10. FIGURE 2 PLOTTING — Stress state
# =============================================================================

STRESS_PROBE_TITLES = {
    "top":          "Cavern top",
    "quarter":      "Cavern quarter-point",
    "mid":          "Cavern mid-point",
    "threequarter": "Cavern three-quarter-point",
    "bottom":       "Cavern bottom",
}


def plot_stress_combined(cases, stress_by_series):
    probe_types = ["top", "quarter", "mid", "threequarter", "bottom"]

    fig, axes = plt.subplots(2, 3, figsize=(20, 11))
    axes = axes.flatten()

    panel_letters = ["A", "B", "C", "D", "E", "F"]

    for i, ptype in enumerate(probe_types):
        ax = axes[i]
        plot_dilatancy_boundaries(ax, show_boundaries=SHOW_DILATANCY)

        for (cav, sc, ps), d in stress_by_series.items():
            p, q = d[ptype]
            color, linestyle = get_case_color_and_style(cav, sc, ps)
            label = get_case_label({"cavern_label": cav, "scenario_preset": sc, "pressure_scenario": ps})

            ax.plot(p, q, linewidth=2.0, color=color, linestyle=linestyle, label=label)
            ax.scatter(p[-1], q[-1], s=30, edgecolors="black", linewidths=0.6,
                       color=color, zorder=5)

        ax.set_title(STRESS_PROBE_TITLES.get(ptype, ptype), fontsize=24, fontweight='bold', pad=6)
        ax.set_xlabel("Mean stress p (MPa)", fontsize=22)
        ax.set_ylabel("Differential stress q (MPa)", fontsize=22)
        ax.tick_params(axis='both', labelsize=20)
        ax.grid(True, alpha=0.3)
        add_panel_label(ax, panel_letters[i])

    axp = axes[5]
    _plotted_pressures = set()
    for c in cases:
        ps = c.get("pressure_scenario")
        if ps in _plotted_pressures:
            continue
        _plotted_pressures.add(ps)
        tH, pMPa = read_pressure_schedule(c["case_path"])
        if tH is not None:
            axp.plot(tH / 24.0, pMPa, linewidth=2.0, color='black',
                     label=PRESSURE_LABELS.get(ps, ps))
    if not _plotted_pressures:
        axp.text(0.5, 0.5, "No pressure_schedule.json found.",
                 ha="center", va="center", transform=axp.transAxes)
        axp.axis("off")
    else:
        axp.set_title("Pressure schedule", fontsize=24, fontweight='bold', pad=6)
        axp.set_xlabel("Time (days)", fontsize=22)
        axp.set_ylabel("Pressure (MPa)", fontsize=22)
        axp.tick_params(axis='both', labelsize=20)
        axp.grid(True, alpha=0.3)
        if len(_plotted_pressures) > 1:
            axp.legend(fontsize=18, frameon=True)
        add_panel_label(axp, panel_letters[5])

    # Build separate legend entries for dilatancy boundaries and cavern shapes
    boundary_handles, boundary_labels = [], []
    shape_handles, shape_labels = [], []
    handles, labels = axes[0].get_legend_handles_labels()
    seen = set()
    for h, l in zip(handles, labels):
        if l in seen:
            continue
        seen.add(l)
        # Dilatancy boundary labels contain known keywords
        if any(kw in l for kw in ["Ratigan", "Spiers", "De Vries"]):
            boundary_handles.append(h)
            boundary_labels.append(l)
        else:
            shape_handles.append(h)
            shape_labels.append(l)

    # Legend above the figure, with explicit section headers so readers can see
    # at a glance which entries are cavern geometries vs. dilatancy criteria.
    all_handles, all_labels = [], []
    if shape_handles:
        all_handles.append(plt.Line2D([], [], color="none"))
        all_labels.append(r"$\bf{Cavern\ geometries}$")
        all_handles.extend(shape_handles)
        all_labels.extend(shape_labels)
    if boundary_handles:
        all_handles.append(plt.Line2D([], [], color="none"))
        all_labels.append(r"$\bf{Dilatancy\ criteria}$")
        all_handles.extend(boundary_handles)
        all_labels.extend(boundary_labels)

    fig.legend(all_handles, all_labels, loc="upper center", bbox_to_anchor=(0.5, 0.99),
               ncol=min(4, len(all_labels)), frameon=True, fontsize=18)

    fig.tight_layout(rect=[0, 0, 1, 0.80])

    outname = f"pq_paths_combined_pressure={SELECT.get('pressure')}_scenario={SELECT.get('scenario')}.png"
    outpath = os.path.join(OUT_DIR, outname.replace(" ", "").replace("'", "").replace("[", "").replace("]", ""))
    fig.savefig(outpath, dpi=DPI)
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close(fig)


def plot_stress_separate(cases, stress_by_series):
    probe_types = ["top", "quarter", "mid", "threequarter", "bottom"]

    panel_letters = ["A", "B", "C", "D", "E", "F"]

    for (cav, sc, ps), d in stress_by_series.items():
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        axes = axes.flatten()

        case_label = get_case_label({"cavern_label": cav, "scenario_preset": sc, "pressure_scenario": ps})
        color, linestyle = get_case_color_and_style(cav, sc, ps)

        for i, ptype in enumerate(probe_types):
            ax = axes[i]
            plot_dilatancy_boundaries(ax, show_boundaries=SHOW_DILATANCY)

            p, q = d[ptype]
            ax.plot(p, q, linewidth=2.5, color=color, linestyle=linestyle, label=case_label)
            ax.scatter(p[-1], q[-1], s=40, edgecolors="black", linewidths=0.8,
                       color=color, zorder=5)

            ax.set_title(f"p-q stress path: {ptype}", fontsize=22, fontweight='bold')
            ax.set_xlabel("Mean stress p (MPa)", fontsize=22)
            ax.set_ylabel("Differential stress q (MPa)", fontsize=22)
            ax.tick_params(axis='both', labelsize=20)
            ax.grid(True, alpha=0.3)
            add_panel_label(ax, panel_letters[i])

        axp = axes[5]
        case_path = None
        for m in cases:
            if m.get("cavern_label") == cav and m.get("scenario_preset") == sc and m.get("pressure_scenario") == ps:
                case_path = m["case_path"]
                break

        if case_path:
            tH, pMPa = read_pressure_schedule(case_path)
            if tH is not None:
                axp.plot(tH / 24.0, pMPa, linewidth=2.0, color='black')
                axp.set_title("Pressure schedule", fontsize=22, fontweight='bold')
                axp.set_xlabel("Time (days)", fontsize=22)
                axp.set_ylabel("Pressure (MPa)", fontsize=22)
                axp.tick_params(axis='both', labelsize=20)
                axp.grid(True, alpha=0.3)
                add_panel_label(axp, panel_letters[5])
            else:
                axp.text(0.5, 0.5, "No pressure data", ha="center", va="center", transform=axp.transAxes)
                axp.axis("off")
        else:
            axp.axis("off")

        # Single legend above the figure, separating geometries from dilatancy criteria.
        b_handles, b_labels = [], []
        s_handles, s_labels = [], []
        handles, labels = axes[0].get_legend_handles_labels()
        seen = set()
        for h, l in zip(handles, labels):
            if l in seen:
                continue
            seen.add(l)
            if any(kw in l for kw in ["Ratigan", "Spiers", "De Vries", "Anhydrite", "Mudstone"]):
                b_handles.append(h); b_labels.append(l)
            else:
                s_handles.append(h); s_labels.append(l)
        all_handles, all_labels = [], []
        if s_handles:
            all_handles.append(plt.Line2D([], [], color="none"))
            all_labels.append(r"$\bf{Cavern\ geometries}$")
            all_handles.extend(s_handles); all_labels.extend(s_labels)
        if b_handles:
            all_handles.append(plt.Line2D([], [], color="none"))
            all_labels.append(r"$\bf{Dilatancy\ criteria}$")
            all_handles.extend(b_handles); all_labels.extend(b_labels)
        if all_handles:
            fig.legend(all_handles, all_labels, loc="upper center",
                       bbox_to_anchor=(0.5, 0.99),
                       ncol=min(4, len(all_labels)), frameon=True, fontsize=18)

        fig.tight_layout(rect=[0, 0, 1, 0.82])

        safe_name = f"{cav}_{sc}_{ps}".replace(" ", "_").replace("/", "_").replace("None", "")
        outname = f"pq_paths_{safe_name}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


def plot_stress_per_cavern(cases, group_fn=None):
    """One stress-state figure per cavern, scenarios/pressures overlaid."""
    probe_types = ["top", "quarter", "mid", "threequarter", "bottom"]
    groups = (group_fn or group_cases_by_cavern)(cases)

    for cav_label, cav_cases in groups.items():
        # Build stress paths for all cases of this cavern
        stress_by_series = {}
        series_meta = {}
        for m in cav_cases:
            folder = m["case_path"]
            try:
                wall_points = load_wall_points(folder)
                probes = auto_generate_probes_from_wall_points(wall_points)
                series_key = get_series_key(m)
                stress_by_series[series_key] = read_stress_paths(folder, probes)
                series_meta[series_key] = m
            except Exception as e:
                print(f"[WARN] Skipping stress state for {m.get('case_name')}: {e}")
        if not stress_by_series:
            print(f"[STRESS STATE] No valid data for {cav_label}")
            continue

        fig, axes = plt.subplots(2, 3, figsize=(20, 11))
        axes = axes.flatten()

        panel_letters = ["A", "B", "C", "D", "E", "F"]

        for i, ptype in enumerate(probe_types):
            ax = axes[i]
            plot_dilatancy_boundaries(ax, show_boundaries=SHOW_DILATANCY)

            for (cav, sc, ps), d in stress_by_series.items():
                p, q = d[ptype]
                color, linestyle = get_case_color_and_style(cav, sc, ps)
                label = get_case_label({"cavern_label": cav, "scenario_preset": sc, "pressure_scenario": ps})

                ax.plot(p, q, linewidth=2.0, color=color, linestyle=linestyle, label=label)
                ax.scatter(p[-1], q[-1], s=30, edgecolors="black", linewidths=0.6,
                           color=color, zorder=5)

            ax.set_title(f"p-q stress path: {ptype}", fontsize=22, fontweight='bold')
            ax.set_xlabel("Mean stress p (MPa)", fontsize=22)
            ax.set_ylabel("Differential stress q (MPa)", fontsize=22)
            ax.tick_params(axis='both', labelsize=20)
            ax.grid(True, alpha=0.3)
            add_panel_label(ax, panel_letters[i])

        # Pressure subplot — overlay all distinct pressure schedules
        axp = axes[5]
        _plotted_pressures = set()
        for c in cav_cases:
            ps = c.get("pressure_scenario")
            if ps in _plotted_pressures:
                continue
            _plotted_pressures.add(ps)
            tH, pMPa = read_pressure_schedule(c["case_path"])
            if tH is not None:
                axp.plot(tH / 24.0, pMPa, linewidth=2.0, color='black',
                         label=PRESSURE_LABELS.get(ps, ps))
        if not _plotted_pressures:
            axp.text(0.5, 0.5, "No pressure_schedule.json found.",
                     ha="center", va="center", transform=axp.transAxes)
            axp.axis("off")
        else:
            axp.set_title("Pressure schedule", fontsize=22, fontweight='bold')
            axp.set_xlabel("Time (days)", fontsize=22)
            axp.set_ylabel("Pressure (MPa)", fontsize=22)
            axp.tick_params(axis='both', labelsize=20)
            axp.grid(True, alpha=0.3)
            if len(_plotted_pressures) > 1:
                axp.legend(fontsize=18, frameon=True)
            add_panel_label(axp, panel_letters[5])

        # Legend above the figure, with section headers that separate cavern
        # geometries / scenarios from dilatancy criteria.
        b_handles, b_labels = [], []
        s_handles, s_labels = [], []
        handles, labels = axes[0].get_legend_handles_labels()
        seen = set()
        for h, l in zip(handles, labels):
            if l in seen:
                continue
            seen.add(l)
            if any(kw in l for kw in ["Ratigan", "Spiers", "De Vries", "Anhydrite", "Mudstone"]):
                b_handles.append(h); b_labels.append(l)
            else:
                s_handles.append(h); s_labels.append(l)
        all_handles, all_labels = [], []
        if s_handles:
            all_handles.append(plt.Line2D([], [], color="none"))
            all_labels.append(r"$\bf{Scenarios}$")
            all_handles.extend(s_handles); all_labels.extend(s_labels)
        if b_handles:
            all_handles.append(plt.Line2D([], [], color="none"))
            all_labels.append(r"$\bf{Dilatancy\ criteria}$")
            all_handles.extend(b_handles); all_labels.extend(b_labels)
        if all_handles:
            fig.legend(all_handles, all_labels, loc="upper center",
                       bbox_to_anchor=(0.5, 0.99),
                       ncol=min(4, len(all_labels)), frameon=True, fontsize=18)

        fig.tight_layout(rect=[0, 0, 1, 0.82])

        safe_cav = cav_label.replace(" ", "_")
        outname = f"pq_paths_{safe_cav}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


# =============================================================================
# 11. FIGURE 3 PLOTTING — FOS
# =============================================================================

def _compute_fos_for_case(c):
    """Load data and compute FOS for a single case. Returns dict with all needed arrays."""
    time_list_pq, p_elems, q_elems = load_p_q(c["case_path"])
    time_list_sig, sig_vals, sig_path = load_sig(c["case_path"])

    n_times = min(len(time_list_pq), len(time_list_sig))
    time_list = time_list_pq[:n_times]
    p_elems = p_elems[:n_times]
    q_elems = q_elems[:n_times]
    sig_vals = sig_vals[:n_times]

    FOS = compute_FOS_field(p_elems, q_elems, sig_vals, compression_positive=True, q_tol_MPa=1e-3)

    if WRITE_FOS_XDMF:
        write_fos_xdmf_from_sig(c["case_path"], sig_path, time_list, FOS)

    fos_min_global = np.nanmin(np.where(np.isfinite(FOS), FOS, np.nan), axis=1)

    geom_msh = path_geom_msh(c["case_path"])
    slices = extract_cavern_wall_cells_slices(geom_msh, cavern_phys_tag=CAVERN_PHYS_TAG,
                                              top_fraction=0.20, bottom_fraction=0.20)

    roof_stats  = fos_series_stats(FOS, slices["roof"])
    mid_stats   = fos_series_stats(FOS, slices["mid"])
    floor_stats = fos_series_stats(FOS, slices["floor"])

    t_days = time_list / DAY

    return {
        "t_days": t_days,
        "fos_min_global": fos_min_global,
        "roof_stats": roof_stats,
        "mid_stats": mid_stats,
        "floor_stats": floor_stats,
    }


def _plot_fos_on_axes(axes, fos_by_case, title_prefix=""):
    """
    Plot FOS on a 2x3 grid: 5 probe subplots + 1 min-FOS summary.

    fos_by_case: list of (label, color, linestyle, t_days, fos_by_probe)
    """
    panel_letters = ["A", "B", "C", "D", "E", "F"]
    probe_axes = []  # track probe subplots for shared y-limits

    for i, ptype in enumerate(PROBE_ORDER):
        ax = axes[i]
        probe_axes.append(ax)
        for label, col, ls, t_days, fos_probes in fos_by_case:
            if ptype not in fos_probes:
                continue
            fos = np.asarray(fos_probes[ptype])
            tx, fy = _downsample_xy(t_days, fos, max_points=FOS_MAX_POINTS)

            if FOS_SHOW_BAND:
                lo, med, hi = _rolling_quantile_band(fy, window=FOS_BAND_WINDOW)
                ax.fill_between(tx, lo, hi, color=col, alpha=0.12, linewidth=0)
                ax.plot(tx, med, linewidth=1.8, linestyle=ls, color=col, label=label)
            else:
                ax.plot(tx, fy, linewidth=1.6, linestyle=ls, color=col, alpha=0.9, label=label)

        ax.axhline(y=1.0, color='red', linestyle=':', linewidth=1.5, alpha=0.8)
        ax.set_title(f"FOS: {ptype}")
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Factor of Safety")
        ax.grid(True, alpha=0.3)
        add_panel_label(ax, panel_letters[i])

    # --- Summary subplot: global field minimum FOS ---
    ax_sum = axes[5]
    for label, col, ls, t_days, fos_probes in fos_by_case:
        if "global_min" in fos_probes:
            min_fos = fos_probes["global_min"]
        else:
            # Fallback: min across probes
            all_fos = np.array([fos_probes[p] for p in PROBE_ORDER if p in fos_probes])
            min_fos = np.min(all_fos, axis=0)
        tx, my = _downsample_xy(t_days, min_fos, max_points=FOS_MAX_POINTS)

        if FOS_SHOW_BAND:
            lo, med, hi = _rolling_quantile_band(my, window=FOS_BAND_WINDOW)
            ax_sum.fill_between(tx, lo, hi, color=col, alpha=0.12, linewidth=0)
            ax_sum.plot(tx, med, linewidth=2.0, linestyle=ls, color=col, label=label)
        else:
            ax_sum.plot(tx, my, linewidth=2.0, linestyle=ls, color=col, label=label)

    # Add labelled failure threshold line on summary (picked up by fig.legend)
    ax_sum.axhline(y=1.0, color='red', linestyle=':', linewidth=1.5, alpha=0.8,
                   label="FOS = 1 (failure threshold)")
    ax_sum.set_title("Min FOS (global field)")
    ax_sum.set_xlabel("Time (days)")
    ax_sum.set_ylabel("Factor of Safety")
    ax_sum.grid(True, alpha=0.3)
    add_panel_label(ax_sum, panel_letters[5])

    # --- Shared y-limits across the 5 probe subplots for easier visual comparison ---
    lo_vals, hi_vals = [], []
    for ax in probe_axes:
        for line in ax.get_lines():
            y = np.asarray(line.get_ydata(), dtype=float)
            y = y[np.isfinite(y)]
            if y.size:
                lo_vals.append(y.min()); hi_vals.append(y.max())
    if lo_vals and hi_vals:
        y_lo = min(lo_vals); y_hi = max(hi_vals)
        # Always include FOS = 1 in view and add a small pad
        y_lo = min(y_lo, 1.0); y_hi = max(y_hi, 1.0)
        pad = 0.05 * (y_hi - y_lo) if y_hi > y_lo else 0.1
        y_lo -= pad; y_hi += pad
        for ax in probe_axes:
            ax.set_ylim(y_lo, y_hi)

    # Add legend at top of figure (works for all callers)
    handles, labels = ax_sum.get_legend_handles_labels()
    if handles:
        fig = ax_sum.get_figure()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.99),
                   ncol=min(4, len(labels)), frameon=True, fontsize=18)


def plot_fos_combined(cases):
    """Single 2x3 figure with all cases overlaid (compare_shapes mode)."""
    fig, axes = plt.subplots(2, 3, figsize=(20, 11))
    axes = axes.flatten()

    fos_by_case = []
    for c in cases:
        cav = c.get("cavern_label")
        sc = c.get("scenario_preset")
        ps = c.get("pressure_scenario")
        col, ls = get_case_color_and_style(cav, sc, ps)
        label = get_case_label(c)

        try:
            t_days, fos_probes = _compute_fos_probes_for_case(c)
            min_val = min(np.min(v) for v in fos_probes.values())
            print(f"    [FOS] {label}: min FOS = {min_val:.3f}")
            fos_by_case.append((label, col, ls, t_days, fos_probes))
        except Exception as e:
            print(f"[SKIP] {cav}/{sc}: {e}")
            continue

    if not fos_by_case:
        plt.close(fig)
        return

    _plot_fos_on_axes(axes, fos_by_case)
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    outname = f"fos_combined_pressure={SELECT.get('pressure')}_scenario={SELECT.get('scenario')}.png"
    outpath = os.path.join(OUT_DIR, outname.replace(" ", ""))
    fig.savefig(outpath, dpi=DPI)
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close(fig)


def plot_fos_separate(cases):
    """One 2x3 figure per case."""
    for c in cases:
        cav = c.get("cavern_label")
        sc = c.get("scenario_preset")
        ps = c.get("pressure_scenario")
        col, _ = get_case_color_and_style(cav, sc, ps)
        label = get_case_label(c)

        try:
            t_days, fos_probes = _compute_fos_probes_for_case(c)
        except Exception as e:
            print(f"[SKIP] {cav}/{sc}: {e}")
            continue

        fig, axes = plt.subplots(2, 3, figsize=(20, 11))
        axes = axes.flatten()
        _plot_fos_on_axes(axes, [(label, col, "-", t_days, fos_probes)])
        fig.tight_layout(rect=[0, 0, 1, 0.93])

        safe_name = c.get("case_name", "unknown").replace(" ", "_")
        outname = f"fos_{safe_name}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


def plot_fos_per_cavern(cases, group_fn=None):
    """One 2x3 FOS figure per cavern, scenarios/pressures overlaid."""
    groups = (group_fn or group_cases_by_cavern)(cases)
    for cav_label, cav_cases in groups.items():
        fig, axes = plt.subplots(2, 3, figsize=(20, 11))
        axes = axes.flatten()

        fos_by_case = []
        for c in cav_cases:
            cav_full = c.get("cavern_label", cav_label)
            sc = c.get("scenario_preset")
            ps = c.get("pressure_scenario")
            col, ls = get_case_color_and_style(cav_full, sc, ps)
            label = get_case_label(c)

            try:
                t_days, fos_probes = _compute_fos_probes_for_case(c)
                min_val = min(np.min(v) for v in fos_probes.values())
                print(f"    [FOS] {label}: min FOS = {min_val:.3f}")
                fos_by_case.append((label, col, ls, t_days, fos_probes))
            except Exception as e:
                print(f"[SKIP] {cav_label}/{sc}: {e}")
                continue

        if not fos_by_case:
            plt.close(fig)
            continue

        _plot_fos_on_axes(axes, fos_by_case)
        fig.tight_layout(rect=[0, 0, 1, 0.93])

        safe_cav = cav_label.replace(" ", "_")
        outname = f"fos_{safe_cav}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


# =============================================================================
# 11b. FIGURE 5 — Global min FOS + pressure profiles side-by-side
# =============================================================================

# Order in which pressure scenarios appear (top to bottom) in the summary figure
PRESSURE_ORDER = ["industry", "transport", "power_generation", "csv"]


def plot_fos_summary(cases):
    """Global min FOS (left) next to stacked pressure profiles (right).

    Left panel:  one subplot with global field-minimum FOS for every case.
    Right panel: one subplot per distinct pressure scenario, stacked vertically.
    """
    import matplotlib.gridspec as gridspec

    # --- Collect FOS data ---
    fos_by_case = []
    for c in cases:
        cav = c.get("cavern_label")
        sc = c.get("scenario_preset")
        ps = c.get("pressure_scenario")
        col, ls = get_case_color_and_style(cav, sc, ps)
        label = get_case_label(c)

        try:
            t_days, fos_probes = _compute_fos_probes_for_case(c)
            fos_by_case.append((label, col, ls, t_days, fos_probes, c))
        except Exception as e:
            print(f"[SKIP] {cav}/{sc}: {e}")
            continue

    if not fos_by_case:
        print("[FOS SUMMARY] No valid cases.")
        return

    # --- Identify distinct pressure scenarios present ---
    seen_pressures = []
    for _, _, _, _, _, c in fos_by_case:
        ps = c.get("pressure_scenario")
        if ps not in seen_pressures:
            seen_pressures.append(ps)
    # Sort by PRESSURE_ORDER, then any remaining at the end
    ordered = [p for p in PRESSURE_ORDER if p in seen_pressures]
    ordered += [p for p in seen_pressures if p not in ordered]
    n_pressures = max(len(ordered), 1)

    # --- Create figure with gridspec ---
    fig = plt.figure(figsize=(20, max(8, 2.8 * n_pressures)))
    gs = gridspec.GridSpec(n_pressures, 2, width_ratios=[3, 2],
                           hspace=0.40, wspace=0.30)

    # Left: global min FOS spanning all rows
    ax_fos = fig.add_subplot(gs[:, 0])

    for label, col, ls, t_days, fos_probes, _ in fos_by_case:
        if "global_min" in fos_probes:
            min_fos = fos_probes["global_min"]
        else:
            all_fos = np.array([fos_probes[p] for p in PROBE_ORDER if p in fos_probes])
            min_fos = np.min(all_fos, axis=0)

        tx, my = _downsample_xy(t_days, min_fos, max_points=FOS_MAX_POINTS)

        if FOS_SHOW_BAND:
            lo, med, hi = _rolling_quantile_band(my, window=FOS_BAND_WINDOW)
            ax_fos.fill_between(tx, lo, hi, color=col, alpha=0.12, linewidth=0)
            ax_fos.plot(tx, med, linewidth=2.0, linestyle=ls, color=col, label=label)
        else:
            ax_fos.plot(tx, my, linewidth=2.0, linestyle=ls, color=col, label=label)

    ax_fos.axhline(y=1.0, color='red', linestyle=':', linewidth=1.5, alpha=0.8)
    ax_fos.set_title("Global min FOS (all cells)", fontsize=22, fontweight='bold')
    ax_fos.set_xlabel("Time (days)")
    ax_fos.set_ylabel("Factor of Safety")
    ax_fos.grid(True, alpha=0.3)
    ax_fos.legend(fontsize=18, frameon=True, loc="upper center",
                   bbox_to_anchor=(0.5, 1.15), ncol=4)

    # Right: one pressure subplot per scenario
    _pressure_plotted = set()
    for row, ps_key in enumerate(ordered):
        ax_p = fig.add_subplot(gs[row, 1])

        # Pick the first case that has this pressure scenario
        for _, _, _, _, _, c in fos_by_case:
            if c.get("pressure_scenario") == ps_key and ps_key not in _pressure_plotted:
                tH, pMPa = read_pressure_schedule(c["case_path"])
                if tH is not None:
                    pcol = PRESSURE_COLORS.get(ps_key, "#333333")
                    ax_p.plot(tH / 24.0, pMPa, linewidth=1.5, color=pcol)
                _pressure_plotted.add(ps_key)
                break

        ax_p.set_title(PRESSURE_LABELS.get(ps_key, ps_key), fontsize=20)
        ax_p.set_ylabel("Pressure (MPa)")
        ax_p.grid(True, alpha=0.3)
        if row == len(ordered) - 1:
            ax_p.set_xlabel("Time (days)")

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    outname = f"fos_summary_pressure={SELECT.get('pressure')}_scenario={SELECT.get('scenario')}.png"
    outpath = os.path.join(OUT_DIR, outname.replace(" ", ""))
    fig.savefig(outpath, dpi=DPI)
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close(fig)


# =============================================================================
# 12. FIGURE 4 PLOTTING — Fracture propagation (always per-case)
# =============================================================================

def plot_cavern_with_probes(ax, wall_points, sample_points, cavern_label=None):
    """Cavern wall plot with five colored region markers (top / quarter / mid /
    threequarter / bottom). The markers indicate which colour on the right-hand
    panel maps to which part of the cavern wall.
    """
    if cavern_label is not None:
        y_wall = z_mesh_to_depth_m(wall_points[:, 2], cavern_label)
        y_label = 'Depth (m)'
    else:
        y_wall = wall_points[:, 2]
        y_label = 'Height z (m)'

    ax.plot(wall_points[:, 0], y_wall, 'k-', linewidth=2, label='Cavern wall')

    z_w = wall_points[:, 2]
    z_min, z_max = float(z_w.min()), float(z_w.max())
    span = max(z_max - z_min, 1e-12)
    region_fracs = {"bottom": 0.0, "threequarter": 0.25, "mid": 0.5,
                    "quarter": 0.75, "top": 1.0}
    for region_name, frac in region_fracs.items():
        idx = int(np.argmin(np.abs(z_w - (z_min + frac * span))))
        color = PROBE_COLORS.get(region_name, 'gray')
        y_pt = z_mesh_to_depth_m(wall_points[idx, 2], cavern_label) \
            if cavern_label is not None else wall_points[idx, 2]
        ax.scatter(wall_points[idx, 0], y_pt, c=color, s=110, marker='o',
                   edgecolors='black', linewidths=0.7, zorder=5)

    ax.set_xlabel('Radial distance (m)')
    ax.set_ylabel(y_label)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    if cavern_label is not None:
        ax.invert_yaxis()


def plot_propagation_depth_global(ax, time_days, region_depths):
    """Plot dilation-zone penetration per region, computed via volumetric scan."""
    region_order = ["top", "quarter", "mid", "threequarter", "bottom"]
    for region_name in region_order:
        if region_name not in region_depths:
            continue
        color = PROBE_COLORS.get(region_name, 'gray')
        ls = PROBE_LINESTYLES.get(region_name, '-')
        lw = PROBE_LINEWIDTHS.get(region_name, 2.0)
        ax.plot(time_days, region_depths[region_name],
                linewidth=lw, color=color, linestyle=ls, label=region_name)

    ax.set_xlabel('Time (days)')
    ax.set_title('Connected dilating zone size')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0.0)
    ax.legend(loc='upper left', fontsize=18)


def plot_propagation_depth(ax, time_days, fos_data, radial_distances, threshold=1.0):
    for probe_name in fos_data.keys():
        color = PROBE_COLORS.get(probe_name, 'gray')

        max_depth = np.zeros(len(time_days))
        for t_idx in range(len(time_days)):
            depth = 0.0
            for dist in radial_distances:
                fos = fos_data[probe_name][dist][t_idx]
                if fos < threshold:
                    depth = dist
            max_depth[t_idx] = depth

        ax.plot(time_days, max_depth, linewidth=2, color=color, label=probe_name)

    ax.set_xlabel('Time (days)')
    ax.set_title('Connected dilating zone size')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0.0)
    ax.legend(loc='upper left', fontsize=18)


def plot_fracture_propagation(case_meta):
    """Create fracture propagation figure for a single case (always per-case)."""
    folder = case_meta["case_path"]
    cav = case_meta.get("cavern_label")
    sc = case_meta.get("scenario_preset")
    ps = case_meta.get("pressure_scenario")
    ps_label = PRESSURE_LABELS.get(ps, ps or "unknown")
    case_label = f"{cav} | {sc} | {ps_label}" if sc else f"{cav} | {ps_label}"

    print(f"\n[PROCESSING fracture propagation] {cav} | {sc} | {ps}")

    try:
        wall_points = load_wall_points(folder)
        time_days, region_depths = compute_global_dilation_penetration(folder, wall_points)

        fig, (ax_cavern, ax_depth) = plt.subplots(1, 2, figsize=(20, 10))

        plot_cavern_with_probes(ax_cavern, wall_points, None, cavern_label=cav)
        plot_propagation_depth_global(ax_depth, time_days, region_depths)
        ax_depth.set_ylabel("Dilatancy zone depth (m)")

        add_panel_label(ax_cavern, "A", fontsize=14)
        add_panel_label(ax_depth, "B", fontsize=14)
        fig.tight_layout()

        safe_name = f"{cav}_{sc}_{ps}".replace(" ", "_").replace("/", "_").replace("None", "none")
        outname = f"fracture_propagation_FOS_{safe_name}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print(f"[SAVED] {outpath}")

        if SHOW:
            plt.show()
        plt.close(fig)

    except Exception as e:
        print(f"[ERROR] Failed to process fracture propagation: {e}")
        import traceback
        traceback.print_exc()


# Groups for the combined dilatancy zone figures
FRACTURE_GROUPS = [
    {
        "label": "group1",
        "shapes": ["Regular", "Direct-circulation", "Reversed-circulation"],
    },
    {
        "label": "group2",
        "shapes": ["String-failure", "Fast-leached", "Tilt"],
    },
]


def _shade_dilatancy_bands(ax, wall_points, cavern_label, region_depths):
    """Shade the salt side of the cavern wall with the five region colours.

    Each band's vertical (z-fraction) extent matches the assignment used in
    ``compute_global_dilation_penetration``; its horizontal extent is the
    maximum contiguous dilatancy depth seen for that region over the whole
    simulation. This makes the cavern plot show, at a glance, what the
    coloured curves on the right panel represent.
    """
    region_order = ["bottom", "threequarter", "mid", "quarter", "top"]
    fracs = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    edges = np.concatenate(([0.0], 0.5 * (fracs[:-1] + fracs[1:]), [1.0]))

    z_w = wall_points[:, 2]
    x_w = wall_points[:, 0]
    z_min, z_max = float(z_w.min()), float(z_w.max())
    span = max(z_max - z_min, 1e-12)

    order = np.argsort(z_w)
    z_sorted = z_w[order]
    x_sorted = x_w[order]
    z_frac_sorted = (z_sorted - z_min) / span

    for i, region in enumerate(region_order):
        if region not in region_depths:
            continue
        depth_series = region_depths[region]
        max_depth = float(np.max(depth_series)) if len(depth_series) else 0.0
        if max_depth <= 0:
            continue
        f_lo, f_hi = edges[i], edges[i + 1]
        m = (z_frac_sorted >= f_lo) & (z_frac_sorted <= f_hi)
        if not m.any():
            continue
        z_band = z_sorted[m]
        x_band = x_sorted[m]
        y_band = z_mesh_to_depth_m(z_band, cavern_label) \
            if cavern_label is not None else z_band
        ax.fill_betweenx(y_band, x_band, x_band + max_depth,
                         color=PROBE_COLORS.get(region, 'gray'),
                         alpha=0.30, linewidth=0, zorder=1)


def plot_fracture_propagation_grouped(frac_cases):
    """Create two grouped dilatancy zone figures, 3 cavern shapes each.

    Each figure has 3 columns (one per shape). Each column has a cavern
    profile subplot (left-ish) and a dilating-zone-size subplot (right-ish),
    tightly packed.
    """
    # Build lookup: cavern_label -> case_meta
    # Strip volume suffix (e.g. " (1,200,000 m³)") so bare names like "Regular" match
    case_by_label = {}
    for c in frac_cases:
        cav = c.get("cavern_label", "")
        if cav not in case_by_label:
            case_by_label[cav] = c
        bare = cav.split(" (")[0] if " (" in cav else cav
        if bare not in case_by_label:
            case_by_label[bare] = c

    for group in FRACTURE_GROUPS:
        shapes = group["shapes"]
        group_label = group["label"]

        # Collect data for available shapes in this group
        group_data = []
        for shape_name in shapes:
            if shape_name not in case_by_label:
                print(f"[SKIP] {shape_name} not found for fracture group")
                continue
            c = case_by_label[shape_name]
            folder = c["case_path"]
            try:
                wall_points = load_wall_points(folder)
                time_days, region_depths = compute_global_dilation_penetration(
                    folder, wall_points
                )
                group_data.append((shape_name, wall_points, None,
                                   time_days, region_depths))
            except Exception as e:
                print(f"[ERROR] {shape_name}: {e}")
                import traceback
                traceback.print_exc()

        if not group_data:
            continue

        n_shapes = len(group_data)
        # Wider cavern subplot so the shaded dilatancy-band overlay fits on the salt side
        fig, axes = plt.subplots(1, n_shapes * 2, figsize=(7.5 * n_shapes, 9),
                                 gridspec_kw={"width_ratios": [1.4, 1.0] * n_shapes,
                                              "wspace": 0.20})
        if n_shapes * 2 == 2:
            axes = [axes[0], axes[1]]

        panel_letters = ["A", "B", "C", "D", "E", "F"]
        title_fs = 14
        axlabel_fs = 14
        tick_fs = 16
        panel_fs = 12
        dep_axes = []
        for idx, (shape_name, wall_pts, _samp_pts, t_days, region_depths) in enumerate(group_data):
            ax_cav = axes[idx * 2]
            ax_dep = axes[idx * 2 + 1]
            dep_axes.append(ax_dep)

            # Look up the cavern_label (with volume tag) so depth conversion works
            cav_full = None
            for key, c in case_by_label.items():
                if key == shape_name or key.startswith(shape_name + " ("):
                    cav_full = c.get("cavern_label", shape_name)
                    if "(" in cav_full:
                        break

            plot_cavern_with_probes(ax_cav, wall_pts, None, cavern_label=cav_full)
            _shade_dilatancy_bands(ax_cav, wall_pts, cav_full, region_depths)
            ax_cav.set_title(shape_name, fontsize=title_fs, fontweight='bold', pad=6)
            ax_cav.set_xlabel("Radial distance (m)", fontsize=axlabel_fs)
            ax_cav.set_ylabel("Depth (m)" if idx == 0 else "", fontsize=axlabel_fs)
            ax_cav.tick_params(axis='both', labelsize=tick_fs)
            leg = ax_cav.get_legend()
            if leg is not None:
                leg.remove()
            add_panel_label(ax_cav, panel_letters[idx * 2], fontsize=panel_fs)

            plot_propagation_depth_global(ax_dep, t_days, region_depths)
            ax_dep.set_title("", fontsize=1)  # clear sub-title
            ax_dep.set_xlabel("Time (days)", fontsize=axlabel_fs)
            ax_dep.set_ylabel("")
            ax_dep.tick_params(axis='both', labelsize=tick_fs)
            leg_dep = ax_dep.get_legend()
            if leg_dep is not None:
                leg_dep.remove()
            add_panel_label(ax_dep, panel_letters[idx * 2 + 1], fontsize=panel_fs)

        if dep_axes:
            for ax in dep_axes:
                ax.set_ylim(0.0, 35.0)

        # Single legend above the row of subplots
        region_order = ["top", "quarter", "mid", "threequarter", "bottom"]
        legend_handles = [
            plt.Line2D([0], [0],
                       color=PROBE_COLORS.get(r, 'gray'),
                       linestyle=PROBE_LINESTYLES.get(r, '-'),
                       linewidth=PROBE_LINEWIDTHS.get(r, 2.0),
                       label=r)
            for r in region_order
        ]
        fig.legend(handles=legend_handles, loc='upper center', ncol=5,
                   fontsize=13, frameon=True, bbox_to_anchor=(0.5, 0.99))

        fig.tight_layout(rect=[0, 0, 1, 0.94])

        outname = f"dilatancy_zone_{group_label}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print(f"[SAVED] {outpath}")

        if SHOW:
            plt.show()
        plt.close(fig)


# =============================================================================
# 12b. FIGURE 6 — Mohr-Coulomb failure (interlayer cases)
# =============================================================================

def _get_interlayer_cell_mask(case_meta, n_cells_xdmf):
    """Return (mask, volumes_xdmf) for interlayer cells.

    mask: boolean array (XDMF cell order) marking interlayer cells.
    volumes_xdmf: float array (XDMF cell order) with per-cell tet volume (m³)
                  for interlayer cells, 0 elsewhere.
    Returns (None, None) if no interlayer cells are found.
    """
    from scipy.spatial import cKDTree

    msh_path = path_geom_msh(case_meta["case_path"])
    mesh, cell_tags, _ = dfx.io.gmshio.read_from_msh(msh_path, MPI.COMM_WORLD, 0)

    # Known interlayer tags from the spike/heterogeneous grids
    INTERLAYER_TAGS = {32, 34}  # Interlayer_1=32, Interlayer_2=34
    il_dfx_cells = np.where(np.isin(cell_tags.values, list(INTERLAYER_TAGS)))[0]
    if len(il_dfx_cells) == 0:
        return None, None

    # Compute dolfinx cell centroids and tet volumes
    tdim = mesh.topology.dim
    c2v = mesh.topology.connectivity(tdim, 0)
    X = mesh.geometry.x
    n_cells = mesh.topology.index_map(tdim).size_local
    dfx_centroids = np.zeros((n_cells, 3))
    dfx_volumes = np.zeros(n_cells)
    for i in range(n_cells):
        verts = X[c2v.links(i)]
        dfx_centroids[i] = verts.mean(axis=0)
        # Tet volume via scalar triple product: |det([v1-v0, v2-v0, v3-v0])| / 6
        if verts.shape[0] == 4:
            e1 = verts[1] - verts[0]
            e2 = verts[2] - verts[0]
            e3 = verts[3] - verts[0]
            dfx_volumes[i] = abs(np.dot(e1, np.cross(e2, e3))) / 6.0

    # Get XDMF centroids from p_elems
    p_path = path_p_xdmf(case_meta["case_path"])
    centroids_xdmf, _, _ = post.read_cell_scalar(p_path)

    # Map dolfinx → XDMF ordering
    tree = cKDTree(centroids_xdmf)
    _, xdmf_idx = tree.query(dfx_centroids[il_dfx_cells])

    il_mask = np.zeros(n_cells_xdmf, dtype=bool)
    il_mask[xdmf_idx] = True
    il_volumes = np.zeros(n_cells_xdmf)
    il_volumes[xdmf_idx] = dfx_volumes[il_dfx_cells]
    print(f"  [MC] {case_meta.get('case_name')}: {il_mask.sum()} interlayer cells, "
          f"total volume {il_volumes.sum():.2f} m³")
    return il_mask, il_volumes


def _compute_mc_failure_timeseries(case_meta, c_MPa=4.0, phi_deg=35.0):
    """Compute MC yield function, failed cell count, and FOS for interlayer cells.

    Uses the raw stress tensor (sig.xdmf) instead of the smoothed p/q fields,
    because the smoother averages thin interlayer cells with surrounding salt
    and washes out the peak stress.

    Returns
    -------
    t_days : ndarray, shape (nt,)
    f_max_il : ndarray, shape (nt,)
        Maximum f value across interlayer cells (NaN if no interlayer).
    frac_failed_vol : ndarray, shape (nt,)
        Fraction of total interlayer (anhydrite) volume with f > 0 at each timestep.
    fos_mc_min_il : ndarray, shape (nt,)
        Minimum MC FOS across interlayer cells (NaN if no interlayer).
    has_interlayer : bool
        Whether this case has interlayer cells.
    total_il_vol : float
        Total interlayer volume (m³); 0 if no interlayer.
    """
    folder = case_meta["case_path"]

    # Load raw stress tensor — unsmoothed, per-element
    t_s, sig33, _ = load_sig(folder)
    t_days = t_s / DAY
    nt, nc = sig33.shape[0], sig33.shape[1]

    il_mask, il_volumes = _get_interlayer_cell_mask(case_meta, nc)
    has_interlayer = il_mask is not None and il_mask.any()

    if not has_interlayer:
        return (t_days, np.full(nt, np.nan), np.full(nt, np.nan),
                np.full(nt, np.nan), False, 0.0)

    # Raw stress for interlayer cells only
    sig_il = sig33[:, il_mask, :, :]  # (nt, n_il, 3, 3)
    vol_il = il_volumes[il_mask]  # (n_il,) per-cell volumes in interlayer order
    total_il_vol = float(vol_il.sum())

    # Compute p and q from raw tensor (SafeInCave: sigma negative in compression)
    I1 = sig_il[:, :, 0, 0] + sig_il[:, :, 1, 1] + sig_il[:, :, 2, 2]
    I2 = (sig_il[:, :, 0, 0] * sig_il[:, :, 1, 1]
        + sig_il[:, :, 1, 1] * sig_il[:, :, 2, 2]
        + sig_il[:, :, 0, 0] * sig_il[:, :, 2, 2]
        - sig_il[:, :, 0, 1]**2
        - sig_il[:, :, 0, 2]**2
        - sig_il[:, :, 1, 2]**2)
    J2 = np.maximum((1.0 / 3.0) * I1**2 - I2, 0.0)

    p_MPa = -I1 / (3.0 * MPA)  # compression-positive
    q_MPa = np.sqrt(3.0 * J2) / MPA

    phi = np.radians(phi_deg)
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    alpha_F = 6.0 * sin_phi / (3.0 - sin_phi)
    k_F = 6.0 * c_MPa * cos_phi / (3.0 - sin_phi)

    # Yield function: f = q - alpha_F * p - k_F  (f > 0 means failure)
    f_il = q_MPa - alpha_F * p_MPa - k_F
    f_max_il = np.nanmax(f_il, axis=1)

    # Volume fraction of anhydrite in MC failure (f > 0)
    failed = (f_il > 0).astype(float)  # (nt, n_il)
    vol_failed = failed @ vol_il        # (nt,) total failed volume per timestep
    frac_failed_vol = vol_failed / total_il_vol if total_il_vol > 0 else np.zeros(nt)

    # MC FOS = q_boundary / q
    q_tol = 1e-3  # MPa
    q_boundary = np.maximum(alpha_F * p_MPa + k_F, 0.0)
    fos_mc_il = np.full_like(q_MPa, np.inf)
    mask_q = q_MPa >= q_tol
    fos_mc_il[mask_q] = q_boundary[mask_q] / q_MPa[mask_q]
    fos_mc_min_il = np.nanmin(
        np.where(np.isfinite(fos_mc_il), fos_mc_il, np.nan), axis=1)

    return t_days, f_max_il, frac_failed_vol, fos_mc_min_il, has_interlayer, total_il_vol


def plot_mc_failure_combined(cases, c_MPa=4.0, phi_deg=35.0):
    """Plot MC failure: yield function f, failed cell count, pressure, and FOS for interlayer cells."""
    fig, axes = plt.subplots(3, 1, figsize=(14, 13), sharex=True,
                             gridspec_kw={"height_ratios": [3, 3, 1]})
    ax_f, ax_n, ax_p = axes

    pressure_plotted = False
    total_vol_by_label = {}  # for footer annotation: label -> total interlayer volume
    for c in cases:
        label = get_case_label(c)
        color, _ = get_case_color_and_style(
            c.get("cavern_label"), c.get("scenario_preset"),
            c.get("pressure_scenario"))
        ls = "-"  # always solid for MC failure plot

        try:
            t_days, f_max_il, frac_failed_vol, fos_mc_min_il, has_il, total_il_vol = \
                _compute_mc_failure_timeseries(c, c_MPa, phi_deg)
        except Exception as e:
            print(f"[WARN] MC failure skipped for {c.get('case_name')}: {e}")
            continue

        if not has_il:
            print(f"  [MC] {c.get('case_name')}: no interlayer — skipped")
            continue

        ax_f.plot(t_days, f_max_il, color=color, linestyle=ls, label=label)
        ax_n.plot(t_days, 100.0 * frac_failed_vol, color=color, linestyle=ls, label=label)
        total_vol_by_label[label] = total_il_vol

        # Plot pressure schedule (once — same for all cases with same pressure)
        if not pressure_plotted:
            tH, pMPa = read_pressure_schedule(c["case_path"])
            if tH is not None:
                ax_p.plot(tH / 24.0, pMPa, color='#555555', linewidth=1.5)
                pressure_plotted = True

    # Yield function panel (interlayer only)
    # Labelled axhline so the failure-threshold explanation lives in the legend
    # (no inline text over the curves).
    ax_f.axhline(0, color='#d62728', linewidth=1.2, linestyle='--', alpha=0.8,
                 label="$f_{MC} = 0$ (failure threshold; $f>0$: failure)")
    ax_f.set_ylabel("Max $f_{MC}$ (MPa)")
    ax_f.set_title(f"Mohr-Coulomb yield function — interlayer cells (c={c_MPa} MPa, $\\varphi$={phi_deg}°)",
                   fontsize=20, fontweight='bold')
    ax_f.legend(loc='best', fontsize=16)

    # Volume-fraction panel (interlayer only)
    ax_n.set_ylabel("Anhydrite in MC failure (vol. %)")
    ax_n.set_ylim(0.0, 10.0)
    ax_n.set_title("Fraction of interlayer (anhydrite) volume in Mohr-Coulomb failure",
                   fontsize=20)
    ax_n.legend(loc='best', fontsize=16)

    # Show total interlayer volume per case under the fraction panel so the
    # reader can convert percent back to absolute m³ if needed.
    if total_vol_by_label:
        summary = "  ".join(f"{lbl}: {v:,.0f} m³" for lbl, v in total_vol_by_label.items())
        ax_n.text(0.99, 0.98, f"Total interlayer volume — {summary}",
                  transform=ax_n.transAxes, fontsize=12, style='italic',
                  va='top', ha='right', color='#444444')

    # Pressure schedule panel
    ax_p.set_ylabel("P (MPa)")
    ax_p.set_xlabel("Time (days)")
    if not pressure_plotted:
        ax_p.text(0.5, 0.5, "No pressure schedule", ha="center", va="center",
                  transform=ax_p.transAxes, fontsize=16)

    fig.tight_layout()
    outpath = os.path.join(OUT_DIR, "mc_failure.png")
    fig.savefig(outpath, dpi=DPI)
    print(f"[SAVED] {outpath}")
    if SHOW:
        plt.show()
    plt.close(fig)


# =============================================================================
# 13. MAIN
# =============================================================================

def main():
    if MPI.COMM_WORLD.rank != 0:
        return

    mode = PLOT_MODE.lower().replace(" ", "_")

    print("=" * 70)
    print("COMBINED PLOT SCRIPT — plot_results.py")
    print("=" * 70)
    print(f"  ROOT:        {ROOT}")
    print(f"  PLOT_MODE:   {PLOT_MODE}")
    print(f"  Caverns:     {SELECT.get('caverns', 'all')}")
    print(f"  Pressure:    {SELECT.get('pressure', 'all')}")
    print(f"  Scenario:    {SELECT.get('scenario', 'all')}")
    print(f"  Contains:    {SELECT.get('case_contains', 'any')}")
    print(f"  Figures:     {[k for k,v in FIGURES.items() if v]}")
    print("=" * 70)

    os.makedirs(OUT_DIR, exist_ok=True)

    all_cases = detect_layout_and_collect_cases(ROOT)
    cases = filter_cases(all_cases, SELECT)
    if not cases:
        print("[DEBUG] Found (first 20):")
        for m in all_cases[:20]:
            print(" -", m.get("cavern_label"), m.get("scenario_preset"),
                  m.get("pressure_scenario"), m.get("case_name"))
        raise RuntimeError(f"No cases matched SELECT={SELECT}")

    # In compare_shapes mode, deduplicate so there is one case per series
    if mode == "compare_shapes":
        cases = pick_one_case_per_series(cases)

    print(f"[INFO] {len(cases)} case(s) after filtering:")
    for m in cases:
        print(f"  - {m.get('cavern_label')} | {m.get('scenario_preset')} | {m.get('case_name')}")

    # --- Figure 1: Convergence ---
    if FIGURES.get("convergence"):
        conv_cases = [c for c in cases if case_has_convergence_files(c["case_path"])]
        if conv_cases:
            print(f"\n[CONVERGENCE] {len(conv_cases)} case(s) with required files")
            if mode == "compare_shapes":
                plot_convergence_combined(conv_cases)
                plot_initial_vs_final_shapes(conv_cases)
            elif mode == "compare_sizes":
                plot_convergence_per_cavern(conv_cases, group_fn=group_cases_by_base_shape)
            elif mode in ("compare_scenarios", "compare_pressures"):
                plot_convergence_per_cavern(conv_cases)
            else:
                print(f"[WARN] Unknown PLOT_MODE '{PLOT_MODE}', using compare_shapes")
                plot_convergence_combined(conv_cases)
                plot_initial_vs_final_shapes(conv_cases)
        else:
            print("[CONVERGENCE] No cases with required files (u.xdmf + geom.msh)")

    # --- Figure 2: Stress state ---
    if FIGURES.get("stress_state"):
        stress_cases = [c for c in cases if case_has_stress_files(c["case_path"])]
        if stress_cases:
            print(f"\n[STRESS STATE] {len(stress_cases)} case(s) with required files")
            if mode == "compare_shapes":
                stress_by_series = {}
                for m in stress_cases:
                    folder = m["case_path"]
                    try:
                        wall_points = load_wall_points(folder)
                        probes = auto_generate_probes_from_wall_points(wall_points)
                        series_key = get_series_key(m)
                        stress_by_series[series_key] = read_stress_paths(folder, probes)
                    except Exception as e:
                        print(f"[WARN] Skipping stress state for {m.get('case_name')}: {e}")
                if stress_by_series:
                    plot_stress_combined(stress_cases, stress_by_series)
                else:
                    print("[STRESS STATE] No valid stress data could be loaded.")
            elif mode == "compare_sizes":
                plot_stress_per_cavern(stress_cases, group_fn=group_cases_by_base_shape)
            elif mode in ("compare_scenarios", "compare_pressures"):
                plot_stress_per_cavern(stress_cases)
            else:
                print(f"[WARN] Unknown PLOT_MODE '{PLOT_MODE}', skipping stress state")
        else:
            print("[STRESS STATE] No cases with required files (u.xdmf + geom.msh + p/q_elems)")

    # --- Figure 3: FOS ---
    if FIGURES.get("fos"):
        fos_cases = [c for c in cases if case_has_fos_files(c["case_path"])]
        if fos_cases:
            print(f"\n[FOS] {len(fos_cases)} case(s) with required files")
            if mode == "compare_shapes":
                plot_fos_combined(fos_cases)
            elif mode == "compare_sizes":
                plot_fos_per_cavern(fos_cases, group_fn=group_cases_by_base_shape)
            elif mode in ("compare_scenarios", "compare_pressures"):
                plot_fos_per_cavern(fos_cases)
            else:
                print(f"[WARN] Unknown PLOT_MODE '{PLOT_MODE}', using compare_shapes")
                plot_fos_combined(fos_cases)
        else:
            print("[FOS] No cases with required files (geom.msh + p/q_elems + sig.xdmf)")

    # --- Figure 5: FOS summary + pressure profiles ---
    if FIGURES.get("fos_summary"):
        fos_sum_cases = [c for c in cases if case_has_fos_files(c["case_path"])]
        if fos_sum_cases:
            print(f"\n[FOS SUMMARY] {len(fos_sum_cases)} case(s) with required files")
            plot_fos_summary(fos_sum_cases)
        else:
            print("[FOS SUMMARY] No cases with required files (geom.msh + p/q_elems + sig.xdmf)")

    # --- Figure 4: Fracture propagation ---
    if FIGURES.get("fracture_propagation"):
        frac_cases = [c for c in cases if case_has_fracture_files(c["case_path"])]
        if frac_cases:
            print(f"\n[FRACTURE PROPAGATION] {len(frac_cases)} case(s) with required files")
            if mode == "compare_shapes":
                plot_fracture_propagation_grouped(frac_cases)
            elif mode == "compare_sizes":
                # One figure per base shape with sizes overlaid
                for base, grp in group_cases_by_base_shape(frac_cases).items():
                    for case_meta in grp:
                        plot_fracture_propagation(case_meta)
            else:
                for case_meta in frac_cases:
                    plot_fracture_propagation(case_meta)
        else:
            print("[FRACTURE PROPAGATION] No cases with required files (u.xdmf + geom.msh + p/q_elems + sig.xdmf)")

    # --- Figure 6: MC failure ---
    if FIGURES.get("mc_failure"):
        mc_cases = [c for c in cases
                    if os.path.isfile(path_sig_xdmf(c["case_path"]))
                    and os.path.isfile(path_geom_msh(c["case_path"]))]
        if mc_cases:
            print(f"\n[MC FAILURE] {len(mc_cases)} case(s) with required files")
            plot_mc_failure_combined(mc_cases)
        else:
            print("[MC FAILURE] No cases with required files (sig.xdmf + geom.msh)")

    print("\n[DONE]")


if __name__ == "__main__":
    main()