import os
import json
import numpy as np
import matplotlib.pyplot as plt
import dolfinx as dfx
from mpi4py import MPI

import safeincave.PostProcessingTools as post
from case_index import detect_layout_and_collect_cases, filter_cases

DAY = 24.0 * 3600.0
MPA = 1e6

# =============================================================================
# USER CONFIGURATION - Edit this section to customize the plot
# =============================================================================

# --- Output folder containing simulation results ---
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/nobian/Simulation/output"

# --- Case selection filters ---
# Set any filter to None to include all values for that parameter
SELECT = {
    "caverns": None,                           # e.g. ["Regular", "Tilted"] or None for all
    "pressure": None,                          # "sinus", "linear", "irregular", "csv_profile", or None
    "scenario": None,                          # e.g. ["full", "desai_only"] or None
    "case_contains": None,                     # substring filter on case name, or None
}

# --- Plot mode ---
# "combined" = all cases on one figure
# "separate" = one figure per case
PLOT_MODE = "combined"

# --- Color coding by cavern shape ---
CAVERN_COLORS = {
    "Asymmetric":          "#1f77b4",   # blue
    "Direct-circulation":  "#ff7f0e",   # orange
    "IrregularFine":       "#d62728",   # red
    "Regular":             "#9467bd",   # purple
    "Reversed-circulation":"#8c564b",   # brown
    "Tilt":                "#e377c2",   # pink
    "Fast-leached":        "#17becf",   # cyan
    "Tube-failure":        "#bcbd22",   # olive
}

SCENARIO_LINESTYLES = {
    "disloc_old_only":        "-",
    "disloc_new_only":        "--",
    "desai_only":             "-.",
    "full_minus_desai":       "-",
    "full":                   "-",
    "full_minus_ps": "-",
    # Munson-Dawson scenarios
    "md_only":                "-",
    "md_steady_only":         "--",
    "full_md":                "-",
    # Interlayer scenarios
    "interlayer":             "-",
    "nointerlayer":           "--",
    None:                     "-",
}

SCENARIO_COLORS = {
    "disloc_old_only":        "#1f77b4",   # blue
    "disloc_new_only":        "#ff7f0e",   # orange
    "desai_only":             "#2ca02c",   # green
    "full_minus_desai":       "#2ca02c",   # green
    "full":                   "#d62728",   # red
    "full_minus_ps": "#d62728",   # red
    # Munson-Dawson scenarios
    "md_only":                "#1f77b4",   # blue
    "md_steady_only":         "#bcbd22",   # olive
    "full_md":                "#e377c2",   # pink
    # Interlayer scenarios
    "interlayer":             "#7f7f7f",   # gray
    "nointerlayer":           "#8c564b",   # brown
    None:                     "#333333",   # dark gray
}
##

# --- Ordering for legend ---
CAVERN_ORDER = ["Asymmetric", "Direct-circulation", "Regular", "Reversed-circulation", "Tilt", "Fast-leached", "Tube-failure", "IrregularFine"]
SCENARIO_ORDER = ["disloc_old_only", "disloc_new_only", "desai_only", "full_minus_desai", "full",
                  "full_minus_ps", "md_only", "md_steady_only", "full_md",
                  "interlayer", "nointerlayer", None]

# --- Other settings ---
CAVERN_PHYS_TAG = 29          # Physical tag for cavern boundary in mesh
OUT_DIR = os.path.join(ROOT, "_figures")
SHOW = True                  # Show plot interactively after saving
DPI = 180

# =============================================================================
# END OF USER CONFIGURATION
# =============================================================================


def get_cavern_color(cavern_label):
    """Get color for a cavern label, with fallback."""
    return CAVERN_COLORS.get(cavern_label, "#333333")


def get_case_color(cavern_label, scenario):
    """Get color: use scenario color if scenario is set, otherwise cavern color."""
    if scenario is not None:
        return SCENARIO_COLORS.get(scenario, "#333333")
    return CAVERN_COLORS.get(cavern_label, "#333333")


def get_scenario_linestyle(scenario):
    """Get linestyle for a scenario, with fallback."""
    return SCENARIO_LINESTYLES.get(scenario, "-")

def path_u_xdmf(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.xdmf")

def path_u_h5(case_folder):
    return os.path.join(case_folder, "operation", "u", "u.h5")

def path_geom_msh(case_folder):
    return os.path.join(case_folder, "operation", "mesh", "geom.msh")

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
        t = np.asarray(data["t_values_s"], float) / 3600.0
        p = np.asarray(data["p_values_Pa"], float) / 1e6
        return t, p
    if "t_values" in data and "p_values" in data:
        t = np.asarray(data["t_values"], float) / 3600.0
        p = np.asarray(data["p_values"], float) / 1e6
        return t, p
    return None, None

def _area_vectors(points_xyz: np.ndarray, tri: np.ndarray):
    p0 = points_xyz[tri[:, 0]]
    p1 = points_xyz[tri[:, 1]]
    p2 = points_xyz[tri[:, 2]]
    return 0.5 * np.cross(p1 - p0, p2 - p0)

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

def print_config_summary():
    """Print configuration summary at startup."""
    print("=" * 60)
    print("CONVERGENCE PLOT CONFIGURATION")
    print("=" * 60)
    print(f"  ROOT:        {ROOT}")
    print(f"  PLOT_MODE:   {PLOT_MODE}")
    print(f"  Caverns:     {SELECT.get('caverns', 'all')}")
    print(f"  Pressure:    {SELECT.get('pressure', 'all')}")
    print(f"  Scenario:    {SELECT.get('scenario', 'all')}")
    print(f"  Contains:    {SELECT.get('case_contains', 'any')}")
    print("=" * 60)


def plot_combined(cases):
    """Plot all cases on a single figure."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 7), sharex=True,
                                   gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12})
    fig.suptitle(f"3D cavern volume convergence | pressure={SELECT.get('pressure')}")

    for c in cases:
        cav = c.get("cavern_label")
        sc = c.get("scenario_preset")
        col = get_case_color(cav, sc)
        ls = get_scenario_linestyle(sc)
        label = f"{cav} | {sc}" if sc is not None else f"{cav}"
        try:
            t_days, conv = compute_convergence_3d_percent(c["case_path"], cavern_phys_tag=CAVERN_PHYS_TAG)
        except Exception as e:
            print(f"[SKIP] {cav}/{sc}: {e}")
            continue
        ax1.plot(t_days, conv, linewidth=2.0, color=col, linestyle=ls, alpha=0.95, label=label)

    ax1.set_ylabel("Convergence (ΔV/V0) (%)")
    ax1.grid(True, alpha=0.3)

    # dedup legend
    h, l = ax1.get_legend_handles_labels()
    uniq = {}
    for hh, ll in zip(h, l):
        if ll not in uniq:
            uniq[ll] = hh
    if uniq:
        ax1.legend(uniq.values(), uniq.keys(), loc="best", fontsize=9, frameon=True)

    # pressure schedule from first case
    tH, pMPa = read_pressure_schedule(cases[0]["case_path"])
    if tH is None:
        ax2.text(0.5, 0.5, "No pressure_schedule.json found.", ha="center", va="center", transform=ax2.transAxes)
    else:
        ax2.plot(tH / 24.0, pMPa, linewidth=1.7)
    ax2.set_ylabel("Pressure (MPa)")
    ax2.set_xlabel("Time (days)")
    ax2.grid(True, alpha=0.3)

    outname = f"convergence_pressure={SELECT.get('pressure')}_scenario={SELECT.get('scenario')}.png"
    outpath = os.path.join(OUT_DIR, outname.replace(" ", ""))
    fig.savefig(outpath, dpi=DPI)
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close(fig)


def plot_separate(cases):
    """Plot each case on its own figure."""
    for c in cases:
        cav = c.get("cavern_label")
        sc = c.get("scenario_preset")
        col = get_case_color(cav, sc)
        label = f"{cav} | {sc}" if sc is not None else f"{cav}"

        try:
            t_days, conv = compute_convergence_3d_percent(c["case_path"], cavern_phys_tag=CAVERN_PHYS_TAG)
        except Exception as e:
            print(f"[SKIP] {cav}/{sc}: {e}")
            continue

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True,
                                       gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12})
        fig.suptitle(f"Convergence: {label}")

        ax1.plot(t_days, conv, linewidth=2.0, color=col, alpha=0.95, label=label)
        ax1.set_ylabel("Convergence (ΔV/V0) (%)")
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc="best", fontsize=9, frameon=True)

        tH, pMPa = read_pressure_schedule(c["case_path"])
        if tH is None:
            ax2.text(0.5, 0.5, "No pressure_schedule.json found.", ha="center", va="center", transform=ax2.transAxes)
        else:
            ax2.plot(tH / 24.0, pMPa, linewidth=1.7)
        ax2.set_ylabel("Pressure (MPa)")
        ax2.set_xlabel("Time (days)")
        ax2.grid(True, alpha=0.3)

        safe_name = c.get("case_name", "unknown").replace(" ", "_")
        outname = f"convergence_{safe_name}.png"
        outpath = os.path.join(OUT_DIR, outname)
        fig.savefig(outpath, dpi=DPI)
        print("[SAVED]", outpath)

        if SHOW:
            plt.show()
        plt.close(fig)


def main():
    print_config_summary()
    os.makedirs(OUT_DIR, exist_ok=True)

    all_cases = detect_layout_and_collect_cases(ROOT)
    # keep only convergence-required
    kept = []
    for m in all_cases:
        cp = m["case_path"]
        ok = (
            os.path.isfile(path_u_xdmf(cp)) and
            os.path.isfile(path_u_h5(cp)) and
            os.path.isfile(path_geom_msh(cp)) and
            is_valid_hdf5(path_u_h5(cp))
        )
        if ok:
            kept.append(m)
    all_cases = kept

    cases = filter_cases(all_cases, SELECT)
    if not cases:
        print("[DEBUG] Found (first 20):")
        for m in all_cases[:20]:
            print(" -", m.get("cavern_label"), m.get("scenario_preset"), m.get("pressure_scenario"), m.get("case_name"))
        raise RuntimeError(f"No cases matched SELECT={SELECT}")

    # one case per (cavern, scenario)
    by_series = {}
    for c in cases:
        key = (c.get("cavern_label"), c.get("scenario_preset"))
        if key not in by_series:
            by_series[key] = c
    cases = list(by_series.values())

    print(f"[INFO] Processing {len(cases)} case(s)...")

    if PLOT_MODE.lower() == "combined":
        plot_combined(cases)
    elif PLOT_MODE.lower() == "separate":
        plot_separate(cases)
    else:
        print(f"[WARNING] Unknown PLOT_MODE '{PLOT_MODE}', using 'combined'")
        plot_combined(cases)


if __name__ == "__main__":
    main()

