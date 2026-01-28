import os
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
import dolfinx as dfx
import json

import safeincave.PostProcessingTools as post

from case_index import detect_layout_and_collect_cases, filter_cases

MPA  = 1e6
HOUR = 3600.0
DAY  = 24.0 * HOUR

# =============================================================================
# USER SELECTION
# =============================================================================
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/3_cavern/output"

SELECT = {
    "pressure": "sinus",
    "caverns": ["Regular"],  # can be None or list
    "scenario": ["disloc_old_only", "disloc_new_only"],  # can be None or list
    "case_contains": None,
}

CAVERN_PHYS_TAG = 29

OUT_DIR = os.path.join(ROOT, "_figures")
SHOW = False
DPI = 180

CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt", "IrregularFine"]
SCENARIO_ORDER = ["disloc_old_only", "disloc_new_only", "desai_only", "full_minus_desai", "full", None]

def build_color_map_for_scenarios(scenarios):
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    scenarios = list(scenarios)
    return {sc: cycle[i % len(cycle)] for i, sc in enumerate(scenarios)}

def build_linestyle_map_for_caverns(cav_labels):
    styles = ["-", "--", "-.", ":"]
    cav_labels = list(cav_labels)
    return {lab: styles[i % len(styles)] for i, lab in enumerate(cav_labels)}

def pick_existing(*paths):
    for p in paths:
        if os.path.isfile(p):
            return p
    return None

def load_sig(case_folder):
    sig_path = pick_existing(
        os.path.join(case_folder, "operation", "sig", "sig.xdmf"),
        os.path.join(case_folder, "operation", "sig.xdmf"),
    )
    if sig_path is None:
        raise FileNotFoundError(f"Missing sig for {case_folder}")
    _, t, sig_vals = post.read_cell_tensor(sig_path)
    return np.asarray(t, float), np.asarray(sig_vals, float)

def q_dil_rd(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    p = np.asarray(p_MPa, dtype=float)
    psi = np.asarray(psi, dtype=float)
    I1 = 3.0 * p
    denom = (np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi))
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)
    num = D1 * (np.abs(I1) / sigma_ref) ** m + T0
    sqrtJ2_dil = num / denom
    return np.sqrt(3.0) * sqrtJ2_dil

def compute_p_q_psi_from_sigma(sig_Pa, compression_positive=True):
    sig = 0.5 * (sig_Pa + np.swapaxes(sig_Pa, -1, -2))
    sig_eff = -sig if compression_positive else sig
    vals = np.linalg.eigvalsh(sig_eff)
    s1, s2, s3 = vals[:, 2], vals[:, 1], vals[:, 0]
    I1 = s1 + s2 + s3
    p_MPa = (I1 / 3.0) / MPA
    mean = I1 / 3.0
    s1d, s2d, s3d = s1 - mean, s2 - mean, s3 - mean
    J2 = (1.0/6.0) * ((s1d - s2d)**2 + (s2d - s3d)**2 + (s3d - s1d)**2)
    J3 = s1d * s2d * s3d
    q_MPa = np.sqrt(3.0 * np.maximum(J2, 0.0)) / MPA
    J2_safe = np.maximum(J2, 1e-30)
    x = (3.0*np.sqrt(3.0)/2.0) * (J3 / (J2_safe**1.5))
    x = np.clip(x, -1.0, 1.0)
    theta = (1.0/3.0) * np.arccos(x)
    psi = theta - np.pi/6.0
    return p_MPa, q_MPa, psi

def compute_FOS(sig_vals, *, compression_positive=True, q_tol_MPa=1e-3):
    nt, nc = sig_vals.shape[0], sig_vals.shape[1]
    FOS = np.full((nt, nc), np.inf, dtype=float)
    for it in range(nt):
        p_MPa, q_MPa, psi = compute_p_q_psi_from_sigma(sig_vals[it], compression_positive=compression_positive)
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

def main():
    if MPI.COMM_WORLD.rank != 0:
        return
    os.makedirs(OUT_DIR, exist_ok=True)

    all_cases = detect_layout_and_collect_cases(ROOT)
    # keep only cases that have sig + mesh for slices
    kept = []
    for m in all_cases:
        case_path = m["case_path"]
        sig_ok = os.path.isfile(os.path.join(case_path, "operation", "sig", "sig.xdmf")) or os.path.isfile(os.path.join(case_path, "operation", "sig.xdmf"))
        msh_ok = os.path.isfile(os.path.join(case_path, "operation", "mesh", "geom.msh"))
        if sig_ok and msh_ok:
            kept.append(m)
    all_cases = kept

    cases = filter_cases(all_cases, SELECT)
    if not cases:
        print("[DEBUG] Found (first 20):")
        for m in all_cases[:20]:
            print(" -", m.get("cavern_label"), m.get("scenario_preset"), m.get("pressure_scenario"), m.get("case_name"))
        raise RuntimeError(f"No cases matched SELECT={SELECT}")

    # One case per (cavern, scenario) to avoid duplicates
    by_series = {}
    for c in cases:
        key = (c.get("cavern_label"), c.get("scenario_preset"))
        if key not in by_series:
            by_series[key] = c
    cases = list(by_series.values())

    cavs = sorted({c.get("cavern_label") for c in cases}, key=lambda x: (CAVERN_ORDER.index(x) if x in CAVERN_ORDER else 999, x))
    scs = sorted({c.get("scenario_preset") for c in cases}, key=lambda x: (SCENARIO_ORDER.index(x) if x in SCENARIO_ORDER else 999, str(x)))

    scenario_colors = build_color_map_for_scenarios(scs)
    cavern_styles = build_linestyle_map_for_caverns(cavs)

    plt.figure(figsize=(14, 7))

    COMPRESSION_POSITIVE = True
    Q_TOL_MPA = 1e-3

    for c in cases:
        cav = c.get("cavern_label")
        sc = c.get("scenario_preset")
        col = scenario_colors.get(sc, "C0")
        ls = cavern_styles.get(cav, "-")
        label = f"{cav} | {sc}" if sc is not None else f"{cav} | (no scenario)"

        time_list, sig_vals = load_sig(c["case_path"])
        FOS = compute_FOS(sig_vals, compression_positive=COMPRESSION_POSITIVE, q_tol_MPa=Q_TOL_MPA)

        fos_min_global = np.nanmin(np.where(np.isfinite(FOS), FOS, np.nan), axis=1)

        geom_msh = os.path.join(c["case_path"], "operation", "mesh", "geom.msh")
        slices = extract_cavern_wall_cells_slices(geom_msh, cavern_phys_tag=CAVERN_PHYS_TAG, top_fraction=0.20, bottom_fraction=0.20)

        roof_stats  = fos_series_stats(FOS, slices["roof"])
        mid_stats   = fos_series_stats(FOS, slices["mid"])
        floor_stats = fos_series_stats(FOS, slices["floor"])

        t_days = time_list / DAY

        plt.plot(t_days, fos_min_global, lw=2.2, color=col, linestyle=ls, label=f"{label} — global min")
        # keep means thinner; same color (scenario) but different alpha
        plt.plot(t_days, roof_stats["mean"],  lw=1.6, color=col, linestyle="--", alpha=0.9, label=f"{label} — roof mean")
        plt.plot(t_days, mid_stats["mean"],   lw=1.6, color=col, linestyle="-.", alpha=0.9, label=f"{label} — mid mean")
        plt.plot(t_days, floor_stats["mean"], lw=1.6, color=col, linestyle=":",  alpha=0.9, label=f"{label} — floor mean")

    plt.axhline(1.0, color="k", lw=1.2)
    plt.xlabel("Time (days)")
    plt.ylabel("FoS (–)")
    plt.grid(True, alpha=0.3)
    plt.title(f"Factor of Safety over time | pressure={SELECT.get('pressure')}")

    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h
    ax.legend(uniq.values(), uniq.keys(), fontsize=8, frameon=True, loc="best")

    plt.tight_layout()

    outname = f"fos_pressure={SELECT.get('pressure')}_scenario={SELECT.get('scenario')}.png"
    outpath = os.path.join(OUT_DIR, outname.replace(" ", ""))
    plt.savefig(outpath, dpi=DPI)
    print("[SAVED]", outpath)

    if SHOW:
        plt.show()
    plt.close()

if __name__ == "__main__":
    main()


