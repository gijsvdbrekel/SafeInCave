#!/usr/bin/env python3
import os
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mpi4py import MPI
import dolfinx as dfx

import safeincave.PostProcessingTools as post
from case_index import detect_layout_and_collect_cases, filter_cases, one_case_per_cavern_label


MPA  = 1e6
HOUR = 3600.0
DAY  = 24.0 * HOUR

CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt", "IrregularFine"]

# =============================================================================
# USER SELECTION
# =============================================================================
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/3_cavern/output"

SELECT = {
    "pressure": "sinus",          # "sinus"/"irregular"/"csv_profile"/"linear"/None
    "scenario": None,             # "full"/"desai_only"/"full_minus_desai"/"disloc_old_only"/"disloc_new_only"/None
    "caverns": None,              # e.g. ["Regular","regular600"] or None
    "n_cycles": None,
    "operation_days": None,
    "case_contains": None,
    "one_case_per_cavern": True,
}

OUTDIR = os.path.join(ROOT, "_plots")
OUTNAME = "fos_time.png"
DPI = 220

CAVERN_PHYS_TAG = 29


def build_color_map(labels):
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    return {lab: cycle[i % len(cycle)] for i, lab in enumerate(labels)}


# ---------- RD model ----------
def q_dil_rd(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    p = np.asarray(p_MPa, dtype=float)
    psi = np.asarray(psi, dtype=float)
    I1 = 3.0 * p  # MPa
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
    s1d = s1 - mean
    s2d = s2 - mean
    s3d = s3 - mean

    J2 = (1.0/6.0) * ((s1d - s2d)**2 + (s2d - s3d)**2 + (s3d - s1d)**2)
    J3 = s1d * s2d * s3d

    q_MPa = np.sqrt(3.0 * np.maximum(J2, 0.0)) / MPA

    J2_safe = np.maximum(J2, 1e-30)
    x = (3.0*np.sqrt(3.0)/2.0) * (J3 / (J2_safe**1.5))
    x = np.clip(x, -1.0, 1.0)
    theta = (1.0/3.0) * np.arccos(x)
    psi = theta - np.pi/6.0

    return p_MPa, q_MPa, psi


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
        raise FileNotFoundError(f"Missing sig.xdmf under {case_folder}")
    _, t, sig_vals = post.read_cell_tensor(sig_path)
    return np.asarray(t, float), np.asarray(sig_vals, float)


def extract_cavern_wall_cells_slices(
    geom_msh_path: str,
    cavern_phys_tag: int = 29,
    top_fraction: float = 0.20,
    bottom_fraction: float = 0.20,
):
    mesh, cell_tags, facet_tags = dfx.io.gmshio.read_from_msh(geom_msh_path, MPI.COMM_SELF, 0)
    if facet_tags is None:
        raise RuntimeError("No facet tags present in geom.msh")

    dim = mesh.topology.dim
    fdim = dim - 1

    facets = facet_tags.find(cavern_phys_tag)
    if facets.size == 0:
        raise RuntimeError(f"No facets with physical tag {cavern_phys_tag}")

    mesh.topology.create_connectivity(fdim, 0)
    f2v = mesh.topology.connectivity(fdim, 0)
    X = mesh.geometry.x

    zc = np.empty(facets.size, dtype=float)
    for i, f in enumerate(facets):
        vs = f2v.links(int(f))
        zc[i] = X[vs, 2].mean()

    z_min = float(zc.min())
    z_max = float(zc.max())
    H = z_max - z_min
    if not np.isfinite(H) or H <= 0.0:
        raise RuntimeError("Invalid cavern height from facets.")

    z_top_thr = z_max - float(np.clip(top_fraction, 0.0, 0.9)) * H
    z_bot_thr = z_min + float(np.clip(bottom_fraction, 0.0, 0.9)) * H

    roof_mask  = zc >= z_top_thr
    floor_mask = zc <= z_bot_thr
    mid_mask   = ~(roof_mask | floor_mask)

    mesh.topology.create_connectivity(fdim, dim)
    f2c = mesh.topology.connectivity(fdim, dim)

    def cells_from_mask(mask: np.ndarray) -> np.ndarray:
        subset = facets[mask]
        cells = set()
        for f in subset:
            for c in f2c.links(int(f)):
                cells.add(int(c))
        return np.array(sorted(cells), dtype=int)

    return {
        "roof":  cells_from_mask(roof_mask),
        "mid":   cells_from_mask(mid_mask),
        "floor": cells_from_mask(floor_mask),
    }


def fos_series_stats(FOS: np.ndarray, cell_idx: np.ndarray):
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


def compute_FOS(time_list, sig_vals, *, compression_positive=True, q_tol_MPa=1e-3):
    nt = sig_vals.shape[0]
    nc = sig_vals.shape[1]
    FOS = np.full((nt, nc), np.inf, dtype=float)

    for it in range(nt):
        p_MPa, q_MPa, psi = compute_p_q_psi_from_sigma(sig_vals[it], compression_positive=compression_positive)
        q_dil = q_dil_rd(p_MPa, psi)
        mask = q_MPa >= float(q_tol_MPa)
        FOS[it, mask] = q_dil[mask] / q_MPa[mask]

    FOS[~np.isfinite(FOS)] = np.inf
    FOS = np.clip(FOS, 0.0, 1e6)
    return FOS


def main():
    if MPI.COMM_WORLD.rank != 0:
        return

    all_cases = detect_layout_and_collect_cases(ROOT)
    cases = filter_cases(all_cases, SELECT)
    if not cases:
        print("[INFO] Found cases (unfiltered):")
        for m in all_cases[:40]:
            print(" -", m["case_name"], "| pressure:", m.get("pressure_scenario"),
                  "| scenario:", m.get("scenario_preset"), "| cavern:", m.get("cavern_label"))
        raise RuntimeError(f"No cases matched SELECT={SELECT}")

    if SELECT.get("one_case_per_cavern", True):
        cases = one_case_per_cavern_label(cases)

    labels_present = []
    for c in cases:
        if c["cavern_label"] not in labels_present:
            labels_present.append(c["cavern_label"])

    labels_sorted = [l for l in CAVERN_ORDER if l in labels_present] + [l for l in labels_present if l not in CAVERN_ORDER]
    color_map = build_color_map(labels_sorted)

    plt.figure(figsize=(14, 7))

    COMPRESSION_POSITIVE = True
    Q_TOL_MPA = 1e-3

    for c in cases:
        lab = c["cavern_label"]
        col = color_map.get(lab, None)

        case_path = c["case_path"]
        print(f"Processing: {c['case_name']}")

        time_list, sig_vals = load_sig(case_path)
        FOS = compute_FOS(time_list, sig_vals, compression_positive=COMPRESSION_POSITIVE, q_tol_MPa=Q_TOL_MPA)

        fos_min_global = np.nanmin(np.where(np.isfinite(FOS), FOS, np.nan), axis=1)

        geom_msh = os.path.join(case_path, "operation", "mesh", "geom.msh")
        slices = extract_cavern_wall_cells_slices(geom_msh, cavern_phys_tag=CAVERN_PHYS_TAG)

        roof_stats  = fos_series_stats(FOS, slices["roof"])
        mid_stats   = fos_series_stats(FOS, slices["mid"])
        floor_stats = fos_series_stats(FOS, slices["floor"])

        t_days = time_list / DAY

        plt.plot(t_days, fos_min_global, lw=2.2, color=col, label=f"{lab} – global min FoS")
        plt.plot(t_days, roof_stats["mean"],  ls="--", lw=2.0, label=f"{lab} – roof mean FoS")
        plt.plot(t_days, mid_stats["mean"],   ls="--", lw=2.0, label=f"{lab} – mid mean FoS")
        plt.plot(t_days, floor_stats["mean"], ls="--", lw=2.0, label=f"{lab} – floor mean FoS")

    plt.axhline(1.0, color="k", lw=1.2)
    plt.xlabel("Time (days)")
    plt.ylabel("FoS (–)")
    plt.grid(True, alpha=0.3)

    ptxt = SELECT.get("pressure", None) or "ALL"
    stxt = SELECT.get("scenario", None) or "ANY"
    plt.title(f"Factor of Safety over time | pressure={ptxt} | scenario={stxt}")

    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h
    ax.legend(uniq.values(), uniq.keys(), fontsize=8, frameon=True, loc="best")

    os.makedirs(OUTDIR, exist_ok=True)
    outpath = os.path.join(OUTDIR, OUTNAME)
    plt.tight_layout()
    plt.savefig(outpath, dpi=DPI, bbox_inches="tight")
    plt.close()
    print("[OK] Saved:", outpath)


if __name__ == "__main__":
    main()


