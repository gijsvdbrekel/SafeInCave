#!/usr/bin/env python3
import os
import json
import numpy as np

# ---------- headless plotting (NO GUI needed) ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mpi4py import MPI
import dolfinx as dfx

import safeincave.PostProcessingTools as post


MPA  = 1e6
HOUR = 3600.0
DAY  = 24.0 * HOUR

# ---------- consistent naming + colors ----------
CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt"]

# === NEW ROOT: folder that directly contains case_* folders ===
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/3_cavern/output"

SELECT = {
    "pressure": "sinus",          # "sinus" / "irregular" / "csv_profile" / "linear" / None (=all)
    # In this output layout there is NO group folder; caverns are inferred from case name (e.g. *_regular600)
    "caverns": ["regular600"],    # e.g. ["regular600", "tilted600"] or None (=all)
    "case_contains": None,        # e.g. "365days" or "8cyc" or "disloc_old_only" or None
}

CAVERN_PHYS_TAG = 29

# ---------- saving ----------
SAVE_DIR_NAME = "_fos_outputs"     # created under ROOT
SAVE_PNG = True
SAVE_PDF = True
DPI = 200


# =============================================================================
# Label helpers (adapted to case folder names like:
# case_disloc_old_only_sinus_8cyc_365days_regular600
# =============================================================================
def infer_cavern_from_case_name(case_name: str) -> str:
    parts = case_name.split("_")
    if len(parts) < 2:
        return case_name
    return parts[-1].lower()  # e.g. "regular600"


def nice_cavern_label(cavern_type: str) -> str:
    low = cavern_type.lower()
    if low.startswith("asymmetric"):
        return "Asymmetric"
    if low.startswith("multichamber"):
        return "Multichamber"
    if low.startswith("teardrop"):
        return "Teardrop"
    if low.startswith("tilt") or low.startswith("tilted"):
        return "Tilt"
    if low.startswith("regular"):
        return "Regular"
    if low.startswith("irregular"):
        return "Irregular"
    return cavern_type


def build_color_map(labels):
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    return {lab: cycle[i % len(cycle)] for i, lab in enumerate(labels)}


# =============================================================================
# RD model
# =============================================================================
def q_dil_rd(p_MPa, psi, D1=0.683, D2=0.512, m=0.75, T0=1.5, sigma_ref=1.0):
    """
    RD dilation boundary q_dil(p, psi) in MPa.
    p_MPa: mean stress p (compression-positive) in MPa
    psi  : Lode angle in [-pi/6, +pi/6]
    """
    p = np.asarray(p_MPa, dtype=float)
    psi = np.asarray(psi, dtype=float)

    I1 = 3.0 * p  # MPa

    denom = (np.sqrt(3.0) * np.cos(psi) - D2 * np.sin(psi))
    denom = np.where(np.abs(denom) < 1e-12, np.sign(denom) * 1e-12, denom)

    num = D1 * (np.abs(I1) / sigma_ref) ** m + T0
    sqrtJ2_dil = num / denom
    q_dil = np.sqrt(3.0) * sqrtJ2_dil
    return q_dil


# =============================================================================
# invariants from sigma
# =============================================================================
def compute_p_q_psi_from_sigma(sig_Pa, compression_positive=True):
    """
    sig_Pa: (ncells, 3, 3) in Pa.
    Returns: p in MPa, q in MPa, psi in [-pi/6, +pi/6]
    """
    sig = 0.5 * (sig_Pa + np.swapaxes(sig_Pa, -1, -2))  # enforce symmetry
    sig_eff = -sig if compression_positive else sig

    vals = np.linalg.eigvalsh(sig_eff)
    s1, s2, s3 = vals[:, 2], vals[:, 1], vals[:, 0]     # s1 >= s2 >= s3

    I1 = s1 + s2 + s3                                   # Pa
    p_MPa = (I1 / 3.0) / MPA                            # MPa

    mean = I1 / 3.0
    s1d = s1 - mean
    s2d = s2 - mean
    s3d = s3 - mean

    J2 = (1.0/6.0) * ((s1d - s2d)**2 + (s2d - s3d)**2 + (s3d - s1d)**2)  # Pa^2
    J3 = s1d * s2d * s3d                                                 # Pa^3

    q_MPa = np.sqrt(3.0 * np.maximum(J2, 0.0)) / MPA

    J2_safe = np.maximum(J2, 1e-30)
    x = (3.0*np.sqrt(3.0)/2.0) * (J3 / (J2_safe**1.5))
    x = np.clip(x, -1.0, 1.0)
    theta = (1.0/3.0) * np.arccos(x)
    psi = theta - np.pi/6.0

    return p_MPa, q_MPa, psi


# =============================================================================
# robust path helpers
# =============================================================================
def pick_existing(*paths):
    for p in paths:
        if os.path.isfile(p):
            return p
    return None


def load_sig(case_folder):
    """
    Run output layout:
      case_folder/operation/sig/sig.xdmf
      case_folder/operation/sig.xdmf          (fallback)
    """
    sig_path = pick_existing(
        os.path.join(case_folder, "operation", "sig", "sig.xdmf"),
        os.path.join(case_folder, "operation", "sig.xdmf"),
    )
    if sig_path is None:
        raise FileNotFoundError(
            f"Missing sig in {case_folder} (checked operation/sig/sig.xdmf and operation/sig.xdmf)"
        )

    _, t, sig_vals = post.read_cell_tensor(sig_path)
    return np.asarray(t, float), np.asarray(sig_vals, float), sig_path


# =============================================================================
# extract cavern wall cells (from geom.msh)
# =============================================================================
def extract_cavern_wall_cells_slices(
    geom_msh_path: str,
    cavern_phys_tag: int = 29,
    mode: str = "fraction",
    top_fraction: float = 0.20,
    bottom_fraction: float = 0.20,
    top_thickness_m: float | None = None,
    bottom_thickness_m: float | None = None,
):
    """
    Returns dict with cell indices for {'roof','mid','floor'} plus z-bounds.
    Slices are defined by facet-centroid z and converted to adjacent volume cells.
    """
    mesh, cell_tags, facet_tags = dfx.io.gmshio.read_from_msh(
        geom_msh_path, MPI.COMM_SELF, 0
    )

    if facet_tags is None:
        raise RuntimeError("No facet tags present in geom.msh")

    dim = mesh.topology.dim
    fdim = dim - 1

    facets = facet_tags.find(cavern_phys_tag)
    if facets.size == 0:
        raise RuntimeError(f"No facets with physical tag {cavern_phys_tag} in {geom_msh_path}")

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
        raise RuntimeError(f"Invalid cavern height from facets: H={H}")

    if mode == "fraction":
        top_fraction = float(np.clip(top_fraction, 0.0, 0.9))
        bottom_fraction = float(np.clip(bottom_fraction, 0.0, 0.9))
        z_top_thr = z_max - top_fraction * H
        z_bot_thr = z_min + bottom_fraction * H
    elif mode == "thickness_m":
        if top_thickness_m is None or bottom_thickness_m is None:
            raise ValueError("Provide top_thickness_m and bottom_thickness_m for mode='thickness_m'")
        z_top_thr = z_max - float(top_thickness_m)
        z_bot_thr = z_min + float(bottom_thickness_m)
    else:
        raise ValueError("mode must be 'fraction' or 'thickness_m'")

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
        "z_bounds": (z_min, z_max),
    }


# =============================================================================
# FoS computations
# =============================================================================
def fos_series_stats(FOS: np.ndarray, cell_idx: np.ndarray):
    """
    FOS: (Nt, Ncells)
    cell_idx: indices into cells (1D)
    Returns dict of time series (Nt,) for min, mean, and 5th percentile.
    """
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
    """
    FoS per time step, per cell:
        FoS = q_dil(p, psi) / q
    q<tol -> FoS=+inf
    """
    nt = sig_vals.shape[0]
    nc = sig_vals.shape[1]
    FOS = np.full((nt, nc), np.inf, dtype=float)

    for it in range(nt):
        p_MPa, q_MPa, psi = compute_p_q_psi_from_sigma(
            sig_vals[it],
            compression_positive=compression_positive
        )
        q_dil = q_dil_rd(p_MPa, psi)
        mask = q_MPa >= float(q_tol_MPa)
        FOS[it, mask] = q_dil[mask] / q_MPa[mask]

    FOS[~np.isfinite(FOS)] = np.inf
    FOS = np.clip(FOS, 0.0, 1e6)
    return FOS


# =============================================================================
# Case indexing for RUN output layout: ROOT/case_*/...
# =============================================================================
def read_case_pressure_scenario(case_path: str):
    """
    In your runs you often store the pressure scheme under:
      pressure_schedule.json["scenario"]  (e.g. "sinus")
    or sometimes under "pressure_scenario".
    """
    pjson = os.path.join(case_path, "pressure_schedule.json")
    if not os.path.isfile(pjson):
        return None
    try:
        with open(pjson, "r") as f:
            data = json.load(f)
    except Exception:
        return None

    if isinstance(data, dict) and "pressure_scenario" in data:
        v = data.get("pressure_scenario")
        return str(v).lower() if v is not None else None

    if isinstance(data, dict) and "scenario" in data:
        v = str(data.get("scenario")).lower()
        if v in ("sinus", "irregular", "csv_profile", "linear"):
            return v

    return None


def collect_cases_flat(ROOT: str, target_pressure: str | None):
    """
    ROOT/case_*/...
    Filter by pressure scenario using pressure_schedule.json.
    """
    cases = []
    for sub in sorted(os.listdir(ROOT)):
        if not sub.lower().startswith("case_"):
            continue

        case_path = os.path.join(ROOT, sub)
        if not os.path.isdir(case_path):
            continue

        pres = read_case_pressure_scenario(case_path)
        if target_pressure is not None:
            if (pres or "").lower() != target_pressure.lower():
                continue

        # require sig and geom
        try:
            load_sig(case_path)
        except Exception:
            continue

        geom_msh = os.path.join(case_path, "operation", "mesh", "geom.msh")
        if not os.path.isfile(geom_msh):
            continue

        cavern_type = infer_cavern_from_case_name(sub)
        label = nice_cavern_label(cavern_type)

        cases.append({
            "label": label,
            "cavern_type": cavern_type,
            "case_name": sub,
            "case_path": case_path,
            "pressure_scenario": pres,
        })

    return cases


# =============================================================================
# Main
# =============================================================================
def main():
    # plotting is single-rank work; avoid multi-rank matplotlib chaos
    if MPI.COMM_WORLD.rank != 0:
        return

    target_pressure = SELECT.get("pressure", None)
    COMPRESSION_POSITIVE = True
    Q_TOL_MPA = 1e-3  # 1 kPa

    cases = collect_cases_flat(ROOT, target_pressure)
    if not cases:
        raise RuntimeError(f"No cases found for pressure='{target_pressure}' with sig under:\n  {ROOT}")

    cav_filter = SELECT.get("caverns", None)
    contains = SELECT.get("case_contains", None)

    filtered = []
    for c in cases:
        if cav_filter is not None:
            cav_filter_l = [x.lower() for x in cav_filter]
            if c["cavern_type"].lower() not in cav_filter_l and c["label"].lower() not in cav_filter_l:
                continue
        if contains is not None and contains.lower() not in c["case_name"].lower():
            continue
        filtered.append(c)

    if not filtered:
        raise RuntimeError(f"No cases left after filters SELECT={SELECT}")

    # Colors by label present
    labels_present = []
    for c in filtered:
        if c["label"] not in labels_present:
            labels_present.append(c["label"])

    labels_sorted = [l for l in CAVERN_ORDER if l in labels_present] + \
                    [l for l in labels_present if l not in CAVERN_ORDER]
    color_map = build_color_map(labels_sorted)

    out_tag = (target_pressure if target_pressure is not None else "ALL")

    # one case per cavern label (to avoid legend spam)
    ONE_CASE_PER_CAVERN = True
    if ONE_CASE_PER_CAVERN:
        by_label = {}
        for c in filtered:
            if c["label"] not in by_label:
                by_label[c["label"]] = c
        filtered = list(by_label.values())

    plt.figure(figsize=(14, 7))

    for c in filtered:
        lab = c["label"]
        col = color_map.get(lab, None)

        print(f"\n=== {c['case_name']} ===")
        time_list, sig_vals, _mesh_source = load_sig(c["case_path"])

        FOS = compute_FOS(
            time_list, sig_vals,
            compression_positive=COMPRESSION_POSITIVE,
            q_tol_MPa=Q_TOL_MPA
        )

        # Global min across all cells (finite-only)
        fos_min_global = np.nanmin(np.where(np.isfinite(FOS), FOS, np.nan), axis=1)

        # Cavern slices
        geom_msh = os.path.join(c["case_path"], "operation", "mesh", "geom.msh")
        slices = extract_cavern_wall_cells_slices(
            geom_msh,
            cavern_phys_tag=CAVERN_PHYS_TAG,
            mode="fraction",
            top_fraction=0.20,
            bottom_fraction=0.20
        )

        # Sanity: indices in range
        ncells = FOS.shape[1]
        for key in ("roof", "mid", "floor"):
            idx = slices[key]
            if idx.size and (idx.min() < 0 or idx.max() >= ncells):
                raise RuntimeError(f"Slice '{key}' has cell index out of range for {lab} ({c['case_name']}).")

        roof_stats  = fos_series_stats(FOS, slices["roof"])
        mid_stats   = fos_series_stats(FOS, slices["mid"])
        floor_stats = fos_series_stats(FOS, slices["floor"])

        t_days = time_list / DAY

        # Plot global min (color = cavern)
        plt.plot(t_days, fos_min_global, lw=2.2, color=col, label=f"{lab} – global min FoS")

        # Plot means with fixed line styles (so roof/mid/floor are comparable)
        plt.plot(t_days, roof_stats["mean"],  ls="--", lw=2.0, alpha=0.9, label=f"{lab} – roof mean FoS")
        plt.plot(t_days, mid_stats["mean"],   ls="--", lw=2.0, alpha=0.7, label=f"{lab} – mid mean FoS")
        plt.plot(t_days, floor_stats["mean"], ls="--", lw=2.0, alpha=0.9, label=f"{lab} – floor mean FoS")

        # Optional conservative lines (5th percentile)
        plt.plot(t_days, roof_stats["p05"],  ls="-.", lw=1.2, color=col, alpha=0.9)
        plt.plot(t_days, mid_stats["p05"],   ls="-.", lw=1.2, color=col, alpha=0.6)
        plt.plot(t_days, floor_stats["p05"], ls="-.", lw=1.2, color=col, alpha=0.9)

    plt.axhline(1.0, color="k", lw=1.2)
    plt.xlabel("Time (days)")
    plt.ylabel("FoS (–)")
    plt.grid(True, alpha=0.3)
    plt.title(f"Factor of Safety over time (pressure = {out_tag})")

    # Deduplicate legend entries
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labels):
        if l not in uniq:
            uniq[l] = h
    ax.legend(uniq.values(), uniq.keys(), fontsize=9, frameon=True, loc="best")

    plt.tight_layout()

    # --- Save instead of show ---
    out_dir = os.path.join(ROOT, SAVE_DIR_NAME)
    os.makedirs(out_dir, exist_ok=True)

    tag_bits = [f"pressure={out_tag}"]
    if SELECT.get("caverns") is not None:
        tag_bits.append("cav=" + "-".join([x.lower() for x in SELECT["caverns"]]))
    if SELECT.get("case_contains") is not None:
        tag_bits.append("contains=" + str(SELECT["case_contains"]).lower())

    base = "fos_time_" + "_".join(tag_bits).replace("/", "_").replace(" ", "")

    saved = []
    if SAVE_PNG:
        png_path = os.path.join(out_dir, base + ".png")
        plt.savefig(png_path, dpi=DPI, bbox_inches="tight")
        saved.append(png_path)
    if SAVE_PDF:
        pdf_path = os.path.join(out_dir, base + ".pdf")
        plt.savefig(pdf_path, bbox_inches="tight")
        saved.append(pdf_path)

    plt.close()

    print("\n[OK] Saved:")
    for p in saved:
        print(" -", p)


if __name__ == "__main__":
    main()

