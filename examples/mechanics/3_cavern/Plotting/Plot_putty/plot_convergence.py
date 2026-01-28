#!/usr/bin/env python3
import os
import json
import numpy as np

# ---- headless plotting (NO GUI needed) ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import safeincave.PostProcessingTools as post

from mpi4py import MPI
import dolfinx as dfx


# Units
DAY = 24.0 * 3600.0
MPA = 1e6

# Label ordering (optional, for nice legends)
CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt", "IrregularFine"]

# Gmsh physical tag for cavern surface
CAVERN_PHYS_TAG = 29

DEBUG_ONE_CASE = False
DEBUG_CASE_PATH = ""


# =============================================================================
# USER SELECTION (edit only this block)
# =============================================================================
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/3_cavern/output"

SELECT = {
    # In this output layout there are no cavern group folders.
    # We infer cavern type from the case folder name (last token), e.g. *_regular600.
    # Use None for all, or a list like ["regular600", "tilted600", "irregular1200"].
    "caverns": None,

    # Pressure scheme from pressure_schedule.json: "sinus", "irregular", "csv_profile", "linear", ...
    # None = all
    "pressure": "sinus",

    # Optional: only keep case folders containing substring
    "case_contains": None,  # e.g. "365days" or "desai_only" or None

    # Plot settings
    "separate_per_pressure": True,
    "cavern_phys_tag": CAVERN_PHYS_TAG,
}

# Saving (no plt.show())
SAVE_DIR_NAME = "_convergence_outputs"  # saved under ROOT
SAVE_PNG = True
SAVE_PDF = True
DPI = 200


# =============================================================================
# Naming helpers (adapted to flat case folder layout)
# =============================================================================
def infer_cavern_type_from_case_name(case_name: str) -> str:
    # case_disloc_old_only_sinus_8cyc_365days_regular600 -> "regular600"
    parts = case_name.split("_")
    if len(parts) < 2:
        return case_name.lower()
    return parts[-1].lower()


def cavern_label_from_cavern_type(cavern_type: str) -> str:
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
    if low.startswith("irregularfine") or low.startswith("irregular_fine"):
        return "IrregularFine"
    if low.startswith("irregular"):
        return "Irregular"
    return cavern_type


def build_color_map(labels):
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    return {lab: cycle[i % len(cycle)] for i, lab in enumerate(labels)}


# =============================================================================
# Folder structure helpers (flat run output)
# ROOT/case_*/operation/u/u.xdmf
# ROOT/case_*/operation/u/u.h5
# ROOT/case_*/operation/mesh/geom.msh
# ROOT/case_*/pressure_schedule.json
# =============================================================================
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


# =============================================================================
# Pressure schedule reader
# =============================================================================
def read_pressure_schedule(case_folder: str):
    """
    Returns (t_hours, p_MPa) or (None, None).
    Supports:
      - new: t_hours + p_MPa
      - new: t_values_s + p_values_Pa
      - old: t_values + p_values   (s/Pa)
    """
    pjson = path_pressure_json(case_folder)
    if not os.path.isfile(pjson):
        return None, None

    with open(pjson, "r") as f:
        data = json.load(f)

    if "t_hours" in data and "p_MPa" in data:
        t = np.asarray(data["t_hours"], dtype=float)
        p = np.asarray(data["p_MPa"], dtype=float)
        return t, p

    if "t_values_s" in data and "p_values_Pa" in data:
        t = np.asarray(data["t_values_s"], dtype=float) / 3600.0
        p = np.asarray(data["p_values_Pa"], dtype=float) / 1e6
        return t, p

    if "t_values" in data and "p_values" in data:
        t = np.asarray(data["t_values"], dtype=float) / 3600.0
        p = np.asarray(data["p_values"], dtype=float) / 1e6
        return t, p

    return None, None


def read_pressure_from_case(case_path: str):
    t_hours, p_mpa = read_pressure_schedule(case_path)
    if t_hours is None:
        return None, None
    return t_hours / 24.0, p_mpa


def read_case_pressure_scenario(case_path: str) -> str | None:
    """
    Returns pressure scenario string ("sinus","irregular","csv_profile",...)
    from pressure_schedule.json.
    """
    pjson = path_pressure_json(case_path)
    if not os.path.isfile(pjson):
        return None

    try:
        with open(pjson, "r") as f:
            data = json.load(f)
    except Exception:
        return None

    # preferred key in some scripts
    if isinstance(data, dict) and "pressure_scenario" in data:
        val = data.get("pressure_scenario")
        return str(val).lower() if val is not None else None

    # your Run.py stores scheme in "scenario"
    if isinstance(data, dict) and "scenario" in data:
        val = str(data.get("scenario")).lower()
        if val in ("sinus", "irregular", "csv_profile", "linear"):
            return val

    return None


# =============================================================================
# Geometry helpers: cavern facets + area vectors
# =============================================================================
def _area_vectors(points_xyz: np.ndarray, tri: np.ndarray):
    """
    Compute area vector (n * area) for each triangle.
    tri: (N,3) vertex indices
    Returns: (N,3) area vectors.
    """
    p0 = points_xyz[tri[:, 0]]
    p1 = points_xyz[tri[:, 1]]
    p2 = points_xyz[tri[:, 2]]
    return 0.5 * np.cross(p1 - p0, p2 - p0)


def extract_cavern_facets_from_msh_dolfinx(geom_msh_path: str, cavern_tag: int = 29):
    """
    Returns:
      X: (Nnodes,3) mesh node coordinates (float)
      facets: (Nfacets,3) vertex indices (int) of triangular facets on the cavern boundary
    """
    # Use COMM_SELF so this works in normal serial plotting runs
    mesh, cell_tags, facet_tags = dfx.io.gmshio.read_from_msh(geom_msh_path, MPI.COMM_SELF, 0)

    if facet_tags is None:
        raise RuntimeError("No facet tags were read from the .msh (facet_tags is None).")

    fdim = mesh.topology.dim - 1
    cavern_facets = facet_tags.find(cavern_tag)
    if cavern_facets.size == 0:
        raise RuntimeError(f"No facets found with physical tag {cavern_tag} in {geom_msh_path}")

    mesh.topology.create_connectivity(fdim, 0)
    f2v = mesh.topology.connectivity(fdim, 0)

    facets = []
    for f in cavern_facets:
        verts = f2v.links(int(f))
        if len(verts) != 3:
            continue
        facets.append(verts)

    facets = np.asarray(facets, dtype=int)
    if facets.size == 0:
        raise RuntimeError(f"Facets for tag {cavern_tag} exist, but none had 3 vertices (unexpected).")

    X = mesh.geometry.x.copy()
    return X, facets


def orient_area_vectors_outward(points: np.ndarray, tris: np.ndarray, area_vecs: np.ndarray):
    """
    Ensure area vectors are consistently oriented outward from the cavern interior.
    Heuristic: use cavern surface centroid as 'center', and flip triangles whose
    normal points toward the center (dot(nA, centroid - center) < 0).
    """
    center = np.mean(points[np.unique(tris.reshape(-1))], axis=0)
    tri_centroids = (points[tris[:, 0]] + points[tris[:, 1]] + points[tris[:, 2]]) / 3.0
    v = tri_centroids - center
    s = np.einsum("ij,ij->i", area_vecs, v)
    flip = s < 0.0
    area_vecs2 = area_vecs.copy()
    area_vecs2[flip] *= -1.0
    return area_vecs2


# =============================================================================
# Convergence computation
# =============================================================================
def compute_convergence_3d_percent(case_path: str, cavern_phys_tag: int = 29):
    """
    3D cavern volume convergence during OPERATION stage, referenced to t=0 of u series.

    Small-strain approximation:
        ΔV(t) ≈ ∫_Γ (u(t) - u(0)) · n dS
    with V0 from same surface:
        V0 = (1/3) ∫_Γ x · n dS

    Returns:
      t_days: (Nt,)
      conv_pct: (Nt,) with conv_pct[0] = 0
    """
    geom_msh = path_geom_msh(case_path)
    u_xdmf = path_u_xdmf(case_path)

    pts_msh, tris_msh = extract_cavern_facets_from_msh_dolfinx(
        geom_msh, cavern_tag=cavern_phys_tag
    )

    area_vecs = _area_vectors(pts_msh, tris_msh)
    area_vecs = orient_area_vectors_outward(pts_msh, tris_msh, area_vecs)

    tri_centroids = (pts_msh[tris_msh[:, 0]] +
                     pts_msh[tris_msh[:, 1]] +
                     pts_msh[tris_msh[:, 2]]) / 3.0
    V0 = (1.0 / 3.0) * np.sum(np.einsum("ij,ij->i", tri_centroids, area_vecs))
    V0 = float(abs(V0))
    if not np.isfinite(V0) or V0 <= 0.0:
        raise RuntimeError(f"Computed invalid V0={V0} from cavern surface in {geom_msh}")

    points_xdmf, time_list, u_field = post.read_node_vector(u_xdmf)  # (Nt,N,3)
    if u_field.ndim != 3 or u_field.shape[2] != 3:
        raise RuntimeError(f"Unexpected u_field shape {u_field.shape} from {u_xdmf}")

    Nt = u_field.shape[0]
    if Nt < 1:
        raise RuntimeError(f"No timesteps found in {u_xdmf}")

    # Map dolfinx mesh nodes -> xdmf nodes
    mapping = post.build_mapping(pts_msh, points_xdmf)

    tri_xdmf = np.vectorize(mapping.__getitem__)(tris_msh).astype(int)
    v0 = tri_xdmf[:, 0]
    v1 = tri_xdmf[:, 1]
    v2 = tri_xdmf[:, 2]

    u_ref = u_field[0]

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
# Index cases: ROOT/case_*/...
# =============================================================================
def collect_cases_all(ROOT: str):
    """
    Collect ALL cases under flat output folder.
    Returns list of dicts:
      {label, cavern_type, case_name, case_path, pressure_scenario}
    """
    cases = []
    for sub in sorted(os.listdir(ROOT)):
        if not sub.lower().startswith("case_"):
            continue

        case_path = os.path.join(ROOT, sub)
        if not os.path.isdir(case_path):
            continue

        # required for 3D convergence method
        u_xdmf = path_u_xdmf(case_path)
        u_h5 = path_u_h5(case_path)
        msh = path_geom_msh(case_path)
        if not (os.path.isfile(u_xdmf) and os.path.isfile(u_h5) and os.path.isfile(msh)):
            continue
        if not is_valid_hdf5(u_h5):
            continue

        cavern_type = infer_cavern_type_from_case_name(sub)
        label = cavern_label_from_cavern_type(cavern_type)
        pres = read_case_pressure_scenario(case_path)

        cases.append({
            "label": label,
            "cavern_type": cavern_type,
            "case_name": sub,
            "case_path": case_path,
            "pressure_scenario": pres,
            "has_pressure_json": os.path.isfile(path_pressure_json(case_path)),
        })
    return cases


def filter_cases(cases, select: dict):
    cavs = select.get("caverns", None)
    pres = select.get("pressure", None)
    contains = select.get("case_contains", None)

    out = []
    for c in cases:
        if cavs is not None:
            cavs_l = [x.lower() for x in cavs]
            if c["cavern_type"].lower() not in cavs_l and c["label"].lower() not in cavs_l:
                continue

        if pres is not None:
            if (c["pressure_scenario"] or "").lower() != pres.lower():
                continue

        if contains is not None and contains.lower() not in c["case_name"].lower():
            continue

        out.append(c)
    return out


# =============================================================================
# Plotting
# =============================================================================
def plot_scheme_from_case_list(cases, *, title, color_map, cavern_phys_tag=CAVERN_PHYS_TAG):
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(13, 7), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12},
    )
    fig.suptitle(title, fontsize=12)

    plotted_any = False
    for c in cases:
        lab = c["label"]
        col = color_map.get(lab, None)

        try:
            t_days, conv = compute_convergence_3d_percent(c["case_path"], cavern_phys_tag=cavern_phys_tag)
        except Exception as e:
            print(f"[SKIP] convergence {c['case_name']}: {e}")
            continue

        ax1.plot(t_days, conv, linewidth=2.0, alpha=0.95, color=col, label=lab)
        plotted_any = True

    if not plotted_any:
        ax1.text(0.5, 0.5, "No convergence curves could be plotted (all cases failed).",
                 ha="center", va="center", transform=ax1.transAxes)

    ax1.set_ylabel("Convergence (ΔV/V0) (%)")
    ax1.grid(True, alpha=0.3)

    # unique legend
    handles, labs = ax1.get_legend_handles_labels()
    uniq = {}
    for h, l in zip(handles, labs):
        uniq[l] = h
    if uniq:
        ax1.legend(uniq.values(), uniq.keys(), loc="best", fontsize=9, frameon=True)

    # --- pressure schedule: take first available in this selection ---
    ref_t, ref_p = None, None
    for c in cases:
        t, p = read_pressure_from_case(c["case_path"])
        if t is not None:
            ref_t, ref_p = t, p
            break

    if ref_t is None:
        ax2.text(0.5, 0.5, "No pressure_schedule.json found in selected cases.",
                 ha="center", va="center", transform=ax2.transAxes)
    else:
        ax2.plot(ref_t, ref_p, linewidth=1.7)

    ax2.set_ylabel("Pressure (MPa)")
    ax2.set_xlabel("Time (days)")
    ax2.grid(True, alpha=0.3)

    return fig


def save_figure(fig, base_name: str):
    out_dir = os.path.join(ROOT, SAVE_DIR_NAME)
    os.makedirs(out_dir, exist_ok=True)

    saved = []
    if SAVE_PNG:
        p = os.path.join(out_dir, base_name + ".png")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
        saved.append(p)
    if SAVE_PDF:
        p = os.path.join(out_dir, base_name + ".pdf")
        fig.savefig(p, bbox_inches="tight")
        saved.append(p)
    return saved


def main():
    if MPI.COMM_WORLD.rank != 0:
        return

    if DEBUG_ONE_CASE:
        if not DEBUG_CASE_PATH:
            raise RuntimeError("Set DEBUG_CASE_PATH when DEBUG_ONE_CASE=True")

        X, facets = extract_cavern_facets_from_msh_dolfinx(
            path_geom_msh(DEBUG_CASE_PATH),
            cavern_tag=CAVERN_PHYS_TAG
        )
        print("DEBUG cavern facets:", facets.shape[0])
        print("DEBUG unique cavern vertices:", len(np.unique(facets.reshape(-1))))
        return

    # 1) Index everything once
    all_cases = collect_cases_all(ROOT)
    if not all_cases:
        raise RuntimeError(f"No cases found under {ROOT}")

    # 2) Apply selection filters
    selected = filter_cases(all_cases, SELECT)
    if not selected:
        raise RuntimeError(f"No cases match SELECT filters: {SELECT}")

    # 3) Build global color map across selected caverns
    labels_present = []
    for c in selected:
        if c["label"] not in labels_present:
            labels_present.append(c["label"])

    labels_sorted = [l for l in CAVERN_ORDER if l in labels_present] + \
                    [l for l in labels_present if l not in CAVERN_ORDER]
    global_cmap = build_color_map(labels_sorted)

    # 4) Decide grouping
    saved_all = []

    if SELECT.get("separate_per_pressure", False):
        groups = {}
        for c in selected:
            key = c["pressure_scenario"] or "unknown_pressure"
            groups.setdefault(key, []).append(c)

        for pres_key, cases_here in sorted(groups.items()):
            fig = plot_scheme_from_case_list(
                cases_here,
                title=f"3D cavern volume convergence – pressure: {pres_key}",
                color_map=global_cmap,
                cavern_phys_tag=SELECT.get("cavern_phys_tag", CAVERN_PHYS_TAG),
            )
            base = f"convergence3d_pressure={pres_key}"
            if SELECT.get("caverns") is not None:
                base += "_cav=" + "-".join([x.lower() for x in SELECT["caverns"]])
            if SELECT.get("case_contains") is not None:
                base += "_contains=" + str(SELECT["case_contains"]).lower()

            saved = save_figure(fig, base)
            saved_all.extend(saved)
            plt.close(fig)

    else:
        fig = plot_scheme_from_case_list(
            selected,
            title="3D cavern volume convergence – selected cases",
            color_map=global_cmap,
            cavern_phys_tag=SELECT.get("cavern_phys_tag", CAVERN_PHYS_TAG),
        )
        base = "convergence3d_selected"
        saved = save_figure(fig, base)
        saved_all.extend(saved)
        plt.close(fig)

    print("\n[OK] Saved:")
    for p in saved_all:
        print(" -", p)


if __name__ == "__main__":
    main()
