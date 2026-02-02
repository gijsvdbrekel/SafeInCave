import os
import json
import numpy as np
import matplotlib.pyplot as plt
import meshio

import safeincave.PostProcessingTools as post



# Units
DAY = 24.0 * 3600.0
MPA = 1e6

# Label ordering (optional, for nice legends)
CAVERN_ORDER = ["Asymmetric", "Irregular", "Multichamber", "Regular", "Teardrop", "Tilt", "IrregularFine"]

# Gmsh physical tag for cavern surface
CAVERN_PHYS_TAG = 29

DEBUG_ONE_CASE = False  # zet True om 1 case te checken
DEBUG_CASE_PATH = ""    # vul in als DEBUG_ONE_CASE=True



# =============================================================================
# USER SELECTION (edit only this block later)
# =============================================================================
ROOT = r"/home/gvandenbrekel/SafeInCave/OutputNobian"

SELECT = {
    # Cavern folders inside OutputNobian: "Regular", "Irregular", "Tilt", ...
    # None = all caverns
    "caverns": ["Regular"],                 # e.g. ["Regular"] or ["Regular","Irregular"] or None = all

    # Pressure scheme from pressure_schedule.json: "sinus", "irregular", "csv_profile", "linear", ...
    # None = all pressure schemes
    "pressure": None,                # e.g. "sinus"

    # Optional: only keep case folders containing substring
    "case_contains": None,           # e.g. "365days" or "desai_only" or None

    # Plot settings
    "separate_per_pressure": True,   # True => 1 figure per pressure scheme found
    "cavern_phys_tag": CAVERN_PHYS_TAG,
}



# ------------------------
# Naming helpers (like your stress-state script)
# ------------------------
def cavern_label_from_group(group_folder: str) -> str:
    low = group_folder.lower()
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
    return group_folder.split("_")[0]





def build_color_map(labels):
    cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    if not cycle:
        cycle = [f"C{i}" for i in range(10)]
    return {lab: cycle[i % len(cycle)] for i, lab in enumerate(labels)}


# ------------------------
# Folder structure helpers
# ROOT/<GROUP>/<CASE>/operation/u/u.xdmf
# ROOT/<GROUP>/<CASE>/operation/u/u.h5
# ROOT/<GROUP>/<CASE>/operation/mesh/geom.msh
# ROOT/<GROUP>/<CASE>/pressure_schedule.json
# ------------------------
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







def read_pressure_from_case(case_path: str):
    t_hours, p_mpa = read_pressure_schedule(case_path)
    if t_hours is None:
        return None, None
    return t_hours / 24.0, p_mpa



def read_pressure_schedule(case_folder: str):
    """
    Returns (t_hours, p_MPa) or (None, None).
    Supports:
      - new: t_hours + p_MPa
      - new: t_values_s + p_values_Pa
      - old: t_values + p_values   (s/Pa)
    """
    pjson = os.path.join(case_folder, "pressure_schedule.json")
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

    # unknown format
    return None, None




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




from mpi4py import MPI
import dolfinx as dfx

def extract_cavern_facets_from_msh_dolfinx(geom_msh_path: str, cavern_tag: int = 29):
    """
    Returns:
      X: (Nnodes,3) mesh node coordinates (float)
      facets: (Nfacets,3) vertex indices (int) of triangular facets on the cavern boundary
    using dolfinx.io.gmshio.read_from_msh to get facet tags.
    """
    mesh, cell_tags, facet_tags = dfx.io.gmshio.read_from_msh(geom_msh_path, MPI.COMM_WORLD, 0)

    if facet_tags is None:
        raise RuntimeError("No facet tags were read from the .msh (facet_tags is None).")

    # Find facet indices with physical tag = cavern_tag
    fdim = mesh.topology.dim - 1
    cavern_facets = facet_tags.find(cavern_tag)
    if cavern_facets.size == 0:
        raise RuntimeError(f"No facets found with physical tag {cavern_tag} (Cavern) in {geom_msh_path}")

    # Need facet -> vertex connectivity
    mesh.topology.create_connectivity(fdim, 0)
    f2v = mesh.topology.connectivity(fdim, 0)

    # Build (Nfacets,3) vertex index array (tet mesh boundary facets are triangles)
    facets = []
    for f in cavern_facets:
        verts = f2v.links(int(f))
        if len(verts) != 3:
            # should not happen for tetra boundary, but keep safe
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
    Heuristic: use global cavern surface centroid as 'center', and flip triangles whose
    normal points toward the center (dot(nA, centroid - center) < 0).
    """
    center = np.mean(points[np.unique(tris.reshape(-1))], axis=0)
    tri_centroids = (points[tris[:, 0]] + points[tris[:, 1]] + points[tris[:, 2]]) / 3.0
    v = tri_centroids - center
    s = np.einsum("ij,ij->i", area_vecs, v)  # dot(nA, centroid-center)
    flip = s < 0.0
    area_vecs2 = area_vecs.copy()
    area_vecs2[flip] *= -1.0
    return area_vecs2


def compute_convergence_3d_percent(case_path: str, cavern_phys_tag: int = 29):
    """
    3D cavern volume convergence during the OPERATION stage, referenced to t=0 of the
    displacement series in operation/u/u.xdmf.

    Uses small-strain approximation:
        ΔV(t) ≈ ∫_Γ (u(t) - u(0)) · n dS
    with V0 estimated from the same cavern surface using:
        V0 = (1/3) ∫_Γ x · n dS

    Returns:
      t_days: (Nt,)
      conv_pct: (Nt,)  where conv_pct[0] == 0 (up to numerical roundoff)
    """
    geom_msh = path_geom_msh(case_path)
    u_xdmf = path_u_xdmf(case_path)

    # --- Cavern boundary facets (triangles) from dolfinx gmsh reader ---
    pts_msh, tris_msh = extract_cavern_facets_from_msh_dolfinx(
        geom_msh, cavern_tag=cavern_phys_tag
    )

    # --- Triangle area vectors (n * area) and consistent outward orientation ---
    area_vecs = _area_vectors(pts_msh, tris_msh)
    area_vecs = orient_area_vectors_outward(pts_msh, tris_msh, area_vecs)

    # --- V0 from divergence theorem: V = (1/3) * sum( x_c · (nA) ) ---
    tri_centroids = (pts_msh[tris_msh[:, 0]] +
                     pts_msh[tris_msh[:, 1]] +
                     pts_msh[tris_msh[:, 2]]) / 3.0
    V0 = (1.0 / 3.0) * np.sum(np.einsum("ij,ij->i", tri_centroids, area_vecs))
    V0 = float(abs(V0))

    if not np.isfinite(V0) or V0 <= 0.0:
        raise RuntimeError(f"Computed invalid V0={V0} from cavern surface in {geom_msh}")

    # --- Read nodal displacement time series from XDMF/HDF5 ---
    points_xdmf, time_list, u_field = post.read_node_vector(u_xdmf)  # u_field: (Nt, N, 3)

    if u_field.ndim != 3 or u_field.shape[2] != 3:
        raise RuntimeError(f"Unexpected u_field shape {u_field.shape} from {u_xdmf}")

    Nt = u_field.shape[0]
    if Nt < 1:
        raise RuntimeError(f"No timesteps found in {u_xdmf}")

    # --- Map dolfinx mesh nodes -> xdmf nodes ---
    mapping = post.build_mapping(pts_msh, points_xdmf)

    # Map triangle vertex indices into xdmf node indices
    tri_xdmf = np.vectorize(mapping.__getitem__)(tris_msh).astype(int)
    v0 = tri_xdmf[:, 0]
    v1 = tri_xdmf[:, 1]
    v2 = tri_xdmf[:, 2]

    # --- Reference displacement at t=0 (operation start) ---
    u_ref = u_field[0]  # (N,3)

    # --- ΔV(t) ≈ ∑ (u_c · nA) with u_c = mean(u_i - u_ref_i) over triangle vertices ---
    dV = np.zeros(Nt, dtype=float)
    for k in range(Nt):
        u0 = u_field[k, v0, :] - u_ref[v0, :]
        u1 = u_field[k, v1, :] - u_ref[v1, :]
        u2 = u_field[k, v2, :] - u_ref[v2, :]
        uc = (u0 + u1 + u2) / 3.0
        dV[k] = np.sum(np.einsum("ij,ij->i", uc, area_vecs))

    # Convergence: positive means volume reduction (closure)
    conv_pct = -100.0 * (dV / V0)

    # Time in days
    t_days = np.asarray(time_list, float) / DAY

    # Force exact 0 at start (numerical tiny noise can exist)
    conv_pct[0] = 0.0

    return t_days, conv_pct




def read_case_pressure_scenario(case_path: str) -> str | None:
    """
    Returns pressure scenario string ("sinus","irregular","csv_profile",...)
    from pressure_schedule.json, without confusing cavern folder names.
    """
    pjson = path_pressure_json(case_path)
    if not os.path.isfile(pjson):
        return None

    try:
        with open(pjson, "r") as f:
            data = json.load(f)
    except Exception:
        return None

    # preferred key in your newer scripts
    if isinstance(data, dict) and "pressure_scenario" in data:
        val = data.get("pressure_scenario")
        return str(val).lower() if val is not None else None

    # fallback: some scripts store it under "scenario"
    if isinstance(data, dict) and "scenario" in data:
        val = data.get("scenario")
        # only accept if it looks like a pressure scheme, not "desai_only" etc.
        if str(val).lower() in ("sinus", "irregular", "csv_profile", "linear"):
            return str(val).lower()

    return None


def collect_cases_all(ROOT: str):
    """
    Collect ALL cases (no scheme filtering via folder name).
    Returns list of dicts:
      {label, group, case_name, case_path, pressure_scenario}
    """
    cases = []
    for group in sorted(os.listdir(ROOT)):
        group_path = os.path.join(ROOT, group)
        if not os.path.isdir(group_path):
            continue
        if group.lower().startswith("pressure_"):
            continue

        label = cavern_label_from_group(group)

        for sub in sorted(os.listdir(group_path)):
            if not sub.lower().startswith("case_"):
                continue

            case_path = os.path.join(group_path, sub)
            if not os.path.isdir(case_path):
                continue

            # required for your 3D convergence method
            u_xdmf = path_u_xdmf(case_path)
            u_h5 = path_u_h5(case_path)
            msh = path_geom_msh(case_path)
            if not (os.path.isfile(u_xdmf) and os.path.isfile(u_h5) and os.path.isfile(msh)):
                continue
            if not is_valid_hdf5(u_h5):
                continue

            pres = read_case_pressure_scenario(case_path)

            cases.append({
                "label": label,
                "group": group,
                "case_name": sub,
                "case_path": case_path,
                "pressure_scenario": pres,  # <- from JSON, not from folder name
                "has_pressure_json": os.path.isfile(path_pressure_json(case_path)),
            })
    return cases


def filter_cases(cases, select: dict):
    cavs = select.get("caverns", None)
    pres = select.get("pressure", None)
    contains = select.get("case_contains", None)

    out = []
    for c in cases:
        if cavs is not None and c["group"] not in cavs and c["label"] not in cavs:
            # allow selecting by raw folder name OR pretty label
            continue
        if pres is not None:
            if (c["pressure_scenario"] or "").lower() != pres.lower():
                continue
        if contains is not None and contains.lower() not in c["case_name"].lower():
            continue
        out.append(c)
    return out

def plot_scheme_from_case_list(cases, *, title, color_map, cavern_phys_tag=CAVERN_PHYS_TAG):
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(13, 7), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12},
    )
    fig.suptitle(title, fontsize=12)

    # --- convergence lines ---
    plotted_any = False
    for c in cases:
        lab = c["label"]
        col = color_map.get(lab, None)

        try:
            t_days, conv = compute_convergence_3d_percent(c["case_path"], cavern_phys_tag=cavern_phys_tag)
        except Exception as e:
            print(f"[SKIP] convergence {c['group']}/{c['case_name']}: {e}")
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


def main():
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

    # 2) Apply selection filters (safe: pressure from JSON)
    selected = filter_cases(all_cases, SELECT)
    if not selected:
        raise RuntimeError(f"No cases match SELECT filters: {SELECT}")

    # 3) Build GLOBAL color map across *selected* caverns
    labels_present = []
    for c in selected:
        if c["label"] not in labels_present:
            labels_present.append(c["label"])

    labels_sorted = [l for l in CAVERN_ORDER if l in labels_present] + \
                    [l for l in labels_present if l not in CAVERN_ORDER]
    global_cmap = build_color_map(labels_sorted)

    # 4) Decide plotting grouping
    if SELECT.get("separate_per_pressure", False):
        # group by pressure_scenario string
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
        plt.show()

    else:
        # single plot with all selected cases together
        fig = plot_scheme_from_case_list(
            selected,
            title="3D cavern volume convergence – selected cases",
            color_map=global_cmap,
            cavern_phys_tag=SELECT.get("cavern_phys_tag", CAVERN_PHYS_TAG),
        )
        plt.show()



if __name__ == "__main__":
    main()
