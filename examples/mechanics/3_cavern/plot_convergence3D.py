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


def scheme_from_case_folder(case_folder_name: str) -> str:
    low = case_folder_name.lower()
    if low.startswith("case_sinus"):
        return "sinus"
    if low.startswith("case_irregular"):
        return "irregular"
    return ""


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


def collect_cases_nested(ROOT, target_scheme: str):
    """
    Returns list of dicts:
      {label, group, case_name, case_path, has_pressure_json}
    Filters by scheme using case folder name: case_sinus... / case_irregular...
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

            case_scheme = scheme_from_case_folder(sub)
            if case_scheme != target_scheme:
                continue

            case_path = os.path.join(group_path, sub)
            if not os.path.isdir(case_path):
                continue

            # required for 3D convergence
            u_xdmf = path_u_xdmf(case_path)
            u_h5 = path_u_h5(case_path)
            msh = path_geom_msh(case_path)

            if not (os.path.isfile(u_xdmf) and os.path.isfile(u_h5) and os.path.isfile(msh)):
                print(f"[SKIP] {group}/{sub} missing u.xdmf/u.h5/geom.msh")
                continue

            if not is_valid_hdf5(u_h5):
                print(f"[SKIP] {group}/{sub} corrupt u.h5 (not HDF5 signature)")
                continue

            cases.append({
                "label": label,
                "group": group,
                "case_name": sub,
                "case_path": case_path,
                "has_pressure_json": os.path.isfile(path_pressure_json(case_path)),
            })

    return cases


def read_pressure_from_case(case_path: str):
    pjson = path_pressure_json(case_path)
    if not os.path.isfile(pjson):
        return None, None
    with open(pjson, "r") as f:
        data = json.load(f)
    t_days = np.asarray(data["t_values"], float) / DAY
    p_mpa = np.asarray(data["p_values"], float) / MPA
    return t_days, p_mpa


# ------------------------
# 3D convergence from cavern surface: V0 and ΔV(t)
# ------------------------
def _iter_surface_cellblocks_with_phys(mesh: meshio.Mesh):
    """
    Yields tuples:
      (cell_type: str, cells: (Nc, nverts) int array, phys_tags: (Nc,) int array)
    for surface cell blocks that have gmsh:physical.
    """
    # Determine how to access cell data for gmsh physical tags
    # meshio typically stores: mesh.cell_data_dict["gmsh:physical"][cell_type] -> tags
    phys = None
    if hasattr(mesh, "cell_data_dict"):
        phys = mesh.cell_data_dict.get("gmsh:physical", None)
    else:
        # older/newer variations – try best effort
        try:
            phys = mesh.cell_data.get("gmsh:physical", None)
        except Exception:
            phys = None

    if phys is None:
        return  # nothing

    # Surface cell types we accept
    def is_surface_type(ct: str) -> bool:
        return isinstance(ct, str) and (ct.startswith("triangle") or ct.startswith("quad"))

    # Iterate through mesh.cells
    # mesh.cells is list[CellBlock] in recent meshio; sometimes dict in older
    if isinstance(mesh.cells, dict):
        for ct, arr in mesh.cells.items():
            if not is_surface_type(ct):
                continue
            tags = phys.get(ct, None)
            if tags is None:
                continue
            yield ct, np.asarray(arr, dtype=int), np.asarray(tags, dtype=int)
    else:
        for cb in mesh.cells:
            ct = getattr(cb, "type", None)
            if not is_surface_type(ct):
                continue
            tags = phys.get(ct, None)
            if tags is None:
                continue
            yield ct, np.asarray(cb.data, dtype=int), np.asarray(tags, dtype=int)


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


def _quad_to_tris(quad: np.ndarray):
    """
    quad: (N,4) -> two triangles per quad: (2N,3)
    """
    t1 = quad[:, [0, 1, 2]]
    t2 = quad[:, [0, 2, 3]]
    return np.vstack([t1, t2])


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





# ------------------------
# Plotting (multi-cavern per scheme, 2-panel)
# ------------------------
def plot_scheme(ROOT: str, scheme: str, *, color_map, cavern_phys_tag: int = CAVERN_PHYS_TAG):
    cases = collect_cases_nested(ROOT, scheme)
    if not cases:
        raise RuntimeError(f"No '{scheme}' cases found under nested folders in {ROOT}.")

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(13, 7), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.12},
    )
    fig.suptitle(f"3D cavern volume convergence – pressure scheme: {scheme}", fontsize=12)

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

    # --- pressure schedule: take first available in this scheme ---
    ref_t, ref_p = None, None
    for c in cases:
        t, p = read_pressure_from_case(c["case_path"])
        if t is not None:
            ref_t, ref_p = t, p
            break

    if ref_t is None:
        ax2.text(0.5, 0.5, "No pressure_schedule.json found for this scheme.",
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

    
    
    ROOT = "/home/gvandenbrekel/SafeInCave/OutputNobian"
    schemes = ("sinus", "irregular")

    # Build GLOBAL color map across all schemes (consistent colors everywhere)
    all_cases = []
    for s in schemes:
        all_cases.extend(collect_cases_nested(ROOT, s))

    all_labels = []
    for c in all_cases:
        if c["label"] not in all_labels:
            all_labels.append(c["label"])

    labels_sorted = [l for l in CAVERN_ORDER if l in all_labels] + \
                    [l for l in all_labels if l not in CAVERN_ORDER]
    global_cmap = build_color_map(labels_sorted)

    # Plot per scheme
    for s in schemes:
        try:
            plot_scheme(ROOT, s, color_map=global_cmap, cavern_phys_tag=CAVERN_PHYS_TAG)
        except Exception as e:
            print(f"[ERROR] {s}: {e}")

    plt.show()


if __name__ == "__main__":
    main()
