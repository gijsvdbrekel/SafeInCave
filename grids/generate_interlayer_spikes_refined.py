"""
Generate refined-interlayer spike cavern grids.

Same geometry as generate_interlayer_spikes.py (regular 1200k capsule + tilted
unleachable interlayer + shadow zone), but with extra mesh refinement applied
to the interlayer volume and its surfaces. Two refinement levels are produced
for both the upper and lower spike variants — four grid folders total:

  cavern_spike_upper_il2x_1200_3D    interlayer mesh size = SIZE_FINE / 2
  cavern_spike_lower_il2x_1200_3D    interlayer mesh size = SIZE_FINE / 2
  cavern_spike_upper_il4x_1200_3D    interlayer mesh size = SIZE_FINE / 4
  cavern_spike_lower_il4x_1200_3D    interlayer mesh size = SIZE_FINE / 4

Physical groups, surface tags, BREP/MSH layout and naming convention are
identical to the baseline spike grids, so run_interlayer.py, plot_results.py
and paraview.py work unchanged for the new grids once the new cavern types
are registered (see run_interlayer.py CAVERN_PARAMS / GRID_FOLDERS and
case_index.py _SPIKE_RE).
"""
import gmsh
import math
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from generate_interlayer_spikes import (
    Lx, Ly, Lz, xc, yc,
    R_cav, H_cav,
    z_bot_tip, z_bot_ring, z_top_ring, z_top_tip,
    DIP_ANGLE_DEG, INTERLAYER_THICKNESS, SHADOW_THICKNESS,
    BAND_LENGTH, BAND_WIDTH,
    UPPER_BAND_Z_AT_CENTER, UPPER_BAND_X_OFFSET,
    LOWER_BAND_Z_AT_CENTER, LOWER_BAND_X_OFFSET,
    SIZE_COARSE, SIZE_FINE,
    create_tilted_slab,
    classify_and_assign_surfaces,
    write_geo_script,
)


REFINEMENT_FACTORS = [2.0, 4.0]
MODES = ["upper", "lower"]
INTERLAYER_VOL_THRESHOLD = 2_000_000  # m^3 — anything smaller is the interlayer


def setup_mesh_refinement_with_interlayer(cavern_surfs, interlayer_vols, il_size):
    """Cavern-wall refinement (as in baseline) PLUS an interlayer-targeted
    Threshold field that drives mesh size to `il_size` inside and just around
    the interlayer slab. Combines all fields with Min so the smallest target
    wins per location.
    """
    all_pts = gmsh.model.getEntities(0)
    gmsh.model.mesh.setSize(all_pts, SIZE_COARSE)

    # Cavern wall point sizes
    refined_pts = set()
    for s in cavern_surfs:
        boundary = gmsh.model.getBoundary([(2, s)], combined=False, oriented=False)
        for b in boundary:
            pts = gmsh.model.getBoundary([b], combined=False, oriented=False)
            for p in pts:
                if p[1] not in refined_pts:
                    refined_pts.add(p[1])
                    gmsh.model.mesh.setSize([p], SIZE_FINE)

    # Interlayer geometry point sizes (set to il_size)
    il_pts = set()
    il_surfs = set()
    for v in interlayer_vols:
        _, surfs = gmsh.model.getAdjacencies(3, v)
        for s in surfs:
            il_surfs.add(s)
            curves = gmsh.model.getBoundary([(2, s)], combined=False, oriented=False)
            for c in curves:
                pts = gmsh.model.getBoundary([c], combined=False, oriented=False)
                for p in pts:
                    if p[1] not in il_pts:
                        il_pts.add(p[1])
                        gmsh.model.mesh.setSize([p], il_size)
    il_surfs = list(il_surfs)
    print(f"  Refined {len(refined_pts)} cavern points to size_fine={SIZE_FINE}")
    print(f"  Refined {len(il_pts)} interlayer points to il_size={il_size}")

    # Cavern wall distance + threshold (same as baseline)
    f_dist_cav = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist_cav, "SurfacesList", cavern_surfs)
    f_thresh_cav = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thresh_cav, "InField", f_dist_cav)
    gmsh.model.mesh.field.setNumber(f_thresh_cav, "SizeMin", SIZE_FINE)
    gmsh.model.mesh.field.setNumber(f_thresh_cav, "SizeMax", SIZE_COARSE)
    gmsh.model.mesh.field.setNumber(f_thresh_cav, "DistMin", 0)
    gmsh.model.mesh.field.setNumber(f_thresh_cav, "DistMax", 100)

    # Cavern bounding box (restrict cavern refinement to the cavern region)
    f_box = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(f_box, "VIn", SIZE_FINE)
    gmsh.model.mesh.field.setNumber(f_box, "VOut", SIZE_COARSE)
    gmsh.model.mesh.field.setNumber(f_box, "XMin", xc - R_cav - 30)
    gmsh.model.mesh.field.setNumber(f_box, "XMax", xc + R_cav + 30)
    gmsh.model.mesh.field.setNumber(f_box, "YMin", yc - R_cav - 30)
    gmsh.model.mesh.field.setNumber(f_box, "YMax", yc + R_cav + 30)
    gmsh.model.mesh.field.setNumber(f_box, "ZMin", z_bot_tip - 20)
    gmsh.model.mesh.field.setNumber(f_box, "ZMax", z_top_tip + 20)
    gmsh.model.mesh.field.setNumber(f_box, "Thickness", 50)

    f_max_cav = gmsh.model.mesh.field.add("Max")
    gmsh.model.mesh.field.setNumbers(f_max_cav, "FieldsList",
                                     [f_thresh_cav, f_box])

    # Interlayer distance + threshold
    if il_surfs:
        f_dist_il = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(f_dist_il, "SurfacesList", il_surfs)
        f_thresh_il = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(f_thresh_il, "InField", f_dist_il)
        gmsh.model.mesh.field.setNumber(f_thresh_il, "SizeMin", il_size)
        # SizeMax must be coarse: gmsh Threshold returns SizeMax for *any*
        # distance above DistMax, not a decay. Setting it to SIZE_FINE would
        # clamp the entire bulk salt at 4.5 m via the Min combine.
        gmsh.model.mesh.field.setNumber(f_thresh_il, "SizeMax", SIZE_COARSE)
        gmsh.model.mesh.field.setNumber(f_thresh_il, "DistMin", 0)
        gmsh.model.mesh.field.setNumber(f_thresh_il, "DistMax", 8.0)

        f_min = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(f_min, "FieldsList",
                                         [f_max_cav, f_thresh_il])
        gmsh.model.mesh.field.setAsBackgroundMesh(f_min)
    else:
        gmsh.model.mesh.field.setAsBackgroundMesh(f_max_cav)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


def generate_refined_mesh(mode, factor, output_dir):
    """Generate one refined-interlayer spike cavern mesh."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    factor_int = int(round(factor))
    gmsh.model.add(f"spike_{mode}_il{factor_int}x")

    occ = gmsh.model.occ

    # Domain box + capsule cavern
    box = occ.addBox(0, 0, 0, Lx, Ly, Lz)
    sphere_bot = occ.addSphere(xc, yc, z_bot_ring, R_cav)
    sphere_top = occ.addSphere(xc, yc, z_top_ring, R_cav)
    cyl = occ.addCylinder(xc, yc, z_bot_ring, 0, 0, H_cav, R_cav)
    capsule_parts = [(3, sphere_bot), (3, sphere_top), (3, cyl)]
    occ.fuse([capsule_parts[0]], capsule_parts[1:],
             removeObject=True, removeTool=True)
    occ.synchronize()

    all_vols = [v[1] for v in gmsh.model.getEntities(3)]
    capsule_vols = [v for v in all_vols if v != box]

    # Build the unleachable barrier (interlayer + shadow zone)
    if mode == "upper":
        z_shift = +SHADOW_THICKNESS / 2.0
        z_at_center = UPPER_BAND_Z_AT_CENTER
        x_offset = UPPER_BAND_X_OFFSET
    elif mode == "lower":
        z_shift = -SHADOW_THICKNESS / 2.0
        z_at_center = LOWER_BAND_Z_AT_CENTER
        x_offset = LOWER_BAND_X_OFFSET
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    barrier_total_thickness = INTERLAYER_THICKNESS + SHADOW_THICKNESS
    barrier = create_tilted_slab(
        z_at_center, x_offset, barrier_total_thickness, DIP_ANGLE_DEG,
        BAND_LENGTH, BAND_WIDTH, z_shift=z_shift
    )

    edges_before = {e[1] for e in gmsh.model.getEntities(1)}
    cut_result = occ.cut(
        [(3, v) for v in capsule_vols], [(3, barrier)],
        removeObject=True, removeTool=True
    )
    occ.synchronize()
    modified_cavern_vols = [v[1] for v in cut_result[0] if v[0] == 3]

    # Fillet the new sharp slab/cavern intersection edges
    FILLET_R = 12.0
    new_cavern_edges = []
    for vol in modified_cavern_vols:
        _, surfs = gmsh.model.getAdjacencies(3, vol)
        for surf in surfs:
            _, e_curves = gmsh.model.getAdjacencies(2, surf)
            for ec in e_curves:
                if ec not in edges_before:
                    new_cavern_edges.append(ec)
    new_cavern_edges = sorted(set(new_cavern_edges))
    if not new_cavern_edges:
        all_cavern_edges = []
        for vol in modified_cavern_vols:
            _, surfs = gmsh.model.getAdjacencies(3, vol)
            for surf in surfs:
                _, e_curves = gmsh.model.getAdjacencies(2, surf)
                all_cavern_edges.extend(e_curves)
        new_cavern_edges = sorted(set(all_cavern_edges))
    if new_cavern_edges:
        try:
            fillet_result = occ.fillet(
                modified_cavern_vols, new_cavern_edges, [FILLET_R],
                removeVolume=True
            )
            occ.synchronize()
            modified_cavern_vols = [v[1] for v in fillet_result if v[0] == 3]
            print(f"  Filleted {len(new_cavern_edges)} knife-edges with R={FILLET_R} m")
        except Exception as exc:
            print(f"  WARNING: fillet failed ({exc})")

    # Cut cavern from domain box
    occ.cut(
        [(3, box)], [(3, v) for v in modified_cavern_vols],
        removeObject=True, removeTool=True
    )
    occ.synchronize()

    # Fragment by the real interlayer band so it becomes a separate volume
    real_band = create_tilted_slab(
        z_at_center, x_offset, INTERLAYER_THICKNESS, DIP_ANGLE_DEG,
        BAND_LENGTH, BAND_WIDTH, z_shift=0.0
    )
    all_vols = [v[1] for v in gmsh.model.getEntities(3)]
    domain_vols = [v for v in all_vols if v != real_band]
    occ.fragment(
        [(3, v) for v in domain_vols], [(3, real_band)],
        removeObject=True, removeTool=True
    )
    occ.synchronize()

    os.makedirs(output_dir, exist_ok=True)
    brep_path = os.path.join(output_dir, "geom.brep")
    gmsh.write(brep_path)
    print(f"  BREP written: {brep_path}")

    # Physical groups (matches baseline tags so downstream code is unchanged)
    all_vols = [v[1] for v in gmsh.model.getEntities(3)]
    all_surfs = [s[1] for s in gmsh.model.getEntities(2)]
    cavern_s = classify_and_assign_surfaces(all_surfs, all_vols)

    interlayer_vols = []
    salt_vols = []
    for v in all_vols:
        mass = occ.getMass(3, v)
        if mass < INTERLAYER_VOL_THRESHOLD:
            interlayer_vols.append(v)
        else:
            salt_vols.append(v)
        print(f"    Vol {v}: mass={mass:.0f} m3")
    if salt_vols:
        gmsh.model.addPhysicalGroup(3, salt_vols, 31, "Salt_bottom")
    if interlayer_vols:
        gmsh.model.addPhysicalGroup(3, interlayer_vols, 32, "Interlayer_1")
    print(f"  Classification: {len(salt_vols)} salt, {len(interlayer_vols)} interlayer")

    # Refinement: cavern wall + interlayer
    il_size = SIZE_FINE / factor
    setup_mesh_refinement_with_interlayer(cavern_s, interlayer_vols, il_size)

    print(f"  Generating 3D mesh (il_size = {il_size:g} m)...")
    gmsh.model.mesh.generate(3)

    msh_path = os.path.join(output_dir, "geom.msh")
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(msh_path)

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    _, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_tets = sum(len(t) for t in elem_tags)
    print(f"  Mesh: {len(node_tags)} nodes, {n_tets} tetrahedra")
    print(f"  Written to: {msh_path}")

    gmsh.finalize()

    label = (f"Cavern with {mode} interlayer spike — "
             f"interlayer mesh refined {factor_int}x (il_size={il_size:g} m)")
    geo_path = os.path.join(output_dir, "geom.geo")
    write_geo_script(geo_path, label, output_dir)


def main():
    print("=" * 70)
    print("REFINED-INTERLAYER SPIKE GRID GENERATOR")
    print("=" * 70)
    for factor in REFINEMENT_FACTORS:
        for mode in MODES:
            tag = f"il{int(round(factor))}x"
            output_dir = os.path.join(SCRIPT_DIR,
                                      f"cavern_spike_{mode}_{tag}_1200_3D")
            print("-" * 70)
            print(f"Generating: {os.path.basename(output_dir)}")
            print(f"  Refinement factor: {factor}x  →  "
                  f"il_size = {SIZE_FINE / factor:g} m")
            print("-" * 70)
            generate_refined_mesh(mode, factor, output_dir)
            print()
    print("=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
