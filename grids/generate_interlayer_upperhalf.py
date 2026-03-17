"""
Generate interlayer_upperhalf cavern grid.

Creates a realistic salt cavern where a single tilted interlayer band crosses
the UPPER HALF of the cavern. The interlayer is not leachable, so:
  - The interlayer itself protrudes into the cavern
  - The salt ABOVE the interlayer (on the outer/far side from the well) was
    SHIELDED from leaching brine and also remains solid
  - This creates a pronounced spike/ledge at the cavern wall

Geometry approach:
  1. Base cavern: regular 1200k capsule (hemisphere + cylinder + hemisphere)
  2. Create a thick asymmetric "barrier" slab = interlayer (3m) + shadow zone
     (~12m of unleached salt above it on the outer side)
  3. Cut the barrier from the capsule → realistic spike shape
     (the curved capsule wall naturally tapers the spike)
  4. Cut the modified cavern from the domain box
  5. Fragment domain by the real 3m interlayer band for material assignment

Physical groups:
  Surfaces: Bottom(27), Top(22), South(23), North(24), East(25), West(26), Cavern(29)
  Volumes:  Salt_bottom(31), Interlayer_1(32)

Usage:
    python generate_interlayer_upperhalf.py
"""
import gmsh
import os
import math
import sys

# ── Domain ─────────────────────────────────────────────────────────────────
Lx, Ly, Lz = 450.0, 450.0, 660.0
xc, yc = Lx / 2, Ly / 2  # 225, 225

# ── Cavern geometry (regular 1200k profile) ───────────────────────────────
h_bottom = 194.951186
R_cav    = 47.996425
H_cav    = 102.270532

z_bot_tip  = h_bottom                          # ~195.0
z_bot_ring = h_bottom + R_cav                  # ~243.0
z_top_ring = h_bottom + R_cav + H_cav          # ~345.2
z_top_tip  = h_bottom + 2 * R_cav + H_cav      # ~393.2
z_cav_mid  = (z_bot_tip + z_top_tip) / 2       # ~294.1

# ── Interlayer parameters ──────────────────────────────────────────────────
DIP_ANGLE_DEG = 65.0
INTERLAYER_THICKNESS = 3.0   # real geological interlayer thickness (m)

# Shadow zone: unleached salt above the interlayer on the outer/far side.
# The brine approaches from below/center and can't leach past the interlayer,
# so the salt on the upper-outer side remains solid.
# This creates the pronounced spike/ledge seen in real caverns.
SHADOW_THICKNESS = 12.0  # meters of unleached salt above interlayer

# Band extent
BAND_LENGTH = 250.0   # along tilted dip direction (m)
BAND_WIDTH  = 450.0   # y-direction (full domain width)

# Upper-half positioning
BAND_Z_AT_CENTER = 350.0
BAND_X_OFFSET    = 20.0

# ── Mesh sizes ─────────────────────────────────────────────────────────────
SIZE_COARSE = 65.0
SIZE_FINE   = 4.5

# ── Output ─────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "cavern_interlayer_upperhalf_1200_3D")


def create_tilted_slab(z_at_center, x_offset, thickness, dip_deg,
                       band_length, band_width, z_shift=0.0):
    """
    Create a tilted slab (rotated box clipped to domain).

    z_shift: offset of slab center from the interlayer midplane BEFORE rotation.
      - z_shift=0 → slab centered on interlayer midplane (symmetric)
      - z_shift>0 → slab shifted upward (toward upper/outer side)

    The rotation pivot is always the interlayer midplane center,
    so the bottom face of the slab stays aligned with the interlayer.
    """
    occ = gmsh.model.occ
    dip_rad = math.radians(dip_deg)

    cx = xc + x_offset
    # Slab center is shifted upward from the interlayer midplane
    z_slab_center = z_at_center + z_shift
    slab = occ.addBox(
        cx - band_length / 2, yc - band_width / 2,
        z_slab_center - thickness / 2,
        band_length, band_width, thickness
    )

    # Rotate around y-axis through the INTERLAYER midplane center
    # (not the slab center), so the interlayer face stays aligned
    occ.rotate([(3, slab)], cx, yc, z_at_center, 0, 1, 0, dip_rad)

    # Clip to domain
    clip_box = occ.addBox(-1.0, -1.0, -1.0, Lx + 2.0, Ly + 2.0, Lz + 2.0)
    out_map = occ.intersect([(3, slab)], [(3, clip_box)],
                             removeObject=True, removeTool=True)
    occ.synchronize()

    result_vols = [v[1] for v in out_map[0] if v[0] == 3]
    assert len(result_vols) >= 1, f"Slab intersection produced {len(result_vols)} volumes"
    return result_vols[0]


def classify_surface(surf_tag, domain_tol=1.0):
    """Classify a boundary surface by its position on the domain boundary."""
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, surf_tag)

    if abs(zmin) < domain_tol and abs(zmax) < domain_tol:
        return "Bottom"
    if abs(zmin - Lz) < domain_tol and abs(zmax - Lz) < domain_tol:
        return "Top"
    if abs(ymin) < domain_tol and abs(ymax) < domain_tol:
        return "South"
    if abs(ymin - Ly) < domain_tol and abs(ymax - Ly) < domain_tol:
        return "North"
    if abs(xmin) < domain_tol and abs(xmax) < domain_tol:
        return "West"
    if abs(xmin - Lx) < domain_tol and abs(xmax - Lx) < domain_tol:
        return "East"
    return None


def is_cavern_surface(surf_tag, margin=10.0):
    """Check if a surface is part of the cavern wall (including spike surfaces)."""
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, surf_tag)
    sx = xmax - xmin
    sy = ymax - ymin

    cx = (xmin + xmax) / 2.0
    cy = (ymin + ymax) / 2.0
    cz = (zmin + zmax) / 2.0

    dist_xy = math.sqrt((cx - xc)**2 + (cy - yc)**2)

    if dist_xy < R_cav + margin and z_bot_tip - margin < cz < z_top_tip + margin:
        if sx < 2 * R_cav + 2 * margin and sy < 2 * R_cav + 2 * margin:
            return True
    return False


def setup_mesh_refinement(cavern_surfs):
    """Set fine mesh at cavern wall with smooth transition to coarse."""
    all_pts = gmsh.model.getEntities(0)
    gmsh.model.mesh.setSize(all_pts, SIZE_COARSE)

    refined_pts = set()
    for s in cavern_surfs:
        boundary = gmsh.model.getBoundary([(2, s)], combined=False, oriented=False)
        for b in boundary:
            pts = gmsh.model.getBoundary([b], combined=False, oriented=False)
            for p in pts:
                if p[1] not in refined_pts:
                    refined_pts.add(p[1])
                    gmsh.model.mesh.setSize([p], SIZE_FINE)

    print(f"  Refined {len(refined_pts)} cavern points to size_fine={SIZE_FINE}")

    if cavern_surfs:
        f_dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", cavern_surfs)

        f_thresh = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(f_thresh, "InField", f_dist)
        gmsh.model.mesh.field.setNumber(f_thresh, "SizeMin", SIZE_FINE)
        gmsh.model.mesh.field.setNumber(f_thresh, "SizeMax", SIZE_COARSE)
        gmsh.model.mesh.field.setNumber(f_thresh, "DistMin", 0)
        gmsh.model.mesh.field.setNumber(f_thresh, "DistMax", 100)

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

        f_min = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(f_min, "FieldsList", [f_thresh, f_box])

        gmsh.model.mesh.field.setAsBackgroundMesh(f_min)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


def generate_mesh():
    """Generate the interlayer_upperhalf mesh."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("interlayer_upperhalf")

    occ = gmsh.model.occ

    # ── 1. Create base shapes ────────────────────────────────────────────────

    # Domain box
    box = occ.addBox(0, 0, 0, Lx, Ly, Lz)

    # Capsule cavern (hemisphere + cylinder + hemisphere)
    sphere_bot = occ.addSphere(xc, yc, z_bot_ring, R_cav)
    sphere_top = occ.addSphere(xc, yc, z_top_ring, R_cav)
    cyl = occ.addCylinder(xc, yc, z_bot_ring, 0, 0, H_cav, R_cav)

    capsule_parts = [(3, sphere_bot), (3, sphere_top), (3, cyl)]
    fuse_result = occ.fuse([capsule_parts[0]], capsule_parts[1:],
                            removeObject=True, removeTool=True)
    occ.synchronize()

    all_vols = [v[1] for v in gmsh.model.getEntities(3)]
    capsule_vols = [v for v in all_vols if v != box]
    print(f"  Capsule volumes: {capsule_vols}")

    # ── 2. Create thick asymmetric barrier ───────────────────────────────────
    # The barrier = interlayer (3m) + shadow zone (12m above interlayer).
    # Total thickness = 15m, shifted so bottom aligns with interlayer bottom.
    #
    # Before rotation the slab is horizontal:
    #   bottom face = z_at_center - INTERLAYER_THICKNESS/2   (interlayer bottom)
    #   top face    = z_at_center + INTERLAYER_THICKNESS/2 + SHADOW_THICKNESS
    #
    # z_shift positions the slab center relative to the interlayer midplane:
    #   slab center = z_at_center + z_shift
    #   z_shift = SHADOW_THICKNESS / 2   (shift upward by half the shadow)

    barrier_total_thickness = INTERLAYER_THICKNESS + SHADOW_THICKNESS
    barrier_z_shift = SHADOW_THICKNESS / 2.0  # shift slab center upward

    print(f"  Creating barrier slab: {INTERLAYER_THICKNESS}m interlayer + "
          f"{SHADOW_THICKNESS}m shadow = {barrier_total_thickness}m total")

    barrier = create_tilted_slab(
        BAND_Z_AT_CENTER, BAND_X_OFFSET,
        barrier_total_thickness, DIP_ANGLE_DEG,
        BAND_LENGTH, BAND_WIDTH,
        z_shift=barrier_z_shift
    )
    print(f"  Barrier slab: vol {barrier}")

    # ── 3. Cut barrier from capsule → modified cavern with spike ─────────────
    # The barrier region (interlayer + shadow zone) stays as solid domain,
    # creating the pronounced spike/ledge at the cavern wall.
    cut_result = occ.cut(
        [(3, v) for v in capsule_vols],
        [(3, barrier)],
        removeObject=True, removeTool=True
    )
    occ.synchronize()

    modified_cavern_vols = [v[1] for v in cut_result[0] if v[0] == 3]
    print(f"  Modified cavern (capsule - barrier): {modified_cavern_vols}")

    # ── 4. Cut modified cavern from box ──────────────────────────────────────
    occ.cut(
        [(3, box)],
        [(3, v) for v in modified_cavern_vols],
        removeObject=True, removeTool=True
    )
    occ.synchronize()

    print(f"  Domain volumes after cavern cut: "
          f"{[v[1] for v in gmsh.model.getEntities(3)]}")

    # ── 5. Create the real 3m interlayer and fragment domain ─────────────────
    # The domain fragmentation uses only the actual geological interlayer (3m),
    # NOT the full barrier. The shadow zone is salt (just unleached salt).
    real_interlayer = create_tilted_slab(
        BAND_Z_AT_CENTER, BAND_X_OFFSET,
        INTERLAYER_THICKNESS, DIP_ANGLE_DEG,
        BAND_LENGTH, BAND_WIDTH,
        z_shift=0.0  # centered on midplane (symmetric 3m)
    )
    print(f"  Real interlayer band (for fragmentation): vol {real_interlayer}")

    all_vols = [v[1] for v in gmsh.model.getEntities(3)]
    domain_vols = [v for v in all_vols if v != real_interlayer]

    occ.fragment(
        [(3, v) for v in domain_vols],
        [(3, real_interlayer)],
        removeObject=True, removeTool=True
    )
    occ.synchronize()

    # ── 6. Save BREP geometry ────────────────────────────────────────────────
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    brep_path = os.path.join(OUTPUT_DIR, "geom.brep")
    gmsh.write(brep_path)
    print(f"  BREP written: {brep_path}")

    # ── 7. Classify and assign physical groups ───────────────────────────────
    all_vols = [v[1] for v in gmsh.model.getEntities(3)]
    all_surfs = [s[1] for s in gmsh.model.getEntities(2)]
    print(f"  Total volumes: {len(all_vols)}, surfaces: {len(all_surfs)}")

    # --- Surface classification ---
    bottom_s, top_s, south_s, north_s, west_s, east_s = [], [], [], [], [], []
    cavern_s = []

    for s in all_surfs:
        cat = classify_surface(s)
        if cat == "Bottom":
            bottom_s.append(s)
        elif cat == "Top":
            top_s.append(s)
        elif cat == "South":
            south_s.append(s)
        elif cat == "North":
            north_s.append(s)
        elif cat == "West":
            west_s.append(s)
        elif cat == "East":
            east_s.append(s)
        elif is_cavern_surface(s):
            cavern_s.append(s)

    for name, tags, pg_id in [
        ("Bottom", bottom_s, 27), ("Top", top_s, 22),
        ("South", south_s, 23), ("North", north_s, 24),
        ("East", east_s, 25), ("West", west_s, 26),
        ("Cavern", cavern_s, 29),
    ]:
        if tags:
            gmsh.model.addPhysicalGroup(2, tags, pg_id, name)

    print(f"  Surfaces: Bottom={len(bottom_s)}, Top={len(top_s)}, "
          f"South={len(south_s)}, North={len(north_s)}, "
          f"West={len(west_s)}, East={len(east_s)}, Cavern={len(cavern_s)}")

    # --- Volume classification ---
    INTERLAYER_VOL_THRESHOLD = 2_000_000  # m3

    vol_data = []
    for v in all_vols:
        mass = occ.getMass(3, v)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(3, v)
        cxv = (xmin + xmax) / 2.0
        czv = (zmin + zmax) / 2.0
        vol_data.append((v, mass, cxv, czv))
        print(f"    Vol {v}: mass={mass:.0f} m3, center=({cxv:.0f}, {czv:.0f})")

    interlayer_vols = []
    salt_vols = []

    for v, mass, cxv, czv in vol_data:
        if mass < INTERLAYER_VOL_THRESHOLD:
            interlayer_vols.append(v)
            print(f"    -> Vol {v}: INTERLAYER_1 (mass={mass:.0f})")
        else:
            salt_vols.append(v)
            print(f"    -> Vol {v}: SALT (mass={mass:.0f})")

    print(f"  Classification: {len(salt_vols)} salt, {len(interlayer_vols)} interlayer")

    if salt_vols:
        gmsh.model.addPhysicalGroup(3, salt_vols, 31, "Salt_bottom")
    if interlayer_vols:
        gmsh.model.addPhysicalGroup(3, interlayer_vols, 32, "Interlayer_1")

    # ── 8. Mesh refinement ───────────────────────────────────────────────────
    setup_mesh_refinement(cavern_s)

    # ── 9. Generate 3D mesh ──────────────────────────────────────────────────
    print("  Generating 3D mesh...")
    gmsh.model.mesh.generate(3)

    # ── 10. Save mesh ────────────────────────────────────────────────────────
    msh_path = os.path.join(OUTPUT_DIR, "geom.msh")
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(msh_path)

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_tets = sum(len(t) for t in elem_tags)
    print(f"  Mesh: {len(node_tags)} nodes, {n_tets} tetrahedra")
    print(f"  Written to: {msh_path}")

    # ── 11. Write .geo script for inspection ─────────────────────────────────
    cavern_pt_tags = set()
    for s in cavern_s:
        boundary = gmsh.model.getBoundary([(2, s)], combined=False, oriented=False)
        for b in boundary:
            pts = gmsh.model.getBoundary([b], combined=False, oriented=False)
            for p in pts:
                cavern_pt_tags.add(p[1])

    geo_path = os.path.join(OUTPUT_DIR, "geom.geo")
    with open(geo_path, "w") as f:
        f.write('// Generated by generate_interlayer_upperhalf.py\n')
        f.write('// Cavern with tilted interlayer spike in upper half\n')
        f.write('// Spike = interlayer (3m) + shadow zone (12m unleached salt above)\n')
        f.write('SetFactory("OpenCASCADE");\n')
        f.write('Merge "geom.brep";\n\n')
        f.write(f'Mesh.CharacteristicLengthMax = {SIZE_COARSE};\n')
        f.write('Mesh.MeshSizeExtendFromBoundary = 0;\n')
        f.write('Mesh.MeshSizeFromPoints = 1;\n')
        f.write('Mesh.MeshSizeFromCurvature = 0;\n\n')
        if cavern_pt_tags:
            pts_str = ", ".join(str(p) for p in sorted(cavern_pt_tags))
            f.write(f'// Fine mesh on cavern geometry points\n')
            f.write(f'MeshSize {{ {pts_str} }} = {SIZE_FINE};\n\n')
        surf_str = ", ".join(str(s) for s in cavern_s)
        f.write('// Distance-based refinement from cavern surfaces\n')
        f.write('Field[1] = Distance;\n')
        f.write(f'Field[1].SurfacesList = {{ {surf_str} }};\n\n')
        f.write('Field[2] = Threshold;\n')
        f.write('Field[2].InField = 1;\n')
        f.write(f'Field[2].SizeMin = {SIZE_FINE};\n')
        f.write(f'Field[2].SizeMax = {SIZE_COARSE};\n')
        f.write('Field[2].DistMin = 0;\n')
        f.write('Field[2].DistMax = 100;\n\n')
        f.write('// Box to restrict refinement near cavern only\n')
        f.write('Field[3] = Box;\n')
        f.write(f'Field[3].VIn = {SIZE_FINE};\n')
        f.write(f'Field[3].VOut = {SIZE_COARSE};\n')
        f.write(f'Field[3].XMin = {xc - R_cav - 30};\n')
        f.write(f'Field[3].XMax = {xc + R_cav + 30};\n')
        f.write(f'Field[3].YMin = {yc - R_cav - 30};\n')
        f.write(f'Field[3].YMax = {yc + R_cav + 30};\n')
        f.write(f'Field[3].ZMin = {z_bot_tip - 20};\n')
        f.write(f'Field[3].ZMax = {z_top_tip + 20};\n')
        f.write('Field[3].Thickness = 50;\n\n')
        f.write('Field[4] = Min;\n')
        f.write('Field[4].FieldsList = {2, 3};\n')
        f.write('Background Field = 4;\n')
    print(f"  .geo script written to: {geo_path}")

    gmsh.finalize()


def main():
    print("=" * 70)
    print("INTERLAYER UPPERHALF GRID GENERATOR")
    print("=" * 70)

    dip_rad = math.radians(DIP_ANGLE_DEG)
    cx_band = xc + BAND_X_OFFSET

    z_at_axis  = BAND_Z_AT_CENTER + (xc - cx_band) * math.tan(dip_rad)
    z_at_right = BAND_Z_AT_CENTER + (xc + R_cav - cx_band) * math.tan(dip_rad)
    z_at_left  = BAND_Z_AT_CENTER + (xc - R_cav - cx_band) * math.tan(dip_rad)

    print(f"  Cavern: regular 1200k")
    print(f"    z_bot_tip={z_bot_tip:.1f}, z_top_tip={z_top_tip:.1f}, "
          f"z_mid={z_cav_mid:.1f}")
    print(f"    R_cav={R_cav:.1f}, H_cav={H_cav:.1f}")
    print(f"  Interlayer band:")
    print(f"    Dip: {DIP_ANGLE_DEG} deg, thickness: {INTERLAYER_THICKNESS} m")
    print(f"    Shadow zone: {SHADOW_THICKNESS} m (unleached salt above interlayer)")
    print(f"    Total barrier: {INTERLAYER_THICKNESS + SHADOW_THICKNESS} m")
    print(f"    Center: x={cx_band:.0f}, z={BAND_Z_AT_CENTER:.0f}")
    print(f"    z at cavern axis  (x={xc:.0f}): {z_at_axis:.1f}")
    print(f"    z at right wall   (x={xc+R_cav:.0f}): {z_at_right:.1f}")
    print(f"    z at left wall    (x={xc-R_cav:.0f}): {z_at_left:.1f}")
    print(f"  Output: {OUTPUT_DIR}")
    print()

    generate_mesh()

    print()
    print("=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
