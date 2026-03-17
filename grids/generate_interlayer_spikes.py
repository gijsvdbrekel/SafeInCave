"""
Generate cavern grids with tilted interlayer spike/ledge features.

The interlayers are unleachable, so:
  - The interlayer itself protrudes into the cavern void
  - The salt on the far side (away from leaching source) was shielded
    from brine and remains solid ("shadow zone")
  - This creates pronounced spike/ledge features at the cavern wall

For the UPPER interlayer: shadow zone is ABOVE (salt above wasn't reached)
For the LOWER interlayer: shadow zone is BELOW (salt below wasn't reached)

Creates THREE grids:
  1. cavern_spike_upper_1200_3D  - upper interlayer only (1 spike)
  2. cavern_spike_lower_1200_3D  - lower interlayer only (1 spike)
  3. cavern_spike_none_1200_3D   - same capsule, no interlayers (control)

Physical groups:
  Surfaces: Bottom(27), Top(22), South(23), North(24), East(25), West(26), Cavern(29)
  Volumes:  Salt(28) for no-interlayer; Salt_bottom(31) + Interlayer_1(32) for upper;
            Salt_bottom(31) + Interlayer_2(34) for lower

Usage:
    python generate_interlayer_spikes.py
"""
import gmsh
import os
import math

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

# Shadow zone: unleached salt on the far side of each interlayer.
# Upper interlayer: shadow above; lower interlayer: shadow below.
SHADOW_THICKNESS = 12.0  # meters of unleached salt

# Band extent
BAND_LENGTH = 250.0   # along tilted dip direction (m)
BAND_WIDTH  = 450.0   # y-direction (full domain width)

# Upper interlayer: enters from upper-right, crosses cavern near the top
UPPER_BAND_Z_AT_CENTER = 360.0
UPPER_BAND_X_OFFSET    = 40.0

# Lower interlayer: enters from lower-left, crosses cavern near the bottom
LOWER_BAND_Z_AT_CENTER = 240.0
LOWER_BAND_X_OFFSET    = -40.0

# ── Mesh sizes ─────────────────────────────────────────────────────────────
SIZE_COARSE = 65.0
SIZE_FINE   = 4.5

# ── Output ─────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def create_tilted_slab(z_at_center, x_offset, thickness, dip_deg,
                       band_length, band_width, z_shift=0.0):
    """
    Create a tilted slab (rotated box clipped to domain).

    z_shift: offset of slab center from the interlayer midplane BEFORE rotation.
      - z_shift=0  -> slab centered on interlayer midplane (symmetric)
      - z_shift>0  -> slab shifted upward (toward upper/outer side)
      - z_shift<0  -> slab shifted downward (toward lower/outer side)
    """
    occ = gmsh.model.occ
    dip_rad = math.radians(dip_deg)

    cx = xc + x_offset
    z_slab_center = z_at_center + z_shift
    slab = occ.addBox(
        cx - band_length / 2, yc - band_width / 2,
        z_slab_center - thickness / 2,
        band_length, band_width, thickness
    )

    # Rotate around y-axis through the INTERLAYER midplane center
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


def find_cavern_surfaces(all_surfs, all_vols):
    """Find cavern surfaces: internal surfaces that border the void.

    A cavern surface belongs to exactly one volume (salt/interlayer on one
    side, void on the other). Domain boundary surfaces also belong to one
    volume but are identified by classify_surface(). Interlayer-salt
    interface surfaces belong to two volumes.
    """
    # Build a set of volumes for fast lookup
    vol_set = set(all_vols)

    # For each surface, count how many of the mesh volumes it borders
    cavern_s = []
    for s in all_surfs:
        # Skip domain boundary surfaces
        if classify_surface(s) is not None:
            continue

        # Get the volumes that use this surface
        up = gmsh.model.getAdjacencies(2, s)
        # getAdjacencies(dim, tag) returns (upward_adjacencies, downward_adjacencies)
        parent_vols = [v for v in up[0] if v in vol_set]

        # Surface borders exactly one volume → other side is the void (cavern)
        if len(parent_vols) == 1:
            cavern_s.append(s)

    return cavern_s


def setup_mesh_refinement(cavern_surfs):
    """Set fine mesh at cavern wall with smooth transition to coarse.

    Uses Max of Distance-threshold and Box fields so that refinement
    only applies close to the cavern wall, not the entire box region.
    """
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

        # Max: only refine where BOTH fields agree (near cavern AND inside box)
        f_max = gmsh.model.mesh.field.add("Max")
        gmsh.model.mesh.field.setNumbers(f_max, "FieldsList", [f_thresh, f_box])

        gmsh.model.mesh.field.setAsBackgroundMesh(f_max)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


def write_geo_script(geo_path, label, output_dir):
    """Write .geo script for visual inspection.

    Re-loads the BREP to get correct surface/point tags, since tags from
    the Python generation session don't match the BREP file numbering.
    Uses the same approach as standard cavern .geo files (e.g. regular_1200k.geo):
    fine mesh size on cavern geometry points, coarse everywhere else.
    """
    # Re-load BREP to get correct tags for the .geo file
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("geo_writer")
    gmsh.model.occ.importShapes(os.path.join(output_dir, "geom.brep"))
    gmsh.model.occ.synchronize()

    brep_vols = [v[1] for v in gmsh.model.getEntities(3)]
    brep_surfs = [s[1] for s in gmsh.model.getEntities(2)]
    vol_set = set(brep_vols)

    # Find cavern surfaces in the BREP (single-volume, non-boundary)
    brep_cavern_surfs = []
    for s in brep_surfs:
        bb = gmsh.model.getBoundingBox(2, s)
        xmin, ymin, zmin, xmax, ymax, zmax = bb
        # Skip domain boundaries
        if (abs(zmin) < 1 and abs(zmax) < 1) or \
           (abs(zmin - Lz) < 1 and abs(zmax - Lz) < 1) or \
           (abs(ymin) < 1 and abs(ymax) < 1) or \
           (abs(ymin - Ly) < 1 and abs(ymax - Ly) < 1) or \
           (abs(xmin) < 1 and abs(xmax) < 1) or \
           (abs(xmin - Lx) < 1 and abs(xmax - Lx) < 1):
            continue
        up = gmsh.model.getAdjacencies(2, s)
        parent_vols = [v for v in up[0] if v in vol_set]
        if len(parent_vols) == 1:
            brep_cavern_surfs.append(s)

    # Collect geometry points on cavern surfaces
    brep_cavern_pts = set()
    for s in brep_cavern_surfs:
        boundary = gmsh.model.getBoundary([(2, s)], combined=False, oriented=False)
        for b in boundary:
            pts = gmsh.model.getBoundary([b], combined=False, oriented=False)
            for p in pts:
                brep_cavern_pts.add(p[1])

    gmsh.finalize()

    # Write .geo
    with open(geo_path, "w") as f:
        f.write(f'// Generated by generate_interlayer_spikes.py\n')
        f.write(f'// {label}\n')
        f.write('SetFactory("OpenCASCADE");\n')
        f.write('Merge "geom.brep";\n\n')
        f.write(f'Mesh.CharacteristicLengthMax = {SIZE_COARSE};\n')
        f.write('Mesh.MeshSizeExtendFromBoundary = 0;\n')
        f.write('Mesh.MeshSizeFromPoints = 1;\n')
        f.write('Mesh.MeshSizeFromCurvature = 0;\n\n')
        if brep_cavern_pts:
            pts_str = ", ".join(str(p) for p in sorted(brep_cavern_pts))
            f.write(f'// Fine mesh on cavern geometry points\n')
            f.write(f'MeshSize {{ {pts_str} }} = {SIZE_FINE};\n\n')
        if brep_cavern_surfs:
            surf_str = ", ".join(str(s) for s in brep_cavern_surfs)
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
            f.write('// Max: only refine where near cavern AND inside box\n')
            f.write('Field[4] = Max;\n')
            f.write('Field[4].FieldsList = {2, 3};\n')
            f.write('Background Field = 4;\n')
    print(f"  .geo script written to: {geo_path}")


def classify_and_assign_surfaces(all_surfs, all_vols):
    """Classify surfaces and assign physical groups. Returns cavern surfaces."""
    bottom_s, top_s, south_s, north_s, west_s, east_s = [], [], [], [], [], []

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

    cavern_s = find_cavern_surfaces(all_surfs, all_vols)

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
    return cavern_s


def generate_mesh(mode, output_dir):
    """
    Generate a spike cavern mesh.

    mode: "upper" - upper interlayer spike only
          "lower" - lower interlayer spike only
          "none"  - no interlayers (same regular capsule, homogeneous salt)
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add(f"spike_{mode}")

    occ = gmsh.model.occ

    # ── 1. Domain box and capsule cavern ─────────────────────────────────────
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
    print(f"  Capsule volumes: {capsule_vols}")

    # ── 2. Create barrier(s) and cut from capsule ────────────────────────────
    barrier_total_thickness = INTERLAYER_THICKNESS + SHADOW_THICKNESS
    barriers = []

    if mode == "upper":
        upper_z_shift = +SHADOW_THICKNESS / 2.0
        barrier = create_tilted_slab(
            UPPER_BAND_Z_AT_CENTER, UPPER_BAND_X_OFFSET,
            barrier_total_thickness, DIP_ANGLE_DEG,
            BAND_LENGTH, BAND_WIDTH,
            z_shift=upper_z_shift
        )
        barriers.append(barrier)
        print(f"  Created UPPER barrier (interlayer + shadow above): vol {barrier}")

    elif mode == "lower":
        lower_z_shift = -SHADOW_THICKNESS / 2.0
        barrier = create_tilted_slab(
            LOWER_BAND_Z_AT_CENTER, LOWER_BAND_X_OFFSET,
            barrier_total_thickness, DIP_ANGLE_DEG,
            BAND_LENGTH, BAND_WIDTH,
            z_shift=lower_z_shift
        )
        barriers.append(barrier)
        print(f"  Created LOWER barrier (interlayer + shadow below): vol {barrier}")

    # mode == "none": no barriers, regular capsule

    if barriers:
        cut_result = occ.cut(
            [(3, v) for v in capsule_vols],
            [(3, b) for b in barriers],
            removeObject=True, removeTool=True
        )
        occ.synchronize()
        modified_cavern_vols = [v[1] for v in cut_result[0] if v[0] == 3]
    else:
        modified_cavern_vols = capsule_vols

    print(f"  Cavern volumes to subtract: {modified_cavern_vols}")

    # ── 3. Cut cavern from domain box ────────────────────────────────────────
    occ.cut(
        [(3, box)],
        [(3, v) for v in modified_cavern_vols],
        removeObject=True, removeTool=True
    )
    occ.synchronize()

    # ── 4. Fragment by interlayer band (if applicable) ───────────────────────
    if mode == "upper":
        real_band = create_tilted_slab(
            UPPER_BAND_Z_AT_CENTER, UPPER_BAND_X_OFFSET,
            INTERLAYER_THICKNESS, DIP_ANGLE_DEG,
            BAND_LENGTH, BAND_WIDTH,
            z_shift=0.0
        )
        all_vols = [v[1] for v in gmsh.model.getEntities(3)]
        domain_vols = [v for v in all_vols if v != real_band]
        occ.fragment(
            [(3, v) for v in domain_vols],
            [(3, real_band)],
            removeObject=True, removeTool=True
        )
        occ.synchronize()

    elif mode == "lower":
        real_band = create_tilted_slab(
            LOWER_BAND_Z_AT_CENTER, LOWER_BAND_X_OFFSET,
            INTERLAYER_THICKNESS, DIP_ANGLE_DEG,
            BAND_LENGTH, BAND_WIDTH,
            z_shift=0.0
        )
        all_vols = [v[1] for v in gmsh.model.getEntities(3)]
        domain_vols = [v for v in all_vols if v != real_band]
        occ.fragment(
            [(3, v) for v in domain_vols],
            [(3, real_band)],
            removeObject=True, removeTool=True
        )
        occ.synchronize()

    # ── 5. Save BREP ─────────────────────────────────────────────────────────
    os.makedirs(output_dir, exist_ok=True)
    brep_path = os.path.join(output_dir, "geom.brep")
    gmsh.write(brep_path)
    print(f"  BREP written: {brep_path}")

    # ── 6. Physical groups ───────────────────────────────────────────────────
    all_vols = [v[1] for v in gmsh.model.getEntities(3)]
    all_surfs = [s[1] for s in gmsh.model.getEntities(2)]
    print(f"  Total volumes: {len(all_vols)}, surfaces: {len(all_surfs)}")

    cavern_s = classify_and_assign_surfaces(all_surfs, all_vols)

    # Volume classification
    if mode == "none":
        # All salt, single physical group
        gmsh.model.addPhysicalGroup(3, all_vols, 28, "Salt")
        print(f"  Volumes: all {len(all_vols)} -> Salt")
    else:
        INTERLAYER_VOL_THRESHOLD = 2_000_000
        dip_rad = math.radians(DIP_ANGLE_DEG)

        vol_data = []
        for v in all_vols:
            mass = occ.getMass(3, v)
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(3, v)
            cxv = (xmin + xmax) / 2.0
            czv = (zmin + zmax) / 2.0
            vol_data.append((v, mass, cxv, czv))
            print(f"    Vol {v}: mass={mass:.0f} m3")

        interlayer_vols = []
        salt_vols = []

        for v, mass, cxv, czv in vol_data:
            if mass < INTERLAYER_VOL_THRESHOLD:
                interlayer_vols.append(v)
            else:
                salt_vols.append(v)

        if salt_vols:
            gmsh.model.addPhysicalGroup(3, salt_vols, 31, "Salt_bottom")
        if interlayer_vols:
            # Always use Interlayer_1 for single-interlayer caverns,
            # so the user only needs to set INTERLAYER_1_MATERIAL in run_interlayer.py
            gmsh.model.addPhysicalGroup(3, interlayer_vols, 32, "Interlayer_1")

        print(f"  Classification: {len(salt_vols)} salt, "
              f"{len(interlayer_vols)} interlayer")

    # ── 7. Mesh refinement and generation ────────────────────────────────────
    setup_mesh_refinement(cavern_s)

    print("  Generating 3D mesh...")
    gmsh.model.mesh.generate(3)

    # ── 8. Save mesh ─────────────────────────────────────────────────────────
    msh_path = os.path.join(output_dir, "geom.msh")
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(msh_path)

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_tets = sum(len(t) for t in elem_tags)
    print(f"  Mesh: {len(node_tags)} nodes, {n_tets} tetrahedra")
    print(f"  Written to: {msh_path}")

    gmsh.finalize()

    # ── 9. Write .geo script (after finalize, re-loads BREP for correct tags)
    labels = {
        "upper": "Cavern with upper tilted interlayer spike",
        "lower": "Cavern with lower tilted interlayer spike",
        "none":  "Regular capsule cavern, no interlayers (control)",
    }
    geo_path = os.path.join(output_dir, "geom.geo")
    write_geo_script(geo_path, labels[mode], output_dir)


def main():
    print("=" * 70)
    print("INTERLAYER SPIKES GRID GENERATOR")
    print("=" * 70)

    dip_rad = math.radians(DIP_ANGLE_DEG)
    barrier_total = INTERLAYER_THICKNESS + SHADOW_THICKNESS

    cx_upper = xc + UPPER_BAND_X_OFFSET
    z_upper_at_axis = UPPER_BAND_Z_AT_CENTER + (xc - cx_upper) * math.tan(dip_rad)

    cx_lower = xc + LOWER_BAND_X_OFFSET
    z_lower_at_axis = LOWER_BAND_Z_AT_CENTER + (xc - cx_lower) * math.tan(dip_rad)

    print(f"  Cavern: regular 1200k")
    print(f"    z_bot_tip={z_bot_tip:.1f}, z_top_tip={z_top_tip:.1f}")
    print(f"    R_cav={R_cav:.1f}, H_cav={H_cav:.1f}")
    print(f"  Upper interlayer: z={UPPER_BAND_Z_AT_CENTER} at x={cx_upper:.0f}, "
          f"z={z_upper_at_axis:.1f} at cavern axis")
    print(f"  Lower interlayer: z={LOWER_BAND_Z_AT_CENTER} at x={cx_lower:.0f}, "
          f"z={z_lower_at_axis:.1f} at cavern axis")
    print(f"  Barrier: {INTERLAYER_THICKNESS}m interlayer + "
          f"{SHADOW_THICKNESS}m shadow = {barrier_total}m")
    print()

    configs = [
        ("upper", os.path.join(SCRIPT_DIR, "cavern_spike_upper_1200_3D")),
        ("lower", os.path.join(SCRIPT_DIR, "cavern_spike_lower_1200_3D")),
        ("none",  os.path.join(SCRIPT_DIR, "cavern_spike_none_1200_3D")),
    ]

    for mode, output_dir in configs:
        print("-" * 70)
        print(f"Generating: {os.path.basename(output_dir)} (mode={mode})")
        print("-" * 70)
        generate_mesh(mode, output_dir)
        print()

    print("=" * 70)
    print("DONE - All 3 grids generated")
    print("=" * 70)


if __name__ == "__main__":
    main()
