"""
Generate heterogeneous cavern grids with TILTED localized interlayer bands.

Creates two grids:
1. cavern_dipping_interlayer_1200_3D  - regular cavern with 2 localized tilted interlayers
2. cavern_dipping_nointerlayer_1200_3D - identical cavern, single homogeneous salt volume

The interlayers are localized bands (not full-domain planes) that cross the
cavern at specific locations, mimicking the pattern seen in real salt domes
(see Cavern_heterogeneity_example.png):
  - Upper band: enters from upper-right, crosses cavern near the top
  - Lower band: enters from lower-left, crosses cavern near the bottom

Based on the regular 1200k cavern profile.

Usage:
    python generate_heterogeneous_tilted.py
"""
import gmsh
import os
import math

# ── Domain ───────────────────────────────────────────────────────────────────
Lx, Ly, Lz = 450.0, 450.0, 660.0
xc, yc = Lx / 2, Ly / 2  # 225, 225

# ── Cavern geometry (regular 1200k) ─────────────────────────────────────────
h_bottom = 194.951186
R_cav    = 47.996425
H_cav    = 102.270532

z_bot_tip  = h_bottom                          # ~195
z_bot_ring = h_bottom + R_cav                  # ~243
z_top_ring = h_bottom + R_cav + H_cav          # ~345
z_top_tip  = h_bottom + 2 * R_cav + H_cav      # ~393
z_cav_mid  = (z_bot_tip + z_top_tip) / 2       # ~294

# ── Interlayer parameters ────────────────────────────────────────────────────
# Dip angle from horizontal (degrees)
DIP_ANGLE_DEG = 65.0

# Interlayer thickness (meters)
INTERLAYER_THICKNESS = 3.0

# Band dimensions: how far each band extends (length along dip, width in y)
BAND_LENGTH = 250.0   # extent along the tilted direction
BAND_WIDTH  = 450.0   # extent in y-direction (full domain width)

# Upper band: enters from upper-right, crosses cavern near the top
# The band center (where it crosses x=xc) is at z ~ upper part of cavern
UPPER_BAND_Z_AT_CENTER = 360.0   # z where upper band crosses x=xc
UPPER_BAND_X_OFFSET    = 40.0    # shift band center to the right of cavern axis

# Lower band: enters from lower-left, crosses cavern near the bottom
LOWER_BAND_Z_AT_CENTER = 240.0   # z where lower band crosses x=xc
LOWER_BAND_X_OFFSET    = -40.0   # shift band center to the left of cavern axis

# ── Mesh sizes ───────────────────────────────────────────────────────────────
SIZE_COARSE = 65.0
SIZE_FINE   = 4.5

# ── Output directories ──────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DIR_INTERLAYER   = os.path.join(SCRIPT_DIR, "cavern_dipping_interlayer_1200_3D")
DIR_NOINTERLAYER = os.path.join(SCRIPT_DIR, "cavern_dipping_nointerlayer_1200_3D")


def create_localized_tilted_band(z_at_center, x_offset, thickness, dip_deg,
                                  band_length, band_width):
    """
    Create a localized tilted interlayer band.

    The band is a tilted slab with finite extent, centered near the cavern.
    It dips in the x-direction at dip_deg from horizontal.

    Args:
        z_at_center: z-coordinate where band crosses x = xc + x_offset
        x_offset: horizontal offset of band center from cavern axis
        thickness: band thickness (meters)
        dip_deg: dip angle from horizontal (degrees)
        band_length: extent along the tilted dip direction (meters)
        band_width: extent in y-direction (meters)
    """
    occ = gmsh.model.occ
    dip_rad = math.radians(dip_deg)

    # Create horizontal slab centered at the band position
    # The slab is band_length long in x, band_width in y, thickness in z
    cx = xc + x_offset
    slab = occ.addBox(
        cx - band_length / 2, yc - band_width / 2, z_at_center - thickness / 2,
        band_length, band_width, thickness
    )

    # Rotate around y-axis through the band center point
    occ.rotate([(3, slab)], cx, yc, z_at_center, 0, 1, 0, dip_rad)

    # Clip to domain box (the band may extend outside after rotation)
    clip_box = occ.addBox(-0.1, -0.1, -0.1, Lx + 0.2, Ly + 0.2, Lz + 0.2)
    out_map = occ.intersect([(3, slab)], [(3, clip_box)],
                             removeObject=True, removeTool=True)
    occ.synchronize()

    result_vols = [v[1] for v in out_map[0] if v[0] == 3]
    assert len(result_vols) >= 1, f"Intersection produced {len(result_vols)} volumes"
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


def is_cavern_surface(surf_tag, margin=5.0):
    """Check if a surface is part of the cavern wall."""
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


def setup_mesh_refinement(cavern_surfs, all_surfs):
    """
    Set fine mesh at the cavern wall with smooth transition to coarse.

    Uses a Distance field from cavern surfaces combined with a Box field
    so that refinement only applies near the cavern, preventing it from
    leaking to the domain boundary via interlayer surfaces.
    """
    # Set all points to coarse by default
    all_pts = gmsh.model.getEntities(0)
    gmsh.model.mesh.setSize(all_pts, SIZE_COARSE)

    # Set fine size on cavern surface points
    refined_pts = set()
    for s in cavern_surfs:
        boundary = gmsh.model.getBoundary([(2, s)], combined=False, oriented=False)
        for b in boundary:
            pts = gmsh.model.getBoundary([b], combined=False, oriented=False)
            for p in pts:
                if p[1] not in refined_pts:
                    refined_pts.add(p[1])
                    gmsh.model.mesh.setSize([p], SIZE_FINE)

    print(f"  Refined {len(refined_pts)} cavern geometry points to size_fine={SIZE_FINE}")

    if cavern_surfs:
        f_dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", cavern_surfs)

        f_thresh = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(f_thresh, "InField", f_dist)
        gmsh.model.mesh.field.setNumber(f_thresh, "SizeMin", SIZE_FINE)
        gmsh.model.mesh.field.setNumber(f_thresh, "SizeMax", SIZE_COARSE)
        gmsh.model.mesh.field.setNumber(f_thresh, "DistMin", 0)
        gmsh.model.mesh.field.setNumber(f_thresh, "DistMax", 100)

        # Box field to restrict refinement to cavern vicinity
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

        # Max of both: refinement only where near cavern AND inside box
        f_max = gmsh.model.mesh.field.add("Max")
        gmsh.model.mesh.field.setNumbers(f_max, "FieldsList", [f_thresh, f_box])

        gmsh.model.mesh.field.setAsBackgroundMesh(f_max)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


def generate_mesh(with_interlayers, output_dir):
    """Generate a mesh with or without tilted interlayer bands."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("heterogeneous_cavern")

    occ = gmsh.model.occ

    # ── 1. Domain box ────────────────────────────────────────────────────────
    box = occ.addBox(0, 0, 0, Lx, Ly, Lz)

    # ── 2. Cavern (capsule: cylinder + two hemispheres) ──────────────────────
    sphere_bot = occ.addSphere(xc, yc, z_bot_ring, R_cav)
    sphere_top = occ.addSphere(xc, yc, z_top_ring, R_cav)
    cyl = occ.addCylinder(xc, yc, z_bot_ring, 0, 0, H_cav, R_cav)

    cavern_parts = [(3, sphere_bot), (3, sphere_top), (3, cyl)]
    occ.fuse([cavern_parts[0]], cavern_parts[1:],
             removeObject=True, removeTool=True)
    occ.synchronize()

    all_vols = [v[1] for v in gmsh.model.getEntities(3)]
    cavern_vols = [v for v in all_vols if v != box]

    occ.cut([(3, box)], [(3, v) for v in cavern_vols],
            removeObject=True, removeTool=True)
    occ.synchronize()

    if with_interlayers:
        # ── 3. Create localized tilted interlayer bands ──────────────────────
        # Upper band: from upper-right, crossing cavern top
        band_upper = create_localized_tilted_band(
            UPPER_BAND_Z_AT_CENTER, UPPER_BAND_X_OFFSET,
            INTERLAYER_THICKNESS, DIP_ANGLE_DEG,
            BAND_LENGTH, BAND_WIDTH
        )

        # Lower band: from lower-left, crossing cavern bottom
        band_lower = create_localized_tilted_band(
            LOWER_BAND_Z_AT_CENTER, LOWER_BAND_X_OFFSET,
            INTERLAYER_THICKNESS, DIP_ANGLE_DEG,
            BAND_LENGTH, BAND_WIDTH
        )
        occ.synchronize()

        # Fragment: split domain-with-hole by the interlayer bands
        all_vols = [v[1] for v in gmsh.model.getEntities(3)]
        domain_vol = [v for v in all_vols if v not in [band_upper, band_lower]]

        occ.fragment([(3, v) for v in domain_vol],
                     [(3, band_upper), (3, band_lower)],
                     removeObject=True, removeTool=True)
        occ.synchronize()

    # ── 4. Save geometry for visual inspection ─────────────────────────────
    os.makedirs(output_dir, exist_ok=True)
    # Save BREP geometry (can be opened in gmsh directly)
    brep_path = os.path.join(output_dir, "geom.brep")
    gmsh.write(brep_path)
    print(f"  BREP geometry written to: {brep_path}")

    # ── 5. Assign physical groups ────────────────────────────────────────────
    all_vols = [v[1] for v in gmsh.model.getEntities(3)]
    all_surfs = [s[1] for s in gmsh.model.getEntities(2)]

    print(f"  Total volumes: {len(all_vols)}, surfaces: {len(all_surfs)}")

    # Classify surfaces
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

    for name, tags, pg_id in [("Bottom", bottom_s, 27), ("Top", top_s, 22),
                               ("South", south_s, 23), ("North", north_s, 24),
                               ("East", east_s, 25), ("West", west_s, 26),
                               ("Cavern", cavern_s, 29)]:
        if tags:
            gmsh.model.addPhysicalGroup(2, tags, pg_id, name)

    # Classify volumes
    if with_interlayers:
        dip_rad = math.radians(DIP_ANGLE_DEG)

        # Use getMass (actual geometric volume) to separate thin bands from salt
        vol_data = []
        for v in all_vols:
            mass = gmsh.model.occ.getMass(3, v)
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(3, v)
            cx = (xmin + xmax) / 2.0
            cz = (zmin + zmax) / 2.0
            vol_data.append((v, mass, cx, cz))
            print(f"    Vol {v}: mass={mass:.0f} m3")

        # Interlayer band volume: ~BAND_LENGTH * BAND_WIDTH * THICKNESS / cos(dip)
        # ~ 250 * 450 * 3 / cos(65) = ~800,000 m3 (before cavern cut)
        # Salt blocks: tens of millions of m3
        INTERLAYER_VOL_THRESHOLD = 2_000_000

        inter1_vols = []  # upper band
        inter2_vols = []  # lower band
        salt_vols = []

        for v, mass, cx, cz in vol_data:
            if mass < INTERLAYER_VOL_THRESHOLD:
                # Determine which band: upper (higher z centroid) or lower
                # Upper band centered at z~360, lower at z~240
                z_upper = UPPER_BAND_Z_AT_CENTER - (cx - xc - UPPER_BAND_X_OFFSET) * math.tan(dip_rad)
                z_lower = LOWER_BAND_Z_AT_CENTER - (cx - xc - LOWER_BAND_X_OFFSET) * math.tan(dip_rad)
                d_upper = abs(cz - z_upper)
                d_lower = abs(cz - z_lower)

                if d_upper < d_lower:
                    inter1_vols.append(v)
                    print(f"    -> Vol {v}: INTERLAYER_1 (upper, mass={mass:.0f})")
                else:
                    inter2_vols.append(v)
                    print(f"    -> Vol {v}: INTERLAYER_2 (lower, mass={mass:.0f})")
            else:
                salt_vols.append(v)

        print(f"  Classification: {len(salt_vols)} salt, "
              f"{len(inter1_vols)} interlayer_1 (upper), "
              f"{len(inter2_vols)} interlayer_2 (lower)")

        # All salt volumes in one group (since bands are localized,
        # there's no clean bottom/middle/top separation)
        if salt_vols:
            gmsh.model.addPhysicalGroup(3, salt_vols, 31, "Salt_bottom")
            # Also need Salt_middle and Salt_top for compatibility with run_interlayer.py
            # but with localized bands the salt is mostly one connected region
            # We put all salt in Salt_bottom for now (run_interlayer matches "salt" in name)

        if inter1_vols:
            gmsh.model.addPhysicalGroup(3, inter1_vols, 32, "Interlayer_1")
        if inter2_vols:
            gmsh.model.addPhysicalGroup(3, inter2_vols, 34, "Interlayer_2")
    else:
        gmsh.model.addPhysicalGroup(3, all_vols, 28, "Salt")

    print(f"  Surfaces: Bottom={len(bottom_s)}, Top={len(top_s)}, "
          f"South={len(south_s)}, North={len(north_s)}, "
          f"West={len(west_s)}, East={len(east_s)}, Cavern={len(cavern_s)}")

    # ── 6. Mesh refinement near cavern (+ interlayer faces near cavern) ────
    setup_mesh_refinement(cavern_s, all_surfs)

    # ── 7. Generate mesh ─────────────────────────────────────────────────────
    gmsh.model.mesh.generate(3)

    # ── 8. Save ──────────────────────────────────────────────────────────────
    output_path = os.path.join(output_dir, "geom.msh")
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(output_path)

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_tets = sum(len(t) for t in elem_tags)
    print(f"  Mesh: {len(node_tags)} nodes, {n_tets} tetrahedra")
    print(f"  Written to: {output_path}")

    # ── 9. Write .geo with mesh fields for visual inspection ─────────────
    # Collect cavern geometry point tags
    cavern_pt_tags = set()
    for s in cavern_s:
        boundary = gmsh.model.getBoundary([(2, s)], combined=False, oriented=False)
        for b in boundary:
            pts = gmsh.model.getBoundary([b], combined=False, oriented=False)
            for p in pts:
                cavern_pt_tags.add(p[1])

    geo_script_path = os.path.join(output_dir, "geom.geo")
    with open(geo_script_path, "w") as f:
        f.write(f'// Generated by generate_heterogeneous_tilted.py\n')
        f.write(f'SetFactory("OpenCASCADE");\n')
        f.write(f'Merge "geom.brep";\n\n')
        # Default coarse size on all points
        f.write(f'Mesh.CharacteristicLengthMax = {SIZE_COARSE};\n')
        f.write(f'Mesh.MeshSizeExtendFromBoundary = 0;\n')
        f.write(f'Mesh.MeshSizeFromPoints = 1;\n')
        f.write(f'Mesh.MeshSizeFromCurvature = 0;\n\n')
        # Fine size on cavern geometry points
        if cavern_pt_tags:
            pts_str = ", ".join(str(p) for p in sorted(cavern_pt_tags))
            f.write(f'// Fine mesh on cavern geometry points\n')
            f.write(f'MeshSize {{ {pts_str} }} = {SIZE_FINE};\n\n')
        # Distance field from cavern surfaces
        surf_str = ", ".join(str(s) for s in cavern_s)
        f.write(f'// Distance-based refinement from cavern surfaces\n')
        f.write(f'Field[1] = Distance;\n')
        f.write(f'Field[1].SurfacesList = {{ {surf_str} }};\n\n')
        f.write(f'Field[2] = Threshold;\n')
        f.write(f'Field[2].InField = 1;\n')
        f.write(f'Field[2].SizeMin = {SIZE_FINE};\n')
        f.write(f'Field[2].SizeMax = {SIZE_COARSE};\n')
        f.write(f'Field[2].DistMin = 0;\n')
        f.write(f'Field[2].DistMax = 100;\n\n')
        # Box field to restrict refinement to cavern vicinity
        f.write(f'// Box to restrict refinement near cavern only\n')
        f.write(f'Field[3] = Box;\n')
        f.write(f'Field[3].VIn = {SIZE_FINE};\n')
        f.write(f'Field[3].VOut = {SIZE_COARSE};\n')
        f.write(f'Field[3].XMin = {xc - R_cav - 30};\n')
        f.write(f'Field[3].XMax = {xc + R_cav + 30};\n')
        f.write(f'Field[3].YMin = {yc - R_cav - 30};\n')
        f.write(f'Field[3].YMax = {yc + R_cav + 30};\n')
        f.write(f'Field[3].ZMin = {z_bot_tip - 20};\n')
        f.write(f'Field[3].ZMax = {z_top_tip + 20};\n')
        f.write(f'Field[3].Thickness = 50;\n\n')
        # Max of both: refinement only where near cavern AND inside box
        f.write(f'Field[4] = Max;\n')
        f.write(f'Field[4].FieldsList = {{2, 3}};\n')
        f.write(f'Background Field = 4;\n')
    print(f"  .geo script written to: {geo_script_path}")

    gmsh.finalize()


def main():
    print("=" * 70)
    print("HETEROGENEOUS TILTED INTERLAYER GRID GENERATOR")
    print("=" * 70)
    print(f"  Cavern: regular 1200k profile")
    print(f"  z_bot_tip = {z_bot_tip:.2f}, z_top_tip = {z_top_tip:.2f}")
    print(f"  Dip angle: {DIP_ANGLE_DEG} deg")
    print(f"  Band length: {BAND_LENGTH} m, width: {BAND_WIDTH} m")
    print(f"  Upper band: z={UPPER_BAND_Z_AT_CENTER}, x_offset={UPPER_BAND_X_OFFSET}")
    print(f"  Lower band: z={LOWER_BAND_Z_AT_CENTER}, x_offset={LOWER_BAND_X_OFFSET}")
    print(f"  Interlayer thickness: {INTERLAYER_THICKNESS} m")
    print()

    print("-" * 70)
    print("Generating: cavern_dipping_nointerlayer_1200_3D (homogeneous)")
    print("-" * 70)
    generate_mesh(with_interlayers=False, output_dir=DIR_NOINTERLAYER)
    print()

    print("-" * 70)
    print("Generating: cavern_dipping_interlayer_1200_3D (heterogeneous)")
    print("-" * 70)
    generate_mesh(with_interlayers=True, output_dir=DIR_INTERLAYER)
    print()

    print("=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
