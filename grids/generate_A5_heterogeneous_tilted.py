"""
Generate A5 cavern grid with TILTED localized interlayer bands.

Overwrites cavern_A5_interlayer_3D/ with a new mesh that has localized tilted
interlayer bands (matching the style of Cavern_heterogeneity_example.png)
instead of the old horizontal interlayers.

A5 cavern profile (from sonar survey October 2023):
  z=145(tip) -> 175(R=35) -> 220(R=45) -> 270(R=30) -> 320(R=38) ->
  370(R=28) -> 420(R=32) -> 460(R=20) -> 490(R=12) -> 515(tip)
  Volume: ~965,531 m³

Usage:
    python generate_A5_heterogeneous_tilted.py
"""
import gmsh
import os
import math

# ── Domain ───────────────────────────────────────────────────────────────────
Lx, Ly, Lz = 450.0, 450.0, 660.0
xc, yc = Lx / 2, Ly / 2  # 225, 225

# ── A5 cavern profile (z, R) pairs ──────────────────────────────────────────
A5_PROFILE = [
    (145.0,  0.0),   # bottom tip
    (175.0, 35.0),   # lower expansion
    (220.0, 45.0),   # maximum lower bulge
    (270.0, 30.0),   # first constriction
    (320.0, 38.0),   # middle bulge
    (370.0, 28.0),   # upper constriction
    (420.0, 32.0),   # upper body
    (460.0, 20.0),   # neck transition
    (490.0, 12.0),   # chimney
    (515.0,  0.0),   # top tip
]

z_bot_tip = A5_PROFILE[0][0]   # 145
z_top_tip = A5_PROFILE[-1][0]  # 515
R_max = max(r for _, r in A5_PROFILE)  # 45

# ── Interlayer parameters ────────────────────────────────────────────────────
DIP_ANGLE_DEG = 65.0
INTERLAYER_THICKNESS = 3.0
BAND_LENGTH = 250.0
BAND_WIDTH  = 450.0

# Upper band: enters from upper-right, crosses cavern in upper region
# (near the original z=370 constriction)
UPPER_BAND_Z_AT_CENTER = 400.0
UPPER_BAND_X_OFFSET    = 40.0

# Lower band: enters from lower-left, crosses cavern in lower region
# (near the original z=270 constriction)
LOWER_BAND_Z_AT_CENTER = 230.0
LOWER_BAND_X_OFFSET    = -40.0

# ── Mesh sizes ───────────────────────────────────────────────────────────────
SIZE_COARSE = 65.0
SIZE_FINE   = 4.5

# ── Output ───────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "cavern_A5_interlayer_3D")


def build_a5_cavern_occ():
    """
    Build A5 cavern shape in OCC by revolving the profile.

    Creates a 2D wire from the (z, R) profile in the xz-plane centered at
    (xc, yc), then revolves it 360° around the vertical axis.
    Returns the cavern volume tag.
    """
    occ = gmsh.model.occ

    # Build 2D profile as a polyline in the xz-plane at y=yc
    # Points go: bottom_tip -> right side up -> top_tip -> back down axis
    profile_pts = []

    # Bottom tip on axis
    profile_pts.append(occ.addPoint(xc, yc, A5_PROFILE[0][0]))

    # Right side going up (positive x offset = radius)
    for z, R in A5_PROFILE[1:-1]:
        profile_pts.append(occ.addPoint(xc + R, yc, z))

    # Top tip on axis
    profile_pts.append(occ.addPoint(xc, yc, A5_PROFILE[-1][0]))

    # Create lines along the profile (right side)
    right_lines = []
    for i in range(len(profile_pts) - 1):
        right_lines.append(occ.addLine(profile_pts[i], profile_pts[i + 1]))

    # Close with a line along the axis (top tip back to bottom tip)
    axis_line = occ.addLine(profile_pts[-1], profile_pts[0])

    # Create wire and surface
    wire = occ.addCurveLoop(right_lines + [axis_line])
    surf = occ.addPlaneSurface([wire])

    # Revolve 360° around vertical axis through (xc, yc)
    # Axis: point (xc, yc, 0), direction (0, 0, 1)
    revolved = occ.revolve([(2, surf)], xc, yc, 0, 0, 0, 1, 2 * math.pi)
    occ.synchronize()

    # Find the 3D volume from the revolution
    vol_tags = [e[1] for e in revolved if e[0] == 3]
    assert len(vol_tags) >= 1, f"Revolution produced {len(vol_tags)} volumes"
    return vol_tags[0]


def create_localized_tilted_band(z_at_center, x_offset, thickness, dip_deg,
                                  band_length, band_width):
    """Create a localized tilted interlayer band, clipped to domain."""
    occ = gmsh.model.occ
    dip_rad = math.radians(dip_deg)

    cx = xc + x_offset
    slab = occ.addBox(
        cx - band_length / 2, yc - band_width / 2, z_at_center - thickness / 2,
        band_length, band_width, thickness
    )
    occ.rotate([(3, slab)], cx, yc, z_at_center, 0, 1, 0, dip_rad)

    clip_box = occ.addBox(-0.1, -0.1, -0.1, Lx + 0.2, Ly + 0.2, Lz + 0.2)
    out_map = occ.intersect([(3, slab)], [(3, clip_box)],
                             removeObject=True, removeTool=True)
    occ.synchronize()

    result_vols = [v[1] for v in out_map[0] if v[0] == 3]
    assert len(result_vols) >= 1
    return result_vols[0]


def classify_surface(surf_tag, domain_tol=1.0):
    """Classify a boundary surface by position on domain boundary."""
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
    if dist_xy < R_max + margin and z_bot_tip - margin < cz < z_top_tip + margin:
        if sx < 2 * R_max + 2 * margin and sy < 2 * R_max + 2 * margin:
            return True
    return False


def setup_mesh_refinement(cavern_surfs, all_surfs):
    """
    Set fine mesh at the cavern wall with smooth transition to coarse.

    Uses a Distance field from cavern surfaces combined with a Restrict field
    so that refinement only applies inside a box around the cavern, preventing
    it from leaking to the domain boundary via interlayer surfaces.
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

    # Distance field from cavern surfaces for mesh grading on the curved walls
    if cavern_surfs:
        f_dist = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", cavern_surfs)

        f_thresh = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(f_thresh, "InField", f_dist)
        gmsh.model.mesh.field.setNumber(f_thresh, "SizeMin", SIZE_FINE)
        gmsh.model.mesh.field.setNumber(f_thresh, "SizeMax", SIZE_COARSE)
        gmsh.model.mesh.field.setNumber(f_thresh, "DistMin", 0)
        gmsh.model.mesh.field.setNumber(f_thresh, "DistMax", 100)

        # Restrict this field to a box around the cavern so it doesn't
        # propagate to the domain boundary through interlayer surfaces
        f_box = gmsh.model.mesh.field.add("Box")
        gmsh.model.mesh.field.setNumber(f_box, "VIn", SIZE_FINE)
        gmsh.model.mesh.field.setNumber(f_box, "VOut", SIZE_COARSE)
        gmsh.model.mesh.field.setNumber(f_box, "XMin", xc - R_max - 30)
        gmsh.model.mesh.field.setNumber(f_box, "XMax", xc + R_max + 30)
        gmsh.model.mesh.field.setNumber(f_box, "YMin", yc - R_max - 30)
        gmsh.model.mesh.field.setNumber(f_box, "YMax", yc + R_max + 30)
        gmsh.model.mesh.field.setNumber(f_box, "ZMin", z_bot_tip - 20)
        gmsh.model.mesh.field.setNumber(f_box, "ZMax", z_top_tip + 20)
        gmsh.model.mesh.field.setNumber(f_box, "Thickness", 50)

        # Take the maximum (coarsest) of threshold and box — this means
        # refinement only where BOTH fields agree (near cavern AND inside box)
        f_max = gmsh.model.mesh.field.add("Max")
        gmsh.model.mesh.field.setNumbers(f_max, "FieldsList", [f_thresh, f_box])

        gmsh.model.mesh.field.setAsBackgroundMesh(f_max)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


def generate_a5_interlayer():
    """Generate A5 cavern mesh with tilted localized interlayer bands."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("A5_interlayer")

    occ = gmsh.model.occ

    # ── 1. Domain box ────────────────────────────────────────────────────────
    box = occ.addBox(0, 0, 0, Lx, Ly, Lz)

    # ── 2. A5 cavern by revolving profile ────────────────────────────────────
    cavern = build_a5_cavern_occ()
    occ.synchronize()

    # Cut cavern from domain
    all_vols = [v[1] for v in gmsh.model.getEntities(3)]
    cavern_vols = [v for v in all_vols if v != box]
    occ.cut([(3, box)], [(3, v) for v in cavern_vols],
            removeObject=True, removeTool=True)
    occ.synchronize()

    # ── 3. Create tilted interlayer bands ────────────────────────────────────
    band_upper = create_localized_tilted_band(
        UPPER_BAND_Z_AT_CENTER, UPPER_BAND_X_OFFSET,
        INTERLAYER_THICKNESS, DIP_ANGLE_DEG,
        BAND_LENGTH, BAND_WIDTH
    )
    band_lower = create_localized_tilted_band(
        LOWER_BAND_Z_AT_CENTER, LOWER_BAND_X_OFFSET,
        INTERLAYER_THICKNESS, DIP_ANGLE_DEG,
        BAND_LENGTH, BAND_WIDTH
    )
    occ.synchronize()

    # Fragment domain by interlayer bands
    all_vols = [v[1] for v in gmsh.model.getEntities(3)]
    domain_vol = [v for v in all_vols if v not in [band_upper, band_lower]]

    occ.fragment([(3, v) for v in domain_vol],
                 [(3, band_upper), (3, band_lower)],
                 removeObject=True, removeTool=True)
    occ.synchronize()

    # ── 4. Save geometry ─────────────────────────────────────────────────────
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    brep_path = os.path.join(OUTPUT_DIR, "geom.brep")
    gmsh.write(brep_path)
    print(f"  BREP written to: {brep_path}")

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

    # Classify volumes by mass
    dip_rad = math.radians(DIP_ANGLE_DEG)
    INTERLAYER_VOL_THRESHOLD = 2_000_000

    inter1_vols = []
    inter2_vols = []
    salt_vols = []

    for v in all_vols:
        mass = gmsh.model.occ.getMass(3, v)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(3, v)
        cx = (xmin + xmax) / 2.0
        cz = (zmin + zmax) / 2.0
        print(f"    Vol {v}: mass={mass:.0f} m3")

        if mass < INTERLAYER_VOL_THRESHOLD:
            z_upper = UPPER_BAND_Z_AT_CENTER - (cx - xc - UPPER_BAND_X_OFFSET) * math.tan(dip_rad)
            z_lower = LOWER_BAND_Z_AT_CENTER - (cx - xc - LOWER_BAND_X_OFFSET) * math.tan(dip_rad)
            d_upper = abs(cz - z_upper)
            d_lower = abs(cz - z_lower)

            if d_upper < d_lower:
                inter1_vols.append(v)
                print(f"    -> INTERLAYER_1 (upper, mass={mass:.0f})")
            else:
                inter2_vols.append(v)
                print(f"    -> INTERLAYER_2 (lower, mass={mass:.0f})")
        else:
            salt_vols.append(v)

    print(f"  Classification: {len(salt_vols)} salt, "
          f"{len(inter1_vols)} interlayer_1, {len(inter2_vols)} interlayer_2")

    if salt_vols:
        gmsh.model.addPhysicalGroup(3, salt_vols, 31, "Salt_bottom")
    if inter1_vols:
        gmsh.model.addPhysicalGroup(3, inter1_vols, 32, "Interlayer_1")
    if inter2_vols:
        gmsh.model.addPhysicalGroup(3, inter2_vols, 34, "Interlayer_2")

    print(f"  Surfaces: Bottom={len(bottom_s)}, Top={len(top_s)}, "
          f"South={len(south_s)}, North={len(north_s)}, "
          f"West={len(west_s)}, East={len(east_s)}, Cavern={len(cavern_s)}")

    # ── 6. Mesh refinement ───────────────────────────────────────────────────
    setup_mesh_refinement(cavern_s, all_surfs)

    # ── 7. Generate mesh ─────────────────────────────────────────────────────
    gmsh.model.mesh.generate(3)

    # ── 8. Save mesh ─────────────────────────────────────────────────────────
    output_path = os.path.join(OUTPUT_DIR, "geom.msh")
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(output_path)

    node_tags, _, _ = gmsh.model.mesh.getNodes()
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_tets = sum(len(t) for t in elem_tags)
    print(f"  Mesh: {len(node_tags)} nodes, {n_tets} tetrahedra")
    print(f"  Written to: {output_path}")

    # ── 9. Write .geo with mesh fields for visual inspection ─────────────────
    # Collect cavern geometry point tags
    cavern_pt_tags = set()
    for s in cavern_s:
        boundary = gmsh.model.getBoundary([(2, s)], combined=False, oriented=False)
        for b in boundary:
            pts = gmsh.model.getBoundary([b], combined=False, oriented=False)
            for p in pts:
                cavern_pt_tags.add(p[1])

    geo_path = os.path.join(OUTPUT_DIR, "geom.geo")
    with open(geo_path, "w") as f:
        f.write('// A5 Cavern with tilted localized interlayer bands\n')
        f.write('// Generated by generate_A5_heterogeneous_tilted.py\n')
        f.write('SetFactory("OpenCASCADE");\n')
        f.write('Merge "geom.brep";\n\n')
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
        f.write(f'Field[3].XMin = {xc - R_max - 30};\n')
        f.write(f'Field[3].XMax = {xc + R_max + 30};\n')
        f.write(f'Field[3].YMin = {yc - R_max - 30};\n')
        f.write(f'Field[3].YMax = {yc + R_max + 30};\n')
        f.write(f'Field[3].ZMin = {z_bot_tip - 20};\n')
        f.write(f'Field[3].ZMax = {z_top_tip + 20};\n')
        f.write(f'Field[3].Thickness = 50;\n\n')
        # Max of both: refinement only where near cavern AND inside box
        f.write(f'Field[4] = Max;\n')
        f.write(f'Field[4].FieldsList = {{2, 3}};\n')
        f.write(f'Background Field = 4;\n')
    print(f"  .geo written to: {geo_path}")

    gmsh.finalize()


def main():
    print("=" * 70)
    print("A5 CAVERN - TILTED INTERLAYER BAND GENERATOR")
    print("=" * 70)
    print(f"  Cavern: A5 Zuidwending (~965k m3)")
    print(f"  z range: {z_bot_tip} - {z_top_tip} ({z_top_tip - z_bot_tip}m)")
    print(f"  R_max: {R_max} m")
    print(f"  Dip angle: {DIP_ANGLE_DEG} deg")
    print(f"  Upper band: z={UPPER_BAND_Z_AT_CENTER}, x_offset={UPPER_BAND_X_OFFSET}")
    print(f"  Lower band: z={LOWER_BAND_Z_AT_CENTER}, x_offset={LOWER_BAND_X_OFFSET}")
    print(f"  Thickness: {INTERLAYER_THICKNESS} m")
    print(f"  Output: {OUTPUT_DIR}")
    print()

    generate_a5_interlayer()

    print()
    print("=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
