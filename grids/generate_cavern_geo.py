"""
Generator script for gmsh .geo cavern files.
Creates axisymmetric caverns defined by (z, R) profiles inside a 450x450x660 box.
First and last levels must have R=0 (tip points).
"""
import math
import os

def frustum_volume(dz, R1, R2):
    return math.pi / 3.0 * dz * (R1**2 + R1*R2 + R2**2)

def profile_volume(levels):
    V = 0
    for i in range(1, len(levels)):
        z1, R1 = levels[i-1]
        z2, R2 = levels[i]
        V += frustum_volume(z2 - z1, R1, R2)
    return V

def scale_profile(levels, s):
    z_vals = [z for z, r in levels]
    z_mid = (min(z_vals) + max(z_vals)) / 2
    return [(z_mid + s * (z - z_mid), s * r) for z, r in levels]

def fit_volume(levels_base, target_volume=1200000.0, tol=10.0):
    s = 1.0
    levels = levels_base
    V = profile_volume(levels)
    for _ in range(100):
        if abs(V - target_volume) < tol:
            break
        s *= (target_volume / V) ** (1.0 / 3.0)
        levels = scale_profile(levels_base, s)
        V = profile_volume(levels)
    return levels, s, V

def generate_geo(name, description, levels, filename):
    """
    Generate a gmsh .geo file for a cavern.
    levels: list of (z, R) pairs. First and last must have R=0 (tips).
    """
    N = len(levels)
    assert levels[0][1] == 0, "First level must be a tip (R=0)"
    assert levels[-1][1] == 0, "Last level must be a tip (R=0)"

    Lx, Ly, Lz = 450.0, 450.0, 660.0
    xc, yc = Lx / 2, Ly / 2
    coarse = 65
    fine = 4.5

    lines = []
    w = lines.append

    w(f"// ============================================================================")
    w(f"// {name.upper()} - SCALED TO 1.2 MILLION m³")
    w(f"// ============================================================================")
    w(f"// {description}")
    w(f"// Target volume:   1.2 Mm³")
    w(f"// Analytical volume: {profile_volume(levels):.0f} m³")
    w(f"// Domain dimensions: {Lx:.0f}×{Ly:.0f}×{Lz:.0f} m")
    w(f"// ============================================================================")
    w(f"")
    w(f"size_coarse = {coarse};")
    w(f"size_fine   = {fine};")
    w(f"")
    w(f"Lz = {Lz:.1f};")
    w(f"Ly = {Ly:.1f};")
    w(f"Lx = {Lx:.1f};")
    w(f"")
    w(f"x_center = Lx/2;  // {xc:.1f}")
    w(f"y_center = Ly/2;  // {yc:.1f}")
    w(f"")

    # ---- OUTER BOX ----
    w(f"//////////////////////////////")
    w(f"// OUTER BOX")
    w(f"//////////////////////////////")
    w(f"")
    w(f"Point(1) = {{0,  0,  0,  size_coarse}};")
    w(f"Point(2) = {{Lx, 0,  0,  size_coarse}};")
    w(f"Point(3) = {{Lx, Ly, 0,  size_coarse}};")
    w(f"Point(4) = {{0,  Ly, 0,  size_coarse}};")
    w(f"")
    w(f"Point(5) = {{0,  0,  Lz, size_coarse}};")
    w(f"Point(6) = {{Lx, 0,  Lz, size_coarse}};")
    w(f"Point(7) = {{Lx, Ly, Lz, size_coarse}};")
    w(f"Point(8) = {{0,  Ly, Lz, size_coarse}};")
    w(f"")

    # ---- CAVERN GEOMETRY POINTS ----
    # Numbering: level i -> base = 100 + 10*i
    # Tip levels: tip=base, center=base+1, px=base+2, mx=base+3, py=base+4, my=base+5
    # Ring levels: center=base, px=base+2, mx=base+3, py=base+4, my=base+5

    def base(i):
        return 100 + 10 * i

    def pt_tip(i):
        return base(i)

    def pt_center(i):
        if i == 0 or i == N - 1:
            return base(i) + 1
        return base(i)

    def pt_px(i):
        return base(i) + 2

    def pt_mx(i):
        return base(i) + 3

    def pt_py(i):
        return base(i) + 4

    def pt_my(i):
        return base(i) + 5

    w(f"//////////////////////////////")
    w(f"// CAVERN GEOMETRY - {N} LEVELS")
    w(f"//////////////////////////////")
    w(f"")

    for i, (z, R) in enumerate(levels):
        is_tip = (i == 0 or i == N - 1)
        if is_tip:
            tip_label = "bottom" if i == 0 else "top"
            # For bottom tip, tip is below center; for top tip, tip is above center
            if i == 0:
                z_tip = z
                z_ring = z + R  # This doesn't work for R=0...
                # Actually for tips, the adjacent level gives the ring z
                # The tip point is at z, and the ring is at the adjacent level's z
                # But we need a center point for the circle arcs of the ring
                z_next = levels[i + 1][0]
                R_next = levels[i + 1][1]
                # Ring z = z + some small offset toward next level
                dz_ring = min(abs(z_next - z) * 0.5, R_next * 0.5)
                z_ring = z + dz_ring
                R_ring = R_next * 0.4  # Small ring near the tip
            else:
                z_prev = levels[i - 1][0]
                R_prev = levels[i - 1][1]
                dz_ring = min(abs(z - z_prev) * 0.5, R_prev * 0.5)
                z_ring = z - dz_ring
                R_ring = R_prev * 0.4

            w(f"// Level {i+1} - {tip_label} tip (z={z:.2f})")
            w(f"Point({pt_tip(i)}) = {{{xc:.6f}, {yc:.6f}, {z:.6f}, size_fine}};  // tip")
            w(f"Point({pt_center(i)}) = {{{xc:.6f}, {yc:.6f}, {z_ring:.6f}, size_coarse}};  // center")
            w(f"Point({pt_px(i)}) = {{{xc + R_ring:.6f}, {yc:.6f}, {z_ring:.6f}, size_fine}};  // +x")
            w(f"Point({pt_mx(i)}) = {{{xc - R_ring:.6f}, {yc:.6f}, {z_ring:.6f}, size_fine}};  // -x")
            w(f"Point({pt_py(i)}) = {{{xc:.6f}, {yc + R_ring:.6f}, {z_ring:.6f}, size_fine}};  // +y")
            w(f"Point({pt_my(i)}) = {{{xc:.6f}, {yc - R_ring:.6f}, {z_ring:.6f}, size_fine}};  // -y")
        else:
            w(f"// Level {i+1} (z={z:.2f}, R={R:.2f})")
            w(f"Point({pt_center(i)}) = {{{xc:.6f}, {yc:.6f}, {z:.6f}, size_coarse}};  // center")
            w(f"Point({pt_px(i)}) = {{{xc + R:.6f}, {yc:.6f}, {z:.6f}, size_fine}};  // +x")
            w(f"Point({pt_mx(i)}) = {{{xc - R:.6f}, {yc:.6f}, {z:.6f}, size_fine}};  // -x")
            w(f"Point({pt_py(i)}) = {{{xc:.6f}, {yc + R:.6f}, {z:.6f}, size_fine}};  // +y")
            w(f"Point({pt_my(i)}) = {{{xc:.6f}, {yc - R:.6f}, {z:.6f}, size_fine}};  // -y")
        w(f"")

    # ---- OUTER BOX EDGES ----
    w(f"//////////////////////////////")
    w(f"// OUTER BOX EDGES")
    w(f"//////////////////////////////")
    w(f"")
    w(f"Line(1)  = {{1,2}};")
    w(f"Line(2)  = {{2,3}};")
    w(f"Line(3)  = {{3,4}};")
    w(f"Line(4)  = {{4,1}};")
    w(f"Line(5)  = {{5,6}};")
    w(f"Line(6)  = {{6,7}};")
    w(f"Line(7)  = {{7,8}};")
    w(f"Line(8)  = {{8,5}};")
    w(f"Line(9)  = {{1,5}};")
    w(f"Line(10) = {{2,6}};")
    w(f"Line(11) = {{3,7}};")
    w(f"Line(12) = {{4,8}};")
    w(f"")

    # ---- CAVERN CONNECTION LINES ----
    # Line IDs: we need a counter. Start at 20.
    line_id = [20]  # mutable counter

    def next_line():
        lid = line_id[0]
        line_id[0] += 1
        return lid

    # Bottom tip lines (from tip to ring)
    w(f"//////////////////////////////")
    w(f"// CAVERN LINES AND ARCS")
    w(f"//////////////////////////////")
    w(f"")

    # Bottom tip → bottom ring
    w(f"// Bottom tip to ring")
    btip_px = next_line()
    w(f"Line({btip_px}) = {{{pt_tip(0)}, {pt_px(0)}}};")
    btip_py = next_line()
    w(f"Line({btip_py}) = {{{pt_tip(0)}, {pt_py(0)}}};")
    btip_mx = next_line()
    w(f"Line({btip_mx}) = {{{pt_tip(0)}, {pt_mx(0)}}};")
    btip_my = next_line()
    w(f"Line({btip_my}) = {{{pt_tip(0)}, {pt_my(0)}}};")
    w(f"")

    # Bottom ring circles
    w(f"// Bottom ring circles")
    bcirc_px_py = next_line()
    w(f"Circle({bcirc_px_py}) = {{{pt_px(0)}, {pt_center(0)}, {pt_py(0)}}};")
    bcirc_py_mx = next_line()
    w(f"Circle({bcirc_py_mx}) = {{{pt_py(0)}, {pt_center(0)}, {pt_mx(0)}}};")
    bcirc_mx_my = next_line()
    w(f"Circle({bcirc_mx_my}) = {{{pt_mx(0)}, {pt_center(0)}, {pt_my(0)}}};")
    bcirc_my_px = next_line()
    w(f"Circle({bcirc_my_px}) = {{{pt_my(0)}, {pt_center(0)}, {pt_px(0)}}};")
    w(f"")

    # Vertical connections between ring levels
    # For levels 0 to N-2, connect ring i to ring i+1
    # (level 0 ring is the bottom tip's ring, level N-1 ring is the top tip's ring)
    vert_px = {}  # vert_px[(i, i+1)] = line_id
    vert_mx = {}
    vert_py = {}
    vert_my = {}

    # Ring levels: 0 (bottom tip ring), 1, 2, ..., N-2, N-1 (top tip ring)
    w(f"// Vertical connections between levels")
    for i in range(N - 1):
        j = i + 1
        vert_px[(i, j)] = next_line()
        w(f"Line({vert_px[(i,j)]}) = {{{pt_px(i)}, {pt_px(j)}}};  // +x level {i+1}->{j+1}")
        vert_mx[(i, j)] = next_line()
        w(f"Line({vert_mx[(i,j)]}) = {{{pt_mx(i)}, {pt_mx(j)}}};  // -x level {i+1}->{j+1}")
        vert_py[(i, j)] = next_line()
        w(f"Line({vert_py[(i,j)]}) = {{{pt_py(i)}, {pt_py(j)}}};  // +y level {i+1}->{j+1}")
        vert_my[(i, j)] = next_line()
        w(f"Line({vert_my[(i,j)]}) = {{{pt_my(i)}, {pt_my(j)}}};  // -y level {i+1}->{j+1}")
    w(f"")

    # Circle arcs at each intermediate level (1 to N-2)
    circ_px_py = {}
    circ_py_mx = {}
    circ_mx_my = {}
    circ_my_px = {}

    w(f"// Circle arcs at intermediate levels")
    for i in range(1, N - 1):
        circ_px_py[i] = next_line()
        w(f"Circle({circ_px_py[i]}) = {{{pt_px(i)}, {pt_center(i)}, {pt_py(i)}}};  // level {i+1}")
        circ_py_mx[i] = next_line()
        w(f"Circle({circ_py_mx[i]}) = {{{pt_py(i)}, {pt_center(i)}, {pt_mx(i)}}};")
        circ_mx_my[i] = next_line()
        w(f"Circle({circ_mx_my[i]}) = {{{pt_mx(i)}, {pt_center(i)}, {pt_my(i)}}};")
        circ_my_px[i] = next_line()
        w(f"Circle({circ_my_px[i]}) = {{{pt_my(i)}, {pt_center(i)}, {pt_px(i)}}};")
    w(f"")

    # Top tip lines and circles
    w(f"// Top tip to ring")
    t = N - 1
    ttip_px = next_line()
    w(f"Line({ttip_px}) = {{{pt_tip(t)}, {pt_px(t)}}};")
    ttip_py = next_line()
    w(f"Line({ttip_py}) = {{{pt_tip(t)}, {pt_py(t)}}};")
    ttip_mx = next_line()
    w(f"Line({ttip_mx}) = {{{pt_tip(t)}, {pt_mx(t)}}};")
    ttip_my = next_line()
    w(f"Line({ttip_my}) = {{{pt_tip(t)}, {pt_my(t)}}};")
    w(f"")

    w(f"// Top ring circles")
    tcirc_px_py = next_line()
    w(f"Circle({tcirc_px_py}) = {{{pt_px(t)}, {pt_center(t)}, {pt_py(t)}}};")
    tcirc_py_mx = next_line()
    w(f"Circle({tcirc_py_mx}) = {{{pt_py(t)}, {pt_center(t)}, {pt_mx(t)}}};")
    tcirc_mx_my = next_line()
    w(f"Circle({tcirc_mx_my}) = {{{pt_mx(t)}, {pt_center(t)}, {pt_my(t)}}};")
    tcirc_my_px = next_line()
    w(f"Circle({tcirc_my_px}) = {{{pt_my(t)}, {pt_center(t)}, {pt_px(t)}}};")
    w(f"")

    # ---- CAVERN SURFACES ----
    surf_id = [201]

    def next_surf():
        sid = surf_id[0]
        surf_id[0] += 1
        return sid

    cavern_surfaces = []

    w(f"//////////////////////////////")
    w(f"// CAVERN SURFACES")
    w(f"//////////////////////////////")
    w(f"")

    # Bottom cap (4 triangular surfaces from tip to bottom ring)
    w(f"// Bottom cap")
    s1 = next_surf()
    w(f"Curve Loop({s1}) = {{{btip_px}, {bcirc_px_py}, -{btip_py}}};")
    w(f"Surface({s1}) = {{{s1}}};")
    cavern_surfaces.append(s1)

    s2 = next_surf()
    w(f"Curve Loop({s2}) = {{{btip_py}, {bcirc_py_mx}, -{btip_mx}}};")
    w(f"Surface({s2}) = {{{s2}}};")
    cavern_surfaces.append(s2)

    s3 = next_surf()
    w(f"Curve Loop({s3}) = {{{btip_mx}, {bcirc_mx_my}, -{btip_my}}};")
    w(f"Surface({s3}) = {{{s3}}};")
    cavern_surfaces.append(s3)

    s4 = next_surf()
    w(f"Curve Loop({s4}) = {{{btip_my}, {bcirc_my_px}, -{btip_px}}};")
    w(f"Surface({s4}) = {{{s4}}};")
    cavern_surfaces.append(s4)
    w(f"")

    # Side walls: bottom ring to level 1
    # Between level 0 (bottom ring) and level 1
    w(f"// Side walls: bottom ring to level 2")
    # quadrant 1: +x to +y
    s = next_surf()
    w(f"Curve Loop({s}) = {{{bcirc_px_py}, {vert_py[(0,1)]}, -{circ_px_py[1]}, -{vert_px[(0,1)]}}};")
    w(f"Surface({s}) = {{{s}}};")
    cavern_surfaces.append(s)

    s = next_surf()
    w(f"Curve Loop({s}) = {{{bcirc_py_mx}, {vert_mx[(0,1)]}, -{circ_py_mx[1]}, -{vert_py[(0,1)]}}};")
    w(f"Surface({s}) = {{{s}}};")
    cavern_surfaces.append(s)

    s = next_surf()
    w(f"Curve Loop({s}) = {{{bcirc_mx_my}, {vert_my[(0,1)]}, -{circ_mx_my[1]}, -{vert_mx[(0,1)]}}};")
    w(f"Surface({s}) = {{{s}}};")
    cavern_surfaces.append(s)

    s = next_surf()
    w(f"Curve Loop({s}) = {{{bcirc_my_px}, {vert_px[(0,1)]}, -{circ_my_px[1]}, -{vert_my[(0,1)]}}};")
    w(f"Surface({s}) = {{{s}}};")
    cavern_surfaces.append(s)
    w(f"")

    # Side walls between intermediate levels (and last intermediate to top tip ring)
    for i in range(1, N - 1):
        j = i + 1
        w(f"// Side walls: level {i+1} to level {j+1}")

        # Use top circles for level i, bottom circles for level j
        if j == N - 1:
            # j is the top tip level - use top ring circles
            cj_px_py = tcirc_px_py
            cj_py_mx = tcirc_py_mx
            cj_mx_my = tcirc_mx_my
            cj_my_px = tcirc_my_px
        else:
            cj_px_py = circ_px_py[j]
            cj_py_mx = circ_py_mx[j]
            cj_mx_my = circ_mx_my[j]
            cj_my_px = circ_my_px[j]

        s = next_surf()
        w(f"Curve Loop({s}) = {{{circ_px_py[i]}, {vert_py[(i,j)]}, -{cj_px_py}, -{vert_px[(i,j)]}}};")
        w(f"Surface({s}) = {{{s}}};")
        cavern_surfaces.append(s)

        s = next_surf()
        w(f"Curve Loop({s}) = {{{circ_py_mx[i]}, {vert_mx[(i,j)]}, -{cj_py_mx}, -{vert_py[(i,j)]}}};")
        w(f"Surface({s}) = {{{s}}};")
        cavern_surfaces.append(s)

        s = next_surf()
        w(f"Curve Loop({s}) = {{{circ_mx_my[i]}, {vert_my[(i,j)]}, -{cj_mx_my}, -{vert_mx[(i,j)]}}};")
        w(f"Surface({s}) = {{{s}}};")
        cavern_surfaces.append(s)

        s = next_surf()
        w(f"Curve Loop({s}) = {{{circ_my_px[i]}, {vert_px[(i,j)]}, -{cj_my_px}, -{vert_my[(i,j)]}}};")
        w(f"Surface({s}) = {{{s}}};")
        cavern_surfaces.append(s)
        w(f"")

    # Top cap (4 triangular surfaces from top ring to top tip)
    w(f"// Top cap")
    s = next_surf()
    w(f"Curve Loop({s}) = {{{ttip_px}, {tcirc_px_py}, -{ttip_py}}};")
    w(f"Surface({s}) = {{{s}}};")
    cavern_surfaces.append(s)

    s = next_surf()
    w(f"Curve Loop({s}) = {{{ttip_py}, {tcirc_py_mx}, -{ttip_mx}}};")
    w(f"Surface({s}) = {{{s}}};")
    cavern_surfaces.append(s)

    s = next_surf()
    w(f"Curve Loop({s}) = {{{ttip_mx}, {tcirc_mx_my}, -{ttip_my}}};")
    w(f"Surface({s}) = {{{s}}};")
    cavern_surfaces.append(s)

    s = next_surf()
    w(f"Curve Loop({s}) = {{{ttip_my}, {tcirc_my_px}, -{ttip_px}}};")
    w(f"Surface({s}) = {{{s}}};")
    cavern_surfaces.append(s)
    w(f"")

    # ---- OUTER BOX SURFACES ----
    w(f"//////////////////////////////")
    w(f"// OUTER BOX SURFACES")
    w(f"//////////////////////////////")
    w(f"")
    w(f"// Bottom (z=0)")
    w(f"Curve Loop(1) = {{1, 2, 3, 4}};")
    w(f"Plane Surface(1) = {{1}};")
    w(f"")
    w(f"// Top (z=Lz)")
    w(f"Curve Loop(2) = {{5, 6, 7, 8}};")
    w(f"Plane Surface(2) = {{2}};")
    w(f"")
    w(f"// South (y=0)")
    w(f"Curve Loop(3) = {{1, 10, -5, -9}};")
    w(f"Plane Surface(3) = {{3}};")
    w(f"")
    w(f"// North (y=Ly)")
    w(f"Curve Loop(4) = {{3, 12, -7, -11}};")
    w(f"Plane Surface(4) = {{4}};")
    w(f"")
    w(f"// West (x=0)")
    w(f"Curve Loop(5) = {{4, 9, -8, -12}};")
    w(f"Plane Surface(5) = {{5}};")
    w(f"")
    w(f"// East (x=Lx)")
    w(f"Curve Loop(6) = {{2, 11, -6, -10}};")
    w(f"Plane Surface(6) = {{6}};")
    w(f"")

    # ---- VOLUME ----
    w(f"//////////////////////////////")
    w(f"// VOLUME")
    w(f"//////////////////////////////")
    w(f"")
    surf_str = ", ".join(str(s) for s in cavern_surfaces)
    w(f"// Cavern surface loop ({len(cavern_surfaces)} surfaces)")
    w(f"Surface Loop(200) = {{{surf_str}}};")
    w(f"")
    w(f"// Outer box surface loop")
    w(f"Surface Loop(201) = {{1, 2, 3, 4, 5, 6}};")
    w(f"")
    w(f"// Salt volume (box minus cavern)")
    w(f"Volume(1) = {{201, 200}};")
    w(f"")

    # ---- PHYSICAL GROUPS ----
    w(f"//////////////////////////////")
    w(f"// PHYSICAL GROUPS")
    w(f"//////////////////////////////")
    w(f"")
    w(f"Physical Surface(\"Bottom\", 27) = {{1}};")
    w(f"Physical Surface(\"Top\",    22) = {{2}};")
    w(f"Physical Surface(\"South\",  23) = {{3}};")
    w(f"Physical Surface(\"North\",  24) = {{4}};")
    w(f"Physical Surface(\"West\",   26) = {{5}};")
    w(f"Physical Surface(\"East\",   25) = {{6}};")
    w(f"")
    w(f"Physical Surface(\"Cavern\", 29) = {{{surf_str}}};")
    w(f"Physical Volume(\"Salt\",    28) = {{1}};")
    w(f"")

    # Wall profile: +x vertical lines from top tip down through all +x connections to bottom tip
    # This traces the cavern profile on the +x side
    wall_lines = [ttip_px]
    for i in range(N - 2, 0, -1):
        wall_lines.append(vert_px[(i, i + 1)] if i < N - 1 else None)
    # Actually: top tip line, then vertical connections from top to bottom, then bottom tip line
    wall_lines = [ttip_px]
    for i in range(N - 2, 0, -1):
        wall_lines.append(vert_px[(i, i + 1)])
    wall_lines.append(vert_px[(0, 1)])
    wall_lines.append(btip_px)
    wall_str = ", ".join(str(l) for l in wall_lines)
    w(f"Physical Curve(\"Wall_profile\", 30) = {{{wall_str}}};")
    w(f"")

    content = "\n".join(lines)
    with open(filename, 'w') as f:
        f.write(content)
    print(f"Written {filename} ({len(lines)} lines, {len(cavern_surfaces)} cavern surfaces)")
    return content


# ==========================================
# FAST-LEACHED CAVERN BASE PROFILE
# ==========================================
# Rough barrel shape with oscillating radii (jagged walls)
# NOTE: Tube-failure cavern uses multi-chamber structure (hand-written .geo)
FL_BASE = [
    (190, 0),      # bottom tip
    (198, 22),     # small bottom
    (208, 38),     # smooth transition
    (218, 52),     # bulge out
    (228, 48),     # smooth transition
    (238, 40),     # neck in
    (248, 49),     # smooth transition
    (258, 55),     # bulge out
    (268, 50),     # smooth transition
    (278, 42),     # neck in
    (288, 49),     # smooth transition
    (298, 53),     # bulge out
    (308, 49),     # smooth transition
    (318, 43),     # neck in
    (328, 47),     # smooth transition
    (338, 50),     # bulge out
    (348, 45),     # smooth transition
    (358, 38),     # neck in
    (366, 43),     # smooth transition
    (375, 46),     # bulge out
    (382, 30),     # smooth taper
    (388, 12),     # narrow top
    (398, 0),      # top tip
]


if __name__ == '__main__':
    grids_dir = os.path.dirname(os.path.abspath(__file__))

    # Fast-leached 1200k
    fl_levels_1200, fl_scale, fl_vol = fit_volume(FL_BASE, target_volume=1200000.0)
    print(f"Fast-leached 1200k: V={fl_vol:.0f} m³, scale={fl_scale:.6f}")
    fl_dir = os.path.join(grids_dir, "cavern_fastleached_1200_3D")
    os.makedirs(fl_dir, exist_ok=True)
    generate_geo(
        "Fast-Leached Cavern",
        "Rough barrel shape with oscillating radii simulating fast/uncontrolled leaching",
        fl_levels_1200,
        os.path.join(fl_dir, "fastleached_1200k.geo")
    )

    # Fast-leached 600k
    fl_levels_600, fl_scale, fl_vol = fit_volume(FL_BASE, target_volume=600000.0)
    print(f"Fast-leached 600k: V={fl_vol:.0f} m³, scale={fl_scale:.6f}")
    fl_dir = os.path.join(grids_dir, "cavern_fastleached_600_3D")
    os.makedirs(fl_dir, exist_ok=True)
    generate_geo(
        "Fast-Leached Cavern - 600k",
        "Rough barrel shape with oscillating radii simulating fast/uncontrolled leaching",
        fl_levels_600,
        os.path.join(fl_dir, "fastleached_full3d.geo")
    )
