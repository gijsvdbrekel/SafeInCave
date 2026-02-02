#!/usr/bin/env python3
"""
Cavern volume calculator using direct tetrahedral formula.

Robust for:
- Old caverns: single Physical Volume("Salt", <tag>) (e.g., tag 28)
- New caverns with interlayers: multiple salt parts (e.g., "Salt_top") and interlayers
  Cavern volume is computed as:
      cavern = bounding_box(domain_solids) - volume(all_solids)
  where "all_solids" are ALL 3D physical volumes present in the model.
"""

import os
import glob
import csv
import numpy as np
import gmsh

DIR = "/home/gvandenbrekel/SafeInCave/VolumeCalculation/Gmsh_600k_interlayer"
OUT_CSV = os.path.join(DIR, "cavern_volumes.csv")


def tetrahedron_volume(p0, p1, p2, p3):
    """Calculate volume of tetrahedron from 4 vertices."""
    v1 = p1 - p0
    v2 = p2 - p0
    v3 = p3 - p0
    det = np.dot(np.cross(v1, v2), v3)
    return abs(det) / 6.0


def list_physical_volumes_3d():
    """Return list of (tag, name) for all 3D physical groups."""
    phys = gmsh.model.getPhysicalGroups(dim=3)
    out = []
    for d, tag in phys:
        if d != 3:
            continue
        name = ""
        try:
            name = gmsh.model.getPhysicalName(d, tag) or ""
        except Exception:
            pass
        out.append((tag, name))
    return out


def find_tags_by_name_contains(substr: str):
    """Find 3D physical group tags whose name contains substr (case-insensitive)."""
    substr = substr.lower()
    tags = []
    for tag, name in list_physical_volumes_3d():
        if substr in (name or "").lower():
            tags.append(tag)
    return tags


def find_all_solid_volume_tags():
    """
    In your .geo files, the solids (salt and possible interlayers) are represented
    by ALL 3D Physical Volumes. This is the most robust assumption across files.

    If a model has NO 3D physical groups, we fail loudly (otherwise you'd compute nonsense).
    """
    tags = [tag for tag, _ in list_physical_volumes_3d()]
    if not tags:
        raise RuntimeError("No 3D Physical Volumes found (dim=3). Cannot define solids.")
    return tags


def get_all_nodes_map():
    """Get node tag -> coordinate mapping once (fast)."""
    all_tags, all_coords, _ = gmsh.model.mesh.getNodes()
    all_tags = np.array(all_tags, dtype=int)
    all_coords = np.array(all_coords, dtype=float).reshape(-1, 3)
    tag_to_idx = {int(tag): idx for idx, tag in enumerate(all_tags)}
    return all_coords, tag_to_idx


def calculate_volume_for_phys_tags(phys_tags) -> float:
    """
    Calculate volume by summing tetrahedral elements for one or more Physical Volume tags.
    Supports 1st order tets (4 nodes) and 2nd order tets (10 nodes; uses corner nodes only).
    """
    if isinstance(phys_tags, int):
        phys_tags = [phys_tags]
    phys_tags = list(phys_tags)

    all_coords, tag_to_idx = get_all_nodes_map()

    total_volume = 0.0
    total_elements = 0

    for phys_tag in phys_tags:
        vol_entities = gmsh.model.getEntitiesForPhysicalGroup(3, phys_tag)

        for entity_tag in vol_entities:
            elem_types, elem_tags_list, node_tags_list = gmsh.model.mesh.getElements(3, entity_tag)

            for elem_type, elem_tags, node_tags in zip(elem_types, elem_tags_list, node_tags_list):
                if len(elem_tags) == 0:
                    continue

                elem_name, _, _, n_nodes, _, _ = gmsh.model.mesh.getElementProperties(elem_type)

                n_elements = len(elem_tags)
                node_tags_array = np.array(node_tags, dtype=int).reshape(n_elements, n_nodes)

                # For 2nd order tets (10 nodes), use only corners (first 4)
                if n_nodes == 10:
                    node_tags_array = node_tags_array[:, :4]
                elif n_nodes != 4:
                    print(f"    Warning: Skipping element type {elem_name} (n_nodes={n_nodes})")
                    continue

                # Sum tet volumes
                for elem_node_tags in node_tags_array:
                    try:
                        indices = [tag_to_idx[int(t)] for t in elem_node_tags[:4]]
                        p0, p1, p2, p3 = all_coords[indices]
                        total_volume += tetrahedron_volume(p0, p1, p2, p3)
                        total_elements += 1
                    except KeyError as e:
                        print(f"    Warning: Node tag {e} not found, skipping element")
                        continue

    print(f"    Computed from {total_elements:,} tetrahedral elements")
    return float(total_volume)


def get_box_volume(dim: int, entity_tags) -> float:
    """Get bounding box volume over a list of entities of given dimension."""
    if not entity_tags:
        return 0.0

    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, entity_tags[0])

    for tag in entity_tags[1:]:
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, zmin = min(xmin, x0), min(ymin, y0), min(zmin, z0)
        xmax, ymax, zmax = max(xmax, x1), max(ymax, y1), max(zmax, z1)

    return (xmax - xmin) * (ymax - ymin) * (zmax - zmin)


def process_file(filepath: str):
    """Process one .geo file."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.Verbosity", 2)

    try:
        basename = os.path.basename(filepath)
        print(f"\n{'='*70}")
        print(f"Processing: {basename}")

        # Load geometry - set tolerances BEFORE opening
        gmsh.option.setNumber("Geometry.Tolerance", 1e-5)
        gmsh.option.setNumber("Geometry.ToleranceBoolean", 1e-5)

        gmsh.open(filepath)

        # Additional geometry fixes
        gmsh.option.setNumber("Geometry.OCCFixDegenerated", 1)
        gmsh.option.setNumber("Geometry.OCCFixSmallEdges", 1)
        gmsh.option.setNumber("Geometry.OCCFixSmallFaces", 1)

        # Sync (geo vs occ)
        try:
            gmsh.model.geo.synchronize()
        except Exception:
            try:
                gmsh.model.occ.synchronize()
            except Exception:
                pass

        # Mesh settings
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 50)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 2)

        # Generate mesh
        print("  Generating mesh...")
        gmsh.model.mesh.generate(3)

        # Stats
        nodes, _, _ = gmsh.model.mesh.getNodes()
        elems = gmsh.model.mesh.getElements(3)
        n_elem = sum(len(t) for t in elems[1])
        print(f"  Mesh: {len(nodes):,} nodes, {n_elem:,} elements")

        # Inspect physical volumes
        phys3d = list_physical_volumes_3d()
        if not phys3d:
            raise RuntimeError("No 3D Physical Volumes found. Check .geo Physical Volume definitions.")

        print("  3D Physical Volumes found:")
        for tag, name in phys3d:
            nm = name if name else "(no name)"
            print(f"    - tag={tag:>3}  name={nm}")

        # Solids are ALL 3D physical volumes (salt + any interlayers)
        solid_tags = find_all_solid_volume_tags()

        # Salt tags: handle both old ("Salt") and new ("Salt_top", etc.)
        salt_tags = find_tags_by_name_contains("salt")

        # If someone forgot to name things but old tag 28 exists, keep that as a last fallback for reporting salt only.
        # (Cavern volume does NOT rely on this.)
        if not salt_tags:
            if any(tag == 28 for tag, _ in phys3d):
                salt_tags = [28]

        print(f"  Solid tags (all 3D phys): {solid_tags}")
        if salt_tags:
            print(f"  Salt tags (name contains 'salt'): {salt_tags}")
        else:
            print("  Salt tags: none found (salt volume will be reported as 0).")

        # Compute volumes
        print("  Computing volumes...")
        v_solid = calculate_volume_for_phys_tags(solid_tags)
        v_salt = calculate_volume_for_phys_tags(salt_tags) if salt_tags else 0.0

        # Domain box from ALL solid entities (robust)
        solid_entities = []
        for t in solid_tags:
            solid_entities.extend(gmsh.model.getEntitiesForPhysicalGroup(3, t))
        v_box = get_box_volume(3, solid_entities)

        # Cavern = bounding box - solids
        v_cavern = v_box - v_solid

        # Hard sanity check
        if v_cavern < 0:
            raise RuntimeError(
                f"Negative cavern volume ({v_cavern:.6e}). "
                "Likely non-closed volume, bad physical groups, or a bounding-box mismatch."
            )

        print("  Results:")
        print(f"    Domain:  {v_box:>15,.0f} m³")
        print(f"    Solids:  {v_solid:>15,.0f} m³")
        print(f"    Salt:    {v_salt:>15,.0f} m³")
        print(f"    Cavern:  {v_cavern:>15,.0f} m³  ({v_cavern/1e6:.3f} Mm³)")
        print(f"{'='*70}")

        # Return cavern, salt, domain (keep same CSV columns as before)
        return v_cavern, v_salt, v_box

    finally:
        gmsh.finalize()


def main():
    """Main function."""
    print(f"\n{'='*70}")
    print("CAVERN VOLUME CALCULATOR")
    print(f"Directory: {DIR}")
    print(f"{'='*70}")

    files = sorted(glob.glob(os.path.join(DIR, "*.geo")))

    if not files:
        print("\nERROR: No .geo files found")
        return

    print(f"\nFound {len(files)} files")

    results = []

    for filepath in files:
        basename = os.path.basename(filepath)

        try:
            v_cav, v_salt, v_box = process_file(filepath)
            results.append([
                basename,
                f"{v_cav:.2f}",
                f"{v_salt:.2f}",
                f"{v_box:.2f}",
                "OK"
            ])
        except Exception as e:
            print(f"\n{'='*70}")
            print(f"FAILED: {basename}")
            print(f"Error: {e}")
            print(f"{'='*70}")

            error_msg = str(e)[:80]
            results.append([basename, "", "", "", f"ERROR: {error_msg}"])

    # Write CSV
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["File", "Cavern_m3", "Salt_m3", "Domain_m3", "Status"])
        w.writerows(results)

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"{'Geometry':<30} {'Cavern (m³)':<15} {'Cavern (Mm³)':<12} {'Status'}")
    print(f"{'-'*70}")

    n_ok = 0
    for row in results:
        name = row[0].replace('_full3d.geo', '').replace('cavern_', '').replace('_', ' ').title()

        if row[1]:
            vcav = float(row[1])
            print(f"{name:<30} {vcav:>14,.0f}  {vcav/1e6:>11.3f}  ✓")
            n_ok += 1
        else:
            print(f"{name:<30} {'FAILED':<15} {'N/A':<12} ✗")

    print(f"{'='*70}")
    print(f"\nSuccess: {n_ok}/{len(results)} files")
    print(f"CSV: {OUT_CSV}\n")


if __name__ == "__main__":
    main()
