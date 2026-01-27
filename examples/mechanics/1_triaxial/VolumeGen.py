#!/usr/bin/env python3
"""
WORKING: Cavern volume calculator using direct tetrahedral formula.
Final version with ALL API issues fixed.
"""
import os
import glob
import csv
import numpy as np
import gmsh

DIR = "/home/gvandenbrekel/SafeInCave/VolumeCalculation/Gmsh_600k_3D"
OUT_CSV = os.path.join(DIR, "cavern_volumes.csv")


def tetrahedron_volume(p0, p1, p2, p3):
    """Calculate volume of tetrahedron from 4 vertices."""
    v1 = p1 - p0
    v2 = p2 - p0
    v3 = p3 - p0
    det = np.dot(np.cross(v1, v2), v3)
    return abs(det) / 6.0


def find_salt_tag():
    """Find Salt physical volume tag."""
    phys = gmsh.model.getPhysicalGroups(dim=3)
    
    for d, tag in phys:
        try:
            name = gmsh.model.getPhysicalName(d, tag)
            if name == "Salt":
                return tag
        except:
            pass
    
    for d, tag in phys:
        if tag == 28:
            return tag
    
    vols = [tag for d, tag in phys if d == 3]
    if len(vols) == 1:
        return vols[0]
    
    raise RuntimeError("Cannot find Salt volume")


def calculate_salt_volume(phys_tag: int) -> float:
    """
    Calculate volume by summing tetrahedral elements.
    FIXED: Correct getNodes() API call.
    """
    vol_entities = gmsh.model.getEntitiesForPhysicalGroup(3, phys_tag)
    
    total_volume = 0.0
    total_elements = 0
    
    for entity_tag in vol_entities:
        # Get all elements
        elem_types, elem_tags_list, node_tags_list = gmsh.model.mesh.getElements(3, entity_tag)
        
        for elem_type, elem_tags, node_tags in zip(elem_types, elem_tags_list, node_tags_list):
            if len(elem_tags) == 0:
                continue
            
            # Get element properties
            elem_name, _, _, n_nodes, _, _ = gmsh.model.mesh.getElementProperties(elem_type)
            
            # Reshape: (n_elements, n_nodes_per_element)
            n_elements = len(elem_tags)
            node_tags_array = np.array(node_tags, dtype=int).reshape(n_elements, n_nodes)
            
            # For 2nd order tets (10 nodes), use only corners (first 4)
            if n_nodes == 10:
                node_tags_array = node_tags_array[:, :4]
            elif n_nodes != 4:
                print(f"    Warning: Skipping element type {elem_name} (n_nodes={n_nodes})")
                continue
            
            # Get ALL mesh nodes at once (more efficient)
            # FIXED API call: getNodes() without arguments gets all nodes
            all_tags, all_coords, _ = gmsh.model.mesh.getNodes()
            all_coords = np.array(all_coords).reshape(-1, 3)
            
            # Create map: node_tag -> coordinate index
            tag_to_idx = {int(tag): idx for idx, tag in enumerate(all_tags)}
            
            # Calculate volume for each tetrahedron
            for elem_node_tags in node_tags_array:
                try:
                    # Get indices for these 4 nodes
                    indices = [tag_to_idx[int(tag)] for tag in elem_node_tags[:4]]
                    p0, p1, p2, p3 = all_coords[indices]
                    
                    # Calculate volume
                    vol = tetrahedron_volume(p0, p1, p2, p3)
                    total_volume += vol
                    total_elements += 1
                except KeyError as e:
                    print(f"    Warning: Node tag {e} not found, skipping element")
                    continue
    
    print(f"    Computed from {total_elements:,} tetrahedral elements")
    return float(total_volume)


def get_box_volume(dim: int, entity_tags) -> float:
    """Get bounding box volume."""
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
        
        # Sync
        try:
            gmsh.model.geo.synchronize()
        except:
            try:
                gmsh.model.occ.synchronize()
            except:
                pass
        
        # Mesh settings
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 50)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 2)
        
        # Generate
        print(f"  Generating mesh...")
        gmsh.model.mesh.generate(3)
        
        # Stats
        nodes, _, _ = gmsh.model.mesh.getNodes()
        elems = gmsh.model.mesh.getElements(3)
        n_elem = sum(len(t) for t in elems[1])
        print(f"  Mesh: {len(nodes):,} nodes, {n_elem:,} elements")
        
        # Find Salt
        salt_tag = find_salt_tag()
        print(f"  Found Salt (tag={salt_tag})")
        
        # Calculate
        print(f"  Computing volumes...")
        v_salt = calculate_salt_volume(salt_tag)
        
        salt_entities = gmsh.model.getEntitiesForPhysicalGroup(3, salt_tag)
        v_box = get_box_volume(3, salt_entities)
        
        v_cavern = v_box - v_salt
        
        print(f"  Results:")
        print(f"    Domain:  {v_box:>15,.0f} m³")
        print(f"    Salt:    {v_salt:>15,.0f} m³")
        print(f"    Cavern:  {v_cavern:>15,.0f} m³  ({v_cavern/1e6:.3f} Mm³)")
        print(f"{'='*70}")
        
        return v_cavern, v_salt, v_box
        
    finally:
        gmsh.finalize()


def main():
    """Main function."""
    print(f"\n{'='*70}")
    print(f"CAVERN VOLUME CALCULATOR")
    print(f"Directory: {DIR}")
    print(f"{'='*70}")
    
    files = sorted(glob.glob(os.path.join(DIR, "*.geo")))
    
    if not files:
        print(f"\nERROR: No .geo files found")
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
    print(f"SUMMARY")
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
            status = row[4][:35] if len(row[4]) > 35 else row[4]
            print(f"{name:<30} {'FAILED':<15} {'N/A':<12} ✗")
    
    print(f"{'='*70}")
    print(f"\nSuccess: {n_ok}/{len(results)} files")
    print(f"CSV: {OUT_CSV}\n")


if __name__ == "__main__":
    main()
