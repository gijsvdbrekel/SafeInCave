"""
plot_results_subsample.py — Same as plot_results.py but only reads every Nth
timestep from the XDMF files, skipping the rest entirely.

Use this when the simulation saved every timestep instead of every 15th,
making the standard plot_results.py too slow.

Usage:
  python plot_results_subsample.py
"""

import meshio as ms
import numpy as np
import safeincave.PostProcessingTools as post

# ── SUBSAMPLE INTERVAL ────────────────────────────────────────────────────────
# Change this to control how many steps to skip.
# STEP=15 means only read steps 0, 15, 30, ... (matching the usual save interval).
STEP = 15


# ── Replacement read functions that skip steps at the I/O level ───────────────

def _read_node_vector_sub(xdmf_field_path):
    reader = ms.xdmf.TimeSeriesReader(xdmf_field_path)
    points, cells = reader.read_points_cells()
    n_nodes = points.shape[0]

    keep = list(range(0, reader.num_steps, STEP))
    vector_field = np.zeros((len(keep), n_nodes, 3))
    time_list = np.zeros(len(keep))

    out_idx = 0
    for k in range(reader.num_steps):
        time, point_data, cell_data = reader.read_data(k)
        if k in keep:
            time_list[out_idx] = time
            field_name = list(point_data.keys())[0]
            vector_field[out_idx, :, :] = point_data[field_name]
            out_idx += 1

    return points, time_list, vector_field


def _read_cell_scalar_sub(xdmf_field_path):
    reader = ms.xdmf.TimeSeriesReader(xdmf_field_path)
    points, cells = reader.read_points_cells()
    n_cells = cells["tetra"].shape[0]

    keep = list(range(0, reader.num_steps, STEP))
    scalar_field = np.zeros((len(keep), n_cells))
    time_list = np.zeros(len(keep))

    out_idx = 0
    for k in range(reader.num_steps):
        time, point_data, cell_data = reader.read_data(k)
        if k in keep:
            time_list[out_idx] = time
            field_name = list(cell_data["tetra"].keys())[0]
            scalar_field[out_idx] = cell_data["tetra"][field_name].flatten()
            out_idx += 1

    return post.compute_cell_centroids(cells["tetra"], points), time_list, scalar_field


def _read_cell_tensor_sub(xdmf_field_path):
    reader = ms.xdmf.TimeSeriesReader(xdmf_field_path)
    points, cells = reader.read_points_cells()
    n_cells = cells["tetra"].shape[0]

    keep = list(range(0, reader.num_steps, STEP))
    tensor_field = np.zeros((len(keep), n_cells, 3, 3))
    time_list = np.zeros(len(keep))

    out_idx = 0
    for k in range(reader.num_steps):
        time, point_data, cell_data = reader.read_data(k)
        if k in keep:
            time_list[out_idx] = time
            field_name = list(cell_data["tetra"].keys())[0]
            tensor_field[out_idx, :, :] = cell_data["tetra"][field_name].reshape((n_cells, 3, 3))
            out_idx += 1

    return post.compute_cell_centroids(cells["tetra"], points), time_list, tensor_field


# ── Monkey-patch before importing plot_results ────────────────────────────────
post.read_node_vector = _read_node_vector_sub
post.read_cell_scalar = _read_cell_scalar_sub
post.read_cell_tensor = _read_cell_tensor_sub

# ── Now import and run the original main ──────────────────────────────────────
from plot_results import main  # noqa: E402

if __name__ == "__main__":
    main()
