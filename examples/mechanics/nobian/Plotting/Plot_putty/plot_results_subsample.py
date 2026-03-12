"""
plot_results_subsample.py — Same as plot_results.py but subsamples every Nth
timestep from the XDMF data before any processing.

Use this when the simulation saved every timestep instead of every 15th,
making the standard plot_results.py too slow.

Usage:
  python plot_results_subsample.py
"""

import numpy as np
import safeincave.PostProcessingTools as post

# ── SUBSAMPLE INTERVAL ────────────────────────────────────────────────────────
# Change this to control how many steps to skip.
# STEP=15 means use steps 0, 15, 30, ... (matching the usual save interval).
STEP = 15

# ── Monkey-patch the three post-processing read functions ─────────────────────

_orig_read_node_vector = post.read_node_vector
_orig_read_cell_scalar = post.read_cell_scalar
_orig_read_cell_tensor = post.read_cell_tensor


def _read_node_vector_sub(path):
    points, time_list, vector_field = _orig_read_node_vector(path)
    return points, time_list[::STEP], vector_field[::STEP]


def _read_cell_scalar_sub(path):
    centroids, time_list, scalar_field = _orig_read_cell_scalar(path)
    return centroids, time_list[::STEP], scalar_field[::STEP]


def _read_cell_tensor_sub(path):
    centroids, time_list, tensor_field = _orig_read_cell_tensor(path)
    return centroids, time_list[::STEP], tensor_field[::STEP]


post.read_node_vector = _read_node_vector_sub
post.read_cell_scalar = _read_cell_scalar_sub
post.read_cell_tensor = _read_cell_tensor_sub

# ── Now import and run the original main ──────────────────────────────────────
from plot_results import main  # noqa: E402

if __name__ == "__main__":
    main()
