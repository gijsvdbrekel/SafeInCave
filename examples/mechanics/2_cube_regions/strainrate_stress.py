import os
import json
import glob
import numpy as np
import matplotlib.pyplot as plt
import meshio

DAY = 24.0 * 3600.0
MPA = 1e6

# -------------------------
# CONFIG: zet deze goed
# -------------------------
ROOT_GROUP = "/home/gvandenbrekel/SafeInCave/OutputNobian/Regular_sinus_600"
CASE_OLD_KEY = "disloc_old_only"
CASE_NEW_KEY = "disloc_new_only"

# Optional: kies welke peak je wil plotten (laatste peak is meestal het meest relevant)
PICK_PEAK = "last"  # "first" or "last"

# -------------------------
# Helpers: find cases
# -------------------------
def find_case_folder(root_group: str, key: str) -> str:
    """Find a case folder under ROOT_GROUP that contains `key` in its name."""
    candidates = []
    for p in sorted(glob.glob(os.path.join(root_group, "case_*"))):
        if os.path.isdir(p) and key in os.path.basename(p).lower():
            candidates.append(p)
    if not candidates:
        raise RuntimeError(f"No case folder found containing '{key}' under {root_group}")
    # if multiple, take last (often newest run); change if you want
    return candidates[-1]

def path_pressure_json(case_folder):
    return os.path.join(case_folder, "pressure_schedule.json")

def path_field_xdmf(case_folder, field):
    return os.path.join(case_folder, "operation", field, f"{field}.xdmf")

def read_pressure_schedule(case_folder: str):
    pjson = path_pressure_json(case_folder)
    if not os.path.isfile(pjson):
        raise RuntimeError(f"Missing pressure_schedule.json in {case_folder}")

    with open(pjson, "r") as f:
        data = json.load(f)

    if "t_hours" in data and "p_MPa" in data:
        t_s = np.asarray(data["t_hours"], float) * 3600.0
        p_mpa = np.asarray(data["p_MPa"], float)
        return t_s, p_mpa

    if "t_values_s" in data and "p_values_Pa" in data:
        t_s = np.asarray(data["t_values_s"], float)
        p_mpa = np.asarray(data["p_values_Pa"], float) / MPA
        return t_s, p_mpa

    # fallback for older names (from your example)
    if "t_values" in data and "p_values" in data:
        t_s = np.asarray(data["t_values"], float)
        p_mpa = np.asarray(data["p_values"], float) / MPA
        return t_s, p_mpa

    raise RuntimeError(f"Unknown pressure JSON format keys: {list(data.keys())}")

def pick_peak_timestep_index(t_s, p_mpa):
    """Return index of a peak pressure timestep (first/last)."""
    p = np.asarray(p_mpa)
    peaks = np.where((p[1:-1] > p[:-2]) & (p[1:-1] > p[2:]))[0] + 1
    if peaks.size == 0:
        # fallback: just take global max
        return int(np.argmax(p))
    if PICK_PEAK == "first":
        return int(peaks[0])
    return int(peaks[-1])

# -------------------------
# XDMF time series reader for cell scalars
# -------------------------
def read_cell_scalar_timeseries(xdmf_path: str):
    """
    Returns:
      times: (Nt,) float seconds (if present; else 0..Nt-1)
      values: list of cell arrays (Nc,) per timestep
    """
    if not os.path.isfile(xdmf_path):
        raise RuntimeError(f"Missing XDMF: {xdmf_path}")

    reader = meshio.xdmf.TimeSeriesReader(xdmf_path)
    points, cells = reader.read_points_cells()

    times = []
    vals = []
    for k in range(reader.num_steps):
        t, point_data, cell_data = reader.read_data(k)
        times.append(float(t) if t is not None else float(k))

        # cell_data is dict: name -> list per cell block
        # We expect exactly one named field, same as folder name.
        # But some writers store it under a generic key; we try robust.
        if not cell_data:
            raise RuntimeError(f"No cell_data in timestep {k} for {xdmf_path}")

        # If multiple keys exist, prefer the one matching basename
        base = os.path.splitext(os.path.basename(xdmf_path))[0].lower()
        key = None
        for cand in cell_data.keys():
            if str(cand).lower() == base:
                key = cand
                break
        if key is None:
            key = list(cell_data.keys())[0]

        # cell_data[key] is list over cell blocks; take first block
        arr_blocks = cell_data[key]
        if isinstance(arr_blocks, list):
            arr = np.asarray(arr_blocks[0], float)
        else:
            arr = np.asarray(arr_blocks, float)

        vals.append(arr)

    return np.asarray(times, float), vals

def bin_median_trend(x, y, nbins=30):
    """Compute median y in bins of x for a trend line."""
    x = np.asarray(x)
    y = np.asarray(y)
    mask = np.isfinite(x) & np.isfinite(y) & (x > -np.inf) & (y > -np.inf)
    x = x[mask]; y = y[mask]
    if x.size < 10:
        return None, None

    edges = np.linspace(np.nanmin(x), np.nanmax(x), nbins + 1)
    xm = []
    ym = []
    for i in range(nbins):
        m = (x >= edges[i]) & (x < edges[i+1])
        if np.count_nonzero(m) < 30:
            continue
        xm.append(0.5 * (edges[i] + edges[i+1]))
        ym.append(np.nanmedian(y[m]))
    return np.asarray(xm), np.asarray(ym)

# -------------------------
# Main plotting function
# -------------------------
def plot_old_vs_new(case_old: str, case_new: str):
    # pressure peak selection based on schedule (assume both have same schedule)
    t_sched_s, p_sched_mpa = read_pressure_schedule(case_old)
    i_peak = pick_peak_timestep_index(t_sched_s, p_sched_mpa)

    # read time series fields
    q_old_t, q_old_list = read_cell_scalar_timeseries(path_field_xdmf(case_old, "q_elems"))
    r_old_t, r_old_list = read_cell_scalar_timeseries(path_field_xdmf(case_old, "eps_ne_rate_eq_disloc"))

    q_new_t, q_new_list = read_cell_scalar_timeseries(path_field_xdmf(case_new, "q_elems"))
    r_new_t, r_new_list = read_cell_scalar_timeseries(path_field_xdmf(case_new, "eps_ne_rate_eq_disloc"))

    # choose timestep k: map schedule peak to nearest available output time
    # (sometimes output is sparse; so nearest match)
    target_t = t_sched_s[i_peak]
    k_old = int(np.argmin(np.abs(q_old_t - target_t)))
    k_new = int(np.argmin(np.abs(q_new_t - target_t)))

    q_old = np.asarray(q_old_list[k_old], float) / MPA
    r_old = np.asarray(r_old_list[k_old], float)

    q_new = np.asarray(q_new_list[k_new], float) / MPA
    r_new = np.asarray(r_new_list[k_new], float)

    # avoid zeros/negatives for log10
    eps = 1e-30
    x_old = np.log10(np.maximum(q_old, eps))
    y_old = np.log10(np.maximum(r_old, eps))

    x_new = np.log10(np.maximum(q_new, eps))
    y_new = np.log10(np.maximum(r_new, eps))

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.set_title("Dislocation creep: log10(strain-rate) vs log10(von Mises stress)\n(peak pressure timestep)")
    ax.scatter(x_old, y_old, s=4, alpha=0.25, label="OLD (CCC) dislocation")
    ax.scatter(x_new, y_new, s=4, alpha=0.25, label="NEW (Herminio) dislocation")

    # add median trend
    xm, ym = bin_median_trend(x_old, y_old, nbins=35)
    if xm is not None:
        ax.plot(xm, ym, linewidth=2.0, label="OLD median trend")
    xm, ym = bin_median_trend(x_new, y_new, nbins=35)
    if xm is not None:
        ax.plot(xm, ym, linewidth=2.0, label="NEW median trend")

    ax.set_xlabel("log10(q_elems [MPa])")
    ax.set_ylabel("log10(epsdot_disloc_eq [1/s])")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    # also show pressure schedule for context
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 3))
    ax2.set_title("Pressure schedule (MPa) – with selected peak")
    ax2.plot(t_sched_s / DAY, p_sched_mpa, linewidth=1.5)
    ax2.scatter([t_sched_s[i_peak] / DAY], [p_sched_mpa[i_peak]], s=50)
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Pressure (MPa)")
    ax2.grid(True, alpha=0.3)

    return fig, fig2

def main():
    case_old = find_case_folder(ROOT_GROUP, CASE_OLD_KEY)
    case_new = find_case_folder(ROOT_GROUP, CASE_NEW_KEY)

    print("OLD case:", case_old)
    print("NEW case:", case_new)

    fig, fig2 = plot_old_vs_new(case_old, case_new)
    plt.show()

if __name__ == "__main__":
    main()
