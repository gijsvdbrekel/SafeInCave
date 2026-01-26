import os
import json
import glob
import numpy as np
import matplotlib.pyplot as plt
import meshio

DAY = 24.0 * 3600.0
MPA = 1e6

# -------------------------
# CONFIG
# -------------------------
ROOT_GROUP = "/home/gvandenbrekel/SafeInCave/OutputNobian/Regular_sinus_600"
CASE_FULL_KEY = "full"
CASE_NODESAI_KEY = "full_minus_desai"

# If your ksp log was written only on output steps, then for strict convergence comparison:
# set your output interval=1 in those runs.
# -------------------------

def find_case_folder(root_group: str, key: str) -> str:
    candidates = []
    for p in sorted(glob.glob(os.path.join(root_group, "case_*"))):
        if os.path.isdir(p) and key in os.path.basename(p).lower():
            candidates.append(p)
    if not candidates:
        raise RuntimeError(f"No case folder found containing '{key}' under {root_group}")
    return candidates[-1]

def path_field_xdmf(case_folder, field):
    return os.path.join(case_folder, "operation", field, f"{field}.xdmf")

def read_cell_scalar_timeseries(xdmf_path: str):
    if not os.path.isfile(xdmf_path):
        raise RuntimeError(f"Missing XDMF: {xdmf_path}")

    reader = meshio.xdmf.TimeSeriesReader(xdmf_path)
    points, cells = reader.read_points_cells()

    times = []
    vals = []
    for k in range(reader.num_steps):
        t, point_data, cell_data = reader.read_data(k)
        times.append(float(t) if t is not None else float(k))

        if not cell_data:
            raise RuntimeError(f"No cell_data in timestep {k} for {xdmf_path}")

        base = os.path.splitext(os.path.basename(xdmf_path))[0].lower()
        key = None
        for cand in cell_data.keys():
            if str(cand).lower() == base:
                key = cand
                break
        if key is None:
            key = list(cell_data.keys())[0]

        arr_blocks = cell_data[key]
        if isinstance(arr_blocks, list):
            arr = np.asarray(arr_blocks[0], float)
        else:
            arr = np.asarray(arr_blocks, float)

        vals.append(arr)

    return np.asarray(times, float), vals

def series_mean(times, vals_list):
    """Turn list of (Nc,) arrays into mean time series."""
    y = np.array([np.nanmean(v) for v in vals_list], dtype=float)
    return np.asarray(times, float), y

def read_ksp_jsonl(case_folder: str):
    p = os.path.join(case_folder, "ksp_operation.jsonl")
    if not os.path.isfile(p):
        raise RuntimeError(f"Missing {p}")

    t = []
    its = []
    rnorm = []
    reason = []
    with open(p, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            d = json.loads(line)
            t.append(float(d.get("t", np.nan)))
            its.append(d.get("ksp_its", None))
            rnorm.append(d.get("ksp_rnorm", None))
            reason.append(d.get("ksp_reason", None))

    t = np.asarray(t, float)
    its = np.asarray([np.nan if v is None else float(v) for v in its], float)
    rnorm = np.asarray([np.nan if v is None else float(v) for v in rnorm], float)
    reason = np.asarray([np.nan if v is None else float(v) for v in reason], float)
    return t, its, rnorm, reason

def plot_desai_effect(case_full: str, case_nodesai: str):
    # --- KSP convergence plot ---
    tF, itsF, rF, reasonF = read_ksp_jsonl(case_full)
    tN, itsN, rN, reasonN = read_ksp_jsonl(case_nodesai)

    fig1, ax = plt.subplots(1, 1, figsize=(11, 4))
    ax.set_title("Convergence: KSP iterations vs time (operation)")
    ax.plot(tF / DAY, itsF, linewidth=1.5, label="FULL (with Desai)")
    ax.plot(tN / DAY, itsN, linewidth=1.5, label="FULL_MINUS_DESAI")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("KSP iterations")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    # --- Stress state over time (mean q_elems) ---
    tqF, qF_list = read_cell_scalar_timeseries(path_field_xdmf(case_full, "q_elems"))
    tqN, qN_list = read_cell_scalar_timeseries(path_field_xdmf(case_nodesai, "q_elems"))
    tqF, qF = series_mean(tqF, qF_list)
    tqN, qN = series_mean(tqN, qN_list)

    fig2, ax2 = plt.subplots(1, 1, figsize=(11, 4))
    ax2.set_title("Stress state: mean q_elems vs time")
    ax2.plot(tqF / DAY, qF / MPA, linewidth=1.5, label="FULL (with Desai)")
    ax2.plot(tqN / DAY, qN / MPA, linewidth=1.5, label="FULL_MINUS_DESAI")
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Mean von Mises (MPa)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="best")

    # --- Strain-rate decomposition (total vs desai component) ---
    ttF, rtotF_list = read_cell_scalar_timeseries(path_field_xdmf(case_full, "eps_ne_rate_eq_total"))
    tdF, rdesF_list = read_cell_scalar_timeseries(path_field_xdmf(case_full, "eps_ne_rate_eq_desai"))

    ttN, rtotN_list = read_cell_scalar_timeseries(path_field_xdmf(case_nodesai, "eps_ne_rate_eq_total"))
    # nodesai case may still have the field written but ~0
    tdN, rdesN_list = read_cell_scalar_timeseries(path_field_xdmf(case_nodesai, "eps_ne_rate_eq_desai"))

    ttF, rtotF = series_mean(ttF, rtotF_list)
    tdF, rdesF = series_mean(tdF, rdesF_list)
    ttN, rtotN = series_mean(ttN, rtotN_list)
    tdN, rdesN = series_mean(tdN, rdesN_list)

    fig3, ax3 = plt.subplots(1, 1, figsize=(11, 4))
    ax3.set_title("Transient influence: mean equivalent strain-rates (log scale)")
    ax3.semilogy(ttF / DAY, np.maximum(rtotF, 1e-30), linewidth=1.5, label="FULL total")
    ax3.semilogy(tdF / DAY, np.maximum(rdesF, 1e-30), linewidth=1.5, label="FULL Desai component")
    ax3.semilogy(ttN / DAY, np.maximum(rtotN, 1e-30), linewidth=1.5, label="NO_DESAI total")
    ax3.semilogy(tdN / DAY, np.maximum(rdesN, 1e-30), linewidth=1.5, label="NO_DESAI Desai component")
    ax3.set_xlabel("Time (days)")
    ax3.set_ylabel("Mean eq strain-rate (1/s)")
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc="best")

    return fig1, fig2, fig3

def main():
    case_full = find_case_folder(ROOT_GROUP, CASE_FULL_KEY)
    case_nodesai = find_case_folder(ROOT_GROUP, CASE_NODESAI_KEY)

    print("FULL case:", case_full)
    print("NO_DESAI case:", case_nodesai)

    plot_desai_effect(case_full, case_nodesai)
    plt.show()

if __name__ == "__main__":
    main()

