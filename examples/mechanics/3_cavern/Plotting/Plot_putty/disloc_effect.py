import os
import json
import numpy as np
import matplotlib.pyplot as plt
import meshio

try:
    from case_index import detect_layout_and_collect_cases, filter_cases
except ModuleNotFoundError:
    from case_Index import detect_layout_and_collect_cases, filter_cases

DAY = 24.0 * 3600.0
MPA = 1e6

# =============================================================================
# USER CONFIG
# =============================================================================
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/3_cavern/output"

SELECT_BASE = {
    "caverns": ["Regular"],     # try also ["regular600"] if needed
    "pressure": "sinus",
    "n_cycles": 8,
    "operation_days": 365,
    "case_contains": None,
}

SCEN_OLD = "disloc_old_only"
SCEN_NEW = "disloc_new_only"

PICK_PEAK = "last"  # "first" or "last"

OUT_DIR = os.path.join(ROOT, "_figures")
SHOW = False
DPI = 180


# =============================================================================
# Helpers
# =============================================================================
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

    if "t_values" in data and "p_values" in data:
        t_s = np.asarray(data["t_values"], float)
        p_mpa = np.asarray(data["p_values"], float) / MPA
        return t_s, p_mpa

    raise RuntimeError(f"Unknown pressure JSON keys: {list(data.keys())}")

def pick_peak_timestep_index(t_s, p_mpa):
    p = np.asarray(p_mpa)
    peaks = np.where((p[1:-1] > p[:-2]) & (p[1:-1] > p[2:]))[0] + 1
    if peaks.size == 0:
        return int(np.argmax(p))
    return int(peaks[0] if PICK_PEAK == "first" else peaks[-1])

def read_cell_scalar_timeseries(xdmf_path: str):
    if not os.path.isfile(xdmf_path):
        raise RuntimeError(f"Missing XDMF: {xdmf_path}")

    reader = meshio.xdmf.TimeSeriesReader(xdmf_path)
    reader.read_points_cells()

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
        arr = np.asarray(arr_blocks[0] if isinstance(arr_blocks, list) else arr_blocks, float)
        vals.append(arr)

    return np.asarray(times, float), vals

def bin_median_trend(x, y, nbins=35):
    x = np.asarray(x); y = np.asarray(y)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if x.size < 100:
        return None, None
    edges = np.linspace(np.min(x), np.max(x), nbins + 1)
    xm, ym = [], []
    for i in range(nbins):
        mi = (x >= edges[i]) & (x < edges[i+1])
        if np.count_nonzero(mi) < 30:
            continue
        xm.append(0.5 * (edges[i] + edges[i+1]))
        ym.append(np.median(y[mi]))
    return np.asarray(xm), np.asarray(ym)

def newest_case(candidates):
    def mtime(m):
        try:
            return os.path.getmtime(m["case_path"])
        except Exception:
            return -1.0
    return sorted(candidates, key=mtime)[-1]

def pick_case_for_scenario(all_cases_meta, scenario_name: str):
    sel = dict(SELECT_BASE)
    sel["scenario"] = scenario_name
    candidates = filter_cases(all_cases_meta, sel)

    needed_fields = ["q_elems", "eps_ne_rate_eq_disloc"]
    keep = []
    for m in candidates:
        cp = m["case_path"]
        ok = os.path.isfile(path_pressure_json(cp))
        ok = ok and all(os.path.isfile(path_field_xdmf(cp, f)) for f in needed_fields)
        if ok:
            keep.append(m)

    if not keep:
        print(f"\n[ERROR] No usable case for scenario='{scenario_name}' with SELECT_BASE={SELECT_BASE}")
        print("[DEBUG] Try setting caverns to ['regular600'] or setting case_contains='regular600'")
        raise RuntimeError("Selection matched nothing. Fix SELECT_BASE or case naming.")

    chosen = newest_case(keep)
    return chosen["case_path"], chosen

def plot_old_vs_new(case_old: str, case_new: str, meta_old: dict):
    tag = f"{meta_old.get('cavern_label')}|p={meta_old.get('pressure_scenario')}|{meta_old.get('n_cycles')}cyc|{meta_old.get('operation_days')}d"
    tag = tag.replace(" ", "")

    t_sched_s, p_sched_mpa = read_pressure_schedule(case_old)
    i_peak = pick_peak_timestep_index(t_sched_s, p_sched_mpa)
    target_t = t_sched_s[i_peak]

    q_old_t, q_old_list = read_cell_scalar_timeseries(path_field_xdmf(case_old, "q_elems"))
    r_old_t, r_old_list = read_cell_scalar_timeseries(path_field_xdmf(case_old, "eps_ne_rate_eq_disloc"))

    q_new_t, q_new_list = read_cell_scalar_timeseries(path_field_xdmf(case_new, "q_elems"))
    r_new_t, r_new_list = read_cell_scalar_timeseries(path_field_xdmf(case_new, "eps_ne_rate_eq_disloc"))

    k_old = int(np.argmin(np.abs(q_old_t - target_t)))
    k_new = int(np.argmin(np.abs(q_new_t - target_t)))

    q_old = np.asarray(q_old_list[k_old], float) / MPA
    r_old = np.asarray(r_old_list[k_old], float)
    q_new = np.asarray(q_new_list[k_new], float) / MPA
    r_new = np.asarray(r_new_list[k_new], float)

    eps = 1e-30
    x_old = np.log10(np.maximum(q_old, eps))
    y_old = np.log10(np.maximum(r_old, eps))
    x_new = np.log10(np.maximum(q_new, eps))
    y_new = np.log10(np.maximum(r_new, eps))

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.set_title(f"Dislocation: log10(strain-rate) vs log10(von Mises)\npeak='{PICK_PEAK}' @ t≈{target_t/DAY:.2f} days — {tag}")
    ax.scatter(x_old, y_old, s=4, alpha=0.25, label="OLD dislocation")
    ax.scatter(x_new, y_new, s=4, alpha=0.25, label="NEW dislocation")

    xm, ym = bin_median_trend(x_old, y_old)
    if xm is not None:
        ax.plot(xm, ym, linewidth=2.0, label="OLD median trend")
    xm, ym = bin_median_trend(x_new, y_new)
    if xm is not None:
        ax.plot(xm, ym, linewidth=2.0, label="NEW median trend")

    ax.set_xlabel("log10(q_elems [MPa])")
    ax.set_ylabel("log10(epsdot_disloc_eq [1/s])")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 3))
    ax2.set_title("Pressure schedule (MPa) – selected peak")
    ax2.plot(t_sched_s / DAY, p_sched_mpa, linewidth=1.5)
    ax2.scatter([t_sched_s[i_peak] / DAY], [p_sched_mpa[i_peak]], s=60)
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Pressure (MPa)")
    ax2.grid(True, alpha=0.3)

    return fig, fig2, tag

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    all_cases = detect_layout_and_collect_cases(ROOT)
    if not all_cases:
        raise RuntimeError(f"No cases found under ROOT={ROOT}")

    case_old, meta_old = pick_case_for_scenario(all_cases, SCEN_OLD)
    case_new, _ = pick_case_for_scenario(all_cases, SCEN_NEW)

    print("OLD case:", case_old)
    print("NEW case:", case_new)

    fig, fig2, tag = plot_old_vs_new(case_old, case_new, meta_old)

    p1 = os.path.join(OUT_DIR, f"disloc_effect_scatter_{tag}_peak={PICK_PEAK}.png")
    p2 = os.path.join(OUT_DIR, f"disloc_effect_pressure_{tag}_peak={PICK_PEAK}.png")

    fig.savefig(p1, dpi=DPI, bbox_inches="tight")
    fig2.savefig(p2, dpi=DPI, bbox_inches="tight")

    print("[SAVED]", p1)
    print("[SAVED]", p2)

    if SHOW:
        plt.show()

    plt.close(fig); plt.close(fig2)

if __name__ == "__main__":
    main()
