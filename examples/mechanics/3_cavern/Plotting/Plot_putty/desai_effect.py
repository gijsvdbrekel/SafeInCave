import os
import json
import numpy as np
import matplotlib.pyplot as plt
import meshio

from case_index import detect_layout_and_collect_cases, filter_cases

DAY = 24.0 * 3600.0
MPA = 1e6

# =============================================================================
# USER CONFIG
# =============================================================================
ROOT = r"/data/home/gbrekel/SafeInCave_new/examples/mechanics/3_cavern/output"

SELECT_BASE = {
    "caverns": ["Regular"],   # or None
    "pressure": "sinus",      # or None
    "n_cycles": 8,            # or None
    "operation_days": 365,    # or None
    "case_contains": "regular600",  # optional extra substring, or None
}

SCEN_FULL = "full"
SCEN_NODESAI = "full_minus_desai"

OUT_DIR = os.path.join(ROOT, "_figures")
SHOW = False
DPI = 180


# =============================================================================
# Helpers
# =============================================================================
def path_field_xdmf(case_folder, field):
    return os.path.join(case_folder, "operation", field, f"{field}.xdmf")

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
        if isinstance(arr_blocks, list):
            arr = np.asarray(arr_blocks[0], float)
        else:
            arr = np.asarray(arr_blocks, float)

        vals.append(arr)

    return np.asarray(times, float), vals

def series_mean(times, vals_list):
    y = np.array([np.nanmean(v) for v in vals_list], dtype=float)
    return np.asarray(times, float), y

def read_ksp_jsonl(case_folder: str):
    p = os.path.join(case_folder, "ksp_operation.jsonl")
    if not os.path.isfile(p):
        raise RuntimeError(f"Missing {p}")

    t, its, rnorm, reason = [], [], [], []
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

def pick_case_for_scenario(all_cases_meta, scenario_name: str):
    sel = dict(SELECT_BASE)
    sel["scenario"] = scenario_name
    candidates = filter_cases(all_cases_meta, sel)

    # require files used in this comparison
    needed_fields = ["q_elems", "eps_ne_rate_eq_total", "eps_ne_rate_eq_desai"]
    keep = []
    for m in candidates:
        cp = m["case_path"]
        ok = os.path.isfile(os.path.join(cp, "ksp_operation.jsonl"))
        ok = ok and all(os.path.isfile(path_field_xdmf(cp, f)) for f in needed_fields)
        if ok:
            keep.append(m)

    if not keep:
        # debug print to show why filtering fails
        print(f"\n[DEBUG] No usable case for scenario='{scenario_name}' with SELECT_BASE={SELECT_BASE}")
        print("[DEBUG] First 20 found cases (label | scen | press | name):")
        for m in all_cases_meta[:20]:
            print(" -", m.get("cavern_label"), m.get("scenario_preset"), m.get("pressure_scenario"), m.get("case_name"))
        raise RuntimeError(f"No usable case found for scenario='{scenario_name}'")

    # pick newest-ish: last in list is usually newest from collector ordering
    return keep[-1]["case_path"], keep[-1]

def plot_desai_effect(case_full: str, case_nodesai: str, meta_full: dict, meta_nodesai: dict):
    tag = f"{meta_full.get('cavern_label')}|p={meta_full.get('pressure_scenario')}|{meta_full.get('n_cycles')}cyc|{meta_full.get('operation_days')}d"
    tag = tag.replace(" ", "")

    # --- KSP iterations vs time ---
    tF, itsF, rF, reasonF = read_ksp_jsonl(case_full)
    tN, itsN, rN, reasonN = read_ksp_jsonl(case_nodesai)

    fig1, ax = plt.subplots(1, 1, figsize=(11, 4))
    ax.set_title(f"KSP iterations vs time (operation) — {tag}")
    ax.plot(tF / DAY, itsF, linewidth=1.8, label="FULL (with Desai)")
    ax.plot(tN / DAY, itsN, linewidth=1.8, label="FULL_MINUS_DESAI")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("KSP iterations")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    # --- mean q_elems vs time ---
    tqF, qF_list = read_cell_scalar_timeseries(path_field_xdmf(case_full, "q_elems"))
    tqN, qN_list = read_cell_scalar_timeseries(path_field_xdmf(case_nodesai, "q_elems"))
    tqF, qF = series_mean(tqF, qF_list)
    tqN, qN = series_mean(tqN, qN_list)

    fig2, ax2 = plt.subplots(1, 1, figsize=(11, 4))
    ax2.set_title(f"Mean q_elems vs time — {tag}")
    ax2.plot(tqF / DAY, qF / MPA, linewidth=1.8, label="FULL (with Desai)")
    ax2.plot(tqN / DAY, qN / MPA, linewidth=1.8, label="FULL_MINUS_DESAI")
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Mean von Mises (MPa)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="best")

    # --- strain-rate decomposition ---
    ttF, rtotF_list = read_cell_scalar_timeseries(path_field_xdmf(case_full, "eps_ne_rate_eq_total"))
    tdF, rdesF_list = read_cell_scalar_timeseries(path_field_xdmf(case_full, "eps_ne_rate_eq_desai"))

    ttN, rtotN_list = read_cell_scalar_timeseries(path_field_xdmf(case_nodesai, "eps_ne_rate_eq_total"))
    tdN, rdesN_list = read_cell_scalar_timeseries(path_field_xdmf(case_nodesai, "eps_ne_rate_eq_desai"))

    ttF, rtotF = series_mean(ttF, rtotF_list)
    tdF, rdesF = series_mean(tdF, rdesF_list)
    ttN, rtotN = series_mean(ttN, rtotN_list)
    tdN, rdesN = series_mean(tdN, rdesN_list)

    fig3, ax3 = plt.subplots(1, 1, figsize=(11, 4))
    ax3.set_title(f"Mean eq strain-rates (log scale) — {tag}")
    ax3.semilogy(ttF / DAY, np.maximum(rtotF, 1e-30), linewidth=1.8, label="FULL total")
    ax3.semilogy(tdF / DAY, np.maximum(rdesF, 1e-30), linewidth=1.8, label="FULL Desai component")
    ax3.semilogy(ttN / DAY, np.maximum(rtotN, 1e-30), linewidth=1.8, label="NO_DESAI total")
    ax3.semilogy(tdN / DAY, np.maximum(rdesN, 1e-30), linewidth=1.8, label="NO_DESAI Desai component")
    ax3.set_xlabel("Time (days)")
    ax3.set_ylabel("Mean eq strain-rate (1/s)")
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc="best")

    return fig1, fig2, fig3, tag

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    all_cases = detect_layout_and_collect_cases(ROOT)

    case_full, meta_full = pick_case_for_scenario(all_cases, SCEN_FULL)
    case_nodesai, meta_nodesai = pick_case_for_scenario(all_cases, SCEN_NODESAI)

    print("FULL case:", case_full)
    print("NO_DESAI case:", case_nodesai)

    fig1, fig2, fig3, tag = plot_desai_effect(case_full, case_nodesai, meta_full, meta_nodesai)

    p1 = os.path.join(OUT_DIR, f"desai_effect_ksp_{tag}.png")
    p2 = os.path.join(OUT_DIR, f"desai_effect_qmean_{tag}.png")
    p3 = os.path.join(OUT_DIR, f"desai_effect_rates_{tag}.png")

    fig1.savefig(p1, dpi=DPI, bbox_inches="tight")
    fig2.savefig(p2, dpi=DPI, bbox_inches="tight")
    fig3.savefig(p3, dpi=DPI, bbox_inches="tight")

    print("[SAVED]", p1)
    print("[SAVED]", p2)
    print("[SAVED]", p3)

    if SHOW:
        plt.show()

    plt.close(fig1)
    plt.close(fig2)
    plt.close(fig3)

if __name__ == "__main__":
    main()
