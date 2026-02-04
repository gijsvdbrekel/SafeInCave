import os
import json
import re
import numpy as np
import meshio

DAY = 24.0 * 3600.0
MPA = 1e6

KNOWN_PRESSURES = {"sinus", "irregular", "csv_profile", "linear"}
KNOWN_SCENARIOS = {"desai_only", "disloc_old_only", "disloc_new_only", "full", "full_minus_desai", "full_md", "md_only", "md_steady_only"}

# matches regular600, irregular1200, tilted600, etc.
_CAVERN_RE = re.compile(r"(regular|irregularfine|irregular|tilted|tilt|teardrop|asymmetric|multichamber)(600|1200)", re.I)

def _safe_read_json(path: str):
    try:
        with open(path, "r") as f:
            return json.load(f)
    except Exception:
        return None

def _infer_case_name_from_path(case_path: str) -> str:
    return os.path.basename(case_path.rstrip("/"))

def _infer_cavern_type_from_case_name(case_name: str) -> str | None:
    m = _CAVERN_RE.search(case_name)
    if not m:
        return None
    head = m.group(1).lower()
    # normalize tilt/tilted
    if head == "tilt":
        head = "tilted"
    return f"{head}{m.group(2)}".lower()

def _nice_cavern_label(cavern_type_or_group: str) -> str:
    low = (cavern_type_or_group or "").lower()
    if low.startswith("asymmetric"):
        return "Asymmetric"
    if low.startswith("multichamber"):
        return "Multichamber"
    if low.startswith("teardrop"):
        return "Teardrop"
    if low.startswith("tilt") or low.startswith("tilted"):
        return "Tilt"
    if low.startswith("regular"):
        return "Regular"
    if low.startswith("irregularfine") or low.startswith("irregular_fine"):
        return "IrregularFine"
    if low.startswith("irregular"):
        return "Irregular"
    return (cavern_type_or_group or "").split("_")[0] if cavern_type_or_group else "Unknown"

def read_case_metadata(case_path: str) -> dict:
    case_name = _infer_case_name_from_path(case_path)
    group = os.path.basename(os.path.dirname(case_path))  # only meaningful in nested layout

    cavern_type = _infer_cavern_type_from_case_name(case_name)
    cavern_label = _nice_cavern_label(cavern_type or group)

    meta = {
        "case_path": case_path,
        "case_name": case_name,
        "group": group,
        "scenario_preset": None,
        "pressure_scenario": None,
        "n_cycles": None,
        "operation_days": None,
        "cavern_type": cavern_type,
        "cavern_label": cavern_label,
        "pressure_json_path": os.path.join(case_path, "pressure_schedule.json"),
    }

    data = _safe_read_json(meta["pressure_json_path"])
    if isinstance(data, dict):
        if "pressure_scenario" in data:
            v = data.get("pressure_scenario")
            meta["pressure_scenario"] = str(v).lower() if v is not None else None

        v = data.get("scenario", None)
        if v is not None:
            vlow = str(v).lower()
            if vlow in KNOWN_SCENARIOS:
                meta["scenario_preset"] = vlow
            elif vlow in KNOWN_PRESSURES and meta["pressure_scenario"] is None:
                meta["pressure_scenario"] = vlow

        meta["n_cycles"] = data.get("n_cycles", meta["n_cycles"])
        meta["operation_days"] = data.get("operation_days", meta["operation_days"])

    # fallback parse from folder name
    name = case_name.lower()

    if meta["pressure_scenario"] is None:
        for p in ("sinus", "irregular", "csv", "linear"):
            if p in name:
                meta["pressure_scenario"] = "csv_profile" if p == "csv" else p
                break

    if meta["scenario_preset"] is None:
        # Sort by length descending to match longer names first (e.g., "full_md" before "full")
        for s in sorted(KNOWN_SCENARIOS, key=len, reverse=True):
            if f"_{s}_" in name or name.startswith(f"case_{s}_"):
                meta["scenario_preset"] = s
                break

    if meta["n_cycles"] is None:
        m = re.search(r"(\d+)\s*cyc", name)
        if m:
            meta["n_cycles"] = int(m.group(1))

    if meta["operation_days"] is None:
        m = re.search(r"(\d+)\s*days", name)
        if m:
            meta["operation_days"] = int(m.group(1))

    return meta

def detect_layout_and_collect_cases(ROOT: str) -> list[dict]:
    ROOT = os.path.abspath(ROOT)
    if not os.path.isdir(ROOT):
        raise FileNotFoundError(f"ROOT not found: {ROOT}")

    items = sorted(os.listdir(ROOT))
    flat_cases = [x for x in items if x.lower().startswith("case_") and os.path.isdir(os.path.join(ROOT, x))]
    out = []

    if flat_cases:
        for c in flat_cases:
            out.append(read_case_metadata(os.path.join(ROOT, c)))
        return out

    for group in items:
        gpath = os.path.join(ROOT, group)
        if not os.path.isdir(gpath):
            continue
        if group.lower().startswith("pressure_"):
            continue
        if group.lower().startswith("_fos_outputs"):
            continue
        for sub in sorted(os.listdir(gpath)):
            if not sub.lower().startswith("case_"):
                continue
            cpath = os.path.join(gpath, sub)
            if os.path.isdir(cpath):
                out.append(read_case_metadata(cpath))
    return out

def _match(val, want):
    """Return True if val matches want.
       want may be: None, a string, or a list/tuple/set of strings."""
    if want is None:
        return True
    val = (val or "").lower()

    # If list/tuple/set → membership check
    if isinstance(want, (list, tuple, set)):
        want = {str(w).lower() for w in want}
        return val in want

    # If single string → equality check
    return val == str(want).lower()

def filter_cases(cases: list[dict], SELECT: dict) -> list[dict]:
    want_pressure = SELECT.get("pressure", None)
    want_scenario = SELECT.get("scenario", None)
    want_caverns  = SELECT.get("caverns", None)
    want_ncyc     = SELECT.get("n_cycles", None)
    want_days     = SELECT.get("operation_days", None)
    contains      = SELECT.get("case_contains", None) or SELECT.get("case_name_contains", None)

    def ok(m: dict) -> bool:
        # pressure
        if not _match(m.get("pressure_scenario"), want_pressure):
            return False

        # scenario
        if not _match(m.get("scenario_preset"), want_scenario):
            return False

        # caverns (special: match label or cavern_type or group)
        if want_caverns is not None:
            W = {str(x).lower() for x in (want_caverns if isinstance(want_caverns,(list,tuple,set)) else [want_caverns])}
            if (m.get("cavern_label","").lower() not in W
                and m.get("cavern_type","").lower() not in W
                and m.get("group","").lower() not in W):
                return False

        # numeric filters
        if want_ncyc is not None and m.get("n_cycles") != want_ncyc:
            return False
        if want_days is not None and m.get("operation_days") != want_days:
            return False

        # substring match
        if contains is not None and str(contains).lower() not in (m.get("case_name") or "").lower():
            return False

        return True

    return [m for m in cases if ok(m)]

def one_case_per_cavern_label(cases: list[dict]) -> list[dict]:
    out = {}
    for c in cases:
        lab = c.get("cavern_label","Unknown")
        if lab not in out:
            out[lab] = c
    return list(out.values())


# =============================================================================
# SHARED PATH HELPERS
# =============================================================================

def path_field_xdmf(case_folder: str, field: str) -> str:
    """Return path to a field's XDMF file in the operation folder."""
    return os.path.join(case_folder, "operation", field, f"{field}.xdmf")


def path_pressure_json(case_folder: str) -> str:
    """Return path to pressure_schedule.json for a case."""
    return os.path.join(case_folder, "pressure_schedule.json")


# =============================================================================
# SHARED DATA READERS
# =============================================================================

def read_pressure_schedule(case_folder: str):
    """
    Read pressure schedule from JSON file.

    Returns
    -------
    t_s : np.ndarray
        Time values in seconds
    p_mpa : np.ndarray
        Pressure values in MPa
    """
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

def read_cell_scalar_timeseries(xdmf_path: str, timesteps: list[int] | None = None):
    """
    Read cell scalar data from XDMF time series.

    Parameters
    ----------
    xdmf_path : str
        Path to XDMF file
    timesteps : list[int] | None
        If provided, only read these specific timestep indices.
        If None, read all timesteps.

    Returns
    -------
    times : np.ndarray
        Time values for each timestep
    vals : list[np.ndarray]
        Cell values for each timestep
    """
    if not os.path.isfile(xdmf_path):
        raise RuntimeError(f"Missing XDMF: {xdmf_path}")

    reader = meshio.xdmf.TimeSeriesReader(xdmf_path)
    reader.read_points_cells()

    # Determine which timesteps to read
    all_steps = range(reader.num_steps)
    if timesteps is not None:
        steps_to_read = [k for k in timesteps if 0 <= k < reader.num_steps]
    else:
        steps_to_read = list(all_steps)

    times = []
    vals = []
    for k in steps_to_read:
        t, point_data, cell_data = reader.read_data(k)
        times.append(float(t) if t is not None else float(k))

        if not cell_data:
            raise RuntimeError(f"No cell_data in timestep {k} for {xdmf_path}")

        # Find matching key by filename
        base = os.path.splitext(os.path.basename(xdmf_path))[0].lower()
        key = None
        for cand in cell_data.keys():
            if str(cand).lower() == base:
                key = cand
                break
        if key is None:
            key = list(cell_data.keys())[0]

        arr_blocks = cell_data[key]
        # Handle different meshio formats: list, dict, or direct array
        if isinstance(arr_blocks, list):
            arr = np.asarray(arr_blocks[0], float)
        elif isinstance(arr_blocks, dict):
            # meshio sometimes returns {'tetra': array, ...} - take first value
            arr = np.asarray(list(arr_blocks.values())[0], float)
        else:
            arr = np.asarray(arr_blocks, float)
        vals.append(arr)

    return np.asarray(times, float), vals

def read_ksp_jsonl(case_folder: str):
    """
    Read KSP solver statistics from ksp_operation.jsonl.

    Returns
    -------
    t : np.ndarray
        Time values
    its : np.ndarray
        KSP iteration counts
    rnorm : np.ndarray
        Residual norms
    reason : np.ndarray
        Convergence reason codes
    """
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


# =============================================================================
# SHARED CASE SELECTION HELPERS
# =============================================================================

def newest_case(candidates: list[dict]) -> dict:
    """Pick the newest case by filesystem modification time."""
    def mtime(m):
        try:
            return os.path.getmtime(m["case_path"])
        except Exception:
            return -1.0
    return sorted(candidates, key=mtime)[-1]

def make_case_tag(meta: dict) -> str:
    """Create a standardized tag string from case metadata."""
    tag = f"{meta.get('cavern_label')}|p={meta.get('pressure_scenario')}|{meta.get('n_cycles')}cyc|{meta.get('operation_days')}d"
    return tag.replace(" ", "")

def debug_print_inventory(all_cases: list[dict], max_lines: int = 40):
    """Print inventory sample for debugging case selection issues."""
    print("\n[DEBUG] Inventory sample:")
    for m in all_cases[:max_lines]:
        print(f" - label={m.get('cavern_label')}, type={m.get('cavern_type')}, "
              f"press={m.get('pressure_scenario')}, scen={m.get('scenario_preset')}, name={m.get('case_name')}")

