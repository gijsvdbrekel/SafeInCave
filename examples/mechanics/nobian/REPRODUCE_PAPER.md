# Reproducing the paper figures

Paper: *Influence of construction-driven geometric deviations on the structural
integrity of salt caverns for underground hydrogen storage* (van den Brekel,
Honório, Hajibeygi, Bakker).

This document lists, figure by figure, which script to run and with which
settings. Everything below refers to paths relative to this folder
(`examples/mechanics/nobian/`).

---

## 0. Before you start

**SafeInCave version.** These scripts were written against **SafeInCave v2.0.0**
(the version in this repo). The paper cites **v3.0.3**
(Zenodo `10.5281/zenodo.20722991`). That is a major version bump, so expect API
changes. The scripts use:

```python
import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
import safeincave.HeatBC   as heatBC
```

and these classes: `Material`, `Spring`, `DislocationCreep`,
`PressureSolutionCreep`, `Viscoelastic`, `MunsonDawsonCreep`,
`ViscoplasticDesai`, `MohrCoulombViscoplastic`, `LinearMomentumMod`.

> **Check `MunsonDawsonCreep` first.** It was implemented as part of this work and
> merged into SafeInCave. All paper results use it (`USE_MUNSON_DAWSON = True`).
> If its constructor signature changed in v3, every run script needs updating.

**Two simulation scripts only:**

| Script | Produces |
|---|---|
| `Simulation/Run.py` | the six homogeneous cavern shapes (Figs. 3, 4, 5, 6, 7 + Table 3) |
| `Simulation/run_interlayer.py` | the two heterogeneous anhydrite cases + control (Fig. 8) |

**Pressure profile file** — `Simulation/drukprofiel_zoutcaverne_2035_8760u.csv`
must stay next to the run scripts. This is the *2035 projected demand* scheme
used for every result in the paper.

**Output folders** are created automatically under `Simulation/output/` and are
named from the settings, e.g.

```
output/case_leaching_linear_csv(99)_1825days_SA_MD_regular1200        <- Run.py
output/case_spike_upper_il4x_csv(99)_1095days_SA_MD                   <- run_interlayer.py
```

The plotting scripts scan `Simulation/output/` and filter on these names, so do
not rename them.

---

## 1. Homogeneous runs — `Simulation/Run.py`

Run this **six times**, changing only `CAVERN_TYPE`. Everything else stays fixed.

```python
CAVERN_TYPE           = "regular"   # then: directcirculation, reversedcirculation,
                                    #       fastleached, tilted, tubefailure
CAVERN_SIZE           = 1200

USE_LEACHING          = True        # <-- leaching ramp, NOT equilibrium
LEACHING_MODE         = "linear"
LEACHING_DAYS         = 91
DEBRINING_DAYS        = 30
LEACHING_END_FRACTION = 0.30

PRESSURE_SCENARIO     = "csv"
CSV_FILE_PATH         = "drukprofiel_zoutcaverne_2035_8760u.csv"
OPERATION_DAYS        = 1825        # 5 years

MATERIAL_SCENARIO     = "A"         # CCC Zuidwending field-calibrated set
USE_MUNSON_DAWSON     = True
```

Naming: `CAVERN_TYPE` + `CAVERN_SIZE` gives the *cavern key* used everywhere
downstream — `regular1200`, `directcirculation1200`, `reversedcirculation1200`,
`fastleached1200`, `tilted1200`, `tubefailure1200`.

> `tubefailure` is what the paper calls **string-failure**. Same thing.

Run from the `Simulation/` folder:

```bash
cd Simulation
python Run.py
```

All six must finish before any of Figs. 4–7 or Table 3 can be made.

---

## 2. Heterogeneous runs — `Simulation/run_interlayer.py`

Run this **three times**, changing only `CAVERN_TYPE`:

```python
CAVERN_TYPE       = "spike_upper_il4x"   # then: spike_lower_il4x, spike_none
INTERLAYER_1_MATERIAL = "anhydrite"
INTERLAYER_MODEL  = "drucker_prager"

USE_LEACHING      = False        # <-- equilibrium start, NOT the leaching ramp
PRESSURE_SCENARIO = "csv"
CSV_FILE_PATH     = "drukprofiel_zoutcaverne_2035_8760u.csv"
OPERATION_DAYS    = 1095         # ~3 years
MATERIAL_SCENARIO = "A"
USE_MUNSON_DAWSON = True
```

**Two things to be careful about here:**

1. **Use the `il4x` grids.** The paper results are on the **four-times refined**
   interlayer mesh: `spike_upper_il4x`, `spike_lower_il4x`, with `spike_none`
   as the homogeneous reference. The `il2x` grids exist only for the
   mesh-refinement study (thesis Appendix D.3) — do **not** use them for the
   paper figures. Note the file currently ships with
   `CAVERN_TYPE = "spike_upper_il2x"` as its default, and the comment block
   above that line does not list the `il2x`/`il4x` options at all, even though
   they are valid. They are defined further down in `VALID_SHAPES`,
   `CAVERN_GEOMETRY` and `GRID_FOLDERS`. Grids live in `grids/`:
   `cavern_spike_upper_il4x_1200_3D`, `cavern_spike_lower_il4x_1200_3D`,
   `cavern_spike_none_1200_3D`.

2. **`USE_LEACHING = False` here, but `True` in `Run.py`.** The heterogeneous
   cases start from a short equilibrium phase instead of the 91-day leaching
   ramp. This is deliberate and is why the paper compares the heterogeneous
   cases only against the dedicated `spike_none` reference, never against the
   `Run.py` regular cavern. Keep it as is, or the comparison is not like-for-like.

---

## 3. Figure-by-figure

All plotting is done from `Plotting/Plot_putty/`.

### Figure 1 — the six cavern geometries in the domain box
### Figure 2 — the two heterogeneous configurations

Not from simulation output — rendered straight from the grids.

```bash
cd ../../../FinalDefense           # repo root /FinalDefense
python render_domain_overview.py       # -> cavern_renders/domain_overview_comparison.png  (Fig. 1)
python render_heterogeneous_caverns.py # -> cavern_renders/heterogeneous_comparison.png    (Fig. 2)
```

Then convert to PDF and drop into the Overleaf `figures/` folder (the paper uses
`.pdf` for all figures).

### Figure 3 — cavern wall profiles + probe locations

```bash
cd Plotting/Plot_putty
python plot_cavern_profiles.py
```

Reads the grids directly, no simulation output needed. Produces the two
`cavern_profiles_group*.png` panels.

### Figure 4 — p–q stress paths

`plot_results.py`, with:

```python
SELECT = {
    "caverns": ["regular1200", "directcirculation1200", "reversedcirculation1200",
                "fastleached1200", "tilted1200", "tubefailure1200"],
    "pressure": ["csv"],
    "scenario": ["MD_A"],
    "n_cycles": None,
    "operation_days": 1825,
    "case_contains": None,
}

PLOT_MODE = "compare_shapes"       # all six shapes overlaid on one figure

FIGURES = {
    "convergence": False,
    "stress_state": True,          # <-- Figure 4
    "fos": False,
    "fracture_propagation": False,
    "fos_summary": False,
    "mc_failure": False,
    "interlayer_comparison": False,
    "interlayer_stress_paths": False,
}

SHOW_DILATANCY = ["ratigan_027", "spiers", "devries_comp", "devries_ext"]
```

```bash
python plot_results.py
```

### Figure 6 — radial extent of the dilatant zone

Same `SELECT` and `PLOT_MODE` as Figure 4. Only the `FIGURES` flags change:

```python
FIGURES = { ..., "fracture_propagation": True, ... }   # all others False
```

This produces the per-region dilatant-zone extent vs. time — the two stacked
panels that became `Dilatancyzone1/2` in the paper.

### Figure 7 — cavern convergence

Same `SELECT` and `PLOT_MODE` again:

```python
FIGURES = { ..., "convergence": True, ... }            # all others False
```

Gives convergence curves for the six shapes plus the applied pressure schedule.

### Table 3 — geometry summary

Not a single script. Assembled from `plot_results.py` output on the same six
cases:

| Table 3 column | Where it comes from |
|---|---|
| Dilatant elements at `p_min` (%) | `"fos": True` — count of cells with FoS < 1 at the minimum-pressure timestep, divided by the number of salt elements within 30 m of the wall |
| Wall regions with FoS < 1 | `"fracture_propagation": True` — which of the five z-regions show a non-zero extent |
| Max. radial extent (m) | `"fracture_propagation": True` |
| Persistent regions | `"fracture_propagation": True` — regions whose extent does **not** decay over the 5 years |
| Convergence 5 yr (%) | `"convergence": True` — final value of each curve |

The near-wall element counts used for the percentages (25 411 for regular up to
32 629 for fast-leached) come from the meshes themselves, not the simulation —
count salt tetrahedra whose centroid is within 30 m of the cavern wall surface
(physical group 29).

### Figure 5 — factor-of-safety fields on the cavern wall

Two steps: export with `paraview.py`, then render in ParaView.

```python
# paraview.py
SELECT = {
    "caverns": ["regular1200", "directcirculation1200", "reversedcirculation1200",
                "fastleached1200", "tilted1200", "tubefailure1200"],
    "pressure": ["csv"],
    "scenario": None,
    "n_cycles": None,
    "operation_days": 1825,
    "case_contains": None,
}
PHASE     = "operation"
TIMESTEPS = [0, "p_min", "p_max", "p_min_top_3", -1]
```

```bash
python paraview.py
```

This writes, per case: `cavern_surface.xdmf`, `cavern_surface_max.xdmf`,
`paraview_volume.xdmf`, `paraview_volume_max.xdmf`, `cross_section.xdmf`,
`convergence.xdmf`.

In ParaView, for each of the six cases:

1. Open `cavern_surface.xdmf` and step to the **`p_min` timestep**
   (the paper uses *p* = 6 MPa, *t* = 134 days).
2. Colour by the **FoS** field.
3. Set the colour scale to a **hard break at FoS = 1**: red/purple below 1
   (dilatant), blue above 1, range roughly 0.2 → 5. The abrupt break is
   deliberate and is called out in the figure caption.
4. Additionally load `paraview_volume.xdmf` and threshold on **FoS < 1** to show
   the dilatant cells in the salt *around* the wall — those are the red specks
   surrounding each cavern in the figure.
5. Use the same camera position and the same colour-scale limits for all six,
   otherwise the panels are not comparable.

### Figure 8 — heterogeneity, per-cell minimum FoS

Needs the three `run_interlayer.py` cases from §2. `paraview.py` is already
configured for this — it currently ships with exactly these settings:

```python
SELECT = {
    "caverns": ["spike_none", "spike_lower_il4x", "spike_upper_il4x"],
    "pressure": None,
    "scenario": None,
    "n_cycles": None,
    "operation_days": 1095,
    "case_contains": None,
}
```

```bash
python paraview.py
```

Use the `*_max.xdmf` outputs — those carry the **per-cell minimum FoS reached at
any time** during the simulation, which is what the figure shows (not a single
timestep). Render the three panels in ParaView with identical camera and colour
scale. The anhydrite interlayer is shown in grey and excluded from the FoS
evaluation.

The related interlayer plots (Mohr–Coulomb failed-volume fraction over the
cycles, and the convergence comparison) come from `plot_results.py` with
`"mc_failure": True` and `"interlayer_comparison": True`.

---

## 4. Suggested order of work

1. Get `Run.py` working on v3 for **one** shape (`regular`) and confirm the
   convergence curve matches the paper's 0.82 % at 1825 days.
2. If that matches, run the other five shapes.
3. Reproduce Figs. 4, 6, 7 and Table 3 from those six.
4. Then do `run_interlayer.py` (three cases) and Fig. 8.
5. Figs. 1, 2, 3 need no simulation output and can be done at any point.

Values to check against the paper as you go:

| Quantity | Paper value |
|---|---|
| Convergence after 1825 d, regular | 0.82 % |
| Convergence after 1825 d, fast-leached | 0.90 % (highest) |
| Convergence after 1825 d, reversed-circulation | 0.65 % (lowest) |
| Dilatant elements at `p_min`, fast-leached | 22.4 % of near-wall |
| Dilatant elements at `p_min`, reversed-circulation | 9.1 % of near-wall |
| Max radial extent of dilatant zone | 30 m (all shapes except reversed-circulation) |
| Convergence after 1100 d, heterogeneous-above / below / reference | 0.26 % / 0.25 % / 0.30 % |

If the first regular-cavern run reproduces 0.82 %, the port is almost certainly
correct and the rest should follow.

---

## 5. Other scripts in these folders

Not needed for the paper, but present:

- `Simulation/Run_pressure_swing.py`, `Run_sensitivity.py`, `MultipleCycles.py`,
  `ScenarioTest.py`, `VariableTimestep.py`, `quick_test.py` — thesis-only studies
  (pressure schemes, daily pressure rate, sensitivity analysis).
- `Plotting/Plot_putty/plot_pressure_swing2.py`,
  `plot_pressure_swing_comparison.py`, `plot_sensitivity.py`,
  `desai_effect.py`, `disloc_effect.py` — plots for those thesis-only studies.
- `Plotting/Plot_putty/case_index.py` — shared helper that scans
  `Simulation/output/` and parses the case-folder names. Both `plot_results.py`
  and `paraview.py` import it. If the output naming ever changes, this is the
  one file to update.
- `plot_results_subsample.py` — same as `plot_results.py` but subsamples large
  XDMF files; useful if memory is tight.
