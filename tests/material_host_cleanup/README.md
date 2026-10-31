# Material host boundary regression checks

Requires Python 3 and gfortran; builds and runs in temporary directories.

```sh
python3 tests/material_host_cleanup/check_soil_column.py
python3 tests/material_host_cleanup/check_diagnostics.py
python3 tests/material_host_cleanup/check_snowmelt_carry.py
python3 tests/material_host_cleanup/check_snow_remap.py
python3 tests/material_host_cleanup/check_flood_input.py
python3 tests/material_host_cleanup/check_flood_evap.py
python3 tests/material_host_cleanup/check_phase_bookkeeping.py
```

- `check_soil_column.py` compiles the real soil-water module, hydraulic functions,
  physical constants and location helper. Only the model namelist and unused MPI
  statistics module are stubbed. Six columns, each run for three steps, exercise
  an unsaturated column, ponded/shallow groundwater, an impermeable middle layer,
  fully impermeable soil, hydraulic root return and higher root demand. Checks
  cover finite outputs, water balance, root/runoff pathway budgets and preservation
  of snow/ice diagnostics while the soil routine resets its own fields.
- `check_diagnostics.py` compiles the production diagnostic type and checks full
  and partial reset, resizing, snow lower bounds including an empty extent,
  and object isolation.
- `check_snow_remap.py --baseline PATH` compiles all four actual standard/SNICAR
  combine/divide routines. The frozen pre-phase-3 host is the optional baseline.
  It compares all water/material/eight-aerosol outputs, including inactive-layer
  sentinels, for 18,720 calls: 1–5 layers, 13 threshold/disappearance scenarios,
  all 16/8 optional-presence combinations, two lower bounds, 0/3 species and
  repeated calls. Complete-store cases also assert material conservation.
- `check_snowmelt_carry.py` compiles the actual new helper extracted from its
  module. Six cases cover thin-snow partitioning, conservation, unrelated material
  providers, independent patches and no registered materials.
- `check_phase_bookkeeping.py` compiles the real lightweight LandWork module
  and extracted interception, water-flux and finish-step helpers with recording
  services. It checks 384 CROP/flood/snow/patch scenarios, argument values, call
  order and scratch lifetime, plus five CaMa capability-guard paths. Optional
  `--baseline PATH` compares callback traces to the frozen pre-phase-2 MAIN
  blocks; the baseline is the phase-1 working tree, not Git HEAD.
- `check_flood_evap.py` checks seven evaporation ledger/guard/unit-conversion
  scenarios. Optional `--baseline PATH` takes the pre-phase-2 CoLMMAIN file.

For a same-compiler bitwise comparison against the pre-refactor interface, place
`git show e766527a:main/HYDRO/MOD_Hydro_SoilWater.F90` in a temporary directory as
`MOD_Hydro_SoilWater.F90`, then pass that directory with `--baseline PATH` to
`check_soil_column.py`. It compares all 684 emitted host/diagnostic real values
as raw bytes. For a later baseline using the diagnostic-object interface, include its frozen
`MOD_MT_HydroDiagnostics.F90` alongside the soil-water module. Ordinary regression
execution does not require Git history.

`check_land_boundary.py --baseline-main PATH` compares the old CoLMMAIN calls
with the current CoLMMAIN → `mt_land_water_step` calls. Supply a file containing
`git show e766527a:main/CoLMMAIN.F90`. It compiles the current wrapper, production
work/diagnostic types and recording stubs with the actual physics interfaces.
All scalar/array values and optional-argument presence must match bitwise across
64 combinations: four snow extents, soil/static-wetland/dynamic-wetland routing,
split snow modes, and flood support on/off. Inputs carry distinct values per
variable and per layer to detect wrong argument forwarding and shifted slices.
The stubs record arguments; they do not run material transport physics.

`check_flood_input.py` exercises the production flood-input helper, including
credit subtraction, dissolved limits, disabled/empty input paths, and rejection
of infiltration without sufficient published credit. Its optional baseline
comparison is separate from the ordinary test run.

These focused checks do not execute the full CoLMMAIN coupling. The separate
`run_coupled.py` runner builds GRID + PFT + BGC + CROP + MaterialTransport in a
private T7920 source snapshot and runs the NEC IsoGSM fixture continuously for
two days and as two one-day segments. It requires declared read-only CASE, RAW,
RT and ISO inputs through run-7920, and writes only private job outputs.
`--baseline` restores the pre-refactor files from `validation-baseline/` after
checking their SHA256 manifest; this directory is prepared from Git HEAD when
the code-only validation snapshot is assembled. `--self-test` checks namelist
replacement without model inputs or Python scientific packages.
An optional declared BUILD input supplies a read-only previous private build.
Unchanged source is verified by content, changed/new sources are marked for
rebuilding, and Make rebuilds their dependencies in the new private snapshot.
Both entry points are always rebuilt to honor the selected runtime flags.

`compare_outputs.py LEFT RIGHT` requires identical nonempty NetCDF file sets,
groups, dimensions, variables, raw data and masking/scaling/time semantics.
It returns nonzero for any difference. Run `check_compare_outputs.py` for its
23 small fixture checks (requires NumPy and netCDF4).

`verify_coupled.py` takes declared BASELINE and CURRENT result directories and
BASELINE_LOG/CURRENT_LOG terminal runner logs. Run it with the server's
scientific Python (`/home/zhwei/miniconda3/bin/python`). It requires successful
segment-exit markers, checks every valid floating-point value is finite,
checks the expected two daily records and material/sediment fields, then
compares both revisions and both continuous/split histories and restart states.
Numerical output files stay on T7920; its result directory contains comparison
summaries and copied model/build logs for retrieval. Job IDs and outcomes are
recorded in local validation artifacts and the pull-request validation summary.

The model runner itself uses only Python's standard library; all numerical
output comparisons are performed by the separate verifier.


## Comprehensive migration batch (phase 4)

Production migration is completed as one batch before tests. Freeze the complete
previous working tree first; do not substitute Git HEAD when earlier refactors
are uncommitted. This session retains it under
`.omx/artifacts/material-host-cleanup-phase4-20261031/baseline`.

With `phase4_baseline` set to that source-tree root, run:

```sh
python3 tests/material_host_cleanup/check_batch_thermal.py --baseline "$phase4_baseline"
python3 tests/material_host_cleanup/check_batch_snow_host.py --baseline "$phase4_baseline"
python3 tests/material_host_cleanup/check_canopy_diagnostics.py --baseline "$phase4_baseline/main/MOD_LeafInterception.F90"
python3 tests/material_host_cleanup/check_batch_bgc_nh4.py --baseline "$phase4_baseline"
python3 tests/material_host_cleanup/check_batch_bgc_surface.py --baseline "$phase4_baseline"
python3 tests/material_host_cleanup/check_batch_routing.py --baseline "$phase4_baseline"
python3 tests/material_host_cleanup/check_runtime_history.py --baseline "$phase4_baseline"
```

These exercise actual phase-change/interception/snow routines, NH4 and clip helpers,
route callback ordering/publication, and history-count compatibility against frozen
host code. Complete methane aggregation bodies must also match their frozen source.
Existing lifetime, soil-column, boundary, physics and real NetCDF/MPI checks remain
part of the unified phase. The history-field generator recognizes the material
sample-count allocation hook and still checks all 525 host fields.

The coupled runner now compiles all three programs (`colm.x`, `mkinidata.x`,
`mksrfdata.x`) before the one continuous/split simulation. `material_naming/run_build.sh`
adds clean PC, bare-IGBP/LULCC, SinglePoint and UNSTRUCTURED builds alongside the
MT-off build. Build flags and exact job evidence are in the phase-4 results report.


## Optional host scaffolding cleanup (phase 5)

The phase-5 baseline is the complete verified phase-4 working tree, retained at
`.omx/artifacts/material-host-cleanup-phase5-20261031/baseline`. Run:

```sh
python3 tests/material_host_cleanup/check_optional_scaffolding.py --baseline "$phase5_baseline"
python3 tests/material_host_cleanup/check_phase_bookkeeping.py
```

The new comparison extracts actual MAIN branch structure and irrigation keyword
mapping, uses mutating water stubs, runs real work/diagnostic storage with the
actual THERMAL entry initializer, and compares old-wrapper/new-direct history
callback arguments. Eight configurations (CROP/Grid on/off, O0/O2) compare 128
water branches, 96 thermal lifecycles and 72 history calls bitwise. Water solvers
and history handlers are recording stubs; actual physics and real NetCDF/MPI
remain covered by the separate existing checks. Work preparation now leaves
thermal observations to THERMAL's initializer; water snapshots and water reset
remain owned by preparation. PhaseChange still supports standalone initialization.

Interface-only factoring also requires equal preprocessed tokens and statement
boundaries with MaterialTransport on/off. Session plans and validation notes stay
local; build evidence is retained under the corresponding local artifact directory.


## Consolidated material stages (phase 6)

Freeze the verified phase-5 tree before changing entry signatures; retained baseline:
`.omx/artifacts/material-host-cleanup-phase6-20261031/baseline`.

```sh
python3 tests/material_host_cleanup/check_stage_consolidation.py --baseline "$phase6_baseline"
python3 tests/material_host_cleanup/check_batch_snow_host.py --baseline "$phase4_baseline" --stage-baseline "$phase6_baseline"
```

The stage comparison compiles actual extracted old/new MAIN calls and LandPhase
routines, real Work/diagnostic storage, canopy partition and snowmelt carry, with
recording soil/evaporation/flood services. It records all physics arguments and
optional presence, callback order, material stores and final work snapshots across
Grid/CROP on/off and O0/O2. It covers 0/3 materials, 0–3 snow layers, soil/wetland,
dynamic wetland on/off and allocated/unallocated flood credit. Large transport
solvers are stubs in this check; this is a boundary equivalence test.

The snow stage comparison executes actual standard/SNICAR topology routines and
compares the new normalized entry with the frozen MAIN sequence, including empty
snow and snow layers collapsing during combine. Legacy land-boundary/scaffolding
checks retain their earlier baselines and adapt to the new owning entry.
Executed checks and build evidence are recorded in local artifacts and the PR summary.
