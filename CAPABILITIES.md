# MKF Capabilities

**Purpose:** authoritative list of what MKF already does, so downstream agents
(Heaviside, PyOpenMagnetics, MVB++, WebFrontend) don't reinvent magnetics
math or guess at APIs. **If something you need isn't here, ask — don't
implement it downstream.**

**Hard rule:** all magnetics math lives in MKF. No `B_sat·N·A_e/L` or
equivalents downstream. Use MKF / PyOM, or add a `MagneticFilter*` here.

Conventions:
- Names below are **C++ symbols**. PyOM names usually mirror them in
  snake_case. The canonical Python surface is in
  `src/PyMKF/PyMKFWrapper.cpp` (see §9) and in
  `/home/alf/OpenMagnetics/PyMKF/PyOpenMagnetics.pyi` (read PyMKF's
  AGENTS.md first — the stub is a signature reference, not a usage guide).
- Methods listed are the ones a downstream caller would actually invoke;
  getters/setters/ctors omitted.

---

## 1. Converter / topology design

Source: `src/converter_models/`. Each topology takes electrical specs and
produces magnetic requirements + (optionally) an ngspice netlist.

**Common interface** (`Topology.h` base): `process_design_requirements()`,
`process_operating_points()`, `run_checks()`, `generate_ngspice_circuit()`.
Per-operating-point diagnostics are exposed as `lastDutyCycle`, `lastMode`,
`lastZvsMargin*`, and `perOpName[]` / `perOpDutyCycle[]` accessors — read
these instead of re-deriving.

Shared analytical kernel: `PwmBridgeSolver` (duty-cycle loss / ZVS margin /
resonant transition time for PSFB and three-level half-bridge).

**Topologies** (all in `src/converter_models/`):

| Topology | Class | Notes |
|---|---|---|
| Buck | `Buck` | Hard-switched step-down |
| Boost | `Boost` | Hard-switched step-up |
| Cuk / SEPIC / Zeta | `Cuk`, `Sepic`, `Zeta` | Fourth-order; Cuk has V1/V2/V3 (uncoupled/coupled/isolated) variants |
| Four-switch buck-boost | `FourSwitchBuckBoost` | |
| Flyback | `Flyback` | Multi-output; CCM/DCM/QRM/BMO |
| Isolated buck / buck-boost | `IsolatedBuck`, `IsolatedBuckBoost` | Flybuck-style |
| Single / two-switch forward | `SingleSwitchForward`, `TwoSwitchForward` | Demag winding; CCM/DCM |
| Active-clamp forward | `ActiveClampForward` | Cclamp across primary (not to GND) |
| Push-pull | `PushPull` | Center-tapped primary, 180° drivers |
| Phase-shifted full / half bridge | `PhaseShiftedFullBridge`, `PhaseShiftedHalfBridge` | ZVS; three rectifier types |
| Asymmetric half bridge | `AsymmetricHalfBridge` | |
| LLC / CLLC / CLLLC / SRC | `Llc`, `Cllc`, `Clllc`, `Src` | Resonant; Nielsen TDA solver |
| DAB | `Dab` | SPS/EPS/DPS/TPS phase-shift modulation |
| Weinberg | `Weinberg` | |
| PFC / Vienna | `PowerFactorCorrection`, `Vienna`, `PfcControllerDesign` | CCM/DCM/CrCM |
| Common-mode choke | `CommonModeChoke` | |
| Differential-mode choke | `DifferentialModeChoke` | |
| Current transformer | `CurrentTransformer` | |

**"Advanced" pattern:** topologies with user-specified `L` / turns-ratio
have a sibling `Advanced<Topology>` overriding `process_design_requirements()`
to emit user values directly instead of auto-sizing.

**Do not re-derive** duty cycle, peak/RMS currents, volt-seconds, or
operating modes downstream — call the topology.

---

## 2. Advisers (selection / ranking)

Source: `src/advisers/`. Per-candidate scoring is delegated to
`MagneticFilter*` classes — adviser orchestration must **never** inline
the math.

| Class | Intent | Key methods |
|---|---|---|
| `CoreAdviser` | Pick a core from a catalog given Inputs | `get_advised_core()`, `calculate_gapping_constraints()`, `optimize_gap_and_turns_binary_search()`, `filter_available_cores_power_application()`, `filter_available_cores_suppression_application()` |
| `CoilAdviser` | Pick wire + section layout | `get_advised_coil()`, `get_advised_sections()`, `score_magnetics()` |
| `WireAdviser` | Pick a wire (round/litz/foil/rect/planar) | `get_advised_wire()`, `get_advised_planar_wire()`, plus filter methods: `filter_by_area_no_parallels()`, `filter_by_effective_resistance()`, `filter_by_skin_losses_density()`, `filter_by_proximity_factor()` |
| `MagneticAdviser` | Top-level: core + coil + scoring | `get_advised_magnetic()`, `get_advised_magnetic_fast()`, `get_advised_magnetic_from_converter()`, `score_magnetics()`, `preview_magnetic()` |
| `CoreCrossReferencer` | Find equivalent cores | (see header) |
| `CoreMaterialCrossReferencer` | Find equivalent core materials | (see header) |

**`MagneticFilter*` (per-candidate scoring, in `src/advisers/`):**
`MagneticFilterSaturation`, `MagneticFilterLosses`, `MagneticFilterArea`,
`MagneticFilterImpedance`, `MagneticFilterEnergyCost`,
`MagneticFilterGeometry`, `MagneticFilterPhysical`,
`MagneticFilterAdvanced`. If a needed score is missing, add a filter
here — don't compute it in the adviser body.

Default sat-margin: **1.2** (Maniktala / production safety), changed from 1.0 on 2026-05-25.

---

## 3. Physical models (the math layer)

Source: `src/physical_models/`. **All magnetics math lives here.** Do not
reimplement these formulas downstream.

| Quantity | Class | Key methods |
|---|---|---|
| Inductance matrix (self + mutual) | `Inductance` | `calculate_inductance_matrix()`, `calculate_leakage_inductance()` |
| Magnetizing inductance | `MagnetizingInductance` | `calculate_inductance_from_number_turns_and_gapping()`, `calculate_gapping_from_inductance()`, `calculate_number_turns_from_gapping_and_inductance()` |
| Leakage inductance | `LeakageInductance` | `calculate_leakage_inductance()`, `calculate_leakage_inductance_all_windings()`, `calculate_leakage_magnetic_field()` |
| Reluctance (gap + ungapped + air-cored) | `Reluctance` | `get_gap_reluctance()`, `get_ungapped_core_reluctance()`, `get_air_cored_reluctance()` (multiple models: Zhang / Muehlethaler / Partridge / Stenglein) |
| Magnetic field (H, B distribution) | `MagneticField` | `calculate_magnetic_field_strength_field()`, `calculate_magnetic_flux()`, `calculate_magnetic_flux_density()`, `factory()` |
| Magnetic energy | `MagneticEnergy` | `get_ungapped_core_maximum_magnetic_energy()`, `calculate_core_maximum_magnetic_energy()`, `calculate_gap_length_by_magnetic_energy()` |
| Stray capacitance | `StrayCapacitance` | `calculate_capacitance()` (Massarini / Duerdoth / Albach / Koch) |
| Impedance (R + jωL, Q, SRF) | `Impedance` | `calculate_impedance()`, `calculate_q_factor()`, `calculate_self_resonant_frequency()`, `calculate_minimum_number_turns()` |
| Initial permeability | `InitialPermeability` | `get_initial_permeability()`, plus `_temperature_dependent` / `_frequency_dependent` variants |
| Amplitude / complex permeability | `AmplitudePermeability`, `ComplexPermeability` | `get_amplitude_permeability()`, `get_complex_permeability()`, `get_hysteresis_loop()` (Roshen) |
| Core losses (Steinmetz / iGSE / Roshen / loss maps) | `CoreLosses` | model factory + loss-per-unit-volume |
| Winding losses (total) | `WindingLosses` | `calculate_losses()`, `calculate_effective_resistance_of_winding()`, `calculate_resistance_matrix()` |
| DC ohmic losses | `WindingOhmicLosses` | `calculate_dc_resistance_per_meter()`, `calculate_effective_resistance_per_meter()`, `calculate_ohmic_losses()` |
| Skin-effect losses | `WindingSkinEffectLosses` | `calculate_skin_depth()`, `calculate_skin_effect_losses()` (Dowell / Wojda) |
| Proximity-effect losses | `WindingProximityEffectLosses` | `calculate_proximity_effect_losses()` (Rossmanith / Wang / Ferreira / Albach) |
| Core temperature (analytical) | `CoreTemperature` | `get_core_temperature()` (Maniktala / Kazimierczuk / TDK / Dixon / Amidon) |
| Full thermal network solver | `Temperature`, `ThermalNode`, `ThermalResistance` | network factory; `calculateConductionResistance()` |
| Resistivity (T-dependent) | `Resistivity` | `get_resistivity()` |

**Two thermal surfaces by design:** `CoreTemperature` (quick analytical) vs
`Temperature` (full network solver). Pick deliberately.

**Saturation current:** use PyOM's `calculate_saturation_current` (delegates
to `Magnetic::calculate_saturation_current` / `MagneticFilterSaturation`).
Never inline `B_sat·N·A_e/L`.

---

## 4. Constructive models (the magnetic itself)

Source: `src/constructive_models/`

| Class | Intent | Key methods |
|---|---|---|
| `Magnetic` | Full magnetic (core + coil + bobbin) | `get_bobbin()`, `get_wires()`, `get_wire()`, `calculate_saturation_current()`, `fits()` |
| `Core` | Geometry + material + gapping + effective params | `resolve_material()`, `resolve_shape()`, `process_gap()`, `get_effective_length()`, `get_effective_area()`, `get_initial_permeability()`, `create_quick_core()` |
| `CorePiece` | Half-core geometry, volumetric loss split | `process()`, `calculate_column_cross_section()`, `calculate_core_part_volumes()`, `calculate_core_loss_fractions()` |
| `Bobbin` | Winding window geometry | `process_data()`, `get_winding_window_dimensions()`, `get_winding_window_area()`, `create_quick_bobbin()` |
| `Coil` | Multi-winding turn layout + inter-layer insulation | (large header — see `Coil.h` for winding indexing + insulation section API) |
| `Wire` | Solid / litz / foil / rect / planar; coating; packing | (large header — top-level via ctors + inherited `MAS::Wire` methods) |
| `Insulation` / `InsulationMaterial` | IEC 60664 clearance / creepage; tape thickness | `extract_available_thicknesses()`, `get_dielectric_strength_by_thickness()`, `calculate_clearance()`, `calculate_creepage_distance()` |
| `NumberTurns` | Turn-count solver from ratio constraints | `increment_number_turns()`, `get_next_number_turns_combination()` |
| `Mas` / `MasMigration` | Top-level container (Inputs + Magnetic + Outputs) | `get_inputs()`, `get_magnetic()`, `get_outputs()` |

---

## 5. Inputs / Outputs / Simulation

Source: `src/processors/`

| Class | Intent | Key methods |
|---|---|---|
| `Inputs` | Build operating points (waveforms, harmonics, DC bias) | `process()`, `get_operating_point()`, `get_primary_excitation()`; statics: `calculate_sampled_waveform()`, `calculate_processed_data()`, `calculate_harmonics_data()`, `prune_harmonics()`, `reflect_waveform()`, `create_quick_operating_point()` |
| `MagneticSimulator` | Full electromagnetic sim | `simulate()`, `calculate_core_losses()`, `calculate_magnetizing_inductance()`, `calculate_leakage_inductance()`, `calculate_winding_losses()` |
| `Sweeper` | Parametric sweeps (freq / temp / bias) | `sweep_impedance_over_frequency()`, `sweep_magnetizing_inductance_over_temperature()`, `sweep_core_losses_over_frequency()`, `sweep_winding_losses_over_frequency()`, `sweep_resistance_matrix_over_frequency()` |
| `Outputs` | Result container (inherits `MAS::Outputs`) | — |
| `NgspiceRunner` | Drive ngspice for transient sims | (see header) |
| `CircuitSimulatorInterface` | Generic circuit-sim interface | (see header) |
| `CircuitSimulatorLtspice` / `Ngspice` / `Nl5` / `Plecs` / `Simba` | Export netlists per format | (see headers) |

---

## 6. Plotting (SVG)

Exposed via PyMKF (see §9): `plot_turns`, `plot_sections_and_turns`,
`plot_field`, `plot_current_density`, `plot_wire`.

Underlying renderer: `Painter` in `src/support/` —
`paint_core()`, `paint_bobbin()`, `paint_coil_turns()`,
`paint_coil_sections()`, `paint_wire()`,
`paint_wire_with_current_density()`, `paint_magnetic_field()`,
`export_svg()`. Don't reimplement geometry rendering downstream.

---

## 7. Support / infra

Source: `src/support/`

| Class | Intent | Key methods |
|---|---|---|
| `DatabaseManager` | Singleton catalog access (cores, materials, shapes, wires, bobbins) | `getInstance()`, `loadAll()`, `loadCores()`, `loadCoreMaterials()`, `loadCoreShapes()`, `loadWires()`, `loadBobbins()`, `findCoreByName()`, `clearAll()`, `reset()` |
| `LibraryContext` | Per-call catalog override + type filtering (RAII) | `loadFromJson()`, `loadFromString()`, `loadFromFile()`, `applyScoped()`, `setMode()` (Merge / Replace) |
| `Cache<T>` | Generic component cache | `load()`, `read()`, `get()`, `clear()`, `references()`, `evict()` |
| `Settings` | Global simulation config (singleton) | painter modes, field grid precision, reluctance / core-temp model selection, margin tape / insulated-wire / layers knobs |
| `CoilMesher` | Field grid for proximity losses | `generate_mesh_inducing_coil()`, `generate_mesh_induced_coil()`, `generate_mesh_induced_grid()` |
| `Logger` | Structured logging | TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF; ConsoleSink, FileSink; module-based filtering |
| `Exceptions.h` | Exception types | — |

---

## 8. Defaults and constants

- `src/Defaults.h` — **every tunable knob with its default**. Check here
  before adding a magic number anywhere.
- `src/Constants.h` — physical constants.
- `src/Definitions.h` — generated MAS schema types (do not hand-edit; come
  from quicktype on the PSMA schema submodules CAS / SAS / RAS).

---

## 9. Python API (PyMKF / PyOpenMagnetics)

`src/PyMKF/PyMKFWrapper.cpp` exposes:

**Module functions:**
- Plotting: `plot_turns`, `plot_sections_and_turns`, `plot_field`,
  `plot_wire`, `plot_current_density`
- Settings: `get_settings`, `set_settings`, `reset_settings`
- Advisers: `get_advised_core(ctx, constraints, weights)`,
  `get_advised_magnetic(ctx, constraints)`

**Classes:**
- `LibraryContext` — `load_from_json`, `load_from_file`, `load_from_string`,
  `clear`, `empty`
- `AdviserConstraints` — type filters (`shape_family`, `core_material_type`,
  `wire_type`); `empty`
- `TypeFilterSet` — `allowed` / `blocked` sets for enum filtering
- `LoadMode` — `Merge` / `Replace`

The PyOpenMagnetics surface in `PyOpenMagnetics.pyi` is broader (all the
adviser / model methods above are exposed). **Read
`PyOpenMagnetics/AGENTS.md` first** — the stub is a signature reference,
not a usage guide (positional args only; some enum values are display
strings; flyback JSON fields are not what the stub suggests).

---

## What MKF does NOT do

- Component cost / BOM optimization beyond `MagneticFilterEnergyCost`.
- PCB layout / mechanical CAD.
- Controller-loop design beyond `PfcControllerDesign`.
- General-purpose SPICE solving — MKF *generates* netlists and *drives*
  ngspice, but isn't a SPICE engine itself.
- Topology selection — MKF designs *for* a chosen topology; picking one
  is upstream (Heaviside).

---

## How to use this from a downstream repo

Add to the downstream repo's `AGENTS.md` (or `CLAUDE.md`):

```markdown
## MKF capabilities
Before implementing any magnetics calculation or assuming a PyOM API,
consult @../MKF/CAPABILITIES.md (path may vary per repo). MKF is
authoritative for all magnetics math.
```
