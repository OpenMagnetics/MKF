# Sprint: Add PLECS Equivalent Component Exporter (#55)

**Issue:** https://github.com/OpenMagnetics/MKF/issues/55
**Status:** WIP — Planning complete, implementation pending
**Date:** 2026-03-10

---

## Background

MKF already exports magnetic equivalent circuits for **SIMBA** (JSON/.jsimba), **NgSpice** (SPICE netlist),
**LTSpice** (SPICE netlist + .asy symbol), and **NL5** (SPICE netlist). Issue #55 requests adding PLECS support.

### PLECS Magnetic Domain — Key Differences

PLECS uses the **permeance-capacitance analogy** (not reluctance-resistance):

| Concept | Reluctance analogy (SIMBA) | Permeance-capacitance analogy (PLECS) |
|---------|---------------------------|---------------------------------------|
| Core/gap element | Reluctance (resistance-like) | Permeance (capacitance-like) |
| Winding coupling | Ideal transformer | Gyrator |
| MMF | Voltage source | Through variable (current-like) |
| Flux | Current | Across variable (voltage-like) |
| Energy storage | Inductors | Capacitors |
| Losses | Resistors in reluctance domain | Magnetic resistance (F²/Rm) |

**Permeance** P = μ₀·μᵣ·A/l = 1/Reluctance

### PLECS Magnetic Component Library

| Component | Parameters | Notes |
|-----------|-----------|-------|
| Constant Permeance | P [Wb/A·t] | Core legs, air gaps |
| Variable Permeance | P(t), dP/dt, Φ (3-element vector) | Saturable/hysteretic cores |
| Saturable Core | BH curve, area, length | Uses Variable Permeance internally |
| Hysteretic Core | Jiles-Atherton params | Uses Variable Permeance internally |
| Winding | N (turns) | Gyrator interface to electrical domain |
| Magnetic Resistance | Rm [A·t·s/Wb] | Core losses (series with permeance) |
| Magnetic Ground | — | Reference node |

### PLECS File Format

PLECS uses `.plecs` files (XML-based). Structure needs reverse-engineering from example files,
but the key approach is generating a PLECS **subsystem** containing magnetic-domain components.

**Alternative:** PLECS supports **Simulink/C-Script blocks** — we could generate a C-Script block
that wraps MKF calculations directly instead of building the magnetic circuit graphically.

---

## Architecture Decision

Follow the existing pattern: create `CircuitSimulatorExporterPlecsModel` inheriting from
`CircuitSimulatorExporterModel`, analogous to the SIMBA exporter but outputting
permeance-capacitance topology instead of reluctance network.

### Mapping from MKF Magnetic to PLECS Components

| MKF concept | SIMBA component | PLECS equivalent |
|-------------|----------------|-----------------|
| Core leg | `Linear Core` (μᵣ, A, l) | `Constant Permeance` (P = μ₀·μᵣ·A/l) |
| Air gap | `Air Gap` (μᵣ=1, A, l) | `Constant Permeance` (P = μ₀·A/l) |
| Winding | `Winding` (N) | `Winding` (N) — gyrator |
| Core losses | R-L ladder network | `Magnetic Resistance` in series |
| AC winding losses | R-L ladder / fracpole | Electrical-domain ladder (same as NgSpice/LTSpice) |
| Leakage inductance | Inductor (electrical) | Inductor (electrical) or leakage permeance |
| Ground | `Magnetic Ground` | `Magnetic Ground` |

---

## Sprint Tasks

### Phase 1: File Format Research & Scaffolding (1 session)

- [ ] **T1.1** Obtain a PLECS example file with magnetic components and reverse-engineer the XML schema
  - Ask user for a `.plecs` file with a simple transformer, OR
  - Use PowerForge export as reference
  - Document XML element names, attributes, coordinate system
- [ ] **T1.2** Add `PLECS` to `CircuitSimulatorExporterModels` enum
- [ ] **T1.3** Create `CircuitSimulatorExporterPlecsModel` class stub in `CircuitSimulatorInterface.h`
  - Inherit from `CircuitSimulatorExporterModel`
  - Declare `export_magnetic_as_subcircuit()` and `export_magnetic_as_symbol()`
- [ ] **T1.4** Add factory case for PLECS in `CircuitSimulatorExporterModel::factory()`

### Phase 2: Core Implementation (1-2 sessions)

- [ ] **T2.1** Implement XML generation helpers:
  - `create_permeance(P, coordinates, name)` — constant permeance element
  - `create_winding(N, coordinates, name)` — gyrator winding
  - `create_magnetic_resistance(Rm, coordinates, name)` — loss element
  - `create_magnetic_ground(coordinates, name)`
  - `create_connection(from, to)` — magnetic domain wiring
- [ ] **T2.2** Implement `export_magnetic_as_subcircuit()`:
  - Convert core geometry → permeance values (P = μ₀·μᵣ·A/l per leg)
  - Convert air gaps → permeance values (P = μ₀·A_gap/l_gap)
  - Place winding gyrators
  - Add electrical-domain AC resistance ladder/fracpole (reuse existing logic)
  - Add leakage inductance
  - Wire magnetic ground
  - Handle multi-column cores (E-core, UR-core, toroidal)
- [ ] **T2.3** Implement `export_magnetic_as_symbol()` (optional, if PLECS supports custom symbols)
- [ ] **T2.4** Add core loss modeling via `Magnetic Resistance`:
  - Convert core loss coefficients to equivalent Rm value(s)
  - Option: frequency-dependent Rm via multiple branches

### Phase 3: Testing (1 session)

- [ ] **T3.1** Create test fixtures:
  - `test_circuitsimulatorexporter_plecs_simple_inductor.json`
  - `test_circuitsimulatorexporter_plecs_transformer.json`
  - `test_circuitsimulatorexporter_plecs_ep_core.json`
  - `test_circuitsimulatorexporter_plecs_toroidal.json`
- [ ] **T3.2** Add test cases in `TestCircuitSimulatorInterface.cpp`:
  - Basic inductor export (single winding, single gap)
  - Transformer export (two windings, no gap)
  - Gapped E-core with multiple columns
  - Toroidal core
  - Verify permeance values match expected P = μ₀·μᵣ·A/l
- [ ] **T3.3** Cross-validate: compare PLECS inductance/reluctance with SIMBA export for same magnetic

### Phase 4: Documentation & Integration (1 session)

- [ ] **T4.1** Add PLECS to docs (algorithms or models section)
- [ ] **T4.2** Add PLECS example to `getting-started/examples.md`
- [ ] **T4.3** Close issue #55

---

## Open Questions

1. **File format:** We need a sample `.plecs` file with magnetic components to confirm XML schema.
   The PLECS manual PDF is behind binary encoding. User may need to provide one.
2. **Saturation modeling:** Should we use `Constant Permeance` (linear) or `Saturable Core`
   (nonlinear BH curve) for the initial implementation?
   - Recommendation: Start with `Constant Permeance` (matches SIMBA's `Linear Core`), add
     saturable core as a follow-up.
3. **Coordinate system:** PLECS uses a graphical schematic layout. Need to understand the
   coordinate/positioning system from a real file.
4. **C-Script alternative:** If XML reverse-engineering proves difficult, we could generate a
   PLECS C-Script block that calls MKF functions directly. This is more portable but less
   visually integrated.

---

## References

- [PLECS Magnetic Domain Overview](https://www.plexim.com/products/plecs/magnetic)
- [Permeance-Capacitance Analogy Paper (IEEE)](https://ieeexplore.ieee.org/document/6251786/)
- [PLECS Magnetic Modeling PDF](https://www.plexim.com/sites/default/files/magnetic_modeling_PLECS.pdf)
- [PLECS User Manual](https://www.plexim.com/sites/default/files/plecsmanual.pdf)
- [Gyrator-Capacitor Model (Wikipedia)](https://en.wikipedia.org/wiki/Gyrator%E2%80%93capacitor_model)
- Existing SIMBA exporter: `src/processors/CircuitSimulatorInterface.cpp:769`
- Existing NgSpice exporter: `src/processors/CircuitSimulatorInterface.cpp:1222`
