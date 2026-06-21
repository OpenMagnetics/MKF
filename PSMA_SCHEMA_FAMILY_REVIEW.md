# PSMA Schema Family Review — CIAS, TAS, RAS, SAS (and how they mingle with MAS)

> Step 1 of the PSMA-family integration. This document records what each new
> schema repo is, how the schemas reference each other, where MAS sits in the
> graph, and what has to happen in MKF to consume them. No C++ codegen wiring is
> done yet — that is the follow-up step.

## 1. What was added

Four new git submodules from the **Power Supply Manufacturers Association**
(`github.com/Power-Supply-Manufacturers-Association`), alongside the existing
`MAS`, `CAS`, `PEAS`:

| Submodule | Pinned to | Role |
|-----------|-----------|------|
| `CIAS` | `main` | Circuit **brick** (`.subckt`): ports + components + connections |
| `TAS`  | `main` | Complete converter **deck**: inputs + topology + outputs + simulation |
| `RAS`  | `main` | **Resistor** component family |
| `SAS`  | `main` | **Semiconductor** component family (MOSFET/diode/IGBT/BJT) |

**All four are pinned to `main`** — the PEAS-family consolidation has already
merged to `main` in every PSMA repo. The pre-existing `PEAS` submodule was also
re-pinned from `feat/peas-family-consolidation` to `origin/main` in this step
(PEAS `main` is +3 commits: CONAS connector branch + `currencyAmount` rename — see §6).
CIAS's remote *default* branch is still the feat branch, but its `main` exists
and is ahead, so CIAS is pinned to `main` explicitly.

## 2. The family architecture

Every repo is a [JSON Schema 2020-12](https://json-schema.org/) model. They form
a strict tier hierarchy, each tier referencing the tier below **inline or by URI**:

```
COAS   converter wrapper          (Converter Agnostic Structure — product level:
 │                                 form-factor, cooling, environmental, around a TAS)
 ▼
TAS    runnable converter deck    (stages wiring CIAS bricks + stimulus + analyses)
 │  stage.circuit  → inline CIAS | URI
 ▼
CIAS   reusable circuit brick     (ports + components + connections)
 │  component.data → inline PEAS | URI
 ▼
PEAS   universal component container  (the abstract base) — 7 discriminator branches
 ├── MAS    magnetics      (this repo's home territory)        [own repo]
 ├── CAS    capacitors                                         [own repo]
 ├── SAS    semiconductors                                     [own repo]
 ├── RAS    resistors / varistors                              [own repo]
 ├── CONAS  connectors                                         [own repo, psma.com/conas]
 ├── controller   (inline PEAS branch today; future "CtAS")
 └── terminal     (inline PEAS branch today; proposed to fold into CONAS — see §8)
```

Two converter-vs-connector acronyms to keep straight: **COAS** = Converter
(wraps TAS), **CONAS** = Connector (a PEAS component family). See §6.4.

**PEAS is the hub.** Every component family is "a valid PEAS document". Each
family carries the same three-section shape — `inputs` (designRequirements +
operatingPoints) / `<component>` / `outputs` — established by PEAS.

## 3. Cross-`$ref` resolution — the key mechanism

Schemas reference each other by **absolute `$id` URI**, never relative path
across repos:

- `MAS` → `https://psma.com/mas/...`
- `CAS` → `https://psma.com/cas/...`
- `SAS` → `https://psma.com/sas/...`
- `RAS` → `https://psma.com/ras/...`
- `PEAS` → `https://psma.com/peas/...`
- `CIAS` → `https://psma.com/cias/...`
- `TAS` → `https://psma.com/tas/...`

`PEAS/schemas/peas.json` is the dispatcher. Its `oneOf` selects a family by which
discriminator property is present, and `$ref`s straight into that family's repo:

```
magnetic       → https://psma.com/mas/magnetic.json
capacitor      → https://psma.com/cas/capacitor.json
semiconductor  → https://psma.com/sas/SAS.json
resistor       → https://psma.com/ras/resistor.json
controller     → ./controller.json   (inline)
terminal       → ./terminal.json     (inline)
```

…and `allOf` `if/then` clauses pin `inputs.designRequirements` to the matching
family's `inputs/designRequirements.json`.

**Measured cross-repo dependency counts** (`$ref`s pointing outside the repo):

| Repo | depends on |
|------|------------|
| CIAS | PEAS (×2), self |
| TAS  | CIAS (×1), PEAS (×1), self (×5) |
| RAS  | PEAS (×42), self (×7) |
| SAS  | PEAS (×81), self (×9) |

RAS and SAS lean **heavily** on `PEAS/schemas/utils.json` for shared
`datasheetInfo` building blocks (dimensionWithTolerance, curve, part/business/
mechanical/thermal blocks) rather than duplicating them — the same move MAS made.

### Transitive closure
Resolving a **TAS** document pulls in *the entire family*:
`TAS + CIAS + PEAS + MAS + CAS + SAS + RAS + controller + terminal`.
Resolving a bare **RAS** or **SAS** document needs only that repo + PEAS.

## 4. How MAS mingles

MAS is **already** a first-class member of this family, not an outsider:

1. **MAS lives under the `psma.com/mas/` `$id` namespace** — identical scheme to
   the PSMA repos, so PEAS's `https://psma.com/mas/magnetic.json` resolves to the
   `MAS` submodule on disk.
2. **MAS now inherits shared types from `PEAS/schemas/utils.json`** (per the
   utils.json header: "MAS now references this file too… magnetics-specific
   construction/physics types remain in MAS"). This is exactly the
   PEAS-integration that landed in MKF commit `1b878c5b`.
3. **PEAS dispatches the `magnetic` branch into MAS** and pins
   `inputs.designRequirements` to `mas/inputs/designRequirements.json`.

So "mingling" is **already the design**: PEAS is the seam, MAS is one spoke,
and CIAS/TAS sit above PEAS treating a magnetic exactly like any other component.
A transformer/coupled-inductor is a *single multi-winding PEAS component* whose
coupling lives inside its MAS model — CIAS deliberately has **no
magnetic-coupling edges** (coupling is never a circuit net).

## 5. How MKF would consume them (the quicktype pattern)

MKF turns schemas into C++ via `quicktype` in `CMakeLists.txt`. Each family
becomes its own namespaced header (`MAS.hpp` in `namespace MAS`, `CAS.hpp` in
`namespace CAS`). For CAS, the build already passes `-S` flags for every PEAS
schema CAS references, so quicktype resolves the `psma.com/...` `$id` URIs
**offline** (no network during CI).

To add the new families later, the same recipe applies — one
`add_custom_command` per family, with `-S` flags covering its transitive
closure:

- `RAS.hpp` / `SAS.hpp`: need RAS/SAS own schemas + all PEAS schemas they `$ref`.
- `CIAS.hpp`: needs CIAS + `peas.json` + `peas/utils.json`.
- `TAS.hpp`: heaviest — needs TAS + CIAS + **the whole PEAS dispatch closure**
  (MAS, CAS, SAS, RAS component + designRequirements schemas) because a TAS
  component can be any family.

This is deferred to the next step on purpose.

## 6. Findings / inconsistencies to raise upstream

Per the project's "surface broken things, don't route around them" rule, these
are flagged rather than silently worked around:

The findings here are **README/doc drift only** — the schemas themselves are
consistent. The wrong READMEs were corrected as part of this step.

1. **SAS schema is correct as-is (not a bug).** PEAS dispatches `semiconductor`
   → `https://psma.com/sas/SAS.json`, and `SAS.json` exposes
   `mosfet / diode / igbt / bjt` as **top-level discriminator properties** —
   exactly one of which is populated per document, and *that property is the
   device type*. There is intentionally **no** `semiconductor` wrapper object and
   **no** separate `semiconductor.json` component file. This is the deliberate SAS
   shape, confirmed by the maintainer.
   *(README fixed:* the SAS README previously claimed a `semiconductor.json` with
   a `oneOf` on `part.deviceType` inside a `semiconductor` object — that file
   never existed and the pattern was wrong. Corrected to match the schema.)

2. **TAS README was a version behind the schema (now fixed).** README described
   TAS v1 (`inputs / components / outputs`, flat netlist). On-disk `TAS.json` is
   **v2** (`inputs / topology / outputs / simulation`, stages instantiating CIAS
   bricks). The schema-overview section of the README was rewritten to v2; the
   database stats / units sections remain valid.

3. **All repos are PEAS-consolidated in `main`.** Initially CIAS/PEAS were pinned
   to `feat/peas-family-consolidation`; the consolidation has since merged to
   `main` in every repo, so all submodules (including the pre-existing PEAS) are
   now pinned to `main` and `.gitmodules` sets `branch = main` for them so
   `git submodule update --remote` tracks main going forward. PEAS `main` adds a
   7th discriminator branch — **CONAS (connector)**, a separate PSMA repo
   `psma.com/conas/` — plus a `cost → currencyAmount` rename in `utils.json`; both
   will need MKF C++ adaptation when PEAS/the new families are wired into the
   quicktype build (follow-up step).

4. **Two distinct `*AS` acronyms — not a collision (resolved).** The PSMA family
   has both **COAS = Converter Agnostic Structure** (a *wrapper around TAS* at the
   converter / product level — owns form-factor, cooling, environmental
   qualification) and **CONAS = Connector Agnostic Structure** (the connector
   component family, `psma.com/conas/`). The `COAS → CONAS` rename commit
   correctly renamed only the *connector* discriminator `$ref`s in `peas.json`;
   the three `COAS` mentions in `PEAS/utils.json` (`mountingStyle`,
   `environmentalQualification`, `coolingMethod`) genuinely refer to COAS=converter
   and are **left as-is, intentionally**.

None of these block step 1 (adding submodules + understanding).

## 7. Status

- [x] Branch created, four submodules added & pinned
- [x] Schema families read and understood
- [x] Cross-`$ref` graph mapped; MAS↔PEAS seam confirmed
- [x] Integration (quicktype) path identified
- [x] All submodules pinned to `main`; `.gitmodules` set to track `main`
- [x] Wrong READMEs fixed and **pushed to `main`** in SAS/TAS/RAS (SAS device
      discriminator, TAS v2 overview, RAS varistor + file lists); submodule pins
      bumped to the new commits
- [ ] **Next steps (follow-up):** wire chosen families into `CMakeLists.txt`
      quicktype (handle PEAS CONAS branch + `currencyAmount`); add
      `MasMigration`-style aliases; consume in advisers/converter models as needed.

## 8. Proposal — merge `terminal` into CONAS

### The overlap
PEAS today has **two** ways to represent a board-boundary conductor:

- the inline **`terminal`** branch — "the abstract board-boundary point" (so every
  CIAS/TAS `connection.endpoint` lands on a *real* pin, not just a label);
- the new **CONAS `connector`** family — "the physical orderable part".

`peas.json` itself admits the redundancy: *"Connector — defined by CONAS.
Coexists with the lightweight 'terminal' branch: terminal is the abstract
board-boundary point, connector is the physical orderable part."*

But the field sets already collide. `terminal.category` enumerates
`header`, `edgeCardFinger`, `boardToWireConnector`, `boardToBoardConnector`,
`rfConnector`, `screwTerminal`, `screwLug`, `bananaJack`, `solderPad` — and all of
those except the bare `solderPad` are exactly CONAS `familyDetails` members
(pin headers, board-to-board, wire-to-board, terminal blocks, card edge, RF,
power). Two schemas, one physical thing → consumers must handle both, and a part
can be modelled two ways.

### Recommendation: **replace `terminal` with CONAS** (single connector family)
Drop the inline `terminal` discriminator branch; let **CONAS be the one family**
for every board-boundary conductor. To preserve terminal's "lightweight abstract
boundary" role, CONAS already has the right shape — its **Tier-1 catalog**
(`datasheetInfo`) is the orderable part, and a **minimal CONAS document**
(just `category` + `currentRating`/`voltageRating`/`contactResistance`, no
geometry/contactSystem, both of which are already optional-but-closed) covers the
abstract-boundary case that `terminal` existed for. The bare `solderPad` (not a
connector) becomes either a minimal CONAS `familyDetails: solderPad` member or a
documented degenerate case.

Concrete edits (all upstream in PSMA, no MKF code yet):
1. **PEAS `peas.json`** — remove the `terminal` `oneOf` branch and its
   `allOf/if-then` designRequirements clause; keep only the `connector` → CONAS
   branch. (CONAS is already wired in.)
2. **PEAS** — delete `schemas/terminal.json`; migrate any field `terminal` had
   that CONAS lacks (e.g. `numPositions`, `wireGaugeRange`, `pitch` — verify
   against CONAS `datasheetInfo`/`utils`).
3. **CONAS** — ensure `familyDetails` covers all former `terminal.category`
   values; add `solderPad`/bare-pad if it must remain representable; confirm a
   minimal (Tier-1-only) doc validates so the abstract-boundary use still works.
4. **CIAS/TAS** — the "every endpoint is a real pin" invariant now resolves
   through a CONAS connector instead of a `terminal`; update prose/examples.

### Alternative (if a zero-data boundary marker must stay distinct)
Keep `terminal` strictly as a **base**: make CONAS `connector` *extend* terminal
(`allOf: [terminal, {datasheetInfo, geometry, contactSystem}]`) so there is one
shared boundary type and CONAS only adds the catalog/3D layers. More backwards-
compatible, but keeps two names for one concept — weaker than a clean merge.

> Note: this merge concerns **CONAS** (connector) only. It is unrelated to
> **COAS** (Converter Agnostic Structure — the wrapper around TAS), whose `COAS`
> mentions in `PEAS/utils.json` are correct and stay as-is.
