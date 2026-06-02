# Real Winding Investigation (U/Z winding, transition-wire squeeze, connection losses)

**Status:** investigation / design — no code changed yet.
**Date:** 2026-06-01.
**Goal:** Implement "real winding" in MKF — model each wire's entrance/exit/connection, let the
user choose how layers are wound (U vs Z) when a section has >1 layer, and (follow-up) account for
the DC + skin + proximity losses of the extra connection/lead copper out to the bobbin termination,
as defined in the design requirements.

This document is the durable record of the investigation behind that feature. Two web-research
workflows were still running when this was first written; their synthesized results are appended in
§7 / §8 once available.

---

## 1. Feature scope (as clarified by the user, in order)

1. **Real geometry.** OM currently places each turn at its *ideal* position (where it would seat if
   no other turns/layers existed). In reality the lead/transition wires cross over layers and steal
   radial space — **asymmetrically**.
2. **Worked PSPS example** (20-turn primary in 2 layers, interleaved P-S-P-S with a 2-layer
   secondary): layer 1 = primary half (top→bottom), its continuation wire stays at the bottom to
   climb to layer 3; layer 2 (1st secondary half) is squeezed by **~2 wire diameters** (a crossing
   wire on top = primary entry, one on the bottom = primary layer1→layer3 continuation); layer 3 =
   primary half 2 (bottom→top), its exit lead ends on **top at a different angular position** than
   the entry; layer 4 (2nd secondary half) is squeezed by **~1** wire diameter, not two.
3. **User-selectable U vs Z winding** when a section has more than one layer, to control parasitic
   capacitance (and the geometry/length consequences). **Narrowed current deliverable:** make the
   selector exist and make the `Coil.cpp` winder honor it.
4. **Connection losses.** DC + skin (+ proximity) of the extra lead/fly-back copper, summed out to
   the bobbin termination, per the design requirements. *(Follow-up; the schema already supports
   the length — see §6.)*

**User constraint:** capacitance is computed from the voltage drop between turns, so ordering and
naming the turns correctly for U vs Z makes the capacitance result fall out automatically — **do not
re-investigate or special-case capacitance.**

---

## 2. Terminology findings (fact-checked web research)

- **"Dragback" is not a standard term.** It appears exactly once, undefined, in a Würth EMC
  transformer-design deck, meaning *"route both ends of an internal wire-wound shield back to a
  bobbin pin"* for automated termination (vs. the manual alternative of burying one shield lead). A
  broad web search for the word in a winding context otherwise returns only soccer / FIFA / fishing.
  → For this feature we should **not** lean on the word "dragback"; we model **lead / transition /
  connection wires**.
- The concept of bringing a **buried inner-layer winding end out to a terminal by folding the lead
  back over the winding and catching it at the bobbin collar** is a real, **patented** technique
  (US8191240B2) — but the patent never names it "dragback." Real mechanism, no canonical name.
- **U vs Z winding** (terminology to be confirmed by the in-flight U/Z research workflow, §8):
  - **U / "back-and-forth" / standard / alternating:** the wire reverses direction each layer. The
    end of layer N is physically adjacent to the start of layer N+1 (short step-up at one side). At
    the turnaround the inter-layer voltage difference is small; at the far end it is ~2× the layer
    voltage. → generally **higher** inter-layer capacitance.
  - **Z / "fly-back" / return-wire:** every layer is wound the same direction; a return/fly-back
    wire carries the conductor back across the layer to restart the next layer at the same side.
    Inter-layer voltage difference is ~1 layer-voltage **everywhere** → generally **lower** and more
    uniform capacitance, at the cost of extra return-wire copper length.

---

## 3. Current MKF winder behavior (verified in code)

Entry: `Coil::wind()` → `wind_by_rectangular_turns()` (rectangular bobbin) or `wind_by_round_turns()`
(toroid). Hierarchy is built sections → layers → turns:

- `wind_by_rectangular_sections()` — `src/constructive_models/Coil.cpp:2890`
- `wind_by_rectangular_layers()` — `Coil.cpp:3835` (builds the Layer list per section; layer count
  from section dims / wire OD; sets each layer's coordinates/dimensions/orientation)
- `wind_by_rectangular_turns()` — `Coil.cpp:4526` (places the individual Turns)

**Key turn-placement facts** (`Coil.cpp:4564–4767`):
- For each conduction layer it picks a **start edge** from `layer.turnsAlignment`
  (CENTERED / INNER_OR_TOP / OUTER_OR_BOTTOM / SPREAD) and a **fixed increment**:
  - OVERLAPPING layer (turns vary in height): start near the top, loop does
    `currentTurnCenterHeight -= currentTurnHeightIncrement` → **top→bottom**.
  - CONTIGUOUS layer (turns vary in width): start near inner edge, loop does
    `currentTurnCenterWidth += currentTurnWidthIncrement` → **inner→outer**.
- **Every layer resets to the same start edge and the same direction.** There is no per-layer
  direction reversal anywhere.
- Turn **names** are assigned sequentially in placement order:
  `"<winding> parallel <p> turn <currentTurnIndex>"`, where `currentTurnIndex[winding][parallel]`
  increments per placed turn and **persists across layers** (`Coil.cpp:4682`, `:4702`). So the turn
  name encodes the electrical sequence = the order turns are pushed.
- `turn.set_orientation(TurnOrientation::CLOCKWISE)` is **hardcoded** (`Coil.cpp:4683`, `:4744`).
- Two winding-style branches: `WIND_BY_CONSECUTIVE_TURNS` (default) and
  `WIND_BY_CONSECUTIVE_PARALLELS` — **both** use the same start-edge + fixed-increment logic and so
  both need the U/Z change.
- `delimit_and_compact()` (`Coil.cpp:6145`) post-fits bounding boxes / slides sections together but
  never inserts gaps for crossing wires.

**Conclusion:** the current winder is implicitly **Z (fly-back)** — the last turn of a layer (e.g.
bottom) is electrically followed by the first turn of the next layer (top), i.e. a return jump — but
it does **not** model the return wire's length or the space it consumes. **U winding does not
exist.** Collision detection exists for toroids only (`get_collision_distances`, `Coil.cpp:5295`);
rectangular windings have none.

---

## 4. How U/Z would be implemented in the winder (design)

The change is small, localized, and backward-compatible.

1. **Determine each layer's index within its section** while iterating layers in
   `wind_by_rectangular_turns` (track previous section name; reset a counter when the section
   changes). Helper `get_layers_by_section(name)` exists at `Coil.cpp:787`.
2. **Decide reversal:** `reverse = (windingOrder == U) && (layerIndexInSection % 2 == 1)`.
   (Z = never reverse = current behavior.)
3. **Generic reversal** that works for every alignment, because the N forward positions are
   `start, start±inc, …, start±(N-1)·inc`. Reversing the electrical order just means starting at the
   last forward position and flipping the increment sign:
   ```
   if (reverse) {
       currentTurnCenterWidth  += (physicalTurnsInLayer - 1) * currentTurnWidthIncrement;
       currentTurnCenterHeight -= (physicalTurnsInLayer - 1) * currentTurnHeightIncrement; // loop uses -=
       currentTurnWidthIncrement  = -currentTurnWidthIncrement;
       currentTurnHeightIncrement = -currentTurnHeightIncrement;
   }
   ```
   Apply just before the turn-placement loop, in **both** winding-style branches.
4. **Turn names stay sequential** (unchanged), so for U the first turn placed in an odd layer lands
   at the edge adjacent to the previous layer's last turn → correct U adjacency. **Capacitance
   recomputes automatically** from the new turn coordinates + voltage-per-turn (no capacitance code
   touched — per the user's instruction).
5. Optionally set `turn.orientation` (the dormant clockwise/counterClockwise field) to reflect the
   reversal, for downstream/visualization — not required for capacitance.

This models the **electrical ordering and turn positions** of U vs Z. Modeling the **return-wire
copper length / radial space** (the squeeze of §1) is the separate, larger geometry task.

---

## 5. Where the U/Z selector should live (decided with user)

- **No U/Z field exists in MAS** (the schema's "flyback" hits are converter topologies; "u" is a
  core/bobbin *shape*). New fields are needed.
- **Resolution hierarchy (user spec, modeled on layer orientation):**
  1. **`Section.windingOrder`** (if present) — per-section override.
  2. else **bobbin winding-window `windingOrder`** — applies to all sections in that window.
  3. else **default = `Z`** (see ⚠ below — user chose `Z` over `U` to avoid changing existing
     results; current implicit behavior is `Z`).
- **New fields to add (schema + quicktype regen of `MAS.hpp`):**
  - `windingOrder` (optional) on **MAS `Section`** — `MAS/schemas/magnetic/coil.json` (`section` def
    ≈369–481).
  - `windingOrder` (optional) on the **bobbin winding window element** —
    `MAS/schemas/magnetic/bobbin.json` (`windingWindowElement`).
  - New enum type `windingOrder` (values TBD — see §9) in `MAS/schemas/utils.json` (or inline).
  - MAS is a submodule; edit the JSON then **regenerate `MAS.hpp` via quicktype** (memory
    `feedback_use_existing_schemas` — never hand-write wrapper classes).
- **Resolution helper** in `Coil` (mirrors `get_turns_alignment(sectionName)` at `Coil.cpp:787`+
  region): `get_winding_order(sectionName)` → section value ?? window value ?? `U`.
- **⚠ Default = Z (decided):** the current implicit behavior is **Z** (every layer same direction).
  Defaulting to **Z** preserves all existing capacitance / turn-coordinate results → **no test
  churn**. The user explicitly chose this over default-`U` once the churn implication was flagged.
  `U` is opt-in per section or per winding window.

---

## 6. MAS schema + MKF loss/capacitance surfaces (verified)

**MAS Turn / Layer / Section (ideal geometry only):**
- **Turn** (`coil.json` ≈110–200, `MAS.hpp:6514`): single `coordinates` (center), `dimensions`,
  `angle` (swept, for partial turns), `rotation`, `orientation` (cw/ccw), `additionalCoordinates`
  (toroids), `length`, `crossSectionalShape`, `coordinateSystem`. **No entry-vs-exit point, no
  angular start/end.**
- **Layer** (`coil.json` ≈201–311): `coordinates`+`dimensions` (bbox), `fillingFactor`,
  `turnsAlignment`, `windingStyle`, `orientation` (contiguous/overlapping), `type`.
- **Section** (`coil.json` ≈369–481): `coordinates`+`dimensions`, `fillingFactor`, `layersAlignment`,
  `layersOrientation`, `windingStyle`, `numberLayers`, `margin`, `group`.

**Connection length — ALREADY EXISTS (no new schema needed for it):**
- `coil.json:645` — `connection.length`: *"Length of the connection, from the exit of the last turn
  to the terminal. Unit: m."* on a partial winding's input/output `connections`.
- C++ `ConnectionElement` at `MAS.hpp:5581` (`get_length()` at ~5603).
- **Never consumed by any loss calc** today.

**Loss pipeline (exists, cleanly extensible) — from loss-surface scout:**
- `WindingLosses::calculate_losses` (`src/physical_models/WindingLosses.cpp:169`) orchestrates ohmic
  + skin + proximity.
- `WindingOhmicLosses::calculate_dc_resistance(Turn, Wire, temp)`
  (`WindingOhmicLosses.cpp:16`) uses `turn.get_length()`; per-winding aggregation at
  `WindingOhmicLosses.cpp:67–97`.
- Skin (`WindingSkinEffectLosses.cpp:209`) and proximity (`WindingProximityEffectLosses.cpp:187`)
  also multiply loss-per-meter by `turn.get_length()`.
- **Injection point for connection loss:** add `ConnectionElement.length` into the per-winding
  effective conductor length in `calculate_ohmic_losses` (and the skin/proximity equivalents).

**Capacitance (exists; do not modify):**
- `StrayCapacitance` (`src/physical_models/StrayCapacitance.h/.cpp`) with Massarini, Duerdoth,
  Albach, Koch, parallel-plate models; `calculate_capacitance(Coil)`, `calculate_voltages_per_turn`,
  `calculate_capacitance_among_turns`. Driven by turn coordinates + per-turn voltage → U/Z falls out
  of correct turn ordering automatically.

**Design requirements:**
- `designRequirements.strayCapacitance` **exists** (`MAS/schemas/inputs/designRequirements.json:41`,
  `dimensionWithTolerance` array) — capacitance targets already have a home.
- `terminalType` exists (PIN/THT/SMT/flyingLead/…); **no lead/connection-length field** in
  DesignRequirements — would need adding only if the termination length is to be a *requirement*
  input (the per-connection `length` already lives on the coil).

**U/Z encoding today:** `Turn.orientation` (cw/ccw) is **dormant** — set to CLOCKWISE, never read by
physics. `Layer.windingStyle` only controls turn-vs-parallel sequencing. No per-layer direction
reversal → real U/Z is not representable without the new field + winder change.

---

## 7. Real-geometry / state-of-the-art research (workflow `wf_c93f693a` — COMPLETE)

_This research is for the **larger "transition-wire squeeze" geometry feature (§1.1–1.2)**, which is
**not** in the current U/Z+connection-loss PR. Captured here for the follow-up. Full output:
`/tmp/claude-1000/.../tasks/weedaiyye.output`._

**Verified build-up anchors (McLyman ch.4, verbatim):**
- Window utilization `Ku = S1·S2·S3·S4 ≈ 0.4` (S1 wire insulation, S2 wire-lay fill, S3 effective
  window, S4 inter-layer/winding insulation; +workmanship folded into 0.4).
- Turns/layer = winding breadth / **insulated** wire OD (Table 4-10, *not* bare copper), derated by a
  wire-lay factor **0.85–0.90** (Table 4-2; random-wound 0.75–0.90 Table 4-3) — tight winding makes
  wound length exceed the wire-diameter prediction by 10–15%.
- **Inter-layer insulation** thickness by gauge (Table 4-5), added per layer boundary.
- **Axial margins** by gauge (Table 4-4), reduce usable breadth.
- **Rectangular-bobbin bowing:** real radial build exceeds ideal by 15–20% → derate **available**
  build to **0.85×** (negligible on round bobbins/tubes).
- **Orthocyclic self-crossover** ≈ **5–10% taller** than regular winding height — **scope caveat:**
  that figure is a *single-coil self-crossover as % of total build*, **not** a per-layer or
  per-transition figure and **not** the interleaved-lead squeeze. Use only as an order-of-magnitude
  sanity check.

**Is there a published algorithm for the squeeze?** **No.** No closed-form for the transition-wire
squeeze or real per-layer build-up exists. The literature gives only the qualitative phenomenology +
aggregate anchors above. The user's "0/1/2 wire-OD squeeze depending on transition wires above/below"
and the per-layer angular advance are **OM's own synthesis** — must be invented and FEA/empirically
validated.

**State of the art:** **no surveyed tool models the squeeze** — PI Expert (single scalar fill, can
exceed 100%), PowerEsim (ideal stacked + packing options), FEMMT (deterministic Square/Hexagonal
grid, abstracts leads), Ansys PExprt (ideal interleaved 2D/3D), Frenetic (constant-cross-section
rectangular section), Würth REDEXPERT (purely measured), PLECS (reluctance, no geometry). Biela/Kolar
recognize the transition **electrically** (C₀/4 vs C₀/3, incomplete-layer `k/z`, six-cap interleave)
but never convert it to a radial-build model. → OM modeling it would be **novel**, *but* this rests
on a bounded survey — a dedicated prior-art search is required before claiming novelty.

**⚠ Fabricated-source warning (caught by adversarial verify):** sub-agents invented patent citations
**US5208571** (actually an NMR field-homogeneity patent — does NOT model build squeeze) and a
**US3989200A** hexagonal-spacing attribution. **Both must be dropped**; re-source the hex-nesting
`√3/2·d` to standard packing/MathWorld and the alternating-direction fact to ludens.cl / Wikipedia.

**Recommended approach for the squeeze feature (future):** a **two-profile, angularly-aware**
extension of the ideal build — (1) keep the ideal grid baseline (insulated OD increment; √3/2·d
nesting only for the first 1–2 orthocyclic layers; margins/insulation/partial-layer-takes-full-slot;
0.85× bowing derate rectangular-only); (2) add a **localized crossover profile**: each transition/lead
wire occupies a migrating angular sector (~60° round / one face rectangular) and **deducts N wire-ODs
from the squeezed neighbor**, N = count of transition wires resting on top/bottom (the 0/1/2 model →
naturally yields layer-2 ≈2d vs layer-4 ≈1d). Reuse the existing `additionalCoordinates` field
(already consumed by StrayCapacitance/LeakageInductance/Temperature). Insertion point: where
`wind_by_*_layers()` sets each layer's `dimensions()[0]`; section aggregation (≈`Coil.cpp:1187–1197`)
already sums layer dims, so a per-layer deduction propagates. Squeeze count / sector width / angular
advance must be **explicit named parameters**, not hardcoded — and missing wire OD / margin /
insulation must **throw**, not default (OM rule). Validate via FEA + the peer-reviewed electrical
hooks (Dowell leakage, Biela/Kolar capacitance / SRF).

---

## 8. U/Z winding + connection-loss research (workflow `wf_eaa20991` — COMPLETE)

Fact-checked, sources re-fetched and verified. Full output:
`/tmp/claude-1000/.../tasks/wolz27pmy.output`.

- **Canonical naming (DTU / IEEE TPEL, "Transformer Winding Arrangements", peer-reviewed):** the four
  multilayer schemes are **A = U-type, B = Z-type**, C = Sectioned, D = Bank.
  - **U-type:** "the next layer starts where the previous layer ended." Turn numbering reverses each
    layer — bottom row `1-2-3-4`, top row `8-7-6-5`.
  - **Z-type:** every layer wound the same direction; "the next layer starts just above the starting
    point of the previous layer," which **requires a return wire**. Bottom `1-2-3-4`, top `5-6-7-8`.
  - → **Directly confirms §3:** current MKF (every layer same start edge/direction, names sequential)
    is **Z-type**. **U-type = reverse the electrical order of alternate layers** (exactly §4).
- **Capacitance tradeoff (Biela & Kolar, IEEE review of stray capacitance — verified verbatim):**
  - Standard/**U**: layer-to-layer voltage rises linearly 0 → 2× the per-layer voltage across the
    layer → equivalent **`C_Layer = C₀/3`**.
  - Flyback/**Z**: inter-layer voltage is **constant = V_Wdg/N_layer**, uniform field → **`C_Layer =
    C₀/4`** (lower). General eq.(21): `C_Layer = (C₀/3)·(V3² + V3·V4 + V4²)/V_LL²`.
  - Confirms: **Z lowers capacitance; capacitance follows purely from turn voltages/ordering** — so
    MKF's existing capacitance calc produces the right answer once turns are ordered/named per U/Z.
    No capacitance code change (matches user instruction).
- **Z's loss penalty (Biela/Kolar — verified verbatim):** "at the end of one layer, the wire of the
  winding has to cross all turns of this and the proximate layer. This leads to a larger overall
  winding length, and to local high-electric field strengths and dielectric losses." → validates the
  connection/return-wire loss work; the extra-length DC `ρL/A·I²` + Dowell-AC model is a sound MKF
  modeling recommendation (engineering extrapolation, not a paper formula).
- **McLyman (Transformer & Inductor Design Handbook, ch.17 — verified verbatim):** uses the literal
  terms **"U Type Winding"** and the **"Foldback winding technique"** (= Z), stating foldback "is
  preferred to the normal U type winding, even though it takes an extra step before starting the next
  layer" and "will also reduce the voltage gradient." Three ways to cut layer capacitance:
  (1) sandwich/interleave, (2) foldback, (3) more insulation (raises leakage). (Toroidal "progressive
  / back-wind" — *wind 5 forward, 4 back* — is a separate technique.)
- **State of the art (PI Expert online help — verified verbatim):** PI Expert exposes a user choice
  between **"C type"** and **"Z type"** windings; "Z-type windings … decrease transformer primary
  capacitance therefore decreasing AC Primary side switching losses," at the cost of "a slight
  increase in leakage inductance and reduced bobbin utilization." → **MKF exposing a U/Z choice
  matches existing production practice; it is not novel.** (PI calls the standard one "C type"; the
  DTU/McLyman literature calls it "U type" — same thing.)

**Naming verdict:** **`U` / `Z` are the standard literal names** in the peer-reviewed and handbook
literature (McLyman writes "U Type"; DTU writes "U-type / Z-type"). Literal `["U","Z"]` enum values
are well-grounded; "foldback" is McLyman's synonym for Z if a descriptive form is preferred.

---

## 9. Open decisions

**Resolved:**
- **Location/resolution** (user): `Section.windingOrder` → bobbin-winding-window `windingOrder` →
  default `U`. New fields on MAS Section + bobbin winding window (§5).
- **Default** (user): **`U`**. ⚠ This flips the current implicit **Z** behavior → existing
  multi-layer capacitance / turn-coordinate results & snapshot tests may change. Confirming the
  test-churn tolerance with the user (Q4 below).
- **Naming** (research): `U` / `Z` are the standard literal terms — strong default for the enum.

**All resolved (user answered 2026-06-01):**
1. **Enum values:** `["U","Z"]` (literal standard terms).
2. **Winder scope:** **rectangular + toroidal** (`wind_by_rectangular_turns` + `wind_by_round_turns`).
   Planar excluded (PCB copper — U/Z non-physical).
3. **PR scope:** U/Z winder **+ connection loss** — also wire `ConnectionElement.length` into
   `WindingOhmicLosses` / skin / proximity in this change.
4. **Default:** **`Z`** (no churn).

## 11. Implementation plan (this PR)

Backward-compat is guaranteed by **two independent defaults**: `windingOrder` default **Z** (= current
ordering) and the new real/ideal setting default **ideal** (= current geometry & losses). With both
defaults, every existing result and test is unchanged.

**A. MAS schema (submodule, then quicktype-regen `MAS.hpp`):**
- New enum type `windingOrder` = `["U","Z"]`.
- Optional `windingOrder` on `Section` (`MAS/schemas/magnetic/coil.json`).
- Optional `windingOrder` on the bobbin winding-window element (`MAS/schemas/magnetic/bobbin.json`).
- Regenerate `MAS.hpp` via the **exact** `quicktype` command in `MAS/scripts/check-mas-hpp.sh`
  (output to `MAS.hpp`); verify with `bash MAS/scripts/check-mas-hpp.sh`. Needs `quicktype` on PATH.
  Do **not** hand-edit generated classes (memory `feedback_use_existing_schemas`).

**B. Coil resolution helper:**
- `Coil::get_winding_order(sectionName)` → `section.windingOrder ?? windingWindow.windingOrder ?? Z`.
  Mirror `get_turns_alignment` (`Coil.cpp:787`+).

**C. Winder honoring U/Z (default Z = current behavior, so Z path is byte-identical):**
- `wind_by_rectangular_turns` (`Coil.cpp:4564`+): track layer-index-within-section; compute
  `reverse = (order==U) && (idxInSection % 2 == 1)`; apply the generic reversal (§4) in **both**
  winding-style branches.
- `wind_by_round_turns` (toroidal): same, reversal in polar/angular terms.
- Turn naming unchanged → capacitance auto-updates (no capacitance code touched).

**D. Real-vs-ideal setting (configurable on the fly, like all settings):**
- Add `_coilUseRealWindingGeometry = false` to `src/support/Settings.h` (`_coil*` block ~41–50) +
  `get_coil_use_real_winding_geometry()` / `set_...` in `Settings.h/.cpp`. **Default = ideal (false)**.
- Gates **both** the connection-space reservation (E) and the connection-loss inclusion (F), so ideal
  mode reproduces current results exactly.

**E. Real-mode connection-space reservation → filling factor:**
- When real: reserve a rectangle per connection/lead wire crossing a section (size = wire OD × span),
  count driven by the winding's `connections` + the inter-layer transitions implied by `windingOrder`
  (Z adds a return-wire crossing per layer boundary; U does not). Add that reserved area to the
  section's **used** area so `set_filling_factor(...)` (`Coil.cpp:3117/3120` rect, `:3513/3517` tor)
  reflects it. Reserved geometry params (sector width, lead OD source) are **explicit named
  parameters**, not hardcoded; missing required wire OD → **throw** (OM rule), never default.
- Ideal mode: numerator unchanged → filling factor identical to today.

**F. Connection loss to termination (gated by real setting):**
- Feed each winding's `ConnectionElement.length` (+ Z return-wire length) into the per-winding
  effective conductor length in `WindingOhmicLosses::calculate_ohmic_losses` (`WindingOhmicLosses.cpp`)
  and skin/proximity (`WindingSkinEffectLosses.cpp:209`, `WindingProximityEffectLosses.cpp:187`).
  Missing-but-required length → **throw**; absent connection = legitimately no contribution. Ideal
  mode: no connection contribution → losses identical to today.

**G. Painter — draw connections / reserved space (debug):**
- Add `paint_coil_connections(Magnetic)` (interface in `Painter.h`, impl in `PainterImpl.cpp`) with
  two-piece + toroidal variants, drawing each connection/reserved region as a rectangle via the
  existing `paint_rectangle(x,y,xdim,ydim,cssClass,...)`, with its own CSS class for styling. Works
  in real mode (actual reserved rects) and can overlay the would-be reserved space in ideal mode for
  debugging.

**H. Tests:**
- U-vs-Z turn ordering/coordinates (2-layer & PSPS) + resolution precedence (section > window > Z).
- **Ideal vs real comparison**: same magnetic wound both ways — assert real filling factor ≥ ideal
  (accounts for connection area) and real losses ≥ ideal (include connection length); assert ideal
  equals the pre-change baseline (no churn).
- Painter smoke test that connection rectangles render (headless SVG output).

**I. Build:** regen MAS, then `ninja -C <build> -j5 <targets>` (memory `feedback_build_with_ninja_j5`).

_Out of scope (follow-up): the full transition-wire **squeeze** asymmetry of §1.1–1.2 / §7 (the
0/1/2-wire-OD per-layer model). This PR reserves connection space at the **section** level, which is
the debuggable MVP that feeds filling-factor + loss comparisons._

---

## 10. Key references

| What | Where |
|---|---|
| Top-level wind dispatch | `src/constructive_models/Coil.cpp:4500` |
| Turn placement (rectangular) | `src/constructive_models/Coil.cpp:4526` (loop 4564–4767) |
| Layer build (rectangular) | `src/constructive_models/Coil.cpp:3835` |
| Section build (rectangular) | `src/constructive_models/Coil.cpp:2890` |
| `delimit_and_compact` | `src/constructive_models/Coil.cpp:6145` |
| `get_layers_by_section` | `src/constructive_models/Coil.cpp:787` |
| Coil settings flags | `src/support/Settings.h:41–50` |
| Ohmic loss / DC R | `src/physical_models/WindingOhmicLosses.cpp:16, 67–97` |
| Skin / proximity loss | `src/physical_models/WindingSkinEffectLosses.cpp:209`, `WindingProximityEffectLosses.cpp:187` |
| Stray capacitance | `src/physical_models/StrayCapacitance.h/.cpp` |
| MAS Section schema | `MAS/schemas/magnetic/coil.json` (`section` ≈369–481) |
| MAS connection length | `MAS/schemas/magnetic/coil.json:645`; `MAS/MAS.hpp:5581` |
| DesignReq strayCapacitance | `MAS/schemas/inputs/designRequirements.json:41` |
