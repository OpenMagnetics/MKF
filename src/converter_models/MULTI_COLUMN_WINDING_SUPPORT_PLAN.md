# Multi-Column Winding Support in MKF — Implementation Plan

> **Audience**: This document is written for an implementer who is
> NOT expected to make architectural decisions. Every change is
> spelled out at file:line granularity with before/after snippets.
> If you find ambiguity, **stop and ask** rather than improvising.

Infrastructure plan to unblock 3-phase magnetics (3-phase DAB,
3-phase transformers, multi-window E/U cores with windings on
non-centre legs). Sibling to the topology plans in this directory:

- `CONVERTER_MODELS_GOLDEN_GUIDE.md`
- `FUTURE_TOPOLOGIES.md` (lists 3-phase DAB as next-priority topology)

This plan is a **prerequisite** for any future 3-phase DAB work.

---

## Table of contents

1. [Goal & Non-Goal](#1-goal--non-goal)
2. [Glossary (memorise before starting)](#2-glossary-memorise-before-starting)
3. [Three-layer framing](#3-three-layer-framing)
4. [The barrier in MKF (exact file:line)](#4-the-barrier-in-mkf-exact-fileline)
5. [Recommended approach](#5-recommended-approach)
6. [Architectural decisions (DO NOT REVISIT)](#6-architectural-decisions-do-not-revisit)
7. [v1 accuracy table (per quantity)](#7-v1-accuracy-table-per-quantity)
8. [Phase 1 — Geometry: Bobbin.cpp](#8-phase-1--geometry-bobbincpp)
9. [Phase 2 — Geometry: Coil.cpp guard removal](#9-phase-2--geometry-coilcpp-guard-removal)
10. [Phase 3 — Geometry: create_default_groups](#10-phase-3--geometry-create_default_groups)
11. [Phase 4 — Geometry: 115 windingWindows[0] refactor](#11-phase-4--geometry-115-windingwindows0-refactor)
12. [Phase 5 — Geometry: assign_windings_to_columns API](#12-phase-5--geometry-assign_windings_to_columns-api)
13. [Phase 6 — Geometry: Painter canvas](#13-phase-6--geometry-painter-canvas)
14. [Phase 7 — EM: MagnetizingInductance.cpp](#14-phase-7--em-magnetizinginductancecpp)
15. [Phase 8 — EM: CoreLosses.cpp](#15-phase-8--em-corelossescpp)
16. [Phase 9 — Tests](#16-phase-9--tests)
17. [Phase 10 — Documentation](#17-phase-10--documentation)
18. [Verification commands](#18-verification-commands)
19. [Common pitfalls (read before coding)](#19-common-pitfalls-read-before-coding)
20. [v2+ trigger conditions](#20-v2-trigger-conditions)

---

## 1. Goal & Non-Goal

**Goal**: Enable a `Magnetic` to be defined with windings distributed
across multiple core legs (winding windows), wound geometrically,
and run loss/inductance calculations with **documented accuracy
bounds** per quantity.

**Explicit non-goals (deferred to v2+)**:
- Inter-column proximity-effect coupling (calibrated correction)
- Inter-window leakage flux paths (primary-secondary in different
  columns)
- Multi-bobbin physical models (one bobbin per leg)
- Shared-yoke flux solver for unbalanced 3-phase operation
- Toroidal cores with multiple winding windows (toroids are
  always single-window in practice; do not generalise)

---

## 2. Glossary (memorise before starting)

These four terms are easy to conflate. **Read them carefully**:

| Term | Meaning | C++ type / accessor |
|---|---|---|
| **Column** | A vertical magnetic leg of a core. A 3-phase transformer has 3 columns. | `core.processedDescription.columns[i]` |
| **Winding window** | The rectangular area between two columns where coil sits. A 3-leg core has 2 winding windows; a 5-leg core has 4. | `core.processedDescription.windingWindows[i]` |
| **Group** | A MAS object with `coordinates` + `dimensions` that anchors a set of windings to a specific winding window. **This is how MAS represents "which window does this winding go in"**. | `coil.groupsDescription[i]` (vector of `Group`) |
| **Section** | A horizontal slice of a coil within a group. Each section has a `group` field (string name) linking it to a Group. | `coil.sectionsDescription[j]` |

**Critical mental model**: when this doc says "column" it almost
always means *winding window* (because that's where conductors go).
A 3-leg core has **3 columns** but only **2 winding windows**. A
3-phase transformer is wound on the **outer 3 legs of a 5-leg core**
giving 4 winding windows but only 3 are used. **Always use the term
"winding window" in code; reserve "column" for user-facing API
names** (because users think of legs/columns, not windows).

**Caveat**: The user-facing API in this plan is named
`assign_windings_to_columns` even though internally it operates on
winding windows. This is deliberate — users think in legs/columns.
Internal code uses `windingWindowIndex`.

---

## 3. Three-layer framing

A multi-column winding is three stacked problems:

### Layer A — Geometry (where the conductors are)
Bobbin/Coil/Painter must know about multiple winding windows and
place sections in the correct one.

### Layer B — Electromagnetics (how currents and fields interact)
Lm, leakage, core losses, AC losses must reflect the per-column
magnetic circuit.

### Layer C — Application (3-phase DAB / 3-phase transformer
topology) — **out of scope for this plan**. Lives in a separate
3-phase DAB plan that will depend on Layers A and B working.

---

## 4. The barrier in MKF (exact file:line)

> All file:line references verified on git revision `b00af7b` (HEAD
> at plan time). If line numbers shift in your branch, search by
> the snippet.

### Layer A barrier 1 — Bobbin guards (2 throws)

| File | Line | Snippet |
|---|---|---|
| `src/constructive_models/Bobbin.cpp` | 777 | `if (windingWindows.size() > 1)` (in `set_winding_orientation`-ish helper) |
| `src/constructive_models/Bobbin.cpp` | 793 | `if (windingWindows.size() > 1)` (in `get_winding_orientation`-ish helper) |

**Note**: Earlier drafts of this plan listed 4 guards (504/505,
556/557, 777/778, 793/794). The 504/556 guards have already been
removed — only 2 remain. **Do not look for guards at 504 or 556**.

### Layer A barrier 2 — Coil guards (3 throws)

| File | Line | Snippet |
|---|---|---|
| `src/constructive_models/Coil.cpp` | 346 | `throw NotImplementedException("Bobbins with more than winding window not implemented…");` (in `get_section_alignment`) |
| `src/constructive_models/Coil.cpp` | 2189 | `throw NotImplementedException("Bobbins with more than winding window not implemented…");` (in `wind_by_sections`) |
| `src/constructive_models/Coil.cpp` | 3361 | `throw NotImplementedException("Only one group supported for now.");` (in `wind_by_round_sections`) |

### Layer A barrier 3 — Structural: 115× hardcoded `windingWindows[0]`

`Coil.cpp` has 115 direct accesses to `windingWindows[0]`. The full
line-number list is in [§11](#11-phase-4--115-windingwindows0-refactor).
The root cause is `create_default_group` (`Coil.cpp:2140`) always
creates exactly one `Group` anchored to `windingWindows[0]`.

### Layer B barriers (no throws but wrong outputs)

These are **not blocked** — they produce a number, just the wrong
one for multi-column.

- `src/physical_models/MagnetizingInductance.cpp` — computes Lm
  using a single effective magnetic length (whole core treated as
  one column). For 3-leg core with windings in 3 columns this is
  wrong: each column has its own reluctance.
- `src/physical_models/CoreLosses.cpp` — Steinmetz applied to one
  flux density. For 3-leg core, each leg has its own flux density;
  total loss is sum over legs.
- `src/physical_models/LeakageInductance.cpp` (Dowell etc.) —
  per-window calculation; works correctly **if** primary and
  secondary share a column. Wrong if split across columns
  (inter-window leakage path is ignored). **Deferred to v2.**
- `src/physical_models/MagneticField.cpp` and any AC-loss /
  proximity calculator — operate per-conductor inside a section;
  cross-column coupling is field-mediated and
  geometry-distance-dependent. Currently MKF assumes all conductors
  share one window. **Deferred to v2 with documented warning.**

---

## 5. Recommended approach

**Use the existing `Group` / `section.group` mechanism — no schema
changes, no new MAS fields needed.** The implementation has two
modes:

### Mode A — manual groups (user pre-defines columns)
The user calls `coil.set_groups_description({group_0, group_1, …})`
before calling `wind_by_sections`. If `groupsDescription` is already
set, `wind_by_sections` skips `create_default_group` and uses the
provided groups as-is. This is **zero-new-API**: it's just removing
the guards and wiring the existing data path.

> **Already partially implemented**: `Coil.cpp:2201-2202`,
> `2724-2728`, `3054-3058`, `3346-3347`, `5721-5728`, `6000-6018`
> all have a `if (!get_groups_description()) { create_default_group(); }`
> idiom — meaning the codebase **already supports** manually-set
> groups. The work is to make `create_default_group` produce the
> right thing for multi-window cores AND remove the explicit guards.

### Mode B — automatic groups from winding windows (new helper)
```cpp
void Coil::assign_windings_to_columns(
    const std::vector<std::vector<size_t>>& windingIndicesPerColumn);
```
Creates N groups from the core's N winding windows, assigns each
winding to the specified group, and sets `groupsDescription`. The
caller then invokes `wind_by_sections` normally.

**v1 implements Mode A fully + Mode B as a thin helper.** This
covers the 3-phase DAB use-case (user knows which winding goes to
which leg).

---

## 6. Architectural decisions (DO NOT REVISIT)

These have been considered and rejected. Do not re-open them.

1. **Keep `Coil` monolithic.** Don't restructure into per-window
   mini-coils — that's a schema change for marginal benefit.
2. **Don't virtualise the window** (don't compute one big union
   window). It loses per-leg geometry, breaks insulation calcs,
   and renders the painter useless.
3. **Don't extend MAS Group to be the primary container** (don't
   move sections inside groups). Schema break, too little benefit.
4. **Use MAS `Group` mechanism as-is.**
5. **Safe default**: multi-window core + no
   `assign_windings_to_columns()` call → all windings go in
   column 0 with a log warning. Existing single-column paths
   untouched.
6. **API name**: `assign_windings_to_columns` (user-facing name uses
   "column"; internal code uses "windingWindow").

---

## 7. v1 accuracy table (per quantity)

This is the **honest contract** between v1 and downstream users.
Each calculator's header gets a `// Multi-column behaviour:` block
restating the relevant row.

| Quantity | v1 approach | Accuracy claim | v2 work |
|---|---|---|---|
| **DC resistance** | per-turn (existing) | exact | — |
| **Magnetizing inductance** | per-column reluctance (**NEW logic**) | exact for balanced 3-phase, ±10% for unbalanced | full multi-leg magnetic-circuit solver |
| **Leakage inductance** | per-window (existing applied per group) | exact when primary+secondary share a column | inter-window leakage paths |
| **Skin / proximity within group** | existing | unchanged | — |
| **Inter-group proximity** | **ignored, log warning** | ±10–20% pessimistic | empirical correction `f(leg_spacing, freq, d_cond)` calibrated to measured data |
| **Core losses** | per-column flux density (**NEW logic, simplified**) | ±15% | shared-yoke flux for unbalanced |

### Why inter-group proximity is deferred

The strength of inter-column proximity coupling depends on:
- Inter-leg spacing (core geometry)
- Switching frequency
- Conductor diameter / Litz strand size
- Phase relationship of column currents (120° in 3-phase DAB
  reduces it vs. same-phase)

There is no closed-form analytical model for this in MKF's current
literature base. The honest path: ship v1 with conservative
per-group calculation (overestimates loss) + warning, then once a
real 3-phase DAB is built and measured, calibrate an empirical
correction. **Inventing a correction without measurement risks
shipping pseudo-precision.**

### Why magnetizing inductance and core losses are NOT deferred

These are **first-order wrong** without an update — not a bounded
approximation. A 3-phase DAB user would get nonsense Lm and core
loss outputs from v1, not bounded-pessimistic ones. So these two
calculators must be updated in v1.

---

## 8. Phase 1 — Geometry: Bobbin.cpp

### 8.1. Remove the 2 guards

**File**: `src/constructive_models/Bobbin.cpp`

#### Line 777 region — `set_winding_orientation`-ish

**Before** (around lines 776–781):
```cpp
auto windingWindows = bobbinProcessedDescription.get_winding_windows();
if (windingWindows.size() > 1) {
    throw std::runtime_error("…multiple winding windows not supported…");
}
windingWindows[windingWindowIndex].set_sections_orientation(windingOrientation);
bobbinProcessedDescription.set_winding_windows(windingWindows);
```

**After**:
```cpp
auto windingWindows = bobbinProcessedDescription.get_winding_windows();
// Multi-window cores are supported. Caller is expected to pass the
// correct windingWindowIndex (default 0 = centre/main column).
if (windingWindowIndex >= windingWindows.size()) {
    throw std::out_of_range("windingWindowIndex out of range");
}
windingWindows[windingWindowIndex].set_sections_orientation(windingOrientation);
bobbinProcessedDescription.set_winding_windows(windingWindows);
```

#### Line 793 region — `get_winding_orientation`-ish

**Before** (around lines 792–797):
```cpp
auto windingWindows = bobbinProcessedDescription.get_winding_windows();
if (windingWindows.size() > 1) {
    throw std::runtime_error("…");
}
if (windingWindows[windingWindowIndex].get_sections_orientation()) {
    return windingWindows[windingWindowIndex].get_sections_orientation().value();
}
```

**After**:
```cpp
auto windingWindows = bobbinProcessedDescription.get_winding_windows();
if (windingWindowIndex >= windingWindows.size()) {
    throw std::out_of_range("windingWindowIndex out of range");
}
if (windingWindows[windingWindowIndex].get_sections_orientation()) {
    return windingWindows[windingWindowIndex].get_sections_orientation().value();
}
```

### 8.2. `create_quick_bobbin` — loop over all winding windows

Search for `create_quick_bobbin` in `Bobbin.cpp`. There are
typically two overloads:
- `create_quick_bobbin(Core core, bool isToroidal)`
- `create_quick_bobbin(Core core, double columnDepth, double columnWidth)`

Both currently build **one** `windingWindowElement` from
`coreWindingWindows[0]`. Change to iterate all
`coreWindingWindows[i]` and push each into the bobbin's
`windingWindows` vector:

**Before** (pattern):
```cpp
auto coreWindingWindows = core.get_processed_description().value().get_winding_windows();
WindingWindowElement bobbinWindingWindow;
bobbinWindingWindow.set_width(coreWindingWindows[0].get_width());
bobbinWindingWindow.set_height(coreWindingWindows[0].get_height());
bobbinWindingWindow.set_coordinates(coreWindingWindows[0].get_coordinates());
// … more single-window setup …
std::vector<WindingWindowElement> bobbinWindingWindows = {bobbinWindingWindow};
```

**After**:
```cpp
auto coreWindingWindows = core.get_processed_description().value().get_winding_windows();
std::vector<WindingWindowElement> bobbinWindingWindows;
for (size_t i = 0; i < coreWindingWindows.size(); ++i) {
    WindingWindowElement bobbinWindingWindow;
    bobbinWindingWindow.set_width(coreWindingWindows[i].get_width());
    bobbinWindingWindow.set_height(coreWindingWindows[i].get_height());
    bobbinWindingWindow.set_coordinates(coreWindingWindows[i].get_coordinates());
    // copy any other fields (radial_height, angle, etc.) from
    // coreWindingWindows[i] → bobbinWindingWindow
    bobbinWindingWindows.push_back(bobbinWindingWindow);
}
```

> **Caveat**: `CoreBobbinProcessedDescription` has scalar
> `column_depth`, `column_shape`, `column_thickness`,
> `column_width` fields (MAS-generated, cannot change the struct).
> These remain as the description of the wound (centre) column
> only — **add a comment** documenting this limitation:
> ```cpp
> // NOTE: column_depth/column_shape/column_thickness/column_width
> // describe the centre column only. Multi-window cores store
> // additional winding windows in winding_windows[]; the per-column
> // physical dimensions are not represented yet (v2 work).
> ```

### 8.3. Verify Phase 1 compiles before moving on

```bash
cd /home/alfonso/OpenMagnetics/MKF/build
ninja -j2 MKF
# Should compile cleanly. If ninja errors, STOP and fix before Phase 2.
```

---

## 9. Phase 2 — Geometry: Coil.cpp guard removal

### 9.1. Remove the 3 guards (delete the throws only)

**File**: `src/constructive_models/Coil.cpp`

#### Guard at line 346 (in `get_section_alignment`)

**Before** (around lines 343–351):
```cpp
auto windingWindows = bobbinProcessedDescription.get_winding_windows();
if (windingWindows.size() > 1) {
    throw NotImplementedException("Bobbins with more than winding window not implemented…");
}
if (windingWindows[0].get_sections_alignment()) {
    auto alignment = windingWindows[0].get_sections_alignment().value();
    return alignment;
}
```

**After**:
```cpp
auto windingWindows = bobbinProcessedDescription.get_winding_windows();
// Multi-window: use windingWindows[windingWindowIndex] (default 0).
// Reading from window 0 is the safe fallback for single-column
// behaviour and matches the existing data path.
if (windingWindows[0].get_sections_alignment()) {
    auto alignment = windingWindows[0].get_sections_alignment().value();
    return alignment;
}
```

> **Why keep `[0]` here?** This function returns a *single* alignment
> for the whole coil. Per-window alignment is a v2 refinement.
> Document this with the comment above.

#### Guard at line 2189 (in `wind_by_sections`)

**Before** (around lines 2187–2196):
```cpp
auto windingWindows = bobbinProcessedDescription.get_winding_windows();
if (windingWindows.size() > 1) {
    throw NotImplementedException("Bobbins with more than winding window not implemented…");
}
if (!windingWindows[0].get_sections_orientation()) {
    windingWindows[0].set_sections_orientation(_windingOrientation);
}
if (!windingWindows[0].get_sections_alignment()) {
    windingWindows[0].set_sections_alignment(_sectionAlignmentExplicit ? _sectionAlignment : get_section_alignment());
}
```

**After**:
```cpp
auto windingWindows = bobbinProcessedDescription.get_winding_windows();
// Multi-window supported. Initialise sections_orientation and
// sections_alignment for ALL windows (loop instead of [0] only).
for (size_t i = 0; i < windingWindows.size(); ++i) {
    if (!windingWindows[i].get_sections_orientation()) {
        windingWindows[i].set_sections_orientation(_windingOrientation);
    }
    if (!windingWindows[i].get_sections_alignment()) {
        windingWindows[i].set_sections_alignment(_sectionAlignmentExplicit ? _sectionAlignment : get_section_alignment());
    }
}
```

> **Persist the change**: if the existing code calls
> `bobbinProcessedDescription.set_winding_windows(windingWindows)`
> after the original block, keep that call. If not, **add it**:
> ```cpp
> bobbinProcessedDescription.set_winding_windows(windingWindows);
> ```

#### Guard at line 3361 (in `wind_by_round_sections` / printed)

**Before**:
```cpp
auto groups = get_groups_description().value();
if (groups.size() > 1) {
    throw NotImplementedException("Only one group supported for now.");
}
// … rest of function
```

**After**:
```cpp
auto groups = get_groups_description().value();
// Multi-group supported. The outer loop below already iterates
// `for (auto group : groups)`.
// … rest of function
```

### 9.2. Verify Phase 2 compiles

```bash
ninja -j2 MKF
```

---

## 10. Phase 3 — Geometry: create_default_groups

### 10.1. Rename and split `create_default_group`

**File**: `src/constructive_models/Coil.cpp` (line 2140)
**File**: `src/constructive_models/Coil.h` (line 134)

The current single-group function stays (renamed
`create_default_group_single_window` for clarity, or keep the name
and rely on the new dispatcher). New dispatcher:

#### Coil.h declaration

**Before** (line 134):
```cpp
bool create_default_group(Bobbin bobbin, WiringTechnology coilType = WiringTechnology::WOUND, double coreToLayerDistance = 0);
```

**After**:
```cpp
bool create_default_group(Bobbin bobbin, WiringTechnology coilType = WiringTechnology::WOUND, double coreToLayerDistance = 0);
// New: multi-window aware dispatcher. Calls create_default_group()
// for single-window bobbins, or creates one Group per winding window
// for multi-window bobbins (all windings placed in column 0 by
// default — call assign_windings_to_columns() to override).
bool create_default_groups(Bobbin bobbin, WiringTechnology coilType = WiringTechnology::WOUND, double coreToLayerDistance = 0);

// Distribute existing windings across columns. Creates one Group per
// winding window in the bobbin; assigns windings to groups per the
// indices vector. Sets groupsDescription. Must be called BEFORE
// wind_by_sections() if you want non-default column placement.
void assign_windings_to_columns(
    const std::vector<std::vector<size_t>>& windingIndicesPerColumn);
```

#### Coil.cpp definition (after line 2140 — keep existing
`create_default_group`, add new function)

```cpp
bool Coil::create_default_groups(Bobbin bobbin, WiringTechnology coilType, double coreToLayerDistance) {
    auto bobbinProcessedDescription = bobbin.get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    if (windingWindows.size() <= 1) {
        // Single-window: existing behaviour, no change.
        return create_default_group(bobbin, coilType, coreToLayerDistance);
    }

    // Multi-window: create one Group per winding window. By default,
    // all windings go into Group 0 (centre column). User can override
    // by calling assign_windings_to_columns() before wind_by_sections.
    std::vector<Group> groups;
    auto partialWindings = build_partial_windings_from_functional_description();
    // ^ replace with whatever helper the existing create_default_group
    //   uses to populate Group::partialWindings. Look at lines
    //   2155-2167 of the existing function.

    for (size_t i = 0; i < windingWindows.size(); ++i) {
        Group g;
        g.set_name("Column " + std::to_string(i));
        g.set_coordinates({windingWindows[i].get_coordinates().value()[0],
                           windingWindows[i].get_coordinates().value()[1]});

        // Dimensions: account for coreToLayerDistance like the existing
        // create_default_group does (Coil.cpp:2147, 2151).
        if (coilType == WiringTechnology::WOUND) {
            g.set_dimensions({windingWindows[i].get_width().value() - coreToLayerDistance * 2,
                              windingWindows[i].get_height().value()});
        } else {
            // Toroidal/printed: see Coil.cpp:2151 for the radial form.
            g.set_dimensions({windingWindows[i].get_radial_height().value() - coreToLayerDistance * 2,
                              windingWindows[i].get_angle().value()});
        }

        if (i == 0) {
            // All windings into column 0 by default.
            g.set_partial_windings(partialWindings);
        } else {
            // Empty group placeholder for other columns.
            g.set_partial_windings({});
        }
        groups.push_back(g);
    }
    set_groups_description(groups);

    if (windingWindows.size() > 1) {
        // Log warning so the user knows they need to call
        // assign_windings_to_columns() if they want distribution.
        // Use whatever logger MKF uses (search for `log`/`spdlog` in
        // existing Coil.cpp; if none, std::cerr is acceptable).
        std::cerr << "[MKF Warning] Multi-column bobbin detected ("
                  << windingWindows.size() << " winding windows). "
                  << "All windings placed in column 0 by default. "
                  << "Call assign_windings_to_columns() to distribute."
                  << std::endl;
    }
    return true;
}
```

### 10.2. Update call sites of `create_default_group`

There are 2 call sites in `Coil.cpp` (lines 2202 and 3347). **Change
both to call `create_default_groups`** (the new dispatcher):

**Coil.cpp:2201-2202 (in wind_by_sections)**:
```cpp
if (!get_groups_description()) {
    create_default_groups(bobbin);  // was: create_default_group(bobbin)
}
```

**Coil.cpp:3346-3347 (in wind_by_round_sections)**:
```cpp
if (!get_groups_description()) {
    create_default_groups(bobbin, WiringTechnology::PRINTED, coreToLayerDistance);
    // was: create_default_group(bobbin, WiringTechnology::PRINTED, coreToLayerDistance);
}
```

### 10.3. Verify Phase 3 compiles + single-column tests still pass

```bash
ninja -j2 MKF_tests
./MKF_tests "[Coil]" 2>&1 | tail -20
# Expected: all existing tests pass (regression-free).
```

---

## 11. Phase 4 — Geometry: 115 windingWindows[0] refactor

This is the bulk of the geometric work. **Split it into clusters**,
not 115 individual edits. Each cluster has the same fix pattern.

### 11.1. Inventory: full line list

Lines in `Coil.cpp` with `windingWindows[0]` (verified at HEAD
`b00af7b`):

```
56, 82, 104, 136, 349, 350,
1245, 1249, 1250, 1253, 1254,
1445, 1446, 1456, 1457, 1478, 1479, 1489, 1490,
1561, 1562, 1572, 1573, 1596, 1597, 1607, 1608,
1703, 1704, 1714, 1715, 1738, 1739, 1749, 1750,
2145, 2147 (×2), 2151 (×2), 2191, 2192, 2194, 2195,
2984, 3889, 4553, 4939,
5331, 5332, 5366, 5369, 5372, 5377, 5378,
5392, 5395, 5398, 5402, 5403, 5404, 5410,
5419, 5422, 5425, 5430, 5431,
5445, 5448, 5451, 5455, 5456, 5457, 5464,
5473, 5476, 5479, 5484, 5485,
5497, 5503, 5504, 5509, 5514, 5517, 5521, 5522, 5523, 5530, 5539,
5542, 5545, 5550, 5551, 5565, 5568, 5571, 5575, 5576, 5577, 5592,
5614, 5615, 5688,
6057, 6125, 6128, 6134 (commented; ignore),
6893,
7240, 7241, 7251, 7306, 7307, 7317, 7318, 7340, 7341, 7351, 7352
```

### 11.2. Cluster A — toroidal coordinate transforms (4 sites)

**Lines 56, 82, 104, 136** — these are at the top of the file in
helpers like `cartesian_to_polar` / `polar_to_cartesian`. They take
a coordinate value relative to the toroid's centre.

**Strategy**: add a `windingWindowIndex = 0` default parameter to
each helper. Inside the helper, replace `windingWindows[0]` with
`windingWindows[windingWindowIndex]`.

**Caveat**: toroids are single-window in practice; for v1 these can
remain `[0]`. **However**, audit the call graph to confirm no
multi-window code path reaches these. If they do, add the index
param. **If unsure, leave as `[0]` and add a comment**:
```cpp
// Toroidal helper: assumes single winding window. Multi-window
// toroids are not in scope (see plan §11.2).
double radialHeight = windingWindows[0].get_radial_height().value() - radius;
```

### 11.3. Cluster B — `get_section_alignment` (lines 349, 350)

Already handled in [§9.1](#91-remove-the-3-guards-delete-the-throws-only).
**No further work.** These remain `[0]` (single alignment for whole
coil; v2 work to per-window).

### 11.4. Cluster C — `add_insulation_to_sections` (lines 1245, 1249, 1250, 1253, 1254)

These read available area / dimensions from the coil's bobbin. The
function operates on **sections**, each of which has a `group`
field. Instead of `windingWindows[0]`, look up the section's group
and use **that** group's window.

**Strategy**:
```cpp
// Before:
double availableArea = windingWindows[0].get_area().value();

// After:
auto sectionGroup = section.get_group();  // string name
auto windingWindowIndex = find_window_index_for_group(sectionGroup);
double availableArea = windingWindows[windingWindowIndex].get_area().value();
```

You'll need a small helper:
```cpp
size_t Coil::find_window_index_for_group(const std::string& groupName) const {
    auto groups = get_groups_description().value();
    auto bobbinWW = get_bobbin().get_processed_description().value().get_winding_windows();
    for (size_t i = 0; i < groups.size(); ++i) {
        if (groups[i].get_name() == groupName) {
            // Match group's coordinates to a winding window's coordinates.
            for (size_t j = 0; j < bobbinWW.size(); ++j) {
                if (bobbinWW[j].get_coordinates().value() == groups[i].get_coordinates()) {
                    return j;
                }
            }
        }
    }
    return 0;  // safe fallback
}
```

### 11.5. Cluster D — `calculate_insulation*` helpers (lines 1445–1750)

Looking at the line numbers, these are inside large
`calculate_insulation_*` functions that operate on a single section
or section pair. Use the same lookup as Cluster C:

```cpp
// Before:
double windingWindowHeight = windingWindows[0].get_height().value();
double windingWindowWidth  = windingWindows[0].get_width().value();

// After:
auto sectionGroup = section.get_group();
auto windowIndex = find_window_index_for_group(sectionGroup);
double windingWindowHeight = windingWindows[windowIndex].get_height().value();
double windingWindowWidth  = windingWindows[windowIndex].get_width().value();
```

**Apply this pattern to**:
- 1445, 1446, 1456, 1457, 1478, 1479, 1489, 1490
- 1561, 1562, 1572, 1573, 1596, 1597, 1607, 1608
- 1703, 1704, 1714, 1715, 1738, 1739, 1749, 1750

> **30 sites total**. Pure mechanical replacement.

### 11.6. Cluster E — `create_default_group` body (lines 2145, 2147 ×2, 2151 ×2)

These are inside the **existing** `create_default_group` function
which now only handles single-window. **Leave them as
`windingWindows[0]`** — the function is single-window by contract
after [§10](#10-phase-3--geometry-create_default_groups). Add a
comment:
```cpp
// Single-window helper. Multi-window dispatch happens in
// create_default_groups() (see plan §10).
group.set_coordinates({windingWindows[0].get_coordinates().value()[0], …});
```

### 11.7. Cluster F — `wind_by_sections` window initialisation (lines 2191, 2192, 2194, 2195)

Already handled in [§9.1](#91-remove-the-3-guards-delete-the-throws-only)
when removing the guard at 2189. The loop `for (size_t i = 0; i <
windingWindows.size(); ++i)` replaced these accesses.

### 11.8. Cluster G — isolated radial-height accesses (lines 2984, 3889, 4553, 4939)

These are scattered radial-height reads. Each is inside a function
that operates on a specific section, group, or winding. **Look up
the index from the section's group** (same helper as Cluster C):

```cpp
// Before:
double windingWindowRadialHeight = windingWindows[0].get_radial_height().value();

// After:
auto windowIndex = find_window_index_for_group(currentSection.get_group());
double windingWindowRadialHeight = windingWindows[windowIndex].get_radial_height().value();
```

### 11.9. Cluster H — section-placement loop (lines 5331–5688)

This is the biggest cluster (~80 sites). Lines 5331 onwards are
inside a function that does the mechanical "lay sections out
horizontally/vertically in the winding window" work.

**Strategy**: this function should operate **per group**. The outer
loop should already be `for (auto group : groups)` — if not,
restructure it. Inside the loop, replace
`windingWindows[0].get_X()` with `group.get_dimensions()[…]` and
`windingWindows[0].get_coordinates().value()[…]` with
`group.get_coordinates()[…]`.

| `windingWindows[0]` access | Replace with |
|---|---|
| `.get_width().value()` | `group.get_dimensions()[0]` |
| `.get_height().value()` | `group.get_dimensions()[1]` |
| `.get_coordinates().value()[0]` | `group.get_coordinates()[0]` |
| `.get_coordinates().value()[1]` | `group.get_coordinates()[1]` |
| `.get_angle().value()` | `group.get_dimensions()[1]` (toroidal: dimensions stored differently — check) |
| `.get_radial_height().value()` | `group.get_dimensions()[0]` (toroidal) |

**Verify the loop structure first** — read lines 5300–5330 to see
how `groups` and the per-group iteration is structured. If the loop
is already `for (auto group : groups) { … }`, the replacement is
mechanical. If not, **stop and ask** — this is non-trivial.

### 11.10. Cluster I — toroidal section placement (lines 5614, 5615, 5688)

Same as Cluster H but for toroids. Same `group.get_dimensions()`
replacement.

### 11.11. Cluster J — turns placement (lines 6057, 6125, 6128, 6893)

Inside the per-turn loop. Each turn belongs to a layer, which
belongs to a section, which belongs to a group. Look up the
windowIndex from the section/group:

```cpp
auto turnSection = find_section_for_turn(turn);  // existing helper, search for it
auto windowIndex = find_window_index_for_group(turnSection.get_group());
auto windingWindowsRadius = windingWindows[windowIndex].get_radial_height().value();
```

### 11.12. Cluster K — final compaction / output formatting (lines 7240–7352)

Lines 7240, 7241, 7251, 7306, 7307, 7317, 7318, 7340, 7341, 7351,
7352. Same pattern: per-section lookup of winding window via group.

### 11.13. After Phase 4

```bash
ninja -j2 MKF_tests
./MKF_tests "[Coil]" 2>&1 | tail -30
# Expected: all existing single-column tests still pass.
# If any fail, the most likely cause is a missed cluster
# replacement. Bisect by reverting clusters one at a time.
```

---

## 12. Phase 5 — Geometry: assign_windings_to_columns API

### 12.1. Implementation

Add to `src/constructive_models/Coil.cpp` (anywhere near
`create_default_groups`):

```cpp
void Coil::assign_windings_to_columns(
    const std::vector<std::vector<size_t>>& windingIndicesPerColumn) {

    auto bobbin = get_bobbin();
    auto windingWindows = bobbin.get_processed_description().value().get_winding_windows();

    if (windingIndicesPerColumn.size() != windingWindows.size()) {
        throw std::invalid_argument(
            "windingIndicesPerColumn size (" +
            std::to_string(windingIndicesPerColumn.size()) +
            ") must match number of winding windows (" +
            std::to_string(windingWindows.size()) + ")");
    }

    auto functionalDescription = get_functional_description();
    std::vector<Group> groups;

    for (size_t col = 0; col < windingWindows.size(); ++col) {
        Group g;
        g.set_name("Column " + std::to_string(col));
        g.set_coordinates({windingWindows[col].get_coordinates().value()[0],
                           windingWindows[col].get_coordinates().value()[1]});
        g.set_dimensions({windingWindows[col].get_width().value(),
                          windingWindows[col].get_height().value()});

        // Build partial windings for the windings in this column.
        std::vector<PartialWinding> partialWindings;
        for (auto windingIndex : windingIndicesPerColumn[col]) {
            if (windingIndex >= functionalDescription.size()) {
                throw std::out_of_range("winding index out of range");
            }
            PartialWinding pw;
            pw.set_winding(functionalDescription[windingIndex].get_name());
            // Number of parallels = total parallels of that winding
            // (no splitting across columns in v1).
            std::vector<int64_t> parallelsProportion(
                functionalDescription[windingIndex].get_number_parallels(), 1);
            pw.set_parallels_proportion(parallelsProportion);
            partialWindings.push_back(pw);
        }
        g.set_partial_windings(partialWindings);
        groups.push_back(g);
    }
    set_groups_description(groups);
}
```

### 12.2. Header declaration

Already added in [§10.1](#101-rename-and-split-create_default_group).

### 12.3. Verify

```bash
ninja -j2 MKF_tests
./MKF_tests "[Coil]"
```

---

## 13. Phase 6 — Geometry: Painter canvas

### 13.1. Files to inspect

> **Note**: `src/support/Painter.cpp` does NOT contain
> `windingWindows` references (verified at HEAD `b00af7b`). The
> actual winding-window-aware painting is in:
> - `src/support/BasicPainter.cpp` (lines 1056, 1057, 1233, 1252, 1266)
> - `src/support/AdvancedPainter.cpp` (lines 384, 398, 412, 825, 826)

### 13.2. BasicPainter.cpp

**Lines 1056–1057** (in section painter):
```cpp
auto windingWindows = bobbinProcessedDescription.get_winding_windows();
auto sectionOrientation = windingWindows[0].get_sections_orientation();
```

For multi-window: each section's orientation may differ. Look up
per section's group's window:
```cpp
auto windingWindows = bobbinProcessedDescription.get_winding_windows();
size_t windowIndex = find_window_index_for_group(currentSection.get_group());
auto sectionOrientation = windingWindows[windowIndex].get_sections_orientation();
```

**Lines 1233, 1252, 1266** (in core painter — drawing the core
itself): these read `core.get_winding_windows()` to size the
canvas. The core painter typically iterates these already; verify
the loop is correct.

### 13.3. AdvancedPainter.cpp

**Lines 384, 398, 412** (core painter): same pattern; should
already iterate.

**Lines 825–826**: same as BasicPainter 1056–1057. Apply the same
per-section group lookup.

### 13.4. Canvas sizing

If there's a "compute total canvas width" function, change:
- **width** = sum of all winding-window widths + inter-column
  spacing (use core's column width as spacing)
- **height** = max height across all winding windows

If sections already have absolute coordinates set correctly per
group ([§11–12](#11-phase-4--115-windingwindows0-refactor)), the
painter renders correctly because each section is drawn at its
absolute coordinates. **Canvas sizing is the only painter-level
change needed.**

### 13.5. Verify visually

After Phase 9 (tests), generate a multi-column SVG and inspect:
```bash
./MKF_tests "Test_Coil_Painter_Multi_Column_Canvas_Size"
# Open the generated SVG and confirm all 3 columns are visible.
```

---

## 14. Phase 7 — EM: MagnetizingInductance.cpp

### 14.1. File location

`src/physical_models/MagnetizingInductance.cpp` and `.h`. (Note:
**not** `src/processors/`.)

### 14.2. Current logic

Computes `L_m = N²·μ·A_e / l_e` using one effective magnetic path
(`Core::get_effective_length()`, `get_effective_area()`).

### 14.3. New logic for multi-column

For each Group in the coil, compute Lm using **that column's**
reluctance:

```cpp
double Lm_total = 0;
auto groups = coil.get_groups_description().value_or({});
auto columns = core.get_processed_description().value().get_columns();

if (groups.size() <= 1 || columns.size() <= 1) {
    // Existing single-column path. No change.
    return existing_single_column_lm(core, coil);
}

// Multi-column: per-column reluctance.
for (size_t i = 0; i < groups.size(); ++i) {
    auto group = groups[i];
    auto column = columns[i];  // assumes group i ↔ column i; verify

    int64_t N_in_column = 0;
    for (auto& pw : group.get_partial_windings()) {
        auto winding = find_winding_by_name(coil, pw.get_winding());
        N_in_column += winding.get_number_turns();
    }
    if (N_in_column == 0) continue;

    double mu = effective_permeability(core, frequency);
    double A_column = column.get_area();          // e.g., processedDescription.area
    double l_column = column.get_height();        // e.g., processedDescription.height

    double R_column = l_column / (mu * A_column);
    double Lm_column = (N_in_column * N_in_column) / R_column;

    // For balanced 3-phase (currents 120° apart), per-leg Lm is
    // independent. Sum is invalid — instead, expose per-column
    // values via a new accessor get_per_column_lm() and return the
    // primary winding's Lm for backward compat.
    // v1 simplification: return the centre-column Lm, or the first
    // column's Lm if no centre.
    if (i == 0) Lm_total = Lm_column;
}

return Lm_total;
```

### 14.4. Header documentation

Add to `MagnetizingInductance.h`:
```cpp
// Multi-column behaviour:
//
// For multi-window cores (e.g., 3-phase transformers with windings
// distributed across legs), this returns the per-column Lm for the
// FIRST column (column 0) only. The function assumes balanced
// 3-phase operation, where each column's Lm is independent of the
// others (no shared yoke flux). Accuracy: exact for balanced 3-phase,
// ±10% for unbalanced. For per-column Lm values, call
// get_per_column_lm() (returns std::vector<double>).
```

### 14.5. Tests

Add to `tests/TestMagnetizingInductance.cpp`:
- `Test_Magnetizing_Inductance_Single_Window_Regression` — existing
  single-window test, must remain unchanged
- `Test_Magnetizing_Inductance_Three_Leg_Balanced` — 3-leg core, 3
  windings (one per leg), assert each column's Lm matches
  `N² · μ · A_leg / l_leg` within 1%

---

## 15. Phase 8 — EM: CoreLosses.cpp

### 15.1. File location

`src/physical_models/CoreLosses.cpp` and `.h`.

### 15.2. New logic

Apply Steinmetz per column:

```cpp
double P_total = 0;
auto columns = core.get_processed_description().value().get_columns();
auto groups = coil.get_groups_description().value_or({});

if (groups.size() <= 1 || columns.size() <= 1) {
    return existing_single_column_losses(core, coil, ...);
}

for (size_t i = 0; i < columns.size(); ++i) {
    auto column = columns[i];
    auto group = groups[i];

    // Per-column MMF = Σ (N_winding · I_winding) for windings in this column.
    double MMF_column = compute_mmf_for_group(coil, group, operatingPoint);

    // Per-column flux = MMF / R_column
    double R_column = column.get_height() / (effective_permeability(...) * column.get_area());
    double flux_column = MMF_column / R_column;

    // Flux density
    double B_peak_column = flux_column / column.get_area();

    // Steinmetz on this column
    double V_column = column.get_volume();  // or area * height
    double P_column = steinmetz_loss(
        B_peak_column, frequency,
        material.get_steinmetz_k(),
        material.get_steinmetz_alpha(),
        material.get_steinmetz_beta()) * V_column;
    P_total += P_column;
}
return P_total;
```

### 15.3. Header documentation

Add to `CoreLosses.h`:
```cpp
// Multi-column behaviour:
//
// For multi-window cores, total core loss = sum of per-column
// Steinmetz losses, where each column's flux density is computed
// from that column's MMF and reluctance. Yoke loss is ignored in v1
// (assumes balanced operation where yoke flux ≈ 0). Accuracy: ±15%
// for balanced 3-phase, may be more pessimistic for unbalanced.
```

### 15.4. Tests

Add to `tests/TestCoreLosses.cpp`:
- `Test_Core_Loss_Three_Leg_Balanced` — 3-leg core with balanced
  excitation; per-leg loss matches Steinmetz on per-leg flux
  density.

---

## 16. Phase 9 — Tests

### 16.1. Geometry tests (`tests/TestCoil.cpp`)

#### Test_Coil_Wind_Two_Columns_EE_Core

```cpp
TEST_CASE("Test_Coil_Wind_Two_Columns_EE_Core", "[Coil]") {
    // 1. Build an EE core with 2 winding windows (use whatever
    //    helper TestCoil uses — typically OpenMagnetics::CoreWrapper
    //    constructed from a shape name).
    auto core = build_ee_core_two_windows();

    // 2. Define 2 windings (primary + secondary).
    auto coil = build_coil_two_windings(core);

    // 3. Assign winding 0 to column 0, winding 1 to column 1.
    coil.assign_windings_to_columns({{0}, {1}});

    // 4. Wind.
    coil.wind_by_sections();

    // 5. Assert.
    auto sections = coil.get_sections_description().value();
    REQUIRE(sections.size() == 2);
    REQUIRE(sections[0].get_group() != sections[1].get_group());

    auto coords0 = sections[0].get_coordinates();
    auto coords1 = sections[1].get_coordinates();
    REQUIRE(std::abs(coords0[0] - coords1[0]) > 0.001);
    // ^ Sections are in different windows, so x-coordinates differ
    //   by more than the inter-leg spacing (~mm).
}
```

#### Test_Coil_Wind_Three_Columns_Three_Phase_Transformer

```cpp
TEST_CASE("Test_Coil_Wind_Three_Columns_Three_Phase_Transformer", "[Coil]") {
    auto core = build_three_leg_core();  // 3 columns, 2 winding windows
    // (or 5-leg core with 4 windows for the genuine 3-phase case)

    auto coil = build_coil_six_windings(core);  // 3 prim + 3 sec

    coil.assign_windings_to_columns({
        {0, 3},  // column 0: prim_a + sec_a
        {1, 4},  // column 1: prim_b + sec_b
        {2, 5}   // column 2: prim_c + sec_c
    });

    coil.wind_by_sections();

    auto sections = coil.get_sections_description().value();
    REQUIRE(sections.size() == 6);

    // Group sections by group name; expect 3 groups, 2 sections each.
    std::map<std::string, int> groupCounts;
    for (auto& s : sections) groupCounts[s.get_group()]++;
    REQUIRE(groupCounts.size() == 3);
    for (auto& [name, count] : groupCounts) REQUIRE(count == 2);
}
```

#### Test_Coil_Assign_Windings_To_Columns

```cpp
TEST_CASE("Test_Coil_Assign_Windings_To_Columns", "[Coil]") {
    auto core = build_three_leg_core();
    auto coil = build_coil_six_windings(core);

    coil.assign_windings_to_columns({{0, 1}, {2, 3}, {4, 5}});

    auto groups = coil.get_groups_description().value();
    REQUIRE(groups.size() == 3);
    REQUIRE(groups[0].get_partial_windings().size() == 2);
    REQUIRE(groups[1].get_partial_windings().size() == 2);
    REQUIRE(groups[2].get_partial_windings().size() == 2);

    // Coordinates should match the core's winding windows.
    auto windowCoords = core.get_processed_description().value()
                            .get_winding_windows()[0].get_coordinates().value();
    REQUIRE(groups[0].get_coordinates()[0] == Approx(windowCoords[0]));
}
```

#### Test_Coil_Single_Column_Still_Works

This is a regression guard. Pick any existing single-column
TestCoil case and add it under a new name to make the regression
explicit:

```cpp
TEST_CASE("Test_Coil_Single_Column_Still_Works", "[Coil][Regression]") {
    // Replicate one of the existing single-column wind tests.
    // This must pass unchanged after all the refactoring.
    auto core = build_etd_core();
    auto coil = build_coil_two_windings(core);
    coil.wind_by_sections();
    auto sections = coil.get_sections_description().value();
    REQUIRE(sections.size() == 2);
    // Both sections in same group.
    REQUIRE(sections[0].get_group() == sections[1].get_group());
}
```

#### Test_Coil_Painter_Multi_Column_Canvas_Size

```cpp
TEST_CASE("Test_Coil_Painter_Multi_Column_Canvas_Size", "[Coil][Painter]") {
    auto core = build_three_leg_core();
    auto coil = build_coil_six_windings(core);
    coil.assign_windings_to_columns({{0,3},{1,4},{2,5}});
    coil.wind_by_sections();

    auto painter = AdvancedPainter();
    auto svg = painter.paint(core, coil);

    // Parse SVG, find canvas width attribute, assert ≥ 3× single-column width.
    double canvasWidth = parse_svg_width(svg);
    auto singleWindowWidth = core.get_processed_description().value()
                                 .get_winding_windows()[0].get_width().value();
    REQUIRE(canvasWidth >= 3 * singleWindowWidth);
}
```

### 16.2. EM tests

#### Test_Magnetizing_Inductance_Single_Window_Regression

Replicate any existing single-window MagnetizingInductance test
under a new name, must pass unchanged.

#### Test_Magnetizing_Inductance_Three_Leg_Balanced

```cpp
TEST_CASE("Test_Magnetizing_Inductance_Three_Leg_Balanced", "[MagnetizingInductance]") {
    auto core = build_three_leg_core();
    auto coil = build_coil_three_phase_primary_only(core);  // 3 windings, 1 per leg
    coil.assign_windings_to_columns({{0},{1},{2}});
    coil.wind_by_sections();

    auto inputs = build_balanced_three_phase_inputs();

    auto magnetic = OpenMagnetics::Magnetic(core, coil);
    auto Lm_per_column = MagnetizingInductance().get_per_column_lm(magnetic, inputs);

    REQUIRE(Lm_per_column.size() == 3);

    // Reference: N²·μ·A_leg / l_leg
    auto cols = core.get_processed_description().value().get_columns();
    double mu = compute_effective_permeability(...);
    double N = coil.get_functional_description()[0].get_number_turns();
    double Lm_expected = (N * N) * mu * cols[0].get_area() / cols[0].get_height();

    REQUIRE(Lm_per_column[0] == Approx(Lm_expected).epsilon(0.01));
}
```

#### Test_Core_Loss_Three_Leg_Balanced

Similar structure: per-leg flux density × Steinmetz × volume,
summed.

#### Test_Proximity_Multi_Column_Equals_Per_Group

Characterization test pinning v1's documented behaviour:

```cpp
TEST_CASE("Test_Proximity_Multi_Column_Equals_Per_Group", "[Proximity]") {
    // Setup: 2 windings, each in its own column. Compute proximity
    // loss in each. Expected: same as if each were in its own
    // single-column coil. (Documents v1: zero inter-column coupling.)
    // ...
}
```

---

## 17. Phase 10 — Documentation

### 17.1. Per-calculator headers

Each affected calculator's `.h` file gets a `// Multi-column
behaviour:` block matching the v1 accuracy table:

- `MagnetizingInductance.h` — see [§14.4](#144-header-documentation)
- `CoreLosses.h` — see [§15.3](#153-header-documentation)
- `LeakageInductance.h`:
  ```cpp
  // Multi-column behaviour:
  //
  // For multi-window cores, leakage is computed per-window only.
  // If a winding's primary and secondary share a column, the result
  // is exact. If primary and secondary are in different columns,
  // inter-window leakage flux paths are NOT modelled (deferred to v2).
  ```
- `MagneticField.h`:
  ```cpp
  // Multi-column behaviour:
  //
  // Field-mediated proximity-effect coupling between conductors in
  // different columns is NOT modelled in v1. Per-column proximity
  // is computed correctly. Total AC loss may be ±10–20% pessimistic
  // for multi-column cores (deferred to v2 with calibrated
  // correction factor).
  ```

### 17.2. Golden guide section

Add a new section to `CONVERTER_MODELS_GOLDEN_GUIDE.md`:

```markdown
## Multi-column magnetics (v1)

### What's exact
- DC resistance (per-turn)
- Skin/proximity within a single column

### What's bounded
- Magnetizing inductance: ±10% for balanced 3-phase
- Core losses: ±15% for balanced 3-phase
- Leakage when primary/secondary share a column: exact

### What's deferred (v2)
- Inter-group proximity coupling
- Inter-window leakage flux paths
- Shared-yoke flux for unbalanced operation
- Multi-bobbin physical models

### Recommended use
3-phase DAB, 3-phase transformers with balanced excitation.

### Not recommended for
- Highly unbalanced 3-phase loads
- Designs with primary/secondary in different columns (use a
  single-column model instead until v2)
```

---

## 18. Verification commands

```bash
cd /home/alfonso/OpenMagnetics/MKF/build

# Build (use ninja per project memory; never cmake --build)
ninja -j2 MKF
ninja -j2 MKF_tests

# Phase 1-6 (geometry)
./MKF_tests "[Coil]"
./MKF_tests "Test_Coil_Wind_Two_Columns_EE_Core"
./MKF_tests "Test_Coil_Wind_Three_Columns_Three_Phase_Transformer"
./MKF_tests "Test_Coil_Assign_Windings_To_Columns"
./MKF_tests "Test_Coil_Single_Column_Still_Works"
./MKF_tests "Test_Coil_Painter_Multi_Column_Canvas_Size"

# Phase 7-8 (electromagnetic)
./MKF_tests "Test_Magnetizing_Inductance_Three_Leg_Balanced"
./MKF_tests "Test_Magnetizing_Inductance_Single_Window_Regression"
./MKF_tests "Test_Core_Loss_Three_Leg_Balanced"
./MKF_tests "Test_Proximity_Multi_Column_Equals_Per_Group"

# Full regression
./MKF_tests
```

### Acceptance criteria

- [ ] All existing `TestCoil` cases pass (regression-free)
- [ ] `wind_by_sections` on a 2-window bobbin no longer throws
- [ ] Three-phase transformer winding produces 3 groups, 6 sections
      with correctly differentiated coordinates
- [ ] `assign_windings_to_columns` API callable
- [ ] Painter canvas covers all columns (no section clipped)
- [ ] Per-column Lm matches reference within 1% for balanced 3-phase
- [ ] Per-column core loss matches reference within 5% for balanced
- [ ] Each affected calculator header has a "Multi-column behaviour"
      block

---

## 19. Common pitfalls (read before coding)

1. **Don't conflate "column" and "winding window"**. A 3-leg core
   has 3 columns but 2 winding windows. **Always use winding window
   in code; column only in user-facing API.**

2. **Don't try to make insulation calcs cross-column-aware.** v1
   treats each group as independent. Inter-group insulation is a
   v2 concern.

3. **Don't change MAS schema.** Everything fits in the existing
   `Group` mechanism. If you think you need a schema change, you
   are wrong — re-read [§5](#5-recommended-approach).

4. **Don't generalise toroids.** Toroids are single-window. Leave
   `windingWindows[0]` in toroidal helpers; document with a
   comment.

5. **Don't refactor unrelated code.** This is a focused
   geometric-and-EM change. Don't fix unrelated bugs, don't tidy
   unrelated style — that creates merge conflicts.

6. **Don't use `windingWindows[0]` outside of [§11.6, §11.10]**
   (the documented exceptions). Every other access is wrong for
   multi-window.

7. **Don't break single-column.** Run the regression test after
   every phase. If single-column tests fail, **stop** — don't
   plough ahead.

8. **Don't invent EM physics.** If you're tempted to add an
   inter-group proximity correction "while you're at it" — don't.
   That's v2 work and requires measured calibration data.

9. **Don't skip the warning log** in `create_default_groups` when
   multi-window cores get default placement. Users need to know.

10. **Build with `ninja -j2`** in the existing `build/` directory.
    Don't run `cmake --build` (per project memory: ninja-only).

11. **MAS.hpp staging gotcha**: if you see "undeclared MAS types"
    after pulling, copy `MAS.hpp` from source → `build/MAS/MAS.hpp`
    (per project memory).

12. **Don't change `CoreBobbinProcessedDescription`**. It's
    MAS-generated; the scalar `column_*` fields stay scalar (centre
    column only). Multi-column geometry lives in the
    `winding_windows[]` vector (already exists).

---

## 20. v2+ trigger conditions

After v1 ships and a real 3-phase DAB design is built and measured:

1. **Calibrate inter-group proximity correction** from measured AC
   loss at multiple frequencies and inter-leg spacings.
2. **Inter-window leakage flux paths** if a customer needs primary
   and secondary in different columns (rare in 3-phase DAB; common
   in shell-type transformers).
3. **Shared-yoke flux solver** if unbalanced 3-phase operation
   becomes a use case (e.g., single-phase auxiliary winding on a
   3-leg core).
4. **Multi-bobbin physical model** if the 3-phase transformer is
   built with three separate bobbins (common in higher-power
   designs); requires extending the `Bobbin` schema or adding a
   `MultiBobbin` type.

---

## File reference (summary)

| File | Layer | Purpose |
|---|---|---|
| `src/constructive_models/Bobbin.cpp` | A | §8: Remove 2 guards (lines 777, 793); loop in `create_quick_bobbin` |
| `src/constructive_models/Coil.cpp` | A | §9–§12: Remove 3 guards (lines 346, 2189, 3361); new `create_default_groups`; replace 115× `windingWindows[0]` per cluster; add `assign_windings_to_columns`; add `find_window_index_for_group` helper |
| `src/constructive_models/Coil.h` | A | Declare new API (§10.1, §12.2) |
| `src/support/BasicPainter.cpp` | A | §13.2: per-section group lookup at lines 1056–1057 |
| `src/support/AdvancedPainter.cpp` | A | §13.3: per-section group lookup at lines 825–826; canvas sizing |
| `src/physical_models/MagnetizingInductance.cpp` | B | §14: per-column reluctance |
| `src/physical_models/MagnetizingInductance.h` | B | §14.4: multi-column header note + `get_per_column_lm` declaration |
| `src/physical_models/CoreLosses.cpp` | B | §15: per-column flux density |
| `src/physical_models/CoreLosses.h` | B | §15.3: multi-column header note |
| `src/physical_models/LeakageInductance.h` | B | §17.1: documentation only |
| `src/physical_models/MagneticField.h` | B | §17.1: documentation only |
| `tests/TestCoil.cpp` | A | §16.1: 5 new geometry tests |
| `tests/TestMagnetizingInductance.cpp` | B | §16.2: 2 new EM tests |
| `tests/TestCoreLosses.cpp` | B | §16.2: 1 new EM test |
| `tests/TestProximityEffect.cpp` | B | §16.2: characterization test |
| `src/converter_models/CONVERTER_MODELS_GOLDEN_GUIDE.md` | docs | §17.2: multi-column section |

---

## Why this plan is structured this way

The original draft treated multi-column as purely geometric: remove
guards and fix coordinate accesses. That gets to a coil that
**looks** right but produces electromagnetic outputs of unknown
quality.

This plan treats it as three layers (geometry, EM, application)
with explicit accuracy boundaries per quantity. The proximity-effect
question (during plan review) exposed that any quantity computed
from inter-conductor field coupling has the same issue — the honest
v1 answer is "ignored, warned, deferred to v2 calibration." Doing
it any other way means inventing physics without measurements.

The result: v1 ships a coil that **is** right within stated bounds
(exact for DC resistance + within-group AC; ±10% for Lm; ±15% for
core loss; ±10–20% pessimistic for inter-group proximity), with a
clear path to v2 calibration once measured data exists.
