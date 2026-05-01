# Multi-Column Coil Customizer — WebFrontend Plan

Companion to `MULTI_COLUMN_WINDING_SUPPORT_PLAN.md` (MKF-side
infrastructure). This plan covers the **WebFrontend** side: how
users select per-winding column placement, what UI/UX makes sense,
and how the data flows down to MKF/WASM.

> Scope: WebFrontend only (NOT WebFrontend2 — per project memory,
> WebFrontend2 is ignored).

---

## Goal

Let a user with a multi-window core (e.g., 3-leg E-core, EI core)
choose which **winding window (column)** each winding goes into,
and have the result render correctly in the magnetic preview.

---

## When does the UI surface a column selector?

The core's processed description has `windingWindows.size() > 1`.
Single-window cores keep the existing flow (no UI change).

Two cases trigger:
1. User picks a multi-window core shape (E, EI, ER variants with 2
   windows; rare 5-leg shapes with 4 windows).
2. User explicitly opts into "advanced coil customizer" mode.

---

## Two UI modes

### Mode A — Simple (default)
A toggle in the coil step: **"Distribute windings across columns"**.
Default OFF → all windings go in column 0 (current behaviour).
When ON → reveals Mode B's column picker.

This keeps the simple case simple.

### Mode B — Advanced Coil Customizer
A new section in the coil step. Layout:

```
┌─────────────────────────────────────────────────┐
│  Column placement                               │
│                                                 │
│  Core: ETD49 (2 winding windows)                │
│                                                 │
│  ┌───────────────────────────────────────────┐  │
│  │  [visual: top-down core slice with        │  │
│  │   windows highlighted, drag-drop targets] │  │
│  └───────────────────────────────────────────┘  │
│                                                 │
│  Winding 1 (Primary, 24T)         [Column 0 ▼] │
│  Winding 2 (Secondary, 6T)        [Column 1 ▼] │
│  Winding 3 (Aux, 3T)              [Column 0 ▼] │
│                                                 │
│  ⚠ Multi-column EM accuracy:                    │
│  • Magnetizing inductance: ±10% (balanced)      │
│  • Core losses: ±15%                            │
│  • Inter-column proximity: ignored (v1)         │
│  [Read more]                                    │
└─────────────────────────────────────────────────┘
```

Components:
- **Visual core slice** — schematic showing the core's columns and
  winding windows. Renders in negative + positive X.
- **Winding-to-column dropdowns** — each winding gets a dropdown
  with one entry per winding window (`Column 0`, `Column 1`, …).
- **Optional drag-drop** — drag winding chips into column drop
  targets in the visual. Nice-to-have, not required for v1 UI.
- **Accuracy disclaimer** — surfaces the v1 EM-accuracy bounds
  from the MKF plan's accuracy table. Don't bury this; users need
  to know.

---

## Data flow

### Frontend → WASM bridge

The MAS Coil schema already supports the Group mechanism — no
schema change. The frontend builds:

```ts
// In the coil object that gets serialised to JSON for WASM:
coil.groupsDescription = [
  {
    name: "Column 0",
    coordinates: [/* matches windingWindows[0].coordinates */],
    dimensions: [/* width, height of window 0 */],
    coordinateSystem: "cartesian",
    sectionsOrientation: "overlapping",
    type: "wound",
    partialWindings: [
      { winding: "Primary", parallelsProportion: [1] },
      { winding: "Aux",     parallelsProportion: [1] },
    ],
  },
  {
    name: "Column 1",
    coordinates: [/* windingWindows[1].coordinates */],
    dimensions: [/* width, height of window 1 */],
    coordinateSystem: "cartesian",
    sectionsOrientation: "overlapping",
    type: "wound",
    partialWindings: [
      { winding: "Secondary", parallelsProportion: [1] },
    ],
  },
];
```

The MKF side's `assign_windings_to_columns` is a convenience
wrapper that builds this same structure programmatically — the
frontend can either:
1. Build `groupsDescription` directly and pass it (recommended; no
   extra WASM call).
2. Pass column-index arrays and call a new WASM-exported function
   `assignWindingsToColumns(coilJson, indicesArray)` that returns
   the updated coil JSON.

**Recommended path**: option 1. Frontend has all the info it needs;
no need to add a new WASM export.

### WASM exports needed
None new. The existing `process_operating_point` (or whatever
takes a Coil JSON) accepts a `groupsDescription` field and MKF
will use it.

---

## Frontend implementation tasks

### Task 1 — Detect multi-window core
Where: the coil-step component (probably
`src/components/Coil/CoilStep.vue` or similar).
Check: `core.processedDescription.windingWindows.length > 1`.
Reveal: the "Distribute windings across columns" toggle.

### Task 2 — Column dropdown UI
Where: a new sub-component, e.g., `WindingColumnPicker.vue`.
Props:
- `windings: Winding[]`
- `windingWindows: WindingWindowElement[]`
- `value: number[]  // column index per winding`
Emits: `input` (the array, when user changes selection).

### Task 3 — Visual core slice
Where: a new sub-component, e.g., `CoreColumnVisualizer.vue`.
Renders the core cross-section using SVG. Highlights each winding
window. Shows winding chips overlaid on their assigned window.

For v1 UI: a simple rectangular schematic is enough. For v2:
hook into the AdvancedPainter's SVG output so the visual matches
what the final render will look like.

### Task 4 — Wire up to coil JSON
On change → rebuild `coil.groupsDescription` → save to coil state.
On WASM call → groups go through.

### Task 5 — Accuracy disclaimer
Where: inline below the column picker.
Source: pull from a static list mirroring MKF's accuracy table.
Surface a "[Read more]" expander with the v1 accuracy table.

### Task 6 — Persistence / load
When loading a saved design with `groupsDescription` already set:
- Detect that columns are already distributed.
- Pre-populate the column dropdowns from the saved
  `groupsDescription` partial-windings → group mapping.
- Show the multi-column toggle as ON.

### Task 7 — Validation
Frontend should validate:
- Every winding is assigned to exactly one column (no
  duplicates, no missing).
- The number of selectable columns matches
  `windingWindows.length`.

If invalid, disable the "Continue" button with a clear error.

---

## UX edge cases

| Case | Handling |
|---|---|
| Single-window core | Hide the column picker entirely. |
| User picks single-column placement on a multi-window core | OK — equivalent to default behaviour. |
| Empty column (no windings) | Allow it. MKF supports empty groups (placeholder). |
| User changes core shape mid-design | Reset column assignments to default (all in column 0); show a small toast "Column assignments reset because new core has different windows". |
| Loading old design without `groupsDescription` | Treat as "all in column 0", toggle OFF. |

---

## Visual rendering

The MKF `AdvancedPainter` (after Phase 6 changes) extends the canvas
symmetrically when `windingWindows.size() > 1`. Sections placed in
column 1 (negative-X coordinates) will render correctly without
frontend changes.

For the **preview** in the coil customizer (before pressing
"Continue"), the frontend can show a lightweight SVG schematic
that mirrors what MKF will render. This avoids round-tripping to
WASM for every drag.

---

## Rollout

1. **Phase F1** (1–2 days) — Detect multi-window, simple toggle +
   per-winding dropdowns. No visual core slice yet. Ship as
   "experimental" feature flag.
2. **Phase F2** (1 week) — Add `CoreColumnVisualizer` SVG, drag-drop
   from winding chips to column drop targets. Polish.
3. **Phase F3** (later) — Hook up live preview by re-running the
   AdvancedPainter on each change.

---

## Out of scope (v2)

- Multi-bobbin physical model (one bobbin per leg) — needs MKF
  schema extension.
- Per-section column assignment (currently all sections of a winding
  are in the same column). Useful for split-primary across columns
  but uncommon.
- Core-shape-specific recommendations (e.g., "for 3-phase, suggest
  3 windings, one per leg").
