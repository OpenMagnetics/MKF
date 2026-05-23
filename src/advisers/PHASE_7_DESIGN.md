# Phase 7 design — cross-referencer unification

**Status:** design proposal — no code changes
**Author:** drafted 2026-05-22 while validating Phase 6
**Goal:** reduce duplication between
`CoreCrossReferencer::MagneticCoreFilter*` (5 classes),
`CoreMaterialCrossReferencer::MagneticCoreFilter*` (7 classes),
and the main `MagneticFilter` hierarchy (~25 classes).

## API shapes — where the three frameworks actually differ

| Aspect | `MagneticFilter` | `CoreCrossReferencer` | `CoreMaterialCrossReferencer` |
|---|---|---|---|
| Subject of evaluation | `Magnetic*` (whole magnetic) | one `Core` | one `CoreMaterial` |
| Required context | `Inputs*` + `std::vector<Outputs>*` | `Inputs` + `referenceCore` | `temperature` + `referenceCoreMaterial` |
| Op signature | `evaluate_magnetic(...) -> std::pair<bool,double>` | `filter_core(list*, ref, ..., weight, limit) -> filtered_list` | `filter_core_materials(list*, ref, temp, weight) -> filtered_list` |
| Iteration model | Caller iterates candidates, filter is per-call stateless | Filter iterates batch internally | Filter iterates batch internally |
| Scoring storage | Returned directly to caller | Embedded raw-pointer maps (`*_scorings`, `*_validScorings`, `*_scoredValues`, `*_filterConfiguration`) | Same shape, different enum |
| Statefulness | Mostly stateless | Holds external map pointers + per-instance enum | Same |
| Returned shape | `(pass: bool, score: double)` per candidate | `vector<pair<Core, double>>` survivors | `vector<pair<CoreMaterial, double>>` survivors |

Two real frameworks live in this code. They look superficially alike
because both produce `(name → score)` tables, but the contracts diverge
in five non-trivial places (subject type, context bundle, iteration
ownership, scoring storage, return shape). Forcing them through a
single API would either widen `MagneticFilter` to a sum-type input (and
break ~25 subclasses), or wrap cross-referencer filters in adapters
that lie about what they evaluate. Neither pays for itself.

## Options considered

### Option A — Generalize `MagneticFilter::evaluate_magnetic` with a variant input

Make the base method take e.g. `std::variant<Magnetic*, Core*, CoreMaterial*>` so all three filter families can subclass it.

- **Pro:** Single `factory(...)` for everything; single mental model.
- **Con:** Breaks every existing `MagneticFilter` subclass signature. ~25 classes touched. The Outputs* third argument is meaningless for cross-referencers. Reference-comparison filters need a second subject (`referenceCore`/`referenceCoreMaterial`) that `MagneticFilter` doesn't carry. **Forces an awkward API on all callers to share a base class that doesn't model the actual contract well.**
- **Verdict:** Reject. The complexity tax outweighs the conceptual benefit.

### Option B — Thin `CoreFilterAdapter` wrappers

Wrap each cross-referencer filter so it implements `evaluate_magnetic`. The adapter:
1. Extracts `Magnetic->get_core()` (or its material) as the subject.
2. Calls the wrapped cross-referencer with a single-element list.
3. Returns the resulting pair as `(true/false, score)`.

- **Pro:** No changes to `MagneticFilter` base. Cross-referencer code stays as-is internally.
- **Con:** Adapter must hold the `referenceCore` / `referenceCoreMaterial` as state, plus the four shared scoring maps. That's exactly the state the wrapped class already has — adapters end up duplicating the wrapped object's identity. Single-element `vector<pair<...>>` round-trip is wasteful per call. **You add an indirection without removing any code.**
- **Verdict:** Reject. Indirection without payoff.

### Option C — Extract a templated `CrossReferencerBase<T>` to collapse the two cross-referencer files

`CoreCrossReferencer` and `CoreMaterialCrossReferencer` are near-clones differing only in subject type (`Core` vs `CoreMaterial`) and the enum they key their score maps on. Templating on:

```cpp
template <typename Subject, typename FilterEnum>
class CrossReferencerBase { /* shared scoring maps, normalize, apply_filters */ };
```

…folds the two files. The Phase 3 F2+F3 audit estimated **~600–750 LOC saved**.

- **Pro:** Concrete payoff. Each cross-referencer-specific filter remains the only thing in its own file. Bug fixes in the shared base apply to both.
- **Con:** C++ template error messages are noisier when something goes wrong. `MagneticCoreFilterCoreLosses` (CCR) and `MagneticCoreFilterVolumetricLosses` (CMCR) have model-management code that doesn't cleanly template (different storage shapes); they stay as specializations. Some adviser-specific filter logic (saturation, permeance, etc.) is intrinsic to the subject type — not generalizable.
- **Verdict:** Accept as the primary Phase 7 win. ~600 LOC saved is real; risk is bounded to the two files being templated.

### Option D — Consolidate `normalize_scoring` only

The same min-max / log-normalize / neutral=0.5 fallback pattern appears in `CoreCrossReferencer`, `CoreMaterialCrossReferencer`, `WireAdviser`, and parts of `CoreAdviser`. Phase 3 F12 already pulled the helper into a shared utility; some sites still re-implement it inline.

- **Pro:** Tiny risk. ~50–150 LOC saved. Bisectable per-site.
- **Con:** Doesn't address the structural duplication between the two cross-referencer hierarchies.
- **Verdict:** Accept as a complementary cleanup that can land before C with no design dependency.

## Recommendation

**Do C + D, skip A and B.**

Concretely:

1. **First, finish D.** Sweep the four cross-referencer / adviser files for any `normalize_scoring`-like inline duplicates that escaped Phase 3 F12. Land as one commit per site, no behaviour change, char tests stay green.
2. **Then C.** Introduce `CrossReferencerBase<Subject, FilterEnum>` in a new header `src/advisers/CrossReferencerBase.h`. Migrate the *scoring storage*, *normalize step*, and *apply_filters dispatch loop* into the template. Keep the per-subject filter classes (`MagneticCoreFilterPermeance`, `MagneticCoreFilterInitialPermeability`, …) as the concrete specializations. The two `CoreCrossReferencer::filter_xxx` / `CoreMaterialCrossReferencer::filter_xxx` driver functions become thin wrappers over `CrossReferencerBase::filter`. Land in chunks: (i) template scaffold + one filter migrated, (ii) remaining 4 CCR filters, (iii) remaining 6 CMCR filters, (iv) delete the dead member functions on the old base.
3. **Do not unify with `MagneticFilter`.** The API mismatch is real, not cosmetic. Document this decision in the new header so a future contributor doesn't reopen it.

## Expected outcome

- ~600 LOC removed from `src/advisers/`.
- One template base (`CrossReferencerBase<T, EnumT>`) replaces two near-identical hand-written bases.
- `MagneticFilter` framework left alone — its 25 subclasses don't pay any tax for this refactor.
- Char tests stay 314 / 55 green throughout (each step is a refactor, not a behaviour change).
- Benchmark unchanged (or marginally better — template instantiation may inline a couple of small helpers).

## Risks and rollback

- **Template error noise** during the migration — bounded to the two cross-referencer files; mitigated by migrating one filter at a time.
- **MagneticCoreFilterCoreLosses / MagneticCoreFilterVolumetricLosses** carry model-management state (`_coreLossesModelNames`, `_coreLossesModels` vector of model+factory pairs) that doesn't template cleanly. Both stay as specializations holding their own state — the template only owns the scoring-map shape and the normalize/dispatch loop.
- **Rollback** at any step is a single git revert; the per-filter migrations are independent.

## Out of scope

- Anything that touches `MagneticFilter.h` or its 25 subclasses.
- Replacing the raw scoring-map pointers (`std::map<...>*`) with something safer — separate refactor, would be its own design doc.
- The `LOGIC-4 NOTE` on the fixed frequency grids in `MagneticCoreFilterCoreLosses` / `MagneticCoreFilterVolumetricLosses` — those are TODO markers for a different fix (sample actual operating frequency), not a unification concern.
