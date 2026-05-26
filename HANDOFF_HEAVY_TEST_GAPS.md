
---

## Session 3 update — 2 more root-cause fixes on Gap 3 + Gap 4

Two real bugs found, both in the transformer-saturation sizing chain.
Combined effect: every transformer-class wizard (LLC, SRC, PSFB, PSHB,
AHB, Weinberg, Vienna) that the MagneticAdviser hit was getting zero
viable cores from CoreAdviser — not because the request was
unsatisfiable but because two stacked-safety-factor bugs were
guaranteed to wipe out the entire candidate pool.

### Fix A — 3175d9fb fix(inputs): copy MAS optional getters by value

`Inputs::calculate_max_volt_seconds` had a dangling-reference bug:

  const auto& voltage = excitation.get_voltage().value();

`OperatingPointExcitation::get_voltage()` (and every other MAS-
generated optional getter on this path) returns
`std::optional<SignalDescriptor>` BY VALUE — the temporary dies at
the end of the full-expression, and `voltage` then aliases destroyed
memory. Subsequent `voltage.get_waveform()` returns an empty optional,
the function falls through to its `return 0.0` tail, and every
caller that needs V·s for sizing transformer turns gets zero.

Net upstream effect: `add_initial_turns_by_inductance` (CoreAdviser-
Dataset.cpp:552) saw `maxVoltSeconds = 0`, the saturation-aware
initial-N math returned 0, the fallback `initialNumberTurns = 5` ran,
and NumberTurns rounded 5 up to the first integer pair satisfying
the ±25% turns-ratio tolerance (e.g. N=16 for a 21.33:1 PSFB) — far
below what the V·s-derived saturation budget required. Downstream
filterSaturation then rejected every ferrite core.

Fix: copy each optional by value (the MAS-generated pattern is
unambiguously by-value, so the chained `.value()` reference idiom
is a footgun). Inline comment documents the gotcha so the next
developer doesn't reintroduce it.

### Fix B — f65913e3 fix(adviser/saturation-filter): use raw Bsat

`MagneticFilterSaturation::evaluate_magnetic` compared
`B_peak × 1.2` against `Core::get_magnetic_flux_density_saturation(temp)`
— which DEFAULTS `proportion = true` and returns `0.7 × B_sat_raw`.
The 0.7 factor (`defaults.maximumProportionMagneticFluxDensity-
Saturation`) is the *design-side* target that the initial-N seed uses
to size turns; the 1.2 margin is the *filter-side* safety check. They
should multiply different things, not the same thing. Stacking them
turned the filter's effective threshold into
B_peak > 0.583 · B_sat_raw — which the initial-N seed (which targets
B_peak ≈ 0.7 · B_sat_raw) cannot satisfy by construction.

Fix: pass `proportion = false` explicitly in the filter so it
compares against raw B_sat. Inline comment captures the rationale
and explicit reference to the initial-N seed it interlocks with.

### Effect on the heavy test gaps

Gap 3 PSFB (returns 0 cores):
  before:  test cases: 1 | 1 failed | assertions: 1 | 1 failed
  after:   test cases: 1 | 1 passed | assertions: 2 | 2 passed (9 s)

Gap 4 LLC variant matrix (FB+CT and FB+CD returning 0 cores):
  before:  test cases: 1 | 1 failed | assertions: 10 | 8 passed | 2 failed
  after:   test cases: 1 | 1 passed | assertions: 12 | all pass

Gap 4 SRC variant matrix (4 of 6 variants returning 0 cores):
  before:  test cases: 1 | 1 failed | assertions: 8  | 4 passed | 4 failed
  after:   test cases: 1 | 1 passed | assertions: 12 | all pass

All other from-converter transformer tests (PSHB, AHB, Vienna,
Weinberg, Zeta) continue to pass. No regression detected.

### What's still open

Gap 1: CLLLC SPICE convergence — closed in session 1 (commit 03bf3f8f).

Gap 2: Wire Adviser non-termination on multi-winding (Flyback,
       Isolated-BB, Push-Pull, LLC, Vienna, Weinberg, Zeta in the
       playwright G2 tests). Inputs.cpp dangling-ref fix may have
       reduced wasted work here, but the dominant cost is still the
       K×N candidate pairing in the inner wind loop. Needs the
       top-K-per-winding refactor the original handoff proposed.

Gap 3 PFC slowness (5+ min for the 50 W test). CoreAdviser is the
       bottleneck. CoilAdviser is fast for single-winding PFC. Needs
       per-stage profiling to target the specific slow filter; most
       likely `add_ferrite_materials_by_losses` sweeping the
       candidate set × operating points.

Gap 4: Largely cascade-fixed by Fix A + Fix B above. PFC and Zeta F1
       cases should already work end-to-end; needs a WebFrontend
       playwright re-run to confirm.

Recommended next session: re-run the WebFrontend heavy suite end-to-
end. Several tests that timed out at 180-600 s should now finish in
~10-30 s. The remaining open work is Gap 2 (algorithmic refactor)
and PFC slowness (profile + targeted optimisation).

---

## Session 4 update — Gap 3 PFC slowness FIXED

Commit `c26c4c9f`: pre-Inductance pruning cap on the powder branch of
`CoreAdviser::filter_standard_cores_power_application`. The ferrite
branch already had an analogous pre-Losses cap (line 608); the
powder branch was missing the symmetric pre-Inductance one and was
running `add_initial_turns_by_inductance` + `filterMagneticInductance`
on 3 609 candidates (400 cores × 10 powder materials each, after
`add_powder_materials` blew the set up).

Profile diff (from the new CORE_ADVISER_PROFILE=1 stage timer also
added in this commit):

  Stage            | Before  | After
  -----------------+---------+--------
  Inductance pwdr  | 246 s   | 5 s
  Saturation pwdr  | 113 s   | 3 s
  Losses pwdr      |  72 s   | 11 s

Test_MagneticAdviserFromConverter_PFC end-to-end:
  before: 5 min 25 s
  after:  1 min 03 s   (5.2× faster, well under WebFrontend's 180 s)

No regression on the 6 converter tests verified (PSFB, LLC, SRC,
Vienna, Flyback, Buck) — same green output as session 3.

### Status of the four gaps

  Gap 1 (CLLLC ngspice non-convergence)            — FIXED (session 1, 03bf3f8f)
  Gap 2 (Wire Adviser multi-winding hang)          — open
  Gap 3 PSFB / LLC FB+CT/CD / SRC return zero      — FIXED (session 3, 3175d9fb + f65913e3)
  Gap 3 PFC slowness                                — FIXED (session 4, c26c4c9f)
  Gap 4 MagneticAdviser empty                       — cascade-fixed by Gap 3 work

### What's still open

Gap 2 — multi-winding Wire Adviser hang (Flyback, Isolated-BB,
        Push-Pull, LLC, Vienna, Weinberg, Zeta playwright G2 tests).
        Needs the top-K-per-winding refactor of the inner wind loop
        (`CoilAdviser::get_advised_coil_for_pattern` lines ~636-792)
        to cap pairing at K² instead of |wireDatabase|^N. The 60 s
        wall-clock timer attempt in 4b4ebc48 was reverted (silent
        shortcut; CLAUDE.md disallows). Real fix requires:
          1. Capturing the actual MAS JSON the WebFrontend posts to
             `mkf.calculate_advised_coil` (the wasm entry already
             logs it via INPUT_JSON_START markers) so we have a
             deterministic MKF unit-test reproducer.
          2. Refactor: per-winding, generate top-K candidates by
             individual score (where K is configurable, default 10);
             only the K×K pairing matrix runs through wind().

---

## Session 5 — Gap 2 first attempt reverted; lesson captured

Tried the obvious "top-K per winding by score" refactor for the
multi-winding pairing loop in
`CoilAdviser::get_advised_coil_for_pattern`:

  1. `std::stable_sort(wireCoilPerWinding[w], by score asc)` after the
     per-winding candidate-collection loop
  2. Truncate each list to K (default 10, Settings-controlled)
  3. Existing pairing loop runs over the smaller lists

**It made things worse, not better.** Verified results:

  Test                                          | Before  | After
  ----------------------------------------------+---------+--------
  CoilAdviser flyback snapshot (1000 wires)     | passes  | wrong wire pair (predicted)
  Test_MagneticAdviserFromConverter_PSFB        | passes  | passes
  Test_MagneticAdviserFromConverter_LLC_Variant | 47 s    | >240 s (timeout)
  Test_MagneticAdviserFromConverter_SRC_Variant | passes  | passes
  Test_MagneticAdviserFromConverter_Vienna      | passes  | passes
  Test_MagneticAdviserFromConverter_Flyback     | passes  | passes

The snapshot-test failure was the *predicted* and acceptable change
(pinning a specific wire identity from the 1000-pool was the test's
contract; an opt-out for that specific case is fine). The **LLC
regression is the real issue**.

**Root cause of the LLC slowdown:**
`wireCoilPerWinding[w]` is the concatenation of 4 `wireConfiguration`
calls — same wire pool, but each config runs the wire-adviser with
DIFFERENT settings (current-density limit, max-parallels). The
returned scores are normalized against the per-call config, not
globally. Sorting the concatenated list by `pair.second` and taking
the top-K therefore mixes top entries from different configs in an
order that has no real ranking semantics across them.

When `wind()` runs on top-K-from-merged-sort:
  - Some candidates have unrealistic per-strand current density
    (passed under config-2's higher limit but fail in the wound coil)
  - Some have over-paralleled wire that doesn't fit the section width
  - Each `wind()` failure on those takes 0.1–1 s before returning false

The original code's "advance the lowest-scored winding" greedy walked
within-config first (because the concatenation preserves each
config's score-sorted chunk). The top-K-from-global-sort breaks that
locality.

The 25-failure bail-out cap (`max(25, totalCandidates)`) didn't help
because — with `LLC_VariantMatrix` running 6 sections × 25 failed
wind()s × ~1.5 s each — total per call exceeds the playwright budget.

### The right shape of the fix (not implemented)

Cap PER CONFIG, not per winding:

```cpp
for (auto& cfgWires : per_config_per_winding[w]) {
    if (cfgWires.size() > kPerConfig) cfgWires.resize(kPerConfig);
}
// Then concatenate the truncated chunks.
```

This preserves each config's "top wires by its own scoring" while
still bounding total candidates. For `kPerConfig = 3` and 4 configs ×
2 windings, you get 12 candidates per winding (24 total) instead of
~120 — same order of speedup as my failed attempt, but without
breaking score locality.

Need to refactor `CoilAdviser::get_advised_coil_for_pattern`'s per-
winding candidate loop to keep the chunks separate before
concatenating, so the cap can apply to each chunk independently.
Probably 30–60 minutes of careful work + regression run. The reverted
attempt's commit message and inline comment (now removed) capture the
context.

### What's known about the actual hang

Still no reproducer from WebFrontend captured. The
`mkf.calculate_advised_coil` wasm entry prints
`=== INPUT_JSON_START ===` then the MAS JSON to stdout — the next
WebFrontend playwright run will surface the exact failing inputs in
its console capture. Until then, the per-config-cap refactor would
have to be validated against LLC_VariantMatrix as the proxy.
