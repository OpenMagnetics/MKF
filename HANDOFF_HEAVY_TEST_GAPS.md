
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
