# Handoff — Fast-Path Ferrite Gapping (Non-Toroidal)

**Date:** 2026-05-27
**Branch:** `main` (5 commits ahead of `origin/main`, plus 1 uncommitted edit)
**Status:** WIP — implementation drafted, fast-path tests pass, broader regression not yet swept.

---

## Problem

`MagneticAdviser::get_advised_magnetic_fast` (src/advisers/MagneticAdviser.cpp)
is missing the gapping step for non-toroidal ferrite cores.

Current pipeline:
1. AreaProduct filter
2. EnergyStored filter
3. **Step 2c** — reject toroidal ferrite for energy-storing topologies
   (added commit `7b1517b5`).
4. `add_initial_turns_by_inductance`
5. `filterSaturation`
6. Wind + losses

Step 4 runs against whatever gap each catalog core arrived with (typically
0 for ferrite blanks). With ungapped high-µ ferrite, `add_initial_turns_by_inductance`
solves L = N²·A_L for tiny N (1–3 turns). Resulting B_peak = L·I_peak/(N·A_e)
saturates instantly. `filterSaturation` (step 5) drops most, but:

- The few survivors are accidentally-OK candidates served ungapped (wrong
  design for any application that stores meaningful flux).
- The slow path (`CoreAdviserPipeline::process_standard_cores`, lines
  573–602) does the right thing: gap non-toroidal ferrite via
  `add_gapping_standard_cores` **before** turns calculation.

User question that triggered this: "wait, we are returning flybacks without
gaps?? why??" — the toroidal-rejection fix in `7b1517b5` plugged the
worst-case (ungappable ferrite toroids) but left the gap-on-discrete-cores
issue.

## What's been done

5 committed fixes on `main`:

| Commit | Subject |
|---|---|
| `c9d7674f` | fix(isolated-buck): seed Cout IC=0 when V_sec,peak < V_out spec |
| `7b1517b5` | fix(magnetic-adviser/fast): reject ferrite toroids for energy-storing topologies |
| `95685f3d` | fix(acf): canonical clamp topology (Cclamp across primary, not to GND) |
| `1d4f6306` | docs(handoff): session 6 — Gap 2 fixed via per-config wire cap |
| `c1c21964` | perf(coil-adviser): cap per-config wire candidates at K=10 (Gap 2) |

**Uncommitted draft** (in working tree on `src/advisers/MagneticAdviser.cpp`):
inserts a new "Step 2d" between the toroidal-ferrite rejection and
`add_initial_turns_by_inductance`. The new step splits candidates by
toroidal vs non-toroidal, calls `coreAdviser.add_gapping_standard_cores`
on the non-toroidal slice, then concatenates them back.

```cpp
// Step 2d: Gap non-toroidal ferrite cores. ...
{
    std::vector<std::pair<Magnetic, double>> nonToroidal;
    std::vector<std::pair<Magnetic, double>> toroidal;
    nonToroidal.reserve(magneticsWithScoring.size());
    for (auto& entry : magneticsWithScoring) {
        if (entry.first.get_core().get_functional_description().get_type()
            == CoreType::TOROIDAL) {
            toroidal.push_back(std::move(entry));
        } else {
            nonToroidal.push_back(std::move(entry));
        }
    }
    if (!nonToroidal.empty()) {
        coreAdviser.add_gapping_standard_cores(&nonToroidal, inputs);
    }
    magneticsWithScoring = std::move(nonToroidal);
    for (auto& entry : toroidal) {
        magneticsWithScoring.push_back(std::move(entry));
    }
}
```

### What's been verified on the uncommitted draft
- Builds clean: `ninja -C build -j5 MKF_tests` → 0 errors.
- Fast-path-tagged tests: `./MKF_tests "[fast]"` → 20 assertions / 3 cases PASS.

### What has NOT been verified
- The broader `[adviser]` suite. **A prior attempt** at this same change
  (in an earlier session) regressed exactly **1 test on top of the
  19 pre-existing failures** on the clean tree. That regression was
  the reason commit `7b1517b5` only added the toroidal rejection and
  explicitly avoided calling `add_gapping_standard_cores`. The exact
  failing test and the root cause were never identified.

## How to reproduce / continue

```bash
cd /home/alf/OpenMagnetics/MKF

# 1. See the uncommitted change
git diff src/advisers/MagneticAdviser.cpp

# 2. Build
ninja -C build -j5 MKF_tests

# 3. Sanity: fast path
./build/MKF_tests "[fast]"          # expect: 20 assertions PASS

# 4. Baseline: stash the change, build clean, run broader suite, record failures
git stash push src/advisers/MagneticAdviser.cpp
ninja -C build -j5 MKF_tests
./build/MKF_tests "[adviser]" 2>&1 | tee /tmp/baseline_adviser.log
grep -E "FAILED|test cases" /tmp/baseline_adviser.log | tail -30

# 5. Apply the change, build, re-run, diff
git stash pop
ninja -C build -j5 MKF_tests
./build/MKF_tests "[adviser]" 2>&1 | tee /tmp/withfix_adviser.log
diff <(grep "FAILED\|test cases" /tmp/baseline_adviser.log) \
     <(grep "FAILED\|test cases" /tmp/withfix_adviser.log)
```

The diff will pinpoint the regressing test. Then:
- Read the test and figure out **why** gapping changes its outcome.
  Candidates: (a) the test uses a non-ferrite non-toroidal core that
  `add_gapping_standard_cores` doesn't handle correctly, (b) the test's
  expected core/gap matches the pre-fix ungapped behavior and needs
  updating, (c) `add_gapping_standard_cores` has a bug that only the
  fast-path-flow exposes (in the slow path the dimensions/fringing
  filters run after gapping and may mask it).
- Heaviside corpus reference (where this issue surfaced):
  `/home/alf/OpenConverters/Heaviside/docs/MKF-HANDOFF.md` §1
  (flyback served `T 38.1/19.05/12.7 N30` ungapped, isat 0.63 A vs
  ipeak 1.8 A). With the toroidal rejection from `7b1517b5` that
  specific failure no longer occurs, but non-toroidal ferrite in
  flyback specs is still served ungapped.

## Open follow-ups (unrelated to this fix)

- 5 commits ahead of `origin/main` — not pushed.
  Push with: `git -c core.sshCommand="/usr/bin/ssh -i ~/.ssh/hephaestus_om" push origin main`
- Other agent's working changes in `src/constructive_models/Coil.cpp`/`.h`
  and untracked `src/converter_models/SPICE_PROBE_HANDOFF.md` left untouched.
- Heaviside corpus re-run (PyMKF rebuild) — Heaviside's job, not MKF's.
