// PtP regression helpers — shared by the per-topology
// Test*ReferenceDesignsPtp.cpp files.
//
// ─────────────────────────────────────────────────────────────────────────
// THE CANONICAL 4-GATE PtP HARNESS PATTERN
// ─────────────────────────────────────────────────────────────────────────
//
// As of phase 8l (May 2026) every implemented converter topology in
// src/converter_models/ has a paired Test<Topology>ReferenceDesignsPtp.cpp
// file that runs three vendor reference designs through a four-gate
// regression: walltime, load consistency, power balance, and primary-
// current shape NRMSE (analytical model vs SPICE switching simulation).
//
// File template (see TestBuckReferenceDesignsPtp.cpp for the canonical
// example, or TestFlybackReferenceDesignsPtp.cpp for the flyback-class
// variant):
//
//   1. RefDesignSpec struct literal per design — Vin, Vout, Iout, Fs,
//      transformer specs, plus the four gate tolerances (tol_walltime,
//      tol_rload_pct, tol_loss_max, tol_nrmse).
//
//   2. build(spec) helper — sets diodeVoltageDrop=0, efficiency=1.0
//      (lossless analytical anchor), wires the operating point with
//      Vout/Iout sized for the desired load (R = Vout_nom/Iout_nom).
//
//   3. consistent_lm(spec) helper — picks Lm so analytical primary-
//      current ripple matches SPICE (the topology-specific formula
//      derived in the corresponding TestTopology<X>.cpp).
//
//   4. run_ptp_gates(spec) — runs analytical pass, then SPICE switching
//      pass with set_num_steady_state_periods(400+) for cap settling,
//      then the four CHECKs:
//
//        Gate 1  walltime              < tol_walltime
//        Gate 2  rload consistency     |Vout_spice/Iout_spice − Rload_nom|
//                                      / Rload_nom < tol_rload_pct/100
//        Gate 3  power balance         flavour depends on topology — see
//                                      "Gate 3 variants" below
//        Gate 4  primary-current NRMSE ptp_nrmse(analytical, SPICE)
//                                      < tol_nrmse
//
//   5. Three TEST_CASE entries with tag set
//      [converter-model][<topo>-topology][refdesign][ptp][slow], one per
//      RefDesignSpec literal.
//
// ─────────────────────────────────────────────────────────────────────────
// GATE 3 VARIANTS — choose by topology
// ─────────────────────────────────────────────────────────────────────────
//
//   A.  Pin/Pout balance, single-output, asymmetric:
//         loss = (pin − pout) / pin;  CHECK(loss >= 0); CHECK(loss <= max);
//       Use for Buck, Boost, ACF, AHB, CLLC, DAB, LLC, PSHB, PSFB, SSF,
//       TSF, Flyback, PushPull. Open-loop SPICE Vout settles within ~5 %
//       of nominal so Pin and Pout are physically consistent and Pin
//       must dominate (no energy creation). Asymmetric bound catches
//       both sign errors and excess loss.
//
//   B.  Pin/Pout balance, multi-output, symmetric:
//         loss = (pin − pout_total) / pin;  CHECK(|loss| <= tol);
//       pout_total sums all secondaries. Use for IsolatedBuckBoost.
//       Open-loop Vout drifts up to 2× from nominal so the load-vs-
//       nominal comparison is meaningless, but eta=1, Vd=0 still implies
//       energy conservation. Symmetric bound permits ±tol because
//       cap-charge dynamics over a single-period extraction window
//       produce small bidirectional residuals.
//
//   C.  Pout(total) vs Pout_nominal — DEPRECATED measurement design:
//       attempted for IsolatedBuck (commit 1719af04) but is fragile
//       because it gates on whether SPICE landed on the analytical
//       operating point, not on physical correctness. Prefer A or B.
//
// ─────────────────────────────────────────────────────────────────────────
// PER-TOPOLOGY GATE-TIGHTENING PROTOCOL
// ─────────────────────────────────────────────────────────────────────────
//
// New PtP files start with loose gates (walltime=30, rload=30 %,
// loss=0.60, NRMSE=0.30..0.55). Procedure to tighten:
//
//   1. Build + run with loose gates: ./MKF_tests "[<topo>][refdesign][ptp]"
//      --durations yes
//
//   2. Note the maximum observed value across the three cases for each
//      gate. Capture from the std::cout banners emitted by
//      run_ptp_gates().
//
//   3. Set tol = ~1.3× max observed for each gate. The 1.3× factor
//      accommodates platform-dependent SPICE timestep variability and
//      mild ngspice version drift; tighter risks flakes, looser misses
//      regressions.
//
//   4. Walltime budget: <10 s/case mandate. If a topology needs more
//      settling periods (flyback-class), reduce inner-loop count first
//      and only relax walltime as a last resort.
//
//   5. Rebuild + run + verify all three cases pass, then commit:
//      "Add canonical PtP harness for 3 <Topology> reference designs"
//      with body explaining gate rationale and observed numbers.
//
//   6. Run full PtP regression to verify nothing else broke:
//      ./MKF_tests "[refdesign][ptp]"
//      Current passing baseline (May 2026): 12 topologies, 48 cases,
//      939 assertions.
//
// ─────────────────────────────────────────────────────────────────────────
// ADDING A NEW TOPOLOGY PtP FILE
// ─────────────────────────────────────────────────────────────────────────
//
//   1. Read src/converter_models/<X>.h for setter API + simulate_and_
//      extract_topology_waveforms / simulate_and_extract_operating_points
//      signatures.
//
//   2. Read tests/TestTopology<X>.cpp for the existing
//      kRefDesign{1,2,3} / RefDesignSpec literals and the
//      consistent_lm formula. Reuse spec values verbatim.
//
//   3. Copy a similar topology's PtP file as the template (Buck/Boost
//      for non-isolated, Flyback/PushPull for isolated, IsolatedBuck-
//      Boost for multi-secondary flyback-class).
//
//   4. cmake reconfigure (CMake GLOB picks up new test files):
//        cd build && cmake . && ninja -j8 MKF_tests
//
//   5. Apply the gate-tightening protocol above.
//
// ─────────────────────────────────────────────────────────────────────────
// HELPER FUNCTIONS PROVIDED BY THIS HEADER
// ─────────────────────────────────────────────────────────────────────────
//
//   ptp_interp        Linear resample of a (t,d) pair onto a uniform
//                     N-point grid spanning [0, t.back()].
//   ptp_nrmse         Mean-subtracted, scale-normalized NRMSE between
//                     two equal-length series with circular phase
//                     alignment over the first 64 shifts.  Output is
//                     shape similarity in [0,1] — invariant to DC
//                     offset and amplitude scale.  Returns 1.0 if
//                     either series has near-zero AC content.
//   ptp_current       Extract the per-winding current waveform data
//                     vector from an OperatingPoint.  Empty vector if
//                     the winding/current/waveform is missing.
//   ptp_current_time  Same, for the time vector.
//
// All four are inline so this header is header-only and safe to include
// from multiple translation units. Place
// `using namespace OpenMagnetics::Testing;` inside the anon namespace
// of the test file (or after `using namespace OpenMagnetics;` at file
// scope) — `OpenMagnetics::` does not auto-import nested namespaces.
//
// ─────────────────────────────────────────────────────────────────────────
// KNOWN LIMITATIONS / DEFERRED WORK
// ─────────────────────────────────────────────────────────────────────────
//
//   * Closed-loop Vout regulator (would let Gate 3 use Pout-vs-nominal
//     for flyback-class converters): attempted in May 2026 but blocked
//     by the IBB SPICE model — Vout settling point is pinned ~22 V
//     regardless of analytical duty cycle, even at D=0.96. Root cause
//     in src/converter_models/IsolatedBuckBoost.cpp generate_ngspice_
//     circuit() is not yet identified; tracked separately.
//
//   * TestTopologyBuck.cpp:30 has its own anon-namespace ptp_interp for
//     Buck-specific shape checks — leave alone.
//
//   * TestPfcReferenceDesignsPtp.cpp uses different helpers
//     (interp / mean_over) suited to PFC envelope analysis — do not
//     migrate to ptp_interp / ptp_nrmse.
//
// ─────────────────────────────────────────────────────────────────────────

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "MAS.hpp"

namespace OpenMagnetics {
namespace Testing {

inline std::vector<double> ptp_interp(const std::vector<double>& t,
                                      const std::vector<double>& d,
                                      int N) {
    std::vector<double> out(N);
    double T = t.back();
    for (int i = 0; i < N; ++i) {
        double ti = T * i / (N - 1);
        int lo = 0;
        for (int k = 0; k + 1 < (int)t.size(); ++k) if (t[k] <= ti) lo = k;
        int hi = std::min(lo + 1, (int)t.size() - 1);
        double dt = t[hi] - t[lo];
        out[i] = (dt < 1e-20) ? d[hi]
                              : d[lo] + (ti - t[lo]) / dt * (d[hi] - d[lo]);
    }
    return out;
}

inline double ptp_nrmse(const std::vector<double>& ref,
                        const std::vector<double>& sim) {
    int N = (int)ref.size();
    double ref_mean = 0.0, sim_mean = 0.0;
    for (int i = 0; i < N; ++i) { ref_mean += ref[i]; sim_mean += sim[i]; }
    ref_mean /= N; sim_mean /= N;
    std::vector<double> r(N), s(N);
    double rAC = 0.0, sAC = 0.0;
    for (int i = 0; i < N; ++i) {
        r[i] = ref[i] - ref_mean; s[i] = sim[i] - sim_mean;
        rAC += r[i] * r[i]; sAC += s[i] * s[i];
    }
    rAC = std::sqrt(rAC / N); sAC = std::sqrt(sAC / N);
    if (rAC < 1e-10 || sAC < 1e-10) return 1.0;
    for (int i = 0; i < N; ++i) { r[i] /= rAC; s[i] /= sAC; }
    int ns = std::min(N, 64);
    double best = std::numeric_limits<double>::max();
    for (int ss = 0; ss < ns; ++ss) {
        int sh = ss * N / ns;
        double ssd = 0.0;
        for (int k = 0; k < N; ++k) {
            double e = r[k] - s[(k + sh) % N];
            ssd += e * e;
        }
        if (ssd < best) best = ssd;
    }
    return std::sqrt(best / N);
}

inline std::vector<double> ptp_current(const MAS::OperatingPoint& op,
                                       std::size_t wi = 0) {
    if (wi >= op.get_excitations_per_winding().size()) return {};
    auto& e = op.get_excitations_per_winding()[wi];
    if (!e.get_current() || !e.get_current()->get_waveform()) return {};
    return e.get_current()->get_waveform()->get_data();
}

inline std::vector<double> ptp_current_time(const MAS::OperatingPoint& op,
                                            std::size_t wi = 0) {
    if (wi >= op.get_excitations_per_winding().size()) return {};
    auto& e = op.get_excitations_per_winding()[wi];
    if (!e.get_current() || !e.get_current()->get_waveform()) return {};
    auto tv = e.get_current()->get_waveform()->get_time();
    return tv ? tv.value() : std::vector<double>{};
}

} // namespace Testing
} // namespace OpenMagnetics
