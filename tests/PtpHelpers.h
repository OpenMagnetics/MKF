// PtP regression helpers — shared by the per-topology
// Test*ReferenceDesignsPtp.cpp files.
//
// Four pure functions previously inline-duplicated in 11 PtP files:
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
// from multiple translation units.

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
