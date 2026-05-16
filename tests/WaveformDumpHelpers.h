// WaveformDumpHelpers.h — side-by-side analytical-vs-SPICE waveform CSVs.
//
// Use when an analytical model and a SPICE simulation disagree by more
// than the DAB-quality NRMSE gate (≤ 0.15) and you need to see the
// actual shape mismatch.  Writes a single CSV file to /tmp containing
// one column for time plus one analytical column and one SPICE column
// per winding, resampled to a common uniform grid spanning the shorter
// of the two simulation horizons.
//
// Typical use:
//
//   #include "WaveformDumpHelpers.h"
//   using namespace OpenMagnetics::Testing;
//   ...
//   auto a = topo.process_operating_points(tr, Lm);
//   auto s = topo.simulate_and_extract_operating_points(tr, Lm);
//   dump_waveforms_csv("psfb_telecom", a[0], s[0]);
//
// The function logs the CSV path to stderr.  Inspect with:
//   gnuplot -p -e "plot '<path>' u 1:2 w l, '' u 1:3 w l"
//   python -c "import pandas as pd; print(pd.read_csv('<path>').describe())"
//
// Inert in CI: gated on the MKF_DUMP_WAVEFORMS env var.  If unset the
// function is a no-op and returns an empty string, so calls left in a
// test file do not slow down regression runs.

#pragma once

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "MAS.hpp"
#include "PtpHelpers.h"

namespace OpenMagnetics {
namespace Testing {

namespace detail {

// Extract a (time, data) pair from a SignalDescriptor's waveform, or {} {}.
inline std::pair<std::vector<double>, std::vector<double>>
sig_wf(const std::optional<MAS::SignalDescriptor>& s) {
    if (!s || !s->get_waveform()) return {{}, {}};
    auto w = *s->get_waveform();
    auto tv = w.get_time();
    return { tv ? tv.value() : std::vector<double>{}, w.get_data() };
}

inline std::vector<double>
resample_or_zeros(const std::vector<double>& t,
                  const std::vector<double>& d, int N) {
    if (t.empty() || d.empty()) return std::vector<double>(N, 0.0);
    return ptp_interp(t, d, N);
}

inline double safe_back(const std::vector<double>& v) {
    return v.empty() ? std::numeric_limits<double>::max() : v.back();
}

} // namespace detail

inline std::string dump_waveforms_csv(const std::string& tag,
                                      const MAS::OperatingPoint& analytical,
                                      const MAS::OperatingPoint& spice,
                                      int N = 1024) {
    const char* enable = std::getenv("MKF_DUMP_WAVEFORMS");
    if (!enable || std::string(enable) == "0" || std::string(enable).empty()) {
        return {};
    }

    auto ts = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::system_clock::now().time_since_epoch()).count();
    std::string path = "/tmp/mkf_wf_" + tag + "_" + std::to_string(ts) + ".csv";

    std::size_t nWind = std::min(analytical.get_excitations_per_winding().size(),
                                 spice.get_excitations_per_winding().size());
    if (nWind == 0) {
        std::cerr << "[dump_waveforms_csv] no windings in operating points; skip\n";
        return {};
    }

    // Per-winding resampled vectors for i, v, im (magnetizing — only meaningful
    // on winding 0).  Common time horizon = shortest non-empty t.back().
    double T = std::numeric_limits<double>::max();

    struct WindingDump {
        std::vector<double> aI, sI, aV, sV;
    };
    std::vector<WindingDump> wd(nWind);
    std::vector<double> aIm, sIm;
    bool haveIm = false;

    for (std::size_t w = 0; w < nWind; ++w) {
        auto& aExc = analytical.get_excitations_per_winding()[w];
        auto& sExc = spice.get_excitations_per_winding()[w];

        auto [aIt, aId] = detail::sig_wf(aExc.get_current());
        auto [sIt, sId] = detail::sig_wf(sExc.get_current());
        auto [aVt, aVd] = detail::sig_wf(aExc.get_voltage());
        auto [sVt, sVd] = detail::sig_wf(sExc.get_voltage());

        T = std::min({T,
                      detail::safe_back(aIt), detail::safe_back(sIt),
                      detail::safe_back(aVt), detail::safe_back(sVt)});

        wd[w].aI = detail::resample_or_zeros(aIt, aId, N);
        wd[w].sI = detail::resample_or_zeros(sIt, sId, N);
        wd[w].aV = detail::resample_or_zeros(aVt, aVd, N);
        wd[w].sV = detail::resample_or_zeros(sVt, sVd, N);

        if (w == 0) {
            auto [aImt, aImd] = detail::sig_wf(aExc.get_magnetizing_current());
            auto [sImt, sImd] = detail::sig_wf(sExc.get_magnetizing_current());
            if (!aImd.empty() || !sImd.empty()) {
                haveIm = true;
                T = std::min({T,
                              detail::safe_back(aImt),
                              detail::safe_back(sImt)});
                aIm = detail::resample_or_zeros(aImt, aImd, N);
                sIm = detail::resample_or_zeros(sImt, sImd, N);
            }
        }
    }
    if (T == std::numeric_limits<double>::max()) T = 1.0;

    std::ofstream out(path);
    if (!out) {
        std::cerr << "[dump_waveforms_csv] cannot open " << path << "\n";
        return {};
    }
    out << "time_s";
    for (std::size_t w = 0; w < nWind; ++w) {
        out << ",analytical_w" << w << "_i_A"
            << ",spice_w"      << w << "_i_A"
            << ",analytical_w" << w << "_v_V"
            << ",spice_w"      << w << "_v_V";
    }
    if (haveIm) {
        out << ",analytical_im_A,spice_im_A";
    }
    out << "\n";
    out.setf(std::ios::scientific);
    out.precision(8);
    for (int k = 0; k < N; ++k) {
        double tk = T * k / (N - 1);
        out << tk;
        for (std::size_t w = 0; w < nWind; ++w) {
            out << "," << wd[w].aI[k] << "," << wd[w].sI[k]
                << "," << wd[w].aV[k] << "," << wd[w].sV[k];
        }
        if (haveIm) out << "," << aIm[k] << "," << sIm[k];
        out << "\n";
    }
    out.close();
    std::cerr << "[dump_waveforms_csv] " << path << " ("
              << N << " samples, " << nWind << " winding(s)"
              << (haveIm ? ", +Im" : "")
              << ", T=" << T << " s)\n";
    return path;
}

} // namespace Testing
} // namespace OpenMagnetics
