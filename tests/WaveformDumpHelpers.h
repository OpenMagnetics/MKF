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

    // Per-winding resampled vectors and a common time horizon.
    double T = std::numeric_limits<double>::max();
    std::vector<std::vector<double>> aI(nWind), sI(nWind);
    for (std::size_t w = 0; w < nWind; ++w) {
        auto aT = ptp_current_time(analytical, w);
        auto aD = ptp_current(analytical, w);
        auto sT = ptp_current_time(spice, w);
        auto sD = ptp_current(spice, w);
        if (aT.empty() || aD.empty() || sT.empty() || sD.empty()) {
            std::cerr << "[dump_waveforms_csv] winding " << w
                      << " missing current waveform; skip column\n";
            aI[w].assign(N, 0.0);
            sI[w].assign(N, 0.0);
            continue;
        }
        T = std::min({T, aT.back(), sT.back()});
        aI[w] = ptp_interp(aT, aD, N);
        sI[w] = ptp_interp(sT, sD, N);
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
            << ",spice_w" << w << "_i_A";
    }
    out << "\n";
    out.setf(std::ios::scientific);
    out.precision(8);
    for (int k = 0; k < N; ++k) {
        double tk = T * k / (N - 1);
        out << tk;
        for (std::size_t w = 0; w < nWind; ++w) {
            out << "," << aI[w][k] << "," << sI[w][k];
        }
        out << "\n";
    }
    out.close();
    std::cerr << "[dump_waveforms_csv] " << path << " ("
              << N << " samples, " << nWind << " winding(s), T="
              << T << " s)\n";
    return path;
}

} // namespace Testing
} // namespace OpenMagnetics
