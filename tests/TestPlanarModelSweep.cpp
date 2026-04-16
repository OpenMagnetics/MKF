// Combinatorial model sweep for planar winding losses.
//
// For each of the 6 planar test MAS files, iterate every combination of
// (MagneticFieldStrengthModels × FringingEffectModels × SkinEffectModels ×
// ProximityEffectModels), compute winding losses at each frequency, compare
// against a FEM reference (planar_fem_results.tsv produced by
// FEM/scripts/benchmark_planar_fem.py), and report which combination
// minimizes error across the sweep.
//
// Tagged [benchmark-planar] so it is opt-in:
//     ./MKF_tests "[benchmark-planar]"
#include <source_location>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include "Models.h"
#include "support/Settings.h"
#include "physical_models/WindingLosses.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/Inputs.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <magic_enum.hpp>

using namespace MAS;
using namespace OpenMagnetics;

namespace TestPlanarModelSweep {

struct PlanarCase {
    std::string name;
    std::string file;
    double temperature;
};

// Same 6 cases as FEM/scripts/benchmark_planar_fem.py.
static const std::vector<PlanarCase> kPlanarCases = {
    {"One_Turn_No_Fringing",
     "Test_Winding_Losses_One_Turn_Planar_Sinusoidal_No_Fringing.json", 22.0},
    {"One_Turn_Fringing",
     "Test_Winding_Losses_One_Turn_Planar_Sinusoidal_Fringing.json", 22.0},
    {"Sixteen_Turns_No_Fringing",
     "Test_Winding_Losses_Sixteen_Turns_Planar_Sinusoidal_No_Fringing.json", 22.0},
    {"Sixteen_Turns_Fringing_Close",
     "Test_Winding_Losses_Sixteen_Turns_Planar_Sinusoidal_Fringing_Close.json", 100.0},
    {"Sixteen_Turns_Fringing_Far",
     "Test_Winding_Losses_Sixteen_Turns_Planar_Sinusoidal_Fringing_Far.json", 100.0},
    {"Sixteen_Turns_Interleaving",
     "Test_Winding_Losses_Sixteen_Turns_Planar_Sinusoidal_No_Fringing_Interleaving.json", 100.0},
};

// (case_name, freq_hz) -> FEM P_cu [W]
using FemMap = std::map<std::pair<std::string, double>, double>;

FemMap load_fem_reference(const std::filesystem::path& tsv) {
    FemMap out;
    std::ifstream f(tsv);
    if (!f) {
        std::cerr << "[sweep] cannot open FEM reference: " << tsv << "\n";
        return out;
    }
    std::string line;
    std::getline(f, line);  // header
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::string name;
        double freq = 0, pcu = 0;
        size_t t1 = line.find('\t');
        size_t t2 = line.find('\t', t1 + 1);
        if (t1 == std::string::npos || t2 == std::string::npos) continue;
        name = line.substr(0, t1);
        freq = std::stod(line.substr(t1 + 1, t2 - t1 - 1));
        pcu = std::stod(line.substr(t2 + 1));
        out[{name, freq}] = pcu;
    }
    return out;
}

double compute_pcu(OpenMagnetics::Magnetic& magnetic,
                   MAS::OperatingPoint& operatingPoint,
                   double temperature) {
    auto losses = OpenMagnetics::WindingLosses().calculate_losses(magnetic, operatingPoint, temperature);
    return losses.get_winding_losses();
}

struct ComboResult {
    MagneticFieldStrengthModels field;
    MagneticFieldStrengthFringingEffectModels fringing;
    WindingSkinEffectLossesModels skin;
    WindingProximityEffectLossesModels prox;
    // Per-frequency ratios (MKF / FEM), aggregated below.
    int n_pts;
    double geom_error;       // mean |log(MKF/FEM)|
    double max_error;        // max |MKF/FEM − 1|
};

TEST_CASE("Planar model combinatorial sweep", "[benchmark-planar][!benchmark]") {
    auto fem_tsv = std::filesystem::path("/home/alf/OpenMagnetics/FEM/scratch/planar_fem_results.tsv");
    auto fem = load_fem_reference(fem_tsv);
    REQUIRE(!fem.empty());

    auto outPath = std::filesystem::path{ std::source_location::current().file_name() }
                       .parent_path().append("..").append("output").append("planar_model_sweep.tsv");
    std::filesystem::create_directories(outPath.parent_path());
    std::ofstream out(outPath);
    out << "case\tfield\tfringing\tskin\tprox\tn_pts\tgeom_error\tmax_error\n";

    // Aggregated across all cases
    std::map<std::tuple<int,int,int,int>, std::vector<double>> combo_log_ratios;
    std::map<std::tuple<int,int,int,int>, double> combo_max_err;

    const auto field_vals = magic_enum::enum_values<MagneticFieldStrengthModels>();
    const auto fringing_vals = magic_enum::enum_values<MagneticFieldStrengthFringingEffectModels>();
    const auto skin_vals = magic_enum::enum_values<WindingSkinEffectLossesModels>();
    const auto prox_vals = magic_enum::enum_values<WindingProximityEffectLossesModels>();

    MagnetizingInductance magIndModel("ZHANG");

    for (const auto& pc : kPlanarCases) {
        auto path = OpenMagneticsTesting::get_test_data_path(
            std::source_location::current(), pc.file);
        auto mas = OpenMagneticsTesting::mas_loader(path);
        auto magnetic = mas.get_magnetic();
        auto inputs = mas.get_inputs();

        // Collect freqs from FEM reference
        std::vector<double> freqs;
        for (const auto& kv : fem) if (kv.first.first == pc.name) freqs.push_back(kv.first.second);
        std::sort(freqs.begin(), freqs.end());
        if (freqs.empty()) {
            std::cerr << "[sweep] no FEM data for " << pc.name << "\n";
            continue;
        }
        std::cerr << "[sweep] " << pc.name << " — " << freqs.size() << " frequencies\n";

        for (auto field : field_vals) for (auto fr : fringing_vals)
        for (auto sk : skin_vals) for (auto px : prox_vals) {
            settings.reset();
            clear_databases();
            settings.set_magnetic_field_strength_model(field);
            settings.set_magnetic_field_strength_fringing_effect_model(fr);
            settings.set_winding_skin_effect_losses_model(sk);
            settings.set_winding_proximity_effect_losses_model(px);
            settings.set_magnetic_field_mirroring_dimension(1);
            settings.set_magnetic_field_include_fringing(true);

            int n_pts = 0;
            double sum_log2 = 0;
            double maxe = 0;
            for (double freq : freqs) {
                double fem_pcu = fem.at({pc.name, freq});
                if (fem_pcu <= 0) continue;  // skip DC / degenerate
                OperatingPoint op = inputs.get_operating_point(0);
                OpenMagnetics::Inputs::scale_time_to_frequency(op, freq, true);
                double Lmag = OpenMagnetics::resolve_dimensional_values(
                    magIndModel.calculate_inductance_from_number_turns_and_gapping(
                        magnetic.get_core(), magnetic.get_coil(), &op)
                        .get_magnetizing_inductance());
                op = OpenMagnetics::Inputs::process_operating_point(op, Lmag);
                double pcu = 0;
                try {
                    pcu = compute_pcu(magnetic, op, pc.temperature);
                } catch (...) {
                    pcu = std::numeric_limits<double>::quiet_NaN();
                }
                if (!std::isfinite(pcu) || pcu <= 0) {
                    maxe = std::numeric_limits<double>::infinity();
                    continue;
                }
                double ratio = pcu / fem_pcu;
                sum_log2 += std::log(ratio) * std::log(ratio);
                maxe = std::max(maxe, std::fabs(ratio - 1.0));
                ++n_pts;
            }
            double geom_err = n_pts ? std::exp(std::sqrt(sum_log2 / n_pts)) - 1.0
                                    : std::numeric_limits<double>::infinity();
            out << pc.name << "\t" << magic_enum::enum_name(field) << "\t"
                << magic_enum::enum_name(fr) << "\t" << magic_enum::enum_name(sk) << "\t"
                << magic_enum::enum_name(px) << "\t" << n_pts << "\t"
                << std::setprecision(6) << geom_err << "\t" << maxe << "\n";
            out.flush();

            auto key = std::make_tuple(static_cast<int>(field), static_cast<int>(fr),
                                        static_cast<int>(sk), static_cast<int>(px));
            if (n_pts > 0) {
                combo_log_ratios[key].push_back(std::log(geom_err + 1.0));
                combo_max_err[key] = std::max(combo_max_err[key], maxe);
            }
        }
    }
    settings.reset();

    // Rank combinations by mean-across-cases geometric error
    struct Ranked {
        std::tuple<int,int,int,int> key;
        double mean_geom_err;
        double worst_max_err;
        int n_cases;
    };
    std::vector<Ranked> ranked;
    for (auto& kv : combo_log_ratios) {
        double s = 0;
        for (double v : kv.second) s += v;
        ranked.push_back({kv.first, std::exp(s / kv.second.size()) - 1.0,
                          combo_max_err[kv.first], (int)kv.second.size()});
    }
    std::sort(ranked.begin(), ranked.end(),
              [](const Ranked& a, const Ranked& b) {
                  if (a.n_cases != b.n_cases) return a.n_cases > b.n_cases;
                  return a.mean_geom_err < b.mean_geom_err;
              });

    std::cerr << "\n===== TOP 15 combinations (lowest mean-geom-error across cases) =====\n";
    std::cerr << std::left << std::setw(4) << "#"
              << std::setw(20) << "field"
              << std::setw(15) << "fringing"
              << std::setw(18) << "skin"
              << std::setw(14) << "prox"
              << std::setw(10) << "cases"
              << std::setw(14) << "geom_err"
              << "max_err\n";
    for (size_t i = 0; i < std::min<size_t>(15, ranked.size()); ++i) {
        auto& r = ranked[i];
        auto field = static_cast<MagneticFieldStrengthModels>(std::get<0>(r.key));
        auto fr = static_cast<MagneticFieldStrengthFringingEffectModels>(std::get<1>(r.key));
        auto sk = static_cast<WindingSkinEffectLossesModels>(std::get<2>(r.key));
        auto px = static_cast<WindingProximityEffectLossesModels>(std::get<3>(r.key));
        std::cerr << std::left << std::setw(4) << i + 1
                  << std::setw(20) << magic_enum::enum_name(field)
                  << std::setw(15) << magic_enum::enum_name(fr)
                  << std::setw(18) << magic_enum::enum_name(sk)
                  << std::setw(14) << magic_enum::enum_name(px)
                  << std::setw(10) << r.n_cases
                  << std::setw(14) << std::setprecision(4) << r.mean_geom_err
                  << r.worst_max_err << "\n";
    }
    std::cerr << "\nDetailed per-case TSV: " << outPath << "\n";
}

}  // namespace TestPlanarModelSweep
