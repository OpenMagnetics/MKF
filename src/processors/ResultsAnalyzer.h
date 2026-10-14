#pragma once
#include "constructive_models/Mas.h"
#include "constructive_models/Magnetic.h"
#include "processors/Inputs.h"
#include "processors/Outputs.h"
#include "physical_models/MagnetizingInductance.h"
#include "json.hpp"
#include <MAS.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

using namespace MAS;

namespace OpenMagnetics {

// ABT #1855 (design decided by Alf, 2026-10-14): MagneticSimulator::simulate does not judge. It gives
// its numbers for every operating point, a saturating one or one that runs at 500 C included, and
// whether those numbers describe a part that works is said here, after the fact, from the simulated
// result alone (inputs + magnetic + outputs). Nothing of it is written into the MAS.
//
// The judgement is made by small per-topic checks. The advisers and the Sweeper reject candidates
// through the same checks (check_saturation), so an adviser and analyze_results give one verdict.

enum class ResultFindingCode {
    // B peak above the material's B_sat at the core's temperature.
    CORE_SATURATES,
    // A winding voltage whose net volt-seconds over the period exceed Inputs::maximumVoltSecondImbalance
    // of its AC flux-linkage swing (or a non-zero constant): MKF simulated its balanced part.
    VOLT_SECOND_IMBALANCE,
    // The core's temperature above a limit the core material states (Curie temperature, maximum
    // recommended operating temperature).
    CORE_TEMPERATURE_ABOVE_LIMIT,
    // The core material states no temperature limit: the core's temperature was not checked.
    CORE_TEMPERATURE_LIMIT_NOT_STATED,
    // The losses drove the core past the last temperature its material's data covers (its Curie
    // temperature, or the end of its measured mu_i(T) table): the loss solve held it there (Alf
    // 2026-10-14), with mu_r = 1 at the Curie temperature and with the last tabulated permeability at
    // the end of a table below Curie, so the core's real temperature is higher than the numbers given.
    CORE_ABOVE_MATERIAL_DATA_LIMIT,
    // The design's hottest point above a limit a winding's wire coating states (its temperature
    // rating, its insulation material's thermal class).
    WINDING_TEMPERATURE_ABOVE_LIMIT,
    // A winding's wire states no temperature limit: its temperature was not checked.
    WINDING_TEMPERATURE_LIMIT_NOT_STATED,
};

std::string to_string(ResultFindingCode code);

// One finding: which operating point (its index in the MAS, and its name when it has one), what
// (code) and MKF's message, readable as it stands.
struct ResultFinding {
    size_t operatingPointIndex;
    std::optional<std::string> operatingPointName;
    ResultFindingCode code;
    std::string message;
};

void to_json(nlohmann::json& j, const ResultFinding& finding);
void to_json(nlohmann::json& j, const std::vector<ResultFinding>& findings);

// The saturation check's verdict: B peak (times the margin) against the material's raw B_sat at the
// temperature, and MKF's message saying so.
struct SaturationVerdict {
    // Either criterion (Alf, 2026-10-14).
    bool saturates;
    // Criterion 1: B peak (times the margin) above B_sat at the temperature.
    bool fluxDensityAboveSaturation;
    // Criterion 2: the permeability model at the DC bias point at or below the vacuum's (unbiased above it).
    bool permeabilityCollapsed;
    // The permeability model's value at the DC bias point, when one was given with a non-zero bias.
    std::optional<double> biasedPermeability;
    double peakFluxDensity;
    double saturationFluxDensity;
    double temperature;
    double margin;
    std::string message;
};

class ResultsAnalyzer {
  public:
    // The saturation check, shared by analyze_results (margin 1, the core's own temperature) and by
    // the advisers' filters and the Sweeper (their margin and temperature). The core saturates
    // (Alf, 2026-10-14) when
    //   1. peakFluxDensity * margin > B_sat(temperature), the material's raw saturation flux density
    //      (Core::get_magnetic_flux_density_saturation, proportion false), OR
    //   2. a DC bias point is given (MagnetizingInductance::get_last_dc_bias_point) and the material's
    //      permeability model there is at or below the vacuum's 1 while it is above 1 without the bias:
    //      the bias collapsed the permeability to mu_r = 1 (main's SaturationException condition until
    //      ABT #1855).
    // Callers that ran the magnetizing model pass its DC point; without one only criterion 1 applies.
    static SaturationVerdict check_saturation(const CoreMaterial& coreMaterial, double peakFluxDensity, double temperature,
                                              double margin = 1.0, std::optional<DcBiasPoint> dcBiasPoint = std::nullopt);
    // |B| peak of a flux-density signal: its processed peak, else the largest magnitude of its
    // waveform. Throws when it has neither.
    static double get_peak_flux_density(const SignalDescriptor& magneticFluxDensity);

    // The per-topic checks of one simulated operating point.
    // Saturation runs the magnetizing model again at the core losses' temperature for the DC bias point of
    // criterion 2; `topology` (designRequirements.topology) decides, as in simulate, which magnetizing
    // current drives the core of a multi-winding choke.
    static std::vector<ResultFinding> analyze_saturation(const Magnetic& magnetic, const OperatingPoint& operatingPoint,
                                                         const Outputs& output, size_t operatingPointIndex,
                                                         std::optional<Topology> topology = std::nullopt);
    static std::vector<ResultFinding> analyze_volt_second_balance(const OperatingPoint& operatingPoint, size_t operatingPointIndex);
    static std::vector<ResultFinding> analyze_temperature(const Magnetic& magnetic, const OperatingPoint& operatingPoint,
                                                          const Outputs& output, size_t operatingPointIndex);

    // Every finding of a simulated MAS (as MagneticSimulator::simulate returns it: one output per
    // operating point), ordered by operating point, then saturation, volt-second balance, temperature.
    // Throws when the MAS has no output for an operating point, or an output lacks what a check reads
    // (the core losses' flux density and temperature, the temperature output's maximum).
    static std::vector<ResultFinding> analyze_results(const Mas& mas);
};

} // namespace OpenMagnetics
