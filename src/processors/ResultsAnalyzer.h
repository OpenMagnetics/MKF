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
    // ABT #1876 (Alf, 2026-10-14): the operating point has no thermal operating point (the thermal network
    // does not converge on its losses, or its winding runs away with the core held at its data limit).
    // simulate gave it its losses at ambient and no temperature; the message is MKF's, with the loop history.
    NO_THERMAL_OPERATING_POINT,
    // ABT #2004 (Alf, 2026-10-22): a winding of rectangular wire is drawn (real winding) around an
    // EDGEWISE bend tighter than the IEC 60317-0-2 Table 6 mandrel its coating was qualified at. The
    // coil is drawn as laid; the coating's integrity at that bend is outside the IEC qualification.
    // One per winding, naming the winding, its drawn centreline radius and the mandrel. A property of
    // the magnetic, not of an operating point: the finding carries no operating point.
    EDGEWISE_BEND_BELOW_COATING_MANDREL,
};

std::string to_string(ResultFindingCode code);

// One finding: which operating point (its index in the MAS, and its name when it has one), what
// (code) and MKF's message, readable as it stands. A finding about the magnetic itself rather than
// one operating point (EDGEWISE_BEND_BELOW_COATING_MANDREL) carries no operating point, and its JSON
// has no operatingPointIndex.
struct ResultFinding {
    std::optional<size_t> operatingPointIndex;
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
    // Criterion 3 (ABT #1942): the B was read on the magnetisation curve (DcBiasPoint::peakMagneticFieldStrength)
    // and its peak field strength (times the margin) is above the material's H_sat at the temperature.
    bool fieldStrengthAboveSaturation;
    // The curve's peak field strength and the material's H_sat it was held against, when criterion 3 applied.
    std::optional<double> peakFieldStrength;
    std::optional<double> saturationFieldStrength;
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
    //      ABT #1855), OR
    //   3. the DC bias point carries a peak field strength (the magnetizing model read B on the material's
    //      magnetisation curve, a unipolar DCM magnetizing current, ABT #1942) and peak H * margin >
    //      H_sat(temperature) (Core::get_magnetic_field_strength_saturation): on the curve B rises only on
    //      the vacuum slope past the knee, so B_sat alone misses a core driven past it (Alf, 2026-10-18).
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
    // ABT #1876: for an output WITHOUT a temperature. The MAS cannot say why a point has none (simulate gives
    // none to a point with no thermal operating point, and none to any point in fast mode), so this solves the
    // point's losses at its operating temperature again (MagneticSimulator::calculate_core_losses, the solve
    // simulate ran) and returns the NO_THERMAL_OPERATING_POINT finding with the solve's message when it finds
    // none. When the solve DOES find a thermal operating point, the missing temperature is not this case (fast
    // mode, or an output edited after simulate) and this THROWS: such an output cannot be analyzed.
    static std::vector<ResultFinding> analyze_missing_temperature(const Magnetic& magnetic, const OperatingPoint& operatingPoint,
                                                                  size_t operatingPointIndex,
                                                                  std::optional<Topology> topology = std::nullopt);

    // ABT #2004: one EDGEWISE_BEND_BELOW_COATING_MANDREL finding per winding whose real-winding layout
    // bends its rectangular wire edgewise tighter than the IEC 60317-0-2 Table 6 mandrel, from the same
    // judgement real winding makes (Coil::judge_corners_against_wire_flexibility). Nothing under ideal
    // winding (the flag off judges no corner). A round or flatwise bend below its limit is refused by
    // real winding and is never a finding: this throws UnwindableBendException for it, as wind does.
    static std::vector<ResultFinding> analyze_winding_bends(const Magnetic& magnetic);

    // Every finding of a simulated MAS (as MagneticSimulator::simulate returns it: one output per
    // operating point), ordered by operating point, then saturation, volt-second balance, temperature;
    // then the findings about the magnetic itself (analyze_winding_bends).
    // An output without a temperature goes to analyze_missing_temperature instead of the temperature checks
    // (ABT #1876): NO_THERMAL_OPERATING_POINT, or a throw for a fast-mode output.
    // Throws when the MAS has no output for an operating point, or an output lacks what a check reads
    // (the core losses' flux density and temperature).
    static std::vector<ResultFinding> analyze_results(const Mas& mas);
};

} // namespace OpenMagnetics
