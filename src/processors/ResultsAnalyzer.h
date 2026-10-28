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
    // of its AC flux-linkage swing (or a non-zero constant), and that the winding's own I*R drop does not
    // explain (analyze_volt_second_balance): MKF simulated its balanced part.
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
    // One per winding, naming the winding, its drawn centreline radius and the mandrel. ABT #2012: the
    // same finding names the winding's terminal leads planned for edgewise bends below that mandrel.
    // A property of the magnetic, not of an operating point: the finding carries no operating point.
    EDGEWISE_BEND_BELOW_COATING_MANDREL,
    // Alf, 2026-10-24: a terminal lead of the real-winding layout is drawn across another turn (its
    // tangent straight crosses that turn's corner closer than their half heights). The design is drawn
    // as routed; the copper collides there. One per lead, naming the winding, the lead, the turn it hits
    // (the worst of them) and the overlap. A property of the magnetic: no operating point.
    LEAD_COLLIDES_WITH_TURNS,
    // Alf, 2026-10-25: a terminal lead of the real-winding layout leaves tangentially on a route MKF does
    // not model in full (a dragback lane riding over its face, a lateral column, a pin run, a revolution
    // that changes layer, a one-turn conductor: TangentLeadCondition). It is drawn on its best route; one
    // per lead, naming the winding, the lead and each condition with its numbers. No operating point.
    LEAD_ROUTE_NOT_MODELLED,
    // Alf, 2026-10-25 (ABT #2071): turn stations of the drawn coil overlap one another (a wind that does
    // not fit its window still produces its turns). One for the coil, naming the overlapping turns, how
    // many places they overlap and the worst overlap. No operating point.
    TURNS_OVERLAP,
    // Alf, 2026-10-25: two terminal leads of the real-winding layout are drawn through each other. One
    // per pair, naming both and their overlap. No operating point.
    LEADS_COLLIDE,
    // Alf, 2026-10-25 (ABT #2071): turn stations of the drawn coil lie outside the winding window they
    // are wound in (a wind that does not fit its window is laid out as wound, overflowing). One per
    // winding, naming each turn outside, how far past the window it reaches and past which side.
    // No operating point.
    TURNS_OUTSIDE_WINDOW,
    // Alf, 2026-10-27 (ABT #2109): whether the turns overlap one another or lie outside the winding
    // window cannot be judged, because a turn stores no cross-section and its wire does not give one
    // either (Coil::turn_outlines_unknown_reason). One for the coil, with the reason, in place of
    // TURNS_OVERLAP, LEAD_COLLIDES_WITH_TURNS, LEADS_COLLIDE and TURNS_OUTSIDE_WINDOW; every other
    // finding is still returned. No operating point.
    TURN_OUTLINES_NOT_EVALUABLE,
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
    // The numbers behind the message, for whoever has to look into it (a solver's iterates); none for most findings.
    std::optional<std::string> details;
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
    // Criterion 4 (Alf, 2026-10-24): the magnetizing model's DC solve found no operating point
    // (DcBiasPoint::noDcOperatingPoint) and took the saturated mu_r = 1.
    bool noDcOperatingPoint;
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
    // Volt-second balance of each winding voltage. A voltage balanced as given is never flagged (ABT #2112:
    // the drop is not taken out of it, which gave an ideal V = 0 at a DC current a mean of -R*I). ABT #2002:
    // with `windingDcResistances` (each winding's R_dc at `resistanceTemperature`) an unbalanced voltage of a
    // winding with a current is judged on v(t) - R_dc i(t), its own resistive drop taken out
    // (Inputs::remove_resistive_drop), when that leaves a smaller mean: no finding when what remains is
    // balanced, the remaining imbalance otherwise. A drop that does not shrink the mean does not explain
    // it, and the voltage is judged as given; so is a voltage without resistances or a current.
    static std::vector<ResultFinding> analyze_volt_second_balance(const OperatingPoint& operatingPoint, size_t operatingPointIndex,
                                                                  std::optional<std::vector<double>> windingDcResistances = std::nullopt,
                                                                  std::optional<double> resistanceTemperature = std::nullopt);
    // The same for a simulated point: the resistances come from the magnetic's coil when it is wound
    // (Inputs::get_winding_dc_resistances), at the temperature of the output's winding losses, or the
    // operating point's temperature when the output has none.
    static std::vector<ResultFinding> analyze_volt_second_balance(const Magnetic& magnetic, const OperatingPoint& operatingPoint,
                                                                  const Outputs& output, size_t operatingPointIndex);
    static std::vector<ResultFinding> analyze_temperature(const Magnetic& magnetic, const OperatingPoint& operatingPoint,
                                                          const Outputs& output, size_t operatingPointIndex);
    // ABT #1876: for an output WITHOUT a temperature. simulate gives none to a point with no thermal operating
    // point, and none to any point in fast mode. ABT #2075: simulate records the reason of the first case on the
    // point's core losses (MagneticSimulator::get_no_thermal_operating_point_reason), and this returns the
    // NO_THERMAL_OPERATING_POINT finding with that message, without solving the point again. An output with no
    // reason recorded is not this case (fast mode, or an output edited after simulate) and this THROWS: such an
    // output cannot be analyzed.
    static std::vector<ResultFinding> analyze_missing_temperature(const OperatingPoint& operatingPoint, const Outputs& output,
                                                                  size_t operatingPointIndex);

    // ABT #2052: for an output WITHOUT core losses. simulate leaves them (and the point's inductance and
    // temperature) unset when the magnetizing model, at the point's ambient temperature, cannot read its flux
    // density on the material's magnetisation curve (SaturationException / MagnetisationCurveExhaustedException:
    // past the curve's end). This runs that model again and returns the CORE_SATURATES finding with its reason;
    // when the model DOES read the point, the missing core losses are not this case and this THROWS.
    static std::vector<ResultFinding> analyze_flux_density_past_magnetisation_curve(const Magnetic& magnetic,
                                                                                   const OperatingPoint& operatingPoint,
                                                                                   size_t operatingPointIndex,
                                                                                   std::optional<Topology> topology = std::nullopt);

    // ABT #2004: one EDGEWISE_BEND_BELOW_COATING_MANDREL finding per winding whose real-winding layout
    // bends its rectangular wire edgewise tighter than the IEC 60317-0-2 Table 6 mandrel, from the same
    // judgement real winding makes (Coil::judge_corners_against_wire_flexibility), naming also the
    // winding's terminal leads planned below it (Coil::judge_leads_against_wire_flexibility). Nothing under ideal
    // winding (the flag off judges no corner). A round or flatwise bend below its limit is refused by
    // real winding and is never a finding: this throws UnwindableBendException for it, as wind does.
    static std::vector<ResultFinding> analyze_winding_bends(const Magnetic& magnetic);

    // One LEAD_COLLIDES_WITH_TURNS finding per terminal lead the real-winding layout draws across
    // another turn (Coil::judge_lead_collisions). Nothing under ideal winding.
    static std::vector<ResultFinding> analyze_lead_collisions(const Magnetic& magnetic);

    // One LEAD_ROUTE_NOT_MODELLED finding per terminal lead drawn on a tangent route MKF does not model
    // in full (Coil::judge_unmodelled_tangent_leads). Nothing under ideal winding.
    static std::vector<ResultFinding> analyze_unmodelled_leads(const Magnetic& magnetic);

    // The copper the drawn design overlaps (Coil::judge_copper_overlaps, the same overlaps the Painter
    // draws in red): one TURNS_OVERLAP finding when turn stations overlap, one LEAD_COLLIDES_WITH_TURNS
    // per radial terminal lead drawn across turns (a tangent lead's are analyze_lead_collisions'), and
    // one LEADS_COLLIDE per pair of terminal leads drawn through each other.
    static std::vector<ResultFinding> analyze_copper_overlaps(const Magnetic& magnetic);

    // One TURNS_OUTSIDE_WINDOW finding per winding with turn stations outside their winding window
    // (Coil::judge_turns_outside_window, the same copper the Painter draws in red).
    static std::vector<ResultFinding> analyze_turns_outside_window(const Magnetic& magnetic);

    // Every finding of a simulated MAS (as MagneticSimulator::simulate returns it: one output per
    // operating point), ordered by operating point, then saturation, volt-second balance, temperature;
    // then the findings about the magnetic itself (analyze_winding_bends, analyze_lead_collisions,
    // analyze_unmodelled_leads, analyze_copper_overlaps, analyze_turns_outside_window).
    // An output without a temperature goes to analyze_missing_temperature instead of the temperature checks
    // (ABT #1876): NO_THERMAL_OPERATING_POINT, or a throw for a fast-mode output. An output without core losses
    // goes to analyze_flux_density_past_magnetisation_curve instead of the saturation and temperature checks
    // (ABT #2052): CORE_SATURATES, or a throw.
    // Throws when the MAS has no output for an operating point, or an output lacks what a check reads
    // (the core losses' flux density and temperature).
    static std::vector<ResultFinding> analyze_results(const Mas& mas);
};

} // namespace OpenMagnetics
