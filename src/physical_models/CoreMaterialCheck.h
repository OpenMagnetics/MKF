#pragma once
// A core material judged on its own, without a design (Alf, 2026-10-28, ABT #1911).
//
// A design that embeds its own copy of a core material can only be judged today by computing the
// whole design: when the design stops for a reason that has nothing to do with the material (a
// thermal overlap, a time-axis mismatch, saturation), the material is never reached. These calls
// judge the material alone.
//
// check_core_material runs the SAME checks that refuse the material inside a design calculation
// (CoreLossesModel::get_temperature_factor, InitialPermeability::get_initial_permeability,
// resolve_core_material_remanence / _coercive_force, Core::get_magnetic_flux_density_saturation,
// resolve_core_material_heat_conductivity, CoreLossesModel::get_steinmetz_fitted_span) over the
// temperatures the material's own data covers, and reports every refusal with the error code and
// the message the design calculation raises. Nothing here re-implements a check: a check changed in
// its model changes here with it.
//
// compare_core_material_volumetric_losses prices two copies of a material with the design's
// core-loss path (CoreLosses::get_core_volumetric_losses) at points taken from the reference
// copy's own data, so a caller can tell whether a copy differs from the library's.
#include "Constants.h"
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <nlohmann/json.hpp>
#include <optional>
#include <string>
#include <vector>

namespace OpenMagnetics {

using json = nlohmann::json;

struct CoreMaterialSpan {
    double minimum;
    double maximum;
};

// One way a design calculation refuses the material.
struct CoreMaterialProblem {
    ErrorCode code;              // the code the design calculation throws (MATERIAL_INVALID_PROPERTY, INVALID_LOSS_DATA, ...)
    std::string property;        // the MAS field the problem lies in: "volumetricLosses", "remanence", ...
    std::string message;         // what is wrong, over what range, and what to do
    std::string engineMessage;   // the exception the design calculation raises there, verbatim (what())
    std::optional<CoreMaterialSpan> temperature;  // C: where inside the judged span the data is refused
    std::optional<CoreMaterialSpan> frequency;    // Hz: the loss range concerned (Steinmetz range)
};

// What the material's own data covers.
struct CoreMaterialCoverage {
    // C: every temperature the material tabulates (initial permeability, saturation, remanence,
    // coercive force, volumetric-loss points).
    CoreMaterialSpan temperature;
    // Hz: the span its Steinmetz ranges are fitted over, and its volumetric-loss points.
    std::optional<CoreMaterialSpan> lossFrequency;
    // T (peak): its volumetric-loss points. A fitted model alone states no flux density span.
    std::optional<CoreMaterialSpan> lossMagneticFluxDensity;
};

struct CoreMaterialCheck {
    std::string material;
    bool valid;                               // no problem found
    std::vector<CoreMaterialProblem> problems;
    CoreMaterialCoverage coverage;
    std::vector<std::string> notChecked;      // "<property>: <why it was not judged>"
};

// One volumetric-loss evaluation: the losses, or the refusal (code and message) of the design path.
struct CoreMaterialLossValue {
    std::optional<double> volumetricLosses;  // W/m3
    std::optional<ErrorCode> code;
    std::optional<std::string> message;
};

struct CoreMaterialLossComparisonPoint {
    double frequency;                 // Hz
    double magneticFluxDensityPeak;   // T
    double temperature;               // C
    CoreMaterialLossValue reference;
    CoreMaterialLossValue other;
    std::optional<double> relativeDifference;  // (other - reference) / reference, when both computed
};

struct CoreMaterialLossComparison {
    std::string reference;
    std::string other;
    std::string pointsSource;  // how the points were taken from the reference's data
    std::vector<CoreMaterialLossComparisonPoint> points;
    std::optional<double> maximumAbsoluteRelativeDifference;  // over the points both copies computed
    size_t pointsRefusedByReference = 0;
    size_t pointsRefusedByOther = 0;
};

// Judge a core material alone. Throws MaterialException(MATERIAL_DATA_MISSING) when the material
// tabulates nothing at any temperature (then there is no span to judge it over).
CoreMaterialCheck check_core_material(const CoreMaterial& material);

// Price both copies at the reference copy's own points: its volumetric-loss points when it carries
// them; otherwise, at each temperature its saturation flux density is tabulated at (below its Curie
// point), sinusoidal flux densities of 1/2, 1/4 and 1/8 of that saturation flux density, at the
// geometric centre of each of its Steinmetz ranges. Throws MaterialException(MATERIAL_DATA_MISSING)
// when the reference has neither loss points nor Steinmetz ranges, or no saturation data.
CoreMaterialLossComparison compare_core_material_volumetric_losses(const CoreMaterial& reference, const CoreMaterial& other);

void to_json(json& j, const CoreMaterialSpan& span);
void to_json(json& j, const CoreMaterialProblem& problem);
void to_json(json& j, const CoreMaterialCoverage& coverage);
void to_json(json& j, const CoreMaterialCheck& check);
void to_json(json& j, const CoreMaterialLossValue& value);
void to_json(json& j, const CoreMaterialLossComparisonPoint& point);
void to_json(json& j, const CoreMaterialLossComparison& comparison);

} // namespace OpenMagnetics
