#include <MAS.hpp>
#include "constructive_models/Magnetic.h"
#include "constructive_models/Mas.h"
#include "processors/Outputs.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/Reluctance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/Temperature.h"
#include "support/Utils.h"
#include "support/Exceptions.h"

#include <cmath>
#include <algorithm>
#include <map>

namespace OpenMagnetics {

Bobbin Magnetic::get_bobbin() {
    return get_mutable_coil().resolve_bobbin();
}

std::vector<Wire> Magnetic::get_wires() {
    return get_mutable_coil().get_wires();
}

std::vector<double> Magnetic::get_turns_ratios() const {
    return get_coil().get_turns_ratios();
}

std::vector<double> Magnetic::get_turns_ratios(MAS::Magnetic magnetic) {
    std::vector<double> turnsRatios;
    if (!magnetic.get_coil()) {
        throw std::runtime_error("Magnetic has no coil; cannot compute turns ratios");
    }
    auto coil = magnetic.get_coil().value();
    for (size_t windingIndex = 1; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        turnsRatios.push_back(double(coil.get_functional_description()[0].get_number_turns()) / coil.get_functional_description()[windingIndex].get_number_turns());
    }
    return turnsRatios;
}

Wire Magnetic::get_wire(size_t windingIndex) {
    return get_mutable_coil().resolve_wire(windingIndex);
}

std::string Magnetic::get_reference() const {
    if (get_manufacturer_info()) {
        if (get_manufacturer_info()->get_reference()) {
            return get_manufacturer_info()->get_reference().value();
        }
    }
    return "Custom component made with OpenMagnetic";
}

std::vector<double> Magnetic::get_maximum_dimensions() {
    if (!_maximumDimensions && !(has_core() && has_coil())) {
        // A datasheet-only catalogue part: its published body size. Ordered as the core's
        // {width, height, depth}, so a footprint is [0] x [2] and the height [1] -- the same
        // reading Area and Height apply to a constructed magnetic.
        if (!has_datasheet_dimensions()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Magnetic '" + get_reference() + "' has no core and coil to take its dimensions from, and its datasheet states no body size (mechanical height plus length and width, or diameter)");
        }
        const auto mechanical = get_manufacturer_info()->get_datasheet_info()->get_mechanical().value();
        double height = resolve_dimensional_values(mechanical.get_height().value());
        if (mechanical.get_length() && mechanical.get_width()) {
            _maximumDimensions = std::vector<double>{resolve_dimensional_values(mechanical.get_length().value()), height, resolve_dimensional_values(mechanical.get_width().value())};
        }
        else {
            double diameter = resolve_dimensional_values(mechanical.get_diameter().value());
            _maximumDimensions = std::vector<double>{diameter, height, diameter};
        }
    }
    if (!_maximumDimensions) {
        auto coreMaximumDimensions = get_mutable_core().get_maximum_dimensions();
        auto coilMaximumDimensions = get_mutable_coil().get_maximum_dimensions();
        _maximumDimensions = {std::max(coreMaximumDimensions[0], coilMaximumDimensions[0]),
                              std::max(coreMaximumDimensions[1], coilMaximumDimensions[1]),
                              std::max(coreMaximumDimensions[2], coilMaximumDimensions[2])};
    }
    return _maximumDimensions.value();
}

bool Magnetic::has_datasheet_dimensions() const {
    if (!get_manufacturer_info() || !get_manufacturer_info()->get_datasheet_info() || !get_manufacturer_info()->get_datasheet_info()->get_mechanical()) {
        return false;
    }
    // The generated getters return optionals BY VALUE: copy, never bind a reference to them.
    const auto mechanical = get_manufacturer_info()->get_datasheet_info()->get_mechanical().value();
    return mechanical.get_height().has_value() && ((mechanical.get_length().has_value() && mechanical.get_width().has_value()) || mechanical.get_diameter().has_value());
}

bool Magnetic::is_datasheet_coupled_inductor() const {
    if (has_core() || has_coil() || !get_manufacturer_info() || !get_manufacturer_info()->get_datasheet_info() || !get_manufacturer_info()->get_datasheet_info()->get_electrical()) {
        return false;
    }
    const auto electricals = get_manufacturer_info()->get_datasheet_info()->get_electrical().value();
    return std::any_of(electricals.begin(), electricals.end(), [](const MagneticDatasheetElectrical& entry) { return entry.get_subtype() == ElectricalSubtype::COUPLED_INDUCTOR; });
}

std::optional<MagneticDatasheetElectrical> Magnetic::get_datasheet_inductor_electrical() const {
    if (!get_manufacturer_info() || !get_manufacturer_info()->get_datasheet_info() || !get_manufacturer_info()->get_datasheet_info()->get_electrical()) {
        return std::nullopt;
    }
    std::optional<MagneticDatasheetElectrical> found;
    const auto electricals = get_manufacturer_info()->get_datasheet_info()->get_electrical().value();
    for (const auto& entry : electricals) {
        if (entry.get_subtype() != ElectricalSubtype::INDUCTOR) {
            continue;
        }
        if (found) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic '" + get_reference() + "' states several single-winding inductor configurations in its datasheet; which one applies cannot be chosen");
        }
        found = entry;
    }
    return found;
}

namespace {

// L at `current` along one measured curve, sorted by current. nullopt beyond the last point.
std::optional<double> inductance_along_curve(const std::vector<std::pair<double, double>>& curve, double current) {
    if (current <= curve.front().first) {
        return curve.front().second;
    }
    for (size_t index = 1; index < curve.size(); ++index) {
        if (current <= curve[index].first) {
            const auto& [currentLow, inductanceLow] = curve[index - 1];
            const auto& [currentHigh, inductanceHigh] = curve[index];
            return inductanceLow + (inductanceHigh - inductanceLow) * (current - currentLow) / (currentHigh - currentLow);
        }
    }
    return std::nullopt;
}

}  // namespace

std::optional<double> Magnetic::calculate_datasheet_inductance(double dcBiasCurrent, double temperature) const {
    auto electrical = get_datasheet_inductor_electrical();
    if (!electrical) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Magnetic '" + get_reference() + "' has no single-winding inductor entry in its datasheet to read an inductance from");
    }
    double current = std::fabs(dcBiasCurrent);
    auto points = electrical->get_inductance_points();
    if (!points || points->empty()) {
        if (!electrical->get_inductance()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Magnetic '" + get_reference() + "' states neither an inductance nor L(I) points in its datasheet");
        }
        return resolve_dimensional_values(electrical->get_inductance().value());
    }

    // One curve per measured temperature. A point without a temperature cannot be placed on
    // any of them, so a datasheet mixing the two is refused rather than read one way or other.
    std::map<double, std::vector<std::pair<double, double>>> curves;
    size_t withoutTemperature = 0;
    for (const auto& point : points.value()) {
        if (!point.get_current()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic '" + get_reference() + "' has a datasheet L(I) point without a current");
        }
        if (!point.get_temperature()) {
            ++withoutTemperature;
        }
        curves[point.get_temperature().value_or(0)].push_back({point.get_current().value(), point.get_inductance()});
    }
    if (withoutTemperature > 0 && withoutTemperature != points->size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic '" + get_reference() + "' mixes datasheet L(I) points with and without a temperature");
    }
    for (auto& [curveTemperature, curve] : curves) {
        std::sort(curve.begin(), curve.end());
    }

    if (curves.size() == 1 || temperature <= curves.begin()->first) {
        return inductance_along_curve(curves.begin()->second, current);
    }
    if (temperature >= curves.rbegin()->first) {
        return inductance_along_curve(curves.rbegin()->second, current);
    }
    auto upper = curves.lower_bound(temperature);
    if (upper->first == temperature) {
        return inductance_along_curve(upper->second, current);
    }
    auto lower = std::prev(upper);
    auto inductanceLow = inductance_along_curve(lower->second, current);
    auto inductanceHigh = inductance_along_curve(upper->second, current);
    if (!inductanceLow || !inductanceHigh) {
        return std::nullopt;
    }
    return inductanceLow.value() + (inductanceHigh.value() - inductanceLow.value()) * (temperature - lower->first) / (upper->first - lower->first);
}

// NOTE (code review L-7): These fits_*_dimension() functions are duplicated in Core.cpp.
// They should be consolidated into a shared utility (e.g., support/DimensionChecker.h).
// See OpenMagnetics Code Review Report for details.
bool fits_one_dimension(std::vector<double> magneticDimensions, double dimension) {
    if (magneticDimensions[2] <= dimension || magneticDimensions[1] <= dimension || magneticDimensions[0] <= dimension) {
        return true;
    }
    return false;
}

bool fits_two_dimensions(std::vector<double> magneticDimensions, double firstDimension, double secondDimension) {
    if ((magneticDimensions[2] <= firstDimension && (magneticDimensions[1] <= secondDimension || magneticDimensions[0] <= secondDimension)) ||
        (magneticDimensions[1] <= firstDimension && (magneticDimensions[2] <= secondDimension || magneticDimensions[0] <= secondDimension)) ||
        (magneticDimensions[0] <= firstDimension && (magneticDimensions[1] <= secondDimension || magneticDimensions[2] <= secondDimension))) {
        return true;
    }
    return false;
}

bool fits_three_dimensions(std::vector<double> magneticDimensions, double firstDimension, double secondDimension, double thirdDimension) {
    if ((magneticDimensions[2] <= firstDimension && magneticDimensions[1] <= secondDimension && magneticDimensions[0] <= thirdDimension) ||
        (magneticDimensions[2] <= firstDimension && magneticDimensions[1] <= thirdDimension && magneticDimensions[0] <= secondDimension) ||
        (magneticDimensions[2] <= secondDimension && magneticDimensions[1] <= firstDimension && magneticDimensions[0] <= thirdDimension) ||
        (magneticDimensions[2] <= secondDimension && magneticDimensions[1] <= thirdDimension && magneticDimensions[0] <= firstDimension) ||
        (magneticDimensions[2] <= thirdDimension && magneticDimensions[1] <= firstDimension && magneticDimensions[0] <= secondDimension) ||
        (magneticDimensions[2] <= thirdDimension && magneticDimensions[1] <= secondDimension && magneticDimensions[0] <= firstDimension)){
        return true;
    }
    return false;
}

bool Magnetic::fits(MaximumDimensions maximumDimensions, bool allowRotation) {

    auto magneticDimensions = get_maximum_dimensions();

    if (!maximumDimensions.get_depth() && !maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        return true;
    }
    else if (maximumDimensions.get_depth() && !maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        auto depth = maximumDimensions.get_depth().value();
        if (allowRotation) {
            return fits_one_dimension(magneticDimensions, depth);
        }
        else {
            return magneticDimensions[2] <= depth;
        }
    }
    else if (!maximumDimensions.get_depth() && maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        auto height = maximumDimensions.get_height().value();
        if (allowRotation) {
            return fits_one_dimension(magneticDimensions, height);
        }
        else {
            return magneticDimensions[1] <= height;
        }
    }
    else if (!maximumDimensions.get_depth() && !maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto width = maximumDimensions.get_width().value();
        if (allowRotation) {
            return fits_one_dimension(magneticDimensions, width);
        }
        else {
            return magneticDimensions[0] <= width;
        }
    }
    else if (maximumDimensions.get_depth() && maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        auto depth = maximumDimensions.get_depth().value();
        auto height = maximumDimensions.get_height().value();
        if (allowRotation) {
            return fits_two_dimensions(magneticDimensions, depth, height);
        }
        else {
            return magneticDimensions[2] <= depth && magneticDimensions[1] <= height;
        }
    }
    else if (!maximumDimensions.get_depth() && maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto width = maximumDimensions.get_width().value();
        auto height = maximumDimensions.get_height().value();
        if (allowRotation) {
            return fits_two_dimensions(magneticDimensions, width, height);
        }
        else {
            return magneticDimensions[0] <= width && magneticDimensions[1] <= height;
        }
    }
    else if (maximumDimensions.get_depth() && !maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto width = maximumDimensions.get_width().value();
        auto depth = maximumDimensions.get_depth().value();
        if (allowRotation) {
            return fits_two_dimensions(magneticDimensions, width, depth);
        }
        else {
            return magneticDimensions[2] <= depth && magneticDimensions[0] <= width;
        }
    }
    else if (maximumDimensions.get_depth() && maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto depth = maximumDimensions.get_depth().value();
        auto height = maximumDimensions.get_height().value();
        auto width = maximumDimensions.get_width().value();
        if (allowRotation) {
            return fits_three_dimensions(magneticDimensions, depth, height, width);
        }
        else {
            return magneticDimensions[2] <= depth && magneticDimensions[1] <= height && magneticDimensions[0] <= width;
        }
    }
    else {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Not sure how this happened");
    }
}

double Magnetic::calculate_saturation_current(OperatingPoint& operatingPoint,
                                              double temperature) {
    // Operating-point variant: same B_sat·N·A_e/L identity, but L is
    // taken at the supplied operating point (DC-bias-rolled-off μ).
    // Use this when comparing against I_peak from the same OP — keeps
    // both numbers on a consistent footing.
    auto magneticFluxDensitySaturation =
        get_mutable_core().get_magnetic_flux_density_saturation(temperature);
    auto numberTurns = get_mutable_coil().get_number_turns(0);
    auto effectiveArea = get_mutable_core().get_effective_area();

    OpenMagnetics::MagnetizingInductance magnetizingInductanceCalc;
    auto inductanceOutput =
        magnetizingInductanceCalc.calculate_inductance_from_number_turns_and_gapping(
            get_mutable_core(), get_mutable_coil(), &operatingPoint);
    auto inductanceNominal =
        inductanceOutput.get_magnetizing_inductance().get_nominal();
    if (!inductanceNominal.has_value() || inductanceNominal.value() <= 0) {
        throw std::runtime_error(
            "calculate_saturation_current(OP): operating-point inductance "
            "is missing or non-positive — cannot derive I_sat.");
    }
    return magneticFluxDensitySaturation * numberTurns * effectiveArea
         / inductanceNominal.value();
}

double Magnetic::calculate_saturation_current(double temperature, bool proportion) {
    // Saturation current of the as-designed magnetic:
    //
    //   I_sat = B_sat(T) · N · A_e / L_actual
    //
    // where L_actual is the inductance the wound + gapped magnetic
    // realises (NOT the bare-core inductance). Earlier this function
    // used `ReluctanceModel::get_core_reluctance(core, μ_initial)`,
    // which returns the *core's intrinsic* reluctance — for an
    // ungapped toroid that happens to equal the operating reluctance,
    // but for a gapped E-core the operating reluctance is dominated
    // by the gap (typically ~10× the core piece's reluctance). The
    // old formula therefore back-solved to a phantom L far above the
    // real inductance and over-reported I_sat by the same factor.
    //
    // Symptom this fixes: on a 48 → 12 V 60 W buck the CoreAdviser
    // picked an EP 17 with N = 6, L_actual ≈ 5.1 µH, ipeak ≈ 10.85 A.
    // The MagneticFilterSaturation filter accepted it correctly
    // (B_pk = 0.21 T, B_sat·1.2 = 0.49 T, margin 1.99×). PyOpenMagnetics'
    // calculate_saturation_current reported 11.79 A — corresponding
    // to a phantom L = 9.4 µH (the EP 17's bare-core L at μ_init,
    // before gap) — so Heaviside's I_sat / I_peak realism check saw
    // 1.086 and falsely FAILed inductor_isat_margin. The honest
    // value is 21.6 A (ratio 1.99 × — matches the filter).
    //
    // Use MagnetizingInductance::calculate_inductance_from_number_turns_and_gapping
    // which solves the gapped magnetic's circuit (core piece + each
    // gap entry + fringing factor) and reports the operating-point
    // L. That L is then inverted into the saturation current via the
    // same B_sat·N·A_e/L identity the realism gate expects.
    auto magneticFluxDensitySaturation =
        get_mutable_core().get_magnetic_flux_density_saturation(temperature, proportion);
    auto numberTurns = get_mutable_coil().get_number_turns(0);
    auto effectiveArea = get_mutable_core().get_effective_area();

    OpenMagnetics::MagnetizingInductance magnetizingInductanceCalc;
    auto inductanceOutput =
        magnetizingInductanceCalc.calculate_inductance_from_number_turns_and_gapping(
            get_mutable_core(), get_mutable_coil());
    auto inductanceNominal =
        inductanceOutput.get_magnetizing_inductance().get_nominal();
    if (!inductanceNominal.has_value() || inductanceNominal.value() <= 0) {
        throw std::runtime_error(
            "calculate_saturation_current: gapped-magnetic inductance "
            "is missing or non-positive — cannot derive I_sat.");
    }
    double inductanceActual = inductanceNominal.value();

    return magneticFluxDensitySaturation * numberTurns * effectiveArea
         / inductanceActual;
}

double Magnetic::calculate_rated_current(double temperatureRise, double ambientTemperature) {
    // Datasheet "Irms" / rated current: the DC current that heats the part by
    // `temperatureRise` K above ambient through the primary winding's ohmic loss
    // (core loss excluded — a DC bias has no flux swing). See the header for the
    // DC-vs-AC convention and references.
    if (temperatureRise <= 0) {
        throw std::invalid_argument("calculate_rated_current: temperatureRise must be positive");
    }

    // The full thermal-network model needs a wound coil to place the turn nodes; wind a
    // local copy so the call is non-mutating.
    Magnetic magnetic = *this;
    if (!magnetic.get_coil().get_turns_description()) {
        magnetic.get_mutable_coil().wind();
    }
    auto coil = magnetic.get_coil();

    // Per-turn DC winding losses for a primary current `dcCurrent` at copper temperature
    // `temperature`, via the MKF ohmic-loss model (R_dc is resolved at `temperature`). The
    // full Temperature model needs this per-turn distribution, not just a total.
    auto ohmicLossesForCurrent = [&](double dcCurrent, double temperature) {
        ProcessedWaveform processed;
        processed.set_offset(dcCurrent);
        processed.set_rms(dcCurrent);
        SignalDescriptor currentSignal;
        currentSignal.set_processed(processed);
        OperatingPointExcitation excitation;
        excitation.set_frequency(0);  // DC bias
        excitation.set_current(currentSignal);
        OperatingPoint dcOperatingPoint;
        dcOperatingPoint.set_excitations_per_winding({excitation});  // primary only; others carry no current
        return WindingOhmicLosses::calculate_ohmic_losses(coil, dcOperatingPoint, temperature);
    };

    // Steady-state hot-spot temperature for a DC current `dcCurrent` in the primary winding,
    // from the full Temperature thermal-network model fed with the per-turn DC ohmic losses.
    // R_dc rises with copper temperature, so loss and temperature are solved together by a
    // short fixed-point iteration (copper tempco is small, so it converges in a few passes).
    auto steadyStateTemperature = [&](double dcCurrent) {
        double windingTemperature = ambientTemperature;
        double maximumTemperature = ambientTemperature;
        for (size_t iteration = 0; iteration < 6; ++iteration) {
            auto losses = ohmicLossesForCurrent(dcCurrent, windingTemperature);

            TemperatureConfig config;
            config.ambientTemperature = ambientTemperature;
            config.coreLosses = 0;            // pure DC bias: no flux swing, no core loss
            config.windingLosses = losses.get_winding_losses();
            config.windingLossesOutput = losses;
            config.plotSchematic = false;     // no SVG side effect per evaluation

            auto result = Temperature(magnetic, config).calculateTemperatures();
            maximumTemperature = result.maximumTemperature;
            if (std::abs(result.averageCoilTemperature - windingTemperature) < 0.05) {
                break;
            }
            windingTemperature = result.averageCoilTemperature;
        }
        return maximumTemperature;
    };

    // The rise is monotonically increasing in current, so bisect. First grow an upper
    // bound until its rise exceeds the target.
    double lowerCurrent = 0;
    double upperCurrent = 1;
    while (steadyStateTemperature(upperCurrent) - ambientTemperature < temperatureRise) {
        upperCurrent *= 2;
        if (upperCurrent > 1e6) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR,
                "calculate_rated_current: could not reach the target temperature rise below 1 MA");
        }
    }

    for (size_t iteration = 0; iteration < 40; ++iteration) {
        double midCurrent = 0.5 * (lowerCurrent + upperCurrent);
        double rise = steadyStateTemperature(midCurrent) - ambientTemperature;
        if (std::abs(rise - temperatureRise) < 0.1) {
            return midCurrent;
        }
        if (rise < temperatureRise) {
            lowerCurrent = midCurrent;
        }
        else {
            upperCurrent = midCurrent;
        }
    }

    return 0.5 * (lowerCurrent + upperCurrent);
}

void to_file(std::filesystem::path filepath, const Magnetic & x) {
    OpenMagnetics::Mas mas;
    auto inputs = OpenMagnetics::get_defaults_inputs(x);
    std::vector<OpenMagnetics::Outputs> outputs = {OpenMagnetics::Outputs()};
    mas.set_magnetic(x);
    mas.set_inputs(inputs);
    mas.set_outputs(outputs);
    json masJson;
    to_json(masJson, mas);

    std::ofstream myfile;
    myfile.open(filepath);
    myfile << masJson;
}

} // namespace OpenMagnetics
