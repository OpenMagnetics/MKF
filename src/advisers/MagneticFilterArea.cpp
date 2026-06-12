#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Insulation.h"
#include "physical_models/CoreLosses.h"
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numbers>
#include <string>

namespace OpenMagnetics {

// Phase 5: extracted from MagneticFilter.cpp.
// This translation unit owns the area-based filters:
//   - MagneticFilterAreaProduct (the heavyweight first-stage core sieve)
//   - MagneticFilterAreaNoParallels (per-section geometric fit, no parallels)
//   - MagneticFilterAreaWithParallels (per-section geometric fit with parallels)
// Shared helpers (DUMMY_SENTINEL_NAME usage, settings, defaults, constants)
// come transitively via MagneticFilter.h + support/Utils.h.

MagneticFilterAreaProduct::MagneticFilterAreaProduct(Inputs inputs) {
    double frequencyReference = 100000;
    SignalDescriptor magneticFluxDensity;
    Processed processed;
    _operatingPointExcitation.set_frequency(frequencyReference);
    processed.set_label(WaveformLabel::SINUSOIDAL);
    processed.set_offset(0);
    processed.set_peak(_magneticFluxDensityReference);
    processed.set_peak_to_peak(2 * _magneticFluxDensityReference);
    magneticFluxDensity.set_processed(processed);
    _operatingPointExcitation.set_magnetic_flux_density(magneticFluxDensity);
    _coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}}));
    _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));

    if (settings.get_core_adviser_include_margin() && inputs.get_design_requirements().get_insulation()) {
        auto clearanceAndCreepageDistance = InsulationCoordinator().calculate_creepage_distance(inputs, true);
        _averageMarginInWindingWindow = clearanceAndCreepageDistance;
    }

    double primaryAreaFactor = 1;
    if (inputs.get_design_requirements().get_turns_ratios().size() > 0) {
        primaryAreaFactor = 0.5;
    }

    _areaProductRequiredPreCalculations.clear();
    // Phase 6 (perf): cache by const-ref to avoid OperatingPoint deep copies.
    const auto& operatingPoints = inputs.get_operating_points();
    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        auto excitation = Inputs::get_primary_excitation(operatingPoints[operatingPointIndex]);
        auto voltageWaveform = excitation.get_voltage().value().get_waveform().value();
        auto currentWaveform = excitation.get_current().value().get_waveform().value();
        double frequency = excitation.get_frequency();
        if (voltageWaveform.get_data().size() != currentWaveform.get_data().size()) {
            size_t maximumNumberPoints = constants.numberPointsSampledWaveforms;
            maximumNumberPoints = std::max(maximumNumberPoints, voltageWaveform.get_data().size());
            maximumNumberPoints = std::max(maximumNumberPoints, currentWaveform.get_data().size());
            voltageWaveform = Inputs::calculate_sampled_waveform(voltageWaveform, frequency, maximumNumberPoints);
            currentWaveform = Inputs::calculate_sampled_waveform(currentWaveform, frequency, maximumNumberPoints);
        }

        std::vector<double> voltageWaveformData = voltageWaveform.get_data();
        std::vector<double> currentWaveformData = currentWaveform.get_data();

        double powerMean = 0;
        for (size_t i = 0; i < voltageWaveformData.size(); ++i)
        {
            powerMean += fabs(voltageWaveformData[i] * currentWaveformData[i]);
        }
        powerMean /= voltageWaveformData.size();

        double switchingFrequency = Inputs::get_switching_frequency(excitation);
        double preCalculation = 0;

        if (inputs.get_wiring_technology() == WiringTechnology::PRINTED) {
            preCalculation = powerMean / (primaryAreaFactor * 2 * switchingFrequency * defaults.maximumCurrentDensityPlanar);
        }
        else {
            preCalculation = powerMean / (primaryAreaFactor * 2 * switchingFrequency * defaults.maximumCurrentDensity);
        }

        if (preCalculation > 1) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "maximumAreaProductRequired cannot be larger than 1 (probably)");
        }
        _areaProductRequiredPreCalculations.push_back(preCalculation);
        if (std::isinf(_areaProductRequiredPreCalculations.back()) || _areaProductRequiredPreCalculations.back() == 0) {
            std::cerr << "powerMean: " << powerMean << std::endl;
            std::cerr << "operatingPointIndex: " << operatingPointIndex << std::endl;
            std::cerr << "primaryAreaFactor: " << primaryAreaFactor << std::endl;
            std::cerr << "switchingFrequency: " << switchingFrequency << std::endl;
            std::cerr << "_areaProductRequiredPreCalculations.back(): " << _areaProductRequiredPreCalculations.back() << std::endl;
            throw NaNResultException("_areaProductRequiredPreCalculations cannot be 0 or NaN");
        }
    }
}

std::pair<bool, double> MagneticFilterAreaProduct::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    const auto& core = magnetic->get_core();

    double bobbinFillingFactor;
    if (core.get_winding_windows().size() == 0)
        return {false, 0.0};
    auto windingWindow = core.get_winding_windows()[0];
    auto windingColumn = core.get_columns()[0];
    if (inputs->get_wiring_technology() == WiringTechnology::PRINTED) {
        bobbinFillingFactor = 1;
    }
    else if (!_bobbinFillingFactors.contains(core.get_shape_name())) {
        if (core.get_functional_description().get_type() != CoreType::TOROIDAL) {
            bobbinFillingFactor = Bobbin::get_filling_factor(windingWindow.get_width().value(), core.get_winding_windows()[0].get_height().value());
        }
        else {
            // For toroids: calculate realistic filling factor based on geometry
            // The inner circumference is smaller than outer, limiting wire packing
            // Manual winding is less efficient than bobbin-based winding
            // Typical toroid fill factors are 0.55-0.70 depending on geometry
            if (windingWindow.get_radial_height()) {
                double radialHeight = windingWindow.get_radial_height().value();
                double outerRadius = core.get_width() / 2;
                double innerRadius = outerRadius - radialHeight;
                // Ratio of inner to outer circumference limits packing
                double circumferenceRatio = (innerRadius > 0) ? (innerRadius / outerRadius) : 0.5;
                // Base filling factor ~0.55, adjusted up to ~0.70 for favorable geometry
                bobbinFillingFactor = 0.55 + 0.15 * circumferenceRatio;
            } else {
                bobbinFillingFactor = 0.6;  // Default for toroids without radial height
            }
        }
        _bobbinFillingFactors[core.get_shape_name()] = bobbinFillingFactor;
    }
    else {
        bobbinFillingFactor = _bobbinFillingFactors[core.get_shape_name()];
    }
    double windingWindowArea = windingWindow.get_area().value();
    if (_averageMarginInWindingWindow > 0) {
        if (core.get_functional_description().get_type() != CoreType::TOROIDAL) {
            if (windingWindow.get_height().value() > windingWindow.get_width().value()) {
                windingWindowArea -= windingWindow.get_width().value() * _averageMarginInWindingWindow;
            }
            else {
                windingWindowArea -= windingWindow.get_height().value() * _averageMarginInWindingWindow;
            }
        }
        else {
            auto radialHeight = windingWindow.get_radial_height().value();
            if (_averageMarginInWindingWindow > radialHeight / 2) {
                return {false, 0.0};
            }
            double marginAngle = wound_distance_to_angle(_averageMarginInWindingWindow, radialHeight / 2);
            if (std::isnan((marginAngle) / 360)) {
                throw NaNResultException("marginAngle: " + std::to_string(marginAngle));
            }
            // The margin blocks marginAngle degrees of the toroidal window; the
            // remaining (360 - marginAngle) degrees are usable.
            windingWindowArea *= (360 - marginAngle) / 360;
        }
    }
    double areaProductCore = windingWindowArea * windingColumn.get_area();
    double maximumAreaProductRequired = 0;


    // Phase 6 (perf): cache operating-points reference; get_operating_point(i)
    // returns by value (deep copy of OperatingPoint). Hot path: this loop is
    // executed once per candidate magnetic (>2000 cores) during AreaProduct
    // filtering.
    const auto& operatingPoints = inputs->get_operating_points();
    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        const auto& operatingPoint = operatingPoints[operatingPointIndex];
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();
        double frequency = Inputs::get_switching_frequency(Inputs::get_primary_excitation(operatingPoint));
        // double switchingFrequency = Inputs::get_switching_frequency(excitation);

        auto skinDepth = _windingSkinEffectLossesModel.calculate_skin_depth("copper", frequency, temperature);  // TODO material hardcoded
        double wireAirFillingFactor = Wire::get_filling_factor_round(2 * skinDepth);
        double windingWindowUtilizationFactor = wireAirFillingFactor * bobbinFillingFactor;
        double magneticFluxDensityPeakAtFrequencyOfReferenceLosses;
        if (core.get_material_name() == DUMMY_SENTINEL_NAME) {
            // Phase 1 fix: explicit branch for the shape-only pre-filter stage.
            // CoreAdviser uses the sentinel material name "Dummy" before material
            // fan-out (see CoreAdviser.cpp:2344, 2407): the actual material is
            // chosen later from the ferrite catalogue. During this stage, no real
            // loss model can be queried, so we sieve on shape using a reference
            // flux density. This is a *named-sentinel* contract, not a silent
            // catch-all fallback.
            magneticFluxDensityPeakAtFrequencyOfReferenceLosses = _magneticFluxDensityReference;
        }
        else {
            try {
                if (!_materialScaledMagneticFluxDensities.contains(core.get_material_name())) {
                    auto coreLossesMethods = core.get_available_core_losses_methods();

                    if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), VolumetricCoreLossesMethodType::STEINMETZ) != coreLossesMethods.end()) {
                        double referenceCoreLosses = _coreLossesModelSteinmetz->get_core_losses(core, _operatingPointExcitation, temperature).get_core_losses();
                        auto aux = _coreLossesModelSteinmetz->get_magnetic_flux_density_from_core_losses(core, frequency, temperature, referenceCoreLosses);
                        magneticFluxDensityPeakAtFrequencyOfReferenceLosses = aux.get_processed().value().get_peak().value();
                    }
                    else {
                        double referenceCoreLosses = _coreLossesModelProprietary->get_core_losses(core, _operatingPointExcitation, temperature).get_core_losses();

                        auto aux = _coreLossesModelProprietary->get_magnetic_flux_density_from_core_losses(core, frequency, temperature, referenceCoreLosses);
                        magneticFluxDensityPeakAtFrequencyOfReferenceLosses = aux.get_processed().value().get_peak().value();
                    }
                    _materialScaledMagneticFluxDensities[core.get_material_name()] = magneticFluxDensityPeakAtFrequencyOfReferenceLosses;
                }
                else {
                    magneticFluxDensityPeakAtFrequencyOfReferenceLosses = _materialScaledMagneticFluxDensities[core.get_material_name()];
                }
            }
            catch (...) {
                // Phase 1 fix: previously silently substituted _magneticFluxDensityReference
                // here, hiding loss-model failures for real materials and letting
                // unsuitable candidates pass area-product check with a fake B. Per the
                // no-silent-fallbacks policy, reject the candidate explicitly when its
                // real-material loss model cannot produce a peak B.
                return {false, 0.0};
            }
        }
        double areaProductRequired = _areaProductRequiredPreCalculations[operatingPointIndex] / (windingWindowUtilizationFactor * magneticFluxDensityPeakAtFrequencyOfReferenceLosses);
        if (std::isnan(magneticFluxDensityPeakAtFrequencyOfReferenceLosses) || magneticFluxDensityPeakAtFrequencyOfReferenceLosses == 0) {
            throw NaNResultException("magneticFluxDensityPeakAtFrequencyOfReferenceLosses cannot be 0 or NaN");
        }
        if (std::isnan(areaProductRequired)) {
            // Phase 1 fix: was a silent `break` that exited the operating-point loop
            // and let whatever maximumAreaProductRequired had accumulated decide
            // validity. NaN here means at least one operating point cannot be
            // characterised for this candidate → reject the candidate.
            return {false, 0.0};
        }
        if (std::isinf(areaProductRequired) || areaProductRequired == 0) {
            throw NaNResultException("areaProductRequired cannot be 0 or NaN");
        }

        maximumAreaProductRequired = std::max(maximumAreaProductRequired, areaProductRequired);
    }
    if (maximumAreaProductRequired > 1) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "maximumAreaProductRequired cannot be larger than 1 (probably)");
    }

    bool valid = areaProductCore >= maximumAreaProductRequired * defaults.coreAdviserThresholdValidity;
    double scoring = fabs(areaProductCore - maximumAreaProductRequired);

    return {valid, scoring};
}

double MagneticFilterAreaProduct::get_estimated_area_product_required(Inputs inputs) {
    double maxAp = 0;
    const auto& operatingPoints = inputs.get_operating_points();
    for (size_t i = 0; i < operatingPoints.size(); ++i) {
        const auto& operatingPoint = operatingPoints[i];
        double temperature = operatingPoint.get_conditions().get_ambient_temperature();
        double frequency = Inputs::get_switching_frequency(Inputs::get_primary_excitation(operatingPoint));
        auto skinDepth = _windingSkinEffectLossesModel.calculate_skin_depth("copper", frequency, temperature);
        double wireAirFillingFactor = Wire::get_filling_factor_round(2 * skinDepth);
        double bobbinFillingFactor = 0.45;
        double kFill = wireAirFillingFactor * bobbinFillingFactor;
        double ap = _areaProductRequiredPreCalculations[i] / (kFill * _magneticFluxDensityReference);
        maxAp = std::max(maxAp, ap);
    }
    return maxAp;
}

MagneticFilterAreaNoParallels::MagneticFilterAreaNoParallels(int maximumNumberParallels) {
    _maximumNumberParallels = maximumNumberParallels;
}

std::pair<bool, double> MagneticFilterAreaNoParallels::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;
    for (auto winding : magnetic->get_coil().get_functional_description()) {
        auto section = magnetic->get_mutable_coil().get_sections_by_winding(winding.get_name())[0];
        auto [auxValid, auxScoring] = evaluate_magnetic(winding, section);
        valid &= auxValid;
        scoring += auxScoring;
    }
    scoring /= magnetic->get_coil().get_functional_description().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterAreaNoParallels::evaluate_magnetic(Winding winding, Section section) {
    auto wire = Coil::resolve_wire(winding);

    if (wire.get_type() == WireType::FOIL && winding.get_number_parallels() * winding.get_number_turns() > _maximumNumberParallels) {
        return {false, 0.0};
    }

    if (!section.get_coordinate_system() || section.get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
        if (wire.get_maximum_outer_width() < section.get_dimensions()[0] && wire.get_maximum_outer_height() < section.get_dimensions()[1]) {
            return {true, 0.0};
        }
        else {
            return {false, 0.0};
        }
    }
    else {
        double wireAngle = wound_distance_to_angle(wire.get_maximum_outer_height(), wire.get_maximum_outer_width());

        if (wire.get_maximum_outer_width() < section.get_dimensions()[0] && wireAngle < section.get_dimensions()[1]) {
            return {true, 0.0};
        }
        else {
            return {false, 0.0};
        }
    }
}

std::pair<bool, double> MagneticFilterAreaWithParallels::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double scoring = 0;

    for (auto winding : magnetic->get_coil().get_functional_description()) {
        auto sections = magnetic->get_mutable_coil().get_sections_by_winding(winding.get_name());
        auto section = sections[0];
        double sectionArea;
        if (!section.get_coordinate_system() || section.get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
            sectionArea = section.get_dimensions()[0] * section.get_dimensions()[1];
        }
        else {
            sectionArea = std::numbers::pi * pow(section.get_dimensions()[0], 2) * section.get_dimensions()[1] / 360;
        }
        auto [auxValid, auxScoring] = evaluate_magnetic(winding, section, sections.size(), sectionArea, false);
        valid &= auxValid;
        scoring += auxScoring;
    }
    scoring /= magnetic->get_coil().get_functional_description().size();

    return {valid, scoring};
}

std::pair<bool, double> MagneticFilterAreaWithParallels::evaluate_magnetic(Winding winding, Section section, double numberSections, double sectionArea, bool allowNotFit) {
    auto wire = Coil::resolve_wire(winding);
    if (!Coil::resolve_wire(winding).get_conducting_area()) {
        throw CoilNotProcessedException("Conducting area is missing");
    }
    auto neededOuterAreaNoCompact = wire.get_maximum_outer_width() * wire.get_maximum_outer_height();

    neededOuterAreaNoCompact *= winding.get_number_parallels() * winding.get_number_turns() / numberSections;

    if (neededOuterAreaNoCompact < sectionArea) {
        // double scoring = (section.get_dimensions()[0] * section.get_dimensions()[1]) - neededOuterAreaNoCompact;
        return {true, 1.0};
    }
    else if (allowNotFit) {
        double extra = (neededOuterAreaNoCompact - sectionArea) / sectionArea;
        if (extra < 0.5) {
            // Higher scoring must mean better: a small overflow scores close to
            // the fitting case (1.0), a large overflow scores lower.
            return {true, 1.0 - extra};
        }
        else {
            return {false, 0.0};
        }
    }
    else {
        return {false, 0.0};
    }
}

} // namespace OpenMagnetics
