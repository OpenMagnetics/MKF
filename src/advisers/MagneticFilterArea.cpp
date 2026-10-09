#include "advisers/MagneticFilter.h"
#include "advisers/CoilAdviser.h"
#include "advisers/WireAdviser.h"
#include "advisers/MagneticFilterInternal.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Insulation.h"
#include "physical_models/CoreLosses.h"
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
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
    ProcessedWaveform processed;
    _operatingPointExcitation.set_frequency(frequencyReference);
    processed.set_label(WaveformLabel::SINUSOIDAL);
    processed.set_offset(0);
    processed.set_peak(_magneticFluxDensityReference);
    processed.set_peak_to_peak(2 * _magneticFluxDensityReference);
    magneticFluxDensity.set_processed(processed);
    _operatingPointExcitation.set_magnetic_flux_density(magneticFluxDensity);
    _coreLossesModelSteinmetz = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Steinmetz"}}));
    _coreLossesModelProprietary = CoreLossesModel::factory(std::map<std::string, std::string>({{"coreLosses", "Proprietary"}}));

    if (settings.get_core_adviser_include_margin() && inputs.has_insulation_coordination_requirements()) {
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
        require_voltage_and_current_for_power(excitation, operatingPointIndex,
                                              operatingPoints[operatingPointIndex],
                                              "sizing the core by area product");
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

double MagneticFilterAreaProduct::get_bobbin_filling_factor(const Core& core, std::optional<WiringTechnology> wiringTechnology) {
    if (core.get_winding_windows().size() == 0) {
        throw CoreNotProcessedException("Bobbin filling factor: core " + core.get_name().value_or("?") + " has no winding window");
    }
    auto windingWindow = core.get_winding_windows()[0];
    if (wiringTechnology && wiringTechnology.value() == WiringTechnology::PRINTED) {
        return 1;
    }
    if (core.get_functional_description().get_type() != CoreType::TOROIDAL) {
        return Bobbin::get_filling_factor(windingWindow.get_width().value(), windingWindow.get_height().value());
    }
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
        return 0.55 + 0.15 * circumferenceRatio;
    }
    return 0.6;  // Default for toroids without radial height
}

std::pair<bool, double> MagneticFilterAreaProduct::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    const auto& core = magnetic->get_core();

    if (core.get_winding_windows().size() == 0)
        return {false, 0.0};
    auto windingWindow = core.get_winding_windows()[0];
    auto windingColumn = core.get_columns()[0];
    double bobbinFillingFactor;
    if (!_bobbinFillingFactors.contains(core.get_shape_name())) {
        bobbinFillingFactor = get_bobbin_filling_factor(core, inputs->get_wiring_technology());
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
    return {wire_fits(wire, winding.get_number_parallels(), winding.get_number_turns(), section), 0.0};
}

bool MagneticFilterAreaNoParallels::wire_fits(Wire& wire, int64_t numberParallels, int64_t numberTurns, const Section& section) const {
    if (wire.get_type() == WireType::FOIL && numberParallels * numberTurns > _maximumNumberParallels) {
        return false;
    }

    if (!section.get_coordinate_system() || section.get_coordinate_system().value() == CoordinateSystem::CARTESIAN) {
        return wire.get_maximum_outer_width() < section.get_dimensions()[0] && wire.get_maximum_outer_height() < section.get_dimensions()[1];
    }
    else {
        double wireAngle = wound_distance_to_angle(wire.get_maximum_outer_height(), wire.get_maximum_outer_width());
        return wire.get_maximum_outer_width() < section.get_dimensions()[0] && wireAngle < section.get_dimensions()[1];
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

// ABT #1446: whether the winding's physical turns can be laid out in whole turns per layer
// inside a rectangular (Cartesian) section, using the same packing rules the layer winder
// (Coil::wind_by_rectangular_layers) applies: turns per layer = floor(section length along the
// layer / conductor size along the layer), layers = ceil(turns / turns per layer), and the
// layers' stacked thickness must fit the section's depth. The bulk outer-area comparison alone
// passed a 0.85 mm litz for 13 turns in a 2.27 x 4.59 mm section (area ratio 0.9), though
// 5 turns per layer x 3 layers needs 2.55 mm of depth: the winder then rejected every wire the
// adviser had offered and the coil adviser answered "No coil found".
static bool fits_in_whole_turns_per_layer(const Wire& wire, const Winding& winding, const Section& section, double numberSections) {
    if (section.get_dimensions().size() < 2) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
            "Section " + section.get_name() + " has no width and height to lay turns out in");
    }
    if (!(numberSections > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
            "Winding " + winding.get_name() + " is split into a non-positive number of sections");
    }
    double sectionWidth = section.get_dimensions()[0];
    double sectionHeight = section.get_dimensions()[1];
    bool layersStackAlongWidth = section.get_layers_orientation() == WindingOrientation::OVERLAPPING;

    double conductorAlongLayer;
    double layerThickness;
    std::optional<double> fixedTurnsPerLayer;
    if (wire.get_type() == WireType::ROUND || wire.get_type() == WireType::LITZ) {
        if (!wire.get_outer_diameter()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Wire " + wire.get_name().value_or("(unnamed)") + " is missing its outer diameter");
        }
        double diameter = resolve_dimensional_values(wire.get_outer_diameter().value());
        conductorAlongLayer = diameter;
        layerThickness = diameter;
    }
    else {
        if (!wire.get_outer_width() || !wire.get_outer_height()) {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Wire " + wire.get_name().value_or("(unnamed)") + " is missing its outer width or height");
        }
        double wireWidth = resolve_dimensional_values(wire.get_outer_width().value());
        double wireHeight = resolve_dimensional_values(wire.get_outer_height().value());
        if (layersStackAlongWidth) {
            conductorAlongLayer = wireHeight;
            layerThickness = wireWidth;
            if (wire.get_type() == WireType::FOIL) {
                fixedTurnsPerLayer = 1;
            }
        }
        else {
            conductorAlongLayer = wireWidth;
            layerThickness = wireHeight;
            if (wire.get_type() == WireType::RECTANGULAR && settings.get_coil_only_one_turn_per_layer_in_contiguous_rectangular()) {
                fixedTurnsPerLayer = 1;
            }
        }
    }
    if (!(conductorAlongLayer > 0) || !(layerThickness > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "Wire " + wire.get_name().value_or("(unnamed)") + " has a non-positive outer dimension");
    }

    double lengthAlongLayer = layersStackAlongWidth ? sectionHeight : sectionWidth;
    double depth = layersStackAlongWidth ? sectionWidth : sectionHeight;
    double turnsPerLayer = fixedTurnsPerLayer ? fixedTurnsPerLayer.value() : std::floor(lengthAlongLayer / conductorAlongLayer);
    if (turnsPerLayer < 1) {
        return false;
    }
    double physicalTurns = std::ceil(static_cast<double>(winding.get_number_turns()) * static_cast<double>(winding.get_number_parallels()) / numberSections);

    // ABT #1699: with real winding geometry on, every layer a parallel enters also carries that
    // parallel's crossing station (ABT #685: a wire making N turns over L layers crosses the
    // window plane N + L times; Coil::wind_inner's RealWindingCrossingBump charges one extra
    // cross-section per layer and parallel). A layer of T positions then holds T - P turns of a
    // P-parallel winding, so N P / S turns need L layers with L (T - P) >= N P / S, which is the
    // fixpoint the bump iterates to ((N/S + L) P <= L T). With T <= P no number of layers closes
    // (each new layer brings as many stations as it has room): the real winder refuses it. The
    // bump charges no station to a foil (ABT #881) or a one-turn winding (the omega), and none in
    // a toroid's round window, whose sections are polar and never reach this function. With the
    // setting off nothing changes: the ideal wind charges no station.
    double stationsPerLayer = 0;
    if (settings.get_coil_use_real_winding_geometry() && wire.get_type() != WireType::FOIL && winding.get_number_turns() > 1) {
        stationsPerLayer = static_cast<double>(winding.get_number_parallels());
    }
    double turnsPerLayerAfterStations = turnsPerLayer - stationsPerLayer;
    if (turnsPerLayerAfterStations < 1) {
        return false;
    }
    double numberLayers = std::ceil(physicalTurns / turnsPerLayerAfterStations);
    return numberLayers * layerThickness <= depth * (1 + 1e-9);
}

std::pair<bool, double> MagneticFilterAreaWithParallels::evaluate_magnetic(Winding winding, Section section, double numberSections, double sectionArea, bool allowNotFit) {
    auto wire = Coil::resolve_wire(winding);
    if (!Coil::resolve_wire(winding).get_conducting_area()) {
        throw CoilNotProcessedException("Conducting area is missing");
    }
    auto neededOuterAreaNoCompact = wire.get_maximum_outer_width() * wire.get_maximum_outer_height();

    neededOuterAreaNoCompact *= winding.get_number_parallels() * winding.get_number_turns() / numberSections;

    // ABT #1446: a strict fit also needs the turns to pack in whole turns per layer. Only
    // rectangular sections are checked: in a toroid's polar section the turns a layer holds
    // depend on that layer's radius (the inner circumference shrinks layer by layer), which the
    // toroidal winder works out turn by turn; a per-layer count here would be a second,
    // diverging model of it. The allowNotFit pass (taken only when no wire fits strictly) keeps
    // its existing bulk-area tolerance. Planar wires are not checked either: a planar winding
    // is laid out per PCB layer by the planar winder (Coil::wind_by_planar_layers), not by the
    // layer winder whose rules fits_in_whole_turns_per_layer mirrors.
    bool packs = true;
    if (!allowNotFit && wire.get_type() != WireType::PLANAR &&
        (!section.get_coordinate_system() || section.get_coordinate_system().value() == CoordinateSystem::CARTESIAN)) {
        packs = fits_in_whole_turns_per_layer(wire, winding, section, numberSections);
    }

    if (neededOuterAreaNoCompact < sectionArea && packs) {
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

std::pair<bool, double> MagneticFilterRealWinding::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, [[maybe_unused]] std::vector<Outputs>* outputs) {
    _lastReason.clear();
    if (!magnetic->has_core() || !magnetic->has_coil()) {
        throw CoilNotProcessedException("Real winding check: the magnetic has no core and coil to rebuild");
    }
    if (magnetic->get_mutable_coil().is_planar()) {
        throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
            "Real winding check: a planar coil has no real-winding model (ABT #492)");
    }
    std::optional<Inputs> realWindingInputs;
    if (inputs) {
        realWindingInputs = *inputs;
    }
    auto refusal = real_winding_refusal(*magnetic, realWindingInputs);
    if (refusal) {
        _lastReason = refusal.value();
        return {false, 0.0};
    }
    return {true, 1.0};
}

bool MagneticFilterWindowCopperCapacity::applies_to(const Inputs& inputs) {
    return inputs.get_wiring_technology() != WiringTechnology::PRINTED;
}

MagneticFilterWindowCopperCapacity::MagneticFilterWindowCopperCapacity(Inputs inputs, double maximumEffectiveCurrentDensity)
    : _maximumEffectiveCurrentDensity(maximumEffectiveCurrentDensity)
{
    if (inputs.get_operating_points().empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Window copper capacity: the inputs have no operating points");
    }
    if (!(maximumEffectiveCurrentDensity > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Window copper capacity: the maximum effective current density must be positive");
    }
    const size_t numberWindings = inputs.get_operating_points()[0].get_excitations_per_winding().size();
    const auto& turnsRatios = inputs.get_design_requirements().get_turns_ratios();
    if (turnsRatios.size() + 1 != numberWindings) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "Window copper capacity: the input has " + std::to_string(numberWindings) + " windings but " +
            std::to_string(turnsRatios.size()) + " turns ratios; the turns of the windings the stand-in coil lacks cannot be derived");
    }
    for (size_t ratioIndex = 0; ratioIndex < turnsRatios.size(); ++ratioIndex) {
        double turnsRatio = resolve_dimensional_values(turnsRatios[ratioIndex], DimensionalValues::NOMINAL);
        if (!(turnsRatio > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "Window copper capacity: turns ratio " + std::to_string(ratioIndex + 1) + " is not positive");
        }
        _turnsRatios.push_back(turnsRatio);
    }
    _temperature = inputs.get_maximum_temperature();
}

std::pair<bool, double> MagneticFilterWindowCopperCapacity::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, [[maybe_unused]] std::vector<Outputs>* outputs) {
    _lastReason.clear();
    const auto& core = magnetic->get_core();
    const std::string coreName = core.get_name().value_or(core.get_shape_name());
    if (core.get_winding_windows().empty() || !core.get_winding_windows()[0].get_area()) {
        throw CoreNotProcessedException("Window copper capacity: core " + coreName + " has no winding window area");
    }
    auto& coil = magnetic->get_mutable_coil();
    const auto& windings = coil.get_functional_description();
    const size_t numberWindings = _turnsRatios.size() + 1;
    if (windings.empty() || windings.size() > numberWindings) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "Window copper capacity: candidate " + coreName + " carries " + std::to_string(windings.size()) +
            " windings for an input with " + std::to_string(numberWindings));
    }
    const double primaryNumberTurns = static_cast<double>(windings[0].get_number_turns());
    if (!(primaryNumberTurns >= 1)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Window copper capacity: candidate " + coreName + " has no seeded turns");
    }

    // The strand: the stand-in's wire, two skin depths at the highest switching frequency.
    Wire strand = coil.resolve_wire(0);
    if (strand.get_type() != WireType::ROUND || !strand.get_conducting_diameter()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "Window copper capacity: the stand-in wire of " + coreName + " is not a round strand with a conducting diameter");
    }
    const double strandConductingArea = strand.calculate_conducting_area();
    const double strandConductingDiameter = resolve_dimensional_values(strand.get_conducting_diameter().value());
    const double bobbinFillingFactor = MagneticFilterAreaProduct::get_bobbin_filling_factor(core, inputs->get_wiring_technology());

    // Per turn, a winding needs strandsPerTurn stand-in strands: the effective current density
    // the current would have in ONE stand-in strand over the maximum. Above one, whole strands
    // (rounded up). At or below one, a single stand-in strand is more copper than the current
    // needs, and the coil stage is free to pick a thinner round wire: the copper is that of the
    // round conductor whose own effective current density is the maximum. At line frequency the stand-in strand is two skin
    // depths at 50 Hz (18.6 mm); counting it whole made a 32 A PFC choke need 272 mm2 of copper
    // per turn instead of 2.7 mm2, and rejected every core of the catalogue.
    // Each winding's copper occupies its area over the round-wire filling factor of its own
    // conductor diameter; the window offers its area times the bobbin filling factor.
    double requiredCopperArea = 0;
    double requiredWindowArea = 0;
    std::string perWinding;
    const bool standIn = std::holds_alternative<std::string>(coil.get_bobbin());
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        double numberTurns;
        // A stand-in (bobbin still a name) seeds the first winding's turns only: the common-mode
        // choke stand-in carries every winding, the others at the placeholder 1 turn.
        if (windingIndex < windings.size() && (windingIndex == 0 || !standIn)) {
            numberTurns = static_cast<double>(windings[windingIndex].get_number_turns());
        }
        else {
            numberTurns = static_cast<double>(std::max<int64_t>(1, std::llround(primaryNumberTurns / _turnsRatios[windingIndex - 1])));
        }
        double strandsPerTurn = 0;
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs->get_operating_points().size(); ++operatingPointIndex) {
            const auto& excitations = inputs->get_operating_points()[operatingPointIndex].get_excitations_per_winding();
            if (windingIndex >= excitations.size() || !excitations[windingIndex].get_current()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA,
                    "Window copper capacity: operating point " + std::to_string(operatingPointIndex) + " has no current for winding " +
                    std::to_string(windingIndex));
            }
            const auto current = excitations[windingIndex].get_current().value();
            double strandLoad = strand.calculate_effective_current_density(current, _temperature) / _maximumEffectiveCurrentDensity;
            if (strandLoad > 1) {
                strandsPerTurn = std::max(strandsPerTurn, std::ceil(strandLoad));
                continue;
            }
            // One strand is more than enough: the copper is that of the round conductor whose own
            // effective current density is the maximum, not the fraction strandLoad of the strand.
            // A thinner conductor carries the current's high-frequency harmonics (a PFC choke's
            // ripple at the switching frequency) no better per unit of copper than the strand, so
            // the fraction overcounts. The search depends on the inputs only: memoised.
            const auto cacheKey = std::make_tuple(strandConductingDiameter, windingIndex, operatingPointIndex);
            auto cached = _conductorAreaCache.find(cacheKey);
            if (cached == _conductorAreaCache.end()) {
                const double conductingArea = get_conducting_area_for_current(current, strandConductingArea, _temperature, _maximumEffectiveCurrentDensity);
                cached = _conductorAreaCache.emplace(cacheKey, conductingArea).first;
            }
            strandsPerTurn = std::max(strandsPerTurn, cached->second / strandConductingArea);
        }
        if (!(strandsPerTurn > 0)) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                "Window copper capacity: winding " + std::to_string(windingIndex) + " of " + coreName + " carries no current to size its copper");
        }
        const double conductorDiameter = strandsPerTurn < 1 ? strandConductingDiameter * std::sqrt(strandsPerTurn) : strandConductingDiameter;
        double windingCopperArea = numberTurns * strandsPerTurn * strandConductingArea;
        requiredCopperArea += windingCopperArea;
        requiredWindowArea += windingCopperArea / Wire::get_filling_factor_round(conductorDiameter);
        perWinding += (windingIndex == 0 ? "" : ", ") + std::to_string(std::llround(numberTurns)) + " turns x " +
                      std::to_string(strandsPerTurn) + " strands";
    }

    const double windingWindowArea = core.get_winding_windows()[0].get_area().value();
    const double usableWindowArea = windingWindowArea * bobbinFillingFactor;
    if (!(usableWindowArea > 0)) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
            "Window copper capacity: core " + coreName + " has no usable winding window");
    }
    const double proportion = requiredWindowArea / usableWindowArea;
    bool valid = proportion <= 1;
    if (!valid) {
        _lastReason = "core " + coreName + ": the windings (" + perWinding + " of " +
                      std::to_string(strandConductingDiameter * 1e3) + " mm strand at <= " +
                      std::to_string(_maximumEffectiveCurrentDensity * 1e-6) + " A/mm2) need " +
                      std::to_string(requiredCopperArea * 1e6) + " mm2 of copper, " +
                      std::to_string(requiredWindowArea * 1e6) + " mm2 of window at their round-wire fill; it offers " +
                      std::to_string(usableWindowArea * 1e6) + " mm2 (" + std::to_string(windingWindowArea * 1e6) +
                      " mm2 x bobbin filling factor " + std::to_string(bobbinFillingFactor) + ")";
    }
    return {valid, proportion};
}

MagneticFilterWireWithinLimits::MagneticFilterWireWithinLimits(double maximumEffectiveCurrentDensity, int maximumNumberParallels)
    : _maximumEffectiveCurrentDensity(maximumEffectiveCurrentDensity), _maximumNumberParallels(maximumNumberParallels)
{
    if (!(maximumEffectiveCurrentDensity > 0) || maximumNumberParallels < 1) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Wire within limits: the current density and parallels limits must be positive");
    }
    _coilAdviser = std::make_shared<CoilAdviser>();
    _wireAdviser = std::make_shared<WireAdviser>();
    _wireAdviser->set_maximum_effective_current_density(_maximumEffectiveCurrentDensity);
    _wireAdviser->set_maximum_number_parallels(_maximumNumberParallels);
    _wireAdviser->set_synthesize_litz(true);
}

std::pair<bool, double> MagneticFilterWireWithinLimits::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, [[maybe_unused]] std::vector<Outputs>* outputs) {
    _lastReason.clear();
    const auto& core = magnetic->get_core();
    const std::string coreName = core.get_name().value_or(core.get_shape_name());
    if (inputs->get_operating_points().empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Wire within limits: the inputs have no operating points");
    }
    const auto& windings = magnetic->get_coil().get_functional_description();
    const size_t numberWindings = inputs->get_operating_points()[0].get_excitations_per_winding().size();
    if (windings.size() != numberWindings) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "Wire within limits: candidate " + coreName + " carries " + std::to_string(windings.size()) +
            " windings for an input with " + std::to_string(numberWindings) + "; complete its coil first");
    }

    // The current each winding's wire is chosen for: the operating point with the highest
    // rms x sqrt(effective frequency), as the coil stage picks it.
    std::vector<SignalDescriptor> windingCurrents;
    for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
        double maximumFigure = -1;
        SignalDescriptor chosen;
        for (const auto& operatingPoint : inputs->get_operating_points()) {
            const auto& current = operatingPoint.get_excitations_per_winding()[windingIndex].get_current();
            if (!current || !current->get_processed() || !current->get_processed()->get_rms() ||
                !current->get_processed()->get_effective_frequency()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA,
                    "Wire within limits: winding " + std::to_string(windingIndex) + " has no processed current (rms, effective frequency)");
            }
            double figure = current->get_processed()->get_rms().value() * sqrt(current->get_processed()->get_effective_frequency().value());
            if (figure > maximumFigure) {
                maximumFigure = figure;
                chosen = current.value();
            }
        }
        windingCurrents.push_back(chosen);
    }
    double temperature = -std::numeric_limits<double>::max();
    for (const auto& operatingPoint : inputs->get_operating_points()) {
        temperature = std::max(temperature, operatingPoint.get_conditions().get_ambient_temperature());
    }

    auto coreType = core.get_functional_description().get_type();
    auto patterns = Coil::get_patterns(*inputs, coreType);
    auto repetitionsOptions = Coil::get_repetitions(*inputs, coreType);
    if (patterns.empty() || repetitionsOptions.empty()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Wire within limits: no winding pattern for candidate " + coreName);
    }
    // Every pattern and repetition the coil stage tries: the core passes if one of them has a
    // wire for every winding.
    std::vector<std::pair<std::vector<size_t>, size_t>> patternsAndRepetitions;
    for (const auto& repetitionsOption : repetitionsOptions) {
        for (const auto& pattern : patterns) {
            patternsAndRepetitions.push_back({pattern, repetitionsOption});
        }
    }

    Mas mas;
    mas.set_inputs(*inputs);
    mas.set_magnetic(*magnetic);

    std::string reasons;
    for (auto [pattern, repetitionsOption] : patternsAndRepetitions) {
        auto [checkedPattern, repetitions] = mas.get_mutable_magnetic().get_mutable_coil().check_pattern_and_repetitions_integrity(pattern, repetitionsOption);
        // As the coil stage does: each combination of solid insulation per wire is a separate
        // try, and a combination that needs margin tape narrows the sections by it.
        std::vector<std::optional<std::vector<WireSolidInsulationRequirements>>> insulationCombinations;
        if (mas.get_mutable_inputs().get_wiring_technology() == WiringTechnology::WOUND) {
            for (auto& combination : InsulationCoordinator::get_solid_insulation_requirements_for_wires(mas.get_mutable_inputs(), checkedPattern, repetitions)) {
                insulationCombinations.push_back(combination);
            }
            if (insulationCombinations.empty()) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Wire within limits: no solid insulation combination for candidate " + coreName);
            }
        }
        else {
            insulationCombinations.push_back(std::nullopt);
        }
        for (auto insulationCombination : insulationCombinations) {
            bool needsMargin = false;
            if (insulationCombination) {
                needsMargin = InsulationCoordinator::needs_margin(insulationCombination.value(), checkedPattern, repetitions,
                                                                  InsulationCoordinator::insulation_class_for_margin(mas.get_mutable_inputs()));
                if (needsMargin) {
                    CoilAdviser::limit_wire_insulation_requirements_for_margin(insulationCombination.value());
                }
            }
            auto sections = _coilAdviser->get_advised_sections(mas, checkedPattern, repetitions, needsMargin);
            if (sections.empty()) {
                reasons += (reasons.empty() ? "" : "; ") + std::string("no section layout");
                continue;
            }
            std::set<size_t> checkedWindings;
            std::string failure;
            for (const auto& section : sections) {
                if (section.get_type() != ElectricalType::CONDUCTION) {
                    continue;
                }
                for (const auto& partialWinding : section.get_partial_windings()) {
                    size_t windingIndex = numberWindings;
                    for (size_t index = 0; index < numberWindings; ++index) {
                        if (windings[index].get_name() == partialWinding.get_winding()) {
                            windingIndex = index;
                            break;
                        }
                    }
                    if (windingIndex == numberWindings) {
                        throw InvalidInputException(ErrorCode::INVALID_INPUT,
                            "Wire within limits: section winding '" + partialWinding.get_winding() + "' is not a winding of " + coreName);
                    }
                    if (!checkedWindings.insert(windingIndex).second) {
                        continue;
                    }
                    // The coil stage's verdict "within the limits": some wire of the catalogue (or
                    // litz synthesised for this current) with the insulation this combination asks
                    // of it fits the section with its parallels.
                    // The wires the coil stage advises from (preferred standard, insulated kept).
                    if (!_catalogueWires) {
                        _catalogueWires = _coilAdviser->get_catalogue_wires();
                    }
                    auto dataset = _wireAdviser->create_dataset(windings[windingIndex], &_catalogueWires.value(), section,
                                                                windingCurrents[windingIndex], temperature);
                    auto fittingAlone = _wireAdviser->filter_by_area_no_parallels(&dataset, section);
                    if (insulationCombination) {
                        if (windingIndex >= insulationCombination->size()) {
                            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                                "Wire within limits: the solid insulation combination has no entry for winding " + std::to_string(windingIndex));
                        }
                        fittingAlone = _wireAdviser->filter_by_solid_insulation_requirements(&fittingAlone, insulationCombination.value()[windingIndex]);
                    }
                    auto fittingWithParallels = _wireAdviser->filter_by_area_with_parallels(&fittingAlone, section, static_cast<double>(repetitions), false);
                    if (fittingWithParallels.empty()) {
                        failure = "winding '" + windings[windingIndex].get_name() + "' (" + std::to_string(windings[windingIndex].get_number_turns()) +
                                  " turns, " + std::to_string(windingCurrents[windingIndex].get_processed()->get_rms().value()) +
                                  " A rms) has no wire within " + std::to_string(_maximumEffectiveCurrentDensity * 1e-6) + " A/mm2 and " +
                                  std::to_string(_maximumNumberParallels) + " parallels" + (needsMargin ? " (margin taped)" : "") +
                                  " that fits its " + std::to_string(section.get_dimensions()[0] * 1e3) + " x " +
                                  std::to_string(section.get_dimensions()[1] * 1e3) + " mm section";
                        break;
                    }
                }
                if (!failure.empty()) {
                    break;
                }
            }
            if (failure.empty()) {
                return {true, 0};
            }
            reasons += (reasons.empty() ? "" : "; ") + failure;
        }
    }
    _lastReason = "core " + coreName + ": " + reasons;
    return {false, 1};
}

} // namespace OpenMagnetics
