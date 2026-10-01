#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"
#include "constructive_models/Bobbin.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/Wire.h"
#include "physical_models/MagneticEnergy.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "support/Exceptions.h"
#include "support/Utils.h"
#include "Constants.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <string>

namespace OpenMagnetics {

MagneticFilterEnergyStored::MagneticFilterEnergyStored(Inputs inputs, std::map<std::string, std::string> models) {
    _models = models;
    _magneticEnergy = MagneticEnergy(models);
    _requiredMagneticEnergy = resolve_dimensional_values(_magneticEnergy.calculate_required_magnetic_energy(inputs));
}

MagneticFilterEnergyStored::MagneticFilterEnergyStored(Inputs inputs) {
    _requiredMagneticEnergy = resolve_dimensional_values(_magneticEnergy.calculate_required_magnetic_energy(inputs));
}

std::pair<bool, double> MagneticFilterEnergyStored::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    bool valid = true;
    double totalStorableMagneticEnergy = 0;
    // Phase 6 (perf): cache operating-points by const-ref.
    const auto& operatingPoints = inputs->get_operating_points();
    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        const auto& operatingPoint = operatingPoints[operatingPointIndex];
        double storableEnergy = _magneticEnergy.calculate_core_maximum_magnetic_energy(magnetic->get_core(), operatingPoint);
        totalStorableMagneticEnergy = std::max(totalStorableMagneticEnergy, storableEnergy);

        if (totalStorableMagneticEnergy >= _requiredMagneticEnergy * defaults.coreAdviserThresholdValidity) {
            if (outputs != nullptr) {
                while (outputs->size() < operatingPointIndex + 1) {
                    outputs->push_back(Outputs());
                }
                MagnetizingInductanceOutput magnetizingInductanceOutput;
                magnetizingInductanceOutput.set_maximum_magnetic_energy_core(storableEnergy);
                magnetizingInductanceOutput.set_method_used(_models["gapReluctance"]);
                magnetizingInductanceOutput.set_origin(ResultOrigin::SIMULATION);
                InductanceOutput inductanceOutput;
                if ((*outputs)[operatingPointIndex].get_inductance()) {
                    inductanceOutput = *(*outputs)[operatingPointIndex].get_inductance();
                }
                inductanceOutput.set_magnetizing_inductance(magnetizingInductanceOutput);
                (*outputs)[operatingPointIndex].set_inductance(inductanceOutput);
            }
        }
        else {
            valid = false;
            break;
        }
    }

    return {valid, totalStorableMagneticEnergy};
}

MagneticFilterEstimatedCost::MagneticFilterEstimatedCost(Inputs inputs) {
    double primaryCurrentRms = 0;
    double frequency = 0;
    double temperature = inputs.get_maximum_temperature();
    // Phase 6 (perf): cache operating-points by const-ref.
    const auto& operatingPoints = inputs.get_operating_points();
    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        const auto& operatingPoint = operatingPoints[operatingPointIndex];
        auto excitation = Inputs::get_primary_excitation(operatingPoint);
        primaryCurrentRms = std::max(primaryCurrentRms, excitation.get_current().value().get_processed().value().get_rms().value());
        frequency = std::max(frequency, Inputs::get_switching_frequency(excitation));
    }

    auto windingSkinEffectLossesModel = WindingSkinEffectLosses();
    _skinDepth = windingSkinEffectLossesModel.calculate_skin_depth("copper", frequency, temperature);  // TODO material hardcoded
    _wireAirFillingFactor = Wire::get_filling_factor_round(2 * _skinDepth);
    double estimatedWireConductingArea = std::numbers::pi * pow(_skinDepth, 2);
    _estimatedWireTotalArea = estimatedWireConductingArea / _wireAirFillingFactor;
    double necessaryWireCopperArea = primaryCurrentRms / defaults.maximumCurrentDensity;
    _estimatedParallels = ceil(necessaryWireCopperArea / estimatedWireConductingArea);

    if (settings.get_core_adviser_include_margin() && inputs.get_design_requirements().get_insulation()) {
        auto clearanceAndCreepageDistance = InsulationCoordinator().calculate_creepage_distance(inputs, true);
        _averageMarginInWindingWindow = clearanceAndCreepageDistance;
    }
}

std::pair<bool, double> MagneticFilterEstimatedCost::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    const auto& core = magnetic->get_core();

    double primaryNumberTurns = magnetic->get_coil().get_functional_description()[0].get_number_turns();
    double estimatedNeededWindingArea = primaryNumberTurns * _estimatedParallels * _estimatedWireTotalArea * (inputs->get_design_requirements().get_turns_ratios().size() + 1);
    WindingWindowElement windingWindow;

    std::string shapeName = core.get_shape_name();
    if (!is_pqi_or_ui_shape(shapeName)) {
        auto bobbin = Bobbin::create_quick_bobbin(core);
        windingWindow = bobbin.get_processed_description().value().get_winding_windows()[0];
    }
    else {
        windingWindow = core.get_winding_windows()[0];

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

    bool valid = windingWindowArea >= estimatedNeededWindingArea * defaults.coreAdviserThresholdValidity;


    double manufacturabilityRelativeCost;
    if (core.get_functional_description().get_type() != CoreType::TOROIDAL) {
        double estimatedNeededLayers = (primaryNumberTurns * _estimatedParallels * (2 * _skinDepth / _wireAirFillingFactor)) / windingWindow.get_height().value();
        manufacturabilityRelativeCost = estimatedNeededLayers;
    }
    else {
        if (_skinDepth >= windingWindow.get_radial_height().value()) {
            return {false, 0.0};
        }
        double layerLength = 2 * std::numbers::pi * (windingWindow.get_radial_height().value() - _skinDepth);
        double estimatedNeededLayers = (primaryNumberTurns * _estimatedParallels * (2 * _skinDepth / _wireAirFillingFactor)) / layerLength;
        if (estimatedNeededLayers < 0) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_INPUT, "estimatedNeededLayers cannot be negative");
        }
        if (estimatedNeededLayers > 1) {
            manufacturabilityRelativeCost = estimatedNeededLayers * 2;
        }
        else {
            manufacturabilityRelativeCost = estimatedNeededLayers;
        }
    }
    if (core.get_functional_description().get_number_stacks() > 1) {
        manufacturabilityRelativeCost *= 2;  // Because we need a custom bobbin
    }

    return {valid, manufacturabilityRelativeCost};
}

CostBasis MagneticFilterCost::default_basis() {
    CostBasis basis;
    // Core price laws, price = c * mass^k, fitted by least squares in log-log (docs/cost/fit_core_prices.py,
    // with its data). Prices are the tier a 1000-piece order pays, from breaks of at least 10 pieces, per set
    // or per toroid:
    //  - DigiKey list prices read through https://www.findchips.com/search/<part number> on 2026-09-30 (exact
    //    part-number rows; four parts cross-checked against DigiKey's own pages, every tier equal);
    //  - NiZn, carbonyl-iron and nanocrystalline toroids also from Farnell, TME and rf-microwave.com, the median
    //    over a part's distributors, EUR and GBP at the ECB reference rate of 2026-09-30;
    //  - Micrometals iron and proprietary-grade toroids from Power Magnetics (UK; highest break 50 pieces) and
    //    Feryster (PL; its breaks are a fixed discount on the 1-piece price, not negotiated volume prices), and
    //    Proterial AMCC amorphous C-core pairs from eu.mouser.com, 2026-10-01.
    // EU shops quote the DigiKey-priced MnZn parts at x1.06 (Gateway) to x1.20 (TME) of DigiKey. The price basis
    // of most TDK and Ferroxcube parts is not stated; Ferroxcube's datasheets that state it sell E and EQ per
    // half, the rest is inferred from the datasheet's mass unit. The mass is what get_mass() returns.
    // "proprietary" names no alloy, so its laws carry the manufacturer they were fitted on.
    // n = parts; spread = median / 90th-percentile factor between a part's price and the law.
    basis.corePrice["ferrite/MnZn"] = {21.03, 0.6407};                    // n=335, 0.08-2504 g, x1.61 / x2.83
    basis.corePrice["ferrite/NiZn"] = {30.61, 0.6973};                    // n=74, 0.05-360 g, x1.37 / x2.49 (toroids)
    basis.corePrice["powder/FeSiAl"] = {14.29, 0.5174};                   // n=63, 0.06-2238 g, x1.33 / x1.95 (Kool Mu)
    basis.corePrice["powder/FeNi"] = {65.05, 0.6985};                     // n=41, 1-678 g, x1.27 / x2.09 (High Flux, Edge)
    basis.corePrice["powder/FeNiMo"] = {105.7, 0.6287};                   // n=18, 0.16-276 g, x1.44 / x2.79 (MPP toroids)
    basis.corePrice["powder/FeSi"] = {23.76, 0.8595};                     // n=27, 6.7-1834 g, x1.17 / x1.43 (XFlux)
    basis.corePrice["powder/carbonylIron"] = {27.87, 0.6799};             // n=72, 0.26-3670 g, x1.24 / x1.97 (toroids)
    basis.corePrice["powder/iron"] = {12.73, 0.6543};                     // n=120, 0.34-5140 g, x1.52 / x2.46 (Micrometals toroids)
    // Weakest law: one shop (Feryster), exponent above 1, x3.34 at the 90th percentile.
    basis.corePrice["powder/proprietary/Micrometals"] = {125.5, 1.2489};  // n=24, 16-620 g, x1.61 / x3.34 (OC/OD/GX toroids)
    basis.corePrice["amorphous/FeSi"] = {47.64, 0.6708};                  // n=20, 150-3670 g, x1.10 / x1.26 (AMCC C-core pairs)
    basis.corePrice["nanocrystalline/proprietary/Vacuumschmelze"] = {60.17, 0.7750};  // n=18, 1.1-167 g, x1.15 / x1.30 (VITROPERM toroids)
    return basis;
}

CostBreakdown MagneticFilterCost::calculate_cost(Magnetic& magnetic) const {
    const CostBasis basis = _basis ? *_basis : default_basis();
    auto& coil = magnetic.get_mutable_coil();
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Missing turns description to evaluate cost filter");
    }
    // Conductor: every turn's length x its wire's conducting area x density, at the copper price, with
    // Sullivan's strand premium for litz.
    Constants constants;
    double conductor = 0;
    double turnsTotal = 0;
    auto wires = coil.get_wires();
    auto turns = coil.get_turns_description().value();
    for (const auto& turn : turns) {
        auto windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        auto wire = wires.at(windingIndex);
        auto materialName = Wire::resolve_material(wire).get_name();
        if (materialName != "copper") {
            throw InvalidInputException(ErrorCode::INVALID_WIRE_DATA, "the COST filter has a price for copper conductors only, not " + materialName);
        }
        double mass = turn.get_length() * wire.calculate_conducting_area() * constants.copperDensity;
        double premium = 1;
        if (wire.get_type() == WireType::LITZ) {
            double d = resolve_dimensional_values(wire.resolve_strand().get_conducting_diameter());
            premium = 1 + basis.litzK1 / std::pow(d, 6) + basis.litzK2 / std::pow(d, 2);
        }
        conductor += mass * basis.copperPricePerKg * premium;
        turnsTotal += 1;
    }
    // Labour: winding time from the turn count, grossed up by winding's share of the labour.
    double minutes = turnsTotal * basis.windingMinutesPerTurn / basis.windingShareOfLabour;
    double labour = minutes / 60 * basis.laborRatePerHour;
    // Core: the price law of its material type and composition at its mass.
    auto& core = magnetic.get_mutable_core();
    auto material = core.resolve_material();
    if (!material.get_material_composition()) {
        return {std::nullopt, conductor, labour, "core material " + material.get_name() + " has no material composition to price"};
    }
    json typeJson = material.get_material();
    json compositionJson = material.get_material_composition().value();
    std::string composition = typeJson.get<std::string>() + "/" + compositionJson.get<std::string>();
    if (compositionJson.get<std::string>() == "proprietary") {
        // A proprietary grade is priced only by a law fitted on its own manufacturer's parts.
        auto manufacturer = material.get_manufacturer_info().get_name();
        composition += "/" + manufacturer;
    }
    auto price = basis.corePrice.find(composition);
    if (price == basis.corePrice.end()) {
        return {std::nullopt, conductor, labour, "the COST filter has no price law for " + composition + " cores (" + material.get_name() + ")"};
    }
    double coreMass = core.get_mass();
    if (std::isnan(coreMass)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_MATERIAL_DATA, "core material " + material.get_name() + " has no density: its mass cannot be priced");
    }
    return {price->second.coefficient * std::pow(coreMass, price->second.exponent), conductor, labour, ""};
}

double CostBreakdown::total() const {
    if (!core) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_MATERIAL_DATA, "unit cost unknown: " + coreUnpricedReason);
    }
    return *core + conductor + labour;
}

std::pair<bool, double> MagneticFilterCost::evaluate_magnetic(Magnetic* magnetic, Inputs* inputs, std::vector<Outputs>* outputs) {
    auto cost = calculate_cost(*magnetic);
    if (!cost.priced()) {
        return {false, 0.0};
    }
    return {true, cost.total()};
}

} // namespace OpenMagnetics
