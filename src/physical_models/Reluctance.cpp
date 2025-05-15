#include "physical_models/Reluctance.h"
#include "physical_models/MagneticEnergy.h"
#include "Defaults.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>

namespace OpenMagnetics {


double ReluctanceModel::get_ungapped_core_reluctance(const CoreWrapper& core, double initialPermeability) {
    auto constants = Constants();
    double absolutePermeability = constants.vacuumPermeability * initialPermeability;
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double effectiveLength = core.get_processed_description()->get_effective_parameters().get_effective_length();

    double reluctanceCore = effectiveLength / (absolutePermeability * effectiveArea);
    return reluctanceCore;
}

double ReluctanceModel::get_air_cored_reluctance(BobbinWrapper bobbin) {
    if (!bobbin.get_processed_description()) {
        throw std::runtime_error("Bobbin not processed");
    }
    double airArea;
    auto bobbinDimensions = bobbin.get_winding_window_dimensions();
    auto wallThickness = bobbin.get_processed_description()->get_wall_thickness();
    double airLength = bobbinDimensions[1] + 2 * wallThickness;
    auto columnShape = bobbin.get_processed_description()->get_column_shape();
    auto columnThickness = bobbin.get_processed_description()->get_column_thickness();
    auto columnWidth = bobbin.get_processed_description()->get_column_width().value();
    auto columnDepth = bobbin.get_processed_description()->get_column_depth();
    if (columnShape == ColumnShape::ROUND) {
        airArea = std::numbers::pi * pow(columnWidth - columnThickness, 2);
    }
    else if (columnShape == ColumnShape::RECTANGULAR) {
        airArea = 4 * (columnWidth - columnThickness) * (columnDepth - columnThickness);
    }
    else {
        throw std::runtime_error("Column shape not implemented yet");
    }

    auto constants = Constants();
    double absolutePermeability = constants.vacuumPermeability;

    double reluctanceAirCore = airLength / (absolutePermeability * airArea);
    return reluctanceAirCore;
}

double ReluctanceModel::get_ungapped_core_reluctance(CoreWrapper core, OperatingPoint* operatingPoint) {
    OpenMagnetics::InitialPermeability initialPermeability;

    auto coreMaterial = core.resolve_material();

    double initialPermeabilityValue;
    if (operatingPoint != nullptr) {
        double temperature = operatingPoint->get_conditions().get_ambient_temperature(); // TODO: Use a future calculated temperature
        _magneticFluxDensitySaturation = core.get_magnetic_flux_density_saturation(temperature, true);
        auto frequency = operatingPoint->get_excitations_per_winding()[0].get_frequency();
        initialPermeabilityValue = initialPermeability.get_initial_permeability(coreMaterial, temperature, std::nullopt, frequency);
    }
    else {
        initialPermeabilityValue = initialPermeability.get_initial_permeability(coreMaterial);
        _magneticFluxDensitySaturation = core.get_magnetic_flux_density_saturation(true);
    }
    return get_ungapped_core_reluctance(core, initialPermeabilityValue);
}


MagnetizingInductanceOutput ReluctanceModel::get_core_reluctance(CoreWrapper core, OperatingPoint* operatingPoint) {
    auto ungappedCoreReluctance = get_ungapped_core_reluctance(core, operatingPoint);
    auto magnetizingInductanceOutput = get_gapping_reluctance(core);

    if (std::isnan(ungappedCoreReluctance)) {
        throw std::runtime_error("Core Reluctance must be a number, not NaN");
    }
    double calculatedReluctance = ungappedCoreReluctance + magnetizingInductanceOutput.get_gapping_reluctance().value();
    if (std::isnan(calculatedReluctance)) {
        throw std::runtime_error("Reluctance must be a number, not NaN");
    }

    if (operatingPoint != nullptr) {
        magnetizingInductanceOutput.set_maximum_magnetic_energy_core(MagneticEnergy::get_ungapped_core_maximum_magnetic_energy(core, operatingPoint));
    }
    magnetizingInductanceOutput.set_core_reluctance(calculatedReluctance);
    magnetizingInductanceOutput.set_ungapped_core_reluctance(ungappedCoreReluctance);

    return magnetizingInductanceOutput;
}

MagnetizingInductanceOutput ReluctanceModel::get_core_reluctance(CoreWrapper core, double initialPermeability) {
    auto ungappedCoreReluctance = get_ungapped_core_reluctance(core, initialPermeability);

    auto magnetizingInductanceOutput = get_gapping_reluctance(core);
    double calculatedReluctance = ungappedCoreReluctance + magnetizingInductanceOutput.get_gapping_reluctance().value();

    magnetizingInductanceOutput.set_core_reluctance(calculatedReluctance);
    magnetizingInductanceOutput.set_ungapped_core_reluctance(ungappedCoreReluctance);

    return magnetizingInductanceOutput;
}

MagnetizingInductanceOutput ReluctanceModel::get_gapping_reluctance(CoreWrapper core) {
    double calculatedReluctance = 0;
    double calculatedCentralReluctance = 0;
    double calculatedLateralReluctance = 0;
    double maximumFringingFactor = 1;
    double maximumStorableMagneticEnergyGapping = 0;
    std::vector<AirGapReluctanceOutput> reluctancePerGap;
    auto gapping = core.get_functional_description().get_gapping();
    if (gapping.size() != 0) {
        // We recompute all gaps in case some is missing coordinates
        for (const auto& gap : gapping) {
            if (!gap.get_coordinates()) {
                core.process_gap();
                gapping = core.get_functional_description().get_gapping();
                break;
            }
        }


        for (const auto& gap : gapping) {
            auto gapReluctance = get_gap_reluctance(gap);
            auto gapColumn = core.find_closest_column_by_coordinates(gap.get_coordinates().value());
            reluctancePerGap.push_back(gapReluctance);
            if (gapColumn.get_type() == OpenMagnetics::ColumnType::LATERAL) {
                calculatedLateralReluctance += 1 / gapReluctance.get_reluctance();
            }
            else {
                calculatedCentralReluctance += gapReluctance.get_reluctance();
            }
            maximumFringingFactor = std::max(maximumFringingFactor, gapReluctance.get_fringing_factor());
            maximumStorableMagneticEnergyGapping += gapReluctance.get_maximum_storable_magnetic_energy();
            if (gapReluctance.get_fringing_factor() < 1) {
                std::cout << "fringing_factor " << gapReluctance.get_fringing_factor() << std::endl;
            }
        }
        calculatedReluctance = calculatedCentralReluctance + 1 / calculatedLateralReluctance;
    }

    MagnetizingInductanceOutput magnetizingInductanceOutput;

    magnetizingInductanceOutput.set_maximum_fringing_factor(maximumFringingFactor);
    magnetizingInductanceOutput.set_maximum_storable_magnetic_energy_gapping(maximumStorableMagneticEnergyGapping);

    magnetizingInductanceOutput.set_gapping_reluctance(calculatedReluctance);
    magnetizingInductanceOutput.set_reluctance_per_gap(reluctancePerGap);
    magnetizingInductanceOutput.set_method_used(methodName);
    magnetizingInductanceOutput.set_origin(ResultOrigin::SIMULATION);

    return magnetizingInductanceOutput;
}

AirGapReluctanceOutput ReluctanceZhangModel::get_gap_reluctance(CoreGap gapInfo) {
    double perimeter = 0;
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    if (!gapInfo.get_shape()) {
        throw std::runtime_error("Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw std::runtime_error("Gap Section Dimensions are not set");
    }
    if (!gapInfo.get_distance_closest_normal_surface()) {
        throw std::runtime_error("Gap Distance Closest Normal Surface is not set");
    }
    auto gapArea = *(gapInfo.get_area());
    auto gapShape = *(gapInfo.get_shape());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    auto distanceClosestNormalSurface = *(gapInfo.get_distance_closest_normal_surface());
    auto reluctanceInternal = gapLength / (constants.vacuumPermeability * gapArea);
    double reluctanceFringing = 0;
    double fringingFactor = 1;
    auto gapSectionWidth = gapSectionDimensions[0];
    auto gapSectionDepth = gapSectionDimensions[1];

    if (gapShape == ColumnShape::ROUND) {
        perimeter = std::numbers::pi * gapSectionWidth;
    }
    else { // TODO: Properly calculate perimeter for all shapes
        perimeter = gapSectionWidth * 2 + gapSectionDepth * 2;
    }

    double reluctance = 0;
    if (gapLength > 0) {
        reluctanceFringing = std::numbers::pi / (constants.vacuumPermeability * perimeter *
                                                  log((2 * distanceClosestNormalSurface + gapLength) / gapLength));

        if (std::isnan(reluctanceInternal) || reluctanceInternal == 0) {
            throw std::runtime_error("reluctanceInternal cannot be 0 or NaN");
        }

        if (std::isnan(reluctanceFringing) || reluctanceFringing == 0) {
            throw std::runtime_error("reluctanceFringing cannot be 0 or NaN");
        }

        reluctance = 1. / (1. / reluctanceInternal + 1. / reluctanceFringing);

        fringingFactor = gapLength / (constants.vacuumPermeability * gapArea * reluctance);
    }
    AirGapReluctanceOutput airGapReluctanceOutput;
    airGapReluctanceOutput.set_maximum_storable_magnetic_energy(get_gap_maximum_storable_energy(gapInfo, fringingFactor));
    airGapReluctanceOutput.set_reluctance(reluctance);
    airGapReluctanceOutput.set_method_used("Zhang");
    airGapReluctanceOutput.set_origin(ResultOrigin::SIMULATION);
    airGapReluctanceOutput.set_fringing_factor(fringingFactor);

    return airGapReluctanceOutput;
};


double ReluctanceMuehlethalerModel::get_basic_reluctance(double l, double w, double h) {
    auto constants = Constants();
    return 1 / constants.vacuumPermeability /
           (w / 2 / l + 2 / std::numbers::pi * (1 + log(std::numbers::pi * h / 4 / l)));
}

double ReluctanceMuehlethalerModel::get_reluctance_type_1(double l, double w, double h) {
    double basicReluctance = get_basic_reluctance(l, w, h);
    return 1 / (1 / (basicReluctance + basicReluctance) + 1 / (basicReluctance + basicReluctance));
}

AirGapReluctanceOutput ReluctanceMuehlethalerModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_shape()) {
        throw std::runtime_error("Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw std::runtime_error("Gap Section Dimensions are not set");
    }
    if (!gapInfo.get_distance_closest_normal_surface()) {
        throw std::runtime_error("Gap Distance Closest Normal Surface is not set");
    }
    auto gapShape = *(gapInfo.get_shape());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    auto distanceClosestNormalSurface = *(gapInfo.get_distance_closest_normal_surface());
    double reluctance;
    double fringingFactor = 1;
    auto gapSectionWidth = gapSectionDimensions[0];
    auto gapSectionDepth = gapSectionDimensions[1];

    if (gapShape == ColumnShape::ROUND) {
        double gammaR = get_reluctance_type_1(gapLength / 2, gapSectionWidth / 2, distanceClosestNormalSurface) /
                         (gapLength / constants.vacuumPermeability / (gapSectionWidth / 2));
        reluctance = pow(gammaR, 2) * gapLength /
                     (constants.vacuumPermeability * std::numbers::pi * pow(gapSectionWidth / 2, 2));
        fringingFactor = 1 / gammaR;
    }
    else {
        double gammaX = get_reluctance_type_1(gapLength / 2, gapSectionWidth, distanceClosestNormalSurface) /
                         (gapLength / constants.vacuumPermeability / gapSectionWidth);
        double gammaY = get_reluctance_type_1(gapLength / 2, gapSectionDepth, distanceClosestNormalSurface) /
                         (gapLength / constants.vacuumPermeability / gapSectionDepth);
        double gamma = gammaX * gammaY;
        reluctance = gamma * gapLength / (constants.vacuumPermeability * gapSectionDepth * gapSectionWidth);
        fringingFactor = 1 / gamma;
    }

    AirGapReluctanceOutput airGapReluctanceOutput;
    airGapReluctanceOutput.set_maximum_storable_magnetic_energy(get_gap_maximum_storable_energy(gapInfo, fringingFactor));
    airGapReluctanceOutput.set_reluctance(reluctance);
    airGapReluctanceOutput.set_method_used("Muehlethaler");
    airGapReluctanceOutput.set_origin(ResultOrigin::SIMULATION);
    airGapReluctanceOutput.set_fringing_factor(fringingFactor);

    return airGapReluctanceOutput;
};

AirGapReluctanceOutput ReluctanceEffectiveAreaModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    if (!gapInfo.get_shape()) {
        throw std::runtime_error("Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw std::runtime_error("Gap Section Dimensions are not set");
    }
    auto gapArea = *(gapInfo.get_area());
    auto gapShape = *(gapInfo.get_shape());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    double reluctance;
    double fringingFactor = 1;
    auto gapSectionWidth = gapSectionDimensions[0];
    auto gapSectionDepth = gapSectionDimensions[1];

    if (gapLength > 0) {
        if (gapShape == ColumnShape::ROUND) {
            fringingFactor = pow(1 + gapLength / gapSectionWidth, 2);
        }
        else {
            fringingFactor = (gapSectionDepth + gapLength) * (gapSectionWidth + gapLength) /
                              (gapSectionDepth * gapSectionWidth);
        }
    }

    reluctance = gapLength / (constants.vacuumPermeability * gapArea * fringingFactor);

    AirGapReluctanceOutput airGapReluctanceOutput;
    airGapReluctanceOutput.set_maximum_storable_magnetic_energy(get_gap_maximum_storable_energy(gapInfo, fringingFactor));
    airGapReluctanceOutput.set_reluctance(reluctance);
    airGapReluctanceOutput.set_method_used("EffectiveArea");
    airGapReluctanceOutput.set_origin(ResultOrigin::SIMULATION);
    airGapReluctanceOutput.set_fringing_factor(fringingFactor);

    return airGapReluctanceOutput;
};

AirGapReluctanceOutput ReluctanceEffectiveLengthModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    if (!gapInfo.get_shape()) {
        throw std::runtime_error("Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw std::runtime_error("Gap Section Dimensions are not set");
    }
    auto gapArea = *(gapInfo.get_area());
    auto gapShape = *(gapInfo.get_shape());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    double reluctance;
    double fringingFactor = 1;
    auto gapSectionWidth = gapSectionDimensions[0];
    auto gapSectionDepth = gapSectionDimensions[1];

    if (gapLength > 0) {
        if (gapShape == ColumnShape::ROUND) {
            fringingFactor = pow(1 + gapLength / gapSectionWidth, 2);
        }
        else {
            fringingFactor = (1 + gapLength / gapSectionDepth) * (1 + gapLength / gapSectionWidth);
        }
    }

    reluctance = gapLength / (constants.vacuumPermeability * gapArea * fringingFactor);

    AirGapReluctanceOutput airGapReluctanceOutput;
    airGapReluctanceOutput.set_maximum_storable_magnetic_energy(get_gap_maximum_storable_energy(gapInfo, fringingFactor));
    airGapReluctanceOutput.set_reluctance(reluctance);
    airGapReluctanceOutput.set_method_used("EffectiveLength");
    airGapReluctanceOutput.set_origin(ResultOrigin::SIMULATION);
    airGapReluctanceOutput.set_fringing_factor(fringingFactor);

    return airGapReluctanceOutput;
};

AirGapReluctanceOutput ReluctancePartridgeModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw std::runtime_error("Gap Section Dimensions are not set");
    }
    if (!gapInfo.get_distance_closest_normal_surface()) {
        throw std::runtime_error("Gap Distance Closest Normal Surface is not set");
    }
    auto gapArea = *(gapInfo.get_area());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    auto distanceClosestNormalSurface = *(gapInfo.get_distance_closest_normal_surface());
    double reluctance;
    double fringingFactor = 1;

    if (gapLength > 0) {
        fringingFactor = 1 + gapLength / sqrt(gapArea) * log(2 * 2 * distanceClosestNormalSurface / gapLength);
    }

    reluctance = gapLength / (constants.vacuumPermeability * gapArea * fringingFactor);

    AirGapReluctanceOutput airGapReluctanceOutput;
    airGapReluctanceOutput.set_maximum_storable_magnetic_energy(get_gap_maximum_storable_energy(gapInfo, fringingFactor));
    airGapReluctanceOutput.set_reluctance(reluctance);
    airGapReluctanceOutput.set_method_used("Partridge");
    airGapReluctanceOutput.set_origin(ResultOrigin::SIMULATION);
    airGapReluctanceOutput.set_fringing_factor(fringingFactor);

    return airGapReluctanceOutput;
};

AirGapReluctanceOutput ReluctanceStengleinModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    if (!gapInfo.get_shape()) {
        throw std::runtime_error("Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw std::runtime_error("Gap Section Dimensions are not set");
    }
    if (!gapInfo.get_distance_closest_normal_surface()) {
        throw std::runtime_error("Gap Distance Closest Normal Surface is not set");
    }
    if (!gapInfo.get_coordinates()) {
        throw std::runtime_error("Gap Corrdinates are not set");
    }
    if (!gapInfo.get_distance_closest_parallel_surface()) {
        throw std::runtime_error("Gap Distance Closest Parallel Surface is not set");
    }
    auto gapArea = *(gapInfo.get_area());
    auto gapCoordinates = *(gapInfo.get_coordinates());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    auto distanceClosestNormalSurface = *(gapInfo.get_distance_closest_normal_surface());
    auto distanceClosestParallelSurface = *(gapInfo.get_distance_closest_parallel_surface());
    double reluctance;
    double fringingFactor = 1;
    auto gapSectionWidth = gapSectionDimensions[0];

    if (gapLength > 0) {
        double c = gapSectionWidth / 2 + distanceClosestParallelSurface;
        double b = gapSectionWidth / 2 + 0.001;
        double l1 = distanceClosestNormalSurface * 2;
        double lg = gapLength;
        double rc = gapSectionWidth / 2;
        double rx = gapSectionWidth / 2;
        double aux1 = 1 + 2. / sqrt(std::numbers::pi) * lg / (2 * rc) * log(2.1 * rx / lg);
        double aux2 = 1. / 6. * (pow(c, 2) + 2 * c * b + pow(b, 2)) / pow(b, 2);

        double gamma = aux1 + (aux2 - aux1) * pow(lg / l1, 2 * std::numbers::pi);

        fringingFactor = alpha(rx, l1, lg) * pow(gapCoordinates[1] / l1, 2) + gamma;
    }

    reluctance = gapLength / (constants.vacuumPermeability * gapArea * fringingFactor);

    AirGapReluctanceOutput airGapReluctanceOutput;
    airGapReluctanceOutput.set_maximum_storable_magnetic_energy(get_gap_maximum_storable_energy(gapInfo, fringingFactor));
    airGapReluctanceOutput.set_reluctance(reluctance);
    airGapReluctanceOutput.set_method_used("Stenglein");
    airGapReluctanceOutput.set_origin(ResultOrigin::SIMULATION);
    airGapReluctanceOutput.set_fringing_factor(fringingFactor);

    return airGapReluctanceOutput;
};

AirGapReluctanceOutput ReluctanceClassicModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    auto gapArea = *(gapInfo.get_area());
    double reluctance;
    double fringingFactor = 1;

    reluctance = gapLength / (constants.vacuumPermeability * gapArea);

    AirGapReluctanceOutput airGapReluctanceOutput;
    airGapReluctanceOutput.set_maximum_storable_magnetic_energy(get_gap_maximum_storable_energy(gapInfo, fringingFactor));
    airGapReluctanceOutput.set_reluctance(reluctance);
    airGapReluctanceOutput.set_method_used("Classic");
    airGapReluctanceOutput.set_origin(ResultOrigin::SIMULATION);
    airGapReluctanceOutput.set_fringing_factor(fringingFactor);

    return airGapReluctanceOutput;
};

AirGapReluctanceOutput ReluctanceBalakrishnanModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw std::runtime_error("Gap Area is not set");
    }
    auto gapArea = *(gapInfo.get_area());
    double reluctance;
    double fringingFactor = 1;
    auto distanceClosestNormalSurface = *(gapInfo.get_distance_closest_normal_surface());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    auto gapSectionDepth = gapSectionDimensions[1];

    reluctance = 1. / (constants.vacuumPermeability *
                       (gapArea / gapLength +
                        2. * gapSectionDepth / std::numbers::pi *
                            (1 + log(std::numbers::pi * distanceClosestNormalSurface / (2 * gapLength)))));

    if (gapLength > 0) {
        fringingFactor = gapLength / (constants.vacuumPermeability * gapArea * reluctance);
    }

    AirGapReluctanceOutput airGapReluctanceOutput;
    airGapReluctanceOutput.set_maximum_storable_magnetic_energy(get_gap_maximum_storable_energy(gapInfo, fringingFactor));
    airGapReluctanceOutput.set_reluctance(reluctance);
    airGapReluctanceOutput.set_method_used("Balakrishnan");
    airGapReluctanceOutput.set_origin(ResultOrigin::SIMULATION);
    airGapReluctanceOutput.set_fringing_factor(fringingFactor);

    return airGapReluctanceOutput;
};

std::shared_ptr<ReluctanceModel> ReluctanceModel::factory(std::map<std::string, std::string> models) {
    return factory(magic_enum::enum_cast<OpenMagnetics::ReluctanceModels>(models["gapReluctance"]).value());
}

std::shared_ptr<ReluctanceModel> ReluctanceModel::factory(ReluctanceModels modelName) {
    if (modelName == ReluctanceModels::ZHANG) {
        return std::make_shared<ReluctanceZhangModel>();
    }
    else if (modelName == ReluctanceModels::PARTRIDGE) {
        return std::make_shared<ReluctancePartridgeModel>();
    }
    else if (modelName == ReluctanceModels::EFFECTIVE_AREA) {
        return std::make_shared<ReluctanceEffectiveAreaModel>();
    }
    else if (modelName == ReluctanceModels::EFFECTIVE_LENGTH) {
        return std::make_shared<ReluctanceEffectiveLengthModel>();
    }
    else if (modelName == ReluctanceModels::MUEHLETHALER) {
        return std::make_shared<ReluctanceMuehlethalerModel>();
    }
    else if (modelName == ReluctanceModels::STENGLEIN) {
        return std::make_shared<ReluctanceStengleinModel>();
    }
    else if (modelName == ReluctanceModels::BALAKRISHNAN) {
        return std::make_shared<ReluctanceBalakrishnanModel>();
    }
    else if (modelName == ReluctanceModels::CLASSIC) {
        return std::make_shared<ReluctanceClassicModel>();
    }

    else
        throw std::runtime_error("Unknown Reluctance mode, available options are: {ZHANG, PARTRIDGE, EFFECTIVE_AREA, "
                                 "EFFECTIVE_LENGTH, MUEHLETHALER, STENGLEIN, BALAKRISHNAN, CLASSIC}");
}
std::shared_ptr<ReluctanceModel> ReluctanceModel::factory() {
    auto defaults = Defaults();
    return factory(defaults.reluctanceModelDefault);
}
} // namespace OpenMagnetics
