#include "physical_models/Reluctance.h"
#include "physical_models/MagneticEnergy.h"
#include "support/Utils.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <streambuf>
#include <vector>
#include "support/Exceptions.h"
#include "support/Logger.h"

namespace OpenMagnetics {


double ReluctanceModel::get_ungapped_core_reluctance(const Core& core, double initialPermeability) {
    auto constants = Constants();
    double absolutePermeability = constants.vacuumPermeability * initialPermeability;
    double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();
    double effectiveLength = core.get_processed_description()->get_effective_parameters().get_effective_length();

    double reluctanceCore = effectiveLength / (absolutePermeability * effectiveArea);
    return reluctanceCore;
}

double ReluctanceModel::get_air_cored_reluctance(Bobbin bobbin) {
    if (!bobbin.get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed");
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
        throw NotImplementedException("Column shape not implemented yet");
    }

    auto constants = Constants();
    double absolutePermeability = constants.vacuumPermeability;

    double reluctanceAirCore = airLength / (absolutePermeability * airArea);
    return reluctanceAirCore;
}

double ReluctanceModel::get_ungapped_core_reluctance(Core core, std::optional<OperatingPoint> operatingPoint) {
    InitialPermeability initialPermeability;

    auto coreMaterial = core.resolve_material();

    double initialPermeabilityValue;
    if (operatingPoint) {
        double temperature = operatingPoint->get_conditions().get_ambient_temperature(); // TODO: Use a future calculated temperature
        _magneticFluxDensitySaturation = core.get_magnetic_flux_density_saturation(temperature, true);
        initialPermeabilityValue = initialPermeability.get_initial_permeability(coreMaterial, operatingPoint.value());
    }
    else {
        _magneticFluxDensitySaturation = core.get_magnetic_flux_density_saturation(true);
        initialPermeabilityValue = initialPermeability.get_initial_permeability(coreMaterial);
    }

    return get_ungapped_core_reluctance(core, initialPermeabilityValue);
}


MagnetizingInductanceOutput ReluctanceModel::get_core_reluctance(Core core, std::optional<OperatingPoint> operatingPoint) {
    auto ungappedCoreReluctance = get_ungapped_core_reluctance(core, operatingPoint);
    auto magnetizingInductanceOutput = get_gapping_reluctance(core);

    if (std::isnan(ungappedCoreReluctance)) {
        throw NaNResultException("Core Reluctance must be a number, not NaN");
    }
    double calculatedReluctance = ungappedCoreReluctance + magnetizingInductanceOutput.get_gapping_reluctance().value();
    if (std::isnan(calculatedReluctance)) {
        throw NaNResultException("Reluctance must be a number, not NaN");
    }

    if (operatingPoint) {
        magnetizingInductanceOutput.set_maximum_magnetic_energy_core(MagneticEnergy::get_ungapped_core_maximum_magnetic_energy(core, operatingPoint));
    }
    magnetizingInductanceOutput.set_core_reluctance(calculatedReluctance);
    magnetizingInductanceOutput.set_ungapped_core_reluctance(ungappedCoreReluctance);

    return magnetizingInductanceOutput;
}

MagnetizingInductanceOutput ReluctanceModel::get_core_reluctance(Core core, double initialPermeability) {
    auto ungappedCoreReluctance = get_ungapped_core_reluctance(core, initialPermeability);

    auto magnetizingInductanceOutput = get_gapping_reluctance(core);
    double calculatedReluctance = ungappedCoreReluctance + magnetizingInductanceOutput.get_gapping_reluctance().value();

    magnetizingInductanceOutput.set_core_reluctance(calculatedReluctance);
    magnetizingInductanceOutput.set_ungapped_core_reluctance(ungappedCoreReluctance);

    return magnetizingInductanceOutput;
}

MagnetizingInductanceOutput ReluctanceModel::get_gapping_reluctance(Core core) {
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
            if (gapColumn.get_type() == ColumnType::LATERAL) {
                calculatedLateralReluctance += 1 / gapReluctance.get_reluctance();
            }
            else {
                calculatedCentralReluctance += gapReluctance.get_reluctance();
            }
            maximumFringingFactor = std::max(maximumFringingFactor, gapReluctance.get_fringing_factor());
            maximumStorableMagneticEnergyGapping += gapReluctance.get_maximum_storable_magnetic_energy();
            if (gapReluctance.get_fringing_factor() < 1) {
                std::cerr << "fringing_factor " << gapReluctance.get_fringing_factor() << std::endl;
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

/**
 * @brief Calculate air gap reluctance using Zhang's improved method for air-gap inductors
 *
 * Reference: X. Zhang, F. Xiao, R. Wang, X. Fan, H. Wang,
 *            "Improved Calculation Method for Inductance Value of the Air-Gap Inductor",
 *            IEEE 1st China International Youth Conference on Electrical Engineering (CIYCEE), 2020.
 *            https://ieeexplore.ieee.org/document/9332553
 *
 * The total reluctance is modeled as a parallel combination of internal reluctance (direct flux path)
 * and fringing reluctance (fringing flux path around the gap edges):
 *
 *   R_g = R_fr || R_in                                                    [Eq. 11]
 *
 * where:
 *   - R_in = d_i / (μ₀ * A_c)  (internal reluctance, uniform field)       [Eq. 9]
 *   - R_fr = π / (μ₀ * C * ln((2*h + d_i) / d_i))  (fringing reluctance)  [Eq. 10]
 *
 * with C being the perimeter of the cross-section and h the winding window height.
 * The fringing field is modeled using an equivalent current source at the gap edge (Fig. 6).
 *
 * The fringing factor F is computed as the ratio of effective to geometric reluctance:
 *   F = l_g / (μ₀ * A_g * R_total)
 *
 * @param gapInfo CoreGap object containing gap geometry parameters
 * @return AirGapReluctanceOutput with calculated reluctance and fringing factor
 */
AirGapReluctanceOutput ReluctanceZhangModel::get_gap_reluctance(CoreGap gapInfo) {
    double perimeter = 0;
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Area is not set");
    }
    if (!gapInfo.get_shape()) {
        throw GapException(ErrorCode::GAP_SHAPE_NOT_SET, "Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Section Dimensions are not set");
    }
    if (!gapInfo.get_distance_closest_normal_surface()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Distance Closest Normal Surface is not set");
    }
    auto gapArea = *(gapInfo.get_area());
    auto gapShape = *(gapInfo.get_shape());
    auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
    auto gapSectionWidth = gapSectionDimensions[0];
    auto gapSectionDepth = gapSectionDimensions[1];
    auto distanceClosestNormalSurface = std::max(*(gapInfo.get_distance_closest_normal_surface()), gapSectionWidth);
    auto reluctanceInternal = gapLength / (constants.vacuumPermeability * gapArea);
    double reluctanceFringing = 0;
    double fringingFactor = 1;

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
            throw NaNResultException("reluctanceInternal cannot be 0 or NaN");
        }

        if (std::isnan(reluctanceFringing) || reluctanceFringing == 0) {
            throw NaNResultException("reluctanceFringing cannot be 0 or NaN");
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


/**
 * @brief Helper function to compute the basic reluctance for a single fringing path
 *
 * Reference: J. Mühlethaler, J.W. Kolar, A. Ecklebe,
 *            "A Novel Approach for 3D Air Gap Reluctance Calculations",
 *            8th International Conference on Power Electronics - ECCE Asia, Jeju, 2011.
 *            https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/10_A_Novel_Approach_ECCEAsia2011_01.pdf
 *            (Referenced by Zhang 2020 as [5])
 *
 * The basic reluctance element accounts for both direct flux path and fringing:
 *
 *   R_basic = 1 / (μ₀ * (w/(2*l) + (2/π) * (1 + ln(π*h / (4*l)))))
 *
 * @param l Half of gap length (l_g/2)
 * @param w Width dimension of the gap cross-section
 * @param h Distance to closest perpendicular surface (winding window height)
 * @return Basic reluctance value in H⁻¹
 */
double ReluctanceMuehlethalerModel::get_basic_reluctance(double l, double w, double h) {
    auto constants = Constants();
    return 1 / constants.vacuumPermeability /
           (w / 2 / l + 2 / std::numbers::pi * (1 + log(std::numbers::pi * h / 4 / l)));
}

/**
 * @brief Compute Type 1 reluctance (full gap surrounded by core on both sides)
 *
 * Reference: J. Mühlethaler, J.W. Kolar, A. Ecklebe,
 *            "A Novel Approach for 3D Air Gap Reluctance Calculations",
 *            8th International Conference on Power Electronics - ECCE Asia, Jeju, 2011.
 *            (Referenced by Zhang 2020 as [5])
 *
 * Combines four basic reluctance elements for a symmetric gap configuration:
 *   R_type1 = 1 / (1/(R_b + R_b) + 1/(R_b + R_b)) = R_b
 *
 * @param l Half of gap length (l_g/2)
 * @param w Width dimension
 * @param h Distance to closest perpendicular surface
 * @return Type 1 reluctance in H⁻¹
 */
double ReluctanceMuehlethalerModel::get_reluctance_type_1(double l, double w, double h) {
    double basicReluctance = get_basic_reluctance(l, w, h);
    return 1 / (1 / (basicReluctance + basicReluctance) + 1 / (basicReluctance + basicReluctance));
}

/**
 * @brief Calculate air gap reluctance using Mühlethaler's 3D approach
 *
 * Reference: J. Mühlethaler, J.W. Kolar, A. Ecklebe,
 *            "A Novel Approach for 3D Air Gap Reluctance Calculations",
 *            8th International Conference on Power Electronics - ECCE Asia, Jeju, 2011.
 *            https://www.pes-publications.ee.ethz.ch/uploads/tx_ethpublications/10_A_Novel_Approach_ECCEAsia2011_01.pdf
 *            (Referenced by Zhang 2020 as [5])
 *
 * This method extends 2D reluctance analysis to 3D by decomposing the air gap into
 * multiple reluctance elements and using superposition. From Zhang Eq. 6-7:
 *
 *   R_x = R_0x / γ_x   where R_0x = d / (μ₀ * w_x)
 *   R_y = R_0y / γ_y   where R_0y = d / (μ₀ * w_y)
 *   R_3D = γ_x * γ_y * d / (μ₀ * w_x * w_y)                               [Zhang Eq. 7]
 *
 * The fringing factor is F = 1/γ where γ = γ_x * γ_y
 *
 * @param gapInfo CoreGap object containing gap geometry parameters
 * @return AirGapReluctanceOutput with calculated reluctance and fringing factor
 */
AirGapReluctanceOutput ReluctanceMuehlethalerModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_shape()) {
        throw GapException(ErrorCode::GAP_SHAPE_NOT_SET, "Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Section Dimensions are not set");
    }
    if (!gapInfo.get_distance_closest_normal_surface()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Distance Closest Normal Surface is not set");
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

/**
 * @brief Calculate air gap reluctance using the effective area method
 *
 * This simplified approach accounts for fringing by increasing the effective cross-sectional
 * area of the gap. The fringing factor expands the area by the gap length in each dimension.
 *
 * For circular cross-section:
 *   A_eff = π * (r + l_g)² = π * r² * (1 + l_g/r)²
 *   F = (1 + l_g/d)²   where d is the column diameter
 *
 * For rectangular cross-section:
 *   A_eff = (w + l_g) * (d + l_g)
 *   F = (1 + l_g/w) * (1 + l_g/d)
 *
 * The effective reluctance is:
 *   R = l_g / (μ₀ * A * F)
 *
 * @param gapInfo CoreGap object containing gap geometry parameters
 * @return AirGapReluctanceOutput with calculated reluctance and fringing factor
 */
AirGapReluctanceOutput ReluctanceEffectiveAreaModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Area is not set");
    }
    if (!gapInfo.get_shape()) {
        throw GapException(ErrorCode::GAP_SHAPE_NOT_SET, "Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Section Dimensions are not set");
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

/**
 * @brief Calculate air gap reluctance using the effective length method
 *
 * Similar to the effective area method, this approach accounts for fringing flux by
 * assuming the flux spreads out from the gap edges. The fringing factor is computed
 * to account for the reduction in effective reluctance due to fringing.
 *
 * For circular cross-section:
 *   F = (1 + l_g/d)²   where d is the column diameter
 *
 * For rectangular cross-section:
 *   F = (1 + l_g/d) * (1 + l_g/w)   where d is depth, w is width
 *
 * The effective reluctance is:
 *   R = l_g / (μ₀ * A * F)
 *
 * Note: This method produces identical results to EffectiveArea for most geometries.
 *
 * @param gapInfo CoreGap object containing gap geometry parameters
 * @return AirGapReluctanceOutput with calculated reluctance and fringing factor
 */
AirGapReluctanceOutput ReluctanceEffectiveLengthModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Area is not set");
    }
    if (!gapInfo.get_shape()) {
        throw GapException(ErrorCode::GAP_SHAPE_NOT_SET, "Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Section Dimensions are not set");
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

/**
 * @brief Calculate air gap reluctance using Partridge's formula
 *
 * Reference: Referenced by Zhang 2020 as [3] (McLyman), Eq. 1:
 *   "Reference [3] provided a correction factor for C-type and E-type core to describe
 *    the influence of the air gap fringing flux on the inductance value."
 *
 * This method uses a logarithmic fringing factor formula that accounts for flux spreading
 * around the gap edges based on the ratio of winding window height to gap length:
 *
 *   F = 1 + (l_g / √A_c) * ln(2 * H_w / d)                               [Zhang Eq. 1]
 *
 * where:
 *   - d = gap length (l_g)
 *   - A_c = cross-sectional area of the core
 *   - H_w = height of the winding window
 *
 * The reluctance is then computed as:
 *   R = l_g / (μ₀ * A * F)
 *
 * @param gapInfo CoreGap object containing gap geometry parameters
 * @return AirGapReluctanceOutput with calculated reluctance and fringing factor
 */
AirGapReluctanceOutput ReluctancePartridgeModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Area is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Section Dimensions are not set");
    }
    if (!gapInfo.get_distance_closest_normal_surface()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Distance Closest Normal Surface is not set");
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

/**
 * @brief Calculate air gap reluctance using Stenglein's method for large air gaps
 *
 * Reference: E. Stenglein, D. Kuebrich, M. Albach, T. Duerbaum,
 *            "The Reluctance of Large Air Gaps in Ferrite Cores",
 *            17th European Conference on Power Electronics and Applications (EPE'16 ECCE Europe), 2016.
 *            https://ieeexplore.ieee.org/document/7695271
 *
 * This method is specifically designed for large air gaps where the gap length is
 * comparable to the winding window dimensions. It accounts for both gap position
 * and proximity to core boundaries.
 *
 * For lM = 0 (centered gap), the ratio Ag/Ac is given by Eq. 12:
 *
 *   γ(lg) = 1 + (2/√π) * (lg/(2*rc)) * ln(2.1*rx/lg) + (aux2 - aux1) * (lg/l1)^(2π)
 *
 * where aux2 = (1/6) * (c² + 2cb + 3b²) / rc²   (from Eq. 11 limit case)
 *
 * For position-dependent gaps, Eq. 13-14 add the α factor:
 *   Ag/Ac = α(lg) * (lM/l1)² + γ(lg)                                      [Eq. 13]
 *
 * with polynomial coefficients from least-squares fitting (Eq. 15-17):
 *   u = 42.7 * rx/l1 - 50.2
 *   v = -55.4 * rx/l1 + 71.6
 *   w = 0.88 * rx/l1 - 0.80
 *   α(lg) = u*(lg/l1)² + v*(lg/l1) + w                                    [Eq. 14]
 *
 * @param gapInfo CoreGap object containing gap geometry parameters
 * @return AirGapReluctanceOutput with calculated reluctance and fringing factor
 */
AirGapReluctanceOutput ReluctanceStengleinModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Area is not set");
    }
    if (!gapInfo.get_shape()) {
        throw GapException(ErrorCode::GAP_SHAPE_NOT_SET, "Gap Shape is not set");
    }
    if (!gapInfo.get_section_dimensions()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Section Dimensions are not set");
    }
    if (!gapInfo.get_distance_closest_normal_surface()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Distance Closest Normal Surface is not set");
    }
    if (!gapInfo.get_coordinates()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Coordinates are not set");
    }
    if (!gapInfo.get_distance_closest_parallel_surface()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Distance Closest Parallel Surface is not set");
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

/**
 * @brief Calculate air gap reluctance using the classic formula (no fringing)
 *
 * Reference: Standard magnetic circuit theory.
 *            https://en.wikipedia.org/wiki/Magnetic_reluctance
 *
 * The classic reluctance formula assumes a uniform magnetic field with no fringing:
 *
 *   R = l_g / (μ₀ * A)
 *
 * where:
 *   - l_g = gap length
 *   - A = gap cross-sectional area
 *   - μ₀ = vacuum permeability (4π × 10⁻⁷ H/m)
 *
 * The fringing factor is fixed at 1 (no fringing compensation).
 *
 * Note: This method typically overestimates reluctance (underestimates inductance)
 * because it ignores the fringing flux that exists in all practical air gaps.
 *
 * @param gapInfo CoreGap object containing gap geometry parameters
 * @return AirGapReluctanceOutput with calculated reluctance and fringing factor = 1
 */
AirGapReluctanceOutput ReluctanceClassicModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Area is not set");
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

/**
 * @brief Calculate air gap reluctance using Balakrishnan's Schwarz-Christoffel method
 *
 * Reference: A. Balakrishnan, W. T. Joines, T. G. Wilson,
 *            "Air-gap reluctance and inductance calculations for magnetic circuits
 *            using a Schwarz-Christoffel transformation",
 *            IEEE Transactions on Power Electronics, vol. 12, no. 4, pp. 654-663, July 1997.
 *            https://ieeexplore.ieee.org/document/602560
 *            (Referenced by Zhang 2020 as [4])
 *
 * Uses conformal mapping (Schwarz-Christoffel transformation) for analytical solution
 * of fringing flux. Per-unit-length reluctance expressions from Section V (Eq. 14-17):
 *
 * For post-plate configuration (Fig. 3a, Eq. 14):
 *   R_a = 1 / (μ₀ * (w/(2*d) + (2/π) * (1 + ln(π*h / (4*d)))))
 *
 * For post-post configuration (Fig. 3b, Eq. 15, used here):
 *   R_b = 1 / (μ₀ * (w/d + (4/π) * (1 + ln(π*h / (4*d)))))
 *
 * The implementation uses Eq. 15 form:
 *   R = 1 / (μ₀ * (A/lg + (2*depth/π) * (1 + ln(π*h / (2*lg)))))
 *
 * The fringing factor F = lg / (μ₀ * A * R)
 *
 * @param gapInfo CoreGap object containing gap geometry parameters
 * @return AirGapReluctanceOutput with calculated reluctance and fringing factor
 */
AirGapReluctanceOutput ReluctanceBalakrishnanModel::get_gap_reluctance(CoreGap gapInfo) {
    auto constants = Constants();
    auto gapLength = gapInfo.get_length();
    if (!gapInfo.get_area()) {
        throw GapException(ErrorCode::GAP_INVALID_DIMENSIONS, "Gap Area is not set");
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
    ReluctanceModels reluctanceModelEnum;
    from_json(models["gapReluctance"], reluctanceModelEnum);
    return factory(reluctanceModelEnum);
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
        throw ModelNotAvailableException("Unknown Reluctance mode, available options are: {ZHANG, PARTRIDGE, EFFECTIVE_AREA, "
                                 "EFFECTIVE_LENGTH, MUEHLETHALER, STENGLEIN, BALAKRISHNAN, CLASSIC}");
}
std::shared_ptr<ReluctanceModel> ReluctanceModel::factory() {
    return factory(defaults.reluctanceModelDefault);
}

double ReluctanceModel::get_gapping_by_fringing_factor(Core core, double fringingFactor) {
    auto centralColumns = core.find_columns_by_type(ColumnType::CENTRAL);
    if (centralColumns.size() == 0) {
        throw CoreNotProcessedException("No columns found in core");
    }
    double gapLength = centralColumns[0].get_height();
    double gapIncrease = gapLength / 2;
    size_t timeout = 100;
    while (true) {
        core.set_gap_length(gapLength);
        auto calculatedFringingFactor = get_core_reluctance(core).get_maximum_fringing_factor().value();
        if ((fabs(calculatedFringingFactor - fringingFactor) / fringingFactor) < 0.001) {
            break;
        }
        if (calculatedFringingFactor < fringingFactor) {
            gapLength += gapIncrease;
            if (gapLength > centralColumns[0].get_height()) {
                return centralColumns[0].get_height() / 2;
            }
        }
        else {
            gapLength -= gapIncrease;
        }
        gapLength = roundFloat(gapLength, 6);
        gapIncrease = std::max(gapIncrease / 2, constants.residualGap);
        timeout--;
        if (timeout <= 0) {
            break;
        }
    }

    return gapLength;
}

} // namespace OpenMagnetics
