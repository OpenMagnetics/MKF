#include "physical_models/MagneticField.h"
#include <numbers>
#include "physical_models/StrayCapacitance.h"
#include "support/Painter.h"
#include "support/CoilMesher.h"
#include "MAS.hpp"
#include "support/Utils.h"
#include "json.hpp"
#ifdef ENABLE_MATPLOTPP
#include <matplot/matplot.h>
#endif
#include <cfloat>
#include <chrono>
#include <thread>
#include <filesystem>
#include <set>
#include <cmath>
#include <iostream>
#include "support/Exceptions.h"


namespace OpenMagnetics {

ComplexField PainterInterface::calculate_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex) {
    if (!operatingPoint.get_excitations_per_winding()[0].get_current()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Current is missing in excitation");
    }
    for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        if (!operatingPoint.get_excitations_per_winding()[windingIndex].get_current()->get_harmonics()) {
            auto current = operatingPoint.get_excitations_per_winding()[windingIndex].get_current().value();
            if (!current.get_waveform()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Waveform is missing from current");
            }
            auto sampledWaveform = Inputs::calculate_sampled_waveform(current.get_waveform().value(), operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            auto harmonics = Inputs::calculate_harmonics_data(sampledWaveform, operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            current.set_harmonics(harmonics);
            if (!current.get_processed()) {
                auto processed = Inputs::calculate_processed_data(harmonics, sampledWaveform, true);
                current.set_processed(processed);
            }
            else {
                if (!current.get_processed()->get_rms()) {
                    auto processed = Inputs::calculate_processed_data(harmonics, sampledWaveform, true);
                    current.set_processed(processed);
                }
            }
            operatingPoint.get_mutable_excitations_per_winding()[windingIndex].set_current(current);
        }
    }

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];

    bool includeFringing = settings.get_painter_include_fringing();
    bool mirroringDimension = settings.get_painter_mirroring_dimension();

    size_t numberPointsX = settings.get_painter_number_points_x();
    size_t numberPointsY = settings.get_painter_number_points_y();
    Field inducedField = CoilMesher::generate_mesh_induced_grid(magnetic, frequency, numberPointsX, numberPointsY, true).first;

    auto modelOverride = settings.get_painter_magnetic_field_strength_model();
    // Use painter override if set, otherwise use the simulation magnetic field strength model
    auto magneticFieldModel = modelOverride.value_or(settings.get_magnetic_field_strength_model());
    // Use the configured fringing effect model from settings
    auto fringingEffectModel = settings.get_magnetic_field_strength_fringing_effect_model();
    MagneticField magneticField(magneticFieldModel, fringingEffectModel);
    settings.set_magnetic_field_include_fringing(includeFringing);
    settings.set_magnetic_field_mirroring_dimension(mirroringDimension);
    ComplexField field;
    {
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField);
        field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];

    }
    auto turns = magnetic.get_coil().get_turns_description().value();

    if (turns[0].get_additional_coordinates()) {
        for (size_t turnIndex = 0; turnIndex < turns.size(); ++turnIndex) {
            if (turns[turnIndex].get_additional_coordinates()) {
                turns[turnIndex].set_coordinates(turns[turnIndex].get_additional_coordinates().value()[0]);
            }
        }
        magnetic.get_mutable_coil().set_turns_description(turns);
        auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField);
        auto additionalField = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
        for (size_t pointIndex = 0; pointIndex < field.get_data().size(); ++pointIndex) {
            field.get_mutable_data()[pointIndex].set_real(field.get_mutable_data()[pointIndex].get_real() + additionalField.get_mutable_data()[pointIndex].get_real());
            field.get_mutable_data()[pointIndex].set_imaginary(field.get_mutable_data()[pointIndex].get_imaginary() + additionalField.get_mutable_data()[pointIndex].get_imaginary());
        }
    }


    return field;
}


// ==================== SDF Primitives ====================

double PainterInterface::sdf_circle(double px, double py, double cx, double cy, double radius) {
    return hypot(px - cx, py - cy) - radius;
}

double PainterInterface::sdf_box(double px, double py, double cx, double cy, double halfW, double halfH) {
    double dx = fabs(px - cx) - halfW;
    double dy = fabs(py - cy) - halfH;
    return hypot(std::max(dx, 0.0), std::max(dy, 0.0)) + std::min(std::max(dx, dy), 0.0);
}

double PainterInterface::sdf_oriented_box(double px, double py, double cx, double cy,
                                           double halfW, double halfH, double angle) {
    double cosA = cos(angle), sinA = sin(angle);
    double lx =  cosA * (px - cx) + sinA * (py - cy);
    double ly = -sinA * (px - cx) + cosA * (py - cy);
    double dx = fabs(lx) - halfW;
    double dy = fabs(ly) - halfH;
    return hypot(std::max(dx, 0.0), std::max(dy, 0.0)) + std::min(std::max(dx, dy), 0.0);
}

double PainterInterface::sdf_turn(double px, double py, const Turn& turn) {
    if (!turn.get_cross_sectional_shape() || !turn.get_dimensions()) {
        return std::numeric_limits<double>::max();
    }
    auto shape = turn.get_cross_sectional_shape().value();
    auto& coords = turn.get_coordinates();
    auto dims = turn.get_dimensions().value();

    if (coords.size() < 2 || dims.empty()) {
        return std::numeric_limits<double>::max();
    }

    if (shape == TurnCrossSectionalShape::ROUND) {
        return sdf_circle(px, py, coords[0], coords[1], dims[0] / 2.0);
    } else {
        double angle = 0.0;
        if (turn.get_rotation()) {
            angle = turn.get_rotation().value() * std::numbers::pi / 180.0;
        }
        double halfWidth = dims[0] / 2.0;
        double halfHeight = dims.size() > 1 ? dims[1] / 2.0 : halfWidth;
        return sdf_oriented_box(px, py, coords[0], coords[1],
                                halfWidth, halfHeight, angle);
    }
}

// ==================== Energy Density Models ====================

double PainterInterface::energy_density_bipolar(double px, double py,
                                                 const Turn& t1, const Turn& t2,
                                                 double voltageDrop, double epsilonEff) {
    auto params = StrayCapacitance::compute_bipolar_params(t1, t2);
    double energy = StrayCapacitance::bipolar_energy_density_at_point(px, py, params, voltageDrop, epsilonEff);
    // Prevent infinity/NaN
    if (!std::isfinite(energy) || energy < 0) return 0;
    if (energy > 1e300) return 1e300;
    return energy;
}

double PainterInterface::energy_density_parallel_plate(double px, double py,
                                                         const Turn& t1, const Turn& t2,
                                                         double voltageDrop, double epsilonEff) {
    auto& coords1 = t1.get_coordinates();
    auto& coords2 = t2.get_coordinates();
    if (coords1.size() < 2 || coords2.size() < 2) return 0;

    double cx1 = coords1[0], cy1 = coords1[1];
    double cx2 = coords2[0], cy2 = coords2[1];
    double centerDist = hypot(cx2 - cx1, cy2 - cy1);
    if (centerDist < 1e-15) return 0;

    double axisAngle = atan2(cy2 - cy1, cx2 - cx1);

    auto getHalfExtent = [](const Turn& t, double ax_angle) -> double {
        if (!t.get_dimensions()) return 0;
        auto d = t.get_dimensions().value();
        if (d.empty()) return 0;

        double turnAngle = 0.0;
        if (t.get_rotation()) {
            turnAngle = t.get_rotation().value() * std::numbers::pi / 180.0;
        }
        double relAngle = ax_angle - turnAngle;
        double hw = d[0] / 2.0;
        double hh = d.size() > 1 ? d[1] / 2.0 : hw;
        return fabs(hw * cos(relAngle)) + fabs(hh * sin(relAngle));
    };

    double halfExt1 = getHalfExtent(t1, axisAngle);
    double halfExt2 = getHalfExtent(t2, axisAngle);
    double gap = centerDist - halfExt1 - halfExt2;
    if (gap < 1e-15) gap = 1e-15;

    double E = voltageDrop / gap;
    double energy = 0.5 * epsilonEff * E * E;
    // Prevent infinity/NaN
    if (!std::isfinite(energy) || energy < 0) return 0;
    if (energy > 1e300) return 1e300;
    return energy;
}

double PainterInterface::energy_density_angled_plates(double px, double py,
                                                       const Turn& t1, const Turn& t2,
                                                       double voltageDrop, double epsilonEff) {
    // Local gap model: use SDF from each turn surface, sum as local gap
    double d1 = sdf_turn(px, py, t1);
    double d2 = sdf_turn(px, py, t2);
    double localGap = d1 + d2;
    if (localGap < 1e-15) localGap = 1e-15;

    double E = voltageDrop / localGap;
    double energy = 0.5 * epsilonEff * E * E;
    // Prevent infinity/NaN
    if (!std::isfinite(energy) || energy < 0) return 0;
    if (energy > 1e300) return 1e300;
    return energy;
}

double PainterInterface::energy_density_round_rect(double px, double py,
                                                     const Turn& roundTurn, const Turn& rectTurn,
                                                     double voltageDrop, double epsilonEff) {
    double dRound = sdf_turn(px, py, roundTurn);
    double dRect  = sdf_turn(px, py, rectTurn);
    double localGap = dRound + dRect;
    if (localGap < 1e-15) localGap = 1e-15;

    double totalDist = dRound + dRect;
    double blendRect = (totalDist > 1e-15) ? dRound / totalDist : 0.5;

    double E_plate = voltageDrop / localGap;
    double u_plate = 0.5 * epsilonEff * E_plate * E_plate;

    // Bipolar correction near the round conductor
    double bipolarBoost = 1.0;
    if (roundTurn.get_dimensions() && !roundTurn.get_dimensions().value().empty()) {
        double r = roundTurn.get_dimensions().value()[0] / 2.0;
        if (r > 1e-15 && dRound < r && dRound > 0) {
            bipolarBoost = localGap / std::max(1e-15, localGap - 0.3 * dRound * dRound / r);
            bipolarBoost = std::max(1.0, std::min(bipolarBoost, 3.0));
        }
    }

    double energy = u_plate * (blendRect + (1.0 - blendRect) * bipolarBoost);
    // Prevent infinity/NaN
    if (!std::isfinite(energy) || energy < 0) return 0;
    if (energy > 1e300) return 1e300;
    return energy;
}

double PainterInterface::compute_energy_density_at_pixel(double px, double py,
                                                          const Turn& t1, const Turn& t2,
                                                          double voltageDrop, double epsilonEff) {
    if (!t1.get_cross_sectional_shape() || !t2.get_cross_sectional_shape()) {
        return 0.0;
    }
    auto shape1 = t1.get_cross_sectional_shape().value();
    auto shape2 = t2.get_cross_sectional_shape().value();
    bool round1 = (shape1 == TurnCrossSectionalShape::ROUND);
    bool round2 = (shape2 == TurnCrossSectionalShape::ROUND);

    if (round1 && round2) {
        return energy_density_bipolar(px, py, t1, t2, voltageDrop, epsilonEff);
    }
    else if (!round1 && !round2) {
        double angle1 = 0.0, angle2 = 0.0;
        if (t1.get_rotation()) {
            angle1 = t1.get_rotation().value();
        }
        if (t2.get_rotation()) {
            angle2 = t2.get_rotation().value();
        }
        double angleDiff = fabs(angle1 - angle2);
        while (angleDiff >= 180.0) angleDiff -= 180.0;

        if (angleDiff < 1.0) {
            return energy_density_parallel_plate(px, py, t1, t2, voltageDrop, epsilonEff);
        } else {
            return energy_density_angled_plates(px, py, t1, t2, voltageDrop, epsilonEff);
        }
    }
    else {
        const Turn& rTurn = round1 ? t1 : t2;
        const Turn& xTurn = round1 ? t2 : t1;
        return energy_density_round_rect(px, py, rTurn, xTurn, voltageDrop, epsilonEff);
    }
}


// ==================== SDF-based Electric Field Calculation ====================

Field PainterInterface::calculate_electric_field_sdf(OperatingPoint operatingPoint,
                                                      Magnetic magnetic,
                                                      size_t harmonicIndex) {
    // Ensure core is processed
    if (!magnetic.get_core().get_processed_description()) {
        magnetic.get_mutable_core().process_data();
    }

    // --- Voltage preprocessing (same as legacy) ---
    if (!operatingPoint.get_excitations_per_winding()[0].get_voltage()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "voltage is missing in excitation");
    }
    for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        if (!operatingPoint.get_excitations_per_winding()[windingIndex].get_voltage()->get_harmonics()) {
            auto voltage = operatingPoint.get_excitations_per_winding()[windingIndex].get_voltage().value();
            if (!voltage.get_waveform()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Waveform is missing from voltage");
            }
            auto sampledWaveform = Inputs::calculate_sampled_waveform(voltage.get_waveform().value(),
                operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            auto harmonics_data = Inputs::calculate_harmonics_data(sampledWaveform,
                operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            voltage.set_harmonics(harmonics_data);
            if (!voltage.get_processed()) {
                auto processed = Inputs::calculate_processed_data(harmonics_data, sampledWaveform, true);
                voltage.set_processed(processed);
            } else if (!voltage.get_processed()->get_rms()) {
                auto processed = Inputs::calculate_processed_data(harmonics_data, sampledWaveform, true);
                voltage.set_processed(processed);
            }
            operatingPoint.get_mutable_excitations_per_winding()[windingIndex].set_voltage(voltage);
        }
    }

    bool includeFringing = settings.get_painter_include_fringing();
    int mirroringDimension = settings.get_painter_mirroring_dimension();

    auto oldCoilMesherInsideTurnsFactor = settings.get_coil_mesher_inside_turns_factor();
    // Use same factor as LEGACY method for consistent mesh generation
    settings.set_coil_mesher_inside_turns_factor(1.2);
    settings.set_magnetic_field_include_fringing(includeFringing);
    settings.set_magnetic_field_mirroring_dimension(mirroringDimension);

    auto strayCapacitanceModel = settings.get_stray_capacitance_model();
    StrayCapacitance strayCapacitance(strayCapacitanceModel);

    auto coil = magnetic.get_coil();
    auto wirePerWinding = coil.get_wires();

    // Create the mesh field (same as legacy: reuse CoilMesher for grid points)
    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_voltage()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];
    size_t numberPointsX = settings.get_painter_number_points_x();
    size_t numberPointsY = settings.get_painter_number_points_y();
    Field inducedField = CoilMesher::generate_mesh_induced_grid(magnetic, frequency, numberPointsX, numberPointsY, false, false).first;

    auto capacitanceOutput = strayCapacitance.calculate_capacitance(coil);
    auto voltageDropAmongTurnsOpt = capacitanceOutput.get_voltage_drop_among_turns();
    if (!voltageDropAmongTurnsOpt) {
        // Return the mesh field even without voltage data
        settings.set_coil_mesher_inside_turns_factor(oldCoilMesherInsideTurnsFactor);
        return inducedField;
    }
    auto voltageDropAmongTurns = voltageDropAmongTurnsOpt.value();

    if (!coil.get_turns_description()) {
        return inducedField;
    }
    auto turns = coil.get_turns_description().value();

    settings.set_coil_mesher_inside_turns_factor(oldCoilMesherInsideTurnsFactor);

    // Clear values to zero - we just want the mesh grid positions
    for (size_t i = 0; i < inducedField.get_data().size(); ++i) {
        inducedField.get_mutable_data()[i].set_value(0);
    }

    // --- Pre-compute turn pair data ---
    struct TurnPairInfo {
        size_t idx1, idx2;
        double voltageDrop;
    };
    std::vector<TurnPairInfo> pairInfos;
    std::set<std::pair<std::string, std::string>> processedPairs;

    for (auto& [firstName, aux] : voltageDropAmongTurns) {
        for (auto& [secondName, vDrop] : aux) {
            auto key = std::make_pair(firstName, secondName);
            auto revKey = std::make_pair(secondName, firstName);
            if (processedPairs.contains(key) || processedPairs.contains(revKey)) continue;
            processedPairs.insert(key);

            size_t idx1 = SIZE_MAX, idx2 = SIZE_MAX;
            for (size_t t = 0; t < turns.size(); ++t) {
                if (turns[t].get_name() == firstName) idx1 = t;
                if (turns[t].get_name() == secondName) idx2 = t;
            }
            if (idx1 == SIZE_MAX || idx2 == SIZE_MAX) continue;
            pairInfos.push_back({idx1, idx2, fabs(vDrop)});
        }
    }

    // Pre-compute bipolar params for round-round pairs (cache for performance)
    std::map<std::pair<size_t, size_t>, StrayCapacitance::BipolarParams> bipolarCache;
    for (auto& pInfo : pairInfos) {
        auto shape1 = turns[pInfo.idx1].get_cross_sectional_shape();
        auto shape2 = turns[pInfo.idx2].get_cross_sectional_shape();
        if (shape1 && shape2 &&
            shape1.value() == TurnCrossSectionalShape::ROUND &&
            shape2.value() == TurnCrossSectionalShape::ROUND) {
            bipolarCache[{pInfo.idx1, pInfo.idx2}] =
                StrayCapacitance::compute_bipolar_params(turns[pInfo.idx1], turns[pInfo.idx2]);
        }
    }

    // Vacuum permittivity as default effective permittivity
    double epsilonEff = 8.854187817e-12;

    // If no pairs found, return the empty field
    if (pairInfos.empty()) {
        settings.set_coil_mesher_inside_turns_factor(oldCoilMesherInsideTurnsFactor);
        return inducedField;
    }

    // --- Per-pixel: SDF Voronoi assignment + physics-based energy density ---
    for (size_t pointIndex = 0; pointIndex < inducedField.get_data().size(); ++pointIndex) {
        auto& datum = inducedField.get_data()[pointIndex];
        auto& point = datum.get_point();
        if (point.size() < 2) {
            inducedField.get_mutable_data()[pointIndex].set_value(0);
            continue;
        }
        double px = point[0];
        double py = point[1];

        // Find 2 nearest turns by SDF
        double dBest = DBL_MAX, dSecond = DBL_MAX;
        size_t iBest = 0, iSecond = 0;

        for (size_t t = 0; t < turns.size(); ++t) {
            double d = sdf_turn(px, py, turns[t]);
            if (d < dBest) {
                dSecond = dBest;  iSecond = iBest;
                dBest = d;        iBest = t;
            } else if (d < dSecond) {
                dSecond = d;  iSecond = t;
            }
        }

        // Only compute for pixels outside all turns
        if (dBest <= 0 || turns.size() < 2) {
            inducedField.get_mutable_data()[pointIndex].set_value(0);
            continue;
        }

        // Find the voltage drop for this (iBest, iSecond) pair
        double voltageDrop = 0;
        bool found = false;
        for (auto& pInfo : pairInfos) {
            if ((pInfo.idx1 == iBest && pInfo.idx2 == iSecond) ||
                (pInfo.idx1 == iSecond && pInfo.idx2 == iBest)) {
                voltageDrop = pInfo.voltageDrop;
                found = true;
                break;
            }
        }

        if (!found || voltageDrop < 1e-15) {
            inducedField.get_mutable_data()[pointIndex].set_value(0);
            continue;
        }

        // Compute physics-accurate energy density based on shape combination
        double u = compute_energy_density_at_pixel(px, py,
            turns[iBest], turns[iSecond], voltageDrop, epsilonEff);

        // Ensure value is valid before setting
        if (!std::isfinite(u) || u < 0) u = 0;
        if (u > 1e300) u = 1e300;

        inducedField.get_mutable_data()[pointIndex].set_value(u);
    }

    return inducedField;
}

Field PainterInterface::calculate_electric_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex, ElectricFieldVisualizationModel model) {
    if (model == ElectricFieldVisualizationModel::SDF_PHYSICS) {
        return calculate_electric_field_sdf(operatingPoint, magnetic, harmonicIndex);
    }
    // === LEGACY implementation below ===
    if (!operatingPoint.get_excitations_per_winding()[0].get_voltage()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "voltage is missing in excitation");
    }
    for (size_t windingIndex = 0; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        if (!operatingPoint.get_excitations_per_winding()[windingIndex].get_voltage()->get_harmonics()) {
            auto voltage = operatingPoint.get_excitations_per_winding()[windingIndex].get_voltage().value();
            if (!voltage.get_waveform()) {
                throw InvalidInputException(ErrorCode::MISSING_DATA, "Waveform is missing from voltage");
            }
            auto sampledWaveform = Inputs::calculate_sampled_waveform(voltage.get_waveform().value(), operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            auto harmonics = Inputs::calculate_harmonics_data(sampledWaveform, operatingPoint.get_excitations_per_winding()[windingIndex].get_frequency());
            voltage.set_harmonics(harmonics);
            if (!voltage.get_processed()) {
                auto processed = Inputs::calculate_processed_data(harmonics, sampledWaveform, true);
                voltage.set_processed(processed);
            }
            else {
                if (!voltage.get_processed()->get_rms()) {
                    auto processed = Inputs::calculate_processed_data(harmonics, sampledWaveform, true);
                    voltage.set_processed(processed);
                }
            }
            operatingPoint.get_mutable_excitations_per_winding()[windingIndex].set_voltage(voltage);
        }
    }

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_voltage()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];

    bool includeFringing = settings.get_painter_include_fringing();
    bool mirroringDimension = settings.get_painter_mirroring_dimension();

    size_t numberPointsX = settings.get_painter_number_points_x();
    size_t numberPointsY = settings.get_painter_number_points_y();
    auto oldCoilMesherInsideTurnsFactor = settings.get_coil_mesher_inside_turns_factor();
    settings.set_coil_mesher_inside_turns_factor(1.2);
    Field inducedField = CoilMesher::generate_mesh_induced_grid(magnetic, frequency, numberPointsX, numberPointsY, false, false).first;
    settings.set_coil_mesher_inside_turns_factor(oldCoilMesherInsideTurnsFactor);

    auto strayCapacitanceModel = settings.get_stray_capacitance_model();
    StrayCapacitance strayCapacitance(strayCapacitanceModel);
    settings.set_magnetic_field_include_fringing(includeFringing);
    settings.set_magnetic_field_mirroring_dimension(mirroringDimension);

    auto [pixelXDimension, pixelYDimension] = Painter::get_pixel_dimensions(magnetic);

    auto coil = magnetic.get_coil();
    auto turns = coil.get_turns_description().value();
    auto wirePerWinding = coil.get_wires();


    auto capacitanceOutput = strayCapacitance.calculate_capacitance(coil);
    auto electricEnergyAmongTurns = capacitanceOutput.get_electric_energy_among_turns().value();

    std::set<std::pair<size_t, size_t>> turnsCombinations;
    for (size_t pointIndex = 0; pointIndex < inducedField.get_data().size(); ++pointIndex) {
        auto inducedFieldPoint = inducedField.get_data()[pointIndex];
        double fieldValue = 0;
        std::set<std::pair<std::string, std::string>> turnsCombinations;

        for (auto [firstTurnName, aux] : electricEnergyAmongTurns) {
            auto firstTurn = coil.get_turn_by_name(firstTurnName);
            for (auto [secondTurnName, energy] : aux) {
                auto key = std::make_pair(firstTurnName, secondTurnName);

                if (turnsCombinations.contains(key) || turnsCombinations.contains(std::make_pair(secondTurnName, firstTurnName))) {
                    continue;
                }
                turnsCombinations.insert(key);

                auto secondTurn = coil.get_turn_by_name(secondTurnName);
                double pixelArea = Painter::get_pixel_area_between_turns(firstTurn.get_coordinates(), firstTurn.get_dimensions().value(), firstTurn.get_cross_sectional_shape().value(), secondTurn.get_coordinates(), secondTurn.get_dimensions().value(), secondTurn.get_cross_sectional_shape().value(), inducedFieldPoint.get_point(), std::max(pixelXDimension, pixelYDimension));
                if (pixelArea > 0) {
                    double area = StrayCapacitance::calculate_area_between_two_turns(firstTurn, secondTurn);
                    double energyDensity = energy / area;
                    fieldValue += energyDensity * pixelArea;
                }
            }
        }
        inducedField.get_mutable_data()[pointIndex].set_value(fieldValue);
    }

    return inducedField;
}

std::shared_ptr<PainterInterface> Painter::factory(bool useAdvancedPainter, std::filesystem::path filepath, bool addProportionForColorBar, bool showTicks) {
#ifdef ENABLE_MATPLOTPP
    if (useAdvancedPainter) {
        return std::make_shared<AdvancedPainter>(filepath, addProportionForColorBar, showTicks);
    }
    else {
        return std::make_shared<BasicPainter>(filepath);
    }
#else
    (void)useAdvancedPainter;
    (void)addProportionForColorBar;
    (void)showTicks;
    return std::make_shared<BasicPainter>(filepath);
#endif
}

void Painter::paint_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex, std::optional<ComplexField> inputField) {
    _painter->paint_magnetic_field(operatingPoint, magnetic, harmonicIndex, inputField);
}

void Painter::paint_electric_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex, std::optional<Field> inputField, ElectricFieldVisualizationModel model, ColorPalette colorPalette) {
    _painter->paint_electric_field(operatingPoint, magnetic, harmonicIndex, inputField, model, colorPalette);
}

void Painter::paint_wire_losses(Magnetic magnetic, std::optional<Outputs> outputs, std::optional<OperatingPoint> operatingPoint, double temperature) {
    _painter->paint_wire_losses(magnetic, outputs, operatingPoint, temperature);
}

std::string Painter::export_svg() {
    return _painter->export_svg();
}

void Painter::export_png() {
    _painter->export_png();
}

void Painter::paint_core(Magnetic magnetic) {
    _painter->paint_core(magnetic);
}

void Painter::paint_bobbin(Magnetic magnetic) {
    _painter->paint_bobbin(magnetic);
}

void Painter::paint_coil_sections(Magnetic magnetic) {
    _painter->paint_coil_sections(magnetic);
}

void Painter::paint_coil_layers(Magnetic magnetic) {
    _painter->paint_coil_layers(magnetic);
}

void Painter::paint_coil_turns(Magnetic magnetic) {
    _painter->paint_coil_turns(magnetic);
}

void Painter::paint_wire(Wire wire) {
    _painter->paint_wire(wire);
}

void Painter::paint_wire_with_current_density(Wire wire, OperatingPoint operatingPoint, size_t windingIndex) {
    _painter->paint_wire_with_current_density(wire, operatingPoint, windingIndex);
}

void Painter::paint_wire_with_current_density(Wire wire, SignalDescriptor current, double frequency, double temperature) {
    _painter->paint_wire_with_current_density(wire, current, frequency, temperature);
}

void Painter::paint_waveform(Waveform waveform) {
    paint_waveform(waveform.get_data(), waveform.get_time());
}

void Painter::paint_waveform(std::vector<double> data, std::optional<std::vector<double>> time) {
    _painter->paint_waveform(data, time);
}

void Painter::paint_curve(Curve2D curve2D, bool logScale) {
    _painter->paint_curve(curve2D, logScale);
}

void Painter::paint_rectangle(double xCoordinate, double yCoordinate, double xDimension, double yDimension) {
    _painter->paint_rectangle(xCoordinate, yCoordinate, xDimension, yDimension);
}

void Painter::paint_circle(double xCoordinate, double yCoordinate, double radius){
    _painter->paint_circle(xCoordinate, yCoordinate, radius);
}

double Painter::get_pixel_area_between_turns(std::vector<double> firstTurnCoordinates, std::vector<double> firstTurnDimensions, TurnCrossSectionalShape firstTurncrossSectionalShape,
                                             std::vector<double> secondTurnCoordinates, std::vector<double> secondTurnDimensions, TurnCrossSectionalShape secondTurncrossSectionalShape,
                                             std::vector<double> pixelCoordinates, double dimension) {
    return dimension * dimension * get_pixel_proportion_between_turns(firstTurnCoordinates, firstTurnDimensions, firstTurncrossSectionalShape, secondTurnCoordinates, secondTurnDimensions, secondTurncrossSectionalShape, pixelCoordinates, dimension);
}

double Painter::get_pixel_proportion_between_turns(std::vector<double> firstTurnCoordinates, std::vector<double> firstTurnDimensions, TurnCrossSectionalShape firstTurncrossSectionalShape,
                                             std::vector<double> secondTurnCoordinates, std::vector<double> secondTurnDimensions, TurnCrossSectionalShape secondTurncrossSectionalShape,
                                             std::vector<double> pixelCoordinates, double dimension) {
    // auto factor = Defaults().overlappingFactorSurroundingTurns;
    auto x1 = firstTurnCoordinates[0];
    auto y1 = firstTurnCoordinates[1];
    auto x2 = secondTurnCoordinates[0];
    auto y2 = secondTurnCoordinates[1];

    if (y2 == y1 && x2 == x1) {
        return 0;
    }


    double firstTurnMaximumDimension = 0;
    double secondTurnMaximumDimension = 0;

    if (firstTurncrossSectionalShape == TurnCrossSectionalShape::RECTANGULAR) {
        firstTurnMaximumDimension = hypot(firstTurnDimensions[0], firstTurnDimensions[1]);
    }
    else {
        firstTurnMaximumDimension = firstTurnDimensions[0];
    }
    if (secondTurncrossSectionalShape == TurnCrossSectionalShape::RECTANGULAR) {
        secondTurnMaximumDimension = hypot(secondTurnDimensions[0], secondTurnDimensions[1]);
    }
    else {
        secondTurnMaximumDimension = secondTurnDimensions[0];
    }

    double semiAverageDimensionOf12 = (firstTurnMaximumDimension + secondTurnMaximumDimension) / 4;

    auto x0 = pixelCoordinates[0];
    auto y0 = pixelCoordinates[1];

    auto distanceFrom0toLine12 = fabs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / sqrt(pow(y2 - y1, 2) + pow(x2 - x1, 2));
    auto distanceFrom0toCenter1 = hypot(x1 - x0, y1 - y0);
    auto distanceFrom0toCenter2 = hypot(x2 - x0, y2 - y0);
    auto distanceFromCenter1toCenter2 = hypot(x2 - x1, y2 - y1);

    if (distanceFrom0toCenter1 > distanceFromCenter1toCenter2 || distanceFrom0toCenter2 > distanceFromCenter1toCenter2) {
        return 0;   
    }

    double proportion;
    if (distanceFrom0toLine12 - dimension / 2 > semiAverageDimensionOf12) {
        proportion = 0;
    }
    else if (distanceFrom0toLine12 + dimension / 2 < semiAverageDimensionOf12) {
        proportion = 1;
    }
    else {
        proportion = (semiAverageDimensionOf12 - (distanceFrom0toLine12 - dimension / 2)) / dimension;
    }

    proportion *= (semiAverageDimensionOf12 - distanceFrom0toLine12) / semiAverageDimensionOf12;

    return proportion;
}

std::pair<double, double> Painter::get_pixel_dimensions(Magnetic magnetic) {

    size_t numberPointsX = settings.get_painter_number_points_x();
    size_t numberPointsY = settings.get_painter_number_points_y();

    auto columns = magnetic.get_mutable_core().get_columns();
    if (columns.empty()) {
        return {0.0, 0.0};
    }
    double coreColumnHeight = columns[0].get_height();

    auto windingWindow = magnetic.get_mutable_core().get_winding_window();
    if (!windingWindow.get_width()) {
        return {0.0, 0.0};
    }
    double coreWindingWindowWidth = windingWindow.get_width().value();

    double pixelXDimension = coreWindingWindowWidth / numberPointsX;
    double pixelYDimension = coreColumnHeight / numberPointsY;

    return {pixelXDimension, pixelYDimension};
}


} // namespace OpenMagnetics
