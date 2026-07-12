#include "physical_models/ReluctanceNetwork.h"
#include "support/Exceptions.h"
#include <cmath>

namespace OpenMagnetics {

ReluctanceNetwork::ReluctanceNetwork(Core core, double ungappedCoreReluctance, std::vector<AirGapReluctanceOutput> reluctancePerGap) {
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("Core is not processed, cannot build its magnetic circuit");
    }
    if (std::isnan(ungappedCoreReluctance) || ungappedCoreReluctance <= 0) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                                    "Ungapped core reluctance must be a positive number, got " + std::to_string(ungappedCoreReluctance));
    }
    _core = core;
    _mainColumnIndex = core.get_main_column_index();
    auto columns = core.get_columns();

    // Split the lumped ungapped reluctance geometrically across the column branches:
    // relative branch reluctance = path length / area, where return branches carry
    // their yoke share on top of the column length. The split is then scaled so the
    // main-column driving-point combination (R_main + parallel of the returns)
    // reproduces the lumped value exactly.
    std::vector<double> coreBranchReluctances(columns.size(), 0);
    if (columns.size() == 1) {
        coreBranchReluctances[0] = ungappedCoreReluctance;
    }
    else {
        double totalColumnLength = 0;
        for (auto& column : columns) {
            totalColumnLength += column.get_height();
        }
        double yokeSharePerSide = std::max((core.get_effective_length() - totalColumnLength) / 2, 0.0);
        std::vector<double> relativeBranchReluctances(columns.size(), 0);
        for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
            double pathLength = columns[columnIndex].get_height();
            if (columnIndex != _mainColumnIndex) {
                pathLength += 2 * yokeSharePerSide;
            }
            if (columns[columnIndex].get_area() <= 0) {
                throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                                            "Column " + std::to_string(columnIndex) + " has no positive area");
            }
            relativeBranchReluctances[columnIndex] = pathLength / columns[columnIndex].get_area();
        }
        double returnAdmittance = 0;
        for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
            if (columnIndex != _mainColumnIndex) {
                returnAdmittance += 1 / relativeBranchReluctances[columnIndex];
            }
        }
        double relativeCombination = relativeBranchReluctances[_mainColumnIndex] + 1 / returnAdmittance;
        double scale = ungappedCoreReluctance / relativeCombination;
        for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
            coreBranchReluctances[columnIndex] = scale * relativeBranchReluctances[columnIndex];
        }
    }

    // Add each gap's reluctance to the branch of the column it sits in. The per-gap
    // reluctances arrive aligned with the core's gapping order.
    auto gapping = _core.get_functional_description().get_gapping();
    if (gapping.size() != reluctancePerGap.size()) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                                    "Per-gap reluctances (" + std::to_string(reluctancePerGap.size()) +
                                        ") do not match the core gapping (" + std::to_string(gapping.size()) + " gaps)");
    }
    for (auto& gap : gapping) {
        if (!gap.get_coordinates()) {
            _core.process_gap();
            gapping = _core.get_functional_description().get_gapping();
            break;
        }
    }
    _columnReluctances = coreBranchReluctances;
    for (size_t gapIndex = 0; gapIndex < gapping.size(); ++gapIndex) {
        auto columnIndex = _core.find_closest_column_index_by_coordinates(gapping[gapIndex].get_coordinates().value());
        _columnReluctances[static_cast<size_t>(columnIndex)] += reluctancePerGap[gapIndex].get_reluctance();
    }
}

std::vector<size_t> ReluctanceNetwork::resolve_winding_column_indexes(Magnetic magnetic) {
    auto core = magnetic.get_core();
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("Core is not processed, cannot resolve winding columns");
    }
    auto windingWindows = core.get_winding_windows();
    size_t mainColumnIndex = core.get_main_column_index();
    std::vector<size_t> columnIndexPerWinding;
    for (auto& winding : magnetic.get_coil().get_functional_description()) {
        auto windowIndex = winding.get_winding_window();
        if (!windowIndex) {
            columnIndexPerWinding.push_back(mainColumnIndex);
            continue;
        }
        if (windowIndex.value() < 0 || static_cast<size_t>(windowIndex.value()) >= windingWindows.size()) {
            throw InvalidInputException(ErrorCode::INVALID_COIL_CONFIGURATION,
                                        "Winding " + winding.get_name() + " references winding window " +
                                            std::to_string(windowIndex.value()) + " but the core has " +
                                            std::to_string(windingWindows.size()) + " winding windows");
        }
        columnIndexPerWinding.push_back(core.get_winding_window_column_index(static_cast<size_t>(windowIndex.value())));
    }
    return columnIndexPerWinding;
}

bool ReluctanceNetwork::has_non_main_placement(Magnetic magnetic) {
    bool anyPlacement = false;
    for (auto& winding : magnetic.get_coil().get_functional_description()) {
        if (winding.get_winding_window()) {
            anyPlacement = true;
            break;
        }
    }
    if (!anyPlacement) {
        return false;
    }
    auto columnIndexPerWinding = resolve_winding_column_indexes(magnetic);
    size_t mainColumnIndex = magnetic.get_core().get_main_column_index();
    for (auto columnIndex : columnIndexPerWinding) {
        if (columnIndex != mainColumnIndex) {
            return true;
        }
    }
    return false;
}

std::vector<std::vector<double>> ReluctanceNetwork::calculate_magnetizing_inductance_matrix(Magnetic magnetic) const {
    auto columnIndexPerWinding = resolve_winding_column_indexes(magnetic);
    auto windings = magnetic.get_coil().get_functional_description();

    std::vector<double> columnAdmittances;
    double totalAdmittance = 0;
    for (auto reluctance : _columnReluctances) {
        columnAdmittances.push_back(1 / reluctance);
        totalAdmittance += 1 / reluctance;
    }

    std::vector<std::vector<double>> inductanceMatrix(windings.size(), std::vector<double>(windings.size(), 0));
    for (size_t i = 0; i < windings.size(); ++i) {
        for (size_t j = 0; j < windings.size(); ++j) {
            double numberTurnsI = windings[i].get_number_turns();
            double numberTurnsJ = windings[j].get_number_turns();
            double sameColumnTerm = columnIndexPerWinding[i] == columnIndexPerWinding[j] ? columnAdmittances[columnIndexPerWinding[i]] : 0;
            inductanceMatrix[i][j] = numberTurnsI * numberTurnsJ *
                                     (sameColumnTerm - columnAdmittances[columnIndexPerWinding[i]] * columnAdmittances[columnIndexPerWinding[j]] / totalAdmittance);
        }
    }
    return inductanceMatrix;
}

std::vector<double> ReluctanceNetwork::calculate_column_magnetic_fluxes(Magnetic magnetic, const std::vector<double>& windingCurrents) const {
    auto columnIndexPerWinding = resolve_winding_column_indexes(magnetic);
    auto windings = magnetic.get_coil().get_functional_description();
    if (windingCurrents.size() != windings.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
                                    "One current per winding is needed: got " + std::to_string(windingCurrents.size()) +
                                        " currents for " + std::to_string(windings.size()) + " windings");
    }

    std::vector<double> columnMagnetomotiveForces(_columnReluctances.size(), 0);
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        columnMagnetomotiveForces[columnIndexPerWinding[windingIndex]] +=
            static_cast<double>(windings[windingIndex].get_number_turns()) * windingCurrents[windingIndex];
    }

    double totalAdmittance = 0;
    double drivenAdmittance = 0;
    for (size_t columnIndex = 0; columnIndex < _columnReluctances.size(); ++columnIndex) {
        totalAdmittance += 1 / _columnReluctances[columnIndex];
        drivenAdmittance += columnMagnetomotiveForces[columnIndex] / _columnReluctances[columnIndex];
    }
    double yokePotentialDifference = drivenAdmittance / totalAdmittance;

    std::vector<double> columnFluxes(_columnReluctances.size(), 0);
    for (size_t columnIndex = 0; columnIndex < _columnReluctances.size(); ++columnIndex) {
        columnFluxes[columnIndex] = (columnMagnetomotiveForces[columnIndex] - yokePotentialDifference) / _columnReluctances[columnIndex];
    }
    return columnFluxes;
}

std::vector<double> ReluctanceNetwork::calculate_column_magnetic_flux_densities(Magnetic magnetic, const std::vector<double>& windingCurrents) const {
    auto columnFluxes = calculate_column_magnetic_fluxes(magnetic, windingCurrents);
    auto columns = _core.get_columns();
    std::vector<double> columnFluxDensities(columnFluxes.size(), 0);
    for (size_t columnIndex = 0; columnIndex < columnFluxes.size(); ++columnIndex) {
        columnFluxDensities[columnIndex] = columnFluxes[columnIndex] / columns[columnIndex].get_area();
    }
    return columnFluxDensities;
}

} // namespace OpenMagnetics
