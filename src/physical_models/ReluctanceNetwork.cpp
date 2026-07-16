#include "physical_models/ReluctanceNetwork.h"
#include "support/Exceptions.h"
#include <algorithm>
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

    // Leg-chain mesh network: the legs ordered by x, one yoke segment (top plus
    // bottom plate) per window between adjacent legs, n-1 mesh loops. The FD
    // validation (ABT #227) showed the yoke segments are NOT negligible: even at
    // high permeability the long span to the FAR leg rivals the leg reluctances, so
    // collapsing the yokes to perfect conductors (plain parallel branches)
    // overweights the far return path. Everything core is scaled so the main-column
    // driving-point reluctance reproduces the lumped IEC value exactly.
    _orderedColumnIndexes.resize(columns.size());
    for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
        _orderedColumnIndexes[columnIndex] = columnIndex;
    }
    std::sort(_orderedColumnIndexes.begin(), _orderedColumnIndexes.end(), [&](size_t left, size_t right) {
        return columns[left].get_coordinates()[0] < columns[right].get_coordinates()[0];
    });

    std::vector<double> coreBranchReluctances(columns.size(), 0);
    _yokeSegmentReluctances.assign(columns.size() > 1 ? columns.size() - 1 : 0, 0);
    if (columns.size() == 1) {
        coreBranchReluctances[0] = ungappedCoreReluctance;
    }
    else {
        std::vector<double> relativeLegReluctances(columns.size(), 0);
        for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
            if (columns[columnIndex].get_area() <= 0) {
                throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                                            "Column " + std::to_string(columnIndex) + " has no positive area");
            }
            relativeLegReluctances[columnIndex] = columns[columnIndex].get_height() / columns[columnIndex].get_area();
        }

        auto windingWindows = core.get_winding_windows();
        if (!windingWindows[0].get_height()) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                                        "Multi-column network needs a rectangular winding window with a height");
        }
        double yokeThickness = (core.get_processed_description().value().get_height() - windingWindows[0].get_height().value()) / 2;
        double yokeArea = yokeThickness * columns[_mainColumnIndex].get_depth();
        if (yokeArea <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core yoke has no positive area");
        }
        std::vector<double> relativeYokeReluctances(columns.size() - 1, 0);
        for (size_t segmentIndex = 0; segmentIndex + 1 < columns.size(); ++segmentIndex) {
            double span = std::abs(columns[_orderedColumnIndexes[segmentIndex + 1]].get_coordinates()[0] -
                                   columns[_orderedColumnIndexes[segmentIndex]].get_coordinates()[0]);
            // Top and bottom plates in series along the loop.
            relativeYokeReluctances[segmentIndex] = 2 * span / yokeArea;
        }

        // Calibrate: the relative network's main-column driving-point reluctance must
        // equal the lumped ungapped value.
        std::vector<double> savedColumnReluctances = _columnReluctances;
        _columnReluctances = relativeLegReluctances;
        _yokeSegmentReluctances = relativeYokeReluctances;
        std::vector<double> unitMagnetomotiveForce(columns.size(), 0);
        unitMagnetomotiveForce[_mainColumnIndex] = 1.0;
        double relativeDrivingReluctance = 1.0 / solve_column_fluxes_for_magnetomotive_forces(unitMagnetomotiveForce)[_mainColumnIndex];
        double scale = ungappedCoreReluctance / relativeDrivingReluctance;
        for (size_t columnIndex = 0; columnIndex < columns.size(); ++columnIndex) {
            coreBranchReluctances[columnIndex] = scale * relativeLegReluctances[columnIndex];
        }
        for (auto& yokeSegmentReluctance : _yokeSegmentReluctances) {
            yokeSegmentReluctance *= scale;
        }
        _columnReluctances = savedColumnReluctances;
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

// Placement precedence mirrors the winder: winding-level windingWindow wins,
// else the winding's wound SECTIONS (via group, already resolved onto each
// section) name the window. Returns nullopt when no placement is recorded.
static std::optional<int64_t> resolve_window_from_sections(const OpenMagnetics::Coil& coil, const std::string& windingName) {
    // Materialize the copy: the generated getter returns the optional BY VALUE,
    // so iterating `get_sections_description().value()` directly dangles.
    auto sectionsDescription = coil.get_sections_description();
    if (!sectionsDescription) {
        return std::nullopt;
    }
    for (const auto& section : sectionsDescription.value()) {
        if (section.get_type() != ElectricalType::CONDUCTION || !section.get_winding_window()) {
            continue;
        }
        for (const auto& partialWinding : section.get_partial_windings()) {
            if (partialWinding.get_winding() == windingName) {
                return section.get_winding_window();
            }
        }
    }
    return std::nullopt;
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
            // The placement may live on the SECTIONS instead of the winding
            // (e.g. hand-authored MAS files, or intents carried by sections
            // only). Without this fallback such a lateral-placed coil got the
            // ideal rank-1 coupling instead of the real flux divider —
            // a silently wrong inductance matrix.
            windowIndex = resolve_window_from_sections(magnetic.get_coil(), winding.get_name());
        }
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
        if (winding.get_winding_window() || resolve_window_from_sections(magnetic.get_coil(), winding.get_name())) {
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

std::vector<double> ReluctanceNetwork::solve_column_fluxes_for_magnetomotive_forces(const std::vector<double>& columnMagnetomotiveForces) const {
    size_t numberColumns = _columnReluctances.size();
    if (numberColumns == 1) {
        return {columnMagnetomotiveForces[0] / _columnReluctances[0]};
    }

    // Mesh analysis over the leg chain (legs x-ordered, all oriented the same way):
    // loop k circulates in the window between ordered legs k and k+1 through the yoke
    // segment k. A[k][k] = R_k + R_{k+1} + Y_k, A[k][k±1] = -R_shared,
    // rhs[k] = MMF_k − MMF_{k+1}; leg flux = φ_k − φ_{k−1}.
    size_t numberLoops = numberColumns - 1;
    std::vector<double> orderedLegReluctances(numberColumns);
    std::vector<double> orderedMagnetomotiveForces(numberColumns);
    for (size_t position = 0; position < numberColumns; ++position) {
        orderedLegReluctances[position] = _columnReluctances[_orderedColumnIndexes[position]];
        orderedMagnetomotiveForces[position] = columnMagnetomotiveForces[_orderedColumnIndexes[position]];
    }

    std::vector<std::vector<double>> meshMatrix(numberLoops, std::vector<double>(numberLoops, 0));
    std::vector<double> rightHandSide(numberLoops, 0);
    for (size_t loopIndex = 0; loopIndex < numberLoops; ++loopIndex) {
        meshMatrix[loopIndex][loopIndex] = orderedLegReluctances[loopIndex] + orderedLegReluctances[loopIndex + 1] +
                                           _yokeSegmentReluctances[loopIndex];
        if (loopIndex > 0) {
            meshMatrix[loopIndex][loopIndex - 1] = -orderedLegReluctances[loopIndex];
        }
        if (loopIndex + 1 < numberLoops) {
            meshMatrix[loopIndex][loopIndex + 1] = -orderedLegReluctances[loopIndex + 1];
        }
        rightHandSide[loopIndex] = orderedMagnetomotiveForces[loopIndex] - orderedMagnetomotiveForces[loopIndex + 1];
    }

    // Gaussian elimination with partial pivoting (the system is tiny: n-1 loops).
    std::vector<double> loopFluxes(numberLoops, 0);
    for (size_t pivotIndex = 0; pivotIndex < numberLoops; ++pivotIndex) {
        size_t bestRow = pivotIndex;
        for (size_t row = pivotIndex + 1; row < numberLoops; ++row) {
            if (std::abs(meshMatrix[row][pivotIndex]) > std::abs(meshMatrix[bestRow][pivotIndex])) {
                bestRow = row;
            }
        }
        std::swap(meshMatrix[pivotIndex], meshMatrix[bestRow]);
        std::swap(rightHandSide[pivotIndex], rightHandSide[bestRow]);
        if (meshMatrix[pivotIndex][pivotIndex] == 0) {
            throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Singular reluctance network");
        }
        for (size_t row = pivotIndex + 1; row < numberLoops; ++row) {
            double factor = meshMatrix[row][pivotIndex] / meshMatrix[pivotIndex][pivotIndex];
            for (size_t column = pivotIndex; column < numberLoops; ++column) {
                meshMatrix[row][column] -= factor * meshMatrix[pivotIndex][column];
            }
            rightHandSide[row] -= factor * rightHandSide[pivotIndex];
        }
    }
    for (size_t rowPlusOne = numberLoops; rowPlusOne > 0; --rowPlusOne) {
        size_t row = rowPlusOne - 1;
        double accumulated = rightHandSide[row];
        for (size_t column = row + 1; column < numberLoops; ++column) {
            accumulated -= meshMatrix[row][column] * loopFluxes[column];
        }
        loopFluxes[row] = accumulated / meshMatrix[row][row];
    }

    std::vector<double> columnFluxes(numberColumns, 0);
    for (size_t position = 0; position < numberColumns; ++position) {
        double flux = 0;
        if (position < numberLoops) {
            flux += loopFluxes[position];
        }
        if (position > 0) {
            flux -= loopFluxes[position - 1];
        }
        columnFluxes[_orderedColumnIndexes[position]] = flux;
    }
    return columnFluxes;
}

std::vector<std::vector<double>> ReluctanceNetwork::calculate_magnetizing_inductance_matrix(Magnetic magnetic) const {
    auto columnIndexPerWinding = resolve_winding_column_indexes(magnetic);
    auto windings = magnetic.get_coil().get_functional_description();

    std::vector<std::vector<double>> inductanceMatrix(windings.size(), std::vector<double>(windings.size(), 0));
    for (size_t j = 0; j < windings.size(); ++j) {
        std::vector<double> columnMagnetomotiveForces(_columnReluctances.size(), 0);
        columnMagnetomotiveForces[columnIndexPerWinding[j]] = static_cast<double>(windings[j].get_number_turns());
        auto columnFluxes = solve_column_fluxes_for_magnetomotive_forces(columnMagnetomotiveForces);
        for (size_t i = 0; i < windings.size(); ++i) {
            inductanceMatrix[i][j] = static_cast<double>(windings[i].get_number_turns()) * columnFluxes[columnIndexPerWinding[i]];
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
    return solve_column_fluxes_for_magnetomotive_forces(columnMagnetomotiveForces);
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
