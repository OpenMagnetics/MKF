#include "constructive_models/Core.h"
#include "constructive_models/Wire.h"
#include "physical_models/Resistivity.h"
#include "physical_models/Reluctance.h"
#include "physical_models/InitialPermeability.h"
#include "spline.h"

#include "json.hpp"
#include "support/Utils.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <streambuf>
#include <vector>

using json = nlohmann::json;

namespace OpenMagnetics {


Core::Core(const json& j, bool includeMaterialData, bool includeProcessedDescription, bool includeGeometricalDescription) {
    _includeMaterialData = includeMaterialData;
    from_json(j, *this);
    
    if (includeProcessedDescription) {
        process_data();
        process_gap();
    }

    if (!get_geometrical_description() && includeGeometricalDescription) {
        auto geometricalDescription = create_geometrical_description();

        set_geometrical_description(geometricalDescription);
    }
}

Core::Core(const MagneticCore core) {
    set_functional_description(core.get_functional_description());

    if (core.get_geometrical_description()) {
        set_geometrical_description(core.get_geometrical_description());
    }
    if (core.get_processed_description()) {
        set_processed_description(core.get_processed_description());
    }
    if (core.get_distributors_info()) {
        set_distributors_info(core.get_distributors_info());
    }
    if (core.get_manufacturer_info()) {
        set_manufacturer_info(core.get_manufacturer_info());
    }
    if (core.get_name()) {
        set_name(core.get_name());
    }
}

Core::Core(const CoreShape shape, std::optional<CoreMaterial> material) {
    get_mutable_functional_description().set_gapping(std::vector<CoreGap>({}));
    get_mutable_functional_description().set_number_stacks(1);
    get_mutable_functional_description().set_shape(shape);
    if (material) {
        get_mutable_functional_description().set_material(material.value());
    }
    else {
        get_mutable_functional_description().set_material("Dummy");
    }
    if (shape.get_magnetic_circuit() == MagneticCircuit::OPEN) {
        get_mutable_functional_description().set_type(CoreType::TWO_PIECE_SET);
    }
    else {
        get_mutable_functional_description().set_type(CoreType::TOROIDAL);
    }

    if (material) {
        set_name(shape.get_name().value() + " " + material.value().get_name());
    }
    else {
        set_name(shape.get_name().value());
    }
}

double Core::get_depth() {
    if (!get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }
    return get_processed_description()->get_depth();
}

double Core::get_height() {
    if (!get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }
    return get_processed_description()->get_height();
}

double Core::get_width() {
    if (!get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }
    return get_processed_description()->get_width();
}

double Core::get_mass() {
    if (get_shape_family() == CoreShapeFamily::T) {
        auto dimensions = flatten_dimensions(resolve_shape().get_dimensions().value());
        double volume = std::numbers::pi * (pow(dimensions["A"] / 2, 2) - pow(dimensions["B"] / 2, 2)) * dimensions["C"];
        return volume * get_density();
    }
    else{
        throw std::runtime_error("get_mass only implemented for toroidal cores for now");
    }
}

double Core::get_effective_length() {
    if (!get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }
    return get_processed_description()->get_effective_parameters().get_effective_length();
}

double Core::get_effective_area() {
    if (!get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }
    return get_processed_description()->get_effective_parameters().get_effective_area();
}

double Core::get_minimum_area() {
    if (!get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }
    return get_processed_description()->get_effective_parameters().get_minimum_area();
}

double Core::get_effective_volume() {
    if (!get_processed_description()) {
        throw std::runtime_error("Core is not processed");
    }
    return get_processed_description()->get_effective_parameters().get_effective_volume();
}

std::string Core::get_reference() {
    if (get_manufacturer_info()) {
        if (get_manufacturer_info()->get_reference()) {
            return get_manufacturer_info()->get_reference().value();
        }
    }
    return "";
}

double interp(std::vector<std::pair<double, double>> data, double temperature) {
    double result;

    if (data.size() == 0) {
        throw std::runtime_error("Data cannot be empty");
    }
    else if (data.size() == 1) {
        result = data[0].second;
    }
    else if (data.size() == 2) {
        result = data[0].second -
        (data[0].first - temperature) *
            (data[0].second - data[1].second) /
            (data[0].first - data[1].first);
    }
    else {
        std::sort(data.begin(), data.end(), [](const std::pair<double, double>& b1, const std::pair<double, double>& b2) {
                return b1.first < b2.first;
        });
        size_t n = data.size();
        std::vector<double> x, y;

        for (size_t i = 0; i < n; i++) {
            if (x.size() == 0 || data[i].first != x.back()) {
                x.push_back(data[i].first);
                y.push_back(data[i].second);
            }
        }
        tk::spline interp(x, y, tk::spline::cspline_hermite, true);
        result = interp(temperature);
    }
    return result;
}

std::optional<std::vector<CoreGeometricalDescriptionElement>> Core::create_geometrical_description() {
    std::vector<CoreGeometricalDescriptionElement> geometricalDescription;
    auto numberStacks = *(get_functional_description().get_number_stacks());
    auto gapping = get_functional_description().get_gapping();

    auto corePiece = CorePiece::factory(std::get<CoreShape>(get_functional_description().get_shape()));
    auto corePieceHeight = corePiece->get_height();
    auto corePieceDepth = corePiece->get_depth();

    CoreGeometricalDescriptionElement spacer;
    std::vector<Machining> machining;
    CoreGeometricalDescriptionElement piece;
    double currentDepth = roundFloat((-corePieceDepth * (numberStacks - 1)) / 2);
    double spacerThickness = 0;

    for (auto& gap : gapping) {
        if (gap.get_type() == GapType::ADDITIVE) {
            spacerThickness = gap.get_length();
        }
        else if (gap.get_type() == GapType::SUBTRACTIVE) {
            Machining aux;
            aux.set_length(gap.get_length());
            if (!gap.get_coordinates()) {
                aux.set_coordinates({0, 0, 0});
            }
            else {
                aux.set_coordinates(gap.get_coordinates().value());
            }
            // else {  // Gap coordinates are centered in the gap, while machining coordinates start at the base of the
            // central column surface, therefore, we must correct them
            //     aux["coordinates"][1] = aux["coordinates"].at(1).get<double>() - gap.get_length() / 2;
            // }
            machining.push_back(aux);
        }
    }

    piece.set_material(resolve_material().get_name());
    piece.set_shape(std::get<CoreShape>(get_functional_description().get_shape()));
    auto topPiece = piece;
    auto bottomPiece = piece;
    switch (get_functional_description().get_type()) {
        case CoreType::TOROIDAL:
            piece.set_type(CoreGeometricalDescriptionElementType::TOROIDAL);
            for (auto i = 0; i < numberStacks; ++i) {
                std::vector<double> coordinates = {0, 0, currentDepth};
                piece.set_coordinates(coordinates);
                piece.set_rotation(std::vector<double>({std::numbers::pi / 2, std::numbers::pi / 2, 0}));

                geometricalDescription.push_back(CoreGeometricalDescriptionElement(piece));

                currentDepth = roundFloat(currentDepth + corePieceDepth);
            }
            break;
        case CoreType::CLOSED_SHAPE:
            piece.set_type(CoreGeometricalDescriptionElementType::CLOSED);
            for (auto i = 0; i < numberStacks; ++i) {
                double currentHeight = roundFloat(corePieceHeight);
                std::vector<double> coordinates = {0, currentHeight, currentDepth};
                piece.set_coordinates(coordinates);
                piece.set_rotation(std::vector<double>({0, 0, 0}));
                piece.set_machining(std::nullopt);
                geometricalDescription.push_back(CoreGeometricalDescriptionElement(piece));
                currentDepth = roundFloat(currentDepth + corePieceDepth);
            }
            break;
        case CoreType::TWO_PIECE_SET:
            topPiece.set_type(CoreGeometricalDescriptionElementType::HALF_SET);
            bottomPiece.set_type(CoreGeometricalDescriptionElementType::HALF_SET);
            for (auto i = 0; i < numberStacks; ++i) {
                double currentHeight = roundFloat(spacerThickness / 2);
                std::vector<Machining> topHalfMachining;
                std::vector<double> coordinates = {0, currentHeight, currentDepth};
                topPiece.set_coordinates(coordinates);
                topPiece.set_rotation(std::vector<double>({std::numbers::pi, std::numbers::pi, 0}));
                for (auto& operating : machining) {
                    if (operating.get_coordinates()[1] >= 0 &&
                        operating.get_coordinates()[1] < operating.get_length() / 2) {
                        Machining brokenDownOperation;
                        brokenDownOperation.set_coordinates(operating.get_coordinates());
                        brokenDownOperation.set_length(operating.get_length() / 2 + operating.get_coordinates()[1]);
                        brokenDownOperation.get_mutable_coordinates()[1] = brokenDownOperation.get_length() / 2;
                        topHalfMachining.push_back(brokenDownOperation);
                    }
                    else if (operating.get_coordinates()[1] < 0 &&
                             roundFloat(operating.get_coordinates()[1] + operating.get_length() / 2, 9) > 0) {
                        Machining brokenDownOperation;
                        brokenDownOperation.set_coordinates(operating.get_coordinates());
                        brokenDownOperation.set_length(operating.get_length() / 2 + operating.get_coordinates()[1]);
                        brokenDownOperation.get_mutable_coordinates()[1] = brokenDownOperation.get_length() / 2;
                        topHalfMachining.push_back(brokenDownOperation);
                    }
                    else if (operating.get_coordinates()[1] > 0) {
                        topHalfMachining.push_back(operating);
                    }
                }

                if (topHalfMachining.size() > 0) {
                    topPiece.set_machining(topHalfMachining);
                }
                geometricalDescription.push_back(topPiece);

                std::vector<Machining> bottomHalfMachining;

                for (auto& operating : machining) {
                    if (operating.get_coordinates()[1] <= 0 &&
                        (-operating.get_coordinates()[1] < operating.get_length() / 2)) {
                        Machining brokenDownOperation;
                        brokenDownOperation.set_coordinates(operating.get_coordinates());
                        brokenDownOperation.set_length(operating.get_length() / 2 - operating.get_coordinates()[1]);
                        brokenDownOperation.get_mutable_coordinates()[1] = -brokenDownOperation.get_length() / 2;
                        bottomHalfMachining.push_back(brokenDownOperation);
                    }
                    else if (operating.get_coordinates()[1] > 0 &&
                             roundFloat(operating.get_coordinates()[1] - operating.get_length() / 2, 9) < 0) {
                        Machining brokenDownOperation;
                        brokenDownOperation.set_coordinates(operating.get_coordinates());
                        brokenDownOperation.set_length(operating.get_length() / 2 - operating.get_coordinates()[1]);
                        brokenDownOperation.get_mutable_coordinates()[1] = -brokenDownOperation.get_length() / 2;
                        bottomHalfMachining.push_back(brokenDownOperation);
                    }
                    else if (operating.get_coordinates()[1] < 0) {
                        bottomHalfMachining.push_back(operating);
                    }
                }

                if (std::get<CoreShape>(get_functional_description().get_shape()).get_family() ==
                        CoreShapeFamily::UR ||
                    std::get<CoreShape>(get_functional_description().get_shape()).get_family() ==
                        CoreShapeFamily::U ||
                    std::get<CoreShape>(get_functional_description().get_shape()).get_family() ==
                        CoreShapeFamily::C) {
                    bottomPiece.set_rotation(std::vector<double>({0, std::numbers::pi, 0}));
                }
                else {
                    bottomPiece.set_rotation(std::vector<double>({0, 0, 0}));
                }

                if (bottomHalfMachining.size() > 0) {
                    bottomPiece.set_machining(bottomHalfMachining);
                }
                currentHeight = -currentHeight;
                coordinates = {0, currentHeight, currentDepth};
                bottomPiece.set_coordinates(coordinates);
                geometricalDescription.push_back(bottomPiece);

                currentDepth = roundFloat(currentDepth + corePieceDepth);
            }

            if (spacerThickness > 0) {
                for (auto& column : corePiece->get_columns()) {
                    auto shape_data = std::get<CoreShape>(get_functional_description().get_shape());
                    if (column.get_type() == ColumnType::LATERAL) {
                        spacer.set_type(CoreGeometricalDescriptionElementType::SPACER);
                        spacer.set_material("plastic");
                        // We cannot use directly column.get_width()
                        auto dimensions = flatten_dimensions(shape_data.get_dimensions().value());
                        double windingWindowWidth;
                        if (dimensions.find("E") == dimensions.end() ||
                            (roundFloat(dimensions["E"]) == 0)) {
                            if (dimensions.find("F") == dimensions.end() ||
                                (roundFloat(dimensions["F"]) == 0)) {
                                windingWindowWidth = dimensions["A"] -
                                                     dimensions["C"] -
                                                     dimensions["H"];
                            }
                            else {
                                windingWindowWidth = dimensions["A"] -
                                                     dimensions["F"] -
                                                     dimensions["H"];
                            }
                        }
                        else {
                            windingWindowWidth = dimensions["E"];
                        }
                        double minimum_column_width;
                        double minimum_column_depth;
                        if ((shape_data.get_family() == CoreShapeFamily::EP ||
                             shape_data.get_family() == CoreShapeFamily::EPX) &&
                            corePiece->get_columns().size() == 2) {
                            minimum_column_width = dimensions["A"];
                        }
                        else if (shape_data.get_family() == CoreShapeFamily::U ||
                                 shape_data.get_family() == CoreShapeFamily::UR ||
                                 shape_data.get_family() == CoreShapeFamily::C) {
                            if (dimensions.find("H") == dimensions.end() ||
                                (roundFloat(dimensions["H"]) == 0)) {
                                minimum_column_width = (dimensions["A"] - windingWindowWidth) / 2;
                            }
                            else {
                                minimum_column_width = dimensions["H"];
                            }
                        }
                        else {
                            minimum_column_width = (dimensions["A"] - windingWindowWidth) / 2;
                        }

                        if ((shape_data.get_family() == CoreShapeFamily::EP ||
                             shape_data.get_family() == CoreShapeFamily::EPX) &&
                            corePiece->get_columns().size() == 2) {
                            minimum_column_depth = column.get_depth();
                        }
                        else if (shape_data.get_family() == CoreShapeFamily::P ||
                                 shape_data.get_family() == CoreShapeFamily::PM) {
                            minimum_column_depth = dimensions["F"];
                        }
                        else if (shape_data.get_family() == CoreShapeFamily::RM) {
                            if (dimensions.find("J") != dimensions.end() &&
                                (roundFloat(dimensions["J"]) != 0)) {
                                minimum_column_depth =
                                    sqrt(2) * dimensions["J"] - dimensions["A"];
                            }
                            else if (dimensions.find("H") != dimensions.end() &&
                                     (roundFloat(dimensions["H"]) != 0)) {
                                minimum_column_depth = dimensions["H"];
                            }
                            else {
                                minimum_column_depth = dimensions["F"];
                            }
                        }
                        else {
                            minimum_column_depth =
                                std::min(dimensions["C"], column.get_depth()) * numberStacks;
                        }
                        minimum_column_width *= (1 + constants.spacerProtudingPercentage);
                        minimum_column_depth *= (1 + constants.spacerProtudingPercentage);
                        double protuding_width = minimum_column_width * constants.spacerProtudingPercentage;
                        double protuding_depth = minimum_column_depth * constants.spacerProtudingPercentage;
                        spacer.set_dimensions(std::vector<double>({minimum_column_width, spacerThickness, minimum_column_depth}));
                        spacer.set_rotation(std::vector<double>({0, 0, 0}));
                        if (column.get_coordinates()[0] == 0) {
                            spacer.set_coordinates({0, column.get_coordinates()[1],
                                                         -dimensions["C"] / 2 +
                                                             minimum_column_depth / 2 - protuding_depth});
                        }
                        else if (column.get_coordinates()[0] < 0) {
                            if (shape_data.get_family() == CoreShapeFamily::U ||
                                shape_data.get_family() == CoreShapeFamily::UR ||
                                shape_data.get_family() == CoreShapeFamily::C) {
                                spacer.set_coordinates(std::vector<double>({column.get_coordinates()[0] - column.get_width() / 2 +
                                                                 minimum_column_width / 2 - protuding_width,
                                                             column.get_coordinates()[1], column.get_coordinates()[2]}));
                            }
                            else {
                                spacer.set_coordinates(std::vector<double>({-dimensions["A"] / 2 +
                                                                 minimum_column_width / 2 - protuding_width,
                                                             column.get_coordinates()[1], column.get_coordinates()[2]}));
                            }
                        }
                        else {
                            if (shape_data.get_family() == CoreShapeFamily::U ||
                                shape_data.get_family() == CoreShapeFamily::UR ||
                                shape_data.get_family() == CoreShapeFamily::C) {
                                spacer.set_coordinates(std::vector<double>({column.get_coordinates()[0] + column.get_width() / 2 -
                                                                 minimum_column_width / 2 + protuding_width,
                                                             column.get_coordinates()[1], column.get_coordinates()[2]}));
                            }
                            else {
                                spacer.set_coordinates(std::vector<double>({dimensions["A"] / 2 -
                                                                 minimum_column_width / 2 + protuding_width,
                                                             column.get_coordinates()[1], column.get_coordinates()[2]}));
                            }
                        }
                        geometricalDescription.push_back(spacer);
                    }
                }
            }
            break;
        case CoreType::PIECE_AND_PLATE:
            // TODO add for toroPIECE_AND_PLATE
            break;
        default:
            throw std::runtime_error(
                "Unknown type of core, options are {TOROIDAL, TWO_PIECE_SET, PIECE_AND_PLATE, CLOSED_SHAPE}");
    }

    return geometricalDescription;
}

std::vector<ColumnElement> Core::find_columns_by_type(ColumnType columnType) {
    std::vector<ColumnElement> foundColumns;
    auto processedDescription = get_processed_description().value();
    for (auto& column : processedDescription.get_columns()) {
        if (column.get_type() == columnType) {
            foundColumns.push_back(column);
        }
    }
    return foundColumns;
}

int Core::find_closest_column_index_by_coordinates(std::vector<double> coordinates) {
    double closestDistance = std::numeric_limits<double>::infinity();
    int closestColumnIndex = -1;
    auto processedDescription = get_processed_description().value();
    auto columns = processedDescription.get_columns();
    for (size_t index = 0; index < columns.size(); ++index) {
        double distance = 0;
        auto columnCoordinates = columns[index].get_coordinates();
        for (size_t i = 0; i < columnCoordinates.size(); ++i) {
            if (i != 1) { // We don't care about how high in the column the gap is, just about its projection, with are
                          // axis X and Z
                distance += fabs(columnCoordinates[i] - coordinates[i]);
            }
        }
        if (distance < closestDistance) {
            closestColumnIndex = index;
            closestDistance = distance;
        }
    }
    return closestColumnIndex;
}

int Core::find_exact_column_index_by_coordinates(std::vector<double> coordinates) {
    auto processedDescription = get_processed_description().value();
    auto columns = processedDescription.get_columns();
    for (size_t index = 0; index < columns.size(); ++index) {
        double distance = 0;
        auto columnCoordinates = columns[index].get_coordinates();

        double maxCoordinate = 1e-6; // Toavoid dividing by 0;
        for (size_t i = 0; i < columnCoordinates.size(); ++i) {
            if (i != 1) { // We don't care about how high in the column the gap is, just about its projection, with are
                          // axis X and Z
                maxCoordinate = std::max(fabs(columnCoordinates[i]), maxCoordinate);
                distance += fabs(columnCoordinates[i] - coordinates[i]);
            }
        }
        double error = distance / maxCoordinate;
        if (error < 0.01) {
            return index;
        }
    }
    return -1;
}

ColumnElement Core::find_closest_column_by_coordinates(std::vector<double> coordinates) {
    double closestDistance = std::numeric_limits<double>::infinity();
    ColumnElement closestColumn;
    auto processedDescription = get_processed_description().value();
    auto columns = processedDescription.get_columns();
    for (auto& column : columns) {
        double distance = 0;
        auto columnCoordinates = column.get_coordinates();
        for (size_t i = 0; i < columnCoordinates.size(); ++i) {
            distance += fabs(columnCoordinates[i] - coordinates[i]);
        }
        if (distance < closestDistance) {
            closestColumn = column;
            closestDistance = distance;
        }
    }
    return closestColumn;
}

std::vector<CoreGap> Core::find_gaps_by_type(GapType gappingType) {
    std::vector<CoreGap> foundGaps;
    for (auto& gap : get_functional_description().get_gapping()) {
        if (gap.get_type() == gappingType) {
            foundGaps.push_back(gap);
        }
    }
    return foundGaps;
}

std::vector<CoreGap> Core::find_gaps_by_column(ColumnElement column) {
    std::vector<CoreGap> foundGaps;
    double columnLeftLimit = column.get_coordinates()[0] - column.get_width() / 2;
    double columnRightLimit = column.get_coordinates()[0] + column.get_width() / 2;
    double columnFrontLimit = column.get_coordinates()[2] + column.get_depth() / 2;
    double columnBackLimit = column.get_coordinates()[2] - column.get_depth() / 2;

    bool gapIsComplete = true;
    for (auto& gap : get_functional_description().get_gapping()) {
        if (!gap.get_coordinates()) {
            gapIsComplete = false;
        }
    }

    if (!gapIsComplete) {
        process_gap();
    }

    for (auto& gap : get_functional_description().get_gapping()) {
        if (gap.get_coordinates().value()[0] >= columnLeftLimit &&
            gap.get_coordinates().value()[0] <= columnRightLimit &&
            gap.get_coordinates().value()[2] <= columnFrontLimit &&
            gap.get_coordinates().value()[2] >= columnBackLimit) {
            foundGaps.push_back(gap);
        }
    }
    return foundGaps;
}

void Core::scale_to_stacks(int64_t numberStacks) {
    auto processedDescription = get_processed_description().value();
    processedDescription.get_mutable_effective_parameters().set_effective_area(
        processedDescription.get_effective_parameters().get_effective_area() * numberStacks);
    processedDescription.get_mutable_effective_parameters().set_minimum_area(
        processedDescription.get_effective_parameters().get_minimum_area() * numberStacks);
    processedDescription.get_mutable_effective_parameters().set_effective_volume(
        processedDescription.get_effective_parameters().get_effective_volume() * numberStacks);
    processedDescription.set_depth(processedDescription.get_depth() * numberStacks);
    for (auto& column : processedDescription.get_mutable_columns()) {
        column.set_area(column.get_area() * numberStacks);
        column.set_depth(column.get_depth() * numberStacks);
    }
    set_processed_description(processedDescription);

    auto functionalDescription = get_functional_description();
    auto gapping = functionalDescription.get_gapping();
    std::vector<CoreGap> scaledGapping;
    for (auto& gapInfo : gapping) {
        if (gapInfo.get_section_dimensions()) {
            auto gapSectionDimensions = *(gapInfo.get_section_dimensions());
            gapSectionDimensions[1] *= numberStacks;
            gapInfo.set_section_dimensions(gapSectionDimensions);
        }
        if (gapInfo.get_area()) {
            gapInfo.set_area(gapInfo.get_area().value() * numberStacks);
        }
        scaledGapping.push_back(gapInfo);
    }
    functionalDescription.set_gapping(scaledGapping);
    set_functional_description(functionalDescription);
}

void Core::set_gap_length(double gapLength) {
    for (size_t i = 0; i < get_functional_description().get_gapping().size(); ++i) {
        if (get_functional_description().get_gapping()[i].get_type() != GapType::RESIDUAL) {
            get_mutable_functional_description().get_mutable_gapping()[i].set_length(gapLength);
        }
    }
    distribute_and_process_gap();
}

bool Core::distribute_and_process_gap() {
    auto constants = Constants();
    std::vector<CoreGap> newGapping;
    auto gapping = get_functional_description().get_gapping();
    double centralColumnGapsHeightOffset;
    double distanceClosestNormalSurface;
    double coreChunkSizePlusGap = 0;
    auto nonResidualGaps = find_gaps_by_type(GapType::SUBTRACTIVE);
    auto additiveGaps = find_gaps_by_type(GapType::ADDITIVE);
    nonResidualGaps.insert(nonResidualGaps.end(), additiveGaps.begin(), additiveGaps.end());
    auto residualGaps = find_gaps_by_type(GapType::RESIDUAL);
    int numberNonResidualGaps = nonResidualGaps.size();
    int numberResidualGaps = residualGaps.size();
    int numberGaps = numberNonResidualGaps + numberResidualGaps;
    auto processedDescription = get_processed_description().value();
    auto columns = processedDescription.get_columns();
    int numberColumns = columns.size();

    // Check if current information is valid
    // if (numberNonResidualGaps + numberResidualGaps < numberColumns && (numberNonResidualGaps + numberResidualGaps) >
    // 0) {
    //     numberNonResidualGaps = 0;
    //     numberResidualGaps = 0;
    //     // if (!(get_functional_description().get_type() == CoreType::TOROIDAL && numberColumns == 1))
    //     {
    //     //     throw std::runtime_error("A TWO_PIECE_SET core cannot have less gaps than columns");
    //     // }
    // }

    if (numberNonResidualGaps == 0 && numberResidualGaps > numberColumns) {
        gapping = std::vector<CoreGap>(gapping.begin(), gapping.end() - (numberResidualGaps - numberColumns));
        get_mutable_functional_description().set_gapping(gapping);
        residualGaps = find_gaps_by_type(GapType::RESIDUAL);
        numberResidualGaps = residualGaps.size();
        numberGaps = numberNonResidualGaps + numberResidualGaps;
    }

    if (numberNonResidualGaps + numberResidualGaps == 0) {
        for (size_t i = 0; i < columns.size(); ++i) {
            CoreGap gap;
            gap.set_type(GapType::RESIDUAL);
            gap.set_length(constants.residualGap);
            gap.set_coordinates(columns[i].get_coordinates());
            gap.set_shape(columns[i].get_shape());
            if (columns[i].get_height() / 2 - constants.residualGap / 2 < 0) {
                return false;
                // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", column of index: " + std::to_string(i));

            }
            gap.set_distance_closest_normal_surface(columns[i].get_height() / 2 - constants.residualGap / 2);
            gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
            gap.set_area(columns[i].get_area());
            gap.set_section_dimensions(std::vector<double>({columns[i].get_width(), columns[i].get_depth()}));
            newGapping.push_back(gap);
        }
    }
    else if (numberNonResidualGaps + numberResidualGaps < numberColumns) {
        for (size_t i = 0; i < columns.size(); ++i) {
            size_t gapIndex = i;
            if (i >= gapping.size()) {
                gapIndex = gapping.size() - 1;
            }
            CoreGap gap;
            gap.set_type(gapping[gapIndex].get_type());
            gap.set_length(gapping[gapIndex].get_length());
            gap.set_coordinates(columns[i].get_coordinates());
            gap.set_shape(columns[i].get_shape());
            if (columns[i].get_height() / 2 - gapping[gapIndex].get_length() / 2 < 0) {
                return false;
                // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", column of index: " + std::to_string(i));

            }
            gap.set_distance_closest_normal_surface(columns[i].get_height() / 2 - gapping[gapIndex].get_length() / 2);
            gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
            gap.set_area(columns[i].get_area());
            gap.set_section_dimensions(std::vector<double>({columns[i].get_width(), columns[i].get_depth()}));
            newGapping.push_back(gap);
        }
    }
    else if ((numberResidualGaps == numberColumns || numberNonResidualGaps == numberColumns) &&
             (numberGaps == numberColumns)) {
        for (size_t i = 0; i < columns.size(); ++i) {
            CoreGap gap;
            gap.set_type(gapping[i].get_type());
            gap.set_length(gapping[i].get_length());
            gap.set_coordinates(columns[i].get_coordinates());
            gap.set_shape(columns[i].get_shape());
            if (columns[i].get_height() / 2 - gapping[i].get_length() / 2 < 0) {
                return false;
                // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", column of index: " + std::to_string(i));

            }
            gap.set_distance_closest_normal_surface(columns[i].get_height() / 2 - gapping[i].get_length() / 2);
            gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
            gap.set_area(columns[i].get_area());
            gap.set_section_dimensions(std::vector<double>({columns[i].get_width(), columns[i].get_depth()}));
            newGapping.push_back(gap);
        }
    }
    else {
        auto lateralColumns = find_columns_by_type(ColumnType::LATERAL);
        auto centralColumns = find_columns_by_type(ColumnType::CENTRAL);

        ColumnElement windingColumn;
        std::vector<ColumnElement> returnColumns;
        if (centralColumns.size() == 0) {
            windingColumn = lateralColumns[0];
            returnColumns = std::vector<ColumnElement>(lateralColumns.begin() + 1, lateralColumns.end());
        }
        else {
            windingColumn = centralColumns[0];
            returnColumns = lateralColumns;
        }

        if (numberGaps == numberColumns) {
            if (windingColumn.get_height() > nonResidualGaps[0].get_length()) {
                centralColumnGapsHeightOffset = roundFloat(nonResidualGaps[0].get_length() / 2);
            }
            else {
                centralColumnGapsHeightOffset = 0;
            }
            distanceClosestNormalSurface = roundFloat(windingColumn.get_height() / 2 - nonResidualGaps[0].get_length() / 2);
        }
        else {
            coreChunkSizePlusGap = roundFloat(windingColumn.get_height() / (nonResidualGaps.size() + 1));
            centralColumnGapsHeightOffset = roundFloat(-coreChunkSizePlusGap * (nonResidualGaps.size() - 1) / 2);
            distanceClosestNormalSurface = roundFloat(coreChunkSizePlusGap - nonResidualGaps[0].get_length() / 2);
        }

        for (size_t i = 0; i < nonResidualGaps.size(); ++i) {
            CoreGap gap;
            gap.set_type(nonResidualGaps[i].get_type());
            gap.set_length(nonResidualGaps[i].get_length());
            gap.set_coordinates(std::vector<double>({windingColumn.get_coordinates()[0],
                                      windingColumn.get_coordinates()[1] + centralColumnGapsHeightOffset,
                                      windingColumn.get_coordinates()[2]}));
            gap.set_shape(windingColumn.get_shape());
            if (distanceClosestNormalSurface < 0) {
                return false;
                // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", non residual gap of index: " + std::to_string(i));

            }
            gap.set_distance_closest_normal_surface(distanceClosestNormalSurface);
            gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
            gap.set_area(windingColumn.get_area());
            gap.set_section_dimensions(std::vector<double>({windingColumn.get_width(), windingColumn.get_depth()}));
            newGapping.push_back(gap);

            centralColumnGapsHeightOffset += roundFloat(windingColumn.get_height() / (nonResidualGaps.size() + 1));
            if (i < nonResidualGaps.size() / 2. - 1) {
                distanceClosestNormalSurface = roundFloat(distanceClosestNormalSurface + coreChunkSizePlusGap);
            }
            else if (i > nonResidualGaps.size() / 2. - 1) {
                distanceClosestNormalSurface = roundFloat(distanceClosestNormalSurface - coreChunkSizePlusGap);
            }
        }

        if (residualGaps.size() < returnColumns.size()) {
            for (size_t i = 0; i < returnColumns.size(); ++i) {
                CoreGap gap;
                gap.set_type(GapType::RESIDUAL);
                gap.set_length(constants.residualGap);
                gap.set_coordinates(returnColumns[i].get_coordinates());
                gap.set_shape(returnColumns[i].get_shape());
                if (returnColumns[i].get_height() / 2 - constants.residualGap / 2 < 0) {
                    return false;
                    // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", return column of index: " + std::to_string(i));

                }
                gap.set_distance_closest_normal_surface(returnColumns[i].get_height() / 2 - constants.residualGap / 2);
                gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
                gap.set_area(returnColumns[i].get_area());
                gap.set_section_dimensions(std::vector<double>({returnColumns[i].get_width(), returnColumns[i].get_depth()}));
                newGapping.push_back(gap);
            }
        }
        else {
            for (size_t i = 0; i < returnColumns.size(); ++i) {
                CoreGap gap;
                gap.set_type(residualGaps[i].get_type());
                gap.set_length(residualGaps[i].get_length());
                gap.set_coordinates(returnColumns[i].get_coordinates());
                gap.set_shape(returnColumns[i].get_shape());
                if (returnColumns[i].get_height() / 2 < 0) {
                    return false;
                    // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", return column of index: " + std::to_string(i));

                }
                gap.set_distance_closest_normal_surface(returnColumns[i].get_height() / 2);
                gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
                gap.set_area(returnColumns[i].get_area());
                gap.set_section_dimensions(std::vector<double>({returnColumns[i].get_width(), returnColumns[i].get_depth()}));
                newGapping.push_back(gap);
            }
        }
    }

    get_mutable_functional_description().set_gapping(newGapping);
    return true;
}

void Core::set_ground_gap(double gapLength) {
    size_t numberColumns = get_columns().size();
    std::vector<CoreGap> gapping;
    {
        CoreGap gap;
        gap.set_type(GapType::SUBTRACTIVE);
        gap.set_length(gapLength);
        gapping.push_back(gap);
    }
    for (size_t i = 0; i < numberColumns - 1; ++i) {
        CoreGap gap;
        gap.set_type(GapType::RESIDUAL);
        gap.set_length(constants.residualGap);
        gapping.push_back(gap);
    }
    get_mutable_functional_description().set_gapping(gapping);
}

void Core::set_distributed_gap(double gapLength, size_t numberGaps) {
    size_t numberColumns = get_columns().size();
    std::vector<CoreGap> gapping;
    for (size_t i = 0; i < numberGaps; ++i) {
        CoreGap gap;
        gap.set_type(GapType::SUBTRACTIVE);
        gap.set_length(gapLength);
        gapping.push_back(gap);
    }
    for (size_t i = 0; i < numberColumns - 1; ++i) {
        CoreGap gap;
        gap.set_type(GapType::RESIDUAL);
        gap.set_length(constants.residualGap);
        gapping.push_back(gap);
    }
    get_mutable_functional_description().set_gapping(gapping);
}

void Core::set_spacer_gap(double gapLength) {
    size_t numberColumns = get_columns().size();
    std::vector<CoreGap> gapping;
    for (size_t i = 0; i < numberColumns; ++i) {
        CoreGap gap;
        gap.set_type(GapType::ADDITIVE);
        gap.set_length(gapLength);
        gapping.push_back(gap);
    }
    get_mutable_functional_description().set_gapping(gapping);
}

void Core::set_residual_gap() {
    size_t numberColumns = get_columns().size();
    std::vector<CoreGap> gapping;
    for (size_t i = 0; i < numberColumns; ++i) {
        CoreGap gap;
        gap.set_type(GapType::RESIDUAL);
        gap.set_length(constants.residualGap);
        gapping.push_back(gap);
    }
    get_mutable_functional_description().set_gapping(gapping);
}


bool Core::is_gapping_missaligned() {
    auto gapping = get_functional_description().get_gapping();
    for (size_t i = 0; i < gapping.size(); ++i) {
        if (!gapping[i].get_coordinates()) {
            return true;
        }

        auto columnIndex = find_exact_column_index_by_coordinates(*gapping[i].get_coordinates());
        if (columnIndex == -1) {
            return true;
        }
    }
    return false;
}

bool Core::is_gap_processed() {
    for (auto& gap : get_gapping()) {
        if (!gap.get_coordinates()) {
            return false;
        }
    }
    return true;
}

bool Core::process_gap() {
    std::vector<CoreGap> newGapping;
    auto gapping = get_functional_description().get_gapping();
    auto family = std::get<CoreShape>(get_functional_description().get_shape()).get_family();
    auto processedDescription = get_processed_description().value();
    auto columns = processedDescription.get_columns();

    if (family == CoreShapeFamily::T && gapping.size() > 0 ) {
        throw std::runtime_error("Toroids cannot be gapped: " + std::to_string(gapping[0].get_length()));
    }

    if (family != CoreShapeFamily::T) {
        if (gapping.size() == 0 || !gapping[0].get_coordinates() || is_gapping_missaligned()) {
            return distribute_and_process_gap();
        }

        for (size_t i = 0; i < gapping.size(); ++i) {
            CoreGap gap;
            auto columnIndex = find_closest_column_index_by_coordinates(*gapping[i].get_coordinates());

            gap.set_type(gapping[i].get_type());
            gap.set_length(gapping[i].get_length());
            gap.set_coordinates(gapping[i].get_coordinates());
            gap.set_shape(columns[columnIndex].get_shape());
            if (roundFloat(columns[columnIndex].get_height() / 2 - fabs((*gapping[i].get_coordinates())[1]) - gapping[i].get_length() / 2) < 0) {
                return false;
                // throw std::runtime_error("distance_closest_normal_surface cannot be negative in shape: " + std::get<CoreShape>(get_functional_description().get_shape()).get_name().value() + ", gap of index: " + std::to_string(i));

            }
            gap.set_distance_closest_normal_surface(roundFloat(columns[columnIndex].get_height() / 2 - fabs((*gapping[i].get_coordinates())[1]) - gapping[i].get_length() / 2));
            gap.set_distance_closest_parallel_surface(processedDescription.get_winding_windows()[0].get_width());
            gap.set_area(columns[columnIndex].get_area());
            gap.set_section_dimensions(std::vector<double>({columns[columnIndex].get_width(), columns[columnIndex].get_depth()}));
            newGapping.push_back(gap);
        }
    }

    get_mutable_functional_description().set_gapping(newGapping);
    return true;
}

CoreMaterial Core::resolve_material() {
    auto material = resolve_material(get_functional_description().get_material());
    get_mutable_functional_description().set_material(material);
    return material;
}

CoreMaterial Core::resolve_material(CoreMaterialDataOrNameUnion coreMaterial) {
    // If the material is a string, we have to load its data from the database
    if (std::holds_alternative<std::string>(coreMaterial)) {
        if (std::get<std::string>(coreMaterial) == "dummy" || std::get<std::string>(coreMaterial) == "Dummy") {
            CoreMaterial coreMaterial;
            coreMaterial.set_name("Dummy");
            return coreMaterial;
        }
        auto coreMaterialData = find_core_material_by_name(std::get<std::string>(coreMaterial));
        coreMaterial = coreMaterialData;
        return coreMaterialData;
    }
    else {
        return std::get<CoreMaterial>(coreMaterial);
    }
}

void Core::set_material(CoreMaterial coreMaterial) {
    get_mutable_functional_description().set_material(coreMaterial);
}

void Core::set_material_initial_permeability(double value) {
    auto coreMaterial = resolve_material();
    PermeabilityPoint permeabilityPoint;
    permeabilityPoint.set_value(value);
    coreMaterial.get_mutable_permeability().set_initial(permeabilityPoint);
}


CoreShape Core::resolve_shape() {
    auto shape = resolve_shape(get_functional_description().get_shape());
    get_mutable_functional_description().set_shape(shape);
    return shape;
}

CoreShape Core::resolve_shape(CoreShapeDataOrNameUnion coreShape) {
    // If the shape is a string, we have to load its data from the database
    if (std::holds_alternative<std::string>(coreShape)) {
        auto coreShapeData = find_core_shape_by_name(std::get<std::string>(coreShape));
        coreShape = coreShapeData;
        return coreShapeData;
    }
    else {
        return std::get<CoreShape>(coreShape);
    }
}

void Core::process_data() {
    // If the shape is a string, we have to load its data from the database
    if (std::holds_alternative<std::string>(get_functional_description().get_shape())) {

        auto shape_data = find_core_shape_by_name(std::get<std::string>(get_functional_description().get_shape()));
        shape_data.set_name(std::get<std::string>(get_functional_description().get_shape()));
        get_mutable_functional_description().set_shape(shape_data);
    }

    // If the material is a string, we have to load its data from the database, unless it is dummy (in order to avoid
    // long loading operatings)
    if (_includeMaterialData && std::holds_alternative<std::string>(get_functional_description().get_material()) &&
        std::get<std::string>(get_functional_description().get_material()) != "dummy" &&
        std::get<std::string>(get_functional_description().get_material()) != "Dummy") {
        auto material_data = find_core_material_by_name(
            std::get<std::string>(get_functional_description().get_material()));
        get_mutable_functional_description().set_material(material_data);
    }

    auto corePiece =
        CorePiece::factory(std::get<CoreShape>(get_functional_description().get_shape()));
    CoreProcessedDescription processedDescription;
    auto coreColumns = corePiece->get_columns();
    auto coreWindingWindow = corePiece->get_winding_window();
    auto coreEffectiveParameters = corePiece->get_partial_effective_parameters();
    switch (get_functional_description().get_type()) {
        case CoreType::TOROIDAL:
        case CoreType::CLOSED_SHAPE:
            processedDescription.set_columns(coreColumns);
            processedDescription.set_effective_parameters(corePiece->get_partial_effective_parameters());
            processedDescription.get_mutable_winding_windows().push_back(corePiece->get_winding_window());
            processedDescription.set_depth(corePiece->get_depth());
            processedDescription.set_height(corePiece->get_height());
            processedDescription.set_width(corePiece->get_width());
            break;

        case CoreType::TWO_PIECE_SET:
            for (auto& column : coreColumns) {
                column.set_height(2 * column.get_height());
            }
            processedDescription.set_columns(coreColumns);

            coreEffectiveParameters.set_effective_length(2 * coreEffectiveParameters.get_effective_length());
            coreEffectiveParameters.set_effective_area(coreEffectiveParameters.get_effective_area());
            coreEffectiveParameters.set_effective_volume(2 * coreEffectiveParameters.get_effective_volume());
            coreEffectiveParameters.set_minimum_area(coreEffectiveParameters.get_minimum_area());
            processedDescription.set_effective_parameters(coreEffectiveParameters);

            coreWindingWindow.set_area(2 * coreWindingWindow.get_area().value());
            coreWindingWindow.set_height(2 * coreWindingWindow.get_height().value());
            processedDescription.get_mutable_winding_windows().push_back(coreWindingWindow);
            processedDescription.set_depth(corePiece->get_depth());
            processedDescription.set_height(corePiece->get_height() * 2);
            processedDescription.set_width(corePiece->get_width());
            break;
        default:
            throw std::runtime_error("Unknown type of core, available options are {TOROIDAL, TWO_PIECE_SET}");
    }
    set_processed_description(processedDescription);
    scale_to_stacks(*(get_functional_description().get_number_stacks()));
}

double Core::get_magnetic_flux_density_saturation(CoreMaterial coreMaterial, double temperature, bool proportion) {
    auto defaults = Defaults();
    auto saturationData = coreMaterial.get_saturation();
    std::vector<std::pair<double, double>> data;

    if (saturationData.size() == 0) {
        return defaults.magneticFluxDensitySaturation;
        // throw std::runtime_error("Missing saturation data in core material");
    }
    for (auto& datum : saturationData) {
        data.push_back(std::pair<double, double>{datum.get_temperature(), datum.get_magnetic_flux_density()});
    }

    double saturationMagneticFluxDensity = interp(data, temperature);

    if (proportion)
        return defaults.maximumProportionMagneticFluxDensitySaturation * saturationMagneticFluxDensity;
    else
        return saturationMagneticFluxDensity;
}

double Core::get_magnetic_flux_density_saturation(double temperature, bool proportion) {
    auto coreMaterial = resolve_material();
    return get_magnetic_flux_density_saturation(coreMaterial, temperature, proportion);
}

double Core::get_magnetic_flux_density_saturation(CoreMaterial coreMaterial, bool proportion) {
    return get_magnetic_flux_density_saturation(coreMaterial, 25, proportion);
}
double Core::get_magnetic_flux_density_saturation(bool proportion) {
    return get_magnetic_flux_density_saturation(25, proportion);
}

double Core::get_magnetic_field_strength_saturation(CoreMaterial coreMaterial, double temperature) {
    auto saturationData = coreMaterial.get_saturation();
    std::vector<std::pair<double, double>> data;

    if (saturationData.size() == 0) {
        return defaults.magneticFluxDensitySaturation;
        // throw std::runtime_error("Missing saturation data in core material");
    }
    for (auto& datum : saturationData) {
        data.push_back(std::pair<double, double>{datum.get_temperature(), datum.get_magnetic_field()});
    }

    double saturationMagneticFieldStrength = interp(data, temperature);

    return saturationMagneticFieldStrength;
}

double Core::get_magnetic_field_strength_saturation(double temperature) {
    auto coreMaterial = resolve_material();
    return get_magnetic_field_strength_saturation(coreMaterial, temperature);
}

double Core::get_remanence(CoreMaterial coreMaterial, double temperature, bool returnZeroIfMissing) {
    if (!coreMaterial.get_remanence()) {
        if (returnZeroIfMissing) {
            return 0;
        }
        else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    auto remanenceData = coreMaterial.get_remanence().value();
    std::vector<std::pair<double, double>> data;

    if (remanenceData.size() == 0) {
        throw std::runtime_error("Missing remanence data in core material");
    }
    for (auto& datum : remanenceData) {
        data.push_back(std::pair<double, double>{datum.get_temperature(), datum.get_magnetic_flux_density()});
    }

    double remanence = interp(data, temperature);

    return remanence;
}

double Core::get_remanence(double temperature) {
    auto coreMaterial = resolve_material();
    return get_remanence(coreMaterial, temperature);
}

double Core::get_curie_temperature(CoreMaterial coreMaterial) {
    if (!coreMaterial.get_curie_temperature()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return coreMaterial.get_curie_temperature().value();
}

double Core::get_curie_temperature() {
    auto coreMaterial = resolve_material();
    return get_curie_temperature(coreMaterial);
}

double Core::get_coercive_force(CoreMaterial coreMaterial, double temperature, bool returnZeroIfMissing) {
    if (!coreMaterial.get_coercive_force()) {
        if (returnZeroIfMissing) {
            return 0;
        }
        else {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    auto coerciveForceData = coreMaterial.get_coercive_force().value();
    std::vector<std::pair<double, double>> data;

    if (coerciveForceData.size() == 0) {
        throw std::runtime_error("Missing coercive force data in core material");
    }
    for (auto& datum : coerciveForceData) {
        data.push_back(std::pair<double, double>{datum.get_temperature(), datum.get_magnetic_field()});
    }

    double coerciveForce = interp(data, temperature);

    return coerciveForce;
}

double Core::get_coercive_force(double temperature) {
    auto coreMaterial = resolve_material();
    return get_coercive_force(coreMaterial, temperature);
}


double Core::get_initial_permeability(double temperature){
    auto coreMaterial = resolve_material();
    return get_initial_permeability(coreMaterial, temperature);
}

double Core::get_initial_permeability(CoreMaterial coreMaterial, double temperature){
    InitialPermeability initialPermeability;
    auto initialPermeabilityValue = initialPermeability.get_initial_permeability(coreMaterial, temperature, std::nullopt, std::nullopt);
    return initialPermeabilityValue;
}

double Core::get_effective_permeability(double temperature){
    auto constants = Constants();
    auto reluctance = get_reluctance(temperature);
    if (!get_processed_description()) {
        process_data();
    }
    auto effectiveLength = get_processed_description().value().get_effective_parameters().get_effective_length();
    auto effectiveArea = get_processed_description().value().get_effective_parameters().get_effective_area();
    double effectivePermeability = effectiveLength / (reluctance * effectiveArea) / constants.vacuumPermeability;
    return effectivePermeability;
}

double Core::get_reluctance(double temperature){
    auto coreMaterial = resolve_material();
    InitialPermeability initialPermeability;
    auto initialPermeabilityValue = initialPermeability.get_initial_permeability(coreMaterial, temperature, std::nullopt, std::nullopt);
    auto reluctanceModel = ReluctanceModel::factory();

    auto magnetizingInductanceOutput = reluctanceModel->get_core_reluctance(*this, initialPermeabilityValue);
    double calculatedReluctance = magnetizingInductanceOutput.get_core_reluctance();
    return calculatedReluctance;
}

double Core::get_resistivity(double temperature){
    auto coreMaterial = resolve_material();
    return get_resistivity(coreMaterial, temperature);
}

double Core::get_resistivity(CoreMaterial coreMaterial, double temperature){
    auto resistivityModel = ResistivityModel::factory(ResistivityModels::CORE_MATERIAL);
    auto resistivity = (*resistivityModel).get_resistivity(coreMaterial, temperature);
    return resistivity;
}

double Core::get_density(CoreMaterial coreMaterial) {
    if (!coreMaterial.get_density()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return coreMaterial.get_density().value();
}

double Core::get_density() {
    auto coreMaterial = resolve_material();
    return get_density(coreMaterial);
}

std::vector<ColumnElement> Core::get_columns() {
    if (get_processed_description()) {
        return get_processed_description().value().get_columns();
    }
    else {
        return std::vector<ColumnElement>();
    }
}

WindingWindowElement Core::get_winding_window(size_t windingWindowIndex) {
    if (get_processed_description()) {
        return get_processed_description().value().get_winding_windows()[windingWindowIndex];
    }
    else {
        return WindingWindowElement();
    }
}

std::vector<WindingWindowElement> Core::get_winding_windows() {
    if (get_processed_description()) {
        return get_processed_description().value().get_winding_windows();
    }
    else {
        return std::vector<WindingWindowElement>();
    }
}

std::string Core::get_material_family() {
    if (resolve_material().get_family())
        return resolve_material().get_family().value();
    else 
        return ""; 
}

CoreShapeFamily Core::get_shape_family() {
    return resolve_shape().get_family();
}

std::string Core::get_shape_name() {
    if (resolve_shape().get_name())
        return resolve_shape().get_name().value();
    else
        return "Custom";
}

int64_t Core::get_number_stacks() {
    if (get_functional_description().get_number_stacks()) {
        return get_functional_description().get_number_stacks().value();
    }
    else {
        return 1;
    }
}

std::string Core::get_material_name() {
    return resolve_material().get_name();
}

std::vector<VolumetricCoreLossesMethodType> Core::get_available_core_losses_methods(){
    auto coreMaterial = resolve_material();
    return get_available_core_losses_methods(coreMaterial);
}

std::vector<VolumetricCoreLossesMethodType> Core::get_available_core_losses_methods(CoreMaterial coreMaterial){
    std::vector<VolumetricCoreLossesMethodType> methods;
    auto volumetricLossesMethodsVariants = coreMaterial.get_volumetric_losses();
    for (auto& volumetricLossesMethodVariant : volumetricLossesMethodsVariants) {
        auto volumetricLossesMethods = volumetricLossesMethodVariant.second;
        for (auto& volumetricLossesMethod : volumetricLossesMethods) {
            if (std::holds_alternative<CoreLossesMethodData>(volumetricLossesMethod)) {
                auto methodData = std::get<CoreLossesMethodData>(volumetricLossesMethod);
                if (std::find(methods.begin(), methods.end(), methodData.get_method()) == methods.end()) {
                    methods.push_back(methodData.get_method());
                }
            }
        }
    }
    return methods;
}


Application Core::resolve_material_application() {
    auto coreMaterial = resolve_material();
    auto application = resolve_material_application(coreMaterial);
    get_mutable_functional_description().set_material(coreMaterial);
    return application;
}

Application Core::resolve_material_application(CoreMaterial& coreMaterial) {
    if (coreMaterial.get_application()) {
        return coreMaterial.get_application().value();
    }
    else {
        coreMaterial.set_application(guess_material_application(coreMaterial));
        return coreMaterial.get_application().value();
    }
}

Application Core::guess_material_application() {
    auto coreMaterial = resolve_material();
    return guess_material_application(coreMaterial);
}

Application Core::guess_material_application(CoreMaterial coreMaterial) {
    for (auto method : get_available_core_losses_methods(coreMaterial)) {
        if (method == VolumetricCoreLossesMethodType::LOSS_FACTOR) {
            if (coreMaterial.get_permeability().get_complex()) {
                return Application::INTERFERENCE_SUPPRESSION;
            }
            return Application::SIGNAL_PROCESSING;
        }
    }
    return Application::POWER;
}

bool Core::check_material_application(Application application) {
    auto coreMaterial = resolve_material();
    return check_material_application(coreMaterial, application);
}

bool Core::check_material_application(CoreMaterial coreMaterial, Application application) {
    if (coreMaterial.get_permeability().get_complex()) {
        if (application == Application::INTERFERENCE_SUPPRESSION) {
            return true;
        }
    }

    for (auto method : get_available_core_losses_methods(coreMaterial)) {
        if (method == VolumetricCoreLossesMethodType::LOSS_FACTOR) {
            if (application == Application::SIGNAL_PROCESSING) {
                return true;
            }
        }
        else {
            if (application == Application::POWER) {
                return true;
            }
        }
    }
    return false;
}

Application Core::guess_material_application(std::string coreMaterialName) {
    auto coreMaterial = find_core_material_by_name(coreMaterialName);
    return guess_material_application(coreMaterial);
}

CoreType Core::get_type() {
    return get_functional_description().get_type();
}

bool fits_one_dimension(CoreProcessedDescription coreProcessedDescription, double dimension) {
    if (coreProcessedDescription.get_depth() <= dimension || coreProcessedDescription.get_height() <= dimension || coreProcessedDescription.get_width() <= dimension) {
        return true;
    }
    return false;
}

bool fits_two_dimensions(CoreProcessedDescription coreProcessedDescription, double firstDimension, double secondDimension) {
    if ((coreProcessedDescription.get_depth() <= firstDimension && (coreProcessedDescription.get_height() <= secondDimension || coreProcessedDescription.get_width() <= secondDimension)) ||
        (coreProcessedDescription.get_height() <= firstDimension && (coreProcessedDescription.get_depth() <= secondDimension || coreProcessedDescription.get_width() <= secondDimension)) ||
        (coreProcessedDescription.get_width() <= firstDimension && (coreProcessedDescription.get_height() <= secondDimension || coreProcessedDescription.get_depth() <= secondDimension))) {
        return true;
    }
    return false;
}

bool fits_three_dimensions(CoreProcessedDescription coreProcessedDescription, double firstDimension, double secondDimension, double thirdDimension) {
    if ((coreProcessedDescription.get_depth() <= firstDimension && coreProcessedDescription.get_height() <= secondDimension && coreProcessedDescription.get_width() <= thirdDimension) ||
        (coreProcessedDescription.get_depth() <= firstDimension && coreProcessedDescription.get_height() <= thirdDimension && coreProcessedDescription.get_width() <= secondDimension) ||
        (coreProcessedDescription.get_depth() <= secondDimension && coreProcessedDescription.get_height() <= firstDimension && coreProcessedDescription.get_width() <= thirdDimension) ||
        (coreProcessedDescription.get_depth() <= secondDimension && coreProcessedDescription.get_height() <= thirdDimension && coreProcessedDescription.get_width() <= firstDimension) ||
        (coreProcessedDescription.get_depth() <= thirdDimension && coreProcessedDescription.get_height() <= firstDimension && coreProcessedDescription.get_width() <= secondDimension) ||
        (coreProcessedDescription.get_depth() <= thirdDimension && coreProcessedDescription.get_height() <= secondDimension && coreProcessedDescription.get_width() <= firstDimension)){
        return true;
    }
    return false;
}

bool Core::fits(MaximumDimensions maximumDimensions, bool allowRotation) {
    if (!get_processed_description()) {
        process_data();
    }

    auto coreProcessedDescription = get_processed_description().value();

    if (!maximumDimensions.get_depth() && !maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        return true;
    }
    else if (maximumDimensions.get_depth() && !maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        auto depth = maximumDimensions.get_depth().value();
        if (allowRotation) {
            return fits_one_dimension(coreProcessedDescription, depth);
        }
        else {
            return coreProcessedDescription.get_depth() <= depth;
        }
    }
    else if (!maximumDimensions.get_depth() && maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        auto height = maximumDimensions.get_height().value();
        if (allowRotation) {
            return fits_one_dimension(coreProcessedDescription, height);
        }
        else {
            return coreProcessedDescription.get_height() <= height;
        }
    }
    else if (!maximumDimensions.get_depth() && !maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto width = maximumDimensions.get_width().value();
        if (allowRotation) {
            return fits_one_dimension(coreProcessedDescription, width);
        }
        else {
            return coreProcessedDescription.get_width() <= width;
        }
    }
    else if (maximumDimensions.get_depth() && maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        auto depth = maximumDimensions.get_depth().value();
        auto height = maximumDimensions.get_height().value();
        if (allowRotation) {
            return fits_two_dimensions(coreProcessedDescription, depth, height);
        }
        else {
            return coreProcessedDescription.get_depth() <= depth && coreProcessedDescription.get_height() <= height;
        }
    }
    else if (!maximumDimensions.get_depth() && maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto width = maximumDimensions.get_width().value();
        auto height = maximumDimensions.get_height().value();
        if (allowRotation) {
            return fits_two_dimensions(coreProcessedDescription, width, height);
        }
        else {
            return coreProcessedDescription.get_width() <= width && coreProcessedDescription.get_height() <= height;
        }
    }
    else if (maximumDimensions.get_depth() && !maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto width = maximumDimensions.get_width().value();
        auto depth = maximumDimensions.get_depth().value();
        if (allowRotation) {
            return fits_two_dimensions(coreProcessedDescription, width, depth);
        }
        else {
            return coreProcessedDescription.get_depth() <= depth && coreProcessedDescription.get_width() <= width;
        }
    }
    else if (maximumDimensions.get_depth() && maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto depth = maximumDimensions.get_depth().value();
        auto height = maximumDimensions.get_height().value();
        auto width = maximumDimensions.get_width().value();
        if (allowRotation) {
            return fits_three_dimensions(coreProcessedDescription, depth, height, width);
        }
        else {
            return coreProcessedDescription.get_depth() <= depth && coreProcessedDescription.get_height() <= height && coreProcessedDescription.get_width() <= width;
        }
    }
    else {
        throw std::runtime_error("Not sure how this happened");
    }
}

std::vector<double> Core::get_maximum_dimensions() {
    if (!get_processed_description()) {
        process_data();
    }

    auto coreProcessedDescription = get_processed_description().value();
    return {coreProcessedDescription.get_width(), coreProcessedDescription.get_height(), coreProcessedDescription.get_depth()};
}

} // namespace OpenMagnetics
