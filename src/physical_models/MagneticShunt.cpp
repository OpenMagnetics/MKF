#include "physical_models/MagneticShunt.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/InitialPermeability.h"
#include "physical_models/LeakageInductance.h"
#include "support/Exceptions.h"
#include "support/Utils.h"
#include "Constants.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <numbers>

namespace OpenMagnetics {

namespace {

// Geometric agreement allowed between a drawn position and a declared gap: one micrometre, or one
// percent of the declared value when that is larger. Positions come from millimetre-scale drawings;
// this only absorbs rounding, never a different design.
double geometric_tolerance(double declaredValue) {
    return std::max(1e-6, 0.01 * std::fabs(declaredValue));
}

double interval_overlap(double a0, double a1, double b0, double b1) {
    return std::max(0.0, std::min(a1, b1) - std::max(a0, b0));
}

std::string shunt_label(const MagneticShunt& shunt, size_t shuntIndex) {
    if (shunt.get_name()) {
        return "shunt '" + shunt.get_name().value() + "'";
    }
    return "shunt " + std::to_string(shuntIndex);
}

// Gap with fringing after Hurley & Wolfle, as used by Li et al. eq. (13): a gap of length g across a
// section a by b behaves as a section (a + g) by (b + g).
double fringed_gap_reluctance(double gapLength, double sectionWidth, double sectionDepth) {
    if (gapLength <= 0) {
        return 0;
    }
    double mu0 = Constants().vacuumPermeability;
    return gapLength / (mu0 * (sectionWidth + gapLength) * (sectionDepth + gapLength));
}

struct SheetBox {
    double x0, x1, y0, y1, z0, z1;
};

SheetBox sheet_box(const MagneticShunt& shunt, const std::string& label) {
    auto coordinates = shunt.get_coordinates();
    auto dimensions = shunt.get_dimensions();
    if (coordinates.size() < 2) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " needs at least x and y coordinates");
    }
    if (dimensions.size() < 3) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " needs width, height and depth");
    }
    double width = dimensions[0];
    double height = dimensions[1];
    double depth = dimensions[2];
    if (width <= 0 || height <= 0 || depth <= 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " has a non-positive dimension");
    }
    double z = coordinates.size() > 2 ? coordinates[2] : 0.0;
    return {coordinates[0] - width / 2, coordinates[0] + width / 2,
            coordinates[1] - height / 2, coordinates[1] + height / 2,
            z - depth / 2, z + depth / 2};
}

struct EColumns {
    ColumnElement main;
    ColumnElement right;
    ColumnElement left;
};

EColumns resolve_e_columns(Core& core) {
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("Magnetic shunt model needs a processed core");
    }
    auto columns = core.get_columns();
    if (columns.size() != 3) {
        throw NotImplementedException("Magnetic shunt model for a core with " + std::to_string(columns.size()) +
                                      " columns (only a central column with two lateral columns is modelled)");
    }
    size_t mainIndex = core.get_main_column_index();
    std::optional<ColumnElement> right;
    std::optional<ColumnElement> left;
    for (size_t index = 0; index < columns.size(); ++index) {
        if (index == mainIndex) {
            continue;
        }
        if (columns[index].get_coordinates()[0] > 0) {
            right = columns[index];
        }
        else {
            left = columns[index];
        }
    }
    auto main = columns[mainIndex];
    if (!right || !left) {
        throw NotImplementedException("Magnetic shunt model for a core whose lateral columns are not on both sides of the main column");
    }
    if (std::fabs(main.get_coordinates()[0]) > 1e-9) {
        throw NotImplementedException("Magnetic shunt model for a core whose main column is not at x = 0");
    }
    for (auto column : {main, right.value(), left.value()}) {
        if (column.get_shape() != ColumnShape::RECTANGULAR) {
            throw NotImplementedException("Magnetic shunt model for a core with non-rectangular columns (a sheet across a round or oblong column is an annulus, not modelled)");
        }
    }
    return {main, right.value(), left.value()};
}

// The gap of `column` that holds the sheet's axial range, or nullopt when there is none.
std::optional<CoreGap> gap_holding_sheet(Core& core, const ColumnElement& column, const SheetBox& box) {
    for (auto gap : core.find_gaps_by_column(column)) {
        if (!gap.get_coordinates()) {
            throw GapNotProcessedException("Magnetic shunt model needs processed gap coordinates");
        }
        double gapCentre = gap.get_coordinates().value()[1];
        double gapBottom = gapCentre - gap.get_length() / 2;
        double gapTop = gapCentre + gap.get_length() / 2;
        double tolerance = 1e-9;
        if (box.y0 >= gapBottom - tolerance && box.y1 <= gapTop + tolerance) {
            return gap;
        }
    }
    return std::nullopt;
}

} // namespace

CoreMaterial MagneticShuntModel::resolve_material(const MagneticShunt& shunt) {
    CoreMaterial material;
    if (std::holds_alternative<std::string>(shunt.get_material())) {
        material = find_core_material_by_name(std::get<std::string>(shunt.get_material()));
    }
    else {
        material = std::get<CoreMaterial>(shunt.get_material());
    }
    auto initial = material.get_permeability().get_initial();
    bool hasInitial = std::holds_alternative<PermeabilityPoint>(initial) ||
                      !std::get<std::vector<PermeabilityPoint>>(initial).empty();
    if (!hasInitial) {
        throw MaterialDataMissingException(material.get_name(), "initial permeability of the magnetic shunt material");
    }
    return material;
}

double MagneticShuntModel::get_relative_permeability(const CoreMaterial& material, double temperature, double frequency) {
    double relativePermeability = InitialPermeability::get_initial_permeability(material, temperature, std::nullopt, frequency);
    if (!(relativePermeability > 1.0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic shunt material " + material.get_name() +
                                    " resolves to a relative permeability of " + std::to_string(relativePermeability) +
                                    ", which carries no leakage flux");
    }
    return relativePermeability;
}

bool MagneticShuntModel::is_leakage_shunt(const MagneticShunt& shunt) {
    return shunt.get_placement() == MagneticShuntPlacement::IN_WINDOW ||
           shunt.get_placement() == MagneticShuntPlacement::BETWEEN_SECTIONS;
}

bool MagneticShuntModel::has_leakage_shunts(const Magnetic& magnetic) {
    if (!magnetic.get_shunts()) {
        return false;
    }
    auto shunts = magnetic.get_shunts().value();
    for (auto& shunt : shunts) {
        if (is_leakage_shunt(shunt)) {
            return true;
        }
    }
    return false;
}

void MagneticShuntModel::check_supported_placements(const Magnetic& magnetic) {
    if (!magnetic.get_shunts()) {
        return;
    }
    auto shunts = magnetic.get_shunts().value();
    for (size_t shuntIndex = 0; shuntIndex < shunts.size(); ++shuntIndex) {
        if (shunts[shuntIndex].get_placement() == MagneticShuntPlacement::OUTSIDE_WINDOW) {
            throw NotImplementedException("Magnetic " + shunt_label(shunts[shuntIndex], shuntIndex) +
                                          " placed outsideWindow (a flux path outside the core)");
        }
    }
}

MagneticShuntNetworkInputs MagneticShuntModel::extract_network_inputs(Magnetic& magnetic, size_t shuntIndex, double temperature, double frequency) {
    if (!magnetic.get_shunts() || shuntIndex >= magnetic.get_shunts()->size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic has no shunt " + std::to_string(shuntIndex));
    }
    auto shunt = magnetic.get_shunts().value()[shuntIndex];
    auto label = shunt_label(shunt, shuntIndex);
    if (!is_leakage_shunt(shunt)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " is not placed in the window (inWindow or betweenSections)");
    }
    auto& core = magnetic.get_mutable_core();
    auto columns = resolve_e_columns(core);
    auto box = sheet_box(shunt, label);

    MagneticShuntNetworkInputs inputs;
    inputs.name = shunt.get_name() ? shunt.get_name().value() : "Shunt " + std::to_string(shuntIndex);
    inputs.placement = shunt.get_placement();
    inputs.thickness = box.y1 - box.y0;
    inputs.depth = box.z1 - box.z0;
    inputs.axialCoordinate = (box.y0 + box.y1) / 2;
    double width = box.x1 - box.x0;
    if (width <= inputs.thickness) {
        throw NotImplementedException("Magnetic " + label + " is not wider than it is tall: a sheet carrying leakage flux along the column axis");
    }

    // Zhang eq. (10): R_c1 = R_c2 = l_e / (2 mu0 mur A_e); one half-core loop is R_c1/2 + R_c2.
    double mu0 = Constants().vacuumPermeability;
    auto effectiveParameters = core.get_processed_description()->get_effective_parameters();
    double coreRelativePermeability = InitialPermeability::get_initial_permeability(core.resolve_material(), temperature, std::nullopt, frequency);
    double coreHalfReluctance = effectiveParameters.get_effective_length() / (2 * mu0 * coreRelativePermeability * effectiveParameters.get_effective_area());
    inputs.coreLoopReluctance = coreHalfReluctance / 2 + coreHalfReluctance;

    double windowHalfHeight = columns.main.get_height() / 2;
    auto gapToColumns = shunt.get_gap_to_columns();

    for (int side : {1, -1}) {
        auto lateral = side > 0 ? columns.right : columns.left;
        // Everything in side coordinates u = side * x, so the window runs from the main column face
        // outwards on both sides.
        double u0 = std::min(side * box.x0, side * box.x1);
        double u1 = std::max(side * box.x0, side * box.x1);
        double innerFace = columns.main.get_width() / 2;
        double outerFace = side * lateral.get_coordinates()[0] - lateral.get_width() / 2;
        double outerEdge = side * lateral.get_coordinates()[0] + lateral.get_width() / 2;
        double spanInWindow = interval_overlap(u0, u1, innerFace, outerFace);
        double innerOverlapWidth = interval_overlap(u0, u1, 0.0, innerFace);
        double outerOverlapWidth = interval_overlap(u0, u1, outerFace, outerEdge);

        auto columnZOverlap = [&](const ColumnElement& column) {
            double columnZ = column.get_coordinates().size() > 2 ? column.get_coordinates()[2] : 0.0;
            return interval_overlap(box.z0, box.z1, columnZ - column.get_depth() / 2, columnZ + column.get_depth() / 2);
        };

        // A sheet inside a column must sit in that column's gap: anywhere else it is inside ferrite.
        auto overlapAir = [&](const ColumnElement& column, const std::string& columnName) -> std::pair<double, double> {
            auto gap = gap_holding_sheet(core, column, box);
            if (!gap) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " runs into the " + columnName +
                                            " column at y = " + std::to_string(inputs.axialCoordinate) +
                                            " m, where that column has no gap to hold it");
            }
            double gapCentre = gap->get_coordinates().value()[1];
            double airAbove = (gapCentre + gap->get_length() / 2) - box.y1;
            double airBelow = box.y0 - (gapCentre - gap->get_length() / 2);
            return {std::max(0.0, airAbove), std::max(0.0, airBelow)};
        };

        if (spanInWindow <= 0) {
            if (innerOverlapWidth > 0 || outerOverlapWidth > 0) {
                // Part of the sheet sits in a column gap on this side without reaching the window:
                // it still has to be held by a gap.
                if (innerOverlapWidth > 0) {
                    overlapAir(columns.main, "central");
                }
                if (outerOverlapWidth > 0) {
                    overlapAir(lateral, "lateral");
                }
            }
            continue;
        }

        if (inputs.axialCoordinate + inputs.thickness / 2 > windowHalfHeight + 1e-9 ||
            inputs.axialCoordinate - inputs.thickness / 2 < -windowHalfHeight - 1e-9) {
            if (innerOverlapWidth <= 0 && outerOverlapWidth <= 0) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " sticks out of the winding window axially");
            }
        }

        MagneticShuntNetworkInputs::Crossing crossing;
        crossing.side = side;
        crossing.spanInWindow = spanInWindow;
        crossing.windowWidth = outerFace - innerFace;
        crossing.depthInWindow = inputs.depth;
        crossing.sheetStart = u0;
        crossing.sheetEnd = u1;

        // Inner (central column) contact.
        if (innerOverlapWidth > 0) {
            if (gapToColumns && gapToColumns->get_inner()) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label +
                                            " declares a gap to the central column but is drawn running into it");
            }
            auto [above, below] = overlapAir(columns.main, "central");
            double depthOverlap = columnZOverlap(columns.main);
            if (depthOverlap <= 0) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " overlaps the central column in x but not in depth");
            }
            crossing.innerOverlapWidth = innerOverlapWidth;
            crossing.innerOverlapDepth = depthOverlap;
            if (std::fabs(above - below) > geometric_tolerance(above)) {
                throw NotImplementedException("Magnetic " + label + " held off-centre in the central column gap (unequal air above and below)");
            }
            crossing.innerOverlapAirPerSide = above;
        }
        else {
            double drawnGap = u0 - innerFace;
            if (!gapToColumns || !gapToColumns->get_inner()) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " stops " + std::to_string(drawnGap) +
                                            " m short of the central column but gapToColumns.inner is not given");
            }
            double declared = gapToColumns->get_inner().value();
            if (std::fabs(drawnGap - declared) > geometric_tolerance(declared)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " gapToColumns.inner is " + std::to_string(declared) +
                                            " m but the sheet is drawn " + std::to_string(drawnGap) + " m from the central column");
            }
            crossing.innerSideGap = declared;
        }

        // Outer (lateral column) contact.
        if (outerOverlapWidth > 0) {
            if (gapToColumns && gapToColumns->get_outer()) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label +
                                            " declares a gap to the lateral column but is drawn running into it");
            }
            auto [above, below] = overlapAir(lateral, "lateral");
            double depthOverlap = columnZOverlap(lateral);
            if (depthOverlap <= 0) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " overlaps a lateral column in x but not in depth");
            }
            crossing.outerOverlapWidth = outerOverlapWidth;
            crossing.outerOverlapDepth = depthOverlap;
            if (std::fabs(above - below) > geometric_tolerance(above)) {
                throw NotImplementedException("Magnetic " + label + " held off-centre in a lateral column gap (unequal air above and below)");
            }
            crossing.outerOverlapAirPerSide = above;
        }
        else {
            double drawnGap = outerFace - u1;
            if (!gapToColumns || !gapToColumns->get_outer()) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " stops " + std::to_string(drawnGap) +
                                            " m short of the lateral column but gapToColumns.outer is not given");
            }
            double declared = gapToColumns->get_outer().value();
            if (std::fabs(drawnGap - declared) > geometric_tolerance(declared)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " gapToColumns.outer is " + std::to_string(declared) +
                                            " m but the sheet is drawn " + std::to_string(drawnGap) + " m from the lateral column");
            }
            crossing.outerSideGap = declared;
        }

        if (shunt.get_segments() && !shunt.get_segments()->empty()) {
            if (crossing.innerOverlapWidth > 0 || crossing.outerOverlapWidth > 0) {
                throw NotImplementedException("Magnetic " + label + ": a segmented sheet running into a column gap");
            }
            double pieces = 0;
            double gaps = 0;
            auto segments = shunt.get_segments().value();
            for (auto& segment : segments) {
                if (!segment.get_length()) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " has a segment without a length");
                }
                pieces += segment.get_length().value();
                double gap = segment.get_gap() ? segment.get_gap().value() : 0.0;
                if (gap < 0) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " has a segment with a negative gap");
                }
                gaps += gap;
                if (gap > 0) {
                    crossing.segmentGaps.push_back(gap);
                }
            }
            if (std::fabs(pieces + gaps - width) > geometric_tolerance(width)) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " segments add up to " + std::to_string(pieces + gaps) +
                                            " m but the sheet is " + std::to_string(width) + " m wide");
            }
            // A sheet that stops short of both columns lies entirely inside one window.
            if (std::fabs(spanInWindow - width) > geometric_tolerance(width)) {
                throw NotImplementedException("Magnetic " + label + ": a segmented sheet spanning more than one window");
            }
        }

        inputs.crossings.push_back(crossing);
    }

    if (inputs.crossings.empty()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " is placed " +
                                    (shunt.get_placement() == MagneticShuntPlacement::IN_WINDOW ? std::string("inWindow") : std::string("betweenSections")) +
                                    " but does not cross a winding window");
    }
    return inputs;
}

double MagneticShuntModel::calculate_enclosed_magnetomotive_force(Magnetic& magnetic, const MagneticShuntNetworkInputs& inputs, size_t sourceIndex, size_t destinationIndex) {
    auto& coil = magnetic.get_mutable_coil();
    if (!coil.get_turns_description()) {
        throw CoilNotProcessedException("Magnetic shunt model needs the turns description");
    }
    size_t numberWindings = coil.get_functional_description().size();
    if (sourceIndex >= numberWindings || destinationIndex >= numberWindings || sourceIndex == destinationIndex) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic shunt model: invalid source/destination winding pair");
    }
    double sourceTurns = static_cast<double>(coil.get_number_turns(sourceIndex));
    double destinationTurns = static_cast<double>(coil.get_number_turns(destinationIndex));

    double sheetBottom = inputs.axialCoordinate - inputs.thickness / 2;
    double sheetTop = inputs.axialCoordinate + inputs.thickness / 2;

    double enclosedAbove = 0;
    double enclosedBelow = 0;
    auto turns = coil.get_turns_description().value();
    for (auto& turn : turns) {
        size_t windingIndex = coil.get_winding_index_by_name(turn.get_winding());
        if (!turn.get_dimensions()) {
            throw CoilNotProcessedException("Turn " + turn.get_name() + " has no dimensions");
        }
        auto dimensions = turn.get_dimensions().value();
        auto coordinates = turn.get_coordinates();
        double turnBottom = coordinates[1] - dimensions[1] / 2;
        double turnTop = coordinates[1] + dimensions[1] / 2;
        // One nanometre of contact is rounding of the drawn stack, not an intersection.
        const double contactTolerance = 1e-9;
        if (interval_overlap(turnBottom, turnTop, sheetBottom, sheetTop) > contactTolerance) {
            // Turns are described in the +x window; the mirrored window holds the same turn.
            double turnStart = std::fabs(coordinates[0]) - dimensions[0] / 2;
            double turnEnd = std::fabs(coordinates[0]) + dimensions[0] / 2;
            for (auto& crossing : inputs.crossings) {
                if (interval_overlap(turnStart, turnEnd, crossing.sheetStart, crossing.sheetEnd) > contactTolerance) {
                    throw InvalidInputException(ErrorCode::INVALID_INPUT, "Turn " + turn.get_name() + " intersects magnetic shunt '" + inputs.name + "'");
                }
            }
        }
        double current = 0;
        double parallels = static_cast<double>(coil.get_number_parallels(windingIndex));
        if (windingIndex == sourceIndex) {
            current = 1.0 / parallels;
        }
        else if (windingIndex == destinationIndex) {
            current = -(sourceTurns / destinationTurns) / parallels;
        }
        if (coordinates[1] > inputs.axialCoordinate) {
            enclosedAbove += current;
        }
        else {
            enclosedBelow += current;
        }
    }
    if (std::fabs(enclosedAbove + enclosedBelow) > 1e-9 * std::max(1.0, std::fabs(enclosedAbove))) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic shunt model: the leakage excitation is not ampere-turn balanced");
    }
    return enclosedAbove;
}

MagneticShuntContribution MagneticShuntModel::solve_network(const MagneticShuntNetworkInputs& inputs, double relativePermeability, double enclosedMagnetomotiveForcePerAmpere) {
    if (!(relativePermeability > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic shunt '" + inputs.name + "' needs a positive relative permeability");
    }
    double mu0 = Constants().vacuumPermeability;
    double thickness = inputs.thickness;
    double force = enclosedMagnetomotiveForcePerAmpere;

    MagneticShuntContribution contribution;
    contribution.name = inputs.name;
    contribution.relativePermeability = relativePermeability;
    contribution.enclosedMagnetomotiveForcePerAmpere = force;

    for (auto& crossing : inputs.crossings) {
        double depth = crossing.depthInWindow;
        double gapsTotal = 0;
        double gapsReluctance = 0;
        for (double gap : crossing.segmentGaps) {
            gapsTotal += gap;
            gapsReluctance += fringed_gap_reluctance(gap, thickness, depth);
        }
        // Zhang eq. (10): R_s2 = b_w / (mu0 mus h l_w), here over the material actually in the window.
        double sheetReluctance = (crossing.spanInWindow - gapsTotal) / (mu0 * relativePermeability * thickness * depth);
        double branch = sheetReluctance + gapsReluctance;
        if (crossing.innerSideGap) {
            branch += fringed_gap_reluctance(crossing.innerSideGap.value(), thickness, depth);
        }
        if (crossing.outerSideGap) {
            branch += fringed_gap_reluctance(crossing.outerSideGap.value(), thickness, depth);
        }

        double loop = inputs.coreLoopReluctance;
        double innerOverlapArea = crossing.innerOverlapWidth * crossing.innerOverlapDepth;
        double outerOverlapArea = crossing.outerOverlapWidth * crossing.outerOverlapDepth;
        if (innerOverlapArea > 0) {
            // Zhang R_s1 = h / (2 mu0 mus b_c l_w): half the sheet thickness across the overlap.
            loop += thickness / (2 * mu0 * relativePermeability * innerOverlapArea);
            loop += fringed_gap_reluctance(crossing.innerOverlapAirPerSide, crossing.innerOverlapWidth, crossing.innerOverlapDepth);
        }
        if (outerOverlapArea > 0) {
            loop += thickness / (2 * mu0 * relativePermeability * outerOverlapArea);
            loop += fringed_gap_reluctance(crossing.outerOverlapAirPerSide, crossing.outerOverlapWidth, crossing.outerOverlapDepth);
        }
        // Air held symmetrically: the loops above and below the sheet are equal.
        double upper = loop;
        double lower = loop;
        double total = branch + upper * lower / (upper + lower);

        MagneticShuntWindowCrossing result;
        result.side = crossing.side;
        result.spanInWindow = crossing.spanInWindow;
        result.depth = depth;
        result.windowWidth = crossing.windowWidth;
        result.shuntBranchReluctance = branch;
        result.upperLoopReluctance = upper;
        result.lowerLoopReluctance = lower;
        result.totalReluctance = total;
        result.innerOverlapArea = innerOverlapArea;
        result.outerOverlapArea = outerOverlapArea;
        contribution.crossings.push_back(result);

        contribution.networkLeakageInductance += force * force / total;
        double airField = force / crossing.windowWidth;
        contribution.displacedAirLeakageInductance += mu0 * airField * airField * thickness * crossing.spanInWindow * depth;

        double flux = force / total;
        double fluxDensity = std::fabs(flux) / (thickness * depth);
        contribution.magneticFluxDensityPerAmpere = std::max(contribution.magneticFluxDensityPerAmpere, fluxDensity);
        contribution.magneticFieldStrengthPerAmpere = std::max(contribution.magneticFieldStrengthPerAmpere,
                                                               fluxDensity / (mu0 * relativePermeability));
    }
    return contribution;
}

void MagneticShuntModel::add_losses(MagneticShuntContribution& contribution, const MagneticShuntNetworkInputs& inputs, const CoreMaterial& material, double frequency) {
    auto complex = material.get_permeability().get_complex();
    if (!complex) {
        return;
    }
    auto imaginary = complex->get_imaginary();
    if (!std::holds_alternative<std::vector<PermeabilityPoint>>(imaginary)) {
        return;
    }
    auto points = std::get<std::vector<PermeabilityPoint>>(imaginary);
    double minimumFrequency = DBL_MAX;
    double maximumFrequency = std::numeric_limits<double>::lowest();
    for (auto& point : points) {
        if (!point.get_frequency()) {
            throw InvalidInputException(ErrorCode::MISSING_DATA, "Complex permeability point without frequency in " + material.get_name());
        }
        minimumFrequency = std::min(minimumFrequency, point.get_frequency().value());
        maximumFrequency = std::max(maximumFrequency, point.get_frequency().value());
    }
    if (points.size() < 2 || frequency < minimumFrequency || frequency > maximumFrequency) {
        return;
    }
    double imaginaryPermeability = ComplexPermeability().get_complex_permeability(material, frequency).second;
    double mu0 = Constants().vacuumPermeability;
    double omega = 2 * std::numbers::pi * frequency;
    // Li eq. (27), over the in-window volume of the sheet. The parts lying in a column gap carry the
    // same flux spread over the overlap area (w d instead of t d), a field (t/w)^2 smaller in energy
    // density; they are not added.
    double losses = 0;
    for (size_t index = 0; index < contribution.crossings.size(); ++index) {
        auto& crossing = contribution.crossings[index];
        double flux = contribution.enclosedMagnetomotiveForcePerAmpere / crossing.totalReluctance;
        double field = flux / (inputs.thickness * crossing.depth) / (mu0 * contribution.relativePermeability);
        double volume = inputs.thickness * crossing.spanInWindow * crossing.depth;
        losses += 0.5 * omega * mu0 * imaginaryPermeability * field * field * volume;
    }
    contribution.lossesPerAmpereSquared = losses;
}

Core MagneticShuntModel::apply_shunts_to_gapping(Magnetic magnetic, double temperature, double frequency,
                                                 std::optional<std::vector<double>> relativePermeabilityPerShunt) {
    auto core = magnetic.get_core();
    if (!magnetic.get_shunts() || magnetic.get_shunts()->empty()) {
        return core;
    }
    check_supported_placements(magnetic);
    auto columns = resolve_e_columns(core);
    auto gapping = core.get_functional_description().get_gapping();
    auto shunts = magnetic.get_shunts().value();
    if (relativePermeabilityPerShunt && relativePermeabilityPerShunt->size() != shunts.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "apply_shunts_to_gapping: " + std::to_string(relativePermeabilityPerShunt->size()) +
                                    " relative permeabilities given for " + std::to_string(shunts.size()) + " shunts");
    }
    for (size_t shuntIndex = 0; shuntIndex < shunts.size(); ++shuntIndex) {
        auto& shunt = shunts[shuntIndex];
        auto label = shunt_label(shunt, shuntIndex);
        auto box = sheet_box(shunt, label);
        double thickness = box.y1 - box.y0;
        std::optional<double> relativePermeability;
        for (auto column : {columns.main, columns.right, columns.left}) {
            double columnX = column.get_coordinates()[0];
            double columnZ = column.get_coordinates().size() > 2 ? column.get_coordinates()[2] : 0.0;
            double overlapX = interval_overlap(box.x0, box.x1, columnX - column.get_width() / 2, columnX + column.get_width() / 2);
            double overlapZ = interval_overlap(box.z0, box.z1, columnZ - column.get_depth() / 2, columnZ + column.get_depth() / 2);
            if (overlapX <= 0 || overlapZ <= 0) {
                continue;
            }
            auto gap = gap_holding_sheet(core, column, box);
            if (!gap) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT, "Magnetic " + label + " overlaps a column at y = " +
                                            std::to_string((box.y0 + box.y1) / 2) + " m, where that column has no gap to hold it");
            }
            if (!relativePermeability) {
                relativePermeability = relativePermeabilityPerShunt ? relativePermeabilityPerShunt.value()[shuntIndex]
                                                                    : get_relative_permeability(resolve_material(shunt), temperature, frequency);
            }
            double coveredFraction = (overlapX * overlapZ) / (column.get_width() * column.get_depth());
            coveredFraction = std::min(1.0, coveredFraction);
            for (auto& storedGap : gapping) {
                if (!storedGap.get_coordinates()) {
                    continue;
                }
                auto storedCoordinates = storedGap.get_coordinates().value();
                auto gapCoordinates = gap->get_coordinates().value();
                bool same = std::fabs(storedCoordinates[0] - gapCoordinates[0]) < 1e-12 &&
                            std::fabs(storedCoordinates[1] - gapCoordinates[1]) < 1e-12 &&
                            std::fabs(storedGap.get_length() - gap->get_length()) < 1e-15;
                if (!same) {
                    continue;
                }
                double gapLength = storedGap.get_length();
                // The covered part is air (g - t) in series with the sheet (t / mus); the uncovered part
                // stays g; the two share the section in parallel.
                double coveredLength = gapLength - thickness + thickness / relativePermeability.value();
                double equivalentLength = coveredFraction >= 1.0
                                              ? coveredLength
                                              : 1.0 / (coveredFraction / coveredLength + (1.0 - coveredFraction) / gapLength);
                storedGap.set_length(equivalentLength);
            }
        }
    }
    core.get_mutable_functional_description().set_gapping(gapping);
    return core;
}

MagneticShuntLeakageResult MagneticShuntModel::assemble_leakage(Magnetic& magnetic, double windingLeakageInductance, double frequency,
                                                                size_t sourceIndex, size_t destinationIndex,
                                                                std::optional<std::vector<double>> relativePermeabilityPerShunt) {
    MagneticShuntLeakageResult result;
    result.windingLeakageInductance = windingLeakageInductance;
    result.leakageInductance = windingLeakageInductance;
    if (!magnetic.get_shunts()) {
        return result;
    }
    check_supported_placements(magnetic);
    double temperature = Defaults().ambientTemperature;
    auto shunts = magnetic.get_shunts().value();
    for (size_t shuntIndex = 0; shuntIndex < shunts.size(); ++shuntIndex) {
        if (!is_leakage_shunt(shunts[shuntIndex])) {
            continue;
        }
        auto inputs = extract_network_inputs(magnetic, shuntIndex, temperature, frequency);
        double force = calculate_enclosed_magnetomotive_force(magnetic, inputs, sourceIndex, destinationIndex);
        MagneticShuntContribution contribution;
        if (relativePermeabilityPerShunt) {
            contribution = solve_network(inputs, relativePermeabilityPerShunt.value()[shuntIndex], force);
        }
        else {
            auto material = resolve_material(shunts[shuntIndex]);
            contribution = solve_network(inputs, get_relative_permeability(material, temperature, frequency), force);
            add_losses(contribution, inputs, material, frequency);
        }
        result.leakageInductance += contribution.networkLeakageInductance - contribution.displacedAirLeakageInductance;
        result.shunts.push_back(contribution);
    }
    return result;
}

std::vector<MagneticShunt> MagneticShuntModel::size_shunt_for_leakage(Magnetic magnetic, double targetLeakageInductance, const CoreMaterial& material,
                                                                      double innerGapToColumn, double outerGapToColumn, double sourceCurrentPeak,
                                                                      double frequency, size_t sourceIndex, size_t destinationIndex) {
    if (magnetic.get_shunts() && !magnetic.get_shunts()->empty()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "size_shunt_for_leakage sizes sheets for a magnetic that has none yet");
    }
    if (!(targetLeakageInductance > 0) || innerGapToColumn < 0 || outerGapToColumn < 0 || !(sourceCurrentPeak > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "size_shunt_for_leakage needs a positive target, non-negative gaps and a positive current");
    }
    auto& coil = magnetic.get_mutable_coil();
    if (!coil.get_sections_description()) {
        throw CoilNotProcessedException("size_shunt_for_leakage needs the sections description");
    }
    auto sourceName = coil.get_functional_description().at(sourceIndex).get_name();
    auto destinationName = coil.get_functional_description().at(destinationIndex).get_name();
    auto sections = coil.get_sections_description().value();
    auto holds = [](const Section& section, const std::string& windingName) {
        for (auto& partialWinding : section.get_partial_windings()) {
            if (partialWinding.get_winding() == windingName) {
                return true;
            }
        }
        return false;
    };

    // The facing pair of source/destination conduction sections with the smallest axial clearance and
    // no conduction section in between.
    double bestClearance = DBL_MAX;
    double interfaceCentre = 0;
    for (auto& upper : sections) {
        if (upper.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        for (auto& lower : sections) {
            if (lower.get_type() != ElectricalType::CONDUCTION) {
                continue;
            }
            bool pair = (holds(upper, sourceName) && holds(lower, destinationName)) ||
                        (holds(upper, destinationName) && holds(lower, sourceName));
            if (!pair) {
                continue;
            }
            double upperBottom = upper.get_coordinates()[1] - upper.get_dimensions()[1] / 2;
            double lowerTop = lower.get_coordinates()[1] + lower.get_dimensions()[1] / 2;
            double clearance = upperBottom - lowerTop;
            double xOverlap = interval_overlap(upper.get_coordinates()[0] - upper.get_dimensions()[0] / 2, upper.get_coordinates()[0] + upper.get_dimensions()[0] / 2,
                                               lower.get_coordinates()[0] - lower.get_dimensions()[0] / 2, lower.get_coordinates()[0] + lower.get_dimensions()[0] / 2);
            if (clearance <= 0 || xOverlap <= 0 || clearance >= bestClearance) {
                continue;
            }
            bool blocked = false;
            for (auto& other : sections) {
                if (other.get_type() != ElectricalType::CONDUCTION) {
                    continue;
                }
                double otherBottom = other.get_coordinates()[1] - other.get_dimensions()[1] / 2;
                double otherTop = other.get_coordinates()[1] + other.get_dimensions()[1] / 2;
                if (interval_overlap(otherBottom, otherTop, lowerTop, upperBottom) > 0) {
                    blocked = true;
                    break;
                }
            }
            if (blocked) {
                continue;
            }
            bestClearance = clearance;
            interfaceCentre = (upperBottom + lowerTop) / 2;
        }
    }
    if (bestClearance == DBL_MAX) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "size_shunt_for_leakage: windings '" + sourceName + "' and '" + destinationName +
                                    "' have no facing sections stacked along the column axis with clearance between them");
    }

    auto& core = magnetic.get_mutable_core();
    auto columns = resolve_e_columns(core);
    double innerFace = columns.main.get_width() / 2;
    double outerFace = columns.right.get_coordinates()[0] - columns.right.get_width() / 2;
    double sheetStart = innerFace + innerGapToColumn;
    double sheetEnd = outerFace - outerGapToColumn;
    if (sheetEnd <= sheetStart) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "size_shunt_for_leakage: the gaps to the columns leave no sheet in the window");
    }
    double depth = columns.main.get_depth();

    auto make_sheets = [&](double thickness) {
        std::vector<MagneticShunt> sheets;
        for (int side : {1, -1}) {
            MagneticShunt sheet;
            sheet.set_name(std::string("Shunt ") + (side > 0 ? "+x" : "-x"));
            sheet.set_placement(MagneticShuntPlacement::BETWEEN_SECTIONS);
            sheet.set_coordinates({side * (sheetStart + sheetEnd) / 2, interfaceCentre, 0});
            sheet.set_dimensions({sheetEnd - sheetStart, thickness, depth});
            MagneticShuntGapToColumns gaps;
            gaps.set_inner(innerGapToColumn);
            gaps.set_outer(outerGapToColumn);
            sheet.set_gap_to_columns(gaps);
            sheet.set_material(material);
            sheets.push_back(sheet);
        }
        return sheets;
    };

    Magnetic withoutShunts = magnetic;
    withoutShunts.set_shunts(std::nullopt);
    auto windingOutput = LeakageInductance().calculate_leakage_inductance(withoutShunts, frequency, sourceIndex, destinationIndex);
    if (windingOutput.get_method_used() != "Energy") {
        throw NotImplementedException("size_shunt_for_leakage on top of the " + windingOutput.get_method_used() + " leakage method");
    }
    double windingLeakageInductance = windingOutput.get_leakage_inductance_per_winding()[0].get_nominal().value();
    if (targetLeakageInductance <= windingLeakageInductance) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "size_shunt_for_leakage: target " + std::to_string(targetLeakageInductance) +
                                    " H is not above the leakage without a shunt (" + std::to_string(windingLeakageInductance) + " H)");
    }

    auto evaluate = [&](double thickness) {
        Magnetic trial = magnetic;
        trial.set_shunts(make_sheets(thickness));
        return assemble_leakage(trial, windingLeakageInductance, frequency, sourceIndex, destinationIndex);
    };

    auto largest = evaluate(bestClearance);
    if (largest.leakageInductance < targetLeakageInductance) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "size_shunt_for_leakage: even a sheet filling the " + std::to_string(bestClearance) +
                                    " m clearance gives " + std::to_string(largest.leakageInductance) + " H, below the target " +
                                    std::to_string(targetLeakageInductance) + " H");
    }
    double low = 0;
    double high = bestClearance;
    MagneticShuntLeakageResult solution = largest;
    for (size_t iteration = 0; iteration < 100; ++iteration) {
        double middle = (low + high) / 2;
        auto trial = evaluate(middle);
        if (trial.leakageInductance < targetLeakageInductance) {
            low = middle;
        }
        else {
            high = middle;
            solution = trial;
        }
        if ((high - low) < 1e-6 * bestClearance) {
            break;
        }
    }

    double saturation = Core::get_magnetic_flux_density_saturation(material, Defaults().ambientTemperature, false);
    for (auto& contribution : solution.shunts) {
        double peakFluxDensity = contribution.magneticFluxDensityPerAmpere * sourceCurrentPeak;
        if (!(peakFluxDensity < 0.7 * saturation)) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "size_shunt_for_leakage: the sized sheet reaches " + std::to_string(peakFluxDensity) +
                                        " T at " + std::to_string(sourceCurrentPeak) + " A, not below 0.7 of the " + material.get_name() +
                                        " saturation (" + std::to_string(saturation) + " T)");
        }
    }
    return make_sheets(high);
}

} // namespace OpenMagnetics
