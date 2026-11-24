#include "physical_models/CoreMagneticFluxDensity.h"
#include "physical_models/CoreLosses.h"
#include "processors/Inputs.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/ReluctanceNetwork.h"
#include "support/Exceptions.h"
#include "support/Utils.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <numbers>

namespace OpenMagnetics {

namespace {

using Point = std::array<double, 2>;
using Polygon = std::vector<Point>;

double polygon_area(const Polygon& polygon) {
    double area = 0;
    for (size_t index = 0; index < polygon.size(); ++index) {
        const auto& current = polygon[index];
        const auto& next = polygon[(index + 1) % polygon.size()];
        area += current[0] * next[1] - next[0] * current[1];
    }
    return std::abs(area) / 2;
}

// Sutherland-Hodgman against the half-plane keep(p) = normal . p <= offset.
Polygon clip_half_plane(const Polygon& polygon, Point normal, double offset) {
    Polygon clipped;
    auto inside = [&](const Point& point) { return normal[0] * point[0] + normal[1] * point[1] <= offset; };
    for (size_t index = 0; index < polygon.size(); ++index) {
        const auto& current = polygon[index];
        const auto& next = polygon[(index + 1) % polygon.size()];
        bool currentInside = inside(current);
        bool nextInside = inside(next);
        if (currentInside) {
            clipped.push_back(current);
        }
        if (currentInside != nextInside) {
            double currentValue = normal[0] * current[0] + normal[1] * current[1] - offset;
            double nextValue = normal[0] * next[0] + normal[1] * next[1] - offset;
            double fraction = currentValue / (currentValue - nextValue);
            clipped.push_back({current[0] + fraction * (next[0] - current[0]), current[1] + fraction * (next[1] - current[1])});
        }
    }
    return clipped;
}

bool point_in_polygon(const Polygon& polygon, double x, double y) {
    bool inside = false;
    for (size_t index = 0, previous = polygon.size() - 1; index < polygon.size(); previous = index++) {
        const auto& current = polygon[index];
        const auto& other = polygon[previous];
        if (((current[1] > y) != (other[1] > y)) &&
            (x < (other[0] - current[0]) * (y - current[1]) / (other[1] - current[1]) + current[0])) {
            inside = !inside;
        }
    }
    return inside;
}

// The Schwarz-Christoffel L channel of right_angle_bend_squares (CorePiece.cpp, ABT #1766): limb of width a along
// -y (x in [-a, 0]), limb of width b along +x (y in [0, b]), window corner at the origin, core corner at (-a, b).
// A point of the bend is named by w = ln(zeta) = u + i theta, theta in (0, pi): theta = 0 is the window wall,
// theta = pi the outer wall; u -> -inf is the far end of limb a, u -> +inf that of limb b.
struct SchwarzChristoffelBend {
    double a;
    double b;
    double q;
    SchwarzChristoffelBend(double widthA, double widthB) : a(widthA), b(widthB), q((widthA / widthB) * (widthA / widthB)) {}

    std::complex<double> s_of(double u, double theta) const {
        std::complex<double> zeta = std::exp(std::complex<double>(u, theta));
        return std::sqrt((zeta - q) / (zeta + 1.0));
    }
    Point z_of(double u, double theta) const {
        auto s = s_of(u, theta);
        auto z = (2.0 / std::numbers::pi) * (b * std::atanh(s) - a * std::atan(b * s / a));
        return {z.real(), z.imag()};
    }
    // |dz/dw|: the length one unit of u (or theta) spans at that point.
    double scale(double u, double theta) const { return b / std::numbers::pi * std::abs(s_of(u, theta)); }
};

}  // namespace

CoreMagneticFluxDensity::CoreMagneticFluxDensity(const Magnetic& magnetic, double temperature, std::optional<double> frequency,
                                                 size_t numberRowsPerBend)
    : _numberRowsPerBend(numberRowsPerBend) {
    if (numberRowsPerBend < 2) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "The core flux density map needs at least 2 rows per bend");
    }
    auto core = magnetic.get_core();
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("needed to map its flux density");
    }
    if (ReluctanceNetwork::has_non_main_placement(magnetic)) {
        throw NotImplementedException("Core flux density map: a winding placed on a column other than the main one is not "
                                      "covered yet; the map drives the flux through the main column");
    }
    _curve = InitialPermeability::get_magnetisation_curve(core.resolve_material(), temperature, frequency);

    auto coreType = core.get_functional_description().get_type();
    if (coreType == CoreType::TOROIDAL) {
        build_toroid(core);
        find_region_twins();
        return;
    }
    auto family = core.get_shape_family();
    bool rectangularLimbFamily = family == CoreShapeFamily::E || family == CoreShapeFamily::PLANAR_E ||
                                 family == CoreShapeFamily::EF || family == CoreShapeFamily::U;
    if (coreType != CoreType::TWO_PIECE_SET || !rectangularLimbFamily) {
        throw NotImplementedException("Core flux density map: only two-piece sets of E, planar E, EF and U shapes and toroids are "
                                      "covered yet, not this " + std::string(magic_enum::enum_name(family)) + " core");
    }
    build_two_piece_set(core);
    find_region_twins();
}

void CoreMagneticFluxDensity::find_region_twins() {
    auto same = [](const Region& left, const Region& right) {
        if (left.depth != right.depth || left.fluxShare != right.fluxShare || left.slices.size() != right.slices.size()) {
            return false;
        }
        for (size_t sliceIndex = 0; sliceIndex < left.slices.size(); ++sliceIndex) {
            const auto& leftPaths = left.slices[sliceIndex].paths;
            const auto& rightPaths = right.slices[sliceIndex].paths;
            if (leftPaths.size() != rightPaths.size()) {
                return false;
            }
            for (size_t pathIndex = 0; pathIndex < leftPaths.size(); ++pathIndex) {
                const auto& a = leftPaths[pathIndex];
                const auto& b = rightPaths[pathIndex];
                if (a.length != b.length || a.width != b.width || a.gapLength != b.gapLength || a.lengthFraction != b.lengthFraction ||
                    a.cellIndexes.empty() != b.cellIndexes.empty()) {
                    return false;
                }
            }
        }
        return true;
    };
    _regionTwins.resize(_regions.size());
    for (size_t regionIndex = 0; regionIndex < _regions.size(); ++regionIndex) {
        _regionTwins[regionIndex] = regionIndex;
        for (size_t candidate = 0; candidate < regionIndex; ++candidate) {
            if (_regionTwins[candidate] == candidate && same(_regions[candidate], _regions[regionIndex])) {
                _regionTwins[regionIndex] = candidate;
                break;
            }
        }
    }
}

void CoreMagneticFluxDensity::add_straight_region(const std::string& name, Point origin, Point alongFlux, double length,
                                                  Point acrossFlux, double width, double depth, double fluxShare) {
    if (!(length > 0) || !(width > 0)) {
        return;
    }
    CoreFluxDensityCell cell;
    cell.region = name;
    cell.polygon = {origin,
                    {origin[0] + alongFlux[0] * length, origin[1] + alongFlux[1] * length},
                    {origin[0] + alongFlux[0] * length + acrossFlux[0] * width, origin[1] + alongFlux[1] * length + acrossFlux[1] * width},
                    {origin[0] + acrossFlux[0] * width, origin[1] + acrossFlux[1] * width}};
    cell.volume = length * width * depth / fluxShare;
    _cells.push_back(cell);
    Slice slice;
    slice.paths.push_back(Path{{_cells.size() - 1}, length, width, 0});
    _regions.push_back(Region{name, {slice}, depth, fluxShare});
}

// The canonical geometry of a bend: every curvilinear square of the grid (row by row, column by column), its scale,
// its full area and what the cuts keep of it.
struct CoreMagneticFluxDensity::BendGeometry {
    struct Square {
        size_t column;
        size_t row;
        double scale;
        double fullArea;
        double clippedArea;
        Polygon clipped;
        bool kept;
    };
    size_t numberColumns;
    std::vector<Square> squares;
};

CoreMagneticFluxDensity::BendGeometry CoreMagneticFluxDensity::compute_bend_geometry(double widthA, double widthB, double cutA, double cutB,
                                                                                    size_t numberRows) {
    SchwarzChristoffelBend bend(widthA, widthB);
    const double step = std::numbers::pi / static_cast<double>(numberRows);
    // The walls themselves are branch cuts of the square root: sample just inside them.
    const double wallOffset = 1e-9;
    auto theta_at = [&](size_t index) {
        return std::clamp(static_cast<double>(index) * step, wallOffset, std::numbers::pi - wallOffset);
    };

    // u range: from beyond the cut on limb a to beyond the cut on limb b (cells past a cut are clipped away).
    const double centreU = std::log(bend.q) / 2;
    double lowestU = centreU;
    while (bend.z_of(lowestU, std::numbers::pi / 2)[1] > -(cutA + widthA)) {
        lowestU -= step;
    }
    double highestU = centreU;
    while (bend.z_of(highestU, std::numbers::pi / 2)[0] < cutB + widthB) {
        highestU += step;
    }
    const size_t numberColumns = static_cast<size_t>(std::ceil((highestU - lowestU) / step));
    const size_t pointsPerEdge = 3;

    BendGeometry geometry;
    geometry.numberColumns = numberColumns;
    for (size_t rowIndex = 0; rowIndex < numberRows; ++rowIndex) {
        double thetaLow = theta_at(rowIndex);
        double thetaHigh = theta_at(rowIndex + 1);
        double thetaCentre = (static_cast<double>(rowIndex) + 0.5) * step;
        for (size_t columnIndex = 0; columnIndex < numberColumns; ++columnIndex) {
            double uLow = lowestU + static_cast<double>(columnIndex) * step;
            double uHigh = uLow + step;
            Polygon canonical;
            // The walls' two corners are singular points of the map: a cell whose wall edge spans one gets the
            // corner itself as a vertex (the window corner at zeta = q on theta = 0, the core's outer corner at
            // zeta = -1 on theta = pi), or its chord would cut the corner off.
            const double windowCornerU = std::log(bend.q);
            for (size_t point = 0; point < pointsPerEdge; ++point) {
                double u = uLow + step * point / pointsPerEdge;
                double next = uLow + step * (point + 1) / pointsPerEdge;
                canonical.push_back(bend.z_of(u, thetaLow));
                if (rowIndex == 0 && windowCornerU > u && windowCornerU < next) {
                    canonical.push_back({0, 0});
                }
            }
            for (size_t point = 0; point < pointsPerEdge; ++point) {
                canonical.push_back(bend.z_of(uHigh, thetaLow + (thetaHigh - thetaLow) * point / pointsPerEdge));
            }
            for (size_t point = 0; point < pointsPerEdge; ++point) {
                double u = uHigh - step * point / pointsPerEdge;
                double next = uHigh - step * (point + 1) / pointsPerEdge;
                canonical.push_back(bend.z_of(u, thetaHigh));
                if (rowIndex + 1 == numberRows && 0 < u && 0 > next) {
                    canonical.push_back({-widthA, widthB});
                }
            }
            for (size_t point = 0; point < pointsPerEdge; ++point) {
                canonical.push_back(bend.z_of(uLow, thetaHigh - (thetaHigh - thetaLow) * point / pointsPerEdge));
            }
            double fullArea = polygon_area(canonical);
            Polygon clipped = clip_half_plane(canonical, {0, -1}, cutA);  // y >= -cutA
            clipped = clip_half_plane(clipped, {1, 0}, cutB);              // x <= cutB
            double scale = bend.scale(uLow + step / 2, thetaCentre);
            BendGeometry::Square square{columnIndex, rowIndex, scale, fullArea, 0, {}, false};
            if (clipped.size() >= 3) {
                double clippedArea = polygon_area(clipped);
                if (clippedArea > 1e-9 * fullArea) {
                    square.clippedArea = clippedArea;
                    square.clipped = std::move(clipped);
                    square.kept = true;
                }
            }
            geometry.squares.push_back(std::move(square));
        }
    }
    return geometry;
}

void CoreMagneticFluxDensity::add_bend_region(const std::string& name, Point corner, bool mirrorX, bool mirrorY, double widthA,
                                              double widthB, double cutA, double cutB, double depth, double fluxShare) {
    if (!(widthA > 0) || !(widthB > 0) || !(cutA >= 0) || !(cutB >= 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                                    "Core flux density map: a bend needs positive limb widths, got " + std::to_string(widthA) + " and " +
                                        std::to_string(widthB));
    }
    const size_t numberRows = _numberRowsPerBend;
    const std::array<double, 4> geometryKey{widthA, widthB, cutA, cutB};
    auto cached = _bendGeometries.find(geometryKey);
    if (cached == _bendGeometries.end()) {
        cached = _bendGeometries.emplace(geometryKey, std::make_shared<const BendGeometry>(
                                                          compute_bend_geometry(widthA, widthB, cutA, cutB, numberRows))).first;
    }
    const BendGeometry& geometry = *cached->second;
    const size_t numberColumns = geometry.numberColumns;
    const double step = std::numbers::pi / static_cast<double>(numberRows);

    // pathsByColumn[column][row]: the bend's curvilinear square there, if the cuts leave any of it.
    std::vector<std::vector<std::optional<Path>>> pathsByColumn(numberColumns, std::vector<std::optional<Path>>(numberRows));
    for (const auto& square : geometry.squares) {
        // The whole curvilinear square is the path; the cut keeps a fraction of it.
        pathsByColumn[square.column][square.row] = Path{{}, square.scale * step, square.scale * step, 0, 0};
        if (!square.kept) {
            continue;
        }
        CoreFluxDensityCell cell;
        cell.region = name;
        for (const auto& point : square.clipped) {
            double x = mirrorX ? -point[0] : point[0];
            double y = mirrorY ? -point[1] : point[1];
            cell.polygon.push_back({corner[0] + x, corner[1] + y});
        }
        cell.volume = square.clippedArea * depth / fluxShare;
        _cells.push_back(cell);
        pathsByColumn[square.column][square.row]->cellIndexes.push_back(_cells.size() - 1);
        pathsByColumn[square.column][square.row]->lengthFraction = std::min(1.0, square.clippedArea / square.fullArea);
    }

    // Slices: one per column of the grid that keeps any cell. A column the end cuts clip is still solved whole, as the
    // band of the bend's field it is: its kept cells take their rows' B, and it adds to the region's MMF only the part
    // of the band the region keeps.
    Region region{name, {}, depth, fluxShare};
    for (size_t columnIndex = 0; columnIndex < numberColumns; ++columnIndex) {
        Slice slice;
        bool keepsCell = false;
        for (auto& path : pathsByColumn[columnIndex]) {
            slice.paths.push_back(path.value());
            keepsCell = keepsCell || !path->cellIndexes.empty();
        }
        if (keepsCell) {
            region.slices.push_back(slice);
        }
    }
    if (region.slices.empty()) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Core flux density map: " + name + " keeps no cell");
    }
    _regions.push_back(region);
}

void CoreMagneticFluxDensity::build_two_piece_set(Core core) {
    auto processedDescription = core.get_processed_description().value();
    auto family = core.get_shape_family();
    const bool uShape = family == CoreShapeFamily::U;
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    auto rightColumn = core.find_closest_column_by_coordinates({processedDescription.get_width() / 2, 0, -processedDescription.get_depth() / 2});
    auto columns = core.get_columns();
    if (uShape ? columns.size() != 2 : columns.size() != 3) {
        throw NotImplementedException("Core flux density map: a " + std::string(magic_enum::enum_name(family)) + " core with " +
                                      std::to_string(columns.size()) + " columns is not covered");
    }
    for (const auto& column : columns) {
        if (column.get_shape() != ColumnShape::RECTANGULAR) {
            throw NotImplementedException("Core flux density map: a core with a non-rectangular column is not covered yet");
        }
        if (column.get_depth() != mainColumn.get_depth() || column.get_height() != mainColumn.get_height()) {
            throw NotImplementedException("Core flux density map: columns of different depth or height are not covered");
        }
    }
    const double depth = mainColumn.get_depth();
    const double windowHeight = mainColumn.get_height();
    const double yokeThickness = (processedDescription.get_height() - windowHeight) / 2;
    if (!(yokeThickness > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: the core has no yoke above its window");
    }
    // The Painter's frame: a symmetric core is drawn from its main column's axis (x = 0) out to +x, the main column
    // as its half; a U is drawn whole from the main leg's outer face.
    const double mainLegWidth = uShape ? mainColumn.get_width() : mainColumn.get_width() / 2;
    const double mainLegInnerFace = mainColumn.get_width() / 2;
    const double showingCoreWidth = uShape ? processedDescription.get_width() - mainColumn.get_width() / 2
                                           : processedDescription.get_width() / 2;
    const double lateralLegWidth = rightColumn.get_width();
    const double lateralLegInnerFace = showingCoreWidth - lateralLegWidth;
    const double windowWidth = lateralLegInnerFace - mainLegInnerFace;
    if (!(windowWidth > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: the core has no winding window");
    }
    if (!uShape) {
        auto leftColumn = core.find_closest_column_by_coordinates({-processedDescription.get_width() / 2, 0, -processedDescription.get_depth() / 2});
        auto rightGaps = core.find_gaps_by_column(rightColumn);
        auto leftGaps = core.find_gaps_by_column(leftColumn);
        double rightGapLength = 0, leftGapLength = 0;
        for (const auto& gap : rightGaps) {
            rightGapLength += gap.get_length();
        }
        for (const auto& gap : leftGaps) {
            leftGapLength += gap.get_length();
        }
        if (leftColumn.get_width() != rightColumn.get_width() || leftGaps.size() != rightGaps.size() ||
            std::abs(leftGapLength - rightGapLength) > 1e-12) {
            throw NotImplementedException("Core flux density map: lateral legs of different width or gapping are not covered yet; "
                                          "the map takes the two halves of the core as mirror images");
        }
    }
    const double fluxShare = uShape ? 1.0 : 0.5;

    // Gaps as intervals of y along each leg, in the Painter's frame (an additive gap separates the two pieces and
    // shifts them by half its length each).
    double topOffset = 0, bottomOffset = 0;
    auto mainGaps = core.find_gaps_by_column(mainColumn);
    bool additive = !mainGaps.empty() && mainGaps.front().get_type() == GapType::ADDITIVE;
    if (additive) {
        topOffset = mainGaps.front().get_length() / 2;
        bottomOffset = -mainGaps.front().get_length() / 2;
    }
    auto gapIntervals = [&](const ColumnElement& column) {
        std::vector<std::pair<double, double>> intervals;
        if (additive) {
            intervals.push_back({bottomOffset, topOffset});
            return intervals;
        }
        for (const auto& gap : core.find_gaps_by_column(column)) {
            if (!gap.get_coordinates()) {
                throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: a gap has no coordinates");
            }
            double centre = gap.get_coordinates().value()[1];
            intervals.push_back({centre - gap.get_length() / 2, centre + gap.get_length() / 2});
        }
        std::sort(intervals.begin(), intervals.end());
        return intervals;
    };
    const double windowTop = topOffset + windowHeight / 2;
    const double windowBottom = bottomOffset - windowHeight / 2;
    const double coreTop = topOffset + processedDescription.get_height() / 2;
    const double coreBottom = bottomOffset - processedDescription.get_height() / 2;

    // A bend's straight arm reaches 1.5 limb widths (its field is uniform there to exp(-1.5 pi), 1 %) or less: to the
    // nearest gap face, or to the middle of the limb shared with the opposite bend.
    const double armReach = 1.5;

    // One leg: the bends at its two ends and the straight ferrite between them, cut by its gaps.
    struct LegCuts {
        double top;
        double bottom;
    };
    auto legCuts = [&](const std::vector<std::pair<double, double>>& gaps, double legWidth) {
        double topReach = armReach * legWidth;
        double bottomReach = armReach * legWidth;
        if (gaps.empty()) {
            topReach = std::min(topReach, (windowTop - windowBottom) / 2);
            bottomReach = std::min(bottomReach, (windowTop - windowBottom) / 2);
        }
        else {
            topReach = std::clamp(windowTop - gaps.back().second, 0.0, topReach);
            bottomReach = std::clamp(gaps.front().first - windowBottom, 0.0, bottomReach);
        }
        return LegCuts{topReach, bottomReach};
    };
    auto addLegStraights = [&](const std::string& name, const std::vector<std::pair<double, double>>& gaps, LegCuts cuts,
                               double xLeft, double width) {
        double segmentTop = windowTop - cuts.top;
        double lowest = windowBottom + cuts.bottom;
        std::vector<std::pair<double, double>> ferrite;
        double cursor = lowest;
        for (const auto& gap : gaps) {
            double gapBottom = std::max(gap.first, lowest);
            double gapTop = std::min(gap.second, segmentTop);
            if (gapTop <= gapBottom) {
                continue;
            }
            ferrite.push_back({cursor, gapBottom});
            cursor = gapTop;
        }
        ferrite.push_back({cursor, segmentTop});
        for (const auto& [bottom, top] : ferrite) {
            add_straight_region(name, {xLeft, bottom}, {0, 1}, top - bottom, {1, 0}, width, depth, fluxShare);
        }
    };

    auto mainGapIntervals = gapIntervals(mainColumn);
    auto lateralGapIntervals = gapIntervals(rightColumn);
    auto mainCuts = legCuts(mainGapIntervals, mainLegWidth);
    auto lateralCuts = legCuts(lateralGapIntervals, lateralLegWidth);
    double mainYokeCut = std::min(armReach * yokeThickness, windowWidth / 2);
    double lateralYokeCut = std::min(armReach * yokeThickness, windowWidth / 2);

    const std::string mainName = uShape ? "main leg" : "central column";
    const std::string lateralName = uShape ? "return leg" : "lateral column";
    add_bend_region(mainName + " / top yoke corner", {mainLegInnerFace, windowTop}, false, false, mainLegWidth, yokeThickness,
                    mainCuts.top, mainYokeCut, depth, fluxShare);
    add_bend_region(lateralName + " / top yoke corner", {lateralLegInnerFace, windowTop}, true, false, lateralLegWidth, yokeThickness,
                    lateralCuts.top, lateralYokeCut, depth, fluxShare);
    add_bend_region(mainName + " / bottom yoke corner", {mainLegInnerFace, windowBottom}, false, true, mainLegWidth, yokeThickness,
                    mainCuts.bottom, mainYokeCut, depth, fluxShare);
    add_bend_region(lateralName + " / bottom yoke corner", {lateralLegInnerFace, windowBottom}, true, true, lateralLegWidth,
                    yokeThickness, lateralCuts.bottom, lateralYokeCut, depth, fluxShare);

    addLegStraights(mainName, mainGapIntervals, mainCuts, mainLegInnerFace - mainLegWidth, mainLegWidth);
    addLegStraights(lateralName, lateralGapIntervals, lateralCuts, lateralLegInnerFace, lateralLegWidth);
    double yokeLength = windowWidth - mainYokeCut - lateralYokeCut;
    add_straight_region("top yoke", {mainLegInnerFace + mainYokeCut, windowTop}, {1, 0}, yokeLength, {0, 1}, yokeThickness, depth,
                        fluxShare);
    add_straight_region("bottom yoke", {mainLegInnerFace + mainYokeCut, coreBottom}, {1, 0}, yokeLength, {0, 1}, yokeThickness,
                        depth, fluxShare);
    (void) coreTop;
}

void CoreMagneticFluxDensity::build_toroid(Core core) {
    auto processedDescription = core.get_processed_description().value();
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    const double outerRadius = processedDescription.get_width() / 2;
    const double radialThickness = mainColumn.get_width();
    const double innerRadius = outerRadius - radialThickness;
    const double depth = mainColumn.get_depth();
    if (!(innerRadius > 0) || !(depth > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: the toroid has no positive bore or height");
    }
    if (mainColumn.get_shape() != ColumnShape::RECTANGULAR) {
        throw NotImplementedException("Core flux density map: a toroid of non-rectangular section is not covered yet");
    }
    double gapLength = 0;
    for (const auto& gap : core.get_functional_description().get_gapping()) {
        gapLength += gap.get_length();
    }
    const size_t numberRings = _numberRowsPerBend;
    const size_t numberSegments = 96;
    const double ringWidth = radialThickness / static_cast<double>(numberRings);
    Region region{"ring", {Slice{}}, depth, 1.0};
    for (size_t ringIndex = 0; ringIndex < numberRings; ++ringIndex) {
        double inner = innerRadius + ringWidth * static_cast<double>(ringIndex);
        double outer = inner + ringWidth;
        double centre = (inner + outer) / 2;
        CoreFluxDensityCell cell;
        cell.region = "ring";
        for (size_t segment = 0; segment <= numberSegments; ++segment) {
            // The closing vertex is the opening one exactly: sin(2 pi) is not 0 in floating point, and a seam vertex off the
            // x axis breaks the crossing count of a point on it.
            double angle = 2 * std::numbers::pi * static_cast<double>(segment % numberSegments) / numberSegments;
            cell.polygon.push_back({outer * std::cos(angle), outer * std::sin(angle)});
        }
        for (size_t segment = 0; segment <= numberSegments; ++segment) {
            double angle = -2 * std::numbers::pi * static_cast<double>(segment % numberSegments) / numberSegments;
            cell.polygon.push_back({inner * std::cos(angle), inner * std::sin(angle)});
        }
        double ringLength = 2 * std::numbers::pi * centre - gapLength;
        if (!(ringLength > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: the toroid's gap is longer than its bore");
        }
        // The exact ring (less its gap), not the drawn polygon's area.
        cell.volume = ringLength * ringWidth * depth;
        _cells.push_back(cell);
        region.slices[0].paths.push_back(Path{{_cells.size() - 1}, ringLength, ringWidth, gapLength});
    }
    _regions.push_back(region);
}

double CoreMagneticFluxDensity::field_strength(double fluxDensity) const {
    // MKF's own inverse of the magnetisation curve, with its policy past the curve's end.
    return InitialPermeability::get_magnetomotive_force_for_flux_density(*_curve, fluxDensity, 1.0, 1.0, 0.0).first;
}

double CoreMagneticFluxDensity::flux_density(double fieldStrength) const {
    if (_curve->paramagnetic) {
        return constants.vacuumPermeability * fieldStrength;
    }
    const auto& fieldStrengths = _curve->fieldStrength;
    const auto& fluxDensities = _curve->fluxDensity;
    if (fieldStrength >= fieldStrengths.back()) {
        // Past the table the curve continues on the vacuum's slope; field_strength applies MKF's policy on whether a
        // flux density there is accepted when the solved cell is read back.
        return fluxDensities.back() + constants.vacuumPermeability * (fieldStrength - fieldStrengths.back());
    }
    size_t upper = std::upper_bound(fieldStrengths.begin(), fieldStrengths.end(), fieldStrength) - fieldStrengths.begin();
    upper = std::max<size_t>(upper, 1);
    size_t lower = upper - 1;
    return fluxDensities[lower] + (fieldStrength - fieldStrengths[lower]) / (fieldStrengths[upper] - fieldStrengths[lower]) *
                                      (fluxDensities[upper] - fluxDensities[lower]);
}

double CoreMagneticFluxDensity::incremental_permeability(double fieldStrength) const {
    if (_curve->paramagnetic) {
        return constants.vacuumPermeability;
    }
    const auto& fieldStrengths = _curve->fieldStrength;
    const auto& fluxDensities = _curve->fluxDensity;
    if (fieldStrength >= fieldStrengths.back()) {
        return constants.vacuumPermeability;
    }
    size_t upper = std::upper_bound(fieldStrengths.begin(), fieldStrengths.end(), fieldStrength) - fieldStrengths.begin();
    upper = std::max<size_t>(upper, 1);
    size_t lower = upper - 1;
    return (fluxDensities[upper] - fluxDensities[lower]) / (fieldStrengths[upper] - fieldStrengths[lower]);
}

double CoreMagneticFluxDensity::path_flux_density(const Path& path, double magnetomotiveForce) const {
    if (path.gapLength == 0) {
        return flux_density(magnetomotiveForce / path.length);
    }
    // A gapped ring: H l + B g / mu0 = F, rising strictly with B; bisection on the ferrite's share of the MMF.
    double low = 0, high = magnetomotiveForce;
    for (size_t iteration = 0; iteration < 200; ++iteration) {
        double ferriteMagnetomotiveForce = (low + high) / 2;
        double fluxDensity = flux_density(ferriteMagnetomotiveForce / path.length);
        double total = ferriteMagnetomotiveForce + fluxDensity / constants.vacuumPermeability * path.gapLength;
        if (total > magnetomotiveForce) {
            high = ferriteMagnetomotiveForce;
        }
        else {
            low = ferriteMagnetomotiveForce;
        }
        if (high - low <= 1e-13 * magnetomotiveForce) {
            break;
        }
    }
    return flux_density((low + high) / 2 / path.length);
}

// The MMF across a slice whose parallel paths carry sliceFlux between them; the paths' B come back in
// pathFluxDensities. Each path's flux rises strictly with the shared MMF: bracketed Newton.
double CoreMagneticFluxDensity::solve_slice(const Slice& slice, double sliceFlux, double depth, std::vector<double>& pathFluxDensities) const {
    pathFluxDensities.assign(slice.paths.size(), 0);
    if (sliceFlux <= 0) {
        return 0;
    }
    auto total_flux = [&](double magnetomotiveForce, double* derivative) {
        double total = 0;
        if (derivative) {
            *derivative = 0;
        }
        for (size_t pathIndex = 0; pathIndex < slice.paths.size(); ++pathIndex) {
            const auto& path = slice.paths[pathIndex];
            double area = path.width * depth;
            double fluxDensity = path_flux_density(path, magnetomotiveForce);
            pathFluxDensities[pathIndex] = fluxDensity;
            total += fluxDensity * area;
            if (derivative) {
                // dPhi/dF of the path: area / (l / mu_inc + g / mu0).
                double fieldStrength = path.gapLength == 0 ? magnetomotiveForce / path.length
                                                           : (magnetomotiveForce - fluxDensity / constants.vacuumPermeability * path.gapLength) / path.length;
                *derivative += area / (path.length / incremental_permeability(fieldStrength) + path.gapLength / constants.vacuumPermeability);
            }
        }
        return total;
    };
    double derivative = 0;
    total_flux(0, &derivative);
    double magnetomotiveForce = sliceFlux / derivative;
    double low = 0;
    double high = magnetomotiveForce;
    while (total_flux(high, nullptr) < sliceFlux) {
        low = high;
        high *= 2;
        if (!std::isfinite(high)) {
            throw CalculationException(ErrorCode::CALCULATION_DIVERGED, "Core flux density map: a slice's MMF does not converge");
        }
    }
    magnetomotiveForce = std::clamp(magnetomotiveForce, low, high);
    for (size_t iteration = 0; iteration < 200; ++iteration) {
        double residual = total_flux(magnetomotiveForce, &derivative) - sliceFlux;
        if (std::abs(residual) <= 1e-10 * sliceFlux) {
            return magnetomotiveForce;
        }
        if (residual > 0) {
            high = magnetomotiveForce;
        }
        else {
            low = magnetomotiveForce;
        }
        double next = magnetomotiveForce - residual / derivative;
        if (!(next > low && next < high)) {
            next = (low + high) / 2;
        }
        magnetomotiveForce = next;
    }
    throw CalculationException(ErrorCode::CALCULATION_DIVERGED, "Core flux density map: a slice's MMF did not converge in 200 iterations");
}

CoreFluxDensityInstant CoreMagneticFluxDensity::calculate_flux_densities(double magneticFlux) const {
    if (!std::isfinite(magneticFlux)) {
        throw InvalidInputException(ErrorCode::CALCULATION_INVALID_INPUT, "Core flux density map: the flux is not a finite number");
    }
    CoreFluxDensityInstant instant;
    instant.magneticFlux = magneticFlux;
    instant.fluxDensities.assign(_cells.size(), 0);
    instant.relativePermeabilities.assign(_cells.size(), 0);
    instant.regionMagnetomotiveForces.assign(_regions.size(), 0);
    const double sign = magneticFlux < 0 ? -1.0 : 1.0;
    const double magnitude = std::abs(magneticFlux);
    const double initialRelativePermeability = incremental_permeability(0) / constants.vacuumPermeability;
    std::vector<double> pathFluxDensities;
    // Every region's solved slices (their MMF and paths' B), kept for the twins that reuse them.
    std::vector<std::vector<std::pair<double, std::vector<double>>>> solvedSlices(_regions.size());
    for (size_t regionIndex = 0; regionIndex < _regions.size(); ++regionIndex) {
        const auto& region = _regions[regionIndex];
        double regionFlux = region.fluxShare * magnitude;
        const size_t twin = _regionTwins[regionIndex];
        for (size_t sliceIndex = 0; sliceIndex < region.slices.size(); ++sliceIndex) {
            const auto& slice = region.slices[sliceIndex];
            double sliceMagnetomotiveForce;
            if (twin == regionIndex) {
                sliceMagnetomotiveForce = solve_slice(slice, regionFlux, region.depth, pathFluxDensities);
                solvedSlices[regionIndex].emplace_back(sliceMagnetomotiveForce, pathFluxDensities);
            }
            else {
                sliceMagnetomotiveForce = solvedSlices[twin][sliceIndex].first;
                pathFluxDensities = solvedSlices[twin][sliceIndex].second;
            }
            // The part of the slice the region keeps: its paths' kept fractions, weighted by the flux each carries.
            double keptFlux = 0;
            double totalFlux = 0;
            for (size_t pathIndex = 0; pathIndex < slice.paths.size(); ++pathIndex) {
                double pathFlux = pathFluxDensities[pathIndex] * slice.paths[pathIndex].width;
                keptFlux += pathFlux * slice.paths[pathIndex].lengthFraction;
                totalFlux += pathFlux;
            }
            if (totalFlux > 0) {
                instant.regionMagnetomotiveForces[regionIndex] += sliceMagnetomotiveForce * keptFlux / totalFlux;
            }
            for (size_t pathIndex = 0; pathIndex < slice.paths.size(); ++pathIndex) {
                double fluxDensity = pathFluxDensities[pathIndex];
                // Read back through MKF's inverse: it applies the policy for a flux density past the curve's end.
                double relativePermeability = fluxDensity > 0 ? fluxDensity / (constants.vacuumPermeability * field_strength(fluxDensity))
                                                              : initialRelativePermeability;
                for (size_t cellIndex : slice.paths[pathIndex].cellIndexes) {
                    instant.fluxDensities[cellIndex] = sign * fluxDensity;
                    instant.relativePermeabilities[cellIndex] = relativePermeability;
                }
            }
        }
    }
    return instant;
}

double CoreMagneticFluxDensity::calculate_flux_density_at_point(double x, double y, double magneticFlux) const {
    for (size_t cellIndex = 0; cellIndex < _cells.size(); ++cellIndex) {
        if (point_in_polygon(_cells[cellIndex].polygon, x, y)) {
            return std::abs(calculate_flux_densities(magneticFlux).fluxDensities[cellIndex]);
        }
    }
    throw InvalidInputException(ErrorCode::INVALID_INPUT, "Core flux density map: the point (" + std::to_string(x) + ", " +
                                                              std::to_string(y) + ") m is not in the core's ferrite");
}

CoreFluxDensityWaveform CoreMagneticFluxDensity::calculate_waveform(const Magnetic& magnetic, const OperatingPoint& operatingPoint,
                                                                   size_t numberRowsPerBend) {
    OperatingPoint workingOperatingPoint = operatingPoint;
    MagnetizingInductance magnetizingInductanceModel;
    auto fluxDensity = magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(magnetic, &workingOperatingPoint).second;
    if (!fluxDensity.get_waveform() || !fluxDensity.get_waveform()->get_time()) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                                   "Core flux density map: the magnetizing model returned a flux density without a sampled waveform");
    }
    const Waveform waveform = fluxDensity.get_waveform().value();
    std::vector<double> vertexData = waveform.get_data();
    std::vector<double> vertexTime = waveform.get_time().value();
    if (vertexData.empty() || vertexData.size() != vertexTime.size()) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Core flux density map: the flux density waveform is empty or has no time for every sample");
    }
    auto excitations = workingOperatingPoint.get_excitations_per_winding();
    if (excitations.empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Core flux density map: the operating point has no excitation");
    }
    const double frequency = excitations[0].get_frequency();
    const double period = 1 / frequency;
    // The map covers exactly one period, as the core-loss models read it. A waveform sampled over [0, T) is closed by
    // its periodic continuation, its first value again at T.
    double span = vertexTime.back() - vertexTime.front();
    if (span > period * (1 + 1e-9)) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Core flux density map: the flux density waveform covers " +
                                   std::to_string(span) + " s, more than the excitation's period " + std::to_string(period) + " s");
    }
    if (span < period * (1 - 1e-9)) {
        vertexTime.push_back(vertexTime.front() + period);
        vertexData.push_back(vertexData.front());
    }
    // A piecewise-linear waveform may come as its vertices only (a triangle is three points): the period is sampled
    // at 64 even instants plus every vertex, so the animation has frames and the true peak is still a sample.
    std::vector<double> time = vertexTime;
    const size_t evenSamples = 64;
    if (vertexTime.size() > 1 && vertexTime.back() > vertexTime.front()) {
        for (size_t sample = 0; sample < evenSamples; ++sample) {
            time.push_back(vertexTime.front() + (vertexTime.back() - vertexTime.front()) * static_cast<double>(sample) / evenSamples);
        }
        std::sort(time.begin(), time.end());
        time.erase(std::unique(time.begin(), time.end(), [&](double left, double right) {
                       return right - left <= 1e-12 * (vertexTime.back() - vertexTime.front());
                   }),
                   time.end());
    }
    std::vector<double> data;
    for (double instant : time) {
        size_t upper = std::lower_bound(vertexTime.begin(), vertexTime.end(), instant) - vertexTime.begin();
        if (upper == 0) {
            data.push_back(vertexData.front());
        }
        else if (upper >= vertexTime.size()) {
            data.push_back(vertexData.back());
        }
        else {
            size_t lower = upper - 1;
            double fraction = (instant - vertexTime[lower]) / (vertexTime[upper] - vertexTime[lower]);
            data.push_back(vertexData[lower] + fraction * (vertexData[upper] - vertexData[lower]));
        }
    }
    double temperature = workingOperatingPoint.get_conditions().get_ambient_temperature();
    CoreMagneticFluxDensity map(magnetic, temperature, frequency, numberRowsPerBend);
    double effectiveArea = magnetic.get_core().get_effective_area();

    CoreFluxDensityWaveform result;
    result.cells = map.get_cells();
    result.time = time;
    result.saturationFluxDensity = map.get_saturation_flux_density();
    result.frequency = frequency;
    result.temperature = temperature;
    double largestFlux = -1;
    // The map is odd in the flux (calculate_flux_densities solves |flux| and applies its sign), so a sample whose |flux|
    // was already solved is that solution, sign applied: the same numbers, bit for bit.
    std::map<double, std::vector<double>> solvedByMagnitude;
    for (size_t sample = 0; sample < data.size(); ++sample) {
        double flux = data[sample] * effectiveArea;
        result.magneticFlux.push_back(flux);
        double magnitude = std::abs(flux);
        auto solved = solvedByMagnitude.find(magnitude);
        if (solved == solvedByMagnitude.end()) {
            solved = solvedByMagnitude.emplace(magnitude, map.calculate_flux_densities(magnitude).fluxDensities).first;
        }
        std::vector<double> fluxDensities = solved->second;
        if (flux < 0) {
            for (auto& fluxDensity : fluxDensities) {
                fluxDensity = -fluxDensity;
            }
        }
        result.fluxDensities.push_back(std::move(fluxDensities));
        if (std::abs(flux) > largestFlux) {
            largestFlux = std::abs(flux);
            result.peakSampleIndex = sample;
        }
    }
    size_t numberCells = map.get_cells().size();
    result.peakFluxDensities.resize(numberCells);
    result.alternatingFluxDensityAmplitudes.resize(numberCells);
    for (size_t cell = 0; cell < numberCells; ++cell) {
        double minimum = std::numeric_limits<double>::max();
        double maximum = std::numeric_limits<double>::lowest();
        for (const auto& sample : result.fluxDensities) {
            minimum = std::min(minimum, sample[cell]);
            maximum = std::max(maximum, sample[cell]);
        }
        result.peakFluxDensities[cell] = std::abs(result.fluxDensities[result.peakSampleIndex][cell]);
        result.alternatingFluxDensityAmplitudes[cell] = (maximum - minimum) / 2;
    }
    return result;
}

CoreLossDensityMap CoreMagneticFluxDensity::calculate_core_losses(const Magnetic& magnetic, const CoreFluxDensityWaveform& waveform) {
    size_t numberCells = waveform.cells.size();
    if (numberCells == 0 || waveform.fluxDensities.size() != waveform.time.size() || waveform.time.size() < 2) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Core loss map: the flux density waveform has no cells or no sampled period");
    }
    if (!(waveform.frequency > 0)) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Core loss map: the flux density waveform has no frequency");
    }
    auto core = magnetic.get_core();
    const auto material = core.resolve_material();
    CoreLosses coreLossesModel;
    CoreLossDensityMap result;
    result.volumetricLosses.resize(numberCells);
    result.coreLosses = 0;
    result.volume = 0;
    // Cells with the same B(t) (a core's mirror-image corners) are priced once: the same excitation gives the same losses.
    std::map<std::vector<double>, double> volumetricLossesByWaveform;
    for (size_t cell = 0; cell < numberCells; ++cell) {
        std::vector<double> data(waveform.time.size());
        for (size_t sample = 0; sample < waveform.time.size(); ++sample) {
            data[sample] = waveform.fluxDensities[sample][cell];
        }
        auto priced = volumetricLossesByWaveform.find(data);
        if (priced != volumetricLossesByWaveform.end()) {
            result.volumetricLosses[cell] = priced->second;
            result.coreLosses += priced->second * waveform.cells[cell].volume;
            result.volume += waveform.cells[cell].volume;
            continue;
        }
        Waveform cellWaveform;
        cellWaveform.set_time(waveform.time);
        cellWaveform.set_data(data);
        SignalDescriptor fluxDensity;
        fluxDensity.set_waveform(cellWaveform);
        fluxDensity.set_harmonics(Inputs::calculate_harmonics_data(cellWaveform, waveform.frequency));
        fluxDensity.set_processed(Inputs::calculate_processed_data(fluxDensity, cellWaveform, false));
        OperatingPointExcitation excitation;
        excitation.set_frequency(waveform.frequency);
        excitation.set_magnetic_flux_density(fluxDensity);
        double volumetricLosses = coreLossesModel.get_core_volumetric_losses(material, excitation, waveform.temperature);
        volumetricLossesByWaveform.emplace(std::move(data), volumetricLosses);
        result.volumetricLosses[cell] = volumetricLosses;
        result.coreLosses += volumetricLosses * waveform.cells[cell].volume;
        result.volume += waveform.cells[cell].volume;
    }
    return result;
}

}  // namespace OpenMagnetics
