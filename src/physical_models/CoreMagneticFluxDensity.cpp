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

namespace {
// The magnetic's core, once its windings are known to sit on the main column, the column the map drives the flux through.
Core main_column_core(const Magnetic& magnetic) {
    if (ReluctanceNetwork::has_non_main_placement(magnetic)) {
        throw NotImplementedException("Core flux density map: a winding placed on a column other than the main one is not "
                                      "covered yet; the map drives the flux through the main column");
    }
    return magnetic.get_core();
}
}  // namespace

CoreMagneticFluxDensity::CoreMagneticFluxDensity(const Magnetic& magnetic, double temperature, std::optional<double> frequency,
                                                 size_t numberRowsPerBend)
    : CoreMagneticFluxDensity(main_column_core(magnetic), temperature, frequency, numberRowsPerBend) {}

CoreMagneticFluxDensity::CoreMagneticFluxDensity(Core core, double temperature, std::optional<double> frequency, size_t numberRowsPerBend)
    : _numberRowsPerBend(numberRowsPerBend) {
    if (numberRowsPerBend < 2) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "The core flux density map needs at least 2 rows per bend");
    }
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("needed to map its flux density");
    }
    _temperature = temperature;
    _frequency = frequency;
    _materials.push_back(core.resolve_material());
    _curves.push_back(InitialPermeability::get_magnetisation_curve(_materials.front(), temperature, frequency));

    auto coreType = core.get_functional_description().get_type();
    const auto shapeFamily = core.get_shape_family();
    if (shapeFamily == CoreShapeFamily::DRUM || shapeFamily == CoreShapeFamily::ROD || shapeFamily == CoreShapeFamily::DRUM_RING ||
        shapeFamily == CoreShapeFamily::DRUM_SEMISHIELDED || shapeFamily == CoreShapeFamily::DRUM_PLATE) {
        build_drum_family(core);
        find_region_twins();
        return;
    }
    if (coreType == CoreType::TOROIDAL) {
        build_toroid(core);
        find_region_twins();
        return;
    }
    auto family = core.get_shape_family();
    if (coreType != CoreType::TWO_PIECE_SET && coreType != CoreType::PIECE_AND_PLATE &&
        !(coreType == CoreType::CLOSED_SHAPE && family == CoreShapeFamily::MOLDED)) {
        throw NotImplementedException("Core flux density map: a " + std::string(magic_enum::enum_name(coreType)) + " " +
                                      std::string(magic_enum::enum_name(family)) + " core is not covered yet");
    }
    if (closed_core_columns_rectangular(core)) {
        // Rectangular columns of one depth: the core's own depth everywhere, its exact 2D geometry.
        build_closed_core(core, core.get_columns().front().get_depth());
        find_region_twins();
        return;
    }
    calibrate_closed_core(core);
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
    cell.materialIndex = _currentMaterialIndex;
    _cells.push_back(cell);
    Slice slice;
    slice.paths.push_back(Path{{_cells.size() - 1}, length, width, 0});
    _regions.push_back(Region{name, {slice}, depth, fluxShare, _currentMaterialIndex});
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
    // Each step moves the limb's coordinate on by about width / numberRows; once zeta = e^u is lost beside q in double
    // precision it stops moving, so a cut that far down a limb is beyond the map's reach: throw, never loop.
    auto beyond_reach = [&](double cut, double width) {
        return CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                                    "Core flux density map: a bend cut " + std::to_string(cut) + " m down a limb " + std::to_string(width) +
                                        " m wide is beyond the conformal map's double-precision reach");
    };
    const double centreU = std::log(bend.q) / 2;
    double lowestU = centreU;
    double lowestY = bend.z_of(lowestU, std::numbers::pi / 2)[1];
    while (lowestY > -(cutA + widthA)) {
        lowestU -= step;
        double next = bend.z_of(lowestU, std::numbers::pi / 2)[1];
        if (!(next < lowestY)) {
            throw beyond_reach(cutA, widthA);
        }
        lowestY = next;
    }
    double highestU = centreU;
    double highestX = bend.z_of(highestU, std::numbers::pi / 2)[0];
    while (highestX < cutB + widthB) {
        highestU += step;
        double next = bend.z_of(highestU, std::numbers::pi / 2)[0];
        if (!(next > highestX)) {
            throw beyond_reach(cutB, widthB);
        }
        highestX = next;
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
                                              double widthB, double cutA, double cutB, double depth, double fluxShare, double limbDepth) {
    if (!(widthA > 0) || !(widthB > 0) || !(cutA >= 0) || !(cutB >= 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA,
                                    "Core flux density map: a bend needs positive limb widths, got " + std::to_string(widthA) + " and " +
                                        std::to_string(widthB));
    }
    if (!(limbDepth > 0) || !(depth > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: " + name + " has no positive depth");
    }
    const size_t numberRows = _numberRowsPerBend;
    // A limb of another depth than the yoke's (a round post's area spread over its drawn width) enters the conformal
    // bend at the yoke's depth with its area-equivalent width; its side of the bend is then stretched back onto the
    // drawn width. Far from the corner its B is exactly its flux over its area.
    const double canonicalWidthA = limbDepth == depth ? widthA : widthA * limbDepth / depth;
    const double stretchA = widthA / canonicalWidthA;
    const std::array<double, 4> geometryKey{canonicalWidthA, widthB, cutA, cutB};
    auto cached = _bendGeometries.find(geometryKey);
    if (cached == _bendGeometries.end()) {
        cached = _bendGeometries.emplace(geometryKey, std::make_shared<const BendGeometry>(
                                                          compute_bend_geometry(canonicalWidthA, widthB, cutA, cutB, numberRows))).first;
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
            double drawnX = point[0] < 0 && canonicalWidthA != widthA ? point[0] * stretchA : point[0];
            double x = mirrorX ? -drawnX : drawnX;
            double y = mirrorY ? -point[1] : point[1];
            cell.polygon.push_back({corner[0] + x, corner[1] + y});
        }
        cell.volume = square.clippedArea * depth / fluxShare;
        cell.materialIndex = _currentMaterialIndex;
        _cells.push_back(cell);
        pathsByColumn[square.column][square.row]->cellIndexes.push_back(_cells.size() - 1);
        pathsByColumn[square.column][square.row]->lengthFraction = std::min(1.0, square.clippedArea / square.fullArea);
    }

    // Slices: one per column of the grid that keeps any cell. A column the end cuts clip is still solved whole, as the
    // band of the bend's field it is: its kept cells take their rows' B, and it adds to the region's MMF only the part
    // of the band the region keeps.
    Region region{name, {}, depth, fluxShare, _currentMaterialIndex};
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

namespace {

// Families the Painter draws whole, from the main leg's outer face: one return leg carries the whole flux.
bool drawn_whole(CoreShapeFamily family) {
    return family == CoreShapeFamily::U || family == CoreShapeFamily::UR || family == CoreShapeFamily::UT ||
           family == CoreShapeFamily::C || family == CoreShapeFamily::UI;
}

}  // namespace

bool CoreMagneticFluxDensity::closed_core_columns_rectangular(Core core) {
    auto columns = core.get_columns();
    for (const auto& column : columns) {
        if (column.get_shape() != ColumnShape::RECTANGULAR || column.get_depth() != columns.front().get_depth()) {
            return false;
        }
    }
    return true;
}

void CoreMagneticFluxDensity::build_closed_core(Core core, double yokeDepth) {
    auto processedDescription = core.get_processed_description().value();
    auto family = core.get_shape_family();
    const bool whole = drawn_whole(family);
    const bool pieceAndPlate = core.get_functional_description().get_type() == CoreType::PIECE_AND_PLATE;
    auto mainColumn = core.find_closest_column_by_coordinates({0, 0, 0});
    auto rightColumn = core.find_closest_column_by_coordinates({processedDescription.get_width() / 2, 0, -processedDescription.get_depth() / 2});
    auto columns = core.get_columns();
    if (columns.size() < 2) {
        throw NotImplementedException("Core flux density map: a " + std::string(magic_enum::enum_name(family)) + " core with " +
                                      std::to_string(columns.size()) + " column has no return path to map");
    }
    for (const auto& column : columns) {
        if (column.get_height() != mainColumn.get_height()) {
            throw NotImplementedException("Core flux density map: columns of different height are not covered");
        }
        if (!(column.get_area() > 0) || !(column.get_width() > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: a column has no positive area or width");
        }
    }
    // Share of the wound column's flux the drawn part carries, and the return legs' area that carries it.
    const double fluxShare = whole ? 1.0 : 0.5;
    double lateralArea = 0;
    for (const auto& column : columns) {
        if (column.get_coordinates() != mainColumn.get_coordinates()) {
            lateralArea += column.get_area();
        }
    }
    lateralArea *= fluxShare;
    if (!whole && columns.size() == 3) {
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
        if (leftColumn.get_width() != rightColumn.get_width() || leftColumn.get_area() != rightColumn.get_area() ||
            leftGaps.size() != rightGaps.size() || std::abs(leftGapLength - rightGapLength) > 1e-12) {
            throw NotImplementedException("Core flux density map: lateral legs of different size or gapping are not covered yet; "
                                          "the map takes the two halves of the core as mirror images");
        }
    }

    // The Painter's frame (paint_two_piece_set_core, paint_piece_and_plate_core): a symmetric core from its main
    // column's axis out to +x, the main column as its half; a U-like core whole from the main leg's outer face. A leg
    // is drawn at its minimum width when it has one (a lateral leg wrapped around a round post).
    const double mainLegWidth = whole ? mainColumn.get_width() : mainColumn.get_width() / 2;
    const double mainLegInnerFace = mainColumn.get_width() / 2;
    const double showingCoreWidth = whole ? processedDescription.get_width() - mainColumn.get_width() / 2
                                          : processedDescription.get_width() / 2;
    double lateralLegWidth = rightColumn.get_minimum_width() ? rightColumn.get_minimum_width().value() : rightColumn.get_width();
    // A moulded body's shell is drawn from the cavity to the body's edge (paint_drum_family_core), not at its column's
    // area-equivalent wall; its area is the column's all the same.
    size_t postMaterial = 0, coverMaterial = 0, baseMaterial = 0;
    if (family == CoreShapeFamily::MOLDED) {
        auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
        lateralLegWidth = (dimensions["A"] - dimensions["E"]) / 2;
        auto regionMaterials = core.resolve_region_materials();
        if (regionMaterials.size() > 1) {
            // CorePieceMolded's regions: post, cover, base (the shell counted with the base).
            std::vector<size_t> indices;
            for (const auto& regionMaterial : regionMaterials) {
                if (!regionMaterial) {
                    throw NotImplementedException("Core flux density map: a moulded body with an air region is not covered yet");
                }
                indices.push_back(add_material(regionMaterial.value()));
            }
            postMaterial = indices[0];
            coverMaterial = indices.size() > 1 ? indices[1] : indices[0];
            baseMaterial = indices.size() > 2 ? indices[2] : coverMaterial;
        }
    }
    const double lateralLegInnerFace = showingCoreWidth - lateralLegWidth;
    const double windowWidth = lateralLegInnerFace - mainLegInnerFace;
    if (!(windowWidth > 0) || !(lateralLegWidth > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: the core has no winding window");
    }
    // Each leg's depth is its area over its drawn width: B there is its flux over its true area.
    // (Rectangular columns of one depth are drawn at their own width: that depth itself.)
    const bool rectangular = closed_core_columns_rectangular(core);
    const double mainDepth = rectangular ? yokeDepth : mainColumn.get_area() * fluxShare / mainLegWidth;
    const double lateralDepth = rectangular ? yokeDepth : lateralArea / lateralLegWidth;
    if (!(yokeDepth > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: the yoke has no positive depth");
    }

    // Gaps as intervals of y along each leg, in the Painter's frame (an additive gap separates the two pieces and
    // shifts them by half its length each).
    const double windowHeight = mainColumn.get_height();
    double topOffset = 0, bottomOffset = 0;
    bool additive = false;
    for (const auto& gap : core.get_functional_description().get_gapping()) {
        if (gap.get_type() == GapType::ADDITIVE) {
            additive = true;
            topOffset = std::max(topOffset, gap.get_length() / 2);
        }
    }
    bottomOffset = -topOffset;
    double legTip = windowHeight / 2;
    double topYokeThickness = (processedDescription.get_height() - windowHeight) / 2;
    double bottomYokeThickness = topYokeThickness;
    // A UT's bar has no letter of its own: MKF's UT (CorePiece get_shape_constants) and the Painter take both yokes
    // (B - D) / 2 around a centred window, the bar meeting the legs at the window's top, as for a two-piece set.
    if (pieceAndPlate && family != CoreShapeFamily::UT) {
        auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
        if (dimensions.find("B2") == dimensions.end() || !(dimensions["B2"] > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: the piece-and-plate core has no plate thickness B2");
        }
        legTip = mainColumn.get_coordinates()[1] + windowHeight / 2;
        topYokeThickness = dimensions["B2"];
        bottomYokeThickness = processedDescription.get_height() - dimensions["B2"] - windowHeight;
    }
    if (!(topYokeThickness > 0) || !(bottomYokeThickness > 0)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: the core has no yoke above or below its window");
    }
    auto gapIntervals = [&](const ColumnElement& column) {
        std::vector<std::pair<double, double>> intervals;
        if (additive) {
            // The spacer sits where the pieces meet: the two-piece set's mid-plane, or a piece's leg tips under its plate.
            const double split = pieceAndPlate ? legTip : 0.0;
            intervals.push_back({split + bottomOffset, split + topOffset});
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
    // A piece and its plate: the window runs from the piece's back to the leg tips, a spacer lifting the plate.
    const double windowTop = pieceAndPlate ? legTip + topOffset : topOffset + windowHeight / 2;
    const double windowBottom = pieceAndPlate ? legTip - windowHeight + bottomOffset : bottomOffset - windowHeight / 2;
    const double coreBottom = windowBottom - bottomYokeThickness;

    // A bend's straight arm reaches 1.5 limb widths (its field is uniform there to exp(-1.5 pi), 1 %) or less: to the
    // nearest gap face, or to the middle of the limb shared with the opposite bend.
    const double armReach = 1.5;
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
                               double xLeft, double width, double legDepth) {
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
            add_straight_region(name, {xLeft, bottom}, {0, 1}, top - bottom, {1, 0}, width, legDepth, fluxShare);
        }
    };

    auto mainGapIntervals = gapIntervals(mainColumn);
    auto lateralGapIntervals = gapIntervals(rightColumn);
    auto mainCuts = legCuts(mainGapIntervals, mainLegWidth);
    auto lateralCuts = legCuts(lateralGapIntervals, lateralLegWidth);
    double topYokeCut = std::min(armReach * topYokeThickness, windowWidth / 2);
    double bottomYokeCut = std::min(armReach * bottomYokeThickness, windowWidth / 2);

    const std::string mainName = whole ? "main leg" : "central column";
    const std::string lateralName = whole ? "return leg" : "lateral column";
    const std::string topName = pieceAndPlate ? "plate" : "top yoke";
    const std::string bottomName = pieceAndPlate ? "back" : "bottom yoke";
    _currentMaterialIndex = coverMaterial;
    add_bend_region(mainName + " / " + topName + " corner", {mainLegInnerFace, windowTop}, false, false, mainLegWidth, topYokeThickness,
                    mainCuts.top, topYokeCut, yokeDepth, fluxShare, mainDepth);
    add_bend_region(lateralName + " / " + topName + " corner", {lateralLegInnerFace, windowTop}, true, false, lateralLegWidth,
                    topYokeThickness, lateralCuts.top, topYokeCut, yokeDepth, fluxShare, lateralDepth);
    _currentMaterialIndex = baseMaterial;
    add_bend_region(mainName + " / " + bottomName + " corner", {mainLegInnerFace, windowBottom}, false, true, mainLegWidth,
                    bottomYokeThickness, mainCuts.bottom, bottomYokeCut, yokeDepth, fluxShare, mainDepth);
    add_bend_region(lateralName + " / " + bottomName + " corner", {lateralLegInnerFace, windowBottom}, true, true, lateralLegWidth,
                    bottomYokeThickness, lateralCuts.bottom, bottomYokeCut, yokeDepth, fluxShare, lateralDepth);

    _currentMaterialIndex = postMaterial;
    addLegStraights(mainName, mainGapIntervals, mainCuts, mainLegInnerFace - mainLegWidth, mainLegWidth, mainDepth);
    _currentMaterialIndex = baseMaterial;
    addLegStraights(lateralName, lateralGapIntervals, lateralCuts, lateralLegInnerFace, lateralLegWidth, lateralDepth);
    _currentMaterialIndex = coverMaterial;
    add_straight_region(topName, {mainLegInnerFace + topYokeCut, windowTop}, {1, 0}, windowWidth - topYokeCut - topYokeCut, {0, 1},
                        topYokeThickness, yokeDepth, fluxShare);
    _currentMaterialIndex = baseMaterial;
    add_straight_region(bottomName, {mainLegInnerFace + bottomYokeCut, coreBottom}, {1, 0}, windowWidth - bottomYokeCut - bottomYokeCut, {0, 1},
                        bottomYokeThickness, yokeDepth, fluxShare);
    _currentMaterialIndex = 0;
}

void CoreMagneticFluxDensity::calibrate_closed_core(Core core) {
    // The plate of a core with round or irregular posts has no 2D depth: it is the one whose linear map has MKF's own
    // reluctance, le / (mu0 mu Ae) (the IEC 60205 constants of the shape's CorePiece), so the map spreads the flux
    // through a core MKF already describes, not through another one.
    const double targetFactor = core.get_effective_length() / core.get_effective_area();
    const double probeFlux = 1e-12;
    auto loopFactor = [&](double yokeDepth) {
        _cells.clear();
        _regions.clear();
        _bendGeometries.clear();
        build_closed_core(core, yokeDepth);
        _regionTwins.assign(_regions.size(), 0);
        for (size_t index = 0; index < _regions.size(); ++index) {
            _regionTwins[index] = index;
        }
        auto instant = calculate_flux_densities(probeFlux);
        // The geometric factor sum(l / A): each region's MMF times its own grade's initial permeability.
        double factor = 0;
        for (size_t regionIndex = 0; regionIndex < _regions.size(); ++regionIndex) {
            factor += instant.regionMagnetomotiveForces[regionIndex] * incremental_permeability(*_curves[_regions[regionIndex].materialIndex], 0);
        }
        return factor / probeFlux;
    };
    // The factor falls strictly as the plate deepens: a bracket grown outwards from the core's own depth by factors of 2
    // (up to 50 either way), then the secant on log(depth), safeguarded by bisection.
    const double processedDepth = core.get_processed_description()->get_depth();
    double low = processedDepth;
    double high = processedDepth;
    double factorLow = loopFactor(low) - targetFactor;
    double factorHigh = factorLow;
    while (!(factorHigh < 0) && high < processedDepth * 50) {
        low = high;
        factorLow = factorHigh;
        high = std::min(high * 2, processedDepth * 50);
        factorHigh = loopFactor(high) - targetFactor;
    }
    while (!(factorLow > 0) && low > processedDepth / 50) {
        high = low;
        factorHigh = factorLow;
        low = std::max(low / 2, processedDepth / 50);
        factorLow = loopFactor(low) - targetFactor;
    }
    if (!(factorLow > 0) || !(factorHigh < 0)) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT,
                                   "Core flux density map: no plate depth between " + std::to_string(processedDepth / 50) + " and " +
                                       std::to_string(processedDepth * 50) +
                                       " m gives this core MKF's effective le / Ae; its legs alone already exceed or fall short of it");
    }
    double logLow = std::log(low), logHigh = std::log(high);
    for (size_t iteration = 0; iteration < 60; ++iteration) {
        double logNext = logLow - factorLow * (logHigh - logLow) / (factorHigh - factorLow);
        if (!(logNext > logLow && logNext < logHigh)) {
            logNext = (logLow + logHigh) / 2;
        }
        double factorNext = loopFactor(std::exp(logNext)) - targetFactor;
        if (std::abs(factorNext) <= 1e-9 * targetFactor) {
            find_region_twins();
            return;
        }
        if (factorNext > 0) {
            logLow = logNext;
            factorLow = factorNext;
        }
        else {
            logHigh = logNext;
            factorHigh = factorNext;
        }
    }
    throw CalculationException(ErrorCode::CALCULATION_DIVERGED, "Core flux density map: the plate depth did not converge");
}

size_t CoreMagneticFluxDensity::add_material(const CoreMaterial& material) {
    for (size_t index = 0; index < _materials.size(); ++index) {
        if (_materials[index].get_name() == material.get_name()) {
            return index;
        }
    }
    _materials.push_back(material);
    _curves.push_back(InitialPermeability::get_magnetisation_curve(material, _temperature, _frequency));
    return _materials.size() - 1;
}

void CoreMagneticFluxDensity::build_drum_family(Core core) {
    auto family = core.get_shape_family();
    auto dimensions = flatten_dimensions(core.resolve_shape().get_dimensions().value());
    auto letter = [&](const std::string& name) {
        auto found = dimensions.find(name);
        if (found == dimensions.end() || !(found->second > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: the " + std::string(magic_enum::enum_name(family)) +
                                                                          " shape has no positive dimension " + name);
        }
        return found->second;
    };
    const double pi = std::numbers::pi;
    // The drawn half section (x = radius, y = height) carries half the flux: a ring of the body at radius rho is
    // pi rho deep there, so B = (flux / 2) / (pi rho t) is the true radial flux density of a disc of thickness t.
    const double fluxShare = 0.5;
    _currentMaterialIndex = 0;

    if (family == CoreShapeFamily::ROD) {
        // A bare cylinder: the whole flux runs through it, uniformly across its section.
        double radius = letter("A") / 2;
        double length = letter("B");
        double area = core.get_columns()[0].get_area();
        add_straight_region("rod", {0, -length / 2}, {0, 1}, length, {1, 0}, radius, area * fluxShare / radius, fluxShare);
        return;
    }

    if (family == CoreShapeFamily::DRUM_PLATE) {
        // A rectangular drum closed by a lid: one loop, post (along x) -> flange -> lid -> flange, each limb its own
        // rectangular section. Drawn whole (x along the post axis, y up), the post at the bottom.
        double length = letter("A");
        double flangeThickness = letter("D");
        double height = letter("B");
        double depth = letter("C");
        double postHeight = letter("F");
        double postDepth = dimensions.count("G") && dimensions["G"] > 0 ? dimensions["G"] : depth;
        double lidThickness = letter("J");
        double lidDepth = dimensions.count("L") && dimensions["L"] > 0 ? dimensions["L"] : depth;
        double windowWidth = length - 2 * flangeThickness;
        double windowHeight = height - postHeight;
        if (!(windowWidth > 0) || !(windowHeight > 0)) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: the drum-and-plate shape has no winding window");
        }
        const double share = 1.0;
        const double armReach = 1.5;
        double legCut = std::min(armReach * flangeThickness, windowHeight / 2);
        double lidCut = std::min(armReach * lidThickness, windowWidth / 2);
        double postCut = std::min(armReach * postHeight, windowWidth / 2);
        double legCutBottom = std::min(armReach * flangeThickness, windowHeight / 2);
        add_bend_region("flange / lid corner", {flangeThickness, height}, false, false, flangeThickness, lidThickness, legCut, lidCut, lidDepth, share, depth);
        add_bend_region("flange / lid corner", {length - flangeThickness, height}, true, false, flangeThickness, lidThickness, legCut, lidCut, lidDepth, share, depth);
        add_bend_region("flange / post corner", {flangeThickness, postHeight}, false, true, flangeThickness, postHeight, legCutBottom, postCut, postDepth, share, depth);
        add_bend_region("flange / post corner", {length - flangeThickness, postHeight}, true, true, flangeThickness, postHeight, legCutBottom, postCut, postDepth, share, depth);
        double legStraight = windowHeight - legCut - legCutBottom;
        add_straight_region("flange", {0, postHeight + legCutBottom}, {0, 1}, legStraight, {1, 0}, flangeThickness, depth, share);
        add_straight_region("flange", {length - flangeThickness, postHeight + legCutBottom}, {0, 1}, legStraight, {1, 0}, flangeThickness, depth, share);
        add_straight_region("lid", {flangeThickness + lidCut, height}, {1, 0}, windowWidth - lidCut - lidCut, {0, 1}, lidThickness, lidDepth, share);
        add_straight_region("post", {flangeThickness + postCut, 0}, {1, 0}, windowWidth - postCut - postCut, {0, 1}, postHeight, postDepth, share);
        return;
    }

    // Drum, shielded drum, semi-shielded drum: a round post between two flange discs (paint_drum_family_core).
    const double flangeRadius = letter("A") / 2;
    const double halfHeight = letter("B") / 2;
    const double postRadius = letter("C") / 2;
    const double topFlange = letter("D");
    const double bottomFlange = letter("F");
    const double postArea = core.get_columns()[0].get_area();
    const double postDepth = postArea * fluxShare / postRadius;
    const double windowTop = halfHeight - topFlange;
    const double windowBottom = -halfHeight + bottomFlange;
    const double windowHeight = windowTop - windowBottom;
    if (!(windowHeight > 0) || !(flangeRadius > postRadius)) {
        throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: the drum has no winding groove");
    }
    // The flux crosses the post and spreads radially through each flange to its rim, the path MKF's drum pieces
    // describe (CorePieceDrum::get_shape_constants): the post uniform, each flange B = flux / (2 pi rho t).
    const double armReach = 1.5;
    const double groove = flangeRadius - postRadius;
    const size_t flangeSlices = 12;
    auto addFlange = [&](const std::string& name, double cornerY, bool below, double thickness) {
        double cut = std::min(armReach * thickness, groove);
        double postCut = std::min(armReach * postRadius, windowHeight / 2);
        add_bend_region("post / " + name + " corner", {postRadius, cornerY}, false, below, postRadius, thickness, postCut, cut,
                        pi * postRadius, fluxShare, postDepth);
        // Beyond the corner, radial rings in a geometric progression to the rim, each pi rho deep at its middle radius.
        double inner = postRadius + cut;
        if (!(flangeRadius > inner)) {
            return postCut;
        }
        double ratio = std::pow(flangeRadius / inner, 1.0 / static_cast<double>(flangeSlices));
        double y = below ? cornerY - thickness : cornerY;
        for (size_t slice = 0; slice < flangeSlices; ++slice) {
            double from = inner * std::pow(ratio, static_cast<double>(slice));
            double to = slice + 1 == flangeSlices ? flangeRadius : inner * std::pow(ratio, static_cast<double>(slice + 1));
            add_straight_region(name, {from, y}, {1, 0}, to - from, {0, 1}, thickness, pi * (from + to) / 2, fluxShare);
        }
        return postCut;
    };
    double topPostCut = addFlange("top flange", windowTop, false, topFlange);
    double bottomPostCut = addFlange("bottom flange", windowBottom, true, bottomFlange);
    add_straight_region("post", {0, windowBottom + bottomPostCut}, {0, 1}, windowHeight - topPostCut - bottomPostCut, {1, 0}, postRadius,
                        postDepth, fluxShare);

    if (family == CoreShapeFamily::DRUM_RING || family == CoreShapeFamily::DRUM_SEMISHIELDED) {
        // The return path: a shield ring (its own grade when it has one), or the magnetic-epoxy shell from the flange rim
        // out, carrying the flux back axially, uniform across its section. The air between flange and ring is the gap
        // MKF models; it is not ferrite and holds no cell.
        double innerRadius, outerRadius, ringHeight = letter("L");
        if (family == CoreShapeFamily::DRUM_RING) {
            innerRadius = letter("K") / 2;
            outerRadius = letter("J") / 2;
            auto ringMaterials = core.resolve_region_materials();
            if (ringMaterials.size() > 1 && ringMaterials[1]) {
                _currentMaterialIndex = add_material(ringMaterials[1].value());
            }
        }
        else {
            innerRadius = flangeRadius;
            outerRadius = letter("J") / 2;
            _currentMaterialIndex = add_material(core.resolve_semishielded_shell_material());
        }
        if (!(outerRadius > innerRadius)) {
            throw InvalidInputException(ErrorCode::INVALID_CORE_DATA, "Core flux density map: the drum's ring or shell has no thickness");
        }
        double ringArea = pi * (outerRadius * outerRadius - innerRadius * innerRadius);
        double ringWidth = outerRadius - innerRadius;
        add_straight_region(family == CoreShapeFamily::DRUM_RING ? "ring" : "shell", {innerRadius, -ringHeight / 2}, {0, 1}, ringHeight,
                            {1, 0}, ringWidth, ringArea * fluxShare / ringWidth, fluxShare);
        _currentMaterialIndex = 0;
    }
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
    Region region{"ring", {Slice{}}, depth, 1.0, _currentMaterialIndex};
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
        cell.materialIndex = _currentMaterialIndex;
        _cells.push_back(cell);
        region.slices[0].paths.push_back(Path{{_cells.size() - 1}, ringLength, ringWidth, gapLength});
    }
    _regions.push_back(region);
}

double CoreMagneticFluxDensity::field_strength(const InitialPermeability::MagnetisationCurve& curve, double fluxDensity) const {
    // MKF's own inverse of the magnetisation curve, with its policy past the curve's end.
    return InitialPermeability::get_magnetomotive_force_for_flux_density(curve, fluxDensity, 1.0, 1.0, 0.0).first;
}

double CoreMagneticFluxDensity::flux_density(const InitialPermeability::MagnetisationCurve& curve, double fieldStrength) const {
    if (curve.paramagnetic) {
        return constants.vacuumPermeability * fieldStrength;
    }
    const auto& fieldStrengths = curve.fieldStrength;
    const auto& fluxDensities = curve.fluxDensity;
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

double CoreMagneticFluxDensity::incremental_permeability(const InitialPermeability::MagnetisationCurve& curve, double fieldStrength) const {
    if (curve.paramagnetic) {
        return constants.vacuumPermeability;
    }
    const auto& fieldStrengths = curve.fieldStrength;
    const auto& fluxDensities = curve.fluxDensity;
    if (fieldStrength >= fieldStrengths.back()) {
        return constants.vacuumPermeability;
    }
    size_t upper = std::upper_bound(fieldStrengths.begin(), fieldStrengths.end(), fieldStrength) - fieldStrengths.begin();
    upper = std::max<size_t>(upper, 1);
    size_t lower = upper - 1;
    return (fluxDensities[upper] - fluxDensities[lower]) / (fieldStrengths[upper] - fieldStrengths[lower]);
}

double CoreMagneticFluxDensity::path_flux_density(const InitialPermeability::MagnetisationCurve& curve, const Path& path, double magnetomotiveForce) const {
    if (path.gapLength == 0) {
        return flux_density(curve, magnetomotiveForce / path.length);
    }
    // A gapped ring: H l + B g / mu0 = F, rising strictly with B; bisection on the ferrite's share of the MMF.
    double low = 0, high = magnetomotiveForce;
    for (size_t iteration = 0; iteration < 200; ++iteration) {
        double ferriteMagnetomotiveForce = (low + high) / 2;
        double fluxDensity = flux_density(curve, ferriteMagnetomotiveForce / path.length);
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
    return flux_density(curve, (low + high) / 2 / path.length);
}

// The MMF across a slice whose parallel paths carry sliceFlux between them; the paths' B come back in
// pathFluxDensities. Each path's flux rises strictly with the shared MMF: bracketed Newton.
double CoreMagneticFluxDensity::solve_slice(const InitialPermeability::MagnetisationCurve& curve, const Slice& slice, double sliceFlux, double depth, std::vector<double>& pathFluxDensities) const {
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
            double fluxDensity = path_flux_density(curve, path, magnetomotiveForce);
            pathFluxDensities[pathIndex] = fluxDensity;
            total += fluxDensity * area;
            if (derivative) {
                // dPhi/dF of the path: area / (l / mu_inc + g / mu0).
                double fieldStrength = path.gapLength == 0 ? magnetomotiveForce / path.length
                                                           : (magnetomotiveForce - fluxDensity / constants.vacuumPermeability * path.gapLength) / path.length;
                *derivative += area / (path.length / incremental_permeability(curve, fieldStrength) + path.gapLength / constants.vacuumPermeability);
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

    std::vector<double> pathFluxDensities;
    // Every region's solved slices (their MMF and paths' B), kept for the twins that reuse them.
    std::vector<std::vector<std::pair<double, std::vector<double>>>> solvedSlices(_regions.size());
    for (size_t regionIndex = 0; regionIndex < _regions.size(); ++regionIndex) {
        const auto& region = _regions[regionIndex];
        const auto& curve = *_curves[region.materialIndex];
        const double initialRelativePermeability = incremental_permeability(curve, 0) / constants.vacuumPermeability;
        double regionFlux = region.fluxShare * magnitude;
        const size_t twin = _regionTwins[regionIndex];
        for (size_t sliceIndex = 0; sliceIndex < region.slices.size(); ++sliceIndex) {
            const auto& slice = region.slices[sliceIndex];
            double sliceMagnetomotiveForce;
            if (twin == regionIndex) {
                sliceMagnetomotiveForce = solve_slice(curve, slice, regionFlux, region.depth, pathFluxDensities);
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
                double relativePermeability = fluxDensity > 0 ? fluxDensity / (constants.vacuumPermeability * field_strength(curve, fluxDensity))
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
    auto excitations = workingOperatingPoint.get_excitations_per_winding();
    if (excitations.empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Core flux density map: the operating point has no excitation");
    }
    const double frequency = excitations[0].get_frequency();
    double temperature = workingOperatingPoint.get_conditions().get_ambient_temperature();
    CoreMagneticFluxDensity map(magnetic, temperature, frequency, numberRowsPerBend);
    return calculate_waveform(map, magnetic.get_core(), fluxDensity.get_waveform().value(), frequency);
}

CoreFluxDensityWaveform CoreMagneticFluxDensity::calculate_waveform(const CoreMagneticFluxDensity& map, const Core& core,
                                                                   const Waveform& fluxDensityWaveform, double frequency) {
    if (!fluxDensityWaveform.get_time()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Core flux density map: the flux density waveform has no time");
    }
    if (!(frequency > 0)) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Core flux density map: the flux density waveform has no frequency");
    }
    std::vector<double> vertexData = fluxDensityWaveform.get_data();
    std::vector<double> vertexTime = fluxDensityWaveform.get_time().value();
    if (vertexData.empty() || vertexData.size() != vertexTime.size()) {
        throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Core flux density map: the flux density waveform is empty or has no time for every sample");
    }
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
    const double temperature = map._temperature;
    // B times the area MKF reads it over is the wound column's flux (the post's area for a family-model core, Ae otherwise).
    double effectiveArea = MagnetizingInductance::get_flux_density_reference_area(core);

    CoreFluxDensityWaveform result;
    result.cells = map.get_cells();
    result.time = time;
    result.saturationFluxDensity = map.get_saturation_flux_density();
    result.materials = map.get_materials();
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

CoreLossDensityMap CoreMagneticFluxDensity::calculate_core_losses(const Magnetic&, const CoreFluxDensityWaveform& waveform) {
    return calculate_core_losses(waveform);
}

CoreLossDensityMap CoreMagneticFluxDensity::calculate_core_losses(const CoreFluxDensityWaveform& waveform) {
    size_t numberCells = waveform.cells.size();
    if (numberCells == 0 || waveform.fluxDensities.size() != waveform.time.size() || waveform.time.size() < 2) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Core loss map: the flux density waveform has no cells or no sampled period");
    }
    if (!(waveform.frequency > 0)) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Core loss map: the flux density waveform has no frequency");
    }
    if (waveform.materials.empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Core loss map: the flux density waveform carries no materials");
    }
    CoreLosses coreLossesModel;
    CoreLossDensityMap result;
    result.volumetricLosses.resize(numberCells);
    result.coreLosses = 0;
    result.volume = 0;
    // Cells with the same B(t) (a core's mirror-image corners) are priced once: the same excitation gives the same losses.
    std::map<std::pair<size_t, std::vector<double>>, double> volumetricLossesByWaveform;
    for (size_t cell = 0; cell < numberCells; ++cell) {
        std::vector<double> data(waveform.time.size());
        for (size_t sample = 0; sample < waveform.time.size(); ++sample) {
            data[sample] = waveform.fluxDensities[sample][cell];
        }
        const size_t materialIndex = waveform.cells[cell].materialIndex;
        if (materialIndex >= waveform.materials.size()) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Core loss map: a cell's material is not among the waveform's");
        }
        const auto& material = waveform.materials[materialIndex];
        // A secondary grade with no loss data (a semi-shielded drum's magnetic epoxy, a moulded body's region) dissipates
        // nothing and is named in unpricedMaterials, as CoreLosses prices those cores; the core's own grade must have it.
        if (materialIndex > 0 && material.get_volumetric_losses().empty()) {
            if (std::find(result.unpricedMaterials.begin(), result.unpricedMaterials.end(), material.get_name()) == result.unpricedMaterials.end()) {
                result.unpricedMaterials.push_back(material.get_name());
            }
            result.volumetricLosses[cell] = 0;
            result.volume += waveform.cells[cell].volume;
            continue;
        }
        auto key = std::make_pair(materialIndex, data);
        auto priced = volumetricLossesByWaveform.find(key);
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
        volumetricLossesByWaveform.emplace(std::move(key), volumetricLosses);
        result.volumetricLosses[cell] = volumetricLosses;
        result.coreLosses += volumetricLosses * waveform.cells[cell].volume;
        result.volume += waveform.cells[cell].volume;
    }
    return result;
}

}  // namespace OpenMagnetics
