#include "support/Utils.h"
#include "constructive_models/Bobbin.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <streambuf>
#include <vector>
#include "spline.h"
#include "support/Exceptions.h"

// tk::spline bobbinFillingFactorInterpWidth;
// tk::spline bobbinFillingFactorInterpHeight;
// tk::spline bobbinWindingWindowProportionInterpWidth;
// tk::spline bobbinWindingWindowProportionInterpHeight;
// double minBobbinWidth;
// double maxBobbinWidth;
// double minBobbinHeight;
// double maxBobbinHeight;
// double minWindingWindowWidth;
// double maxWindingWindowWidth;
// double minWindingWindowHeight;
// double maxWindingWindowHeight;

namespace OpenMagnetics {

class BobbinEDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(OpenMagnetics::Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin.get_functional_description()->get_dimensions());
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_shape(ColumnShape::RECTANGULAR);
            processedDescription.set_column_thickness(dimensions["s1"]);
            processedDescription.set_wall_thickness(dimensions["s2"]);
            WindingWindowElement windingWindowElement;
            // ABT #107: coordinates[0] is the winding-window CENTER (consumers do
            // left = coords[0] - width/2). The window starts at the central-column
            // surface (inner edge = f/2 + column thickness s1) and spans `windowWidth`,
            // so the centre is innerEdge + windowWidth/2 — matching create_quick_bobbin.
            double windowWidth = (dimensions["e"] - dimensions["f"] - 2 * dimensions["s1"]) / 2;
            double innerEdge = dimensions["f"] / 2 + dimensions["s1"];
            std::vector<double> coordinates({innerEdge + windowWidth / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["l2"] - 2 * dimensions["s2"]);
            windingWindowElement.set_width(windowWidth);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            // ABT #685 (Alf, 2026-08-16): store the column width — innerEdge IS the
            // column surface half-extent including the bobbin wall, the same convention
            // create_quick_bobbin uses (core half-width + columnThickness). Every family
            // processor computed it and threw it away, and paint_bobbin (among others)
            // then died on the empty optional for every database bobbin.
            processedDescription.set_column_width(innerEdge);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinRmDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(OpenMagnetics::Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin.get_functional_description()->get_dimensions());
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_shape(ColumnShape::ROUND);
            processedDescription.set_column_thickness((dimensions["D2"] - dimensions["D3"]) / 2);
            processedDescription.set_wall_thickness(dimensions["H5"]);
            WindingWindowElement windingWindowElement;
            // ABT #107: coordinates[0] is the winding-window CENTER. Inner edge is the
            // column surface radius D2/2; window spans windowWidth outward.
            double windowWidth = (dimensions["D1"] - dimensions["D2"]) / 2;
            double innerEdge = dimensions["D2"] / 2;
            std::vector<double> coordinates({innerEdge + windowWidth / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["H2"] - dimensions["H4"] - dimensions["H5"]);
            windingWindowElement.set_width(windowWidth);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            // ABT #685 (Alf, 2026-08-16): store the column width — innerEdge IS the
            // column surface half-extent including the bobbin wall, the same convention
            // create_quick_bobbin uses (core half-width + columnThickness). Every family
            // processor computed it and threw it away, and paint_bobbin (among others)
            // then died on the empty optional for every database bobbin.
            processedDescription.set_column_width(innerEdge);
            // Round column: the depth equals the width (both are the tube radius).
            processedDescription.set_column_depth(innerEdge);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinEpDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(OpenMagnetics::Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin.get_functional_description()->get_dimensions());
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_shape(ColumnShape::ROUND);
            processedDescription.set_column_thickness((dimensions["d2"] - dimensions["d3"]) / 2);
            processedDescription.set_wall_thickness(dimensions["s"]);
            WindingWindowElement windingWindowElement;
            // ABT #107: coordinates[0] is the winding-window CENTER (inner edge d2/2 + half width).
            double windowWidth = (dimensions["d1"] - dimensions["d2"]) / 2;
            double innerEdge = dimensions["d2"] / 2;
            std::vector<double> coordinates({innerEdge + windowWidth / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["h"] - 2 * dimensions["s"]);
            windingWindowElement.set_width(windowWidth);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            // ABT #685 (Alf, 2026-08-16): store the column width — innerEdge IS the
            // column surface half-extent including the bobbin wall, the same convention
            // create_quick_bobbin uses (core half-width + columnThickness). Every family
            // processor computed it and threw it away, and paint_bobbin (among others)
            // then died on the empty optional for every database bobbin.
            processedDescription.set_column_width(innerEdge);
            // Round column: the depth equals the width (both are the tube radius).
            processedDescription.set_column_depth(innerEdge);

            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinEtdDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(OpenMagnetics::Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin.get_functional_description()->get_dimensions());
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_shape(ColumnShape::ROUND);
            processedDescription.set_column_thickness((dimensions["d2"] - dimensions["d3"]) / 2);
            processedDescription.set_wall_thickness((dimensions["h1"] - dimensions["h2"]) / 2);
            WindingWindowElement windingWindowElement;
            // ABT #107: coordinates[0] is the winding-window CENTER. This branch
            // previously stored the FULL diameter d2 (twice the inner-edge radius d2/2,
            // placing the window centre beyond its own outer edge). Inner edge is the
            // column surface radius d2/2; centre is d2/2 + windowWidth/2.
            double windowWidth = (dimensions["d1"] - dimensions["d2"]) / 2;
            double innerEdge = dimensions["d2"] / 2;
            std::vector<double> coordinates({innerEdge + windowWidth / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["h2"]);
            windingWindowElement.set_width(windowWidth);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            // ABT #685 (Alf, 2026-08-16): store the column width — innerEdge IS the
            // column surface half-extent including the bobbin wall, the same convention
            // create_quick_bobbin uses (core half-width + columnThickness). Every family
            // processor computed it and threw it away, and paint_bobbin (among others)
            // then died on the empty optional for every database bobbin.
            processedDescription.set_column_width(innerEdge);
            // Round column: the depth equals the width (both are the tube radius).
            processedDescription.set_column_depth(innerEdge);

            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinPmDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(OpenMagnetics::Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin.get_functional_description()->get_dimensions());
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_shape(ColumnShape::ROUND);
            processedDescription.set_column_thickness((dimensions["d2"] - dimensions["d3"]) / 2);
            processedDescription.set_wall_thickness(dimensions["s1"]);
            WindingWindowElement windingWindowElement;
            // ABT #107: coordinates[0] is the winding-window CENTER (inner edge d2/2 + half width).
            double windowWidth = (dimensions["d1"] - dimensions["d2"]) / 2;
            double innerEdge = dimensions["d2"] / 2;
            std::vector<double> coordinates({innerEdge + windowWidth / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["h"] - dimensions["s1"] - dimensions["s2"]);
            windingWindowElement.set_width(windowWidth);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            // ABT #685 (Alf, 2026-08-16): store the column width — innerEdge IS the
            // column surface half-extent including the bobbin wall, the same convention
            // create_quick_bobbin uses (core half-width + columnThickness). Every family
            // processor computed it and threw it away, and paint_bobbin (among others)
            // then died on the empty optional for every database bobbin.
            processedDescription.set_column_width(innerEdge);
            // Round column: the depth equals the width (both are the tube radius).
            processedDescription.set_column_depth(innerEdge);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinPqDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(OpenMagnetics::Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin.get_functional_description()->get_dimensions());
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_shape(ColumnShape::ROUND);
            processedDescription.set_column_thickness((dimensions["D2"] - dimensions["D3"]) / 2);
            processedDescription.set_wall_thickness((dimensions["H1"] - dimensions["H2"]) / 2);
            WindingWindowElement windingWindowElement;
            // ABT #107: coordinates[0] is the winding-window CENTER. Previously stored the
            // FULL diameter D2 (twice the inner-edge radius). Centre is D2/2 + windowWidth/2.
            double windowWidth = (dimensions["D1"] - dimensions["D2"]) / 2;
            double innerEdge = dimensions["D2"] / 2;
            std::vector<double> coordinates({innerEdge + windowWidth / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["H2"]);
            windingWindowElement.set_width(windowWidth);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            // ABT #685 (Alf, 2026-08-16): store the column width — innerEdge IS the
            // column surface half-extent including the bobbin wall, the same convention
            // create_quick_bobbin uses (core half-width + columnThickness). Every family
            // processor computed it and threw it away, and paint_bobbin (among others)
            // then died on the empty optional for every database bobbin.
            processedDescription.set_column_width(innerEdge);
            // Round column: the depth equals the width (both are the tube radius).
            processedDescription.set_column_depth(innerEdge);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinEcDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(OpenMagnetics::Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin.get_functional_description()->get_dimensions());
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_shape(ColumnShape::ROUND);
            processedDescription.set_column_thickness((dimensions["D2"] - dimensions["D3"]) / 2);
            processedDescription.set_wall_thickness((dimensions["H1"] - dimensions["H2"]) / 2);
            WindingWindowElement windingWindowElement;
            // ABT #107: coordinates[0] is the winding-window CENTER. Previously stored the
            // FULL diameter D2 (twice the inner-edge radius). Centre is D2/2 + windowWidth/2.
            double windowWidth = (dimensions["D1"] - dimensions["D2"]) / 2;
            double innerEdge = dimensions["D2"] / 2;
            std::vector<double> coordinates({innerEdge + windowWidth / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["H2"]);
            windingWindowElement.set_width(windowWidth);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            // ABT #685 (Alf, 2026-08-16): store the column width — innerEdge IS the
            // column surface half-extent including the bobbin wall, the same convention
            // create_quick_bobbin uses (core half-width + columnThickness). Every family
            // processor computed it and threw it away, and paint_bobbin (among others)
            // then died on the empty optional for every database bobbin.
            processedDescription.set_column_width(innerEdge);
            // Round column: the depth equals the width (both are the tube radius).
            processedDescription.set_column_depth(innerEdge);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinEfdDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(OpenMagnetics::Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin.get_functional_description()->get_dimensions());
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_shape(ColumnShape::RECTANGULAR);
            processedDescription.set_column_thickness(dimensions["S1"]);
            processedDescription.set_wall_thickness(dimensions["S2"]);
            WindingWindowElement windingWindowElement;
            // ABT #107: coordinates[0] is the winding-window CENTER (inner edge f1/2 + S1, + half width).
            double windowWidth = (dimensions["e"] - dimensions["f1"] - 2 * dimensions["S1"]) / 2;
            double innerEdge = dimensions["f1"] / 2 + dimensions["S1"];
            std::vector<double> coordinates({innerEdge + windowWidth / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["d"] - 2 * dimensions["S2"]);
            windingWindowElement.set_width(windowWidth);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            // ABT #685 (Alf, 2026-08-16): store the column width — innerEdge IS the
            // column surface half-extent including the bobbin wall, the same convention
            // create_quick_bobbin uses (core half-width + columnThickness). Every family
            // processor computed it and threw it away, and paint_bobbin (among others)
            // then died on the empty optional for every database bobbin.
            processedDescription.set_column_width(innerEdge);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinTDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(OpenMagnetics::Bobbin bobbin) {
            // Toroidal "virtual" bobbin: the winding is held directly on the core
            // ring with no physical former. Column and wall thicknesses are 0.
            // Dimensions A (outer diameter), B (inner diameter), C (height) match
            // the ring-core shape dimensions (see CorePieceT::process_winding_window).
            auto dimensions = flatten_dimensions(bobbin.get_functional_description()->get_dimensions());
            CoreBobbinProcessedDescription processedDescription;
            double columnWidth = (dimensions["A"] - dimensions["B"]) / 2;
            processedDescription.set_column_shape(ColumnShape::RECTANGULAR);
            processedDescription.set_column_thickness(0);
            processedDescription.set_wall_thickness(0);
            processedDescription.set_column_depth(dimensions["C"] / 2);
            processedDescription.set_column_width(columnWidth / 2);
            WindingWindowElement windingWindowElement;
            windingWindowElement.set_shape(WindingWindowShape::ROUND);
            windingWindowElement.set_radial_height(dimensions["B"] / 2);
            windingWindowElement.set_angle(360);
            windingWindowElement.set_area(std::numbers::pi * pow(dimensions["B"] / 2, 2));
            windingWindowElement.set_coordinates(std::vector<double>({dimensions["B"] / 2, 0, 0}));
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

std::shared_ptr<BobbinDataProcessor> BobbinDataProcessor::factory(Bobbin bobbin) {

    // ABT #763: same undefined behaviour ABT #631 removed from get_winding_window_shape,
    // still live here. get_functional_description() returns the optional BY VALUE; when it
    // is disengaged, `->get_family()` reads the optional's *uninitialised* storage as a
    // BobbinFamily. That is not an exception, it is UB: the byte that comes back is
    // whatever the stack happened to hold, so the same input decides its own fate at
    // runtime. Observed on PyOpenMagnetics 1.7.0 with a bobbin json that carries no
    // functionalDescription (e.g. a caller passing json.dumps(row) — a json *string* —
    // instead of the object): most runs land on a value matching no enumerator and report
    // the misleading "Unknown bobbin family"; some land on E/ER and go on to read the
    // equally uninitialised dimensions map, yielding either a SEGFAULT or a silent
    // 0.0 x 0.0 winding window. Ask the question before dereferencing.
    if (!bobbin.get_functional_description()) {
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
            "Bobbin has no functionalDescription, so its family is unknown and it cannot be "
            "processed. A bobbin to be processed must be a json OBJECT carrying "
            "functionalDescription.family and functionalDescription.dimensions.");
    }

    auto family = bobbin.get_functional_description()->get_family();
    if (family == BobbinFamily::E) {
        return std::make_shared<BobbinEDataProcessor>();
    }
    else if (family == BobbinFamily::RM) {
        return std::make_shared<BobbinRmDataProcessor>();
    }
    else if (family == BobbinFamily::EP) {
        return std::make_shared<BobbinEpDataProcessor>();
    }
    else if (family == BobbinFamily::ETD) {
        return std::make_shared<BobbinEtdDataProcessor>();
    }
    else if (family == BobbinFamily::PM) {
        return std::make_shared<BobbinPmDataProcessor>();
    }
    else if (family == BobbinFamily::PQ) {
        return std::make_shared<BobbinPqDataProcessor>();
    }
    else if (family == BobbinFamily::EC) {
        return std::make_shared<BobbinEcDataProcessor>();
    }
    else if (family == BobbinFamily::EFD) {
        return std::make_shared<BobbinEfdDataProcessor>();
    }
    else if (family == BobbinFamily::T) {
        return std::make_shared<BobbinTDataProcessor>();
    }
    // ER, EL share E-style geometry (round/elliptical centre column on a
    // rectangular winding window), so the BobbinEDataProcessor is the right
    // dimensions interpreter. P, U are different geometries but are treated
    // as E-like here as a non-blocking fallback so the adviser doesn't reject
    // entire core families.
    else if (family == BobbinFamily::ER ||
             family == BobbinFamily::EL ||
             family == BobbinFamily::P  ||
             family == BobbinFamily::U) {
        return std::make_shared<BobbinEDataProcessor>();
    }
    else
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
            "Unknown bobbin family (enumerator value " + std::to_string(static_cast<int>(family)) +
            "), available options are: {E, EC, EFD, EL, EP, ER, ETD, P, PM, PQ, RM, T, U}");
}

void load_interpolators() {
    if (bobbinDatabase.empty() ||
        bobbinFillingFactorInterpWidth.get_x().size() == 0 || 
        bobbinFillingFactorInterpHeight.get_x().size() == 0 || 
        bobbinWindingWindowProportionInterpWidth.get_x().size() == 0 || 
        bobbinWindingWindowProportionInterpHeight.get_x().size() == 0) {
        // ABT #113: only load when the catalog is actually missing. This used
        // to call load_bobbins() unconditionally whenever an interpolator was
        // cold, RELOADING the shared bobbin catalog — a mutation, which throws
        // while databases are frozen (each thread's interpolators are
        // thread_local and start cold; the shared catalog stays loaded).
        if (bobbinDatabase.empty()) {
            load_bobbins();
        }

        struct AuxFillingFactorWidth
        {
            double windingWindowWidth, fillingFactor;
        };
        struct AuxFillingFactorHeight
        {
            double windingWindowHeight, fillingFactor;
        };
        std::vector<AuxFillingFactorWidth> auxFillingFactorWidth;
        std::vector<AuxFillingFactorHeight> auxFillingFactorHeight;

        struct AuxWindingWindowWidth
        {
            double windingWindowWidth, windingWindowProportion;
        };
        struct AuxWindingWindowHeight
        {
            double windingWindowHeight, windingWindowProportion;
        };
        std::vector<AuxWindingWindowWidth> auxWindingWindowWidth;
        std::vector<AuxWindingWindowHeight> auxWindingWindowHeight;


        minBobbinWallThickness = std::numeric_limits<double>::infinity();
        minBobbinColumnThickness = std::numeric_limits<double>::infinity();

        // ABT #631: rows this scan cannot use, reported once at the end instead of
        // vanishing into a bare `continue` — 34 of MAS's 504 bobbins point at core
        // shapes MAS does not ship (10 "EI …", 24 "M …"), and nobody knew because the
        // skip was silent.
        std::vector<std::string> unusableBobbins;

        for (auto& datum : bobbinDatabase) {
            // ABT #631: these skips must NOT be spelled as throw-and-catch. Wherever the
            // engine is compiled with exception catching disabled — the Emscripten
            // default, and how MVB++'s WASM module builds MKF (not one of its 92 objects
            // carries a __cxa_begin_catch) — the catch below is deleted and the first
            // unresolvable row escapes this loop, taking down every design whose bobbin
            // or shape is given by NAME in the browser. Ask whether the row is usable.
            if (!datum.second.get_functional_description()) {
                unusableBobbins.push_back(datum.first + " (no functionalDescription)");
                continue;
            }
            auto coreShapeName = datum.second.get_functional_description()->get_shape();
            auto coreShapeOrNone = try_find_core_shape_by_name(coreShapeName);
            if (!coreShapeOrNone) {
                unusableBobbins.push_back(datum.first + " -> core shape '" + coreShapeName + "' is not in the shape database");
                continue;
            }
            if (!datum.second.get_processed_description() ||
                datum.second.get_processed_description()->get_winding_windows().empty()) {
                unusableBobbins.push_back(datum.first + " (no processed winding window)");
                continue;
            }
            auto bobbinWindingWindow = datum.second.get_processed_description()->get_winding_windows()[0];
            if (!bobbinWindingWindow.get_area() || !bobbinWindingWindow.get_width() || !bobbinWindingWindow.get_height()) {
                unusableBobbins.push_back(datum.first + " (winding window has no area/width/height)");
                continue;
            }
            try {
                auto coreShape = coreShapeOrNone.value();
                auto corePiece = CorePiece::factory(coreShape);
                auto coreWindingWindow = corePiece->get_winding_window();
                if (!coreWindingWindow.get_area() || !coreWindingWindow.get_width() || !coreWindingWindow.get_height()) {
                    // Same reason as the checks above: a bad_optional_access here would be a
                    // throw where a skip is meant, and would escape in an exceptionless build.
                    unusableBobbins.push_back(datum.first + " -> core shape '" + coreShapeName + "' has no processed winding window");
                    continue;
                }

                auto bobbinWindingWindowArea = datum.second.get_processed_description()->get_winding_windows()[0].get_area().value();
                auto coreShapeWindingWindowArea = corePiece->get_winding_window().get_area().value() * 2; // Because if we are using a bobbin we have a two piece set
                double bobbinFillingFactor = bobbinWindingWindowArea / coreShapeWindingWindowArea;
                double bobbinWindingWindowWidth = datum.second.get_processed_description()->get_winding_windows()[0].get_width().value();
                double bobbinWindingWindowHeight = datum.second.get_processed_description()->get_winding_windows()[0].get_height().value();
                double coreWindingWindowWidth = corePiece->get_winding_window().get_width().value();
                double coreWindingWindowHeight = corePiece->get_winding_window().get_height().value() * 2; // Because if we are using a bobbin we have a two piece set
                double bobbinWindingWindowWidthProportion = bobbinWindingWindowWidth / coreWindingWindowWidth;
                double bobbinWindingWindowHeightProportion = bobbinWindingWindowHeight / coreWindingWindowHeight;

                // Track minimum real-world wall/column thicknesses. Prefer the bobbin's
                // own processedDescription values (already populated by the family-specific
                // BobbinDataProcessor); fall back to (core - bobbin) leftover otherwise.
                auto bobbinPd = datum.second.get_processed_description();
                double sampleWallThickness = bobbinPd->get_wall_thickness();
                if (!(sampleWallThickness > 0)) {
                    sampleWallThickness = (coreWindingWindowHeight - bobbinWindingWindowHeight) / 2;
                }
                double sampleColumnThickness = bobbinPd->get_column_thickness();
                if (!(sampleColumnThickness > 0)) {
                    sampleColumnThickness = coreWindingWindowWidth - bobbinWindingWindowWidth;
                }
                if (sampleWallThickness > 0) {
                    minBobbinWallThickness = std::min(minBobbinWallThickness, sampleWallThickness);
                }
                if (sampleColumnThickness > 0) {
                    minBobbinColumnThickness = std::min(minBobbinColumnThickness, sampleColumnThickness);
                }
                AuxFillingFactorWidth bobbinAuxFillingFactorWidth = { bobbinWindingWindowWidth, bobbinFillingFactor };
                AuxFillingFactorHeight bobbinAuxFillingFactorHeight = { bobbinWindingWindowHeight, bobbinFillingFactor };
                auxFillingFactorWidth.push_back(bobbinAuxFillingFactorWidth);
                auxFillingFactorHeight.push_back(bobbinAuxFillingFactorHeight);
                AuxWindingWindowWidth coreAuxWindingWindowWidth = { coreWindingWindowWidth, bobbinWindingWindowWidthProportion };
                AuxWindingWindowHeight coreAuxWindingWindowHeight = { coreWindingWindowHeight, bobbinWindingWindowHeightProportion };
                auxWindingWindowWidth.push_back(coreAuxWindingWindowWidth);
                auxWindingWindowHeight.push_back(coreAuxWindingWindowHeight);
            }
            catch (const std::exception &e)
            {
                // Backstop for the one step left that can still only report failure by
                // throwing (CorePiece::factory, for a family it does not implement).
                unusableBobbins.push_back(datum.first + " (" + std::string(e.what()) + ")");
                continue;
            }
        }

        if (!unusableBobbins.empty()) {
            std::string entry = std::to_string(unusableBobbins.size()) + " of " +
                                std::to_string(bobbinDatabase.size()) +
                                " catalogue bobbins could not be used to fit the bobbin interpolators:";
            for (const auto& unusableBobbin : unusableBobbins) {
                entry += "\n  - " + unusableBobbin;
            }
            logEntry(entry, "Bobbin", 1);
        }

        {
            size_t n = auxFillingFactorWidth.size();
            std::vector<double> x, y;

            std::sort(auxFillingFactorWidth.begin(), auxFillingFactorWidth.end(), [](const AuxFillingFactorWidth& b1, const AuxFillingFactorWidth& b2) {
                return b1.windingWindowWidth < b2.windingWindowWidth;
            });
            minBobbinWidth = auxFillingFactorWidth[0].windingWindowWidth;
            maxBobbinWidth = auxFillingFactorWidth[n - 1].windingWindowWidth;

            for (size_t i = 0; i < n; i++) {
                if (x.size() == 0 || fabs(auxFillingFactorWidth[i].windingWindowWidth - x.back()) > 1e-9) {
                    x.push_back(auxFillingFactorWidth[i].windingWindowWidth);
                    y.push_back(auxFillingFactorWidth[i].fillingFactor);
                }
            }

            bobbinFillingFactorInterpWidth = tk::spline(x, y, tk::spline::cspline_hermite, true);

        }
        {
            size_t n = auxFillingFactorHeight.size();
            std::vector<double> x, y;

            std::sort(auxFillingFactorHeight.begin(), auxFillingFactorHeight.end(), [](const AuxFillingFactorHeight& b1, const AuxFillingFactorHeight& b2) {
                return b1.windingWindowHeight < b2.windingWindowHeight;
            });
            minBobbinHeight = auxFillingFactorHeight[0].windingWindowHeight;
            maxBobbinHeight = auxFillingFactorHeight[n - 1].windingWindowHeight;

            for (size_t i = 0; i < n; i++) {
                if (x.size() == 0 || fabs(auxFillingFactorHeight[i].windingWindowHeight - x.back()) > 1e-9) {
                    x.push_back(auxFillingFactorHeight[i].windingWindowHeight);
                    y.push_back(auxFillingFactorHeight[i].fillingFactor);
                }
            }

            bobbinFillingFactorInterpHeight = tk::spline(x, y, tk::spline::cspline_hermite, true);
        }

        {
            size_t n = auxWindingWindowWidth.size();
            std::vector<double> x, y;

            std::sort(auxWindingWindowWidth.begin(), auxWindingWindowWidth.end(), [](const AuxWindingWindowWidth& b1, const AuxWindingWindowWidth& b2) {
                return b1.windingWindowWidth < b2.windingWindowWidth;
            });
            minWindingWindowWidth = auxWindingWindowWidth[0].windingWindowWidth;
            maxWindingWindowWidth = auxWindingWindowWidth[n - 1].windingWindowWidth;

            for (size_t i = 0; i < n; i++) {
                if (x.size() == 0 || fabs(auxWindingWindowWidth[i].windingWindowWidth - x.back()) > 1e-9) {
                    x.push_back(auxWindingWindowWidth[i].windingWindowWidth);
                    y.push_back(auxWindingWindowWidth[i].windingWindowProportion);
                }
            }

            bobbinWindingWindowProportionInterpWidth = tk::spline(x, y, tk::spline::linear, false);
        }
        {
            size_t n = auxWindingWindowHeight.size();
            std::vector<double> x, y;

            std::sort(auxWindingWindowHeight.begin(), auxWindingWindowHeight.end(), [](const AuxWindingWindowHeight& b1, const AuxWindingWindowHeight& b2) {
                return b1.windingWindowHeight < b2.windingWindowHeight;
            });
            minWindingWindowHeight = auxWindingWindowHeight[0].windingWindowHeight;
            maxWindingWindowHeight = auxWindingWindowHeight[n - 1].windingWindowHeight;

            for (size_t i = 0; i < n; i++) {
                if (x.size() == 0 || fabs(auxWindingWindowHeight[i].windingWindowHeight - x.back()) > 1e-9) {
                    x.push_back(auxWindingWindowHeight[i].windingWindowHeight);
                    y.push_back(auxWindingWindowHeight[i].windingWindowProportion);
                }
            }

            bobbinWindingWindowProportionInterpHeight = tk::spline(x, y, tk::spline::linear, false);
        }
    }
}
double Bobbin::get_filling_factor(double windingWindowWidth, double windingWindowHeight){
    load_interpolators();

    windingWindowWidth = std::max(windingWindowWidth, minBobbinWidth);
    windingWindowWidth = std::min(windingWindowWidth, maxBobbinWidth);

    double fillingFactorWidth = bobbinFillingFactorInterpWidth(windingWindowWidth);

    windingWindowHeight = std::max(windingWindowHeight, minBobbinHeight);
    windingWindowHeight = std::min(windingWindowHeight, maxBobbinHeight);

    double fillingFactorHeight = bobbinFillingFactorInterpHeight(windingWindowHeight);

    return (fillingFactorWidth + fillingFactorHeight) / 2;
}

std::vector<double> Bobbin::get_winding_window_dimensions(double coreWindingWindowWidth, double coreWindingWindowHeight){

    load_interpolators();

    double coreWindingWindowWidthForInterpolator = coreWindingWindowWidth;
    coreWindingWindowWidthForInterpolator = std::max(coreWindingWindowWidthForInterpolator, minWindingWindowWidth);
    coreWindingWindowWidthForInterpolator = std::min(coreWindingWindowWidthForInterpolator, maxWindingWindowWidth);
    double bobbinWindingWindowWidthProportion = bobbinWindingWindowProportionInterpWidth(coreWindingWindowWidthForInterpolator);
    // The proportion is bobbinWindow/coreWindow, which is physically ≤ 1.
    // Spline extrapolation outside the sample range can produce values > 1
    // (or < 0); clamp to keep the result physical.
    bobbinWindingWindowWidthProportion = std::clamp(bobbinWindingWindowWidthProportion, 0.0, 0.999);
    double bobbinWindingWindowWidth = bobbinWindingWindowWidthProportion * coreWindingWindowWidth;

    double coreWindingWindowHeightForInterpolator = coreWindingWindowHeight;
    coreWindingWindowHeightForInterpolator = std::max(coreWindingWindowHeightForInterpolator, minWindingWindowHeight);
    coreWindingWindowHeightForInterpolator = std::min(coreWindingWindowHeightForInterpolator, maxWindingWindowHeight);
    double bobbinWindingWindowHeighProportion = bobbinWindingWindowProportionInterpHeight(coreWindingWindowHeightForInterpolator);
    bobbinWindingWindowHeighProportion = std::clamp(bobbinWindingWindowHeighProportion, 0.0, 0.999);
    double bobbinWindingWindowHeight = bobbinWindingWindowHeighProportion * coreWindingWindowHeight;

    auto minimumThickness = std::min(coreWindingWindowWidth - bobbinWindingWindowWidth, (coreWindingWindowHeight - bobbinWindingWindowHeight) / 2);
    double maximumDisproportion = 1.2;  // hardcoded

    if ((coreWindingWindowWidth - bobbinWindingWindowWidth) > minimumThickness * maximumDisproportion) {
        bobbinWindingWindowWidth = coreWindingWindowWidth - minimumThickness * maximumDisproportion;
    }

    if ((coreWindingWindowHeight - bobbinWindingWindowHeight) / 2 > minimumThickness * maximumDisproportion) {
        bobbinWindingWindowHeight = coreWindingWindowHeight - minimumThickness * maximumDisproportion * 2;
    }

    return {bobbinWindingWindowWidth, bobbinWindingWindowHeight};
}

Bobbin Bobbin::create_quick_bobbin(double windingWindowHeight, double windingWindowWidth, ColumnShape shape) {
    CoreBobbinProcessedDescription coreBobbinProcessedDescription;
    WindingWindowElement windingWindowElement;

    windingWindowElement.set_height(windingWindowHeight);
    windingWindowElement.set_width(windingWindowWidth);
    windingWindowElement.set_area(windingWindowHeight * windingWindowWidth);

    windingWindowElement.set_coordinates(std::vector<double>({windingWindowWidth, 0, 0}));
    coreBobbinProcessedDescription.set_winding_windows(std::vector<WindingWindowElement>({windingWindowElement}));
    coreBobbinProcessedDescription.set_wall_thickness(0.001);
    coreBobbinProcessedDescription.set_column_thickness(0.001);
    coreBobbinProcessedDescription.set_column_shape(shape);
    coreBobbinProcessedDescription.set_column_depth(windingWindowWidth / 2);
    coreBobbinProcessedDescription.set_column_width(windingWindowWidth / 2);

    Bobbin bobbin;
    bobbin.set_processed_description(coreBobbinProcessedDescription);
    return bobbin;
}

Bobbin Bobbin::create_quick_bobbin(Core core, bool nullDimensions) {
    if (!core.get_processed_description()) {
        core.process_data();
    }

    // Multi-window cores supported. Use the centre/main window (index 0)
    // to derive default wall/column thicknesses; create_quick_bobbin
    // (Core, double, double) below will iterate over all windows.
    auto coreWindingWindow = core.get_processed_description()->get_winding_windows()[0];

    WindingWindowShape bobbinWindingWindowShape;
    if (core.get_shape_family() == CoreShapeFamily::T) {
        bobbinWindingWindowShape = WindingWindowShape::ROUND;
    }
    else {
        bobbinWindingWindowShape = WindingWindowShape::RECTANGULAR;
    }

    CoreBobbinProcessedDescription coreBobbinProcessedDescription;
    WindingWindowElement windingWindowElement;

    double bobbinColumnThickness = 0;
    double bobbinWallThickness = 0;

    if (!nullDimensions && bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {

        std::vector<double> bobbinWindingWindowDimensions = get_winding_window_dimensions(coreWindingWindow.get_width().value(), coreWindingWindow.get_height().value());
        bobbinColumnThickness = coreWindingWindow.get_width().value() - bobbinWindingWindowDimensions[0];
        bobbinWallThickness = (coreWindingWindow.get_height().value() - bobbinWindingWindowDimensions[1]) / 2;
        // Floor the thicknesses at the smallest real bobbin in the database. The
        // proportion interpolator clamps at 0.999 when the core's window is outside
        // its training range, which produces ~µm walls that are physically impossible
        // (real injection-molded bobbins are at least ~0.3 mm). Use the database min
        // instead, then re-derive the bobbin window dimensions to stay consistent.
        if (std::isfinite(minBobbinWallThickness) && bobbinWallThickness < minBobbinWallThickness) {
            bobbinWallThickness = minBobbinWallThickness;
        }
        if (std::isfinite(minBobbinColumnThickness) && bobbinColumnThickness < minBobbinColumnThickness) {
            bobbinColumnThickness = minBobbinColumnThickness;
        }
        if (bobbinWallThickness <= 0) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "bobbinWallThickness cannot be negative or 0: " + std::to_string(bobbinWallThickness));
        }
    }
    return create_quick_bobbin(core, bobbinWallThickness, bobbinColumnThickness);
}

Bobbin Bobbin::create_quick_bobbin(Core core, double thickness) {
    return create_quick_bobbin(core, thickness, thickness);
}

Bobbin Bobbin::create_quick_bobbin(Core core, double wallThickness, double columnThickness) {
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException("Core has not been processed yet");
    }

    auto coreWindingWindows = core.get_processed_description()->get_winding_windows();

    WindingWindowShape bobbinWindingWindowShape;
    if (core.get_shape_family() == CoreShapeFamily::T) {
        bobbinWindingWindowShape = WindingWindowShape::ROUND;
    }
    else {
        bobbinWindingWindowShape = WindingWindowShape::RECTANGULAR;
    }

    // A catalogue coating lines the toroid on ALL sides. Core::process_data already shrinks the
    // winding-window radial height (the bore) by the coating; the OUTER wrap, however, is derived
    // from the bobbin column, so fold the coating into the toroid column thickness here. It flows
    // into column_width (= coreColumn.width/2 + columnThickness) below, so the outer passes wrap
    // the coated OD for EVERY wind path — autocomplete and the coil-only wind alike, since the
    // latter reuses this resolved bobbin. get_air_cored_reluctance uses (column_width -
    // column_thickness), so the coating cancels there and the magnetic path stays bare ferrite;
    // the residual effect is on winding geometry (turn-to-core distance, mean turn length), which
    // is physically correct. Only an explicit coating on a toroid; a no-op otherwise.
    if (bobbinWindingWindowShape == WindingWindowShape::ROUND && core.get_functional_description().get_coating()) {
        columnThickness += core.get_coating_thickness();
    }

    auto coreCentralColumn = core.get_processed_description()->get_columns()[0];
    CoreBobbinProcessedDescription coreBobbinProcessedDescription;
    std::vector<WindingWindowElement> bobbinWindingWindows;
    bool anyFallback = false;

    for (size_t windowIndex = 0; windowIndex < coreWindingWindows.size(); ++windowIndex) {
        auto& coreWindingWindow = coreWindingWindows[windowIndex];
        WindingWindowElement windingWindowElement;
        std::vector<double> bobbinWindingWindowDimensions;

        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            bobbinWindingWindowDimensions = {std::max(0.0, coreWindingWindow.get_width().value() - columnThickness), std::max(0.0, coreWindingWindow.get_height().value() - wallThickness * 2)};
        }
        else {
            bobbinWindingWindowDimensions = {coreWindingWindow.get_radial_height().value(), coreWindingWindow.get_angle().value()};
        }

        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            if ((bobbinWindingWindowDimensions[0] < 0) || (bobbinWindingWindowDimensions[0] > 1) || (bobbinWindingWindowDimensions[1] < 0) || (bobbinWindingWindowDimensions[1] > 1)) {
                windingWindowElement.set_width(coreWindingWindow.get_width().value());
                windingWindowElement.set_height(coreWindingWindow.get_height().value());
                if (windowIndex == 0) {
                    windingWindowElement.set_coordinates(std::vector<double>({coreCentralColumn.get_width() / 2, 0, 0}));
                }
                else {
                    // Multi-window: mirror from the core's window coordinates so each
                    // additional window is placed symmetrically. The sign of the core
                    // window's x-coordinate determines the side.
                    auto coreCoords = coreWindingWindow.get_coordinates().value();
                    double sign = (coreCoords[0] < 0) ? -1.0 : 1.0;
                    windingWindowElement.set_coordinates(std::vector<double>({sign * coreCentralColumn.get_width() / 2, 0, 0}));
                }
                anyFallback = true;
            }
            else {
                windingWindowElement.set_width(bobbinWindingWindowDimensions[0]);
                windingWindowElement.set_height(bobbinWindingWindowDimensions[1]);
                windingWindowElement.set_area(bobbinWindingWindowDimensions[0] * bobbinWindingWindowDimensions[1]);
                if (windowIndex == 0) {
                    windingWindowElement.set_coordinates(std::vector<double>({coreCentralColumn.get_width() / 2 + columnThickness + bobbinWindingWindowDimensions[0] / 2, 0, 0}));
                }
                else {
                    // Mirror across the centre column for additional windows.
                    auto coreCoords = coreWindingWindow.get_coordinates().value();
                    double sign = (coreCoords[0] < 0) ? -1.0 : 1.0;
                    windingWindowElement.set_coordinates(std::vector<double>({sign * (coreCentralColumn.get_width() / 2 + columnThickness + bobbinWindingWindowDimensions[0] / 2), 0, 0}));
                }
            }
        }
        else {
            windingWindowElement.set_radial_height(bobbinWindingWindowDimensions[0]);
            windingWindowElement.set_angle(bobbinWindingWindowDimensions[1]);
            windingWindowElement.set_area(std::numbers::pi * pow(bobbinWindingWindowDimensions[0], 2) * bobbinWindingWindowDimensions[1] / 360);
            windingWindowElement.set_coordinates(std::vector<double>({bobbinWindingWindowDimensions[0], 0, 0}));
        }
        windingWindowElement.set_shape(bobbinWindingWindowShape);
        // Multi-column winding support: carry the window->column edge from the core
        // window so coil placement can resolve which column this bobbin window wraps.
        if (coreWindingWindow.get_column()) {
            windingWindowElement.set_column(coreWindingWindow.get_column());
        }
        // Carry a windingOrder set on the core's window (U serpentine vs Z dragback) so
        // Coil::get_winding_order finds it on the autocompleted bobbin (ABT #352).
        if (coreWindingWindow.get_winding_order()) {
            windingWindowElement.set_winding_order(coreWindingWindow.get_winding_order());
        }

        if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
            if ((windingWindowElement.get_width().value() < 0) || (windingWindowElement.get_width().value() > 1)) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened in section bobbin first : " + std::to_string(bobbinWindingWindowDimensions[0]));
            }
            if ((windingWindowElement.get_height().value() < 0) || (windingWindowElement.get_height().value() > 1)) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened in section bobbin second : " + std::to_string(bobbinWindingWindowDimensions[1]));
            }
        }
        else {
            if ((windingWindowElement.get_radial_height().value() < 0) || (windingWindowElement.get_radial_height().value() > 1)) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened in section bobbin first : " + std::to_string(bobbinWindingWindowDimensions[0]));
            }
            if ((windingWindowElement.get_angle().value() < 0) || (windingWindowElement.get_angle().value() > 360)) {
                throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened in section bobbin second : " + std::to_string(bobbinWindingWindowDimensions[1]));
            }
        }
        bobbinWindingWindows.push_back(windingWindowElement);
    }

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        if (anyFallback) {
            coreBobbinProcessedDescription.set_wall_thickness(0);
            coreBobbinProcessedDescription.set_column_thickness(0);
        }
        else {
            coreBobbinProcessedDescription.set_wall_thickness(wallThickness);
            coreBobbinProcessedDescription.set_column_thickness(columnThickness);
        }
    }
    else {
        coreBobbinProcessedDescription.set_wall_thickness(0);
        // Toroid: no bobbin plastic, but carry the coating folded into columnThickness above so
        // the coated-OD outer wrap is preserved (and cancels out of get_air_cored_reluctance).
        coreBobbinProcessedDescription.set_column_thickness(columnThickness);

        // THE EDGE THE WIRE IS PULLED OVER. A turn on a toroid is a closed loop around the ring
        // cross-section, so it bends at that section's edges and nowhere else -- there is no
        // moulded former corner here, and the injection-moulding fallback in
        // get_column_corner_radius() (half the wall thickness) is the wrong rule for a part that
        // has no wall. MAS names this case directly: cornerRadius is "the radius of the core
        // cross-section edges (or of its coating, when coated), which is what the wire is pulled
        // over".
        //
        // What that radius is, from the sources: NOTHING publishes the bare ferrite edge as a
        // number. IEC 62317-12 dimensions a ring core as A/B/C and admits the chamfer only
        // through the effective height; the MMPA/IMA "Standard Specification for Ferrite Toroid
        // Cores" says only that "the toroid corners shall not be sharp or rough" (4.2.4-4.2.5);
        // Fair-Rite states its toroids are "supplied burnished to break sharp edges"; TDK gives
        // the treatment by size tier (small = edges rounded by tumbling, medium/large = chamfer)
        // and no dimension. The one quantity that IS dimensioned on the winding surface is the
        // COATING: Ferroxcube draws PA11 at ~0.3 mm on its TN ring cores, TDK specifies epoxy
        // < 0.4 mm and parylene at 12.7 or 25 um, and Magnetics' coated-vs-uncoated limits imply
        // ~0.3-0.4 mm of epoxy buildup per surface.
        //
        // So the jacket is what we can source, and it is also what the wire actually touches: a
        // conformal coating cannot reproduce an edge sharper than itself, and the tumbled or
        // chamfered ferrite underneath only adds to it. Core::get_toroid_edge_radius() carries
        // that resolution -- the drawn ring-core value, floored at this part's own coating and
        // clamped to what the ring section can geometrically carry.
        //
        // It is deliberately NOT get_coating_thickness(). That datum is the dielectric path
        // normal to the flat faces, which is what StrayCapacitance integrates over and what the
        // outer-wrap fold above uses; this one is the curvature at the corner. Reusing a single
        // number for both would tie the windability verdict to a capacitance constant, and moving
        // either would silently move the other.
        //
        // Unlike that OUTER-WRAP fold -- which takes only an EXPLICIT coating, because inflating
        // the wound OD from a jacket the catalogue never claimed would move mean turn length and
        // reluctance -- the edge radius resolves for every toroid, defaults and all: a toroid
        // ships jacketed even when the catalogue omits it, and a bare-ferrite toroid wound with
        // bare wire is not a real part.
        coreBobbinProcessedDescription.set_column_corner_radius(core.get_toroid_edge_radius());
    }

    // NOTE: column_depth/column_shape/column_width describe the centre/main
    // column only. Multi-window cores store per-window geometry in
    // winding_windows[]; per-column physical dimensions are not represented
    // (v2 work).
    coreBobbinProcessedDescription.set_winding_windows(bobbinWindingWindows);
    coreBobbinProcessedDescription.set_column_shape(coreCentralColumn.get_shape());
    coreBobbinProcessedDescription.set_column_depth(coreCentralColumn.get_depth() / 2 + columnThickness);
    coreBobbinProcessedDescription.set_column_width(coreCentralColumn.get_width() / 2 + columnThickness);
    coreBobbinProcessedDescription.set_coordinates(std::vector<double>({0, 0, 0}));

    Bobbin bobbin;
    bobbin.set_processed_description(coreBobbinProcessedDescription);
    auto windingWindowCoordinates = bobbin.get_winding_window_coordinates(0);
    return bobbin;
}

std::vector<double> Bobbin::get_winding_window_dimensions(size_t windingWindowIndex) {
    if (get_winding_window_shape(windingWindowIndex) == WindingWindowShape::RECTANGULAR) {
        double width = get_processed_description()->get_winding_windows()[windingWindowIndex].get_width().value();
        double height = get_processed_description()->get_winding_windows()[windingWindowIndex].get_height().value();
        return {width, height};
    }
    else {
        double radialHeight = get_processed_description()->get_winding_windows()[windingWindowIndex].get_radial_height().value();
        double angle = get_processed_description()->get_winding_windows()[windingWindowIndex].get_angle().value();
        return {radialHeight, angle};
    }
}


double Bobbin::get_winding_window_area(size_t windingWindowIndex) {
    if (windingWindowIndex >= get_processed_description()->get_winding_windows().size()) {
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA, "Winding window does not exist");
    }
    auto windingWindow = get_processed_description()->get_winding_windows()[windingWindowIndex];
    if (get_processed_description()->get_winding_windows()[windingWindowIndex].get_area()) {
        return windingWindow.get_area().value();
    }
    else {
        if (get_winding_window_shape(windingWindowIndex) == WindingWindowShape::RECTANGULAR) {
            double width = windingWindow.get_width().value();
            double height = windingWindow.get_height().value();
            return width * height;
        }
        else {
            double radialHeight = windingWindow.get_radial_height().value();
            double angle = windingWindow.get_angle().value();
            return std::numbers::pi * pow(radialHeight, 2) * angle / 360;
        }
    }
}

std::vector<double> Bobbin::get_winding_window_coordinates(size_t windingWindowIndex) {
    return get_processed_description()->get_winding_windows()[windingWindowIndex].get_coordinates().value();
}

std::pair<double, double> Bobbin::get_column_and_wall_thickness(size_t windingWindowIndex) {
    if (!get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed");
    }
    auto bobbinProcessedDescription = get_processed_description().value();

    double columnThickness = bobbinProcessedDescription.get_column_thickness();
    double wallThickness = bobbinProcessedDescription.get_wall_thickness();
    return {columnThickness, wallThickness};
}

WindingOrientation Bobbin::get_winding_window_sections_orientation(size_t windingWindowIndex) {
    if (windingWindowIndex >= get_processed_description()->get_winding_windows().size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Invalid windingWindowIndex: " + std::to_string(windingWindowIndex) + ", bobbin only has" + std::to_string(get_processed_description()->get_winding_windows().size()) + " winding windows.");
    }
    if (!get_processed_description()->get_winding_windows()[windingWindowIndex].get_sections_orientation()) {
        if (get_winding_window_shape() == WindingWindowShape::ROUND) {
            return defaults.defaultRoundWindowSectionsOrientation;
        }
        else {
            return defaults.defaultRectangularWindowSectionsOrientation;
        }
    }
    return get_processed_description()->get_winding_windows()[windingWindowIndex].get_sections_orientation().value();
}

CoilAlignment Bobbin::get_winding_window_sections_alignment(size_t windingWindowIndex) {
    if (windingWindowIndex >= get_processed_description()->get_winding_windows().size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Invalid windingWindowIndex: " + std::to_string(windingWindowIndex) + ", bobbin only has" + std::to_string(get_processed_description()->get_winding_windows().size()) + " winding windows.");
    }
    if (!get_processed_description()->get_winding_windows()[windingWindowIndex].get_sections_alignment()) {
        if (get_winding_window_shape() == WindingWindowShape::ROUND) {
            return defaults.defaultRoundWindowSectionsAlignment;
        }
        else {
            return defaults.defaultRectangularWindowSectionsAlignment;
        }
    }
    return get_processed_description()->get_winding_windows()[windingWindowIndex].get_sections_alignment().value();
}

WindingWindowShape Bobbin::get_winding_window_shape(size_t windingWindowIndex) {
    if (get_processed_description() && windingWindowIndex >= get_processed_description()->get_winding_windows().size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Invalid windingWindowIndex: " + std::to_string(windingWindowIndex) + ", bobbin only has" + std::to_string(get_processed_description()->get_winding_windows().size()) + " winding windows.");
    }
    if (!get_processed_description()) {
        // ABT #631: find_bobbin_by_name returns a bobbin with NEITHER description for the
        // documented placeholder names ("basic"/"Basic"/"Dummy"/"None"), and Coil::resolve_bobbin
        // hands that straight to callers. Dereferencing the disengaged optional here is undefined
        // behaviour, not an exception: it reads whatever the returned-by-value optional's storage
        // happens to contain as a std::string. Say what is actually wrong instead.
        if (!get_functional_description()) {
            throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                "Bobbin has neither a processedDescription nor a functionalDescription: it is an "
                "unresolved placeholder ('basic'/'Dummy'/'None') and must be resolved against the "
                "core (magnetic_autocomplete, or Bobbin::create_quick_bobbin) before its winding "
                "window can be queried.");
        }
        auto coreShapeName = get_functional_description()->get_shape();
        auto coreShape = find_core_shape_by_name(coreShapeName);
        if (coreShape.get_family() == CoreShapeFamily::T) {
            return WindingWindowShape::ROUND;
        }
        else {
            return WindingWindowShape::RECTANGULAR;
        }
    }
    if (!get_processed_description()->get_winding_windows()[windingWindowIndex].get_shape()) {
        return WindingWindowShape::RECTANGULAR;
    }
    return get_processed_description()->get_winding_windows()[windingWindowIndex].get_shape().value();
}

void Bobbin::process_data() {
    // ABT #763: guard before the factory too. factory() takes its Bobbin BY VALUE, so a
    // guard there protects its own copy only if it is reached; this is the public entry
    // point every binding calls, and it must not be able to hand a description-less bobbin
    // any further down.
    if (!get_functional_description()) {
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
            "Bobbin has no functionalDescription, so there is nothing to process. "
            "process_data expects a bobbin json OBJECT with functionalDescription.family and "
            "functionalDescription.dimensions; a json string, an array, a number or an empty "
            "object all arrive here as an empty bobbin.");
    }

    auto processor = BobbinDataProcessor::factory(*this);
    auto processedDescription = (*processor).process_data(*this);

    // ABT #763 / ABT #634: a processed bobbin whose winding window has no area, or an area
    // of zero, is never a valid answer — it is the signature of a family whose processor
    // read dimension keys the bobbin does not declare (flatten_dimensions returns 0 for a
    // missing key without complaining; that is exactly how the 10 "Bobbin EI …" rows
    // declaring family "etd" while carrying the E dimension set processed to zero area,
    // silently, for as long as they did). Refuse to return a plausible zero.
    if (processedDescription.get_winding_windows().empty()) {
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
            "Processing bobbin '" + (get_name() ? get_name().value() : std::string("<unnamed>")) +
            "' produced no winding window at all.");
    }
    auto windingWindow = processedDescription.get_winding_windows()[0];
    if (!windingWindow.get_area() || !(windingWindow.get_area().value() > 0)) {
        // get_functional_description() returns the optional BY VALUE; holding the
        // description in a named local is what keeps the dimensions map alive while it is
        // iterated (binding the range-for to `...->get_dimensions()` directly walks a map
        // owned by a temporary that is already gone — gcc's -Wdangling-pointer catches it).
        auto functionalDescription = get_functional_description().value();
        std::string declaredDimensions;
        for (const auto& [key, _] : functionalDescription.get_dimensions()) {
            if (!declaredDimensions.empty()) {
                declaredDimensions += ", ";
            }
            declaredDimensions += key;
        }
        throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
            "Processing bobbin '" + (get_name() ? get_name().value() : std::string("<unnamed>")) +
            "' as family '" + to_string(functionalDescription.get_family()) +
            "' produced a winding window of zero area. The family's processor did not find the "
            "dimensions it reads; the bobbin declares {" + declaredDimensions + "}. Either the "
            "declared family does not match the declared dimension set, or a dimension is missing.");
    }

    set_processed_description(processedDescription);
}

bool Bobbin::check_if_fits(double dimension, bool isHorizontalOrRadial, size_t windingWindowIndex) {
    if (get_winding_window_shape(windingWindowIndex) == WindingWindowShape::RECTANGULAR) {
        auto windingWindowDimensions = get_winding_window_dimensions(windingWindowIndex);
        if (isHorizontalOrRadial) {
            return dimension < windingWindowDimensions[0];
        }
        else {
            return dimension < windingWindowDimensions[1];
        }
    }
    else if (get_winding_window_shape(windingWindowIndex) == WindingWindowShape::ROUND) {
        auto windingWindowDimensions = get_winding_window_dimensions(windingWindowIndex);
        if (isHorizontalOrRadial) {
            return dimension < windingWindowDimensions[0];
        }
        else {
            auto halfPerimeter = std::numbers::pi * windingWindowDimensions[0];
            return dimension < halfPerimeter;
        }
    }
    else {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Unsupported winding window shape: " + to_string(get_winding_window_shape(windingWindowIndex)));
    }
}



void Bobbin::set_winding_orientation(WindingOrientation windingOrientation, size_t windingWindowIndex) {
    if (!get_processed_description()) {
        throw CoilNotProcessedException("Bobbin has not been processed yet");
    }

    auto bobbinProcessedDescription = get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    if (windingWindowIndex >= windingWindows.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "windingWindowIndex out of range: " + std::to_string(windingWindowIndex));
    }
    windingWindows[windingWindowIndex].set_sections_orientation(windingOrientation);
    bobbinProcessedDescription.set_winding_windows(windingWindows);
    set_processed_description(bobbinProcessedDescription);
}


std::optional<WindingOrientation> Bobbin::get_winding_orientation(size_t windingWindowIndex) {
    if (!get_processed_description()) {
        return std::nullopt;
    }

    auto bobbinProcessedDescription = get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    if (windingWindowIndex >= windingWindows.size()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "windingWindowIndex out of range: " + std::to_string(windingWindowIndex));
    }
    if (windingWindows[windingWindowIndex].get_sections_orientation()) {
        return windingWindows[windingWindowIndex].get_sections_orientation().value();
    }

    return std::nullopt;
}

std::vector<double> Bobbin::get_maximum_dimensions() {
    if (!get_processed_description()) {
        process_data();
    }

    auto bobbinProcessedDescription = get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();

    double width = 0;
    double height = 0;
    double depth = 0;

    auto windingWindowDimensions = get_winding_window_dimensions();
    if (get_winding_window_shape(0) == WindingWindowShape::RECTANGULAR) {
        width = 2 * (bobbinProcessedDescription.get_column_width().value() + windingWindowDimensions[0]);
        height = 2 * bobbinProcessedDescription.get_wall_thickness() + windingWindowDimensions[1];
        depth = 2 * (bobbinProcessedDescription.get_column_depth() + windingWindowDimensions[0]);
    }
    else {
    }

    return {width, height, depth};
}

// ============================================================================
// Thermal Surface Area Calculations
// ============================================================================

double Bobbin::get_column_right_face_area(double coreDepth, size_t windingWindowIndex) {
    if (!get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed");
    }
    
    auto shape = get_processed_description()->get_column_shape();
    double windingWindowHeight = get_winding_window_height(windingWindowIndex);
    
    if (shape == ColumnShape::RECTANGULAR) {
        // For rectangular: height * depth
        return windingWindowHeight * coreDepth;
    }
    else {
        // For round: curved surface area = 2 * pi * r * h (for half cylinder facing winding)
        double columnRadius = get_column_width();  // column_width is radius for round columns
        return std::numbers::pi * columnRadius * windingWindowHeight;
    }
}

double Bobbin::get_column_top_face_area(double coreDepth, size_t windingWindowIndex) {
    if (!get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed");
    }
    
    auto shape = get_processed_description()->get_column_shape();
    double wallThickness = get_processed_description()->get_wall_thickness();
    
    if (shape == ColumnShape::RECTANGULAR) {
        // For rectangular: thickness * depth
        return wallThickness * coreDepth;
    }
    else {
        // For round: circular segment area = pi * r² (for top surface)
        double columnRadius = get_column_width();
        return std::numbers::pi * columnRadius * columnRadius;
    }
}

double Bobbin::get_column_bottom_face_area(double coreDepth, size_t windingWindowIndex) {
    // Same as top face for symmetrical bobbins
    return get_column_top_face_area(coreDepth, windingWindowIndex);
}

double Bobbin::get_yoke_interior_face_area(double coreDepth, bool isTopYoke, size_t windingWindowIndex) {
    if (!get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed");
    }
    
    auto shape = get_processed_description()->get_column_shape();
    double windingWindowWidth = get_winding_window_width(windingWindowIndex);
    double columnWidth = get_column_width();
    
    if (shape == ColumnShape::RECTANGULAR) {
        // For rectangular: width from column edge to outer edge * depth
        // This is the face facing the winding window (bottom of top yoke / top of bottom yoke)
        return windingWindowWidth * coreDepth;
    }
    else {
        // For round: annular sector area
        // Inner radius = column_radius, outer radius = column_radius + winding_window_radial_height
        double innerRadius = columnWidth;
        double outerRadius = columnWidth + windingWindowWidth;  // windingWindowWidth is radial height for round
        // Full annulus area = pi * (R² - r²), but we use half (one side of bobbin)
        return 0.5 * std::numbers::pi * (outerRadius * outerRadius - innerRadius * innerRadius);
    }
}

double Bobbin::get_yoke_exterior_face_area(double coreDepth, bool isTopYoke, size_t windingWindowIndex) {
    if (!get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed");
    }
    
    auto shape = get_processed_description()->get_column_shape();
    double windingWindowWidth = get_winding_window_width(windingWindowIndex);
    double columnWidth = get_column_width();
    double wallThickness = get_processed_description()->get_wall_thickness();
    
    if (shape == ColumnShape::RECTANGULAR) {
        // For rectangular: exterior face is (winding_window_width + wall_thickness) * depth
        // This accounts for the yoke extending beyond the winding window
        return (windingWindowWidth + wallThickness) * coreDepth;
    }
    else {
        // For round: similar to interior but including the wall thickness extension
        double innerRadius = columnWidth;
        double outerRadius = columnWidth + windingWindowWidth + wallThickness;
        return 0.5 * std::numbers::pi * (outerRadius * outerRadius - innerRadius * innerRadius);
    }
}

double Bobbin::get_yoke_right_face_area(double wallThickness, double coreDepth, size_t windingWindowIndex) {
    // The right face of the yoke is always: wallThickness * coreDepth / 2
    // (half depth because we model one side of the bobbin)
    return wallThickness * coreDepth / 2.0;
}

double Bobbin::get_winding_window_height(size_t windingWindowIndex) {
    if (get_winding_window_shape(windingWindowIndex) == WindingWindowShape::RECTANGULAR) {
        return get_processed_description()->get_winding_windows()[windingWindowIndex].get_height().value();
    }
    else {
        // For round, we need to calculate the arc length or use the radial height
        // For thermal purposes, we use the winding window radial height
        return get_processed_description()->get_winding_windows()[windingWindowIndex].get_radial_height().value();
    }
}

double Bobbin::get_winding_window_width(size_t windingWindowIndex) {
    if (get_winding_window_shape(windingWindowIndex) == WindingWindowShape::RECTANGULAR) {
        return get_processed_description()->get_winding_windows()[windingWindowIndex].get_width().value();
    }
    else {
        // For round, width is the radial height
        return get_processed_description()->get_winding_windows()[windingWindowIndex].get_radial_height().value();
    }
}

double Bobbin::get_column_width() {
    if (!get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed");
    }
    return get_processed_description()->get_column_width().value();
}

double Bobbin::get_column_depth() {
    if (!get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed");
    }
    return get_processed_description()->get_column_depth();
}

double Bobbin::get_column_corner_radius() {
    if (!get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed");
    }
    auto processedDescription = get_processed_description().value();

    switch (processedDescription.get_column_shape()) {
        case ColumnShape::ROUND:
            // column_width IS the radius for a round column (see get_column_right_face_area),
            // and a cylinder is nothing but corner, so that radius is the bend radius.
            if (!processedDescription.get_column_width()) {
                throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                                            "Bobbin: a round column has no column width, so its "
                                            "radius -- which is what a turn bends around -- is "
                                            "not defined");
            }
            return processedDescription.get_column_width().value();
        case ColumnShape::OBLONG:
            // A stadium's ends are semicircles of half the depth, and column_depth is already
            // the half-dimension. This has always been derivable and was simply never used.
            return processedDescription.get_column_depth();
        case ColumnShape::RECTANGULAR:
        case ColumnShape::IRREGULAR:
            if (processedDescription.get_column_corner_radius()) {
                return processedDescription.get_column_corner_radius().value();
            }
            // No datum: a moulded bobbin still cannot have a sharp corner. The injection-moulding
            // rule puts the inside radius at half the wall thickness, which is a sourced design
            // rule rather than a tuned constant, and it is the same treatment the wall thickness
            // itself gets in create_quick_bobbin when a real bobbin is not available. A catalogue
            // bobbin should carry the radius off its drawing instead.
            return 0.5 * processedDescription.get_column_thickness();
    }
    throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                                "Bobbin: unknown column shape, cannot resolve its corner radius");
}

double Bobbin::get_column_corner_half_angle() {
    if (!get_processed_description()) {
        throw CoilNotProcessedException("Bobbin not processed");
    }
    switch (get_processed_description()->get_column_shape()) {
        case ColumnShape::RECTANGULAR:
        case ColumnShape::IRREGULAR:
            // Two perpendicular faces meet at each corner of the racetrack.
            return 0.25 * std::numbers::pi;
        case ColumnShape::ROUND:
        case ColumnShape::OBLONG:
            // Nothing meets at an angle: the wire never leaves the curve, so there is no
            // lift-off term to apply.
            return 0.5 * std::numbers::pi;
    }
    throw InvalidInputException(ErrorCode::INVALID_BOBBIN_DATA,
                                "Bobbin: unknown column shape, cannot resolve its corner angle");
}


} // namespace OpenMagnetics
