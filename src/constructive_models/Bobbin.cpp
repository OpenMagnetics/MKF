#include "support/Utils.h"
#include "constructive_models/Bobbin.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <magic_enum.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
#include "spline.h"

tk::spline bobbinFillingFactorInterpWidth;
tk::spline bobbinFillingFactorInterpHeight;
tk::spline bobbinWindingWindowProportionInterpWidth;
tk::spline bobbinWindingWindowProportionInterpHeight;
double minBobbinWidth;
double maxBobbinWidth;
double minBobbinHeight;
double maxBobbinHeight;
double minWindingWindowWidth;
double maxWindingWindowWidth;
double minWindingWindowHeight;
double maxWindingWindowHeight;

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
            std::vector<double> coordinates({dimensions["f"] / 2 + dimensions["s1"], 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["l2"] - 2 * dimensions["s2"]);
            windingWindowElement.set_width((dimensions["e"] - dimensions["f"] - 2 * dimensions["s1"]) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
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
            std::vector<double> coordinates({dimensions["D2"] / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["H2"] - dimensions["H4"] - dimensions["H5"]);
            windingWindowElement.set_width((dimensions["D1"] - dimensions["D2"]) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
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
            std::vector<double> coordinates({dimensions["d2"] / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["h"] - 2 * dimensions["s"]);
            windingWindowElement.set_width((dimensions["d1"] - dimensions["d2"]) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);

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
            std::vector<double> coordinates({dimensions["d2"], 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["h2"]);
            windingWindowElement.set_width((dimensions["d1"] - dimensions["d2"]) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);

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
            std::vector<double> coordinates({dimensions["d2"] / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["h"] - dimensions["s1"] - dimensions["s2"]);
            windingWindowElement.set_width((dimensions["d1"] - dimensions["d2"]) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
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
            std::vector<double> coordinates({dimensions["D2"], 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["H2"]);
            windingWindowElement.set_width((dimensions["D1"] - dimensions["D2"]) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
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
            std::vector<double> coordinates({dimensions["D2"], 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["H2"]);
            windingWindowElement.set_width((dimensions["D1"] - dimensions["D2"]) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
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
            std::vector<double> coordinates({dimensions["f1"] / 2 + dimensions["S1"], 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(dimensions["d"] - 2 * dimensions["S2"]);
            windingWindowElement.set_width((dimensions["e"] - dimensions["f1"] - 2 * dimensions["S1"]) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

std::shared_ptr<BobbinDataProcessor> BobbinDataProcessor::factory(Bobbin bobbin) {

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
    else
        throw std::runtime_error("Unknown bobbin family, available options are: {E, EC, EFD, EP, ETD, PM, PQ, RM}");
}

void load_interpolators() {
    if (bobbinDatabase.empty() ||
        bobbinFillingFactorInterpWidth.get_x().size() == 0 || 
        bobbinFillingFactorInterpHeight.get_x().size() == 0 || 
        bobbinWindingWindowProportionInterpWidth.get_x().size() == 0 || 
        bobbinWindingWindowProportionInterpHeight.get_x().size() == 0) {
        load_bobbins();

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


        for (auto& datum : bobbinDatabase) {
            try {
                auto coreShapeName = datum.second.get_functional_description()->get_shape();
                auto coreShape = find_core_shape_by_name(coreShapeName);
                auto corePiece = CorePiece::factory(coreShape);

                auto bobbinWindingWindowArea = datum.second.get_processed_description()->get_winding_windows()[0].get_area().value();
                auto coreShapeWindingWindowArea = corePiece->get_winding_window().get_area().value() * 2; // Because if we are using a bobbin we have a two piece set
                double bobbinFillingFactor = bobbinWindingWindowArea / coreShapeWindingWindowArea;
                double bobbinWindingWindowWidth = datum.second.get_processed_description()->get_winding_windows()[0].get_width().value();
                double bobbinWindingWindowHeight = datum.second.get_processed_description()->get_winding_windows()[0].get_height().value();
                double coreWindingWindowWidth = corePiece->get_winding_window().get_width().value();
                double coreWindingWindowHeight = corePiece->get_winding_window().get_height().value() * 2; // Because if we are using a bobbin we have a two piece set
                double bobbinWindingWindowWidthProportion = bobbinWindingWindowWidth / coreWindingWindowWidth;
                double bobbinWindingWindowHeightProportion = bobbinWindingWindowHeight / coreWindingWindowHeight;
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
                continue;
            }
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
    double bobbinWindingWindowWidth = bobbinWindingWindowWidthProportion * coreWindingWindowWidth;

    double coreWindingWindowHeightForInterpolator = coreWindingWindowHeight;
    coreWindingWindowHeightForInterpolator = std::max(coreWindingWindowHeightForInterpolator, minWindingWindowHeight);
    coreWindingWindowHeightForInterpolator = std::min(coreWindingWindowHeightForInterpolator, maxWindingWindowHeight);
    double bobbinWindingWindowHeighProportion = bobbinWindingWindowProportionInterpHeight(coreWindingWindowHeightForInterpolator);
    double bobbinWindingWindowHeight = bobbinWindingWindowHeighProportion * coreWindingWindowHeight;


    if (bobbinWindingWindowHeight > coreWindingWindowHeight) {
        // Not proud, but the library is missbehaving...
        // bobbinWindingWindowHeight = coreWindingWindowHeight - 0.001;
        throw std::runtime_error("bobbinWindingWindowHeight cannot be greater than coreWindingWindowHeight: " + std::to_string(bobbinWindingWindowHeight) + " < " + std::to_string(coreWindingWindowHeight));
    }

    if (bobbinWindingWindowWidth > coreWindingWindowWidth) {
        // Not proud, but the library is missbehaving...
        // bobbinWindingWindowWidth = coreWindingWindowWidth - 0.001;
        throw std::runtime_error("bobbinWindingWindowWidth cannot be greater than coreWindingWindowWidth: " + std::to_string(bobbinWindingWindowWidth) + " < " + std::to_string(coreWindingWindowWidth));
    }

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

Bobbin Bobbin::create_quick_bobbin(double windingWindowHeight, double windingWindowWidth) {
    CoreBobbinProcessedDescription coreBobbinProcessedDescription;
    WindingWindowElement windingWindowElement;

    windingWindowElement.set_height(windingWindowHeight);
    windingWindowElement.set_width(windingWindowWidth);
    windingWindowElement.set_area(windingWindowHeight * windingWindowWidth);

    windingWindowElement.set_coordinates(std::vector<double>({windingWindowWidth, 0, 0}));
    coreBobbinProcessedDescription.set_winding_windows(std::vector<WindingWindowElement>({windingWindowElement}));
    coreBobbinProcessedDescription.set_wall_thickness(0.001);
    coreBobbinProcessedDescription.set_column_thickness(0.001);
    coreBobbinProcessedDescription.set_column_shape(ColumnShape::ROUND);
    coreBobbinProcessedDescription.set_column_depth(windingWindowWidth / 2);
    coreBobbinProcessedDescription.set_column_width(windingWindowWidth / 2);

    Bobbin bobbin;
    bobbin.set_processed_description(coreBobbinProcessedDescription);
    return bobbin;
}

Bobbin Bobbin::create_quick_bobbin(Core core, bool nullDimensions) {
    if (!core.get_processed_description()) {
        throw std::runtime_error("Core has not been processed yet");
    }

    if (core.get_processed_description()->get_winding_windows().size() > 1) {
        throw std::runtime_error("More than one winding window not supported yet");
    }

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
        if (bobbinWallThickness <= 0) {
            throw std::runtime_error("bobbinWallThickness cannot be negative or 0: " + std::to_string(bobbinWallThickness));
        }
    }
    return create_quick_bobbin(core, bobbinWallThickness, bobbinColumnThickness);
}

Bobbin Bobbin::create_quick_bobbin(Core core, double thickness) {
    return create_quick_bobbin(core, thickness, thickness);
}

Bobbin Bobbin::create_quick_bobbin(Core core, double wallThickness, double columnThickness) {
    if (!core.get_processed_description()) {
        throw std::runtime_error("Core has not been processed yet");
    }

    if (core.get_processed_description()->get_winding_windows().size() > 1) {
        throw std::runtime_error("More than one winding window not supported yet");
    }

    auto coreWindingWindow = core.get_processed_description()->get_winding_windows()[0];

    WindingWindowShape bobbinWindingWindowShape;
    if (core.get_shape_family() == CoreShapeFamily::T) {
        bobbinWindingWindowShape = WindingWindowShape::ROUND;
    }
    else {
        bobbinWindingWindowShape = WindingWindowShape::RECTANGULAR;
    }

    std::vector<double> bobbinWindingWindowDimensions;
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        bobbinWindingWindowDimensions = {std::max(0.0, coreWindingWindow.get_width().value() - columnThickness), std::max(0.0, coreWindingWindow.get_height().value() - wallThickness * 2)};
    }
    else {
        bobbinWindingWindowDimensions = {coreWindingWindow.get_radial_height().value(), coreWindingWindow.get_angle().value()};
    }

    CoreBobbinProcessedDescription coreBobbinProcessedDescription;
    WindingWindowElement windingWindowElement;

    auto coreCentralColumn = core.get_processed_description()->get_columns()[0];
    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        if ((bobbinWindingWindowDimensions[0] < 0) || (bobbinWindingWindowDimensions[0] > 1) || (bobbinWindingWindowDimensions[1] < 0) || (bobbinWindingWindowDimensions[1] > 1)) {
            windingWindowElement.set_width(coreWindingWindow.get_width().value());
            windingWindowElement.set_height(coreWindingWindow.get_height().value());
            windingWindowElement.set_coordinates(std::vector<double>({coreCentralColumn.get_width() / 2, 0, 0}));
            coreBobbinProcessedDescription.set_wall_thickness(0);
            coreBobbinProcessedDescription.set_column_thickness(0);
        }
        else {
            windingWindowElement.set_width(bobbinWindingWindowDimensions[0]);
            windingWindowElement.set_height(bobbinWindingWindowDimensions[1]);
            windingWindowElement.set_area(bobbinWindingWindowDimensions[0] * bobbinWindingWindowDimensions[1]);
            windingWindowElement.set_coordinates(std::vector<double>({coreCentralColumn.get_width() / 2 + columnThickness + bobbinWindingWindowDimensions[0] / 2, 0, 0}));
            coreBobbinProcessedDescription.set_wall_thickness(wallThickness);
            coreBobbinProcessedDescription.set_column_thickness(columnThickness);
        }
    }
    else {
        windingWindowElement.set_radial_height(bobbinWindingWindowDimensions[0]);
        windingWindowElement.set_angle(bobbinWindingWindowDimensions[1]);
        windingWindowElement.set_area(std::numbers::pi * pow(bobbinWindingWindowDimensions[0], 2) * bobbinWindingWindowDimensions[1] / 360);
        windingWindowElement.set_coordinates(std::vector<double>({bobbinWindingWindowDimensions[0], 0, 0}));
        coreBobbinProcessedDescription.set_wall_thickness(0);
        coreBobbinProcessedDescription.set_column_thickness(0);

    }
    windingWindowElement.set_shape(bobbinWindingWindowShape);
    coreBobbinProcessedDescription.set_winding_windows(std::vector<WindingWindowElement>({windingWindowElement}));
    coreBobbinProcessedDescription.set_column_shape(coreCentralColumn.get_shape());
    coreBobbinProcessedDescription.set_column_depth(coreCentralColumn.get_depth() / 2 + columnThickness);
    coreBobbinProcessedDescription.set_column_width(coreCentralColumn.get_width() / 2 + columnThickness);
    coreBobbinProcessedDescription.set_coordinates(std::vector<double>({0, 0, 0}));

    if (bobbinWindingWindowShape == WindingWindowShape::RECTANGULAR) {
        if ((windingWindowElement.get_width().value() < 0) || (windingWindowElement.get_width().value() > 1)) {
            throw std::runtime_error("Something wrong happened in section bobbin first : " + std::to_string(bobbinWindingWindowDimensions[0]));
        }
        if ((windingWindowElement.get_height().value() < 0) || (windingWindowElement.get_height().value() > 1)) {
            throw std::runtime_error("Something wrong happened in section bobbin second : " + std::to_string(bobbinWindingWindowDimensions[1]));
        }
    }
    else {
        if ((windingWindowElement.get_radial_height().value() < 0) || (windingWindowElement.get_radial_height().value() > 1)) {
            throw std::runtime_error("Something wrong happened in section bobbin first : " + std::to_string(bobbinWindingWindowDimensions[0]));
        }
        if ((windingWindowElement.get_angle().value() < 0) || (windingWindowElement.get_angle().value() > 360)) {
            throw std::runtime_error("Something wrong happened in section bobbin second : " + std::to_string(bobbinWindingWindowDimensions[1]));
        }
    }
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
        throw std::runtime_error("Winding window does not exist");
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

WindingOrientation Bobbin::get_winding_window_sections_orientation(size_t windingWindowIndex) {
    if (windingWindowIndex >= get_processed_description()->get_winding_windows().size()) {
        throw std::runtime_error("Invalid windingWindowIndex: " + std::to_string(windingWindowIndex) + ", bobbin only has" + std::to_string(get_processed_description()->get_winding_windows().size()) + " winding windows.");
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
        throw std::runtime_error("Invalid windingWindowIndex: " + std::to_string(windingWindowIndex) + ", bobbin only has" + std::to_string(get_processed_description()->get_winding_windows().size()) + " winding windows.");
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
    if (windingWindowIndex >= get_processed_description()->get_winding_windows().size()) {
        throw std::runtime_error("Invalid windingWindowIndex: " + std::to_string(windingWindowIndex) + ", bobbin only has" + std::to_string(get_processed_description()->get_winding_windows().size()) + " winding windows.");
    }
    if (!get_processed_description()) {
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
    auto processor = BobbinDataProcessor::factory(*this);
    set_processed_description((*processor).process_data(*this));
}

bool Bobbin::check_if_fits(double dimension, bool isHorizontalOrRadial, size_t windingWindowIndex) {
    if (get_winding_window_shape() == WindingWindowShape::RECTANGULAR) {
        auto windingWindowDimensions = get_winding_window_dimensions();
        if (isHorizontalOrRadial) {
            return dimension < windingWindowDimensions[0];
        }
        else {
            return dimension < windingWindowDimensions[1];
        }
    }
    else if (get_winding_window_shape() == WindingWindowShape::ROUND) {
        auto windingWindowDimensions = get_winding_window_dimensions();
        if (isHorizontalOrRadial) {
            return dimension < windingWindowDimensions[0];
        }
        else {
            auto halfPerimeter = std::numbers::pi * windingWindowDimensions[0];
            return dimension < halfPerimeter;
        }
    }
    else {
        throw std::runtime_error("Unsupported winding window shape: " + std::string{magic_enum::enum_name(get_winding_window_shape())});
    }
}



void Bobbin::set_winding_orientation(WindingOrientation windingOrientation, size_t windingWindowIndex) {
    if (!get_processed_description()) {
        throw std::runtime_error("Boobbin has not been processed yet");
    }

    auto bobbinProcessedDescription = get_processed_description().value();
    auto windingWindows = bobbinProcessedDescription.get_winding_windows();
    if (windingWindows.size() > 1) {
        throw std::runtime_error("Bobbins with more than winding window not implemented yet");
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
    if (windingWindows.size() > 1) {
        throw std::runtime_error("Bobbins with more than winding window not implemented yet");
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


} // namespace OpenMagnetics
