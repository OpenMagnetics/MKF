#include "Utils.h"
#include "CoreWrapper.h"
#include "BobbinWrapper.h"
#include "json.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <magic_enum.hpp>
#include <nlohmann/json-schema.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
#include <libInterpolate/Interpolate.hpp>
#include "../tests/TestingUtils.h"

using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;

_1D::LinearInterpolator<double> bobbinFillingFactorInterpWidth;
_1D::LinearInterpolator<double> bobbinFillingFactorInterpHeight;
_1D::LinearInterpolator<double> bobbinWindingWindowInterpWidth;
_1D::LinearInterpolator<double> bobbinWindingWindowInterpHeight;
double minBobbinWidth;
double maxBobbinWidth;
double minBobbinHeight;
double maxBobbinHeight;
double minWindingWindowWidth;
double maxWindingWindowWidth;
double minWindingWindowHeight;
double maxWindingWindowHeight;

namespace OpenMagnetics {

std::map<std::string, Dimension> flatten_dimensions(Bobbin bobbin) {
    std::map<std::string, Dimension> dimensions = bobbin.get_functional_description().value().get_dimensions();
    std::map<std::string, Dimension> flattenedDimensions;
    for (auto& dimension : dimensions) {
        double value = resolve_dimensional_values<OpenMagnetics::DimensionalValues::NOMINAL>(dimension.second);
        flattenedDimensions[dimension.first] = value;
    }
    return flattenedDimensions;
}

class BobbinEDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin);
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_thickness(std::get<double>(dimensions["s1"]));
            processedDescription.set_wall_thickness(std::get<double>(dimensions["s2"]));
            WindingWindowElement windingWindowElement;
            std::vector<double> coordinates({std::get<double>(dimensions["f"]) / 2 + std::get<double>(dimensions["s1"]), 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(std::get<double>(dimensions["l2"]) - 2 * std::get<double>(dimensions["s2"]));
            windingWindowElement.set_width((std::get<double>(dimensions["e"]) - std::get<double>(dimensions["f"]) - 2 * std::get<double>(dimensions["s1"])) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinRmDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin);
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_thickness((std::get<double>(dimensions["D2"]) - std::get<double>(dimensions["D3"])) / 2);
            processedDescription.set_wall_thickness(std::get<double>(dimensions["H5"]));
            WindingWindowElement windingWindowElement;
            std::vector<double> coordinates({std::get<double>(dimensions["D2"]) / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(std::get<double>(dimensions["H2"]) - std::get<double>(dimensions["H4"]) - std::get<double>(dimensions["H5"]));
            windingWindowElement.set_width((std::get<double>(dimensions["D1"]) - std::get<double>(dimensions["D2"])) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinEpDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin);
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_thickness((std::get<double>(dimensions["d2"]) - std::get<double>(dimensions["d3"])) / 2);
            processedDescription.set_wall_thickness(std::get<double>(dimensions["s"]));
            WindingWindowElement windingWindowElement;
            std::vector<double> coordinates({std::get<double>(dimensions["d2"]) / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(std::get<double>(dimensions["h"]) - 2 * std::get<double>(dimensions["s"]));
            windingWindowElement.set_width((std::get<double>(dimensions["d1"]) - std::get<double>(dimensions["d2"])) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);

            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinEtdDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin);
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_thickness((std::get<double>(dimensions["d2"]) - std::get<double>(dimensions["d3"])) / 2);
            processedDescription.set_wall_thickness((std::get<double>(dimensions["h1"]) - std::get<double>(dimensions["h2"])) / 2);
            WindingWindowElement windingWindowElement;
            std::vector<double> coordinates({std::get<double>(dimensions["d2"]), 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(std::get<double>(dimensions["h2"]));
            windingWindowElement.set_width((std::get<double>(dimensions["d1"]) - std::get<double>(dimensions["d2"])) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);

            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinPmDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin);
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_thickness((std::get<double>(dimensions["d2"]) - std::get<double>(dimensions["d3"])) / 2);
            processedDescription.set_wall_thickness(std::get<double>(dimensions["s1"]));
            WindingWindowElement windingWindowElement;
            std::vector<double> coordinates({std::get<double>(dimensions["d2"]) / 2, 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(std::get<double>(dimensions["h"]) - std::get<double>(dimensions["s1"]) - std::get<double>(dimensions["s2"]));
            windingWindowElement.set_width((std::get<double>(dimensions["d1"]) - std::get<double>(dimensions["d2"])) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinPqDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin);
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_thickness((std::get<double>(dimensions["D2"]) - std::get<double>(dimensions["D3"])) / 2);
            processedDescription.set_wall_thickness((std::get<double>(dimensions["H1"]) - std::get<double>(dimensions["H2"])) / 2);
            WindingWindowElement windingWindowElement;
            std::vector<double> coordinates({std::get<double>(dimensions["D2"]), 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(std::get<double>(dimensions["H2"]));
            windingWindowElement.set_width((std::get<double>(dimensions["D1"]) - std::get<double>(dimensions["D2"])) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinEcDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin);
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_thickness((std::get<double>(dimensions["D2"]) - std::get<double>(dimensions["D3"])) / 2);
            processedDescription.set_wall_thickness((std::get<double>(dimensions["H1"]) - std::get<double>(dimensions["H2"])) / 2);
            WindingWindowElement windingWindowElement;
            std::vector<double> coordinates({std::get<double>(dimensions["D2"]), 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(std::get<double>(dimensions["H2"]));
            windingWindowElement.set_width((std::get<double>(dimensions["D1"]) - std::get<double>(dimensions["D2"])) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

class BobbinEfdDataProcessor : public BobbinDataProcessor{
    public:
        CoreBobbinProcessedDescription process_data(Bobbin bobbin) {
            auto dimensions = flatten_dimensions(bobbin);
            CoreBobbinProcessedDescription processedDescription;
            processedDescription.set_column_thickness(std::get<double>(dimensions["S1"]));
            processedDescription.set_wall_thickness(std::get<double>(dimensions["S2"]));
            WindingWindowElement windingWindowElement;
            std::vector<double> coordinates({std::get<double>(dimensions["f1"]) / 2 + std::get<double>(dimensions["S1"]), 0});
            windingWindowElement.set_coordinates(coordinates);
            windingWindowElement.set_height(std::get<double>(dimensions["d"]) - 2 * std::get<double>(dimensions["S2"]));
            windingWindowElement.set_width((std::get<double>(dimensions["e"]) - std::get<double>(dimensions["f1"]) - 2 * std::get<double>(dimensions["S1"])) / 2);
            windingWindowElement.set_area(windingWindowElement.get_height().value() * windingWindowElement.get_width().value());
            processedDescription.get_mutable_winding_windows().push_back(windingWindowElement);
            processedDescription.set_coordinates(std::vector<double>({0, 0, 0}));
            return processedDescription;
        }
};

std::shared_ptr<BobbinDataProcessor> BobbinDataProcessor::factory(Bobbin bobbin) {

    auto family = bobbin.get_functional_description().value().get_family();
    if (family == BobbinFamily::E) {
        std::shared_ptr<BobbinDataProcessor> bobbinDataGenerator(new BobbinEDataProcessor);
        return bobbinDataGenerator;
    }
    else if (family == BobbinFamily::RM) {
        std::shared_ptr<BobbinDataProcessor> bobbinDataGenerator(new BobbinRmDataProcessor);
        return bobbinDataGenerator;
    }
    else if (family == BobbinFamily::EP) {
        std::shared_ptr<BobbinDataProcessor> bobbinDataGenerator(new BobbinEpDataProcessor);
        return bobbinDataGenerator;
    }
    else if (family == BobbinFamily::ETD) {
        std::shared_ptr<BobbinDataProcessor> bobbinDataGenerator(new BobbinEtdDataProcessor);
        return bobbinDataGenerator;
    }
    else if (family == BobbinFamily::PM) {
        std::shared_ptr<BobbinDataProcessor> bobbinDataGenerator(new BobbinPmDataProcessor);
        return bobbinDataGenerator;
    }
    else if (family == BobbinFamily::PQ) {
        std::shared_ptr<BobbinDataProcessor> bobbinDataGenerator(new BobbinPqDataProcessor);
        return bobbinDataGenerator;
    }
    else if (family == BobbinFamily::EC) {
        std::shared_ptr<BobbinDataProcessor> bobbinDataGenerator(new BobbinEcDataProcessor);
        return bobbinDataGenerator;
    }
    else if (family == BobbinFamily::EFD) {
        std::shared_ptr<BobbinDataProcessor> bobbinDataGenerator(new BobbinEfdDataProcessor);
        return bobbinDataGenerator;
    }
    else
        throw std::runtime_error("Unknown bobbin family, available options are: {E, EC, EFD, EP, ETD, PM, PQ, RM}");
}

void load_interpolators() {
    if (bobbinDatabase.empty() ||
        bobbinFillingFactorInterpWidth.getXData().size() == 0 || 
        bobbinFillingFactorInterpHeight.getXData().size() == 0 || 
        bobbinWindingWindowInterpWidth.getXData().size() == 0 || 
        bobbinWindingWindowInterpHeight.getXData().size() == 0) {
        load_databases(true);

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
            double windingWindowWidth, WindingWindow;
        };
        struct AuxWindingWindowHeight
        {
            double windingWindowHeight, WindingWindow;
        };
        std::vector<AuxWindingWindowWidth> auxWindingWindowWidth;
        std::vector<AuxWindingWindowHeight> auxWindingWindowHeight;


        for (auto& datum : bobbinDatabase) {
            auto coreShapeName = datum.second.get_functional_description().value().get_shape();
            auto coreShape = find_core_shape_by_name(coreShapeName);
            auto corePiece = OpenMagnetics::CorePiece::factory(coreShape);

            auto bobbinWindingWindowArea = datum.second.get_processed_description().value().get_winding_windows()[0].get_area().value();
            auto coreShapeWindingWindowArea = corePiece->get_winding_window().get_area().value() * 2; // Because if we are using a bobbin we have a two piece set
            double bobbinFillingFactor = bobbinWindingWindowArea / coreShapeWindingWindowArea;
            double bobbinWindingWindowWidth = datum.second.get_processed_description().value().get_winding_windows()[0].get_width().value();
            double bobbinWindingWindowHeight = datum.second.get_processed_description().value().get_winding_windows()[0].get_height().value();
            double coreWindingWindowWidth = corePiece->get_winding_window().get_width().value();
            double coreWindingWindowHeight = corePiece->get_winding_window().get_height().value() * 2; // Because if we are using a bobbin we have a two piece set
            AuxFillingFactorWidth bobbinAuxFillingFactorWidth = { bobbinWindingWindowWidth, bobbinFillingFactor };
            AuxFillingFactorHeight bobbinAuxFillingFactorHeight = { bobbinWindingWindowHeight, bobbinFillingFactor };
            auxFillingFactorWidth.push_back(bobbinAuxFillingFactorWidth);
            auxFillingFactorHeight.push_back(bobbinAuxFillingFactorHeight);
            AuxWindingWindowWidth coreAuxWindingWindowWidth = { coreWindingWindowWidth, bobbinWindingWindowWidth };
            AuxWindingWindowHeight coreAuxWindingWindowHeight = { coreWindingWindowHeight, bobbinWindingWindowHeight };
            auxWindingWindowWidth.push_back(coreAuxWindingWindowWidth);
            auxWindingWindowHeight.push_back(coreAuxWindingWindowHeight);
        }

        {
            size_t n = auxFillingFactorWidth.size();
            _1D::LinearInterpolator<double>::VectorType xx(n), yy(n);

            std::sort(auxFillingFactorWidth.begin(), auxFillingFactorWidth.end(), [](const AuxFillingFactorWidth& b1, const AuxFillingFactorWidth& b2) {
                return b1.windingWindowWidth < b2.windingWindowWidth;
            });
            minBobbinWidth = auxFillingFactorWidth[0].windingWindowWidth;
            maxBobbinWidth = auxFillingFactorWidth[n - 1].windingWindowWidth;

            for (size_t i = 0; i < n; i++) {
                xx(i) = auxFillingFactorWidth[i].windingWindowWidth;
                yy(i) = auxFillingFactorWidth[i].fillingFactor;
            }

            bobbinFillingFactorInterpWidth.setData(xx, yy);
        }
        {
            size_t n = auxFillingFactorWidth.size();
            _1D::LinearInterpolator<double>::VectorType xx(n), yy(n);

            std::sort(auxFillingFactorHeight.begin(), auxFillingFactorHeight.end(), [](const AuxFillingFactorHeight& b1, const AuxFillingFactorHeight& b2) {
                return b1.windingWindowHeight < b2.windingWindowHeight;
            });
            minBobbinHeight = auxFillingFactorHeight[0].windingWindowHeight;
            maxBobbinHeight = auxFillingFactorHeight[n - 1].windingWindowHeight;

            for (size_t i = 0; i < n; i++) {
                xx(i) = auxFillingFactorHeight[i].windingWindowHeight;
                yy(i) = auxFillingFactorHeight[i].fillingFactor;
            }

            bobbinFillingFactorInterpHeight.setData(xx, yy);
        }

        {
            size_t n = auxWindingWindowWidth.size();
            _1D::LinearInterpolator<double>::VectorType xx(n), yy(n);

            std::sort(auxWindingWindowWidth.begin(), auxWindingWindowWidth.end(), [](const AuxWindingWindowWidth& b1, const AuxWindingWindowWidth& b2) {
                return b1.windingWindowWidth < b2.windingWindowWidth;
            });
            minWindingWindowWidth = auxWindingWindowWidth[0].windingWindowWidth;
            maxWindingWindowWidth = auxWindingWindowWidth[n - 1].windingWindowWidth;

            for (size_t i = 0; i < n; i++) {
                xx(i) = auxWindingWindowWidth[i].windingWindowWidth;
                yy(i) = auxWindingWindowWidth[i].WindingWindow;
            }

            bobbinWindingWindowInterpWidth.setData(xx, yy);
        }
        {
            size_t n = auxWindingWindowWidth.size();
            _1D::LinearInterpolator<double>::VectorType xx(n), yy(n);

            std::sort(auxWindingWindowHeight.begin(), auxWindingWindowHeight.end(), [](const AuxWindingWindowHeight& b1, const AuxWindingWindowHeight& b2) {
                return b1.windingWindowHeight < b2.windingWindowHeight;
            });
            minWindingWindowHeight = auxWindingWindowHeight[0].windingWindowHeight;
            maxWindingWindowHeight = auxWindingWindowHeight[n - 1].windingWindowHeight;

            for (size_t i = 0; i < n; i++) {
                xx(i) = auxWindingWindowHeight[i].windingWindowHeight;
                yy(i) = auxWindingWindowHeight[i].WindingWindow;
            }

            bobbinWindingWindowInterpHeight.setData(xx, yy);
        }
    }
}
double BobbinWrapper::get_filling_factor(double windingWindowWidth, double windingWindowHeight){
    load_interpolators();

    windingWindowWidth = std::max(windingWindowWidth, minBobbinWidth);
    windingWindowWidth = std::min(windingWindowWidth, maxBobbinWidth);

    double fillingFactorWidth = bobbinFillingFactorInterpWidth(windingWindowWidth);

    windingWindowHeight = std::max(windingWindowHeight, minBobbinHeight);
    windingWindowHeight = std::min(windingWindowHeight, maxBobbinHeight);

    double fillingFactorHeight = bobbinFillingFactorInterpHeight(windingWindowHeight);

    return (fillingFactorWidth + fillingFactorHeight) / 2;
}

std::pair<double, double> BobbinWrapper::get_winding_window_dimensions(double coreWindingWindowWidth, double coreWindingWindowHeight){

    load_interpolators();
    double bobbinWindingWindowWidth;
    if (coreWindingWindowWidth > maxWindingWindowWidth) {
        double bobbinMaxWindingWindowWidth = bobbinWindingWindowInterpWidth(maxWindingWindowWidth);
        bobbinWindingWindowWidth = coreWindingWindowWidth - (maxWindingWindowWidth - bobbinMaxWindingWindowWidth);
    }
    else if (coreWindingWindowWidth < minWindingWindowWidth) {
        double bobbinMinWindingWindowWidth = bobbinWindingWindowInterpWidth(minWindingWindowWidth);
        bobbinWindingWindowWidth = std::max(coreWindingWindowWidth / 2, coreWindingWindowWidth - (minWindingWindowWidth - bobbinMinWindingWindowWidth));
    }
    else {
        bobbinWindingWindowWidth = bobbinWindingWindowInterpWidth(coreWindingWindowWidth);
    }

    double bobbinWindingWindowHeight;
    if (coreWindingWindowHeight > maxWindingWindowHeight) {
        double bobbinMaxWindingWindowHeight = bobbinWindingWindowInterpHeight(maxWindingWindowHeight);
        bobbinWindingWindowHeight = coreWindingWindowHeight - (maxWindingWindowHeight - bobbinMaxWindingWindowHeight);
    }
    else if (coreWindingWindowHeight < minWindingWindowHeight) {
        double bobbinMinWindingWindowHeight = bobbinWindingWindowInterpHeight(minWindingWindowHeight);
        bobbinWindingWindowHeight = std::max(coreWindingWindowHeight / 2, coreWindingWindowHeight - (minWindingWindowHeight - bobbinMinWindingWindowHeight));
    }
    else {
        bobbinWindingWindowHeight = bobbinWindingWindowInterpHeight(coreWindingWindowHeight);
    }

    return std::pair<double, double>({bobbinWindingWindowWidth, bobbinWindingWindowHeight});
}

BobbinWrapper BobbinWrapper::create_quick_bobbin(double windingWindowHeight, double windingWindowWidth) {
    json bobbinJson;

    CoreBobbinProcessedDescription coreBobbinProcessedDescription;
    WindingWindowElement windingWindowElement;

    windingWindowElement.set_height(windingWindowHeight);
    windingWindowElement.set_width(windingWindowWidth);
    windingWindowElement.set_coordinates(std::vector<double>({windingWindowWidth, 0, 0}));
    coreBobbinProcessedDescription.set_winding_windows(std::vector<WindingWindowElement>({windingWindowElement}));
    coreBobbinProcessedDescription.set_wall_thickness(0.001);
    coreBobbinProcessedDescription.set_column_thickness(0.001);
    coreBobbinProcessedDescription.set_column_shape(ColumnShape::ROUND);
    coreBobbinProcessedDescription.set_column_depth(windingWindowWidth / 2);
    coreBobbinProcessedDescription.set_column_width(windingWindowWidth / 2);

    BobbinWrapper bobbin;
    bobbin.set_processed_description(coreBobbinProcessedDescription);
    return bobbin;
}

BobbinWrapper BobbinWrapper::create_quick_bobbin(MagneticCore core) {
    json bobbinJson;

    CoreWrapper coreWrapper(core);

    if (coreWrapper.get_processed_description().value().get_winding_windows().size() > 1) {
        throw std::runtime_error("More than one winding window not supported yet");
    }

    auto coreWindingWindow = coreWrapper.get_processed_description().value().get_winding_windows()[0];
    auto coreCentralColumn = coreWrapper.get_processed_description().value().get_columns()[0];
    auto aux = get_winding_window_dimensions(coreWindingWindow.get_width().value(), coreWindingWindow.get_height().value());

    CoreBobbinProcessedDescription coreBobbinProcessedDescription;
    WindingWindowElement windingWindowElement;

    double bobbinColumnThickness = coreWindingWindow.get_width().value() - aux.first;
    double bobbinWallThickness = (coreWindingWindow.get_height().value() - aux.second) / 2;

    windingWindowElement.set_width(aux.first);
    windingWindowElement.set_height(aux.second);
    windingWindowElement.set_coordinates(std::vector<double>({coreCentralColumn.get_width() / 2 + bobbinColumnThickness + aux.first / 2, 0, 0}));
    coreBobbinProcessedDescription.set_winding_windows(std::vector<WindingWindowElement>({windingWindowElement}));
    coreBobbinProcessedDescription.set_wall_thickness(bobbinWallThickness);
    coreBobbinProcessedDescription.set_column_thickness(bobbinColumnThickness);
    coreBobbinProcessedDescription.set_column_shape(ColumnShape::ROUND);
    coreBobbinProcessedDescription.set_column_depth(coreCentralColumn.get_depth() / 2 + bobbinColumnThickness);
    coreBobbinProcessedDescription.set_column_width(coreCentralColumn.get_width() / 2 + bobbinColumnThickness);
    coreBobbinProcessedDescription.set_coordinates(std::vector<double>({0, 0, 0}));

    BobbinWrapper bobbin;
    bobbin.set_processed_description(coreBobbinProcessedDescription);
    return bobbin;
}

} // namespace OpenMagnetics
