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
double minBobbinWidth;
double maxBobbinWidth;
double minBobbinHeight;
double maxBobbinHeight;

namespace OpenMagnetics {

std::map<std::string, Dimension> flatten_dimensions(Bobbin bobbin) {
    std::map<std::string, Dimension> dimensions = bobbin.get_functional_description().get_dimensions();
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
            return processedDescription;
        }
};

std::shared_ptr<BobbinDataProcessor> BobbinDataProcessor::factory(Bobbin bobbin) {

    auto family = bobbin.get_functional_description().get_family();
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

double BobbinWrapper::get_filling_factor(double windingWindowWidth, double windingWindowHeight){
    bool loadInterpolateros = false;
    if (bobbinDatabase.empty()) {
        load_databases(true);

        std::vector<double> windingWindowWidths;
        std::vector<double> windingWindowHeights;
        std::vector<double> fillingFactors;
        struct AuxWidth
        {
            double windingWindowWidth, fillingFactor;
        };
        struct AuxHeight
        {
            double windingWindowHeight, fillingFactor;
        };
        std::vector<AuxWidth> auxWidth;
        std::vector<AuxHeight> auxHeight;


        for (auto& datum : bobbinDatabase) {
            auto coreShapeName = datum.second.get_functional_description().get_shape();
            auto coreShape = find_core_shape_by_name(coreShapeName);
            auto corePiece = OpenMagnetics::CorePiece::factory(coreShape);

            auto bobbinWindingWindowArea = datum.second.get_processed_description().value().get_winding_windows()[0].get_area().value();
            auto coreShapeWindingWindowArea = corePiece->get_winding_window().get_area().value() * 2; // Because if we are using a bobbin we have a two piece set
            double bobbingFillingFactor = bobbinWindingWindowArea / coreShapeWindingWindowArea;
            double bobbinWindingWindowWidth = datum.second.get_processed_description().value().get_winding_windows()[0].get_width().value();
            double bobbinWindingWindowHeight = datum.second.get_processed_description().value().get_winding_windows()[0].get_height().value();
            AuxWidth bobbingAuxWidth = { bobbinWindingWindowWidth, bobbingFillingFactor };
            AuxHeight bobbingAuxHeight = { bobbinWindingWindowHeight, bobbingFillingFactor };
            auxWidth.push_back(bobbingAuxWidth);
            auxHeight.push_back(bobbingAuxHeight);
        }

        {
            size_t n = auxWidth.size();
            _1D::LinearInterpolator<double>::VectorType xx(n), yy(n);

            std::sort(auxWidth.begin(), auxWidth.end(), [](const AuxWidth& b1, const AuxWidth& b2) {
                return b1.windingWindowWidth < b2.windingWindowWidth;
            });
            minBobbinWidth = auxWidth[0].windingWindowWidth;
            maxBobbinWidth = auxWidth[n - 1].windingWindowWidth;

            for (size_t i = 0; i < n; i++) {
                xx(i) = auxWidth[i].windingWindowWidth;
                yy(i) = auxWidth[i].fillingFactor;
            }

            bobbinFillingFactorInterpWidth.setData(xx, yy);
        }
        {
            size_t n = auxWidth.size();
            _1D::LinearInterpolator<double>::VectorType xx(n), yy(n);

            std::sort(auxHeight.begin(), auxHeight.end(), [](const AuxHeight& b1, const AuxHeight& b2) {
                return b1.windingWindowHeight < b2.windingWindowHeight;
            });
            minBobbinHeight = auxHeight[0].windingWindowHeight;
            maxBobbinHeight = auxHeight[n - 1].windingWindowHeight;

            for (size_t i = 0; i < n; i++) {
                xx(i) = auxHeight[i].windingWindowHeight;
                yy(i) = auxHeight[i].fillingFactor;
            }

            bobbinFillingFactorInterpHeight.setData(xx, yy);
        }
    }

    windingWindowWidth = std::max(windingWindowWidth, minBobbinWidth);
    windingWindowWidth = std::min(windingWindowWidth, maxBobbinWidth);

    double fillingFactorWidth = bobbinFillingFactorInterpWidth(windingWindowWidth);

    windingWindowHeight = std::max(windingWindowHeight, minBobbinHeight);
    windingWindowHeight = std::min(windingWindowHeight, maxBobbinHeight);

    double fillingFactorHeight = bobbinFillingFactorInterpHeight(windingWindowHeight);

    return (fillingFactorWidth + fillingFactorHeight) / 2;
}
} // namespace OpenMagnetics
