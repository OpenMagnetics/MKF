#include "WireWrapper.h"
#include "WindingSkinEffectLosses.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <vector>
#include <algorithm>
#include "Constants.h"
#include "Utils.h"
#include "spline.h"

std::map<std::string, tk::spline> wireOuterDimensionInterps;
std::map<std::string, tk::spline> wireFillingFactorInterps;
std::map<std::string, tk::spline> wirePackingFactorInterps;
std::map<std::string, tk::spline> wireConductingAreaInterps;
std::map<std::string, double> minWireConductingDimensions;
std::map<std::string, double> maxWireConductingDimensions;
std::map<std::string, int64_t> minLitzWireNumberConductors;
std::map<std::string, int64_t> maxLitzWireNumberConductors;

namespace OpenMagnetics {
    std::optional<InsulationWireCoating> WireWrapper::resolve_coating(WireWrapper wire) {
        // If the coating is a string, we have to load its data from the database
        if (!wire.get_coating()) {
            return std::nullopt;
        }
        if (std::holds_alternative<std::string>(wire.get_coating().value())) {
            throw std::runtime_error("Coating database not implemented yet");
        }
        else {
            return std::get<InsulationWireCoating>(wire.get_coating().value());
        }

    }
    std::optional<InsulationWireCoating> WireWrapper::resolve_coating() {
        // If the coating is a string, we have to load its data from the database
        if (!get_coating()) {
            return std::nullopt;
        }
        if (std::holds_alternative<std::string>(get_coating().value())) {
            throw std::runtime_error("Coating database not implemented yet");
        }
        else {
            return std::get<InsulationWireCoating>(get_coating().value());
        }

    }

    WireWrapper WireWrapper::resolve_strand(WireWrapper wire) { 
        if (!wire.get_strand())
            throw std::runtime_error("Litz wire is missing strand information");

        // If the strand is a string, we have to load its data from the database
        if (std::holds_alternative<std::string>(wire.get_strand().value())) {
            auto strand = find_wire_by_name(std::get<std::string>(wire.get_strand().value()));

            return strand;
        }
        else {
            return WireWrapper(std::get<WireRound>(wire.get_strand().value()));
        }

    }

    WireWrapper WireWrapper::resolve_strand() { 
        if (!get_strand())
            throw std::runtime_error("Litz wire is missing strand information");

        // If the strand is a string, we have to load its data from the database
        if (std::holds_alternative<std::string>(get_strand().value())) {
            auto strand = find_wire_by_name(std::get<std::string>(get_strand().value()));

            return strand;
        }
        else {
            return WireWrapper(std::get<WireRound>(get_strand().value()));
        }

    }

    WireMaterial WireWrapper::resolve_material() { 
        if (get_type() == WireType::LITZ) {
            auto strand = resolve_strand();
            return strand.resolve_material();
        }

        if (!get_material())
            throw std::runtime_error("Wire is missing material information");

        auto material = get_material().value();
        // If the material is a string, we have to load its data from the database
        if (std::holds_alternative<std::string>(material)) {
            auto materialData = find_wire_material_by_name(std::get<std::string>(material));

            return materialData;
        }
        else {
            return std::get<WireMaterial>(material);
        }
    }

    WireMaterial WireWrapper::resolve_material(WireWrapper wire) { 
        if (wire.get_type() == WireType::LITZ) {
            auto strand = wire.resolve_strand();
            return strand.resolve_material();
        }
        if (!wire.get_material())
            throw std::runtime_error("Wire is missing material information");

        auto material = wire.get_material().value();
        // If the material is a string, we have to load its data from the database
        if (std::holds_alternative<std::string>(material)) {
            auto materialData = find_wire_material_by_name(std::get<std::string>(material));

            return materialData;
        }
        else {
            return std::get<WireMaterial>(material);
        }
    }

    void create_interpolators(std::optional<double> conductingDiameter,
                              std::optional<double> conductingWidth,
                              std::optional<double> conductingHeight,
                              int numberConductors,
                              std::optional<int> grade,
                              std::optional<int> numberLayers,
                              std::optional<double> thicknessLayers,
                              std::optional<InsulationWireCoatingType> insulationWireCoatingType,
                              std::optional<WireStandard> standard,
                              WireType wireType,
                              bool includeAirInCell,
                              std::string key) {
        struct InterpolatorDatum
        {
            double wireConductingDimension, wireOuterDimension, wireFillingFactor, wirePackingFactor;
        };

        std::vector<InterpolatorDatum> interpolatorData;

        for (auto& datum : wireDatabase) {
            if (datum.second.get_type() != wireType) {
                continue;
            }
            auto coatingOption = WireWrapper::resolve_coating(datum.second);
            if (coatingOption) {
                auto coating = coatingOption.value();
                if (standard && !datum.second.get_standard()) {
                    continue;
                }
                if (numberLayers && !coating.get_number_layers()) {
                    continue;
                }
                if (thicknessLayers && !coating.get_thickness_layers()) {
                    continue;
                }
                bool numberConductorsPortionOfAnd = true;

                if (numberConductors > 1) {

                    numberConductorsPortionOfAnd = datum.second.get_number_conductors().value() > 1;
                }


                bool gradePortionOfAnd = !grade;
                int strandGrade = 0;
                if (grade) {
                    if (wireType == WireType::LITZ) {
                        auto strandWire = WireWrapper::resolve_strand(datum.second);
                        auto strandCoatingOption = WireWrapper::resolve_coating(strandWire);
                        auto strandCoating = strandCoatingOption.value();
                        if (grade && !strandCoating.get_grade()) {
                            continue;
                        }
                        strandGrade = strandCoating.get_grade().value();
                        gradePortionOfAnd = strandGrade == grade.value();
                    }
                    else {
                        if (coating.get_grade()) {
                            gradePortionOfAnd = coating.get_grade().value() == grade.value();
                        }
                    }
                }

                if (gradePortionOfAnd && 
                    numberConductorsPortionOfAnd && 
                    (!numberLayers || (numberLayers && coating.get_number_layers().value() == numberLayers.value())) && 
                    (!insulationWireCoatingType || (insulationWireCoatingType && coating.get_type().value() == insulationWireCoatingType.value())) && 
                    (!thicknessLayers || (thicknessLayers && coating.get_thickness_layers().value() == thicknessLayers.value())) && 
                    (!standard || (standard && datum.second.get_standard().value() == standard.value()))) {

                    double wireOuterDimension;
                    if (conductingDiameter) {
                        wireOuterDimension = resolve_dimensional_values(datum.second.get_outer_diameter().value()); 
                    }
                    else if (conductingWidth) {
                        wireOuterDimension = resolve_dimensional_values(datum.second.get_outer_width().value()); 
                    }
                    else if (conductingHeight) {
                        wireOuterDimension = resolve_dimensional_values(datum.second.get_outer_height().value()); 
                    }
                    else {
                        throw std::runtime_error("Missing wire dimension");
                    }
                    double wireConductingDimension;
                    double wirePackingFactor;
                    if (wireType == WireType::LITZ) {
                        auto strandWire = WireWrapper::resolve_strand(datum.second);
                        wireConductingDimension = resolve_dimensional_values(strandWire.get_conducting_diameter().value());

                        auto outerStrandDiameter = resolve_dimensional_values(strandWire.get_outer_diameter().value());
                        if (wireConductingDimension != conductingDiameter.value())
                            wirePackingFactor = 0;
                        else if (coating.get_thickness_layers())
                            wirePackingFactor = (wireOuterDimension - 2 * coating.get_number_layers().value() * coating.get_thickness_layers().value()) / (sqrt(datum.second.get_number_conductors().value()) * outerStrandDiameter);
                        else
                            wirePackingFactor = (wireOuterDimension) / (sqrt(datum.second.get_number_conductors().value()) * outerStrandDiameter);

                    }
                    else {
                        wirePackingFactor = 1;
                        if (conductingDiameter) {
                            wireConductingDimension = resolve_dimensional_values(datum.second.get_conducting_diameter().value()); 
                        }
                        else if (conductingWidth) {
                            wireConductingDimension = resolve_dimensional_values(datum.second.get_conducting_width().value()); 
                        }
                        else if (conductingHeight) {
                            wireConductingDimension = resolve_dimensional_values(datum.second.get_conducting_height().value()); 
                        }
                        else {
                            throw std::runtime_error("Missing wire dimension");
                        }
                    }
                    double wireNumberConductors = datum.second.get_number_conductors().value(); 
                    double outerArea;

                    if (conductingDiameter) {
                        if (includeAirInCell) {
                            outerArea = pow(wireOuterDimension, 2);
                        } else {
                            outerArea = std::numbers::pi * pow(wireOuterDimension / 2, 2);
                        }
                    }
                    else {
                        if (datum.second.get_conducting_area())
                            outerArea = resolve_dimensional_values(datum.second.get_conducting_area().value());
                        else {
                            outerArea = resolve_dimensional_values(datum.second.get_outer_width().value()) * resolve_dimensional_values(datum.second.get_outer_height().value());
                        }
                    }
                    double conductingArea = std::numbers::pi * pow(wireConductingDimension / 2, 2) * wireNumberConductors;
                    double wireFillingFactor = conductingArea / outerArea;
                    InterpolatorDatum interpolatorDatum = { wireConductingDimension, wireOuterDimension, wireFillingFactor, wirePackingFactor };

                    interpolatorData.push_back(interpolatorDatum);
                }
            }
        }
        if (interpolatorData.size() == 0) {
            throw std::runtime_error("No wires with that specification");
        }

        size_t n = interpolatorData.size();
        std::vector<double> x, xPackingFactor, yFillingFactor, yOuterDimension, yPackingFactor;
        std::sort(interpolatorData.begin(), interpolatorData.end(), [](const InterpolatorDatum& b1, const InterpolatorDatum& b2) {
            return b1.wireConductingDimension < b2.wireConductingDimension;
        });

        minWireConductingDimensions[key] = interpolatorData[0].wireConductingDimension;
        maxWireConductingDimensions[key] = interpolatorData[n - 1].wireConductingDimension;

        for (size_t i = 0; i < n; i++) {
            if (x.size() == 0 || interpolatorData[i].wireConductingDimension != x.back()) {
                x.push_back(interpolatorData[i].wireConductingDimension);
                yFillingFactor.push_back(interpolatorData[i].wireFillingFactor);
                yOuterDimension.push_back(interpolatorData[i].wireOuterDimension);
                if (wireType == WireType::LITZ && interpolatorData[i].wirePackingFactor > 0) {
                    xPackingFactor.push_back(interpolatorData[i].wireConductingDimension);
                    yPackingFactor.push_back(interpolatorData[i].wirePackingFactor);
                }
            }
        }
        tk::spline interpFillingFactor(x, yFillingFactor, tk::spline::cspline_hermite, true);
        tk::spline interpOuterDimension(x, yOuterDimension, tk::spline::cspline_hermite, true);
        tk::spline interpPackingFactor;
        if (wireType == WireType::LITZ && xPackingFactor.size() > 0) {
            interpPackingFactor = tk::spline(xPackingFactor, yPackingFactor, tk::spline::cspline_hermite, true);
        }
        wireFillingFactorInterps[key] = interpFillingFactor;
        wireOuterDimensionInterps[key] = interpOuterDimension;
        // if (wireType == WireType::LITZ && xPackingFactor.size() > 0) {
        //     wirePackingFactorInterps[key] = interpPackingFactor;
        // }
    }

    void create_packing_factor_interpolator(std::optional<int> grade,
                                            std::optional<int> numberLayers,
                                            std::optional<double> thicknessLayers,
                                            std::optional<InsulationWireCoatingType> insulationWireCoatingType,
                                            std::optional<WireStandard> standard,
                                            std::string key) {
        struct InterpolatorDatum
        {
            double wireNumberConductors, wirePackingFactor;
        };

        std::vector<InterpolatorDatum> interpolatorData;

        for (auto& datum : wireDatabase) {
            if (datum.second.get_type() != WireType::LITZ) {
                continue;
            }
            auto coatingOption = WireWrapper::resolve_coating(datum.second);
            if (coatingOption) {
                auto coating = coatingOption.value();
                if (standard && !datum.second.get_standard()) {
                    continue;
                }
                if (numberLayers && !coating.get_number_layers()) {
                    continue;
                }
                if (thicknessLayers && !coating.get_thickness_layers()) {
                    continue;
                }

                bool gradePortionOfAnd = !grade;
                int strandGrade = 0;
                if (grade) {
                    auto strandWire = WireWrapper::resolve_strand(datum.second);
                    auto strandCoatingOption = WireWrapper::resolve_coating(strandWire);
                    auto strandCoating = strandCoatingOption.value();
                    if (grade && !strandCoating.get_grade()) {
                        continue;
                    }
                    strandGrade = strandCoating.get_grade().value();
                    gradePortionOfAnd = strandGrade == grade.value();
                }
                if (gradePortionOfAnd && 
                    (!numberLayers || (numberLayers && coating.get_number_layers().value() == numberLayers.value())) && 
                    (!insulationWireCoatingType || (insulationWireCoatingType && coating.get_type().value() == insulationWireCoatingType.value())) && 
                    (!thicknessLayers || (thicknessLayers && coating.get_thickness_layers().value() == thicknessLayers.value())) && 
                    (!standard || (standard && datum.second.get_standard().value() == standard.value()))) {

                    double wireOuterDimension = resolve_dimensional_values(datum.second.get_outer_diameter().value()); 
                    double wirePackingFactor;

                    auto strandWire = WireWrapper::resolve_strand(datum.second);
                    double wireNumberConductors = datum.second.get_number_conductors().value(); 

                    auto outerStrandDiameter = resolve_dimensional_values(strandWire.get_outer_diameter().value());
                    if (coating.get_thickness_layers())
                        wirePackingFactor = (wireOuterDimension - 2 * coating.get_number_layers().value() * coating.get_thickness_layers().value()) / (sqrt(wireNumberConductors) * outerStrandDiameter);
                    else
                        wirePackingFactor = (wireOuterDimension) / (sqrt(wireNumberConductors) * outerStrandDiameter);

                    InterpolatorDatum interpolatorDatum = { wireNumberConductors, wirePackingFactor };
                    interpolatorData.push_back(interpolatorDatum);
                }
            }
        }
        if (interpolatorData.size() == 0) {
            throw std::runtime_error("No wires with that specification");
        }

        size_t n = interpolatorData.size();
        std::vector<double> x, y;
        std::sort(interpolatorData.begin(), interpolatorData.end(), [](const InterpolatorDatum& b1, const InterpolatorDatum& b2) {
            return b1.wireNumberConductors < b2.wireNumberConductors;
        });

        minLitzWireNumberConductors[key] = interpolatorData[0].wireNumberConductors;
        maxLitzWireNumberConductors[key] = interpolatorData[n - 1].wireNumberConductors;

        for (size_t i = 0; i < n; i++) {
            if (x.size() == 0 || interpolatorData[i].wireNumberConductors != x.back()) {
                x.push_back(interpolatorData[i].wireNumberConductors);
                y.push_back(interpolatorData[i].wirePackingFactor);
            }
        }
        tk::spline interpPackingFactor;
        if (x.size() > 0) {
            interpPackingFactor = tk::spline(x, y, tk::spline::cspline_hermite, true);
            wirePackingFactorInterps[key] = interpPackingFactor;
        }
    }

    void create_conducting_area_interpolator(std::optional<WireStandard> standard,
                                             std::string key) {
        struct InterpolatorDatum
        {
            double wireTheoreticalConductingArea, wireRealConductingArea;
        };

        std::vector<InterpolatorDatum> interpolatorData;

        for (auto& datum : wireDatabase) {
            if (datum.second.get_type() != WireType::RECTANGULAR) {
                continue;
            }
            auto coatingOption = WireWrapper::resolve_coating(datum.second);
            if (coatingOption) {
                auto coating = coatingOption.value();
                if (standard && !datum.second.get_standard()) {
                    continue;
                }

                if (!standard || (standard && datum.second.get_standard().value() == standard.value())) {

                    double wireRealConductingArea = resolve_dimensional_values(datum.second.get_conducting_area().value()); 
                    double wireTheoreticalConductingArea = resolve_dimensional_values(datum.second.get_conducting_width().value()) * resolve_dimensional_values(datum.second.get_conducting_height().value()); 

                    InterpolatorDatum interpolatorDatum = { wireTheoreticalConductingArea, wireRealConductingArea };
                    interpolatorData.push_back(interpolatorDatum);
                }
            }
        }
        if (interpolatorData.size() == 0) {
            throw std::runtime_error("No wires with that specification");
        }

        size_t n = interpolatorData.size();
        std::vector<double> x, y;
        std::sort(interpolatorData.begin(), interpolatorData.end(), [](const InterpolatorDatum& b1, const InterpolatorDatum& b2) {
            return b1.wireTheoreticalConductingArea < b2.wireTheoreticalConductingArea;
        });

        for (size_t i = 0; i < n; i++) {
            if (x.size() == 0 || interpolatorData[i].wireTheoreticalConductingArea != x.back()) {
                x.push_back(interpolatorData[i].wireTheoreticalConductingArea);
                y.push_back(interpolatorData[i].wireRealConductingArea);
            }
        }
        tk::spline interpPackingFactor;
        if (x.size() > 0) {
            interpPackingFactor = tk::spline(x, y, tk::spline::cspline_hermite, true);
            wireConductingAreaInterps[key] = interpPackingFactor;
        }
    }

    double get_filling_factor(std::optional<double> conductingDiameter,
                              std::optional<double> conductingWidth,
                              std::optional<double> conductingHeight,
                              int numberConductors,
                              std::optional<int> grade,
                              std::optional<int> numberLayers,
                              std::optional<double> thicknessLayers,
                              std::optional<InsulationWireCoatingType> insulationWireCoatingType,
                              std::optional<WireStandard> standard,
                              WireType wireType,
                              bool includeAirInCell,
                              std::string key) {
        if (wireDatabase.empty()) {
            load_databases(true);
        }

        if (!wireFillingFactorInterps.contains(key)) {
            create_interpolators(conductingDiameter,
                                 conductingWidth,
                                 conductingHeight,
                                 numberConductors,
                                 grade,
                                 numberLayers,
                                 thicknessLayers,
                                 insulationWireCoatingType,
                                 standard,
                                 wireType,
                                 includeAirInCell,
                                 key);
        }

        double wireConductingDimension;
        if (conductingDiameter) {
            wireConductingDimension = conductingDiameter.value();
        }
        else if (conductingWidth) {
            wireConductingDimension = conductingWidth.value();
        }
        else if (conductingHeight) {
            wireConductingDimension = conductingHeight.value();
        }
        else {
            throw std::runtime_error("Missing wire dimension");
        }
        if (wireConductingDimension < minWireConductingDimensions[key])
        wireConductingDimension = std::max(wireConductingDimension, minWireConductingDimensions[key]);
        wireConductingDimension = std::min(wireConductingDimension, maxWireConductingDimensions[key]);
        double fillingFactor = wireFillingFactorInterps[key](wireConductingDimension);
        return fillingFactor;
    }

    double get_outer_dimension(std::optional<double> conductingDiameter,
                               std::optional<double> conductingWidth,
                               std::optional<double> conductingHeight,
                               int numberConductors,
                               std::optional<int> grade,
                               std::optional<int> numberLayers,
                               std::optional<double> thicknessLayers,
                               std::optional<InsulationWireCoatingType> insulationWireCoatingType,
                               std::optional<WireStandard> standard,
                               WireType wireType,
                               std::string key) {
        if (wireDatabase.empty()) {
            load_databases(true);
        }

        if (!wireFillingFactorInterps.contains(key)) {
            create_interpolators(conductingDiameter,
                                 conductingWidth,
                                 conductingHeight,
                                 numberConductors,
                                 grade,
                                 numberLayers,
                                 thicknessLayers,
                                 insulationWireCoatingType,
                                 standard,
                                 wireType,
                                 false,
                                 key);
        }

        double wireConductingDimension;
        if (conductingDiameter) {
            wireConductingDimension = conductingDiameter.value();
        }
        else if (conductingWidth) {
            wireConductingDimension = conductingWidth.value();
        }
        else if (conductingHeight) {
            wireConductingDimension = conductingHeight.value();
        }
        else {
            throw std::runtime_error("Missing wire dimension");
        }
        wireConductingDimension = std::max(wireConductingDimension, minWireConductingDimensions[key]);
        wireConductingDimension = std::min(wireConductingDimension, maxWireConductingDimensions[key]);
        double outerDimension = wireOuterDimensionInterps[key](wireConductingDimension);
        return outerDimension;
    }

    double get_packing_factor(int numberConductors,
                               std::optional<int> grade,
                               std::optional<int> numberLayers,
                               std::optional<double> thicknessLayers,
                               std::optional<InsulationWireCoatingType> insulationWireCoatingType,
                               std::optional<WireStandard> standard,
                               std::string key) {
        if (wireDatabase.empty()) {
            load_databases(true);
        }

        if (!wirePackingFactorInterps.contains(key)) {
            create_packing_factor_interpolator(grade,
                                 numberLayers,
                                 thicknessLayers,
                                 insulationWireCoatingType,
                                 standard,
                                 key);
        }
        if (!wirePackingFactorInterps.contains(key)) {
            return 0;
        }
        else {

            auto wireNumberConductors = numberConductors;
            wireNumberConductors = std::max(wireNumberConductors, int(minLitzWireNumberConductors[key]));
            wireNumberConductors = std::min(wireNumberConductors, int(maxLitzWireNumberConductors[key]));
            double packingFactor = wirePackingFactorInterps[key](wireNumberConductors);
            return packingFactor;
        }
    }

    double get_conducting_area_rectangular_from_interpolator(double conductingWidth,
                                                             double conductingHeight,
                                                             std::optional<WireStandard> standard,
                                                             std::string key) {
        if (wireDatabase.empty()) {
            load_databases(true);
        }

        if (!wireConductingAreaInterps.contains(key)) {
            create_conducting_area_interpolator(standard,
                                                key);
        }
        if (!wireConductingAreaInterps.contains(key)) {
            return 0;
        }
        else {

            double wireTheoreticalConductingArea = conductingWidth * conductingHeight;
            double conductingArea = wireConductingAreaInterps[key](wireTheoreticalConductingArea);
            return conductingArea;
        }
    }

    double get_packing_factor_from_standard(WireStandard standard, double numberConductors) {
        if (standard == WireStandard::IEC_60317) {
            // Accoding to standard IEC 60317 - 11
            if (numberConductors < 12) {
                return 1.25;
            }
            else if (numberConductors < 16) {
                return 1.26;
            }
            else if (numberConductors < 20) {
                return 1.27;
            }
            else {
                return 1.28;
            }
        }
        else {
            // Accoding to Rubadue, page 25 of https://www.psma.com/sites/default/files/uploads/files/Litz%20Wire%20Practical%20Design%20Considerations%20for%20Todays%20HF%20Applications%20Jensen%2C%20Rubadue.pdf
            return 1.155;
        }
    }

    double get_serving_thickness_from_standard(int numberLayers, double outerDiameter) {
        if (numberLayers == 0) {
            return 0;
        }
        else if (numberLayers <= 2) {
            if (outerDiameter < 0.00045)
                return 0.000035 * numberLayers;
            else if (outerDiameter < 0.0006)
                return 0.00004 * numberLayers;
            else if (outerDiameter < 0.001)
                return 0.00007 * numberLayers;
            else
                return 0.00008 * numberLayers;
        }
        else {
            throw std::runtime_error("Unsupported number of layers in litz serving");
        }

    }

    // Thought for enamelled round wires
    double WireWrapper::get_filling_factor_round(double conductingDiameter, int grade, WireStandard standard, bool includeAirInCell) {

        auto wireType = WireType::ROUND;

        std::string standardString = " " + static_cast<std::string>(magic_enum::enum_name(standard));
        std::string wireTypeString = static_cast<std::string>(magic_enum::enum_name(wireType));
        std::string key = wireTypeString + " enamelled " + std::to_string(grade) + standardString;
        return get_filling_factor(conductingDiameter,
                                  std::nullopt,
                                  std::nullopt,
                                  1,
                                  grade,
                                  std::nullopt,
                                  std::nullopt,
                                  InsulationWireCoatingType::ENAMELLED,
                                  standard,
                                  wireType,
                                  includeAirInCell,
                                  key);

    }

    // Thought for enamelled round wires
    double WireWrapper::get_outer_diameter_round(double conductingDiameter, int grade, WireStandard standard) {
        auto wireType = WireType::ROUND;

        std::string standardString = " " + static_cast<std::string>(magic_enum::enum_name(standard));
        std::string wireTypeString = static_cast<std::string>(magic_enum::enum_name(wireType));
        std::string key = wireTypeString + " enamelled " + std::to_string(grade) + standardString;
        return get_outer_dimension(conductingDiameter,
                                   std::nullopt,
                                   std::nullopt,
                                   1,
                                   grade,
                                   std::nullopt,
                                   std::nullopt,
                                   InsulationWireCoatingType::ENAMELLED,
                                   standard,
                                   wireType,
                                   key);
    }

    // Thought for insulated round wires
    double WireWrapper::get_filling_factor_round(double conductingDiameter, int numberLayers, double thicknessLayers, WireStandard standard, bool includeAirInCell) {
        auto wireType = WireType::ROUND;
        if (numberLayers == 0 || numberLayers > 3) {
            throw std::runtime_error("Unsupported number of layers");
        }

        if (thicknessLayers != 5.08e-05 && thicknessLayers != 7.62e-05 && thicknessLayers != 3.81e-05 && thicknessLayers != 2.54e-05 && thicknessLayers != 1.27e-04) {
            throw std::runtime_error("Unsupported layer thickness: " + std::to_string(thicknessLayers));
        }

        std::string standardString = " " + static_cast<std::string>(magic_enum::enum_name(standard));
        std::string wireTypeString = static_cast<std::string>(magic_enum::enum_name(wireType));
        std::string key = wireTypeString + " insulated " + std::to_string(numberLayers) + " layers " + std::to_string(thicknessLayers * 1e6) + " µm " + standardString;
        return get_filling_factor(conductingDiameter,
                                  std::nullopt,
                                  std::nullopt,
                                  1,
                                  std::nullopt,
                                  numberLayers,
                                  thicknessLayers,
                                  InsulationWireCoatingType::INSULATED,
                                  standard,
                                  wireType,
                                  includeAirInCell,
                                  key);
    }

    // Thought for insulated round wires
    double WireWrapper::get_outer_diameter_round(double conductingDiameter, int numberLayers, double thicknessLayers, WireStandard standard) {
        auto wireType = WireType::ROUND;
        if (numberLayers == 0 || numberLayers > 3) {
            throw std::runtime_error("Unsupported number of layers");
        }

        if (thicknessLayers != 5.08e-05 && thicknessLayers != 7.62e-05 && thicknessLayers != 3.81e-05 && thicknessLayers != 2.54e-05 && thicknessLayers != 1.27e-04) {
            throw std::runtime_error("Unsupported layer thickness: " + std::to_string(thicknessLayers));
        }

        std::string standardString = " " + static_cast<std::string>(magic_enum::enum_name(standard));
        std::string wireTypeString = static_cast<std::string>(magic_enum::enum_name(wireType));
        std::string key = wireTypeString + " insulated " + std::to_string(numberLayers) + " layers " + std::to_string(thicknessLayers * 1e6) + " µm " + standardString;
        return get_outer_dimension(conductingDiameter,
                                   std::nullopt,
                                   std::nullopt,
                                   1,
                                   std::nullopt,
                                   numberLayers,
                                   thicknessLayers,
                                   InsulationWireCoatingType::INSULATED,
                                   standard,
                                   wireType,
                                   key);
    }

    // Thought for enamelled litz wires with or without serving
    double WireWrapper::get_filling_factor_served_litz(double conductingDiameter, int numberConductors, int grade, int numberLayers, WireStandard standard, bool includeAirInCell) {
        auto wireType = WireType::LITZ;
        std::string standardString = " " + static_cast<std::string>(magic_enum::enum_name(standard));
        std::string wireTypeString = static_cast<std::string>(magic_enum::enum_name(wireType));
        std::string key = wireTypeString + " enamelled " + std::to_string(grade) + standardString;

        return get_filling_factor(conductingDiameter,
                                  std::nullopt,
                                  std::nullopt,
                                  numberConductors,
                                  grade,
                                  numberLayers,
                                  std::nullopt,
                                  InsulationWireCoatingType::SERVED,
                                  standard,
                                  wireType,
                                  includeAirInCell,
                                  key);
    }

    double WireWrapper::get_outer_diameter_served_litz(double conductingDiameter, int numberConductors, int grade, int numberLayers, WireStandard standard) {
        double packingFactor = get_packing_factor_from_standard(standard, numberConductors);
        double outerStrandDiameter = get_outer_diameter_round(conductingDiameter, grade, standard);

        double outerDiameter = packingFactor * sqrt(numberConductors) * outerStrandDiameter;

        double servingThickness = get_serving_thickness_from_standard(numberLayers, outerDiameter);

        outerDiameter += servingThickness;

        return outerDiameter;
    }
    // Thought for insulated litz wires
    double WireWrapper::get_filling_factor_insulated_litz(double conductingDiameter, int numberConductors, int numberLayers, double thicknessLayers, int grade, WireStandard standard, bool includeAirInCell) {
        auto wireType = WireType::LITZ;

        std::string standardString = " " + static_cast<std::string>(magic_enum::enum_name(standard));
        std::string wireTypeString = static_cast<std::string>(magic_enum::enum_name(wireType));
        std::string key = wireTypeString + " insulated " + std::to_string(numberLayers) + " layers " + std::to_string(thicknessLayers * 1e6) + " µm " + standardString;

        return get_filling_factor(conductingDiameter,
                                  std::nullopt,
                                  std::nullopt,
                                  numberConductors,
                                  grade,
                                  numberLayers,
                                  thicknessLayers,
                                  InsulationWireCoatingType::INSULATED,
                                  standard,
                                  wireType,
                                  includeAirInCell,
                                  key);
    }

    // Thought for insulated litz wires
    double WireWrapper::get_outer_diameter_insulated_litz(double conductingDiameter, int numberConductors, int numberLayers, double thicknessLayers, int grade, WireStandard standard) {
        auto wireType = WireType::LITZ;

        std::string standardString = " " + static_cast<std::string>(magic_enum::enum_name(standard));
        std::string wireTypeString = static_cast<std::string>(magic_enum::enum_name(wireType));
        std::string key = wireTypeString + " insulated " + std::to_string(numberLayers) + " layers " + std::to_string(thicknessLayers * 1e6) + " µm " + standardString;

        auto packingFactor = get_packing_factor(numberConductors,
                                                grade,
                                                numberLayers,
                                                thicknessLayers,
                                                InsulationWireCoatingType::INSULATED,
                                                standard,
                                                key);
        if (packingFactor == 0) {
            packingFactor = get_packing_factor_from_standard(standard, numberConductors);
        }

        double outerStrandDiameter = get_outer_diameter_round(conductingDiameter, grade, standard);

        double unservedOuterDiameter = packingFactor * sqrt(numberConductors) * outerStrandDiameter;
        return unservedOuterDiameter + 2 * thicknessLayers * numberLayers;
    }

    // Thought for enamelled reactangular wires
    double WireWrapper::get_filling_factor_rectangular(double conductingWidth, double conductingHeight, int grade, WireStandard standard) {
        double realConductingArea = get_conducting_area_rectangular(conductingWidth, conductingHeight, standard);
        double outerWidth = get_outer_width_rectangular(conductingWidth, grade, standard);
        double outerHeight = get_outer_height_rectangular(conductingHeight, grade, standard);
        double outerArea = outerWidth * outerHeight;

        return realConductingArea / outerArea;
    }

    // Thought for enamelled reactangular wires
    double WireWrapper::get_outer_width_rectangular(double conductingWidth, int grade, WireStandard standard) {
        auto wireType = WireType::RECTANGULAR;

        std::string standardString = " " + static_cast<std::string>(magic_enum::enum_name(standard));
        std::string wireTypeString = static_cast<std::string>(magic_enum::enum_name(wireType));
        std::string key = wireTypeString + " enamelled width " + std::to_string(grade) + standardString;
        return get_outer_dimension(std::nullopt,
                                   conductingWidth,
                                   std::nullopt,
                                   1,
                                   grade,
                                   std::nullopt,
                                   std::nullopt,
                                   InsulationWireCoatingType::ENAMELLED,
                                   standard,
                                   wireType,
                                   key);
    }

    // Thought for enamelled reactangular wires
    double WireWrapper::get_outer_height_rectangular(double conductingHeight, int grade, WireStandard standard) {
        auto wireType = WireType::RECTANGULAR;

        std::string standardString = " " + static_cast<std::string>(magic_enum::enum_name(standard));
        std::string wireTypeString = static_cast<std::string>(magic_enum::enum_name(wireType));
        std::string key = wireTypeString + " enamelled height " + std::to_string(grade) + standardString;
        return get_outer_dimension(std::nullopt,
                                   std::nullopt,
                                   conductingHeight,
                                   1,
                                   grade,
                                   std::nullopt,
                                   std::nullopt,
                                   InsulationWireCoatingType::ENAMELLED,
                                   standard,
                                   wireType,
                                   key);
    }

    // Thought for enamelled reactangular wires
    double WireWrapper::get_conducting_area_rectangular(double conductingWidth, double conductingHeight, WireStandard standard) {
        auto wireType = WireType::RECTANGULAR;

        std::string standardString = " " + static_cast<std::string>(magic_enum::enum_name(standard));
        std::string wireTypeString = static_cast<std::string>(magic_enum::enum_name(wireType));
        std::string key = wireTypeString + " enamelled " + standardString;
        return get_conducting_area_rectangular_from_interpolator(conductingWidth, conductingHeight, standard, key);
    }

    int WireWrapper::get_equivalent_insulation_layers(double voltageToInsulate) {
        if (!resolve_coating()) {
            return 0;
        }
        auto coating = resolve_coating().value();
        // 0.85 according to IEC 61558
        // https://www.elektrisola.com/en/Products/Fully-Insulated-Wire/Breakdown-Voltage
        auto voltageCoveredByWire = coating.get_breakdown_voltage().value() * 0.85;
        int timesVoltageisCovered = floor(voltageCoveredByWire / voltageToInsulate);
        int numberLayers;
        if (coating.get_number_layers()) {
            numberLayers = coating.get_number_layers().value();
        }
        else {
            numberLayers = 1;
        }
        return std::min(numberLayers, timesVoltageisCovered);
    }

    double WireWrapper::calculate_conducting_area() {
        if (!get_number_conductors()) {
            if (get_type() == WireType::LITZ) {
                throw std::runtime_error("Missing number of conductors for wire");
            }
            else {
                set_number_conductors(1);
            }
        }
        switch (get_type()) {
            case WireType::LITZ:
                {
                    auto strand = resolve_strand();
                    if (!strand.get_conducting_diameter()) {
                        throw std::runtime_error("Missing conducting diameter in litz strand");
                    }
                    auto conductingDiameter = resolve_dimensional_values(strand.get_conducting_diameter().value());
                    return std::numbers::pi * pow(conductingDiameter / 2, 2) * get_number_conductors().value();
                }
            case WireType::ROUND:
                {
                    if (!get_conducting_diameter()) {
                        throw std::runtime_error("Missing conducting diameter in round wire");
                    }
                    auto conductingDiameter = resolve_dimensional_values(get_conducting_diameter().value());
                    return std::numbers::pi * pow(conductingDiameter / 2, 2) * get_number_conductors().value();
                }
            case WireType::RECTANGULAR:
                {
                    if (!get_conducting_width()) {
                        throw std::runtime_error("Missing conducting width in rectangular wire");
                    }
                    if (!get_conducting_height()) {
                        throw std::runtime_error("Missing conducting height in rectangular wire");
                    }
                    auto conductingWidth = resolve_dimensional_values(get_conducting_width().value());
                    auto conductingHeight = resolve_dimensional_values(get_conducting_height().value());
                    return get_conducting_area_rectangular(conductingWidth, conductingHeight) * get_number_conductors().value();
                }
            case WireType::FOIL:
                {
                    if (!get_conducting_width()) {
                        throw std::runtime_error("Missing conducting width in foil wire");
                    }
                    if (!get_conducting_height()) {
                        throw std::runtime_error("Missing conducting height in foil wire");
                    }
                    auto conductingWidth = resolve_dimensional_values(get_conducting_width().value());
                    auto conductingHeight = resolve_dimensional_values(get_conducting_height().value()) * get_number_conductors().value();
                    return conductingWidth * conductingHeight;
                }
            default:
                throw std::runtime_error("Unknow type of wire");
        }
    }

    double WireWrapper::calculate_effective_current_density(OperatingPointExcitation excitation, double temperature) {
        if (!excitation.get_current()) {
            throw std::runtime_error("Current is missing");
        }
        return calculate_effective_current_density(excitation.get_current().value(), temperature);
    }

    double WireWrapper::calculate_effective_current_density(SignalDescriptor current, double temperature) {
        if (!current.get_processed()) {
            throw std::runtime_error("Current is not processed");
        }
        if (!current.get_processed()->get_effective_frequency()) {
            throw std::runtime_error("Current is missing effective frequency");
        }
        if (!current.get_processed()->get_rms()) {
            throw std::runtime_error("Current is missing RMS");
        }
        double effectiveFrequency = current.get_processed()->get_effective_frequency().value();
        double rms = current.get_processed()->get_rms().value();
        return calculate_effective_current_density(rms, effectiveFrequency, temperature);
    }

    double WireWrapper::calculate_effective_conducting_area(double frequency, double temperature) {
        auto material = resolve_material();
        auto skinDepth = WindingSkinEffectLosses::calculate_skin_depth(material,  frequency, temperature);
        double conductingSmallestDimension;
        double effectiveConductingArea;
        switch (get_type()) {
            case WireType::LITZ: 
                {
                    auto strand = resolve_strand();
                    if (!strand.get_conducting_diameter()) {
                        throw std::runtime_error("Missing conducting diameter in litz strand");
                    }
                    conductingSmallestDimension = resolve_dimensional_values(strand.get_conducting_diameter().value());
                    break;
                }
            case WireType::ROUND: 
                {
                    if (!get_conducting_diameter()) {
                        throw std::runtime_error("Missing conducting diameter in round wire");
                    }
                    conductingSmallestDimension = resolve_dimensional_values(get_conducting_diameter().value());
                    break;
                }
            case WireType::RECTANGULAR:
            case WireType::FOIL:
                {
                    if (!get_conducting_width()) {
                        throw std::runtime_error("Missing conducting width in foil wire");
                    }
                    if (!get_conducting_height()) {
                        throw std::runtime_error("Missing conducting height in foil wire");
                    }
                    conductingSmallestDimension = std::min(resolve_dimensional_values(get_conducting_width().value()), resolve_dimensional_values(get_conducting_height().value()));
                    break;
                }
            default:
                throw std::runtime_error("Unknow type of wire");
        }

        auto conductingArea = calculate_conducting_area();
        if (!get_conducting_area()) {
            DimensionWithTolerance dimensionWithTolerance;
            dimensionWithTolerance.set_nominal(conductingArea);
            set_conducting_area(dimensionWithTolerance);
        }
        double nonConductingArea = 0;
        if (skinDepth < conductingSmallestDimension / 2) {
            switch (get_type()) {
                case WireType::LITZ:
                    {
                        auto strand = resolve_strand();
                        nonConductingArea = std::numbers::pi * pow(resolve_dimensional_values(strand.get_conducting_diameter().value()) / 2 - skinDepth, 2) * get_number_conductors().value();
                        break;
                    }
                case WireType::ROUND:
                    {
                        nonConductingArea = std::numbers::pi * pow(resolve_dimensional_values(get_conducting_diameter().value()) / 2 - skinDepth, 2) * get_number_conductors().value();
                        break;
                    }
                case WireType::RECTANGULAR:
                case WireType::FOIL:
                    {
                        nonConductingArea = (resolve_dimensional_values(get_conducting_width().value()) - skinDepth * 2) * (resolve_dimensional_values(get_conducting_height().value()) - skinDepth * 2) * get_number_conductors().value();
                        auto scaleDueToRoundCorners = conductingArea / (resolve_dimensional_values(get_conducting_width().value()) * resolve_dimensional_values(get_conducting_height().value()));
                        nonConductingArea *= scaleDueToRoundCorners;
                        break;
                    }
                default:
                    throw std::runtime_error("Unknow type of wire");
            }

            effectiveConductingArea = conductingArea - nonConductingArea;
        }
        else {
            effectiveConductingArea = conductingArea;
        }

        if (effectiveConductingArea < 0) {
            throw std::runtime_error("effectiveConductingArea cannot be negative: " + std::to_string(conductingArea) + ", " + std::to_string(nonConductingArea));
        }

        return effectiveConductingArea;
    }

    double WireWrapper::calculate_effective_current_density(double rms, double frequency, double temperature) {
        return rms / calculate_effective_conducting_area(frequency, temperature);
    }

    int WireWrapper::calculate_number_parallels_needed(InputsWrapper inputs, WireWrapper& wire, double maximumEffectiveCurrentDensity, size_t windingIndex) {
        int maximumNumberParallels = 0;
        for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex){
            double temperature = inputs.get_operating_points()[operatingPointIndex].get_conditions().get_ambient_temperature();
            auto excitation = inputs.get_winding_excitation(operatingPointIndex, windingIndex);
            auto effectiveCurrentDensity = wire.calculate_effective_current_density(excitation, temperature);
            maximumNumberParallels = std::max(maximumNumberParallels, int(ceil(effectiveCurrentDensity / maximumEffectiveCurrentDensity)));
        }

        return maximumNumberParallels;
    }

    int WireWrapper::calculate_number_parallels_needed(double rms, double effectiveFrequency, double temperature, WireWrapper& wire, double maximumEffectiveCurrentDensity) {
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(rms, effectiveFrequency, temperature);
        int numberParallels = int(ceil(effectiveCurrentDensity / maximumEffectiveCurrentDensity));
        return numberParallels;
    }

    int WireWrapper::calculate_number_parallels_needed(OperatingPointExcitation excitation, double temperature, WireWrapper& wire, double maximumEffectiveCurrentDensity) {
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(excitation, temperature);
        int numberParallels = int(ceil(effectiveCurrentDensity / maximumEffectiveCurrentDensity));
        return numberParallels;
    }

    int WireWrapper::calculate_number_parallels_needed(SignalDescriptor current, double temperature, WireWrapper& wire, double maximumEffectiveCurrentDensity) {
        auto effectiveCurrentDensity = wire.calculate_effective_current_density(current, temperature);
        int numberParallels = int(ceil(effectiveCurrentDensity / maximumEffectiveCurrentDensity));
        return numberParallels;
    }

    double WireWrapper::get_maximum_outer_width() {
        switch (get_type()) {
            case WireType::LITZ:
            case WireType::ROUND:
                if (get_outer_diameter())
                    return resolve_dimensional_values(get_outer_diameter().value());
                else 
                    return get_maximum_conducting_width();
            case WireType::RECTANGULAR:
            case WireType::FOIL:
                if (get_outer_width())
                    return resolve_dimensional_values(get_outer_width().value());
                else 
                    return resolve_dimensional_values(get_conducting_width().value());
            default:
                throw std::runtime_error("Unknow type of wire");
        }
    }

    double WireWrapper::get_maximum_outer_height() {
        switch (get_type()) {
            case WireType::LITZ:
            case WireType::ROUND:
                if (get_outer_diameter())
                    return resolve_dimensional_values(get_outer_diameter().value());
                else 
                    return get_maximum_conducting_height();
            case WireType::RECTANGULAR:
            case WireType::FOIL:
                if (get_outer_height())
                    return resolve_dimensional_values(get_outer_height().value());
                else 
                    return resolve_dimensional_values(get_conducting_height().value());
            default:
                throw std::runtime_error("Unknow type of wire");
        }
    }

    double WireWrapper::get_maximum_conducting_width() {
        switch (get_type()) {
            case WireType::LITZ:
                {
                    auto strand = resolve_strand();
                    if (!strand.get_conducting_diameter()) {
                        throw std::runtime_error("Missing conducting diameter in litz strand");
                    }
                    return resolve_dimensional_values(strand.get_conducting_diameter().value());
                }
            case WireType::ROUND:
                return resolve_dimensional_values(get_conducting_diameter().value());
            case WireType::RECTANGULAR:
            case WireType::FOIL:
                return resolve_dimensional_values(get_conducting_width().value());
            default:
                throw std::runtime_error("Unknow type of wire");
        }
    }

    double WireWrapper::get_maximum_conducting_height() {
        switch (get_type()) {
            case WireType::LITZ:
                {
                    auto strand = resolve_strand();
                    if (!strand.get_conducting_diameter()) {
                        throw std::runtime_error("Missing conducting diameter in litz strand");
                    }
                    return resolve_dimensional_values(strand.get_conducting_diameter().value());
                }
            case WireType::ROUND:
                return resolve_dimensional_values(get_conducting_diameter().value());
            case WireType::RECTANGULAR:
            case WireType::FOIL:
                return resolve_dimensional_values(get_conducting_height().value());
            default:
                throw std::runtime_error("Unknow type of wire");
        }
    }

    double WireWrapper::get_minimum_conducting_dimension() {
        switch (get_type()) {
            case WireType::LITZ:
                {
                    auto strand = resolve_strand();
                    if (!strand.get_conducting_diameter()) {
                        throw std::runtime_error("Missing conducting diameter in litz strand");
                    }
                    return resolve_dimensional_values(strand.get_conducting_diameter().value());
                }
            case WireType::ROUND:
                    return resolve_dimensional_values(get_conducting_diameter().value());
            case WireType::RECTANGULAR:
                    return resolve_dimensional_values(get_conducting_width().value());
            case WireType::FOIL:
                    return resolve_dimensional_values(get_conducting_height().value());
            default:
                throw std::runtime_error("Unknow type of wire");
        }
    }

    void WireWrapper::cut_foil_wire_to_section(Section section) {
        auto constants = Constants();
        if (get_type() != WireType::FOIL) {
            throw std::runtime_error("Method only valid for Foil wire");
        }
        DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_maximum(section.get_dimensions()[1] * (1 - constants.foilToSectionMargin));
        set_conducting_height(dimensionWithTolerance);
        set_outer_width(get_conducting_width());
        set_outer_height(get_conducting_height());
    }


} // namespace OpenMagnetics
