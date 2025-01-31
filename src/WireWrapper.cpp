#include "WireWrapper.h"
#include "WindingSkinEffectLosses.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <vector>
#include <algorithm>
#include "Constants.h"
#include "Defaults.h"
#include "Utils.h"
#include "spline.h"

std::map<std::string, tk::spline> wireCoatingThicknessProportionInterps;
std::map<std::string, tk::spline> wireFillingFactorInterps;
std::map<std::string, tk::spline> wirePackingFactorInterps;
std::map<std::string, tk::spline> wireConductingAreaProportionInterps;
std::map<std::string, double> minWireConductingDimensions;
std::map<std::string, double> maxWireConductingDimensions;
std::map<std::string, int64_t> minLitzWireNumberConductors;
std::map<std::string, int64_t> maxLitzWireNumberConductors;

namespace OpenMagnetics {

    WireRound WireWrapper::convert_from_wire_to_strand(WireWrapper wire) {
        WireRound strand;
        strand.set_type(wire.get_type());
        if (wire.get_conducting_diameter()) 
            strand.set_conducting_diameter(wire.get_conducting_diameter().value());
        if (wire.get_material()) 
            strand.set_material(wire.get_material().value());
        if (wire.get_outer_diameter()) 
            strand.set_outer_diameter(wire.get_outer_diameter().value());
        if (wire.get_coating()) 
            strand.set_coating(wire.get_coating().value());
        if (wire.get_conducting_area()) 
            strand.set_conducting_area(wire.get_conducting_area().value());
        if (wire.get_manufacturer_info()) 
            strand.set_manufacturer_info(wire.get_manufacturer_info().value());
        if (wire.get_name()) 
            strand.set_name(wire.get_name().value());
        if (wire.get_number_conductors()) 
            strand.set_number_conductors(wire.get_number_conductors().value());
        if (wire.get_standard()) 
            strand.set_standard(wire.get_standard().value());
        if (wire.get_standard_name()) 
            strand.set_standard_name(wire.get_standard_name().value());
        return strand;
    }

    std::optional<InsulationWireCoating> WireWrapper::resolve_coating(const WireWrapper& wire) {
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

    WireRound wire_to_wire_round(const WireWrapper& wire) {
        WireRound wireRound;
        wireRound.set_type(wire.get_type());
        if (wire.get_coating())
            wireRound.set_coating(wire.get_coating().value());
        wireRound.set_conducting_diameter(wire.get_conducting_diameter().value());
        if (wire.get_manufacturer_info())
            wireRound.set_manufacturer_info(wire.get_manufacturer_info().value());
        if (wire.get_material())
            wireRound.set_material(wire.get_material().value());
        if (wire.get_name())
            wireRound.set_name(wire.get_name().value());
        if (wire.get_number_conductors())
            wireRound.set_number_conductors(wire.get_number_conductors().value());
        if (wire.get_outer_diameter())
            wireRound.set_outer_diameter(wire.get_outer_diameter().value());
        if (wire.get_standard())
            wireRound.set_standard(wire.get_standard().value());
        if (wire.get_standard_name())
            wireRound.set_standard_name(wire.get_standard_name().value());
        if (wire.get_conducting_area())
            wireRound.set_conducting_area(wire.get_conducting_area().value());

        return wireRound;
    }

    WireRound WireWrapper::resolve_strand(const WireWrapper& wire) {
        if (!wire.get_strand())
            throw std::runtime_error("Litz wire is missing strand information");

        // If the strand is a string, we have to load its data from the database
        if (std::holds_alternative<std::string>(wire.get_strand().value())) {
            auto strand = find_wire_by_name(std::get<std::string>(wire.get_strand().value()));

            return wire_to_wire_round(strand);
        }
        else {
            return std::get<WireRound>(wire.get_strand().value());
        }

    }

    WireRound WireWrapper::resolve_strand() {
        if (!get_strand())
            throw std::runtime_error("Litz wire is missing strand information");

        // If the strand is a string, we have to load its data from the database
        if (std::holds_alternative<std::string>(get_strand().value())) {
            auto strand = find_wire_by_name(std::get<std::string>(get_strand().value()));

            return wire_to_wire_round(strand);
        }
        else {
            return std::get<WireRound>(get_strand().value());
        }

    }

    WireMaterial WireWrapper::resolve_material() { 
        if (get_type() == WireType::LITZ) {
            auto strand = resolve_strand();
            return resolve_material(strand);
        }

        if (!get_material()) {
            set_material(Defaults().defaultConductorMaterial);
            // throw std::runtime_error("Wire is missing material information");
        }

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
            return resolve_material(strand);
        }
        if (!wire.get_material()) {
            wire.set_material(Defaults().defaultConductorMaterial);
            // throw std::runtime_error("Wire is missing material information");
        }

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


    WireMaterial WireWrapper::resolve_material(WireRound wire) { 
        if (!wire.get_material()) {
            wire.set_material(Defaults().defaultConductorMaterial);
            // throw std::runtime_error("Wire is missing material information");
        }

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
            double wireConductingDimension, wireCoatingThicknessProportion, wireFillingFactor, wirePackingFactor;
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
                        wireConductingDimension = resolve_dimensional_values(strandWire.get_conducting_diameter());

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
                    double wireCoatingThickness = (wireOuterDimension - wireConductingDimension) / 2;
                    double wireCoatingThicknessProportion = wireCoatingThickness / wireConductingDimension;
                    InterpolatorDatum interpolatorDatum = { wireConductingDimension, wireCoatingThicknessProportion, wireFillingFactor, wirePackingFactor };

                    interpolatorData.push_back(interpolatorDatum);
                }
            }
        }
        if (interpolatorData.size() == 0) {
            // throw std::runtime_error("No wires with that specification");
            return;
        }

        size_t n = interpolatorData.size();
        std::vector<double> x, xPackingFactor, yFillingFactor, yCoatingThicknessProportion, yPackingFactor;
        std::sort(interpolatorData.begin(), interpolatorData.end(), [](const InterpolatorDatum& b1, const InterpolatorDatum& b2) {
            return b1.wireConductingDimension < b2.wireConductingDimension;
        });

        minWireConductingDimensions[key] = interpolatorData[0].wireConductingDimension;
        maxWireConductingDimensions[key] = interpolatorData[n - 1].wireConductingDimension;

        for (size_t i = 0; i < n; i++) {
            if (x.size() == 0 || interpolatorData[i].wireConductingDimension != x.back()) {
                x.push_back(interpolatorData[i].wireConductingDimension);
                yFillingFactor.push_back(interpolatorData[i].wireFillingFactor);
                yCoatingThicknessProportion.push_back(interpolatorData[i].wireCoatingThicknessProportion);
                if (wireType == WireType::LITZ && interpolatorData[i].wirePackingFactor > 0) {
                    xPackingFactor.push_back(interpolatorData[i].wireConductingDimension);
                    yPackingFactor.push_back(interpolatorData[i].wirePackingFactor);
                }
            }
        }
        tk::spline interpFillingFactor(x, yFillingFactor, tk::spline::cspline_hermite, true);
        tk::spline interpCoatingThicknessProportion(x, yCoatingThicknessProportion, tk::spline::cspline_hermite, true);
        tk::spline interpPackingFactor;
        if (wireType == WireType::LITZ && xPackingFactor.size() > 0) {
            interpPackingFactor = tk::spline(xPackingFactor, yPackingFactor, tk::spline::cspline_hermite, true);
        }
        wireFillingFactorInterps[key] = interpFillingFactor;
        wireCoatingThicknessProportionInterps[key] = interpCoatingThicknessProportion;
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
            // throw std::runtime_error("No wires with that specification");
            return;
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
            // throw std::runtime_error("No wires with that specification");
            return;
        }

        size_t n = interpolatorData.size();
        std::vector<double> x, y;
        std::sort(interpolatorData.begin(), interpolatorData.end(), [](const InterpolatorDatum& b1, const InterpolatorDatum& b2) {
            return b1.wireTheoreticalConductingArea < b2.wireTheoreticalConductingArea;
        });

        for (size_t i = 0; i < n; i++) {
            if (x.size() == 0 || interpolatorData[i].wireTheoreticalConductingArea != x.back()) {
                x.push_back(interpolatorData[i].wireTheoreticalConductingArea);
                y.push_back(interpolatorData[i].wireRealConductingArea / interpolatorData[i].wireTheoreticalConductingArea);
            }
        }
        tk::spline interpConductingAreaProportion;
        if (x.size() > 0) {
            interpConductingAreaProportion = tk::spline(x, y, tk::spline::cspline_hermite, true);
            wireConductingAreaProportionInterps[key] = interpConductingAreaProportion;
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
            load_wires();
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
            load_wires();
        }

        if (standard && (wireType == WireType::RECTANGULAR || wireType == WireType::FOIL || wireType == WireType::PLANAR)) {
            standard = std::nullopt;
        } 

        if (!wireCoatingThicknessProportionInterps.contains(key)) {
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
        double wireConductingDimensionForInterpolator = wireConductingDimension;
        wireConductingDimensionForInterpolator = std::max(wireConductingDimensionForInterpolator, minWireConductingDimensions[key]);
        wireConductingDimensionForInterpolator = std::min(wireConductingDimensionForInterpolator, maxWireConductingDimensions[key]);
        double coatingThickness = wireCoatingThicknessProportionInterps[key](wireConductingDimensionForInterpolator) * wireConductingDimension;
        double outerDimension = coatingThickness * 2 + wireConductingDimension;
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
            load_wires();
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
            load_wires();
        }

        if (!wireConductingAreaProportionInterps.contains(key)) {
            create_conducting_area_interpolator(standard,
                                                key);
        }
        double wireTheoreticalConductingArea = conductingWidth * conductingHeight;

        if (!wireConductingAreaProportionInterps.contains(key)) {
            return wireTheoreticalConductingArea;
        }
        else {

            double wireTheoreticalConductingArea = conductingWidth * conductingHeight;
            double conductingArea = wireConductingAreaProportionInterps[key](wireTheoreticalConductingArea) * wireTheoreticalConductingArea;
            if (conductingArea <= 0 || conductingArea < wireTheoreticalConductingArea / 2) {
                return wireTheoreticalConductingArea;
            }
            return conductingArea;
        }
    }

    double get_packing_factor_from_standard(WireStandard standard, double numberConductors) {
        if (standard == WireStandard::IEC_60317) {
            // Accoding to standard IEC 60317 - 11
            if (numberConductors == 2) {
                return sqrt(2);
            }
            else if (numberConductors < 12) {
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

    double WireWrapper::get_serving_thickness_from_standard(int numberLayers, double outerDiameter) {
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
        double outerDiameter = get_outer_diameter_bare_litz(conductingDiameter, numberConductors, grade, standard);

        double servingThickness = get_serving_thickness_from_standard(numberLayers, outerDiameter);

        outerDiameter += servingThickness * 2;

        return outerDiameter;
    }

    double WireWrapper::get_outer_diameter_bare_litz(double conductingDiameter, int numberConductors, int grade, WireStandard standard) {
        double packingFactor = get_packing_factor_from_standard(standard, numberConductors);
        double outerStrandDiameter = get_outer_diameter_round(conductingDiameter, grade, standard);

        double outerDiameter = packingFactor * sqrt(numberConductors) * outerStrandDiameter;

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
        if (packingFactor < 1 / 0.9069) {
            // if packing if smaller than what is physically possible (https://en.wikipedia.org/wiki/Circle_packing) we use the formula. 
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
                    auto conductingDiameter = resolve_dimensional_values(strand.get_conducting_diameter());
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
            case WireType::PLANAR:
                {
                    if (!get_conducting_width()) {
                        throw std::runtime_error("Missing conducting width in planar wire");
                    }
                    if (!get_conducting_height()) {
                        throw std::runtime_error("Missing conducting height in planar wire");
                    }
                    auto conductingWidth = resolve_dimensional_values(get_conducting_width().value());
                    auto conductingHeight = resolve_dimensional_values(get_conducting_height().value()) * get_number_conductors().value();
                    return conductingWidth * conductingHeight;
                }
            default:
                throw std::runtime_error("Unknow type of wire");
        }
    }

    std::vector<std::vector<double>> WireWrapper::calculate_current_density_distribution(SignalDescriptor current, double frequency, double temperature, size_t numberPoints) {
        auto material = resolve_material();
        if (!current.get_processed()) {
            throw std::runtime_error("Current is not processed");
        }

        if (!current.get_processed()->get_rms()) {
            throw std::runtime_error("Current is missing RMS");
        }
        auto skinDepth = WindingSkinEffectLosses::calculate_skin_depth(material, frequency, temperature);
        std::complex<double> waveNumber(1, -1);
        waveNumber /= skinDepth;

        if (get_type() == WireType::LITZ) {
            throw std::runtime_error("For LITZ use ROUND in the strands");
        }
        if (get_type() == WireType::ROUND) {
            std::vector<double> distributionRadial;
            auto conductingRadius = get_minimum_conducting_dimension() / 2;
            auto currentRms = current.get_processed()->get_rms().value();

            for (size_t pointIndex = 0; pointIndex < numberPoints; ++pointIndex) {
                double radius = std::max(conductingRadius / 1000, conductingRadius * pointIndex / (numberPoints - 1));
                double currentDensityPoint = abs(waveNumber * currentRms / (2 * std::numbers::pi * conductingRadius) * bessel_first_kind(0, waveNumber * radius) / bessel_first_kind(1, waveNumber * conductingRadius));
                distributionRadial.push_back(currentDensityPoint);
            }

            return {distributionRadial};
        }
        if (get_type() == WireType::RECTANGULAR || get_type() == WireType::FOIL || get_type() == WireType::PLANAR) {
            std::vector<double> distributionHorizontal;
            std::vector<double> distributionVertical;
            auto conductingSemiWidth = resolve_dimensional_values(get_conducting_width().value()) / 2;
            auto conductingSemiHeight = resolve_dimensional_values(get_conducting_height().value()) / 2;
            auto currentRms = current.get_processed()->get_rms().value();

            for (size_t pointIndex = 0; pointIndex < numberPoints; ++pointIndex) {
                double horizontalRadius = std::max(conductingSemiWidth / 1000, conductingSemiWidth * pointIndex / (numberPoints - 1));
                double verticalRadius = std::max(conductingSemiHeight / 1000, conductingSemiHeight * pointIndex / (numberPoints - 1));
                double horizontalCurrentDensityPoint = abs(waveNumber * currentRms / (2 * std::numbers::pi * conductingSemiWidth) * bessel_first_kind(0, waveNumber * horizontalRadius) / bessel_first_kind(1, waveNumber * conductingSemiWidth));
                double verticalCurrentDensityPoint = abs(waveNumber * currentRms / (2 * std::numbers::pi * conductingSemiHeight) * bessel_first_kind(0, waveNumber * verticalRadius) / bessel_first_kind(1, waveNumber * conductingSemiHeight));
                distributionHorizontal.push_back(horizontalCurrentDensityPoint);
                distributionVertical.push_back(verticalCurrentDensityPoint);
            }

            return {distributionHorizontal, distributionVertical};
        }
        else {
            throw std::runtime_error("Unknonw wire type");
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
                    conductingSmallestDimension = resolve_dimensional_values(strand.get_conducting_diameter());
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
            case WireType::PLANAR:
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
                        nonConductingArea = std::numbers::pi * pow(resolve_dimensional_values(strand.get_conducting_diameter()) / 2 - skinDepth, 2) * get_number_conductors().value();
                        break;
                    }
                case WireType::ROUND:
                    {
                        nonConductingArea = std::numbers::pi * pow(resolve_dimensional_values(get_conducting_diameter().value()) / 2 - skinDepth, 2) * get_number_conductors().value();
                        break;
                    }
                case WireType::PLANAR:
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
            case WireType::PLANAR:
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
            case WireType::PLANAR:
            case WireType::RECTANGULAR:
            case WireType::FOIL:
                if (get_outer_height())
                    return resolve_dimensional_values(get_outer_height().value());
                else if (get_conducting_height())
                    return resolve_dimensional_values(get_conducting_height().value());
                else
                    return 0;
            default:
                throw std::runtime_error("Unknow type of wire");
        }
    }


    double WireWrapper::get_maximum_outer_dimension() {
        return std::max(get_maximum_outer_width(), get_maximum_outer_height());
    }

    double WireWrapper::get_maximum_conducting_width() {
        switch (get_type()) {
            case WireType::LITZ:
                {
                    auto strand = resolve_strand();
                    return resolve_dimensional_values(strand.get_conducting_diameter());
                }
            case WireType::ROUND:
                return resolve_dimensional_values(get_conducting_diameter().value());
            case WireType::PLANAR:
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
                    return resolve_dimensional_values(strand.get_conducting_diameter());
                }
            case WireType::ROUND:
                return resolve_dimensional_values(get_conducting_diameter().value());
            case WireType::PLANAR:
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
                    return resolve_dimensional_values(strand.get_conducting_diameter());
                }
            case WireType::ROUND:
                    return resolve_dimensional_values(get_conducting_diameter().value());
            case WireType::PLANAR:
            case WireType::RECTANGULAR:
                    return resolve_dimensional_values(get_conducting_height().value());
            case WireType::FOIL:
                    return resolve_dimensional_values(get_conducting_width().value());
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

    void WireWrapper::cut_planar_wire_to_section(Section section) {
        auto constants = Constants();
        if (get_type() != WireType::PLANAR) {
            throw std::runtime_error("Method only valid for Planar wire");
        }
        DimensionWithTolerance dimensionWithTolerance;
        dimensionWithTolerance.set_maximum(section.get_dimensions()[0] * (1 - constants.planarToSectionMargin));
        set_conducting_width(dimensionWithTolerance);
        set_outer_width(get_conducting_width());
        set_outer_height(get_conducting_height());
    }

    double WireWrapper::get_relative_cost() {
        double cost = 0;
        if (get_type() == WireType::LITZ) {
            cost++;
            if (get_number_conductors().value() > 100) {
                cost++;
            }
            if (get_number_conductors().value() > 500) {
                cost++;
            }
            if (get_number_conductors().value() > 1000) {
                cost++;
            }
        }
        auto coatingMaybe = resolve_coating();
        if (coatingMaybe) {
            auto coating = coatingMaybe.value();
    
            if (coating.get_type() == InsulationWireCoatingType::INSULATED) {
                cost += coating.get_number_layers().value();
    
                // Due to lead time
                if (coating.get_number_layers().value() == 1 || coating.get_number_layers().value() == 2) {
                    cost += 2;
                }
            }
    
            if (coating.get_type() == InsulationWireCoatingType::ENAMELLED) {
                cost += coating.get_grade().value();
            }
        }

        return cost;
    }


    double WireWrapper::get_coating_thickness() {
        return get_coating_thickness(*this);
    }

    double WireWrapper::get_coating_dielectric_strength() {
        return get_coating_dielectric_strength(*this);
    }

    double WireWrapper::get_coating_thickness(WireWrapper wire) {
        auto coating = resolve_coating(wire);

        if (coating->get_thickness()) {
            return resolve_dimensional_values(coating->get_thickness().value());
        }

        if (coating->get_thickness_layers() && coating->get_number_layers()) {
            return coating->get_thickness_layers().value() * coating->get_number_layers().value();
        }

        auto maximumOuterWidth = wire.get_maximum_outer_width();
        auto maximumOuterHeight = wire.get_maximum_outer_height();
        auto maximumConductingWidth = wire.get_maximum_conducting_width();
        auto maximumConductingHeight = wire.get_maximum_conducting_height();

        auto coatingThicknessWidth = (maximumOuterWidth - maximumConductingWidth) / 2;
        auto coatingThicknessHeight = (maximumOuterHeight - maximumConductingHeight) / 2;

        return std::min(coatingThicknessWidth, coatingThicknessHeight);
    }

    double WireWrapper::get_coating_dielectric_strength(WireWrapper wire) {
        auto coatingInsulationMaterial = resolve_coating_insulation_material(wire);
        auto coatingThickness = get_coating_thickness(wire);
        return coatingInsulationMaterial.get_dielectric_strength_by_thickness(coatingThickness);
    }

    double WireWrapper::get_coating_relative_permittivity() {
        return get_coating_relative_permittivity(*this);
    }

    double WireWrapper::get_coating_relative_permittivity(WireWrapper wire) {
        auto coatingInsulationMaterial = resolve_coating_insulation_material(wire);
        if (!coatingInsulationMaterial.get_relative_permittivity())
            throw std::runtime_error("Coating insulation material is missing dielectric constant");
        return coatingInsulationMaterial.get_relative_permittivity().value();
    }

    InsulationMaterialWrapper WireWrapper::resolve_coating_insulation_material(){
        return resolve_coating_insulation_material(*this);
    }

    InsulationMaterialWrapper WireWrapper::resolve_coating_insulation_material(WireWrapper wire) {
        auto coating = resolve_coating(wire);

        if (!coating->get_material())
        {
            if (coating->get_type().value() == InsulationWireCoatingType::ENAMELLED) {
                coating->set_material(Defaults().defaultEnamelledInsulationMaterial);
            }
            else {
                throw std::runtime_error("Coating is missing material information");
            }
        }

        auto insulationMaterial = coating->get_material().value();
        // If the material is a string, we have to load its data from the database
        if (std::holds_alternative<std::string>(insulationMaterial)) {
            auto insulationMaterialData = find_insulation_material_by_name(std::get<std::string>(insulationMaterial));

            return insulationMaterialData;
        }
        else {
            return std::get<InsulationMaterial>(insulationMaterial);
        }
    }

    InsulationMaterialWrapper WireWrapper::resolve_coating_insulation_material(WireRound wire) {
        auto coating = resolve_coating(wire);

        if (!coating->get_material())
        {
            if (coating->get_type().value() == InsulationWireCoatingType::ENAMELLED) {
                coating->set_material(Defaults().defaultEnamelledInsulationMaterial);
            }
            else {
                throw std::runtime_error("Coating is missing material information");
            }
        }

        auto insulationMaterial = coating->get_material().value();
        // If the material is a string, we have to load its data from the database
        if (std::holds_alternative<std::string>(insulationMaterial)) {
            auto insulationMaterialData = find_insulation_material_by_name(std::get<std::string>(insulationMaterial));

            return insulationMaterialData;
        }
        else {
            return std::get<InsulationMaterial>(insulationMaterial);
        }
    }

    std::string WireWrapper::encode_coating_label() {
        return encode_coating_label(*this);
    }

    std::string WireWrapper::encode_coating_label(WireWrapper wire) {
        auto coating = resolve_coating(wire);
        if (!coating) {
            return "Bare";
        }
        std::string text = "";
        if (coating->get_type().value() == InsulationWireCoatingType::INSULATED) {
            if (!coating->get_number_layers()) {
                throw std::invalid_argument("Missing number of layers in insulated wire");
            }
            if (coating->get_number_layers().value() == 1) {
                text += "SIW";
            }
            else if (coating->get_number_layers().value() == 2) {
                text += "DIW";
            }
            else if (coating->get_number_layers().value() == 3) {
                text += "TIW";
            }

            if (coating->get_temperature_rating()) {
                std::ostringstream auxStream;
                auxStream << std::fixed << std::setprecision(0) << coating->get_temperature_rating().value();
                text += ", TR " + auxStream.str();
            }

            if (coating->get_breakdown_voltage()) {
                if (coating->get_breakdown_voltage().value() > 1000) {
                    std::ostringstream auxStream;
                    auxStream << std::fixed << std::setprecision(1) << (coating->get_breakdown_voltage().value() / 1000);
                    text += ", BV " + auxStream.str() + " kV";
                }
                else {
                    std::ostringstream auxStream;
                    auxStream << std::fixed << std::setprecision(0) << coating->get_breakdown_voltage().value();
                    text += ", BV " + auxStream.str() + " V";
                }
            }
        }
        else if (coating->get_type().value() == InsulationWireCoatingType::ENAMELLED) {
            if (!coating->get_grade()) {
                throw std::invalid_argument("Missing grade in enamelled wire");
            }
            text += "Enamel";
            WireStandard standard;
            if (wire.get_standard()) {
                standard = wire.get_standard().value();
            }
            else {
                standard = WireStandard::IEC_60317;
            }

            if (coating->get_grade().value() == 1) {
                if (standard == WireStandard::NEMA_MW_1000_C) {
                    text += " single build";
                }
                else {
                    text += " grade 1";
                }
            }
            else if (coating->get_grade().value() == 2) {
                if (standard == WireStandard::NEMA_MW_1000_C) {
                    text += " heavy build";
                }
                else {
                    text += " grade 2";
                }
            }
            else if (coating->get_grade().value() == 3) {
                if (standard == WireStandard::NEMA_MW_1000_C) {
                    text += " triple build";
                }
                else {
                    text += " grade 3";
                }
            }
            else {
                text += " FIW" + std::to_string(coating->get_grade().value());
            }
        }
        else if (coating->get_type().value() == InsulationWireCoatingType::SERVED) {
            if (!coating->get_number_layers()) {
                throw std::invalid_argument("Missing number of layers in served wire: " + wire.get_name().value());
            }
            if (coating->get_number_layers().value() == 1) {
                text += "Single served";
            }
            else if (coating->get_number_layers().value() == 2) {
                text += "Double served";
            }
        }
        else if (coating->get_type().value() == InsulationWireCoatingType::BARE) {
            if (wire.get_type() == WireType::LITZ) {
                text += "Unserved";
            }
            else {
                text += "Bare";
            }
        }
        else {
            throw std::invalid_argument("Coating type not implemented yet: " + std::string{magic_enum::enum_name(coating->get_type().value())});
        }

        return text;
    }

    std::optional<InsulationWireCoating> WireWrapper::decode_coating_label(std::string label) {
        std::optional<InsulationWireCoating> coating;

        if ((label.find("Bare") != std::string::npos) || (label.find("Unserved") != std::string::npos)) {
            coating->set_type(InsulationWireCoatingType::BARE);
        }
        if (label.find("SIW") != std::string::npos) {
            coating->set_type(InsulationWireCoatingType::INSULATED);
            coating->set_number_layers(1);
        }
        if (label.find("DIW") != std::string::npos) {
            coating->set_type(InsulationWireCoatingType::INSULATED);
            coating->set_number_layers(2);
        }
        if (label.find("TIW") != std::string::npos) {
            coating->set_type(InsulationWireCoatingType::INSULATED);
            coating->set_number_layers(3);
        }

        if (label.find("TR 155") != std::string::npos) {
            coating->set_temperature_rating(155);
        }
        if (label.find("TR 180") != std::string::npos) {
            coating->set_temperature_rating(180);
        }

        if (label.find(", BV ") != std::string::npos) {
            auto aux = split(label,", BV ")[1];
            if (aux.find(" kV") != std::string::npos) {
                auto breakdownVoltage = std::stod(split(aux," kV")[0]);
                coating->set_breakdown_voltage(breakdownVoltage * 1000);
            }
            if (aux.find(" V") != std::string::npos) {
                auto breakdownVoltage = std::stod(split(aux," V")[0]);
                coating->set_breakdown_voltage(breakdownVoltage);
            }
        }

        if (label.find("Enamel") != std::string::npos) {
            coating->set_type(InsulationWireCoatingType::ENAMELLED);
            if ((label.find("grade 1") != std::string::npos) || (label.find("single build") != std::string::npos)) {
                coating->set_grade(1);
            }
            else if ((label.find("grade 2") != std::string::npos) || (label.find("heavy build") != std::string::npos)) {
                coating->set_grade(2);
            }
            else if ((label.find("grade 3") != std::string::npos) || (label.find("triple build") != std::string::npos)) {
                coating->set_grade(3);
            }
            else {
                auto aux = split(label,"FIW")[1];
                auto grade = std::stoi(aux);
                coating->set_grade(grade);
            }
        }
        if (label.find("Single served") != std::string::npos) {
            coating->set_type(InsulationWireCoatingType::SERVED);
            coating->set_number_layers(1);
        }
        if (label.find("Double served") != std::string::npos) {
            coating->set_type(InsulationWireCoatingType::SERVED);
            coating->set_number_layers(2);
        }

        return coating;
    }
    

    WireWrapper WireWrapper::get_wire_for_frequency(double effectiveFrequency, double temperature, bool exact) {
        auto skinDepth = WindingSkinEffectLosses::calculate_skin_depth("copper", effectiveFrequency, temperature);
        double wireConductingDiameter = skinDepth * 2;
        if (exact) {
            WireWrapper wire;
            wire.set_nominal_value_conducting_diameter(wireConductingDiameter);

            wire.set_nominal_value_outer_diameter(get_outer_diameter_round(wireConductingDiameter, 1, WireStandard::IEC_60317));
            wire.set_material("copper");
            wire.set_type(WireType::ROUND);
            return wire;
        }
        else {
            WireWrapper wire = find_wire_by_dimension(wireConductingDiameter, WireType::ROUND, WireStandard::IEC_60317, false);
            return wire;
        }
    }
    

    WireWrapper WireWrapper::get_wire_for_conducting_area(double conductingArea, double temperature, bool exact) {
        double wireConductingDiameter = sqrt(4 * conductingArea / std::numbers::pi);
        if (exact) {
            WireWrapper wire;
            wire.set_nominal_value_conducting_diameter(wireConductingDiameter);
            wire.set_nominal_value_outer_diameter(get_outer_diameter_round(wireConductingDiameter, 1, WireStandard::IEC_60317));
            wire.set_material("copper");
            wire.set_type(WireType::ROUND);
            return wire;
        }
        else {
            WireWrapper wire = find_wire_by_dimension(wireConductingDiameter, WireType::ROUND, WireStandard::IEC_60317, false);
            return wire;
        }
    }


    WireWrapper WireWrapper::get_equivalent_wire(WireWrapper oldWire, WireType newWireType, double effectiveFrequency, double temperature) {

        WireStandard standard;
        if (oldWire.get_standard()) {
            standard = oldWire.get_standard().value();
        }
        else {
            standard = WireStandard::IEC_60317;
        }

        switch (newWireType) {
            case WireType::LITZ:
                {
                    double strandConductingDiameter = WindingSkinEffectLosses::calculate_skin_depth("copper", effectiveFrequency, temperature);
                    double strandConductingArea = std::numbers::pi * pow(strandConductingDiameter / 2, 2);
                    WireWrapper strand;
                    size_t numberConductors;

                    switch (oldWire.get_type()) {
                        case WireType::LITZ:
                            return oldWire;
                            break;
                        default:
                            double oldConductingArea = oldWire.calculate_conducting_area();
                            strand = OpenMagnetics::find_wire_by_dimension(strandConductingDiameter, OpenMagnetics::WireType::ROUND, standard);
                            numberConductors = round(oldConductingArea / strandConductingArea);
                            break;
                    }


                    InsulationWireCoating newCoating;
                    auto oldCoating = oldWire.resolve_coating();

                    // Served coating by default
                    newCoating.set_type(InsulationWireCoatingType::SERVED);
                    newCoating.set_number_layers(1);

                    // We get the old coating only if it's insulated
                    if (oldCoating) {
                        if (oldCoating->get_type() == InsulationWireCoatingType::INSULATED) {
                            newCoating = oldCoating.value();
                        }
                    }

                    WireWrapper newWire;
                    newWire.set_type(WireType::LITZ);
                    newWire.set_strand(convert_from_wire_to_strand(strand));
                    newWire.set_coating(newCoating);
                    newWire.set_number_conductors(numberConductors);

                    if (newCoating.get_type() == InsulationWireCoatingType::SERVED) {
                        double outerDiameter = get_outer_diameter_served_litz(strandConductingDiameter, numberConductors, 1, 1, standard);
                        newWire.set_nominal_value_outer_diameter(outerDiameter);
                    }
                    else if (newCoating.get_type() == InsulationWireCoatingType::INSULATED) {
                        if (!newCoating.get_number_layers()) {
                            throw std::runtime_error("Missing number of layer in insulated wire coating");
                        }
                        if (!newCoating.get_thickness_layers()) {
                            throw std::runtime_error("Missing layer thickness in insulated wire coating");
                        }
                        double outerDiameter = get_outer_diameter_insulated_litz(strandConductingDiameter, numberConductors, newCoating.get_number_layers().value(), newCoating.get_thickness_layers().value(), 1, standard);
                        newWire.set_nominal_value_outer_diameter(outerDiameter);
                    }

                    return newWire;
                }
            case WireType::ROUND:
                {
                    double conductingDiameter;

                    switch (oldWire.get_type()) {
                        case WireType::ROUND:
                            return oldWire;
                        case WireType::LITZ:
                            {
                                double oldConductingArea = oldWire.calculate_conducting_area();
                                conductingDiameter = sqrt(oldConductingArea / std::numbers::pi) * 2;
                                break;
                            }
                        case WireType::PLANAR:
                        case WireType::RECTANGULAR:
                        case WireType::FOIL:
                            conductingDiameter = oldWire.get_minimum_conducting_dimension();
                            break;
                        default:
                            throw std::runtime_error("Unknown wire type");
                    }
                    auto newWire = OpenMagnetics::find_wire_by_dimension(conductingDiameter, OpenMagnetics::WireType::ROUND, standard);
                    auto oldCoating = oldWire.resolve_coating();
                    newWire.set_type(WireType::ROUND);
                    newWire.set_coating(oldCoating);
                    return newWire;
                }
            case WireType::PLANAR:
                {
                    double conductingHeight;

                    switch (oldWire.get_type()) {
                        case WireType::PLANAR:
                            return oldWire;
                        case WireType::LITZ:
                            {
                                double oldConductingArea = oldWire.calculate_conducting_area();
                                conductingHeight = sqrt(oldConductingArea / std::numbers::pi) * 2;
                                break;
                            }
                        case WireType::RECTANGULAR:
                        case WireType::ROUND:
                        case WireType::FOIL:
                            conductingHeight = oldWire.get_minimum_conducting_dimension();
                            break;
                        default:
                            throw std::runtime_error("Unknown wire type");
                    }
                    auto newWire = OpenMagnetics::find_wire_by_dimension(conductingHeight, OpenMagnetics::WireType::PLANAR);

                    // No coating by default

                    newWire.set_type(WireType::PLANAR);
                    newWire.set_nominal_value_outer_height(resolve_dimensional_values(newWire.get_conducting_height().value()));
                    return newWire;
                }
            case WireType::RECTANGULAR:
                {
                    double conductingHeight;

                    switch (oldWire.get_type()) {
                        case WireType::RECTANGULAR:
                            return oldWire;
                        case WireType::LITZ:
                            {
                                double oldConductingArea = oldWire.calculate_conducting_area();
                                conductingHeight = sqrt(oldConductingArea / std::numbers::pi) * 2;
                                break;
                            }
                        case WireType::PLANAR:
                        case WireType::ROUND:
                        case WireType::FOIL:
                            conductingHeight = oldWire.get_minimum_conducting_dimension();
                            break;
                        default:
                            throw std::runtime_error("Unknown wire type");
                    }
                    auto newWire = OpenMagnetics::find_wire_by_dimension(conductingHeight, OpenMagnetics::WireType::RECTANGULAR);

                    InsulationWireCoating newCoating;
                    auto oldCoating = oldWire.resolve_coating();

                    // Enamelled coating by default
                    newCoating.set_type(InsulationWireCoatingType::ENAMELLED);
                    newCoating.set_grade(1);

                    // We get the old coating only if it's insulated
                    if (oldCoating) {
                        if (oldCoating->get_type() == InsulationWireCoatingType::ENAMELLED) {
                            newCoating = oldCoating.value();
                        }
                    }
                    newWire.set_coating(newCoating);
                    auto outerWidth = get_outer_width_rectangular(resolve_dimensional_values(newWire.get_conducting_width().value()), newCoating.get_grade().value());
                    auto outerHeight = get_outer_height_rectangular(resolve_dimensional_values(newWire.get_conducting_height().value()), newCoating.get_grade().value());

                    newWire.set_type(WireType::RECTANGULAR);
                    newWire.set_coating(newCoating);
                    newWire.set_nominal_value_outer_width(outerWidth);
                    newWire.set_nominal_value_outer_height(outerHeight);
                    return newWire;
                }
            case WireType::FOIL:
                {
                    double conductingWidth;

                    switch (oldWire.get_type()) {
                        case WireType::FOIL:
                            return oldWire;
                        case WireType::LITZ:
                            {
                                double oldConductingArea = oldWire.calculate_conducting_area();
                                conductingWidth = sqrt(oldConductingArea / std::numbers::pi) * 2;
                                break;
                            }
                        case WireType::ROUND:
                        case WireType::RECTANGULAR:
                            conductingWidth = oldWire.get_minimum_conducting_dimension();
                            break;
                        default:
                            throw std::runtime_error("Unknown wire type");
                    }
                    auto newWire = OpenMagnetics::find_wire_by_dimension(conductingWidth, OpenMagnetics::WireType::FOIL);

                    // No coating by default

                    newWire.set_type(WireType::FOIL);
                    newWire.set_nominal_value_outer_width(resolve_dimensional_values(newWire.get_conducting_width().value()));
                    return newWire;
                }
            default:
                throw std::runtime_error("Unknow type of wire");
        }
    }




} // namespace OpenMagnetics
