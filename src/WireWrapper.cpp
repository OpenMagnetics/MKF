#include <cmath>
#include <filesystem>
#include <fstream>
#include <vector>
#include <algorithm>
#include "Constants.h"
#include "Utils.h"
#include "WireWrapper.h"
#include "spline.h"

std::map<std::string, tk::spline> wireFillingFactorInterps;
std::map<std::string, double> minWireConductingWidths;
std::map<std::string, double> maxWireConductingWidths;

namespace OpenMagnetics {
    std::optional<InsulationWireCoating> WireWrapper::get_coating(WireS wire) {
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

    WireWrapper WireWrapper::get_strand(WireS wire) {
        // If the strand is a string, we have to load its data from the database
        if (std::holds_alternative<std::string>(wire.get_strand().value())) {
            throw std::runtime_error("Litz database not implemented yet");
        }
        else {
            return WireWrapper(std::get<WireSolid>(wire.get_strand().value()));
        }

    }

    double WireWrapper::get_filling_factor(double conductingDiameter, int grade, WireStandard standard, bool includeAirInCell){
        if (wireDatabase.empty()) {
            load_databases(true);
        }

        std::string standarString = static_cast<std::string>(magic_enum::enum_name(standard));
        std::string key = std::to_string(grade) + " " + standarString;

        if (!wireFillingFactorInterps.contains(key)) {
            std::vector<double> wireConductingDiameters;
            std::vector<double> fillingFactors;

            for (auto& datum : wireDatabase) {
                if (get_coating(datum.second)) {
                    auto coating = get_coating(datum.second).value();
                    if (coating.get_grade().value() == grade && datum.second.get_standard().value() == standard) {
                        double wireOuterDiameter = resolve_dimensional_values(datum.second.get_outer_diameter().value()); 
                        double wireConductingDiameter = resolve_dimensional_values(datum.second.get_conducting_diameter().value()); 
                        double wireNumberConductors = datum.second.get_number_conductors().value(); 
                        double outerArea;
                        if (includeAirInCell) {
                            outerArea = pow(wireOuterDiameter, 2);
                        } else {
                            outerArea = std::numbers::pi * pow(wireOuterDiameter / 2, 2);
                        }
                        double conductiveArea = std::numbers::pi * pow(wireConductingDiameter / 2, 2) * wireNumberConductors;
                        double wireFillingFactor = conductiveArea / outerArea;
                        wireConductingDiameters.push_back(wireConductingDiameter);
                        fillingFactors.push_back(wireFillingFactor);
                    }
                }
                else {
                    double wireConductingDiameter = resolve_dimensional_values(datum.second.get_conducting_diameter().value()); 
                    wireConductingDiameters.push_back(wireConductingDiameter);
                    fillingFactors.push_back(1);
                }
            }

            size_t n = wireConductingDiameters.size();
            std::vector<double> x, y;

            if (standard == WireStandard::NEMA_MW_1000_C) {
                // Invert the vector so interpolation works fine
                std::reverse(wireConductingDiameters.begin(), wireConductingDiameters.end());
                std::reverse(fillingFactors.begin(), fillingFactors.end());
            }
            minWireConductingWidths[key] = wireConductingDiameters[0];
            maxWireConductingWidths[key] = wireConductingDiameters[n - 1];

            for (size_t i = 0; i < n; i++) {
                if (x.size() == 0 || wireConductingDiameters[i] != x.back()) {
                    x.push_back(wireConductingDiameters[i]);
                    y.push_back(fillingFactors[i]);
                }
            }
            tk::spline interp(x, y, tk::spline::cspline_hermite, true);
            wireFillingFactorInterps[key] = interp;
        }


        conductingDiameter = std::max(conductingDiameter, minWireConductingWidths[key]);
        conductingDiameter = std::min(conductingDiameter, maxWireConductingWidths[key]);

        double fillingFactor = wireFillingFactorInterps[key](conductingDiameter);
        return fillingFactor;
    }

    int WireWrapper::get_equivalent_insulation_layers(double voltageToInsulate) {
        if (!get_coating(*this)) {
            return 0;
        }
        auto coating = get_coating(*this).value();
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

} // namespace OpenMagnetics
