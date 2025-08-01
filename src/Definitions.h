#pragma once
#include <MAS.hpp>
#include "json.hpp"

using json = nlohmann::json;

namespace OpenMagnetics {

class missing_material_data_exception : public std::exception {
private:
    std::string message;

public:
    missing_material_data_exception(std::string msg)
        : message(msg)
    {
    }

    const char* what() const throw()
    {
        return message.c_str();
    }
};

enum class DimensionalValues : int {
    MAXIMUM,
    NOMINAL,
    MINIMUM
};

enum class GappingType : int {
    GROUND,
    SPACER,
    RESIDUAL,
    DISTRIBUTED
};

enum class OrderedIsolationSide : int { 
    PRIMARY,
    SECONDARY,
    TERTIARY,
    QUATERNARY,
    QUINARY,
    SENARY,
    SEPTENARY,
    OCTONARY,
    NONARY,
    DENARY,
    UNDENARY,
    DUODENARY
};

enum class PainterModes : int {
    CONTOUR,
    QUIVER,
    SCATTER,
};


class Curve2D {

    std::vector<double> xPoints;
    std::vector<double> yPoints;
    std::string title;

    public: 
        Curve2D() = default;
        virtual ~Curve2D() = default;
        Curve2D(std::vector<double> xPoints, std::vector<double> yPoints, std::string title) : xPoints(xPoints), yPoints(yPoints), title(title){}

        const std::vector<double> & get_x_points() const { return xPoints; }
        std::vector<double> & get_mutable_x_points() { return xPoints; }
        void set_x_points(const std::vector<double> & value) { this->xPoints = value; }

        const std::vector<double> & get_y_points() const { return yPoints; }
        std::vector<double> & get_mutable_y_points() { return yPoints; }
        void set_y_points(const std::vector<double> & value) { this->yPoints = value; }

        const std::string & get_title() const { return title; }
        std::string & get_mutable_title() { return title; }
        void set_title(const std::string & value) { this->title = value; }
};
                
void from_json(const json & j, Curve2D & x);
void to_json(json & j, const Curve2D & x);

inline void from_json(const json & j, Curve2D& x) {
    x.set_x_points(j.at("xPoints").get<std::vector<double>>());
    x.set_y_points(j.at("yPoints").get<std::vector<double>>());
    x.set_title(j.at("title").get<std::string>());
}

inline void to_json(json & j, const Curve2D & x) {
    j = json::object();
    j["xPoints"] = x.get_x_points();
    j["yPoints"] = x.get_y_points();
    j["title"] = x.get_title();
}

enum class MagneticFilters : int {
    AREA_PRODUCT,
    ENERGY_STORED,
    ESTIMATED_COST,
    COST,
    CORE_AND_DC_LOSSES,
    LOSSES,
    DIMENSIONS,
    CORE_MINIMUM_IMPEDANCE,
    AREA_NO_PARALLELS,
    AREA_WITH_PARALLELS,
    EFFECTIVE_RESISTANCE,
    PROXIMITY_FACTOR,
    SOLID_INSULATION_REQUIREMENTS,
    TURNS_RATIOS,
    MAXIMUM_DIMENSIONS,
    SATURATION,
    DC_CURRENT_DENSITY,
    EFFECTIVE_CURRENT_DENSITY,
    IMPEDANCE,
    MAGNETIZING_INDUCTANCE,
    FRINGING_FACTOR,
    VOLUME,
    AREA,
    HEIGHT,
    TEMPERATURE_RISE,
    LOSSES_TIMES_VOLUME,
    VOLUME_TIMES_TEMPERATURE_RISE,
    LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE
};

class MagneticFilterOperation {
    public:
    MagneticFilterOperation(MagneticFilters filter, bool invert, bool log, double weight) : filter(filter), invert(invert), log(log), weight(weight) {};
    MagneticFilterOperation(MagneticFilters filter, bool invert, bool log, bool strictlyRequired, double weight) : filter(filter), invert(invert), log(log), strictlyRequired(strictlyRequired), weight(weight) {};
    MagneticFilterOperation() {};
    virtual ~MagneticFilterOperation() = default;

    private:
    MagneticFilters filter = MagneticFilters::DIMENSIONS;
    bool invert = true;
    bool log = false;
    bool strictlyRequired = false;
    double weight = 1;

    public:

    MagneticFilters get_filter() const { return filter; }
    void set_filter(const MagneticFilters & value) { this->filter = value; }

    bool get_invert() const { return invert; }
    void set_invert(const bool & value) { this->invert = value; }

    bool get_log() const { return log; }
    void set_log(const bool & value) { this->log = value; }

    bool get_strictly_required() const { return strictlyRequired; }
    void set_strictly_required(const bool & value) { this->strictlyRequired = value; }

    double get_weight() const { return weight; }
    void set_weight(const double & value) { this->weight = value; }
};


void from_json(const json & j, MagneticFilters & x);
void to_json(json & j, const MagneticFilters & x);
void from_json(const json& j, MagneticFilterOperation& x);
void to_json(json& j, const MagneticFilterOperation& x);
void from_json(const json& j, std::vector<MagneticFilterOperation>& v);
void to_json(json& j, const std::vector<MagneticFilterOperation>& v);


inline void from_json(const json & j, MagneticFilters & x) {
    if (j == "Area Product") x = MagneticFilters::AREA_PRODUCT;
    else if (j == "Energy Stored") x = MagneticFilters::ENERGY_STORED;
    else if (j == "Cost") x = MagneticFilters::COST;
    else if (j == "Estimated Cost") x = MagneticFilters::ESTIMATED_COST;
    else if (j == "Core And DC Losses") x = MagneticFilters::CORE_AND_DC_LOSSES;
    else if (j == "Losses") x = MagneticFilters::LOSSES;
    else if (j == "Dimensions") x = MagneticFilters::DIMENSIONS;
    else if (j == "Core Minimum Impedance") x = MagneticFilters::CORE_MINIMUM_IMPEDANCE;
    else if (j == "Area No Parallels") x = MagneticFilters::AREA_NO_PARALLELS;
    else if (j == "Area With Parallels") x = MagneticFilters::AREA_WITH_PARALLELS;
    else if (j == "Effective Resistance") x = MagneticFilters::EFFECTIVE_RESISTANCE;
    else if (j == "Proximity Factor") x = MagneticFilters::PROXIMITY_FACTOR;
    else if (j == "Solid Insulation Requirements") x = MagneticFilters::SOLID_INSULATION_REQUIREMENTS;
    else if (j == "Turns Ratios") x = MagneticFilters::TURNS_RATIOS;
    else if (j == "Maximum Dimensions") x = MagneticFilters::MAXIMUM_DIMENSIONS;
    else if (j == "Saturation") x = MagneticFilters::SATURATION;
    else if (j == "Dc Current Density") x = MagneticFilters::DC_CURRENT_DENSITY;
    else if (j == "Effective Current Density") x = MagneticFilters::EFFECTIVE_CURRENT_DENSITY;
    else if (j == "Impedance") x = MagneticFilters::IMPEDANCE;
    else if (j == "Magnetizing Inductance") x = MagneticFilters::MAGNETIZING_INDUCTANCE;
    else if (j == "Fringing Factor") x = MagneticFilters::FRINGING_FACTOR;
    else if (j == "Volume") x = MagneticFilters::VOLUME;
    else if (j == "Area") x = MagneticFilters::AREA;
    else if (j == "Height") x = MagneticFilters::HEIGHT;
    else if (j == "Temperature Rise") x = MagneticFilters::TEMPERATURE_RISE;
    else if (j == "Losses Times Volume") x = MagneticFilters::LOSSES_TIMES_VOLUME;
    else if (j == "Volume Times Temperature_Rise") x = MagneticFilters::VOLUME_TIMES_TEMPERATURE_RISE;
    else if (j == "Losses Times Volume Times Temperature_Rise") x = MagneticFilters::LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE;
    else { throw std::runtime_error("Input JSON does not conform to MagneticFilters schema!"); }
}

inline void to_json(json & j, const MagneticFilters & x) {
    switch (x) {
        case MagneticFilters::AREA_PRODUCT: j = "Area Product"; break;
        case MagneticFilters::ENERGY_STORED: j = "Energy Stored"; break;
        case MagneticFilters::COST: j = "Cost"; break;
        case MagneticFilters::ESTIMATED_COST: j = "Estimated Cost"; break;
        case MagneticFilters::CORE_AND_DC_LOSSES: j = "Core And DC Losses"; break;
        case MagneticFilters::LOSSES: j = "Losses"; break;
        case MagneticFilters::DIMENSIONS: j = "Dimensions"; break;
        case MagneticFilters::CORE_MINIMUM_IMPEDANCE: j = "Core Minimum Impedance"; break;
        case MagneticFilters::AREA_NO_PARALLELS: j = "Area No Parallels"; break;
        case MagneticFilters::AREA_WITH_PARALLELS: j = "Area With Parallels"; break;
        case MagneticFilters::EFFECTIVE_RESISTANCE: j = "Effective Resistance"; break;
        case MagneticFilters::PROXIMITY_FACTOR: j = "Proximity Factor"; break;
        case MagneticFilters::SOLID_INSULATION_REQUIREMENTS: j = "Solid Insulation Requirements"; break;
        case MagneticFilters::TURNS_RATIOS: j = "Turns Ratios"; break;
        case MagneticFilters::MAXIMUM_DIMENSIONS: j = "Maximum Dimensions"; break;
        case MagneticFilters::SATURATION: j = "Saturation"; break;
        case MagneticFilters::DC_CURRENT_DENSITY: j = "Dc Current Density"; break;
        case MagneticFilters::EFFECTIVE_CURRENT_DENSITY: j = "Effective Current Density"; break;
        case MagneticFilters::IMPEDANCE: j = "Impedance"; break;
        case MagneticFilters::MAGNETIZING_INDUCTANCE: j = "Magnetizing Inductance"; break;
        case MagneticFilters::FRINGING_FACTOR: j = "Fringing Factor"; break;
        case MagneticFilters::VOLUME: j = "Volume"; break;
        case MagneticFilters::AREA: j = "Area"; break;
        case MagneticFilters::HEIGHT: j = "Height"; break;
        case MagneticFilters::TEMPERATURE_RISE: j = "Temperature Rise"; break;
        case MagneticFilters::LOSSES_TIMES_VOLUME: j = "Losses Times Volume"; break;
        case MagneticFilters::VOLUME_TIMES_TEMPERATURE_RISE: j = "Volume Times Temperature_Rise"; break;
        case MagneticFilters::LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE: j = "Losses Times Volume Times Temperature_Rise"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"[object Object]\": " + std::to_string(static_cast<int>(x)));
    }
}

inline void from_json(const json& j, MagneticFilterOperation& x) {
    x.set_filter(j.at("filter").get<MagneticFilters>());
    x.set_invert(j.at("invert").get<bool>());
    x.set_log(j.at("log").get<bool>());
    x.set_strictly_required(j.at("strictlyRequired").get<bool>());
    x.set_weight(j.at("weight").get<double>());
}

inline void to_json(json& j, const MagneticFilterOperation& x) {
    j = json::object();
    j["filter"] = x.get_filter();
    j["invert"] = x.get_invert();
    j["log"] = x.get_log();
    j["strictlyRequired"] = x.get_strictly_required();
    j["weight"] = x.get_weight();
}

inline void from_json(const json& j, std::vector<MagneticFilterOperation>& v) {
    for (auto e : j) {
        MagneticFilterOperation x(
            e.at("filter").get<MagneticFilters>(),
            e.at("invert").get<bool>(),
            e.at("log").get<bool>(),
            e.at("strictlyRequired").get<bool>(),
            e.at("weight").get<double>()
        );
        v.push_back(x);
    }
}

inline void to_json(json& j, const std::vector<MagneticFilterOperation>& v) {
    j = json::array();
    for (auto x : v) {
        json e;
        e["filter"] = x.get_filter();
        e["invert"] = x.get_invert();
        e["log"] = x.get_log();
        e["strictlyRequired"] = x.get_strictly_required();
        e["weight"] = x.get_weight();
        j.push_back(e);
    }
}

double resolve_dimensional_values(MAS::Dimension dimensionValue, DimensionalValues preferredValue = DimensionalValues::NOMINAL);
inline double resolve_dimensional_values(MAS::Dimension dimensionValue, DimensionalValues preferredValue) {
    double doubleValue = 0;
    if (std::holds_alternative<MAS::DimensionWithTolerance>(dimensionValue)) {
        switch (preferredValue) {
            case DimensionalValues::MAXIMUM:
                if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().value();
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_nominal().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_nominal().value();
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().value();
                break;
            case DimensionalValues::NOMINAL:
                if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_nominal().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_nominal().value();
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().has_value() &&
                         std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().has_value())
                    doubleValue =
                        (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().value() +
                         std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().value()) /
                        2;
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().value();
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().value();
                break;
            case DimensionalValues::MINIMUM:
                if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_minimum().value();
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_nominal().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_nominal().value();
                else if (std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().has_value())
                    doubleValue = std::get<MAS::DimensionWithTolerance>(dimensionValue).get_maximum().value();
                break;
            default:
                throw std::runtime_error("Unknown type of dimension, options are {MAXIMUM, NOMINAL, MINIMUM}");
        }
    }
    else if (std::holds_alternative<double>(dimensionValue)) {
        doubleValue = std::get<double>(dimensionValue);
    }
    else {
        throw std::runtime_error("Unknown variant in dimensionValue, holding variant");
    }
    return doubleValue;
}

}