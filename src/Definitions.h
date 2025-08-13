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

void from_json(const json & j, DimensionalValues & x);
void to_json(json & j, const DimensionalValues & x);

inline void from_json(const json & j, DimensionalValues & x) {
    if (j == "maximum") x = DimensionalValues::MAXIMUM;
    else if (j == "nominal") x = DimensionalValues::NOMINAL;
    else if (j == "minimum") x = DimensionalValues::MINIMUM;
    else { throw std::runtime_error("Input JSON does not conform to schema!"); }
}

inline void to_json(json & j, const DimensionalValues & x) {
    switch (x) {
        case DimensionalValues::MAXIMUM: j = "maximum"; break;
        case DimensionalValues::NOMINAL: j = "nominal"; break;
        case DimensionalValues::MINIMUM: j = "minimum"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"DimensionalValues\": " + std::to_string(static_cast<int>(x)));
    }
}

enum class GappingType : int {
    GROUND,
    SPACER,
    RESIDUAL,
    DISTRIBUTED
};

void from_json(const json & j, GappingType & x);
void to_json(json & j, const GappingType & x);

inline void from_json(const json & j, GappingType & x) {
    if (j == "ground") x = GappingType::GROUND;
    else if (j == "spacer") x = GappingType::SPACER;
    else if (j == "residual") x = GappingType::RESIDUAL;
    else if (j == "distributed") x = GappingType::DISTRIBUTED;
    else { throw std::runtime_error("Input JSON does not conform to schema!"); }
}

inline void to_json(json & j, const GappingType & x) {
    switch (x) {
        case GappingType::GROUND: j = "ground"; break;
        case GappingType::SPACER: j = "spacer"; break;
        case GappingType::RESIDUAL: j = "residual"; break;
        case GappingType::DISTRIBUTED: j = "distributed"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"GappingType\": " + std::to_string(static_cast<int>(x)));
    }
}

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

void from_json(const json & j, OrderedIsolationSide & x);
void to_json(json & j, const OrderedIsolationSide & x);

inline void from_json(const json & j, OrderedIsolationSide & x) {
    if (j == "primary") x = OrderedIsolationSide::PRIMARY;
    else if (j == "secondary") x = OrderedIsolationSide::SECONDARY;
    else if (j == "tertiary") x = OrderedIsolationSide::TERTIARY;
    if (j == "quaternary") x = OrderedIsolationSide::QUATERNARY;
    else if (j == "quinary") x = OrderedIsolationSide::QUINARY;
    else if (j == "senary") x = OrderedIsolationSide::SENARY;
    if (j == "septenary") x = OrderedIsolationSide::SEPTENARY;
    else if (j == "octonary") x = OrderedIsolationSide::OCTONARY;
    else if (j == "nonary") x = OrderedIsolationSide::NONARY;
    if (j == "denary") x = OrderedIsolationSide::DENARY;
    else if (j == "undenary") x = OrderedIsolationSide::UNDENARY;
    else if (j == "duodenary") x = OrderedIsolationSide::DUODENARY;
    else { throw std::runtime_error("Input JSON does not conform to schema!"); }
}

inline void to_json(json & j, const OrderedIsolationSide & x) {
    switch (x) {
        case OrderedIsolationSide::PRIMARY: j = "primary"; break;
        case OrderedIsolationSide::SECONDARY: j = "secondary"; break;
        case OrderedIsolationSide::TERTIARY: j = "tertiary"; break;
        case OrderedIsolationSide::QUATERNARY: j = "quaternary"; break;
        case OrderedIsolationSide::QUINARY: j = "quinary"; break;
        case OrderedIsolationSide::SENARY: j = "senary"; break;
        case OrderedIsolationSide::SEPTENARY: j = "septenary"; break;
        case OrderedIsolationSide::OCTONARY: j = "octonary"; break;
        case OrderedIsolationSide::NONARY: j = "nonary"; break;
        case OrderedIsolationSide::DENARY: j = "denary"; break;
        case OrderedIsolationSide::UNDENARY: j = "undenary"; break;
        case OrderedIsolationSide::DUODENARY: j = "duodenary"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"OrderedIsolationSide\": " + std::to_string(static_cast<int>(x)));
    }
}


enum class PainterModes : int {
    CONTOUR,
    QUIVER,
    SCATTER,
};

void from_json(const json & j, PainterModes & x);
void to_json(json & j, const PainterModes & x);

inline void from_json(const json & j, PainterModes & x) {
    if (j == "contour") x = PainterModes::CONTOUR;
    else if (j == "quiver") x = PainterModes::QUIVER;
    else if (j == "scatter") x = PainterModes::SCATTER;
    else { throw std::runtime_error("Input JSON does not conform to schema!"); }
}

inline void to_json(json & j, const PainterModes & x) {
    switch (x) {
        case PainterModes::CONTOUR: j = "contour"; break;
        case PainterModes::QUIVER: j = "quiver"; break;
        case PainterModes::SCATTER: j = "scatter"; break;
        default: throw std::runtime_error("Unexpected value in enumeration \"PainterModes\": " + std::to_string(static_cast<int>(x)));
    }
}


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
    CORE_DC_AND_SKIN_LOSSES,
    LOSSES,
    LOSSES_NO_PROXIMITY,
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
    SKIN_LOSSES_DENSITY,
    VOLUME,
    AREA,
    HEIGHT,
    TEMPERATURE_RISE,
    VOLUME_TIMES_TEMPERATURE_RISE,
    LOSSES_TIMES_VOLUME,
    LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE,
    LOSSES_NO_PROXIMITY_TIMES_VOLUME,
    LOSSES_NO_PROXIMITY_TIMES_VOLUME_TIMES_TEMPERATURE_RISE,
    MAGNETOMOTIVE_FORCE
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
    else if (j == "Core DC And Skin Losses") x = MagneticFilters::CORE_DC_AND_SKIN_LOSSES;
    else if (j == "Losses") x = MagneticFilters::LOSSES;
    else if (j == "Losses No Proximity") x = MagneticFilters::LOSSES_NO_PROXIMITY;
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
    else if (j == "Skin Losses Density") x = MagneticFilters::SKIN_LOSSES_DENSITY;
    else if (j == "Volume") x = MagneticFilters::VOLUME;
    else if (j == "Area") x = MagneticFilters::AREA;
    else if (j == "Height") x = MagneticFilters::HEIGHT;
    else if (j == "Temperature Rise") x = MagneticFilters::TEMPERATURE_RISE;
    else if (j == "Losses Times Volume") x = MagneticFilters::LOSSES_TIMES_VOLUME;
    else if (j == "Volume Times Temperature Rise") x = MagneticFilters::VOLUME_TIMES_TEMPERATURE_RISE;
    else if (j == "Losses Times Volume Times Temperature Rise") x = MagneticFilters::LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE;
    else if (j == "Losses No Proximity Times Volume") x = MagneticFilters::LOSSES_NO_PROXIMITY_TIMES_VOLUME;
    else if (j == "Losses No Proximity Times Volume Times Temperature Rise") x = MagneticFilters::LOSSES_NO_PROXIMITY_TIMES_VOLUME_TIMES_TEMPERATURE_RISE;
    else if (j == "MagnetomotiveForce") x = MagneticFilters::MAGNETOMOTIVE_FORCE;
    else { throw std::runtime_error("Input JSON does not conform to MagneticFilters schema!"); }
}

inline void to_json(json & j, const MagneticFilters & x) {
    switch (x) {
        case MagneticFilters::AREA_PRODUCT: j = "Area Product"; break;
        case MagneticFilters::ENERGY_STORED: j = "Energy Stored"; break;
        case MagneticFilters::COST: j = "Cost"; break;
        case MagneticFilters::ESTIMATED_COST: j = "Estimated Cost"; break;
        case MagneticFilters::CORE_AND_DC_LOSSES: j = "Core And DC Losses"; break;
        case MagneticFilters::CORE_DC_AND_SKIN_LOSSES: j = "Core DC And Skin Losses"; break;
        case MagneticFilters::LOSSES: j = "Losses"; break;
        case MagneticFilters::LOSSES_NO_PROXIMITY: j = "Losses No Proximity"; break;
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
        case MagneticFilters::SKIN_LOSSES_DENSITY: j = "Skin Losses Density"; break;
        case MagneticFilters::VOLUME: j = "Volume"; break;
        case MagneticFilters::AREA: j = "Area"; break;
        case MagneticFilters::HEIGHT: j = "Height"; break;
        case MagneticFilters::TEMPERATURE_RISE: j = "Temperature Rise"; break;
        case MagneticFilters::LOSSES_TIMES_VOLUME: j = "Losses Times Volume"; break;
        case MagneticFilters::VOLUME_TIMES_TEMPERATURE_RISE: j = "Volume Times Temperature Rise"; break;
        case MagneticFilters::LOSSES_TIMES_VOLUME_TIMES_TEMPERATURE_RISE: j = "Losses Times Volume Times Temperature Rise"; break;
        case MagneticFilters::LOSSES_NO_PROXIMITY_TIMES_VOLUME: j = "Losses No Proximity Times Volume"; break;
        case MagneticFilters::LOSSES_NO_PROXIMITY_TIMES_VOLUME_TIMES_TEMPERATURE_RISE: j = "Losses No Proximity Times Volume Times Temperature Rise"; break;
        case MagneticFilters::MAGNETOMOTIVE_FORCE: j = "MagnetomotiveForce"; break;
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