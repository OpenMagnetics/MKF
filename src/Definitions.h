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

}