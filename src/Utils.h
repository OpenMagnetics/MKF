#pragma once
#include "WireWrapper.h"
#include "CoreWrapper.h"
#include "BobbinWrapper.h"
#include "InsulationMaterialWrapper.h"

#include "Constants.h"

#include <MAS.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>
#include <complex>
#include "spline.h"
#include <exception>

extern std::vector<OpenMagnetics::CoreWrapper> coreDatabase;
extern std::map<std::string, OpenMagnetics::CoreMaterial> coreMaterialDatabase;
extern std::map<std::string, OpenMagnetics::CoreShape> coreShapeDatabase;
extern std::map<std::string, OpenMagnetics::WireWrapper> wireDatabase;
extern std::map<std::string, OpenMagnetics::BobbinWrapper> bobbinDatabase;
extern std::map<std::string, OpenMagnetics::InsulationMaterialWrapper> insulationMaterialDatabase;
extern std::map<std::string, OpenMagnetics::WireMaterial> wireMaterialDatabase;
extern OpenMagnetics::Constants constants;

extern tk::spline bobbinFillingFactorInterpWidth;
extern tk::spline bobbinFillingFactorInterpHeight;
extern tk::spline bobbinWindingWindowInterpWidth;
extern tk::spline bobbinWindingWindowInterpHeight;
extern std::map<std::string, tk::spline> wireFillingFactorInterps;
extern std::map<std::string, tk::spline> wirePackingFactorInterps;
extern std::map<std::string, tk::spline> wireOuterDimensionInterps;
extern std::map<std::string, tk::spline> wireConductingAreaInterps;
extern double minBobbinWidth;
extern double maxBobbinWidth;
extern double minBobbinHeight;
extern double maxBobbinHeight;
extern double minWindingWindowWidth;
extern double maxWindingWindowWidth;
extern double minWindingWindowHeight;
extern double maxWindingWindowHeight;

extern std::map<std::string, tk::spline> complexPermeabilityRealInterps;
extern std::map<std::string, tk::spline> complexPermeabilityImaginaryInterps;

extern std::map<std::string, tk::spline> lossFactorInterps;

extern bool _addInternalData;

extern std::map<std::string, double> minWireConductingDimensions;
extern std::map<std::string, double> maxWireConductingDimensions;
extern std::map<std::string, int64_t> minLitzWireNumberConductors;
extern std::map<std::string, int64_t> maxLitzWireNumberConductors;


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

double resolve_dimensional_values(OpenMagnetics::Dimension dimensionValue, DimensionalValues preferredValue = DimensionalValues::NOMINAL);
bool check_requirement(DimensionWithTolerance requirement, double value);
OpenMagnetics::CoreWrapper find_core_by_name(std::string name);
OpenMagnetics::CoreMaterial find_core_material_by_name(std::string name);
OpenMagnetics::CoreShape find_core_shape_by_name(std::string name);
OpenMagnetics::WireWrapper find_wire_by_name(std::string name);
OpenMagnetics::WireWrapper find_wire_by_dimension(double dimension, std::optional<WireType> wireType=std::nullopt, std::optional<WireStandard> wireStandard=std::nullopt);
OpenMagnetics::BobbinWrapper find_bobbin_by_name(std::string name);
OpenMagnetics::InsulationMaterialWrapper find_insulation_material_by_name(std::string name);
OpenMagnetics::WireMaterial find_wire_material_by_name(std::string name);
OpenMagnetics::CoreShape find_core_shape_by_winding_window_perimeter(double desiredPerimeter);


void clear_loaded_cores();
void clear_databases();
void load_cores(bool includeToroidalCores=true, bool useOnlyCoresInStock=false, bool includeConcentricCores=true);
void load_core_materials();
void load_core_shapes(bool withAliases=true);
void load_wires();
void load_bobbins();
void load_insulation_materials();
void load_wire_materials();
void load_databases(json data, bool withAliases=true, bool addInternalData=true);

std::vector<std::string> get_shape_names(std::optional<std::string> manufacturer);
std::vector<std::string> get_material_names(std::optional<std::string> manufacturer);
std::vector<std::string> get_shape_names();
std::vector<std::string> get_wire_names();
std::vector<std::string> get_bobbin_names();
std::vector<std::string> get_insulation_material_names();
std::vector<std::string> get_wire_material_names();

std::vector<OpenMagnetics::CoreMaterial> get_materials(std::optional<std::string> manufacturer);
std::vector<OpenMagnetics::CoreShape> get_shapes(bool includeToroidal = true);
std::vector<OpenMagnetics::WireWrapper> get_wires(std::optional<WireType> wireType=std::nullopt, std::optional<WireStandard> wireStandard=std::nullopt);
std::vector<OpenMagnetics::BobbinWrapper> get_bobbins();
std::vector<OpenMagnetics::InsulationMaterialWrapper> get_insulation_materials();
std::vector<OpenMagnetics::WireMaterial> get_wire_materials();

template<int decimals> double roundFloat(double value);

bool is_size_power_of_2(std::vector<double> data);
size_t round_up_size_to_power_of_2(std::vector<double> data);
size_t round_up_size_to_power_of_2(size_t size);

std::vector<size_t> get_main_current_harmonic_indexes(OperatingPoint operatingPoint, double windingLossesHarmonicAmplitudeThreshold);
std::vector<size_t> get_main_harmonic_indexes(Harmonics harmonics, double windingLossesHarmonicAmplitudeThreshold);

double roundFloat(double value, int64_t decimals);
double ceilFloat(double value, size_t decimals);
double floorFloat(double value, size_t decimals);

CoreShape flatten_dimensions(CoreShape shape);
std::map<std::string, double> flatten_dimensions(std::map<std::string, Dimension> dimensions);

double try_get_duty_cycle(Waveform waveform, double frequency);
std::complex<double> modified_bessel_first_kind(double order, std::complex<double> x);
std::complex<double> bessel_first_kind(double order, std::complex<double> z);
double kelvin_function_real(double order, double x);
double kelvin_function_imaginary(double order, double x);
double derivative_kelvin_function_real(double order, double x);
double derivative_kelvin_function_imaginary(double order, double x);
double comp_ellint_1(double x);
double comp_ellint_2(double x);

std::string to_title_case(std::string text);

double wound_distance_to_angle(double distance, double radius);
double angle_to_wound_distance(double angle, double radius);


bool check_collisions(std::map<std::string, std::vector<double>> dimensionsByName, std::map<std::string, std::vector<double>> coordinatesByName, bool roundTurn = false);
IsolationSide get_isolation_side_from_index(size_t index);

std::vector<std::string> split(std::string s, std::string delimiter);
std::vector<double> linear_spaced_array(double a, double b, size_t N);
std::vector<double> logarithmic_spaced_array(double a, double b, size_t N);

double decibels_to_amplitude(double decibels);
double amplitude_to_decibels(double amplitude);

std::string fix_filename(std::string filename);

} // namespace OpenMagnetics