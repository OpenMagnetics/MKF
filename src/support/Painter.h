#pragma once

#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include <MAS.hpp>
#include <matplot/matplot.h>

using namespace MAS;

namespace OpenMagnetics {

class Painter{
    protected:
        double _fixedScale = 30000;
        double _scale = 30000;
        double _fontSize = 10;
        double _extraDimension = 1;
        double _offsetForColorBar = 0;
        bool _paintingFullMagnetic = false;
        bool _addProportionForColorBar = true;
        int _mirroringDimension = defaults.magneticFieldMirroringDimension;
        std::filesystem::path _filepath;
        std::map<std::string, std::string> _postProcessingChanges;
        std::map<std::string, std::string> _postProcessingColors;
        std::vector<std::string> _postProcessingDefs;
        uint32_t _currentMapIndex = 0x000001;

        // PainterModes _mode;
        // bool _logarithmicScale = false;
        // bool _includeFringing = true;


    public:

        Painter(std::filesystem::path filepath, bool addProportionForColorBar = false, bool showTicks = false){
            _addProportionForColorBar = addProportionForColorBar;
            _filepath = filepath;
            matplot::gcf()->quiet_mode(true);
            matplot::cla();
            if (!showTicks) {
                matplot::xticks({});
                matplot::yticks({});
            }
            matplot::hold(matplot::off);
        };
        virtual ~Painter() = default;

    void increment_current_map_index();
    // std::string get_current_post_processing_index(std::optional<std::string> changes = std::nullopt, std::optional<std::string> colors = std::nullopt);

    ComplexField calculate_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1);
    ComplexField calculate_magnetic_field_additional_coordinates(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex);
    void paint_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1, std::optional<ComplexField> inputField = std::nullopt);


    std::vector<double> get_image_size(Magnetic magnetic);
    void set_image_size(Wire wire);

    void export_svg();
    void export_png();
    std::string paint_rectangle(std::vector<double> cornerData, bool fill=true, double strokeWidth=0);

    void paint_core(Magnetic magnetic);
    void paint_toroidal_core(Magnetic magnetic);
    void paint_two_piece_set_core(Magnetic magnetic);

    void paint_bobbin(Magnetic magnetic);
    void paint_two_piece_set_bobbin(Magnetic magnetic);

    void paint_coil_sections(Magnetic magnetic);
    void paint_two_piece_set_winding_sections(Magnetic magnetic);
    void paint_toroidal_winding_sections(Magnetic magnetic);

    void paint_coil_layers(Magnetic magnetic);
    void paint_two_piece_set_winding_layers(Magnetic magnetic);
    void paint_toroidal_winding_layers(Magnetic magnetic);

    void paint_wire(Wire wire);
    void paint_round_wire(double xCoordinate, double yCoordinate, Wire wire);
    void paint_litz_wire(double xCoordinate, double yCoordinate, Wire wire);
    void paint_rectangular_wire(double xCoordinate, double yCoordinate, Wire wire, double angle=0, std::vector<double> center=std::vector<double>{});
    void paint_coil_turns(Magnetic magnetic);
    void paint_two_piece_set_winding_turns(Magnetic magnetic);
    void paint_toroidal_winding_turns(Magnetic magnetic);
    void paint_wire_with_current_density(Wire wire, OperatingPoint operatingPoint, size_t windingIndex = 0);
    void paint_wire_with_current_density(Wire wire, SignalDescriptor current, double frequency, double temperature=defaults.ambientTemperature);

    void paint_two_piece_set_margin(Magnetic magnetic);
    void paint_toroidal_margin(Magnetic magnetic);

    void calculate_extra_margin_for_toroidal_cores(Magnetic magnetic);

    void paint_waveform(Waveform waveform);
    void paint_waveform(std::vector<double> waveform, std::optional<std::vector<double>> time = std::nullopt);

    void paint_curve(Curve2D curve2D, bool logScale = false);

    void paint_background(Magnetic magnetic);
};
} // namespace OpenMagnetics
