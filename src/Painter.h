#pragma once

#include "Defaults.h"
#include "CoreWrapper.h"
#include "MagneticWrapper.h"
#include "Utils.h"
#include <MAS.hpp>
#include <matplot/matplot.h>

namespace OpenMagnetics {

class Painter{
    public:
        enum class PainterModes : int {
            CONTOUR,
            QUIVER,
        };

    protected:
        double _scale = 30000;
        size_t _numberPointsX = 25;
        size_t _numberPointsY = 50;
        std::filesystem::path _filepath;
        PainterModes _mode;
        bool _logarithmicScale = false;
        bool _includeFringing = true;
        std::optional<double> _maximumValueColorbar = std::nullopt;
        std::optional<double> _minimumValueColorbar = std::nullopt;

        std::string _colorFerrite = "0x8F7b7c7d";
        std::string _colorBobbin = "0x8F1b1b1b";
        std::string _colorCopper = "0x8Fb87333";
        std::string _colorInsulation = "0x18E37E00";
 
        int _mirroringDimension = Defaults().magneticFieldMirroringDimension;



    public:

        Painter(std::filesystem::path filepath, PainterModes mode = PainterModes::QUIVER){
            _filepath = filepath;
            _mode = mode;
            matplot::gcf()->quiet_mode(true);
            matplot::cla();
            matplot::hold(matplot::off);
            matplot::xticks({});
            matplot::yticks({});
        };
        virtual ~Painter() = default;

    void set_number_points_x(double numberPoints) {
        _numberPointsX = numberPoints;
    }
    void set_number_points_y(double numberPoints) {
        _numberPointsY = numberPoints;
    }
    void set_mode(PainterModes mode) {
        _mode = mode;
    }
    void set_logarithmic_scale(bool mode) {
        _logarithmicScale = mode;
    }
    void set_maximum_scale_value(std::optional<double> value) {
        _maximumValueColorbar = value;
    }
    void set_minimum_scale_value(std::optional<double> value) {
        _minimumValueColorbar = value;
    }
    void set_fringing_effect(bool value) {
        _includeFringing = value;
    }
    void set_mirroring_dimension(int mirroringDimension) {
        _mirroringDimension = mirroringDimension;
    }
    
    ComplexField calculate_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex = 1);
    void paint_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex = 1, std::optional<ComplexField> inputField = std::nullopt);

    void export_svg();
    void export_png();

    void paint_core(const MagneticWrapper& magnetic);
    void paint_two_piece_set_core(CoreWrapper core);
    void paint_bobbin(const MagneticWrapper& magnetic);
    void paint_two_piece_set_bobbin(MagneticWrapper magnetic);
    void paint_coil_sections(const MagneticWrapper& magnetic);
    void paint_two_piece_set_winding_sections(const MagneticWrapper& magnetic);
    void paint_coil_layers(const MagneticWrapper& magnetic);
    void paint_two_piece_set_winding_layers(const MagneticWrapper& magnetic);
    void paint_coil_turns(const MagneticWrapper& magnetic);
    void paint_two_piece_set_winding_turns(const MagneticWrapper& magnetic);


};
} // namespace OpenMagnetics
