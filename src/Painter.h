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
        int _mirroringDimension = Defaults().magneticFieldMirroringDimension;
        std::filesystem::path _filepath;
        std::map<std::string, std::string> _postProcessingChanges;
        std::map<std::string, std::string> _postProcessingColors;
        uint32_t _currentMapIndex = 0;

        // PainterModes _mode;
        // bool _logarithmicScale = false;
        // bool _includeFringing = true;



    public:

        Painter(std::filesystem::path filepath){
            _filepath = filepath;
            matplot::gcf()->quiet_mode(true);
            matplot::cla();
            matplot::hold(matplot::off);
            matplot::xticks({});
            matplot::yticks({});
        };
        virtual ~Painter() = default;
    
    ComplexField calculate_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex = 1);
    void paint_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex = 1, std::optional<ComplexField> inputField = std::nullopt);

    void export_svg();
    void export_png();

    void paint_core(MagneticWrapper magnetic);
    void paint_toroidal_core(CoreWrapper core);
    void paint_two_piece_set_core(CoreWrapper core);

    void paint_bobbin(MagneticWrapper magnetic);
    void paint_two_piece_set_bobbin(MagneticWrapper magnetic);

    void paint_coil_sections(MagneticWrapper magnetic);
    void paint_two_piece_set_winding_sections(MagneticWrapper magnetic);
    void paint_toroidal_winding_sections(MagneticWrapper magnetic);

    void paint_coil_layers(MagneticWrapper magnetic);
    void paint_two_piece_set_winding_layers(MagneticWrapper magnetic);
    void paint_toroidal_winding_layers(MagneticWrapper magnetic);

    void paint_coil_turns(MagneticWrapper magnetic);
    void paint_two_piece_set_winding_turns(MagneticWrapper magnetic);
    void paint_toroidal_winding_turns(MagneticWrapper magnetic);

    void paint_two_piece_set_margin(MagneticWrapper magnetic);
    void paint_toroidal_margin(MagneticWrapper magnetic);


};
} // namespace OpenMagnetics
