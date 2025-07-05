#pragma once

#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include <MAS.hpp>
#include "svg.hpp"
#include <matplot/matplot.h>

using namespace MAS;

namespace OpenMagnetics {

class PainterInterface {
    private:
    protected:

        struct CoatingInfo {
            double strokeWidth;
            size_t numberLines;
            double lineRadiusIncrease;
            std::string coatingColor;
        };

        CoatingInfo process_coating(double insulationThickness, InsulationWireCoating coating) {
            InsulationWireCoatingType insulationWireCoatingType = coating.get_type().value();
            size_t numberLines = 0;
            double strokeWidth = 0;
            double lineRadiusIncrease = 0;
            auto coatingColor = settings->get_painter_color_insulation();

            switch (insulationWireCoatingType) {
                case InsulationWireCoatingType::BARE:
                    break;
                case InsulationWireCoatingType::ENAMELLED:
                    if (!coating.get_grade()) {
                        throw std::runtime_error("Enamelled wire missing grade");
                    }
                    numberLines = coating.get_grade().value() + 1;
                    lineRadiusIncrease = insulationThickness / coating.get_grade().value() * 2;
                    coatingColor = settings->get_painter_color_enamel();
                    break;
                case InsulationWireCoatingType::SERVED:
                case InsulationWireCoatingType::INSULATED:
                    {
                        if (!coating.get_number_layers()) {
                            throw std::runtime_error("Insulated wire missing number layers");
                        }
                        numberLines = coating.get_number_layers().value() + 1;
                        lineRadiusIncrease = insulationThickness / coating.get_number_layers().value() * 2;
                        if (insulationWireCoatingType == InsulationWireCoatingType::SERVED) {
                            coatingColor = settings->get_painter_color_silk();
                        }
                        else {
                            coatingColor = settings->get_painter_color_fep();
                            if (coating.get_material()) {
                                std::string coatingMaterial = std::get<std::string>(coating.get_material().value());
                                if (coatingMaterial == "PFA") {
                                    coatingColor = settings->get_painter_color_pfa();
                                }
                                else if (coatingMaterial == "FEP") {
                                    coatingColor = settings->get_painter_color_fep();
                                }
                                else if (coatingMaterial == "ETFE") {
                                    coatingColor = settings->get_painter_color_etfe();
                                }
                                else if (coatingMaterial == "TCA") {
                                    coatingColor = settings->get_painter_color_tca();
                                }
                                else {
                                    throw std::runtime_error("Unknown insulated wire material");
                                }

                            }
                        }
                    }
                    break;
                default:
                    throw std::runtime_error("Coating type plot not implemented yet");
            }
            strokeWidth = insulationThickness / 10 / numberLines;
            CoatingInfo coatingInfo;
            coatingInfo.strokeWidth = strokeWidth;
            coatingInfo.numberLines = numberLines;
            coatingInfo.lineRadiusIncrease = lineRadiusIncrease;
            coatingInfo.coatingColor = coatingColor;

            return coatingInfo;
        }

    public:
        virtual void paint_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1, std::optional<ComplexField> inputField = std::nullopt) = 0;
        virtual std::string export_svg() = 0;
        virtual void export_png() = 0;
        virtual void paint_core(Magnetic magnetic) = 0;
        virtual void paint_bobbin(Magnetic magnetic) = 0;
        virtual void paint_coil_sections(Magnetic magnetic) = 0;
        virtual void paint_coil_layers(Magnetic magnetic) = 0;
        virtual void paint_wire(Wire wire) = 0;
        virtual void paint_coil_turns(Magnetic magnetic) = 0;
        virtual void paint_wire_with_current_density(Wire wire, OperatingPoint operatingPoint, size_t windingIndex = 0) = 0;
        virtual void paint_wire_with_current_density(Wire wire, SignalDescriptor current, double frequency, double temperature=defaults.ambientTemperature) = 0;
        virtual void paint_waveform(Waveform waveform) = 0;
        virtual void paint_waveform(std::vector<double> waveform, std::optional<std::vector<double>> time = std::nullopt) = 0;
        virtual void paint_curve(Curve2D curve2D, bool logScale = false) = 0;
};

class BasicPainter : public PainterInterface {
    private:
        std::filesystem::path _filepath = ".";
        std::filesystem::path _filename = "_.svg";
        double _extraDimension = 1;
        double _opacity = 1;
        double _imageHeight = 0;
        double _imageWidth = 0;
        double _scale = constants.coilPainterScale;

        std::vector<double> get_image_size(Magnetic magnetic);
        void set_opacity(double opacity) {
            _opacity = opacity;
        }
        void paint_two_piece_set_core(Core core);
        void paint_two_piece_set_bobbin(Magnetic magnetic);
        void paint_two_piece_set_coil_sections(Magnetic magnetic);
        void paint_two_piece_set_coil_layers(Magnetic magnetic);
        void paint_two_piece_set_coil_turns(Magnetic magnetic);
        void paint_two_piece_set_margin(Magnetic magnetic);

        void set_image_size(Wire wire);
        void paint_round_wire(double xCoordinate, double yCoordinate, Wire wire);
        void paint_litz_wire(double xCoordinate, double yCoordinate, Wire wire);
        void paint_rectangle(double xCoordinate, double yCoordinate, double xDimension, double yDimension, std::string cssClassName);
        void paint_rectangular_wire(double xCoordinate, double yCoordinate, Wire wire, double angle, std::vector<double> center);

    public:
        SVG::SVG* _root;
        BasicPainter(){};
        BasicPainter(std::filesystem::path filepath){
            _filename = filepath.filename();
            _filepath = filepath.remove_filename();
            _root = new SVG::SVG();
                
            _root->style(".ferrite").set_attr("fill", std::regex_replace(std::string(settings->get_painter_color_ferrite()), std::regex("0x"), "#"));
            _root->style(".bobbin").set_attr("fill", std::regex_replace(std::string(settings->get_painter_color_bobbin()), std::regex("0x"), "#"));
            _root->style(".margin").set_attr("fill", std::regex_replace(std::string(settings->get_painter_color_margin()), std::regex("0x"), "#"));
            _root->style(".spacer").set_attr("fill", std::regex_replace(std::string(settings->get_painter_color_spacer()), std::regex("0x"), "#"));
            _root->style(".copper").set_attr("fill", std::regex_replace(std::string(settings->get_painter_color_copper()), std::regex("0x"), "#"));
            _root->style(".insulation").set_attr("fill", std::regex_replace(std::string(settings->get_painter_color_insulation()), std::regex("0x"), "#"));
            _root->style(".current_density").set_attr("fill", std::regex_replace(std::string(settings->get_painter_color_current_density()), std::regex("0x"), "#"));
            _root->style(".text").set_attr("fill", std::regex_replace(std::string(settings->get_painter_color_text()), std::regex("0x"), "#"));
            _root->style(".white").set_attr("fill", "#ffffff");
        };
        virtual ~BasicPainter() = default;

    void paint_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1, std::optional<ComplexField> inputField = std::nullopt) {
        throw std::runtime_error("Not implemented in basic painter");
    }
    std::string export_svg();
    void export_png() {
        throw std::runtime_error("Not implemented in basic painter");
    }
    void paint_core(Magnetic magnetic);
    void paint_bobbin(Magnetic magnetic);
    void paint_coil_sections(Magnetic magnetic);
    void paint_coil_layers(Magnetic magnetic);
    void paint_wire(Wire wire);
    void paint_coil_turns(Magnetic magnetic);
    void paint_wire_with_current_density(Wire wire, OperatingPoint operatingPoint, size_t windingIndex = 0) {
        throw std::runtime_error("Not implemented in basic painter");
    }
    void paint_wire_with_current_density(Wire wire, SignalDescriptor current, double frequency, double temperature=defaults.ambientTemperature) {
        throw std::runtime_error("Not implemented in basic painter");
    }
    void paint_waveform(Waveform waveform) {
        throw std::runtime_error("Not implemented in basic painter");
    }
    void paint_waveform(std::vector<double> waveform, std::optional<std::vector<double>> time = std::nullopt) {
        throw std::runtime_error("Not implemented in basic painter");
    }
    void paint_curve(Curve2D curve2D, bool logScale = false) {
        throw std::runtime_error("Not implemented in basic painter");
    }
};

class AdvancedPainter : public PainterInterface {
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

        ComplexField calculate_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1);
        ComplexField calculate_magnetic_field_additional_coordinates(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex);
        std::vector<double> get_image_size(Magnetic magnetic);
        void set_image_size(Wire wire);
        std::string paint_rectangle(std::vector<double> cornerData, bool fill=true, double strokeWidth=0);

        void increment_current_map_index();
        void paint_toroidal_core(Magnetic magnetic);
        void paint_two_piece_set_core(Magnetic magnetic);
        void paint_two_piece_set_bobbin(Magnetic magnetic);
        void paint_two_piece_set_winding_sections(Magnetic magnetic);
        void paint_toroidal_winding_sections(Magnetic magnetic);
        void paint_two_piece_set_winding_layers(Magnetic magnetic);
        void paint_toroidal_winding_layers(Magnetic magnetic);
        void paint_round_wire(double xCoordinate, double yCoordinate, Wire wire);
        void paint_litz_wire(double xCoordinate, double yCoordinate, Wire wire);
        void paint_rectangular_wire(double xCoordinate, double yCoordinate, Wire wire, double angle=0, std::vector<double> center=std::vector<double>{});
        void paint_two_piece_set_winding_turns(Magnetic magnetic);
        void paint_toroidal_winding_turns(Magnetic magnetic);
        void paint_two_piece_set_margin(Magnetic magnetic);
        void paint_toroidal_margin(Magnetic magnetic);
        void calculate_extra_margin_for_toroidal_cores(Magnetic magnetic);
        void paint_background(Magnetic magnetic);


    public:

        AdvancedPainter(){};
        AdvancedPainter(std::filesystem::path filepath, bool addProportionForColorBar = false, bool showTicks = false){
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
        virtual ~AdvancedPainter() = default;

    void paint_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1, std::optional<ComplexField> inputField = std::nullopt);

    std::string export_svg();
    void export_png();

    void paint_core(Magnetic magnetic);

    void paint_bobbin(Magnetic magnetic);

    void paint_coil_sections(Magnetic magnetic);

    void paint_coil_layers(Magnetic magnetic);

    void paint_wire(Wire wire);
    void paint_coil_turns(Magnetic magnetic);
    void paint_wire_with_current_density(Wire wire, OperatingPoint operatingPoint, size_t windingIndex = 0);
    void paint_wire_with_current_density(Wire wire, SignalDescriptor current, double frequency, double temperature=defaults.ambientTemperature);


    void paint_waveform(Waveform waveform);
    void paint_waveform(std::vector<double> waveform, std::optional<std::vector<double>> time = std::nullopt);

    void paint_curve(Curve2D curve2D, bool logScale = false);
};

class Painter{
    public:
        std::shared_ptr<PainterInterface> _painter;

        static std::shared_ptr<PainterInterface> factory(bool useAdvancedPainter, std::filesystem::path filepath, bool addProportionForColorBar = false, bool showTicks = false);
        Painter(std::filesystem::path filepath, bool addProportionForColorBar = false, bool showTicks = false, bool useAdvancedPainter = true){
            _painter = factory(useAdvancedPainter, filepath, addProportionForColorBar, showTicks);
        };
        virtual ~Painter() = default;

    void paint_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1, std::optional<ComplexField> inputField = std::nullopt);

    std::string export_svg();
    void export_png();

    void paint_core(Magnetic magnetic);

    void paint_bobbin(Magnetic magnetic);

    void paint_coil_sections(Magnetic magnetic);

    void paint_coil_layers(Magnetic magnetic);

    void paint_wire(Wire wire);
    void paint_coil_turns(Magnetic magnetic);
    void paint_wire_with_current_density(Wire wire, OperatingPoint operatingPoint, size_t windingIndex = 0);
    void paint_wire_with_current_density(Wire wire, SignalDescriptor current, double frequency, double temperature=defaults.ambientTemperature);

    void paint_waveform(Waveform waveform);
    void paint_waveform(std::vector<double> data, std::optional<std::vector<double>> time = std::nullopt);

    void paint_curve(Curve2D curve2D, bool logScale = false);

};

} // namespace OpenMagnetics
