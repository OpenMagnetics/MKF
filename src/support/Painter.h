#pragma once

#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "support/Utils.h"
#include <MAS.hpp>
#include "svg.hpp"
#ifdef ENABLE_MATPLOTPP
#include <matplot/matplot.h>
#endif
#include "support/Exceptions.h"

using namespace MAS;

namespace OpenMagnetics {

class PainterInterface {
    private:
    protected:

        // Convert uint to hex string
        std::string uint_to_hex(const uint32_t color, std::string prefix = "0x") {
            std::stringstream stream;
            stream << prefix << std::setfill ('0') << std::setw(6) << std::hex << color;
            return stream.str();
        }

        std::string key_to_rgb_color(uint32_t i) {
            std::stringstream stream;
            stream << "rgb(" << std::setfill (' ') << std::setw(3) << ((i >> 16) & 0xFF) << ", " << std::setfill (' ') << std::setw(3) << ((i >> 8) & 0xFF) << ", " << std::setfill (' ') << std::setw(3) << ((i >> 0) & 0xFF) << ")";
            return stream.str();    
        }

        std::string key_to_rgb_color(std::string s) {
            if (s[0] != '0') {
                return s;
            }
            auto i = stoi(s, 0, 16);
            std::stringstream stream;
            stream << "rgb(" << std::setfill (' ') << std::setw(3) << ((i >> 16) & 0xFF) << ", " << std::setfill (' ') << std::setw(3) << ((i >> 8) & 0xFF) << ", " << std::setfill (' ') << std::setw(3) << ((i >> 0) & 0xFF) << ")";
            return stream.str();
        }

        // Convert hex color string (#RRGGBB or RRGGBB) to uint32_t
        uint32_t hex_to_uint(const std::string& hex) {
            std::string hexStr = hex;
            
            // Remove leading '#' if present
            if (!hexStr.empty() && hexStr[0] == '#') {
                hexStr = hexStr.substr(1);
            }
            
            // Handle short format (3 digits) - expand to 6 digits
            if (hexStr.length() == 3) {
                hexStr = std::string() + hexStr[0] + hexStr[0] + 
                         hexStr[1] + hexStr[1] + 
                         hexStr[2] + hexStr[2];
            }
            
            // Ensure we have a valid hex string
            if (hexStr.length() < 6) {
                return 0x000000; // Return black as default
            }
            
            // Convert hex string to RGB
            uint32_t colorValue = std::stoul(hexStr, nullptr, 16);

            return colorValue;
        }

        // Clamp value between min and max
        double clamp(double value, double minVal, double maxVal) {
            if (value < minVal) return minVal;
            if (value > maxVal) return maxVal;
            return value;
        }

        uint32_t get_uint_color_from_ratio(double ratio) {
            // We want 4 regions (red->yellow, yellow->green, green->cyan, cyan->blue)
            int normalized = int(ratio * 256 * 4);
            if (normalized > 1023) normalized = 1023;
            
            int region = normalized / 256;
            int x = normalized % 256;
            
            uint8_t r = 0, g = 0, b = 0;
            
            switch (region)
            {
            case 0: // Red to Yellow: Red stays 255, Green increases from 0 to 255
                r = 255;
                g = x;
                b = 0;
                break;
            case 1: // Yellow to Green: Green stays 255, Red decreases from 255 to 0
                r = 255 - x;
                g = 255;
                b = 0;
                break;
            case 2: // Green to Cyan: Green stays 255, Blue increases from 0 to 255
                r = 0;
                g = 255;
                b = x;
                break;
            case 3: // Cyan to Blue: Blue stays 255, Green decreases from 255 to 0
                r = 0;
                g = 255 - x;
                b = 255;
                break;
            }
            
            return r + (g << 8) + (b << 16);
        }

        struct CoatingInfo {
            double strokeWidth;
            size_t numberLines;
            double lineRadiusIncrease;
            std::string coatingColor;
        };

        ComplexField calculate_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1);
        Field calculate_electric_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1);

        CoatingInfo process_coating(double insulationThickness, InsulationWireCoating coating) {
            InsulationWireCoatingType insulationWireCoatingType = coating.get_type().value();
            size_t numberLines = 0;
            double strokeWidth = 0;
            double lineRadiusIncrease = 0;
            auto coatingColor = settings.get_painter_color_insulation();

            switch (insulationWireCoatingType) {
                case InsulationWireCoatingType::BARE:
                    break;
                case InsulationWireCoatingType::ENAMELLED:
                    if (!coating.get_grade()) {
                        throw std::runtime_error("Enamelled wire missing grade");
                    }
                    numberLines = coating.get_grade().value() + 1;
                    lineRadiusIncrease = insulationThickness / coating.get_grade().value() * 2;
                    coatingColor = settings.get_painter_color_enamel();
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
                            coatingColor = settings.get_painter_color_silk();
                        }
                        else {
                            coatingColor = settings.get_painter_color_fep();
                            if (coating.get_material()) {
                                std::string coatingMaterial = Wire::resolve_coating_insulation_material(coating).get_name();
                                if (coatingMaterial == "PFA") {
                                    coatingColor = settings.get_painter_color_pfa();
                                }
                                else if (coatingMaterial == "FEP") {
                                    coatingColor = settings.get_painter_color_fep();
                                }
                                else if (coatingMaterial == "ETFE") {
                                    coatingColor = settings.get_painter_color_etfe();
                                }
                                else if (coatingMaterial == "TCA") {
                                    coatingColor = settings.get_painter_color_tca();
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
        virtual void paint_electric_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1, std::optional<Field> inputField = std::nullopt) = 0;
        virtual void paint_wire_losses(Magnetic magnetic, std::optional<Outputs> outputs = std::nullopt, std::optional<OperatingPoint> operatingPoint = std::nullopt, double temperature=defaults.ambientTemperature) = 0;
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
        virtual void paint_rectangle(double xCoordinate, double yCoordinate, double xDimension, double yDimension) = 0;
        virtual void paint_circle(double xCoordinate, double yCoordinate, double radius) = 0;
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
        bool _fieldPainted = false;

        std::vector<double> get_image_size(Magnetic magnetic);
        void set_opacity(double opacity) {
            _opacity = opacity;
        }
        void paint_two_piece_set_core(Core core);
        void paint_toroidal_core(Core core);
        void paint_two_piece_set_bobbin(Magnetic magnetic);
        void paint_two_piece_set_coil_sections(Magnetic magnetic);
        void paint_toroidal_coil_sections(Magnetic magnetic);
        void paint_two_piece_set_coil_layers(Magnetic magnetic);
        void paint_toroidal_coil_layers(Magnetic magnetic);
        void paint_two_piece_set_coil_turns(Magnetic magnetic);
        void paint_toroidal_coil_turns(Magnetic magnetic);
        void paint_two_piece_set_margin(Magnetic magnetic);
        void paint_toroidal_margin(Magnetic magnetic);

        void set_image_size(Wire wire);
        void set_image_size(Magnetic magnetic);
        void paint_round_wire(double xCoordinate, double yCoordinate, Wire wire, std::optional<std::string> label = std::nullopt);
        void paint_litz_wire(double xCoordinate, double yCoordinate, Wire wire, std::optional<std::string> label = std::nullopt);
        void paint_rectangular_wire(double xCoordinate, double yCoordinate, Wire wire, double angle, std::vector<double> center, std::optional<std::string> label = std::nullopt);
        void paint_rectangle(double xCoordinate, double yCoordinate, double xDimension, double yDimension, std::string cssClassName, SVG::Group* group = nullptr, double angle = 0, std::vector<double> center = {0, 0}, std::optional<std::string> label = std::nullopt);
        void paint_circle(double xCoordinate, double yCoordinate, double radius, std::string cssClassName, SVG::Group* group = nullptr, double fillAngle=360, double angle = 0, std::vector<double> center = {0, 0}, std::optional<std::string> label = std::nullopt);
        std::string get_color(double minimumValue, double maximumValue, std::string minimumColor, std::string maximumColor, double value);
        void paint_field_point(double xCoordinate, double yCoordinate, double xDimension, double yDimension, std::string color, std::string label);
        void paint_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex, std::optional<ComplexField> inputField);
        void paint_electric_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex, std::optional<Field> inputField);
        void paint_wire_losses(Magnetic magnetic, std::optional<Outputs> outputs = std::nullopt, std::optional<OperatingPoint> operatingPoint = std::nullopt, double temperature=defaults.ambientTemperature);

    public:
        SVG::SVG _root;
        BasicPainter(){};
        BasicPainter(std::filesystem::path filepath){
            _filename = filepath.filename();
            _filepath = filepath.remove_filename();
            _root = SVG::SVG();
                
            _root.style(".ferrite").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_ferrite()), std::regex("0x"), "#"));
            _root.style(".bobbin").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_bobbin()), std::regex("0x"), "#"));
            _root.style(".bobbin_translucent").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_bobbin()), std::regex("0x"), "#")).set_attr("opacity", 0.5);
            _root.style(".margin").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_margin()), std::regex("0x"), "#"));
            _root.style(".margin_translucent").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_margin()), std::regex("0x"), "#")).set_attr("opacity", 0.5);
            _root.style(".spacer").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_spacer()), std::regex("0x"), "#"));
            _root.style(".copper").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_copper()), std::regex("0x"), "#"));
            _root.style(".copper_translucent").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_copper()), std::regex("0x"), "#")).set_attr("opacity", 0.5);
            _root.style(".insulation").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_insulation()), std::regex("0x"), "#"));
            _root.style(".insulation_translucent").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_insulation()), std::regex("0x"), "#")).set_attr("opacity", 0.5);
            _root.style(".fr4").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_fr4()), std::regex("0x"), "#"));
            _root.style(".fr4_translucent").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_fr4()), std::regex("0x"), "#")).set_attr("opacity", 0.5);
            _root.style(".current_density").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_current_density()), std::regex("0x"), "#"));
            _root.style(".text").set_attr("fill", std::regex_replace(std::string(settings.get_painter_color_text()), std::regex("0x"), "#"));
            _root.style(".white").set_attr("fill", "#ffffff");
            _root.style(".point").set_attr("fill", "#ff0000");
        };
        virtual ~BasicPainter() = default;


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
    void paint_rectangle(double xCoordinate, double yCoordinate, double xDimension, double yDimension);
    void paint_circle(double xCoordinate, double yCoordinate, double radius);
    
    /**
     * @brief Paint waveforms from an OperatingPoint
     * 
     * Extracts voltage and current waveforms from all windings in the operating point
     * and renders them as SVG.
     * 
     * @param operatingPoint The operating point containing winding excitations
     * @param title Optional title for the plot
     * @param width SVG width in pixels
     * @param height SVG height in pixels
     * @return SVG string
     */
    std::string paint_operating_point_waveforms(
        const OperatingPoint& operatingPoint,
        const std::string& title = "Operating Point Waveforms",
        double width = 1200,
        double height = 800);
    
    /**
     * @brief Paint a single waveform as SVG
     * 
     * @param waveform The waveform to plot
     * @param name Name/label for the waveform
     * @param color Line color (hex format, e.g., "#FF0000")
     * @param xOffset X offset within the SVG
     * @param yOffset Y offset within the SVG
     * @param plotWidth Width of the plot area
     * @param plotHeight Height of the plot area
     */
    void paint_waveform_svg(
        const Waveform& waveform,
        const std::string& name,
        const std::string& color,
        double xOffset,
        double yOffset,
        double plotWidth,
        double plotHeight);
};

#ifdef ENABLE_MATPLOTPP
class AdvancedPainter : public PainterInterface {
    protected:
        double _fixedScale = constants.coilPainterScale;
        double _scale = constants.coilPainterScale;
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
        void paint_rectangle(double xCoordinate, double yCoordinate, double xDimension, double yDimension) {
            throw std::runtime_error("Not implemented");
        }
        void paint_circle(double xCoordinate, double yCoordinate, double radius) {
            throw std::runtime_error("Not implemented");
        }


    public:

        AdvancedPainter(){};
        AdvancedPainter(std::filesystem::path filepath, bool addProportionForColorBar = false, bool showTicks = false){
            _addProportionForColorBar = addProportionForColorBar;
            _filepath = filepath;
            matplot::gcf()->quiet_mode(true);
            matplot::cla();
            // Reset axis scales to linear (in case previous operations used loglog)
            auto ax = matplot::gca();
            ax->x_axis().scale(matplot::axis_type::axis_scale::linear);
            ax->y_axis().scale(matplot::axis_type::axis_scale::linear);
            if (!showTicks) {
                matplot::xticks({});
                matplot::yticks({});
            }
            matplot::hold(matplot::off);
        };
        virtual ~AdvancedPainter() = default;

    void paint_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1, std::optional<ComplexField> inputField = std::nullopt);
    void paint_electric_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1, std::optional<Field> inputField = std::nullopt) {
        throw std::runtime_error("Not implemented");
    }
    void paint_wire_losses(Magnetic magnetic, std::optional<Outputs> outputs = std::nullopt, std::optional<OperatingPoint> operatingPoint = std::nullopt, double temperature=defaults.ambientTemperature) {
        throw std::runtime_error("Not implemented");
    }

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
#endif  // ENABLE_MATPLOTPP

class Painter{
    public:
        std::shared_ptr<PainterInterface> _painter;

        static std::shared_ptr<PainterInterface> factory(bool useAdvancedPainter, std::filesystem::path filepath, bool addProportionForColorBar = false, bool showTicks = false);
        Painter(std::filesystem::path filepath, bool addProportionForColorBar = false, bool showTicks = false, bool useAdvancedPainter = false){
#ifdef ENABLE_MATPLOTPP
            if (addProportionForColorBar || showTicks) {
                _painter = factory(true, filepath, addProportionForColorBar, showTicks);
            }
            else {
                _painter = factory(useAdvancedPainter, filepath, addProportionForColorBar, showTicks);
            }
#else
            (void)useAdvancedPainter;
            (void)addProportionForColorBar;
            (void)showTicks;
            _painter = factory(false, filepath, false, false);
#endif
        };
        virtual ~Painter() = default;

    void paint_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1, std::optional<ComplexField> inputField = std::nullopt);
    void paint_electric_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex = 1, std::optional<Field> inputField = std::nullopt);
    void paint_wire_losses(Magnetic magnetic, std::optional<Outputs> outputs = std::nullopt, std::optional<OperatingPoint> operatingPoint = std::nullopt, double temperature=defaults.ambientTemperature);
    static double get_pixel_proportion_between_turns(std::vector<double> firstTurnCoordinates, std::vector<double> firstTurnDimensions, TurnCrossSectionalShape firstTurncrossSectionalShape, std::vector<double> secondTurnCoordinates, std::vector<double> secondTurnDimensions, TurnCrossSectionalShape secondTurncrossSectionalShape, std::vector<double> pixelCoordinates, double dimension);
    static double get_pixel_area_between_turns(std::vector<double> firstTurnCoordinates, std::vector<double> firstTurnDimensions, TurnCrossSectionalShape firstTurncrossSectionalShape, std::vector<double> secondTurnCoordinates, std::vector<double> secondTurnDimensions, TurnCrossSectionalShape secondTurncrossSectionalShape, std::vector<double> pixelCoordinates, double dimension);
    static std::pair<double, double> get_pixel_dimensions(Magnetic magnetic);

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
    void paint_rectangle(double xCoordinate, double yCoordinate, double xDimension, double yDimension);
    void paint_circle(double xCoordinate, double yCoordinate, double radius);

    void paint_curve(Curve2D curve2D, bool logScale = false);

};

} // namespace OpenMagnetics