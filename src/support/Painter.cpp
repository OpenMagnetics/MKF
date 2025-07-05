#include "physical_models/MagneticField.h"
#include "support/Painter.h"
#include "support/CoilMesher.h"
#include "MAS.hpp"
#include "support/Utils.h"
#include "json.hpp"
#include <matplot/matplot.h>
#include <cfloat>
#include <chrono>
#include <thread>
#include <filesystem>

namespace OpenMagnetics {

std::shared_ptr<PainterInterface> Painter::factory(bool useAdvancedPainter, std::filesystem::path filepath, bool addProportionForColorBar, bool showTicks) {
    if (useAdvancedPainter) {
        return std::make_shared<AdvancedPainter>(filepath, addProportionForColorBar, showTicks);
    }
    else {
        return std::make_shared<BasicPainter>(filepath);
    }
}

void Painter::paint_magnetic_field(OperatingPoint operatingPoint, Magnetic magnetic, size_t harmonicIndex, std::optional<ComplexField> inputField) {
    _painter->paint_magnetic_field(operatingPoint, magnetic, harmonicIndex, inputField);
}

std::string Painter::export_svg() {
    return _painter->export_svg();
}

void Painter::export_png() {
    _painter->export_png();
}

void Painter::paint_core(Magnetic magnetic) {
    _painter->paint_core(magnetic);
}

void Painter::paint_bobbin(Magnetic magnetic) {
    _painter->paint_bobbin(magnetic);
}

void Painter::paint_coil_sections(Magnetic magnetic) {
    _painter->paint_coil_sections(magnetic);
}

void Painter::paint_coil_layers(Magnetic magnetic) {
    _painter->paint_coil_layers(magnetic);
}

void Painter::paint_coil_turns(Magnetic magnetic) {
    _painter->paint_coil_turns(magnetic);
}

void Painter::paint_wire(Wire wire) {
    _painter->paint_wire(wire);
}

void Painter::paint_wire_with_current_density(Wire wire, OperatingPoint operatingPoint, size_t windingIndex) {
    _painter->paint_wire_with_current_density(wire, operatingPoint, windingIndex);
}

void Painter::paint_wire_with_current_density(Wire wire, SignalDescriptor current, double frequency, double temperature) {
    _painter->paint_wire_with_current_density(wire, current, frequency, temperature);
}

void Painter::paint_waveform(Waveform waveform) {
    paint_waveform(waveform.get_data(), waveform.get_time());
}

void Painter::paint_waveform(std::vector<double> data, std::optional<std::vector<double>> time) {
    _painter->paint_waveform(data, time);
}

void Painter::paint_curve(Curve2D curve2D, bool logScale) {
    _painter->paint_curve(curve2D, logScale);
}

} // namespace OpenMagnetics
