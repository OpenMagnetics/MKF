#include <MAS.hpp>
#include "constructive_models/Magnetic.h"
#include "constructive_models/Mas.h"
#include "processors/Outputs.h"
#include "physical_models/Reluctance.h"
#include "support/Utils.h"
#include "support/Exceptions.h"

namespace OpenMagnetics {

Bobbin Magnetic::get_bobbin() {
    return get_mutable_coil().resolve_bobbin();
}

std::vector<Wire> Magnetic::get_wires() {
    return get_mutable_coil().get_wires();
}

std::vector<double> Magnetic::get_turns_ratios() {
    return get_mutable_coil().get_turns_ratios();
}

std::vector<double> Magnetic::get_turns_ratios(MAS::Magnetic magnetic) {
    std::vector<double> turnsRatios;
    for (size_t windingIndex = 1; windingIndex < magnetic.get_coil().get_functional_description().size(); ++windingIndex) {
        turnsRatios.push_back(double(magnetic.get_coil().get_functional_description()[0].get_number_turns()) / magnetic.get_coil().get_functional_description()[windingIndex].get_number_turns());
    }
    return turnsRatios;
}

Wire Magnetic::get_wire(size_t windingIndex) {
    return get_mutable_coil().resolve_wire(windingIndex);
}

std::string Magnetic::get_reference() {
    if (get_manufacturer_info()) {
        if (get_manufacturer_info()->get_reference()) {
            return get_manufacturer_info()->get_reference().value();
        }
    }
    return "Custom component made with OpenMagnetic";
}

std::vector<double> Magnetic::get_maximum_dimensions() {
    if (!_maximumDimensions) {
        auto coreMaximumDimensions = get_mutable_core().get_maximum_dimensions();
        auto coilMaximumDimensions = get_mutable_coil().get_maximum_dimensions();
        _maximumDimensions = {std::max(coreMaximumDimensions[0], coilMaximumDimensions[0]),
                              std::max(coreMaximumDimensions[1], coilMaximumDimensions[1]),
                              std::max(coreMaximumDimensions[2], coilMaximumDimensions[2])};
    }
    return _maximumDimensions.value();
}

bool fits_one_dimension(std::vector<double> magneticDimensions, double dimension) {
    if (magneticDimensions[2] <= dimension || magneticDimensions[1] <= dimension || magneticDimensions[0] <= dimension) {
        return true;
    }
    return false;
}

bool fits_two_dimensions(std::vector<double> magneticDimensions, double firstDimension, double secondDimension) {
    if ((magneticDimensions[2] <= firstDimension && (magneticDimensions[1] <= secondDimension || magneticDimensions[0] <= secondDimension)) ||
        (magneticDimensions[1] <= firstDimension && (magneticDimensions[2] <= secondDimension || magneticDimensions[0] <= secondDimension)) ||
        (magneticDimensions[0] <= firstDimension && (magneticDimensions[1] <= secondDimension || magneticDimensions[2] <= secondDimension))) {
        return true;
    }
    return false;
}

bool fits_three_dimensions(std::vector<double> magneticDimensions, double firstDimension, double secondDimension, double thirdDimension) {
    if ((magneticDimensions[2] <= firstDimension && magneticDimensions[1] <= secondDimension && magneticDimensions[0] <= thirdDimension) ||
        (magneticDimensions[2] <= firstDimension && magneticDimensions[1] <= thirdDimension && magneticDimensions[0] <= secondDimension) ||
        (magneticDimensions[2] <= secondDimension && magneticDimensions[1] <= firstDimension && magneticDimensions[0] <= thirdDimension) ||
        (magneticDimensions[2] <= secondDimension && magneticDimensions[1] <= thirdDimension && magneticDimensions[0] <= firstDimension) ||
        (magneticDimensions[2] <= thirdDimension && magneticDimensions[1] <= firstDimension && magneticDimensions[0] <= secondDimension) ||
        (magneticDimensions[2] <= thirdDimension && magneticDimensions[1] <= secondDimension && magneticDimensions[0] <= firstDimension)){
        return true;
    }
    return false;
}

bool Magnetic::fits(MaximumDimensions maximumDimensions, bool allowRotation) {

    auto magneticDimensions = get_maximum_dimensions();

    if (!maximumDimensions.get_depth() && !maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        return true;
    }
    else if (maximumDimensions.get_depth() && !maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        auto depth = maximumDimensions.get_depth().value();
        if (allowRotation) {
            return fits_one_dimension(magneticDimensions, depth);
        }
        else {
            return magneticDimensions[2] <= depth;
        }
    }
    else if (!maximumDimensions.get_depth() && maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        auto height = maximumDimensions.get_height().value();
        if (allowRotation) {
            return fits_one_dimension(magneticDimensions, height);
        }
        else {
            return magneticDimensions[1] <= height;
        }
    }
    else if (!maximumDimensions.get_depth() && !maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto width = maximumDimensions.get_width().value();
        if (allowRotation) {
            return fits_one_dimension(magneticDimensions, width);
        }
        else {
            return magneticDimensions[0] <= width;
        }
    }
    else if (maximumDimensions.get_depth() && maximumDimensions.get_height() && !maximumDimensions.get_width()) {
        auto depth = maximumDimensions.get_depth().value();
        auto height = maximumDimensions.get_height().value();
        if (allowRotation) {
            return fits_two_dimensions(magneticDimensions, depth, height);
        }
        else {
            return magneticDimensions[2] <= depth && magneticDimensions[1] <= height;
        }
    }
    else if (!maximumDimensions.get_depth() && maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto width = maximumDimensions.get_width().value();
        auto height = maximumDimensions.get_height().value();
        if (allowRotation) {
            return fits_two_dimensions(magneticDimensions, width, height);
        }
        else {
            return magneticDimensions[0] <= width && magneticDimensions[1] <= height;
        }
    }
    else if (maximumDimensions.get_depth() && !maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto width = maximumDimensions.get_width().value();
        auto depth = maximumDimensions.get_depth().value();
        if (allowRotation) {
            return fits_two_dimensions(magneticDimensions, width, depth);
        }
        else {
            return magneticDimensions[2] <= depth && magneticDimensions[0] <= width;
        }
    }
    else if (maximumDimensions.get_depth() && maximumDimensions.get_height() && maximumDimensions.get_width()) {
        auto depth = maximumDimensions.get_depth().value();
        auto height = maximumDimensions.get_height().value();
        auto width = maximumDimensions.get_width().value();
        if (allowRotation) {
            return fits_three_dimensions(magneticDimensions, depth, height, width);
        }
        else {
            return magneticDimensions[2] <= depth && magneticDimensions[1] <= height && magneticDimensions[0] <= width;
        }
    }
    else {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Not sure how this happened");
    }
}

double Magnetic::calculate_saturation_current(double temperature) {
    auto magneticFluxDensitySaturation = get_mutable_core().get_magnetic_flux_density_saturation();
    auto numberTurns = get_mutable_coil().get_number_turns(0);
    auto effectiveArea = get_mutable_core().get_effective_area();
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory();
    auto initialPermeability = get_mutable_core().get_initial_permeability(temperature);
    auto reluctance = reluctanceModel->get_core_reluctance(core, initialPermeability).get_core_reluctance();

    double saturationCurrent = magneticFluxDensitySaturation * effectiveArea * reluctance / numberTurns;
    return saturationCurrent;
}

void to_file(std::filesystem::path filepath, const Magnetic & x) {
    OpenMagnetics::Mas mas;
    auto inputs = OpenMagnetics::get_defaults_inputs(x);
    std::vector<OpenMagnetics::Outputs> outputs = {OpenMagnetics::Outputs()};
    mas.set_magnetic(x);
    mas.set_inputs(inputs);
    mas.set_outputs(outputs);
    json masJson;
    to_json(masJson, mas);

    std::ofstream myfile;
    myfile.open(filepath);
    myfile << masJson;
}

} // namespace OpenMagnetics
