#include "MagneticField.h"
#include "LeakageInductance.h"
#include "MAS.hpp"
#include "Utils.h"
#include "json.hpp"
#include <matplot/matplot.h>
#include <cfloat>
#include <../tests/TestingUtils.h>

namespace OpenMagnetics {

ComplexField LeakageInductance::calculate_magnetic_field(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex) {
    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()){
        throw std::runtime_error("Bobbin is not processed");
    }
    if (!bobbin.get_processed_description()->get_column_width()){
        throw std::runtime_error("Bobbin is missing column width");
    }
    double bobbinWidthStart = bobbin.get_processed_description()->get_winding_windows()[0].get_coordinates().value()[0] - bobbin.get_processed_description()->get_winding_windows()[0].get_width().value() / 2;
    double bobbinWidth = bobbin.get_processed_description()->get_winding_windows()[0].get_width().value();

    double coreColumnWidth = magnetic.get_mutable_core().get_columns()[0].get_width();
    double coreColumnHeight = magnetic.get_mutable_core().get_columns()[0].get_height();

    auto harmonics = operatingPoint.get_excitations_per_winding()[0].get_current()->get_harmonics().value();
    auto frequency = harmonics.get_frequencies()[harmonicIndex];

    std::vector<double> bobbinPointsX = matplot::linspace(coreColumnWidth / 2, bobbinWidthStart + bobbinWidth, _numberPointsX);
    std::vector<double> bobbinPointsY = matplot::linspace(-coreColumnHeight / 2, coreColumnHeight / 2, _numberPointsY);
    std::vector<OpenMagnetics::FieldPoint> points;
    for (size_t j = 0; j < bobbinPointsY.size(); ++j) {
        for (size_t i = 0; i < bobbinPointsX.size(); ++i) {
            OpenMagnetics::FieldPoint fieldPoint;
            fieldPoint.set_point(std::vector<double>{bobbinPointsX[i], bobbinPointsY[j]});
            points.push_back(fieldPoint);
        }
    }

    Field inducedField;
    inducedField.set_data(points);
    inducedField.set_frequency(frequency);

    OpenMagnetics::MagneticField magneticField;
    magneticField.set_fringing_effect(_includeFringing);
    magneticField.set_mirroring_dimension(_mirroringDimension);
    auto windingWindowMagneticStrengthFieldOutput = magneticField.calculate_magnetic_field_strength_field(operatingPoint, magnetic, inducedField);

    auto field = windingWindowMagneticStrengthFieldOutput.get_field_per_frequency()[0];
    return field;
}

double LeakageInductance::calculate_leakage_inductance(OperatingPoint operatingPoint, MagneticWrapper magnetic, size_t harmonicIndex, std::optional<ComplexField> inputField) {
    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_processed_description()){
        throw std::runtime_error("Bobbin is not processed");
    }
    if (!bobbin.get_processed_description()->get_column_width()){
        throw std::runtime_error("Bobbin is missing column width");
    }
    double bobbinColumnWidth = bobbin.get_processed_description()->get_column_width().value();
    double bobbinColumnDepth = bobbin.get_processed_description()->get_column_depth();
    auto coreColumnShape = magnetic.get_mutable_core().get_columns()[0].get_shape();
    double coreWindingWindowWidth = magnetic.get_mutable_core().get_winding_windows()[0].get_width().value();
    double coreWindingWindowHeight = magnetic.get_mutable_core().get_winding_windows()[0].get_height().value();

    double currentRms = operatingPoint.get_excitations_per_winding()[0].get_current()->get_processed()->get_rms().value();

    ComplexField field;
    if (inputField) {
        field = inputField.value();
    }
    else {
        field = calculate_magnetic_field(operatingPoint, magnetic, harmonicIndex);
    }

    double dx = coreWindingWindowWidth / _numberPointsX;
    double dy = coreWindingWindowHeight / _numberPointsY;
    double dA = dx * dy;
    // std::cout << "dx: " << dx << std::endl;
    // std::cout << "dy: " << dy << std::endl;
    // std::cout << "dA: " << dA << std::endl;

    double vacuumPermeability = Constants().vacuumPermeability;

    double energy = 0;
    for (auto datum : field.get_data()){
        double length;

        if (coreColumnShape == ColumnShape::ROUND) {
            length = 2 * std::numbers::pi * datum.get_point()[0];
        }
        else {
            length = 2 * std::numbers::pi * (datum.get_point()[0] - bobbinColumnWidth) + 2 * 2 * bobbinColumnWidth + 2 * 2 * bobbinColumnDepth;
        }
        // std::cout << "round length: " << (2 * std::numbers::pi * datum.get_point()[0]) << std::endl;
        // std::cout << "squared length: " << (2 * std::numbers::pi * (datum.get_point()[0] - bobbinColumnWidth) + 2 * 2 * bobbinColumnWidth + 2 * 2 * bobbinColumnDepth) << std::endl;
        // std::cout << "length: " << length << std::endl;
        // std::cout << "length: " << length << std::endl;
        double HeSquared = pow(datum.get_real(), 2) + pow(datum.get_imaginary(), 2);
        energy += 0.5 * vacuumPermeability * HeSquared * dA * length;
    }
    // std::cout << "energy" << std::endl;
    // std::cout << energy << std::endl;

    // std::cout << "energy" << std::endl;
    // std::cout << energy << std::endl;

    double leakageInductance = 2.0 / pow(currentRms, 2) * energy;
    // leakageInductance = leakageInductance / pow(magnetic.get_coil().get_functional_description()[0].get_number_parallels(), 2);
    return leakageInductance;
}

} // namespace OpenMagnetics
