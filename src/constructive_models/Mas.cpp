#include <MAS.hpp>
#include "constructive_models/Mas.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/Utils.h"

namespace OpenMagnetics {


void from_file(std::filesystem::path filepath, Mas & x) {
    std::ifstream f(filepath);
    std::string data((std::istreambuf_iterator<char>(f)),
                     std::istreambuf_iterator<char>());
    auto masJson = json::parse(data);
    auto inputsJson = masJson["inputs"];
    auto magneticJson = masJson["magnetic"];
    auto outputsJson = masJson["outputs"];
    for (auto layerJson : magneticJson["coil"]["layersDescription"]) {
        auto layer = MAS::Layer(layerJson);
    }
    auto coil = OpenMagnetics::Coil(magneticJson["coil"]);
    auto magnetic = OpenMagnetics::Magnetic(magneticJson);
    std::vector<OpenMagnetics::Outputs> outputs;
    if (outputsJson != nullptr) {
        for (auto outputJson : outputsJson) {
            outputs.push_back(OpenMagnetics::Outputs(outputJson));
        }
    }
    OpenMagnetics::Inputs inputs;
    std::vector<double> magnetizingInductancePerPoint;
    for (auto output : outputs) {
        if (output.get_magnetizing_inductance()) {
            auto magnetizingInductance = resolve_dimensional_values(output.get_magnetizing_inductance()->get_magnetizing_inductance());
            magnetizingInductancePerPoint.push_back(magnetizingInductance);
        }
    }

    if (magnetizingInductancePerPoint.size() > 0) {
        inputs = OpenMagnetics::Inputs(inputsJson, true, magnetizingInductancePerPoint);

    }
    else {
        try {
            MagnetizingInductance magnetizingInductanceModel;
            double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_core(), magnetic.get_coil()).get_magnetizing_inductance().get_nominal().value();
            inputs = OpenMagnetics::Inputs(inputsJson, true, magnetizingInductance);
        }
        catch (const std::exception &e)
        {
            (void)e; // Suppress unused variable warning
            inputs = OpenMagnetics::Inputs(inputsJson, true);
        }
    }

    x.set_inputs(inputs);
    x.set_magnetic(magnetic);
    x.set_outputs(outputs);
}

} // namespace OpenMagnetics
