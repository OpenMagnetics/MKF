#include "constructive_models/MasMigration.h"
#include <MAS.hpp>
#include "constructive_models/Mas.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/Utils.h"

namespace OpenMagnetics {

// ABT #1876 (Alf, 2026-10-14): an output's parts that were not computed are ABSENT, never null. MAS.hpp writes a
// disengaged optional as null, which the outputs schema does not accept; a point with no thermal operating point
// (or a fast-mode simulation) has no "temperature" key at all.
void to_json(json & j, const Mas & x) {
    j = json::object();
    j["inputs"] = x.get_inputs();
    j["magnetic"] = x.get_magnetic();
    json outputs = json::array();
    for (const auto& output : x.get_outputs()) {
        json outputJson = output;
        for (auto entry = outputJson.begin(); entry != outputJson.end();) {
            if (entry->is_null()) {
                entry = outputJson.erase(entry);
            }
            else {
                ++entry;
            }
        }
        outputs.push_back(outputJson);
    }
    j["outputs"] = outputs;
}


void from_file(std::filesystem::path filepath, Mas & x) {
    std::ifstream f(filepath);
    std::string data((std::istreambuf_iterator<char>(f)),
                     std::istreambuf_iterator<char>());
    auto masJson = json::parse(data);
    OpenMagnetics::compat::migrate_pre_1_0(masJson);
    auto inputsJson = masJson["inputs"];
    auto magneticJson = masJson["magnetic"];
    auto outputsJson = masJson["outputs"];
    for (auto layerJson : magneticJson["coil"]["layersDescription"]) {
        auto layer = MAS::Layer(layerJson);
    }
    // The coil is built (and wound) by the Magnetic, after it has its core. A coil wound here on
    // its own did not know the core its bobbin's bore follows, and was thrown away anyway.
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
        if (output.get_inductance()) {
            auto magnetizingInductance = resolve_dimensional_values(output.get_inductance()->get_magnetizing_inductance().get_magnetizing_inductance());
            magnetizingInductancePerPoint.push_back(magnetizingInductance);
        }
    }

    if (magnetizingInductancePerPoint.size() > 0) {
        inputs = OpenMagnetics::Inputs(inputsJson, true, magnetizingInductancePerPoint);

    }
    else {
        // A file whose magnetizing inductance cannot be computed is broken input;
        // loading it silently without the inductance hid real core/coil defects.
        MagnetizingInductance magnetizingInductanceModel;
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance().get_nominal().value();
        inputs = OpenMagnetics::Inputs(inputsJson, true, magnetizingInductance);
    }

    x.set_inputs(inputs);
    x.set_magnetic(magnetic);
    x.set_outputs(outputs);
}

} // namespace OpenMagnetics
