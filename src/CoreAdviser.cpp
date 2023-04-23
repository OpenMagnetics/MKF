#include "CoreAdviser.h"
#include "MagnetizingInductance.h"
#include "MagneticEnergy.h"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <execution>


namespace OpenMagnetics {

std::vector<Magnetic> CoreAdviser::MagneticCoreFilterAreaProduct::filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models) {
    return unfilteredMagnetics;
}

std::vector<Magnetic> CoreAdviser::MagneticCoreFilterEnergyStored::filter_magnetics(std::vector<Magnetic> unfilteredMagnetics, InputsWrapper inputs, std::map<std::string, std::string> models) {
    OpenMagnetics::MagneticEnergy magneticEnergy(models);
    std::vector<Magnetic> filteredMagnetics;
    double requiredMagneticEnergy = magneticEnergy.required_magnetic_energy(inputs).get_nominal().value();
    auto operationPoint = inputs.get_operation_point(0);

    for (auto& magnetic : unfilteredMagnetics){  
        double totalStorableMagneticEnergy = magneticEnergy.get_core_maximum_magnetic_energy(static_cast<CoreWrapper>(magnetic.get_core()), &operationPoint);
        if (totalStorableMagneticEnergy >= requiredMagneticEnergy) {
            filteredMagnetics.push_back(magnetic);
        }
    }

    return filteredMagnetics;
}


WindingWrapper get_dummy_winding() {
    json primaryWindingJson = json();
    primaryWindingJson["isolationSide"] = IsolationSide::PRIMARY;
    primaryWindingJson["name"] = "primary";
    primaryWindingJson["numberParallels"] = 1;
    primaryWindingJson["numberTurns"] = 1;
    primaryWindingJson["wire"] = "Dummy";
    WindingFunctionalDescription primaryWindingFunctionalDescription(primaryWindingJson);

    WindingWrapper winding;
    winding.set_functional_description({primaryWindingFunctionalDescription});
    return winding;
}

CoreWrapper CoreAdviser::get_advised_core(InputsWrapper inputs) {
    std::string file_path = __FILE__;
    auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/PyMKF/inventory.ndjson");
    std::ifstream ndjsonFile(inventory_path);
    std::string jsonLine;
    std::vector<CoreWrapper> cores;
    while (std::getline(ndjsonFile, jsonLine)) {
        json jf = json::parse(jsonLine);
        CoreWrapper core(jf, false, true, false);
        cores.push_back(core);
    }

    std::vector<Magnetic> magnetics;
    MagnetizingInductance magnetizing_inductance(
        std::map<std::string, std::string>({{"gapReluctance", "ZHANG"}}));
    WindingWrapper winding = get_dummy_winding();

    Magnetic magnetic;
    magnetic.set_winding(std::move(winding));
    for (auto& core : cores){  
        magnetic.set_core(std::move(core));
        magnetics.push_back(std::move(magnetic));
    }

    magnetics = MagneticCoreFilterAreaProduct::filter_magnetics(magnetics, inputs, _models);
    magnetics = MagneticCoreFilterEnergyStored::filter_magnetics(magnetics, inputs, _models);

    for (size_t i = 0; i < magnetics.size(); ++i){
        double numberTurns = magnetizing_inductance.get_number_turns_from_gapping_and_inductance(cores[i], &inputs);
        magnetics[i].get_mutable_winding().get_mutable_functional_description()[0].set_number_turns(numberTurns);
    }


    std::cout << cores.size() << std::endl;
    std::cout << magnetics.size() << std::endl;

    return cores[0];
}

} // namespace OpenMagnetics
