#include "InputsWrapper.h"
#include "MagneticAdviser.h"
#include "Defaults.h"
#include "../tests/TestingUtils.h"


namespace OpenMagnetics {


    std::vector<std::pair<MasWrapper, double>> MagneticAdviser::get_advised_magnetic(InputsWrapper inputs, size_t maximumNumberResults) {
        std::vector<WireWrapper> wires;
        std::vector<CoreWrapper> cores;
        std::string file_path = __FILE__;

        {
            auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/../../MAS/data/wires.ndjson");
            std::ifstream ndjsonFile(inventory_path);
            std::string jsonLine;
            while (std::getline(ndjsonFile, jsonLine)) {
                json jf = json::parse(jsonLine);
                WireWrapper wire(jf);
                if ((_includeFoil || wire.get_type() != WireType::FOIL) &&
                    (_includeRectangular || wire.get_type() != WireType::RECTANGULAR) &&
                    (_includeLitz || wire.get_type() != WireType::LITZ) &&
                    (_includeRound || wire.get_type() != WireType::ROUND)) {
                    wires.push_back(wire);
                }
            }
        }
        {
            auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/../../MAS/data/cores_stock.ndjson");
            std::ifstream ndjsonFile(inventory_path);
            std::string jsonLine;
            while (std::getline(ndjsonFile, jsonLine)) {
                json jf = json::parse(jsonLine);
                CoreWrapper core(jf, false, true, false);
                if (_includeToroids || core.get_type() != CoreType::TOROIDAL) {
                    cores.push_back(core);
                }
            }
        }
        return get_advised_magnetic(&cores, &wires, inputs, maximumNumberResults);
    }

    std::vector<std::pair<MasWrapper, double>> MagneticAdviser::get_advised_magnetic(std::vector<CoreWrapper>* cores, std::vector<WireWrapper>* wires, InputsWrapper inputs, size_t maximumNumberResults){
        std::vector<std::pair<MasWrapper, double>> magnetics;
        OpenMagnetics::CoreAdviser coreAdviser(false);
        OpenMagnetics::CoilAdviser coilAdviser;
        coilAdviser.set_interleaving_level(1);  // TODO add more combinations

        auto masMagneticsWithCore = coreAdviser.get_advised_core(inputs, cores, 1);
        if (masMagneticsWithCore.size() > 0) {
            for (auto& masMagneticWithCore : masMagneticsWithCore) {
                auto masMagneticsWithCoreAndCoil = coilAdviser.get_advised_coil(wires, masMagneticWithCore.first, 1);
                std::move(masMagneticsWithCoreAndCoil.begin(), masMagneticsWithCoreAndCoil.end(), std::back_inserter(magnetics));
            }
        }

        return magnetics;

    }
} // namespace OpenMagnetics
