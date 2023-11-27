#include "InputsWrapper.h"
#include "MagneticAdviser.h"
#include "MagneticSimulator.h"
#include "Painter.h"
#include "Defaults.h"
#include "../tests/TestingUtils.h"


namespace OpenMagnetics {


    std::vector<MasWrapper> MagneticAdviser::get_advised_magnetic(InputsWrapper inputs, size_t maximumNumberResults) {
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

    std::vector<MasWrapper> MagneticAdviser::get_advised_magnetic(std::vector<CoreWrapper>* cores, std::vector<WireWrapper>* wires, InputsWrapper inputs, size_t maximumNumberResults){
        std::vector<MasWrapper> masData;
        CoreAdviser coreAdviser(false);
        coreAdviser.set_unique_core_shapes(true);
        CoilAdviser coilAdviser;
        MagneticSimulator magneticSimulator;

        std::vector<int> interleavingCombinations = {1, 2};
        if (inputs.get_design_requirements().get_turns_ratios().size() == 0) {
            interleavingCombinations = {1};
        }

        auto masMagneticsWithCore = coreAdviser.get_advised_core(inputs, cores, maximumNumberResults * 2);
        for (auto& masMagneticWithCore : masMagneticsWithCore) {
            std::cout << masMagneticWithCore.first.get_mutable_magnetic().get_mutable_core().get_shape_name() << std::endl;
            for (auto interleavingCombination : interleavingCombinations) {
                coilAdviser.set_interleaving_level(interleavingCombination);
                auto masMagneticsWithCoreAndCoil = coilAdviser.get_advised_coil(wires, masMagneticWithCore.first, int(maximumNumberResults / 2));
        std::cout << "masMagneticsWithCoreAndCoil.size()" << std::endl;
        std::cout << masMagneticsWithCoreAndCoil.size() << std::endl;
                for (auto masMagnetic : masMagneticsWithCoreAndCoil) {

                    masMagnetic.first = magneticSimulator.simulate(masMagnetic.first);

                    masData.push_back(masMagnetic.first);
                    if (masData.size() == maximumNumberResults) {
                        break;
                    }
                }
                if (masData.size() == maximumNumberResults) {
                    break;
                }
            }
            if (masData.size() == maximumNumberResults) {
                break;
            }
        }
        std::cout << "masData.size()" << std::endl;
        std::cout << masData.size() << std::endl;

        sort(masData.begin(), masData.end(), [](MasWrapper& b1, MasWrapper& b2) {
            return b1.get_outputs()[0].get_core_losses().value().get_core_losses() + b1.get_outputs()[0].get_winding_losses().value().get_winding_losses() <
             b2.get_outputs()[0].get_core_losses().value().get_core_losses() + b2.get_outputs()[0].get_winding_losses().value().get_winding_losses();
        });

        for (auto masMagnetic : masData) {

            std::cout << "losses" << std::endl;
            std::cout << (masMagnetic.get_outputs()[0].get_core_losses().value().get_core_losses() + masMagnetic.get_outputs()[0].get_winding_losses().value().get_winding_losses()) << std::endl;
        }

        if (masData.size() > maximumNumberResults) {
            masData = std::vector<MasWrapper>(masData.begin(), masData.end() - (masData.size() - maximumNumberResults));
        }

        return masData;

    }

void preview_magnetic(MasWrapper mas) {
    std::string text = "";
    text += "Core shape: " + mas.get_mutable_magnetic().get_mutable_core().get_shape_name() + "\n";
    text += "Core material: " + mas.get_mutable_magnetic().get_mutable_core().get_material_name() + "\n";
    text += "Core gap: " + std::to_string(mas.get_mutable_magnetic().get_mutable_core().get_functional_description().get_gapping()[0].get_length()) + "\n";
    for (size_t windingIndex = 0; windingIndex < mas.get_mutable_magnetic().get_mutable_coil().get_functional_description().size(); ++windingIndex) {
        auto winding = mas.get_mutable_magnetic().get_mutable_coil().get_functional_description()[windingIndex];
        auto wire = OpenMagnetics::CoilWrapper::resolve_wire(winding);
        text += "Winding: " + winding.get_name() + "\n";
        text += "\tNumber Turns: " + std::to_string(winding.get_number_turns()) + "\n";
        text += "\tNumber Parallels: " + std::to_string(winding.get_number_parallels()) + "\n";
        text += "\tWire: " + std::string(magic_enum::enum_name(wire.get_type())) + " " + std::string(magic_enum::enum_name(wire.get_standard().value())) + " " + wire.get_name().value() + "\n";
    }

    for (size_t operatingPointIndex = 0; operatingPointIndex < mas.get_outputs().size(); ++operatingPointIndex) {
        auto output = mas.get_outputs()[operatingPointIndex];
        text += "Operating Point: " + std::to_string(operatingPointIndex + 1) + "\n";
        text += "\tCore losses: " + std::to_string(output.get_core_losses().value().get_core_losses()) + "\n";
        text += "\tCore temperature: " + std::to_string(output.get_core_losses().value().get_temperature().value()) + "\n";
        text += "\tWinding losses: " + std::to_string(output.get_winding_losses().value().get_winding_losses()) + "\n";
        for (size_t windingIndex = 0; windingIndex < output.get_winding_losses().value().get_winding_losses_per_winding().value().size(); ++windingIndex) {
            auto winding = mas.get_mutable_magnetic().get_mutable_coil().get_functional_description()[windingIndex];
            auto windingLosses = output.get_winding_losses().value().get_winding_losses_per_winding().value()[windingIndex];
            text += "\t\tLosses for winding: " + winding.get_name() + "\n";
            double skinEffectLosses = 0;
            double proximityEffectLosses = 0;
            for (size_t i = 0; i < windingLosses.get_skin_effect_losses().value().get_losses_per_harmonic().size(); ++i)
            {
                skinEffectLosses += windingLosses.get_skin_effect_losses().value().get_losses_per_harmonic()[i];
            }
            for (size_t i = 0; i < windingLosses.get_proximity_effect_losses().value().get_losses_per_harmonic().size(); ++i)
            {
                proximityEffectLosses += windingLosses.get_proximity_effect_losses().value().get_losses_per_harmonic()[i];
            }


            text += "\t\t\tOhmic losses: " + std::to_string(windingLosses.get_ohmic_losses().value().get_losses()) + "\n";
            text += "\t\t\tSkin effect losses: " + std::to_string(skinEffectLosses) + "\n";
            text += "\t\t\tProximity effect losses: " + std::to_string(proximityEffectLosses) + "\n";
        }
    }
    std::cout << text << std::endl;
}

} // namespace OpenMagnetics

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        throw std::invalid_argument("Missing inputs file");
    }
    else {
        int numberMagnetics = 1;
        std::filesystem::path inputFilepath = argv[1];
        std::string filePath = __FILE__;
        std::filesystem::path outputFilepath = filePath.substr(0, filePath.rfind("/")).append("/../output/");
        if (argc >= 3) {
            if (is_number(argv[2])) {
                numberMagnetics = std::stoi(argv[2]);
            }
            else {
                outputFilepath = argv[2];
            }
        }
        if (argc >= 4) {
            if (!is_number(argv[2]) && is_number(argv[3])) {
                numberMagnetics = std::stoi(argv[3]);
            }
        }

        std::ifstream f(inputFilepath);
        json masJson;
        std::string str;
        if(f) {
            std::ostringstream ss;
            ss << f.rdbuf(); // reading data
            masJson = json::parse(ss.str());
        }

        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "Group IIIa") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "Group IIIA";
        }
        if (masJson["inputs"]["designRequirements"]["insulation"]["cti"] == "Group IIIb") {
            masJson["inputs"]["designRequirements"]["insulation"]["cti"] = "Group IIIB";
        }

        OpenMagnetics::InputsWrapper inputs(masJson["inputs"]);

        auto outputFolder = inputFilepath.parent_path();

        std::cout << "numberMagnetics: " << numberMagnetics << std::endl; 
        std::cout << "inputFilepath: " << inputFilepath << std::endl;
        std::cout << "outputFilepath: " << outputFilepath << std::endl; 
        std::cout << "filename(): " << inputFilepath.filename() << std::endl; 
        // std::cout << "str: " << str << std::endl; 

        OpenMagnetics::MagneticAdviser MagneticAdviser;
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, numberMagnetics);
        for (size_t i = 0; i < masMagnetics.size(); ++i){
            preview_magnetic(masMagnetics[i]);

            std::filesystem::path outputFilename = outputFilepath;
            outputFilename += inputFilepath.filename();
            std::cout << "outputFilename: " << outputFilename << std::endl; 
            outputFilename += "_design_" + std::to_string(i) + ".json";
            std::cout << "outputFilename: " << outputFilename << std::endl; 
            std::ofstream outputFile(outputFilename);
            
            json output;
            to_json(output, masMagnetics[i]);
            outputFile << output;
            
            outputFile.close();

            outputFilename.replace_extension("svg");
            OpenMagnetics::Painter painter(outputFilename, OpenMagnetics::Painter::PainterModes::CONTOUR);

            painter.set_number_points_x(20);
            painter.set_number_points_y(20);
            painter.set_fringing_effect(false);
            painter.set_mirroring_dimension(0);
            painter.paint_magnetic_field(masMagnetics[i].get_mutable_inputs().get_operating_point(0), masMagnetics[i].get_magnetic());
            painter.paint_core(masMagnetics[i].get_magnetic());
            painter.paint_bobbin(masMagnetics[i].get_magnetic());
            painter.paint_coil_turns(masMagnetics[i].get_magnetic());
            painter.export_svg();
        }


    }
}