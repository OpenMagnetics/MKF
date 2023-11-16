#include "Painter.h"
#include "MagneticAdviser.h"
#include "InputsWrapper.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <vector>

SUITE(MagneticAdviser) {

    TEST(MagneticAdviser) {
        srand (time(NULL));

        std::vector<double> turnsRatios;

        std::vector<int64_t> numberTurns = {24, 78, 76};
        std::vector<int64_t> numberParallels = {1, 1, 1};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 100e-6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::TRIANGULAR;
        double peakToPeak = 100;
        double dutyCycle = 0.5;
        double dcCurrent = 0;


        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point_only_current(frequency,
                                                                                              magnetizingInductance,
                                                                                              temperature,
                                                                                              waveShape,
                                                                                              peakToPeak,
                                                                                              dutyCycle,
                                                                                              dcCurrent,
                                                                                              turnsRatios);

        OpenMagnetics::MasWrapper masMagnetic;
        inputs.process_waveforms();

        OpenMagnetics::MagneticAdviser MagneticAdviser;
        // MagneticAdviser.set_interleaving_level(interleavingLevel);
        auto masMagnetics = MagneticAdviser.get_advised_magnetic(inputs, 1);
        std::cout << masMagnetics.size() << std::endl;
        std::cout << masMagnetics[0].first.get_mutable_magnetic().get_mutable_core().get_shape_name() << std::endl;
        // std::cout << masMagnetics[0].first.get_magnetic().get_coil().get_manufacturer_info().value().get_reference().value() << std::endl;

        if (masMagnetics.size() > 0) {
            auto masMagnetic = masMagnetics[0].first;
            OpenMagneticsTesting::check_turns_description(masMagnetic.get_magnetic().get_coil());
            std::string filePath = __FILE__;
            auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");
            auto outFile = outputFilePath;
            std::string filename = "MagneticAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            OpenMagnetics::Painter painter(outFile);

            painter.paint_core(masMagnetic.get_magnetic());
            painter.paint_bobbin(masMagnetic.get_magnetic());
            painter.paint_coil_turns(masMagnetic.get_magnetic());
            painter.export_svg();
        }

    }
}