#include "Painter.h"
#include "CoreAdviser.h"
#include "CoilAdviser.h"
#include "InputsWrapper.h"
#include "TestingUtils.h"

#include <UnitTest++.h>
#include <vector>

SUITE(CoilAdviser) {


    TEST(Test_CoilAdviser_Random) {
        srand (time(NULL));

        int count = 10;
        while (count > 0) {
            std::vector<double> turnsRatios;
            int64_t numberStacks = 1;

            std::vector<int64_t> numberTurns;
            std::vector<int64_t> numberParallels;
            int64_t numberPhysicalTurns = std::numeric_limits<int64_t>::max();
            // for (size_t windingIndex = 0; windingIndex < std::rand() % 2 + 1UL; ++windingIndex)
            for (size_t windingIndex = 0; windingIndex < std::rand() % 4 + 1UL; ++windingIndex)
            {
                numberTurns.push_back(std::rand() % 100 + 1L);
                numberParallels.push_back(1);
                numberPhysicalTurns = std::min(numberPhysicalTurns, numberTurns.back() * numberParallels.back());
            }
            for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
                turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
            }

            double frequency = std::rand() % 1000000 + 10000L;
            double magnetizingInductance = double(std::rand() % 10000) *  1e-6;
            double temperature = 25;
            OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
            double peakToPeak = std::rand() % 30 + 1L;
            double dutyCycle = double(std::rand() % 100) / 100;
            double dcCurrent = 0;
            if (numberTurns.size() > 1) {
                dcCurrent = std::rand() % 30;
            }
            int interleavingLevel = std::rand() % 5 + 1;
            interleavingLevel = std::min(int(numberPhysicalTurns), interleavingLevel);
            auto coreShapeNames = OpenMagnetics::get_shape_names(false);
            std::string coreShapeName;
            OpenMagnetics::MagneticWrapper magnetic;

            double gapLength = 1;
            double columnHeight = 1;
            int numberGaps = 1;
            int gappingTypeIndex = std::rand() % 4;
            OpenMagnetics::GappingType gappingType = magic_enum::enum_cast<OpenMagnetics::GappingType>(gappingTypeIndex).value();
            if (gappingType == OpenMagnetics::GappingType::DISTRIBUTED) {
                numberGaps = std::rand() % 3;
                numberGaps *= 2;
                numberGaps += 3;
            }

            while (true) {
                coreShapeName = coreShapeNames[rand() % coreShapeNames.size()];
                auto shape = OpenMagnetics::find_core_shape_by_name(coreShapeName);
                if (shape.get_family() == OpenMagnetics::CoreShapeFamily::PQI || 
                    shape.get_family() == OpenMagnetics::CoreShapeFamily::UI || 
                    shape.get_family() == OpenMagnetics::CoreShapeFamily::EI) {
                    continue;
                }

                do {
                    auto shape = OpenMagnetics::find_core_shape_by_name(coreShapeName);
                    auto corePiece = OpenMagnetics::CorePiece::factory(shape);
                    gapLength = double(rand() % 10000 + 1.0) / 1000000;
                    columnHeight = corePiece->get_columns()[0].get_height();
                }
                while (columnHeight < (gapLength * numberGaps));

                std::vector<OpenMagnetics::CoreGap> gapping;
                switch (gappingType) {
                    case OpenMagnetics::GappingType::GRINDED:
                        gapping = OpenMagneticsTesting::get_grinded_gap(gapLength);
                        break;
                    case OpenMagnetics::GappingType::SPACER:
                        gapping = OpenMagneticsTesting::get_spacer_gap(gapLength);
                        break;
                    case OpenMagnetics::GappingType::RESIDUAL:
                        gapping = OpenMagneticsTesting::get_residual_gap();
                        break;
                    case OpenMagnetics::GappingType::DISTRIBUTED:
                        gapping = OpenMagneticsTesting::get_distributed_gap(gapLength, numberGaps);
                        break;
                }
                
                magnetic = OpenMagneticsTesting::get_quick_magnetic(coreShapeName, gapping, numberTurns, numberStacks, "3C91");
                break;
            }

            int turnsAlignmentIndex = std::rand() % 4;
            int windingOrientationIndex = (std::rand() % 2) * 2;  // To avoid 1, which is RADIAL
            int layersOrientationIndex = (std::rand() % 2) * 2;  // To avoid 1, which is RADIAL
            OpenMagnetics::CoilAlignment turnsAlignment = magic_enum::enum_cast<OpenMagnetics::CoilAlignment>(turnsAlignmentIndex).value();
            OpenMagnetics::WindingOrientation windingOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(windingOrientationIndex).value();
            OpenMagnetics::WindingOrientation layersOrientation = magic_enum::enum_cast<OpenMagnetics::WindingOrientation>(layersOrientationIndex).value();


            // for (size_t windingIndex = 0; windingIndex < numberTurns.size(); ++windingIndex) {
            //     std::cout << "numberTurns: " << numberTurns[windingIndex] << std::endl;
            // }
            // std::cout << "frequency: " << frequency << std::endl;
            // std::cout << "peakToPeak: " << peakToPeak << std::endl;
            // std::cout << "magnetizingInductance: " << magnetizingInductance << std::endl;
            // std::cout << "dutyCycle: " << dutyCycle << std::endl;
            // std::cout << "dcCurrent: " << dcCurrent << std::endl;
            // std::cout << "interleavingLevel: " << interleavingLevel << std::endl;
            // std::cout << "coreShapeName: " << coreShapeName << std::endl;
            // std::cout << "gappingTypeIndex: " << gappingTypeIndex << std::endl;
            // std::cout << "gapLength: " << gapLength << std::endl;
            // std::cout << "numberGaps: " << numberGaps << std::endl;
            // std::cout << "magnetic.get_mutable_core().get_shape_family(): " << magic_enum::enum_name(magnetic.get_mutable_core().get_shape_family()) << std::endl;
            // std::cout << "turnsAlignmentIndex: " << turnsAlignmentIndex << std::endl;
            // std::cout << "windingOrientationIndex: " << windingOrientationIndex << std::endl;
            // std::cout << "layersOrientationIndex: " << layersOrientationIndex << std::endl;

            magnetic.get_mutable_coil().set_turns_alignment(turnsAlignment);
            magnetic.get_mutable_coil().set_winding_orientation(windingOrientation);
            magnetic.get_mutable_coil().set_layers_orientation(layersOrientation);

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
            masMagnetic.set_inputs(inputs);
            masMagnetic.set_magnetic(magnetic);

            OpenMagnetics::CoilAdviser coilAdviser;
            coilAdviser.set_interleaving_level(interleavingLevel);
            auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

            if (masMagneticsWithCoil.size() > 0) {
                auto masMagneticWithCoil = masMagneticsWithCoil[0].first;
                if (!masMagneticWithCoil.get_magnetic().get_coil().get_turns_description()) {
                    continue;
                }
                OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
                std::string filePath = __FILE__;
                auto randNumber =  std::rand();
                // auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");
                // {
                //     auto outFile = outputFilePath;
                //     std::string filename = "Test_CoilAdviser" + std::to_string(randNumber) + "_Quiver.svg";
                //     outFile.append(filename);
                //     OpenMagnetics::Painter painter(outFile, OpenMagnetics::Painter::PainterModes::QUIVER);
                //     painter.set_number_points_x(50);
                //     painter.set_number_points_y(50);
                //     painter.set_logarithmic_scale(false);
                //     painter.set_fringing_effect(true);
                //     painter.set_maximum_scale_value(std::nullopt);
                //     painter.set_minimum_scale_value(std::nullopt);
                //     painter.paint_magnetic_field(inputs.get_operating_point(0), masMagneticWithCoil.get_magnetic());
                //     painter.paint_core(masMagneticWithCoil.get_magnetic());
                //     painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
                //     painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
                //     painter.export_svg();
                // }
                // {
                //     auto outFile = outputFilePath;
                //     std::string filename = "Test_CoilAdviser" + std::to_string(randNumber) + "_Quiver_Log.svg";
                //     outFile.append(filename);
                //     OpenMagnetics::Painter painter(outFile, OpenMagnetics::Painter::PainterModes::QUIVER);
                //     painter.set_number_points_x(50);
                //     painter.set_number_points_y(50);
                //     painter.set_logarithmic_scale(true);
                //     painter.set_fringing_effect(true);
                //     painter.set_maximum_scale_value(std::nullopt);
                //     painter.set_minimum_scale_value(std::nullopt);
                //     painter.paint_magnetic_field(inputs.get_operating_point(0), masMagneticWithCoil.get_magnetic());
                //     painter.paint_core(masMagneticWithCoil.get_magnetic());
                //     painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
                //     painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
                //     painter.export_svg();
                // }
                // {
                //     auto outFile = outputFilePath;
                //     std::string filename = "Test_CoilAdviser" + std::to_string(randNumber) + "_Contour.svg";
                //     outFile.append(filename);
                //     OpenMagnetics::Painter painter(outFile, OpenMagnetics::Painter::PainterModes::CONTOUR);
                //     painter.set_number_points_x(50);
                //     painter.set_number_points_y(50);
                //     painter.set_logarithmic_scale(false);
                //     painter.set_fringing_effect(true);
                //     painter.set_maximum_scale_value(std::nullopt);
                //     painter.set_minimum_scale_value(std::nullopt);
                //     painter.paint_magnetic_field(inputs.get_operating_point(0), masMagneticWithCoil.get_magnetic());
                //     painter.paint_core(masMagneticWithCoil.get_magnetic());
                //     painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
                //     painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
                //     painter.export_svg();
                // }
                // {
                //     auto outFile = outputFilePath;
                //     std::string filename = "Test_CoilAdviser" + std::to_string(randNumber) + "_Contour_Log.svg";
                //     outFile.append(filename);
                //     OpenMagnetics::Painter painter(outFile, OpenMagnetics::Painter::PainterModes::CONTOUR);
                //     painter.set_number_points_x(50);
                //     painter.set_number_points_y(50);
                //     painter.set_logarithmic_scale(true);
                //     painter.set_fringing_effect(true);
                //     painter.set_maximum_scale_value(std::nullopt);
                //     painter.set_minimum_scale_value(std::nullopt);
                //     painter.paint_magnetic_field(inputs.get_operating_point(0), masMagneticWithCoil.get_magnetic());
                //     painter.paint_core(masMagneticWithCoil.get_magnetic());
                //     painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
                //     painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
                //     painter.export_svg();
                // }
                count--;
            }
        }

    }

    TEST(Test_CoilAdviser_Random_0) {
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.003);
        std::vector<double> turnsRatios;
        int64_t numberStacks = 1;


        double magnetizingInductance = 10e-6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        std::vector<int64_t> numberTurns = {82, 5};
        std::vector<int64_t> numberParallels = {1, 1};
        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }
        double frequency = 675590;
        double peakToPeak = 26;

        auto magnetic = OpenMagneticsTesting::get_quick_magnetic("EP 20", gapping, numberTurns, numberStacks, "3C91");
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
        masMagnetic.set_inputs(inputs);
        masMagnetic.set_magnetic(magnetic);

        OpenMagnetics::CoilAdviser coilAdviser;
        auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

        if (masMagneticsWithCoil.size() > 0) {
            auto masMagneticWithCoil = masMagneticsWithCoil[0].first;
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
        }
    }
    TEST(Test_CoilAdviser_Random_1) {
        std::vector<double> turnsRatios;
        int64_t numberStacks = 1;


        double magnetizingInductance = 10e-6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double dutyCycle = 0.5;
        double dcCurrent = 0;

        std::vector<int64_t> numberTurns = {16, 34};
        std::vector<int64_t> numberParallels = {1, 1};
        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }
        double frequency = 811022;
        double peakToPeak = 1;
        int interleavingLevel = 3;

        OpenMagnetics::MagneticWrapper magnetic;
        std::string coreShapeName  = "P 7.35X3.6";
        double gapLength = 1;
        double columnHeight = 1;
        do {
            auto shape = OpenMagnetics::find_core_shape_by_name(coreShapeName);
            auto corePiece = OpenMagnetics::CorePiece::factory(shape);
            gapLength = double(rand() % 10000 + 1.0) / 1000000;
            columnHeight = corePiece->get_columns()[0].get_height();
        }
        while (columnHeight < gapLength);
        auto gapping = OpenMagneticsTesting::get_grinded_gap(gapLength);
        magnetic = OpenMagneticsTesting::get_quick_magnetic(coreShapeName, gapping, numberTurns, numberStacks, "3C91");
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
        masMagnetic.set_inputs(inputs);
        masMagnetic.set_magnetic(magnetic);

        OpenMagnetics::CoilAdviser coilAdviser;
        coilAdviser.set_interleaving_level(interleavingLevel);
        auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

        if (masMagneticsWithCoil.size() > 0) {
            auto masMagneticWithCoil = masMagneticsWithCoil[0].first;
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            std::string filePath = __FILE__;
            auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            // OpenMagnetics::Painter painter(outFile);

            // painter.paint_core(masMagneticWithCoil.get_magnetic());
            // painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
            // painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
            // painter.export_svg();
        }
    }

    TEST(Test_CoilAdviser_Random_2) {
        srand (time(NULL));

        std::vector<double> turnsRatios;
        int64_t numberStacks = 1;

        std::vector<int64_t> numberTurns = {24, 78, 76};
        std::vector<int64_t> numberParallels = {1, 1, 1};

        for (size_t windingIndex = 1; windingIndex < numberTurns.size(); ++windingIndex) {
            turnsRatios.push_back(double(numberTurns[0]) / numberTurns[windingIndex]);
        }

        double frequency = 507026;
        double magnetizingInductance = 10e-6;
        double temperature = 25;
        OpenMagnetics::WaveformLabel waveShape = OpenMagnetics::WaveformLabel::SINUSOIDAL;
        double peakToPeak = 10;
        double dutyCycle = 0.5;
        double dcCurrent = 0;
        int interleavingLevel = 4;
        auto coreShapeNames = OpenMagnetics::get_shape_names(false);
        std::string coreShapeName;
        OpenMagnetics::MagneticWrapper magnetic;

        auto gapping = OpenMagneticsTesting::get_residual_gap();
        magnetic = OpenMagneticsTesting::get_quick_magnetic("ETD 19", gapping, numberTurns, numberStacks, "3C91");


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
        masMagnetic.set_inputs(inputs);
        masMagnetic.set_magnetic(magnetic);

        OpenMagnetics::CoilAdviser coilAdviser;
        coilAdviser.set_interleaving_level(interleavingLevel);
        auto masMagneticsWithCoil = coilAdviser.get_advised_coil(masMagnetic, 2);

        if (masMagneticsWithCoil.size() > 0) {
            auto masMagneticWithCoil = masMagneticsWithCoil[0].first;
            OpenMagneticsTesting::check_turns_description(masMagneticWithCoil.get_magnetic().get_coil());
            std::string filePath = __FILE__;
            auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");
            auto outFile = outputFilePath;
            std::string filename = "Test_CoilAdviser" + std::to_string(std::rand()) + ".svg";
            outFile.append(filename);
            // OpenMagnetics::Painter painter(outFile);

            // painter.paint_core(masMagneticWithCoil.get_magnetic());
            // painter.paint_bobbin(masMagneticWithCoil.get_magnetic());
            // painter.paint_coil_turns(masMagneticWithCoil.get_magnetic());
            // painter.export_svg();
        }

    }
}