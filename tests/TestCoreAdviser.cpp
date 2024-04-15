#include "Settings.h"
#include "CoreAdviser.h"
#include "Utils.h"
#include "InputsWrapper.h"

#include <UnitTest++.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>


SUITE(CoreAdviser) {
    auto settings = OpenMagnetics::Settings::GetInstance();

    void prepare_test_parameters(double dcCurrent, double ambientTemperature, double frequency, std::vector<double> turnsRatios,
                                 double desiredMagnetizingInductance, OpenMagnetics::InputsWrapper& inputs,
                                 double peakToPeak = 20, double dutyCycle = 0.5) {
        inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(
            frequency, desiredMagnetizingInductance, ambientTemperature, OpenMagnetics::WaveformLabel::SINUSOIDAL,
            peakToPeak, dutyCycle, dcCurrent, turnsRatios);
    }

    std::vector<OpenMagnetics::CoreWrapper> load_test_data() {
        std::string file_path = __FILE__;
        auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/testData/test_cores.ndjson");
        std::ifstream ndjsonFile(inventory_path);
        std::string jsonLine;
        std::vector<OpenMagnetics::CoreWrapper> cores;
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            try {

                OpenMagnetics::CoreWrapper core(jf, false, true, false);
                cores.push_back(core);
            }
            catch(const std::runtime_error& re)
            {
                continue;
                // std::cerr << "Error occurred: " << ex.what() << std::endl;
            }
        }

        return cores;
    }

    TEST(Test_All_Cores) {
        OpenMagnetics::clear_databases();
        double voltagePeakToPeak = 600;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &cores);


        CHECK(masMagnetics.size() == 1);
        double bestScoring = masMagnetics[0].second;
        for (size_t i = 0; i < masMagnetics.size(); ++i)
        {
            CHECK(masMagnetics[i].second <= bestScoring);
        }

        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_name() == "T 18/9.0/7.1 - Kool Mµ Hƒ 40 - Ungapped");
        settings->reset();
    }


    TEST(Test_All_Cores_Load_Internally_Only_Stock) {
        OpenMagnetics::clear_databases();
        double voltagePeakToPeak = 600;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;


        settings->set_use_only_cores_in_stock(true);
        OpenMagnetics::clear_loaded_cores();

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights);


        CHECK(coreDatabase.size() < 3000);
        CHECK(masMagnetics.size() == 1);
        settings->reset();
    }

    TEST(Test_All_Cores_Load_Internally_All) {
        OpenMagnetics::clear_databases();
        double voltagePeakToPeak = 600;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;


        settings->set_use_only_cores_in_stock(false);
        OpenMagnetics::clear_loaded_cores();

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights);


        CHECK(coreDatabase.size() > 3000);
        CHECK(masMagnetics.size() == 1);
        settings->reset();
    }

    TEST(Test_All_Cores_Two_Chosen_Ones) {
        OpenMagnetics::clear_databases();
        double voltagePeakToPeak = 600;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;

        settings->set_use_toroidal_cores(true);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &cores, 2);

        CHECK(masMagnetics.size() == 2);
        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_name() == "T 18/9.0/7.1 - Kool Mµ Hƒ 40 - Ungapped");
        CHECK(masMagnetics[1].first.get_magnetic().get_core().get_name() == "EP 20 - 3C91 - Gapped 0.382 mm");
        settings->reset();
    }

    TEST(Test_No_Toroids_High_Power) {
        OpenMagnetics::clear_databases();
        double voltagePeakToPeak = 6000;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &cores);

        CHECK(masMagnetics.size() == 1);
        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_name() == "E 70/33/32 - 95 - Distributed gapped 1.36 mm");
        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_functional_description().get_number_stacks() == 1);
        settings->reset();
    }


    TEST(Test_Only_Toroids_High_Power) {
        OpenMagnetics::clear_databases();
        double voltagePeakToPeak = 600;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 1e-3;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        settings->set_use_toroidal_cores(true);
        settings->set_use_concentric_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, 10);

        // CHECK(masMagnetics.size() == 1);
        for (auto [mas, scoring] : masMagnetics) {
            CHECK(mas.get_mutable_magnetic().get_mutable_core().get_shape_family() == OpenMagnetics::CoreShapeFamily::T);
        }
        settings->reset();
    }


    TEST(Test_No_Toroids_High_Power_High_Frequency) {
        OpenMagnetics::clear_databases();
        double voltagePeakToPeak = 600000;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 500000;
        double desiredMagnetizingInductance = 10e-3;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &cores);

        CHECK(masMagnetics.size() == 1);

        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_name() == "E 70/33/32 - N87 - Gapped 4.0 mm");
        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_functional_description().get_number_stacks() == 3);
        settings->reset();
    }

    TEST(Test_No_Toroids_Low_Power) {
        OpenMagnetics::clear_databases();
        double voltagePeakToPeak = 60;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &cores);

        CHECK(masMagnetics.size() == 1);

        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_name() == "EFD 10/5/3 - 3C95 - Gapped 0.13999999999999999 mm");
        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_functional_description().get_number_stacks() == 1);
        settings->reset();
    }

    TEST(Test_No_Toroids_Low_Power_Low_Losses) {
        OpenMagnetics::clear_databases();
        double voltagePeakToPeak = 60;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &cores);

        CHECK(masMagnetics.size() == 1);

        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_name() == "E 14/3.5/5 - 3C97 - Gapped 0.27999999999999997 mm");
        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_functional_description().get_number_stacks() == 1);
        settings->reset();
    }

    TEST(Test_No_Toroids_Low_Power_Low_Losses_No_Care_About_Size) {
        OpenMagnetics::clear_databases();
        double voltagePeakToPeak = 60;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 0.1;

        OpenMagnetics::OperatingPoint operatingPoint;
        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &cores);

        CHECK(masMagnetics.size() == 1);

        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_name() == "RM 12 - 3C97 - Gapped 0.27 mm");
        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_functional_description().get_number_stacks() == 1);
        settings->reset();
    }

    TEST(Test_No_Toroids_Redo_Culling) {
        OpenMagnetics::clear_databases();
        double voltagePeakToPeak = 6000;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-3;
        std::vector<double> turnsRatios = {};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 0.1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &cores);

        CHECK(masMagnetics.size() == 1);

        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_name() == "E 25.4/10/7 - N27 - Gapped 0.5 mm");
        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_functional_description().get_number_stacks() == 1);
        settings->reset();
    }

    TEST(Test_No_Toroids_Two_Windings) {
        OpenMagnetics::clear_databases();
        double voltagePeakToPeak = 600;
        double dcCurrent = 0;
        double ambientTemperature = 25;
        double frequency = 100000;
        double desiredMagnetizingInductance = 10e-5;
        std::vector<double> turnsRatios = {0.1};
        OpenMagnetics::InputsWrapper inputs;

        prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &cores, 2);

        CHECK(masMagnetics.size() == 2);

        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_name() == "PQ 26/20 - 3C95 - Gapped 0.365 mm");

        CHECK(masMagnetics[1].first.get_magnetic().get_core().get_name() == "PQ 26/20 - 3C97 - Gapped 0.365 mm");
        settings->reset();
    }

    TEST(Test_No_Toroids_Two_Points_High_Power_Low_Power) {
        OpenMagnetics::clear_databases();
        OpenMagnetics::InputsWrapper inputs;
        std::vector<double> turnsRatios = {1};
        {
            double voltagePeakToPeak = 6000;
            double dcCurrent = 0;
            double ambientTemperature = 25;
            double frequency = 100000;
            double desiredMagnetizingInductance = 10e-5;

            prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);
        }
        auto highPowerOperatingPoint = inputs.get_operating_point(0);


        {
            double voltagePeakToPeak = 60;
            double dcCurrent = 0;
            double ambientTemperature = 25;
            double frequency = 100000;
            double desiredMagnetizingInductance = 10e-5;
            prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);
        }
        inputs.get_mutable_operating_points().push_back(highPowerOperatingPoint);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        OpenMagnetics::OperatingPoint operatingPoint;
        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &cores);

        CHECK(masMagnetics.size() == 1);


        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_name() == "E 100/60/28 - 3C95 - Gapped 2.0 mm");
        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_functional_description().get_number_stacks() == 1);

        CHECK(masMagnetics[0].first.get_outputs().size() == 2);
        CHECK(masMagnetics[0].first.get_magnetic().get_coil().get_functional_description().size() == 2);

        settings->reset();
    }

    TEST(Test_Two_Points_Equal) {
        OpenMagnetics::clear_databases();
        OpenMagnetics::InputsWrapper inputs;
        std::vector<double> turnsRatios = {};
        {
            double voltagePeakToPeak = 6000;
            double dcCurrent = 0;
            double ambientTemperature = 25;
            double frequency = 100000;
            double desiredMagnetizingInductance = 10e-5;

            prepare_test_parameters(dcCurrent, ambientTemperature, frequency, turnsRatios, desiredMagnetizingInductance, inputs, voltagePeakToPeak);
        }
        auto operatingPoint = inputs.get_operating_point(0);
        inputs.get_mutable_operating_points().push_back(operatingPoint);

        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::ENERGY_STORED] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::COST] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::EFFICIENCY] = 1;
        weights[OpenMagnetics::CoreAdviser::CoreAdviserFilters::DIMENSIONS] = 1;

        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto cores = load_test_data();
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &cores, 1);

        CHECK(masMagnetics.size() == 1);


        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_name() == "E 70/33/32 - 95 - Distributed gapped 1.36 mm");
        CHECK(masMagnetics[0].first.get_magnetic().get_core().get_functional_description().get_number_stacks() == 1);
        auto scorings = coreAdviser.get_scorings();

        settings->reset();
    }

    TEST(Test_CoreAdviser_Web_0) {
        OpenMagnetics::clear_databases();
        OpenMagnetics::InputsWrapper inputs = json::parse(R"({"designRequirements":{"insulation":null,"leakageInductance":null,"magnetizingInductance":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1e-05},"market":null,"maximumDimensions":null,"maximumWeight":null,"name":null,"operatingTemperature":null,"strayCapacitance":null,"terminalType":null,"topology":null,"turnsRatios":[]},"operatingPoints":[{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":42.0,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[2.6732088764802597e-15,3.8218284740811588,1.352034707786142,0.4253304394765631,1.0937048636780878e-15,0.1536120213245562,0.15119524913566243,0.07875264997066994,7.241115656602231e-16,0.047948399955734654,0.05513623287662578,0.03235759328251896,4.090564786249015e-16,0.023392967261331783,0.028681570420819397,0.017771008577657642,2.4672126293013086e-16,0.014016396066471176,0.017807158306087076,0.011386475028370351,1.714592334083309e-16,0.009474348163803117,0.012316262725232965,0.008041926341760328,1.1212044817922766e-16,0.006942375707716124,0.009173279000770123,0.006081161484994412,2.111598525624007e-17,0.00539520003414077,0.007217893995285857,0.004841101358177077,2.7775084398586677e-16,0.004388239086960148,0.005929248667932901,0.004014508967799147,9.486404471810101e-17,0.0037036469670079972,0.005045721135369522,0.0034434847131706546,7.352254058261595e-17,0.0032247860657257685,0.004424650744237962,0.0030404540187279,7.970041183113523e-17,0.002884979494813953,0.003983384372505874,0.0027540516952750064,7.285867231955118e-17,0.0026442785821523935,0.0036719561660329553,0.002552983887311817,1.6777299785974844e-16,0.002478058285030729,0.0034594522045109273,0.002417849611852019,9.246617731272175e-17,0.0023710817878942196,0.0033268345221317414,0.0023367953000158204,1.8921291035788944e-16,0.002314304311501306,0.0032630646073174008,0.002303167019442115],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":131904.75985510723,"average":2.609024107869118e-15,"dutyCycle":0.24999999999999975,"effectiveFrequency":131904.75985510723,"label":"Triangular","offset":-4.440892098500626e-16,"peak":4.999999999999999,"peakToPeak":10.0,"rms":2.8876908907560668,"thd":0.3765598442298381},"waveform":{"ancillaryLabel":null,"data":[-5.0,5.0,-5.0],"numberPeriods":null,"time":[0.0,2.4271844660194174e-06,9.70873786407767e-06]}},"frequency":103000.0,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":{"harmonics":{"amplitudes":[0.006798274035336369,0.6948836814766789,0.24493767435427322,0.07659321568695032,0.0011851486650484215,0.02713826890868042,0.026338258337802543,0.013524003490593392,0.001185148665048509,0.007934821893277648,0.00885178693142459,0.005121342057388876,0.0011851486650485484,0.003521567475518845,0.004036010998935581,0.0025392588578627244,0.001185148665048555,0.0019074529966875186,0.002056445686787696,0.0014922821772143813,0.0011851486650485722,0.00122003158014941,0.001056889452408924,0.0010459731468581636,0.0011851486650485724,0.0009399767435062919,0.0004847443569363152,0.0008801360125169896,0.0011851486650485668,0.0008503270906284321,0.00012878834865636887,0.0008391415445710964,0.0011851486650485696,0.0008389428243523104,0.0001057951454006002,0.0008448591470041428,0.001185148665048572,0.000853888382471694,0.000266631468075954,0.0008642251603076917,0.0011851486650485711,0.0008748056010587637,0.0003796903985003957,0.0008850169925296278,0.0011851486650485696,0.0008945181777242711,0.0004600180147137569,0.000903129917240823,0.0011851486650485716,0.000910768242017427,0.0005167100425910805,0.0009174039220630434,0.00118514866504857,0.0009230377659956728,0.0005553940161691154,0.000927685539989642,0.0011851486650485783,0.0009313687576957952,0.0005795355845633968,0.0009341090721450906,0.001185148665048575,0.0009359248919968734,0.0005911441855611932,0.0009368293872760598],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":133347.9143873257,"average":-0.006798274035336401,"dutyCycle":0.25,"effectiveFrequency":133342.3181567079,"label":"Custom","offset":0.0020903409525253602,"peak":0.9028033263894203,"peakToPeak":1.8014259708737899,"rms":0.5247313861450679,"thd":0.3745218047109},"waveform":{"ancillaryLabel":null,"data":[-0.8417355085620355,-0.8606978872028123,0.9028033263894203,-0.8986226444843696,-0.8417355085620355],"numberPeriods":null,"time":[0.0,7.58495145631068e-08,2.427184466019418e-06,9.632888349514574e-06,9.708737864077681e-06]}},"name":"Primary winding excitation","voltage":{"harmonics":{"amplitudes":[0.078125,4.390192203074813,3.1805418164042503,1.6082999269932894,0.15625,0.7853068032814723,1.0533519383460939,0.747228683429235,0.15625,0.38152854838936,0.6237849662140766,0.5097284946091541,0.15625,0.22395874930175205,0.4366894957016387,0.3969081532881607,0.15625,0.1387657119725448,0.3303628683669758,0.32998229455569084,0.15625,0.08454863298834173,0.2606873758724238,0.2849610019773425,0.15625,0.046391828940503435,0.21067873648230084,0.25205943526842706,0.15625,0.017597888807606342,0.1723953087083574,0.22653191551027138,0.15625,0.005294556301253328,0.14161674515924202,0.20579082584016142,0.15625,0.02426175514926864,0.11588289785500591,0.18829805684603593,0.15625,0.040521304735682134,0.09365264588780098,0.17307466014248524,0.15625,0.05487174035875039,0.07390074623301954,0.1594602594300079,0.15625,0.06786649949748992,0.05590714395539467,0.1469848919856514,0.15625,0.07990999700987583,0.03913858752989208,0.13529576856110245,0.15625,0.09131439305868466,0.023177498052867507,0.1241125131778793,0.15625,0.10233553928853857,0.007676070276481761,0.11319769974196572],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":711072.940860665,"average":-0.078125,"dutyCycle":0.24218749999999975,"effectiveFrequency":711013.8104343135,"label":"Rectangular","offset":0.0,"peak":7.5,"peakToPeak":10.0,"rms":4.284784125250653,"thd":0.9507084995254541},"waveform":{"ancillaryLabel":null,"data":[-2.5,7.5,7.5,-2.5,-2.5],"numberPeriods":null,"time":[0.0,0.0,2.4271844660194174e-06,2.4271844660194174e-06,9.70873786407767e-06]}}}],"name":"Operating Point No. 2"}]})");
        auto cores = load_test_data();
        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto masMagnetics = coreAdviser.get_advised_core(inputs, &cores, 1);

        settings->reset();
    }

    TEST(Test_CoreAdviser_Web_1) {
        OpenMagnetics::clear_databases();
        std::string inputsString = R"({"designRequirements":{"insulation":null,"leakageInductance":null,"magnetizingInductance":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1e-05},"market":null,"maximumDimensions":null,"maximumWeight":null,"name":null,"operatingTemperature":null,"strayCapacitance":null,"terminalType":null,"topology":null,"turnsRatios":[]},"operatingPoints":[{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":42.0,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[2.6732088764802597e-15,3.8218284740811588,1.352034707786142,0.4253304394765631,1.0937048636780878e-15,0.1536120213245562,0.15119524913566243,0.07875264997066994,7.241115656602231e-16,0.047948399955734654,0.05513623287662578,0.03235759328251896,4.090564786249015e-16,0.023392967261331783,0.028681570420819397,0.017771008577657642,2.4672126293013086e-16,0.014016396066471176,0.017807158306087076,0.011386475028370351,1.714592334083309e-16,0.009474348163803117,0.012316262725232965,0.008041926341760328,1.1212044817922766e-16,0.006942375707716124,0.009173279000770123,0.006081161484994412,2.111598525624007e-17,0.00539520003414077,0.007217893995285857,0.004841101358177077,2.7775084398586677e-16,0.004388239086960148,0.005929248667932901,0.004014508967799147,9.486404471810101e-17,0.0037036469670079972,0.005045721135369522,0.0034434847131706546,7.352254058261595e-17,0.0032247860657257685,0.004424650744237962,0.0030404540187279,7.970041183113523e-17,0.002884979494813953,0.003983384372505874,0.0027540516952750064,7.285867231955118e-17,0.0026442785821523935,0.0036719561660329553,0.002552983887311817,1.6777299785974844e-16,0.002478058285030729,0.0034594522045109273,0.002417849611852019,9.246617731272175e-17,0.0023710817878942196,0.0033268345221317414,0.0023367953000158204,1.8921291035788944e-16,0.002314304311501306,0.0032630646073174008,0.002303167019442115],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":131904.75985510723,"average":2.609024107869118e-15,"dutyCycle":0.24999999999999975,"effectiveFrequency":131904.75985510723,"label":"Triangular","offset":-4.440892098500626e-16,"peak":4.999999999999999,"peakToPeak":10.0,"rms":2.8876908907560668,"thd":0.3765598442298381},"waveform":{"ancillaryLabel":null,"data":[-5.0,5.0,-5.0],"numberPeriods":null,"time":[0.0,2.4271844660194174e-06,9.70873786407767e-06]}},"frequency":103000.0,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":{"harmonics":{"amplitudes":[0.006798274035336369,0.6948836814766789,0.24493767435427322,0.07659321568695032,0.0011851486650484215,0.02713826890868042,0.026338258337802543,0.013524003490593392,0.001185148665048509,0.007934821893277648,0.00885178693142459,0.005121342057388876,0.0011851486650485484,0.003521567475518845,0.004036010998935581,0.0025392588578627244,0.001185148665048555,0.0019074529966875186,0.002056445686787696,0.0014922821772143813,0.0011851486650485722,0.00122003158014941,0.001056889452408924,0.0010459731468581636,0.0011851486650485724,0.0009399767435062919,0.0004847443569363152,0.0008801360125169896,0.0011851486650485668,0.0008503270906284321,0.00012878834865636887,0.0008391415445710964,0.0011851486650485696,0.0008389428243523104,0.0001057951454006002,0.0008448591470041428,0.001185148665048572,0.000853888382471694,0.000266631468075954,0.0008642251603076917,0.0011851486650485711,0.0008748056010587637,0.0003796903985003957,0.0008850169925296278,0.0011851486650485696,0.0008945181777242711,0.0004600180147137569,0.000903129917240823,0.0011851486650485716,0.000910768242017427,0.0005167100425910805,0.0009174039220630434,0.00118514866504857,0.0009230377659956728,0.0005553940161691154,0.000927685539989642,0.0011851486650485783,0.0009313687576957952,0.0005795355845633968,0.0009341090721450906,0.001185148665048575,0.0009359248919968734,0.0005911441855611932,0.0009368293872760598],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":133347.9143873257,"average":-0.006798274035336401,"dutyCycle":0.25,"effectiveFrequency":133342.3181567079,"label":"Custom","offset":0.0020903409525253602,"peak":0.9028033263894203,"peakToPeak":1.8014259708737899,"rms":0.5247313861450679,"thd":0.3745218047109},"waveform":{"ancillaryLabel":null,"data":[-0.8417355085620355,-0.8606978872028123,0.9028033263894203,-0.8986226444843696,-0.8417355085620355],"numberPeriods":null,"time":[0.0,7.58495145631068e-08,2.427184466019418e-06,9.632888349514574e-06,9.708737864077681e-06]}},"name":"Primary winding excitation","voltage":{"harmonics":{"amplitudes":[0.078125,4.390192203074813,3.1805418164042503,1.6082999269932894,0.15625,0.7853068032814723,1.0533519383460939,0.747228683429235,0.15625,0.38152854838936,0.6237849662140766,0.5097284946091541,0.15625,0.22395874930175205,0.4366894957016387,0.3969081532881607,0.15625,0.1387657119725448,0.3303628683669758,0.32998229455569084,0.15625,0.08454863298834173,0.2606873758724238,0.2849610019773425,0.15625,0.046391828940503435,0.21067873648230084,0.25205943526842706,0.15625,0.017597888807606342,0.1723953087083574,0.22653191551027138,0.15625,0.005294556301253328,0.14161674515924202,0.20579082584016142,0.15625,0.02426175514926864,0.11588289785500591,0.18829805684603593,0.15625,0.040521304735682134,0.09365264588780098,0.17307466014248524,0.15625,0.05487174035875039,0.07390074623301954,0.1594602594300079,0.15625,0.06786649949748992,0.05590714395539467,0.1469848919856514,0.15625,0.07990999700987583,0.03913858752989208,0.13529576856110245,0.15625,0.09131439305868466,0.023177498052867507,0.1241125131778793,0.15625,0.10233553928853857,0.007676070276481761,0.11319769974196572],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":711072.940860665,"average":-0.078125,"dutyCycle":0.24218749999999975,"effectiveFrequency":711013.8104343135,"label":"Rectangular","offset":0.0,"peak":7.5,"peakToPeak":10.0,"rms":4.284784125250653,"thd":0.9507084995254541},"waveform":{"ancillaryLabel":null,"data":[-2.5,7.5,7.5,-2.5,-2.5],"numberPeriods":null,"time":[0.0,0.0,2.4271844660194174e-06,2.4271844660194174e-06,9.70873786407767e-06]}}}],"name":"Operating Point No. 2"}]})";
        std::string file_path = __FILE__;
        auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/testData/cores_stock.ndjson");

        std::ifstream ndjsonFile(inventory_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            coresJson.push_back(jf);
            OpenMagnetics::CoreWrapper core(jf);
        }

        std::string coresString = coresJson.dump(4);

        OpenMagnetics::InputsWrapper inputs(json::parse(inputsString));
        std::vector<OpenMagnetics::CoreWrapper> cores = json::parse(coresString);

        std::vector<OpenMagnetics::CoreWrapper> coresInInventory;
        for (auto& core : cores) {
            if (core.get_distributors_info()) {
                if (core.get_distributors_info().value().size() > 0) {
                    coresInInventory.push_back(core);
                }
            }
        }

        OpenMagnetics::CoreAdviser coreAdviser;
        auto masMagnetics = coreAdviser.get_advised_core(inputs, &coresInInventory, 1000);
        auto log = coreAdviser.read_log();
        auto scores = coreAdviser.get_scorings();
        json result;
        json resultScores = scores;

        json resultCores = json::array();
        for (auto& masMagnetic : masMagnetics) {
            json aux;
            to_json(aux, masMagnetic.first.get_magnetic());
            to_json(aux, masMagnetic.first.get_outputs()[0].get_core_losses().value());
            to_json(aux, masMagnetic.first.get_outputs()[0].get_winding_losses().value());
        if (masMagnetic.first.get_outputs()[0].get_magnetizing_inductance()) {
            to_json(aux, masMagnetic.first.get_outputs()[0].get_magnetizing_inductance().value());
        }
            resultCores.push_back(aux);
        }
        result["cores"] = resultCores;
        result["log"] = log;
        result["scores"] = scores;

        result.dump(4);

        settings->reset();
    }

    TEST(Test_CoreAdviser_Web_2) {
        OpenMagnetics::clear_databases();
        std::string inputsString = R"({"designRequirements":{"insulation":null,"leakageInductance":null,"magnetizingInductance":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1e-05},"market":null,"maximumDimensions":null,"maximumWeight":null,"name":null,"operatingTemperature":null,"strayCapacitance":null,"terminalType":null,"topology":null,"turnsRatios":[]},"operatingPoints":[{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":42.0,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[2.6732088764802597e-15,3.8218284740811588,1.352034707786142,0.4253304394765631,1.0937048636780878e-15,0.1536120213245562,0.15119524913566243,0.07875264997066994,7.241115656602231e-16,0.047948399955734654,0.05513623287662578,0.03235759328251896,4.090564786249015e-16,0.023392967261331783,0.028681570420819397,0.017771008577657642,2.4672126293013086e-16,0.014016396066471176,0.017807158306087076,0.011386475028370351,1.714592334083309e-16,0.009474348163803117,0.012316262725232965,0.008041926341760328,1.1212044817922766e-16,0.006942375707716124,0.009173279000770123,0.006081161484994412,2.111598525624007e-17,0.00539520003414077,0.007217893995285857,0.004841101358177077,2.7775084398586677e-16,0.004388239086960148,0.005929248667932901,0.004014508967799147,9.486404471810101e-17,0.0037036469670079972,0.005045721135369522,0.0034434847131706546,7.352254058261595e-17,0.0032247860657257685,0.004424650744237962,0.0030404540187279,7.970041183113523e-17,0.002884979494813953,0.003983384372505874,0.0027540516952750064,7.285867231955118e-17,0.0026442785821523935,0.0036719561660329553,0.002552983887311817,1.6777299785974844e-16,0.002478058285030729,0.0034594522045109273,0.002417849611852019,9.246617731272175e-17,0.0023710817878942196,0.0033268345221317414,0.0023367953000158204,1.8921291035788944e-16,0.002314304311501306,0.0032630646073174008,0.002303167019442115],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":131904.75985510723,"average":2.609024107869118e-15,"dutyCycle":0.24999999999999975,"effectiveFrequency":131904.75985510723,"label":"Triangular","offset":-4.440892098500626e-16,"peak":4.999999999999999,"peakToPeak":10.0,"rms":2.8876908907560668,"thd":0.3765598442298381},"waveform":{"ancillaryLabel":null,"data":[-5.0,5.0,-5.0],"numberPeriods":null,"time":[0.0,2.4271844660194174e-06,9.70873786407767e-06]}},"frequency":103000.0,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":{"harmonics":{"amplitudes":[0.006798274035336369,0.6948836814766789,0.24493767435427322,0.07659321568695032,0.0011851486650484215,0.02713826890868042,0.026338258337802543,0.013524003490593392,0.001185148665048509,0.007934821893277648,0.00885178693142459,0.005121342057388876,0.0011851486650485484,0.003521567475518845,0.004036010998935581,0.0025392588578627244,0.001185148665048555,0.0019074529966875186,0.002056445686787696,0.0014922821772143813,0.0011851486650485722,0.00122003158014941,0.001056889452408924,0.0010459731468581636,0.0011851486650485724,0.0009399767435062919,0.0004847443569363152,0.0008801360125169896,0.0011851486650485668,0.0008503270906284321,0.00012878834865636887,0.0008391415445710964,0.0011851486650485696,0.0008389428243523104,0.0001057951454006002,0.0008448591470041428,0.001185148665048572,0.000853888382471694,0.000266631468075954,0.0008642251603076917,0.0011851486650485711,0.0008748056010587637,0.0003796903985003957,0.0008850169925296278,0.0011851486650485696,0.0008945181777242711,0.0004600180147137569,0.000903129917240823,0.0011851486650485716,0.000910768242017427,0.0005167100425910805,0.0009174039220630434,0.00118514866504857,0.0009230377659956728,0.0005553940161691154,0.000927685539989642,0.0011851486650485783,0.0009313687576957952,0.0005795355845633968,0.0009341090721450906,0.001185148665048575,0.0009359248919968734,0.0005911441855611932,0.0009368293872760598],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":133347.9143873257,"average":-0.006798274035336401,"dutyCycle":0.25,"effectiveFrequency":133342.3181567079,"label":"Custom","offset":0.0020903409525253602,"peak":0.9028033263894203,"peakToPeak":1.8014259708737899,"rms":0.5247313861450679,"thd":0.3745218047109},"waveform":{"ancillaryLabel":null,"data":[-0.8417355085620355,-0.8606978872028123,0.9028033263894203,-0.8986226444843696,-0.8417355085620355],"numberPeriods":null,"time":[0.0,7.58495145631068e-08,2.427184466019418e-06,9.632888349514574e-06,9.708737864077681e-06]}},"name":"Primary winding excitation","voltage":{"harmonics":{"amplitudes":[0.078125,4.390192203074813,3.1805418164042503,1.6082999269932894,0.15625,0.7853068032814723,1.0533519383460939,0.747228683429235,0.15625,0.38152854838936,0.6237849662140766,0.5097284946091541,0.15625,0.22395874930175205,0.4366894957016387,0.3969081532881607,0.15625,0.1387657119725448,0.3303628683669758,0.32998229455569084,0.15625,0.08454863298834173,0.2606873758724238,0.2849610019773425,0.15625,0.046391828940503435,0.21067873648230084,0.25205943526842706,0.15625,0.017597888807606342,0.1723953087083574,0.22653191551027138,0.15625,0.005294556301253328,0.14161674515924202,0.20579082584016142,0.15625,0.02426175514926864,0.11588289785500591,0.18829805684603593,0.15625,0.040521304735682134,0.09365264588780098,0.17307466014248524,0.15625,0.05487174035875039,0.07390074623301954,0.1594602594300079,0.15625,0.06786649949748992,0.05590714395539467,0.1469848919856514,0.15625,0.07990999700987583,0.03913858752989208,0.13529576856110245,0.15625,0.09131439305868466,0.023177498052867507,0.1241125131778793,0.15625,0.10233553928853857,0.007676070276481761,0.11319769974196572],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":711072.940860665,"average":-0.078125,"dutyCycle":0.24218749999999975,"effectiveFrequency":711013.8104343135,"label":"Rectangular","offset":0.0,"peak":7.5,"peakToPeak":10.0,"rms":4.284784125250653,"thd":0.9507084995254541},"waveform":{"ancillaryLabel":null,"data":[-2.5,7.5,7.5,-2.5,-2.5],"numberPeriods":null,"time":[0.0,0.0,2.4271844660194174e-06,2.4271844660194174e-06,9.70873786407767e-06]}}}],"name":"Operating Point No. 2"}]})";
        std::string file_path = __FILE__;
        auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/testData/cores_stock.ndjson");

        std::ifstream ndjsonFile(inventory_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            coresJson.push_back(jf);
        }

        std::string coresString = coresJson.dump(4);

        OpenMagnetics::InputsWrapper inputs(json::parse(inputsString));
        std::vector<OpenMagnetics::CoreWrapper> cores = json::parse(coresString);

        std::vector<OpenMagnetics::CoreWrapper> coresInInventory;
        for (auto& core : cores) {
            if (core.get_distributors_info()) {
                if (core.get_distributors_info().value().size() > 0) {
                    coresInInventory.push_back(core);
                }
            }
        }
        OpenMagnetics::CoreAdviser coreAdviser;
        auto masMagnetics = coreAdviser.get_advised_core(inputs, &coresInInventory, 100, 500);
        auto log = coreAdviser.read_log();
        auto scores = coreAdviser.get_scorings();
        json result;
        json resultScores = scores;

        json resultCores = json::array();
        for (auto& masMagnetic : masMagnetics) {
            json aux;
            to_json(aux, masMagnetic.first.get_magnetic());
            to_json(aux, masMagnetic.first.get_outputs()[0].get_core_losses().value());
            to_json(aux, masMagnetic.first.get_outputs()[0].get_winding_losses().value());
            if(masMagnetic.first.get_outputs()[0].get_magnetizing_inductance()) {
                to_json(aux, masMagnetic.first.get_outputs()[0].get_magnetizing_inductance().value());
            }
            resultCores.push_back(aux);
        }
        result["cores"] = resultCores;
        result["log"] = log;
        result["scores"] = scores;

        result.dump(4);

        settings->reset();
    }

    TEST(Test_CoreAdviser_Web_3) {
        OpenMagnetics::clear_databases();
        std::string inputsString = R"({"designRequirements":{"insulation":null,"leakageInductance":null,"magnetizingInductance":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1e-05},"market":null,"maximumDimensions":null,"maximumWeight":null,"name":null,"operatingTemperature":null,"strayCapacitance":null,"terminalType":null,"topology":null,"turnsRatios":[]},"operatingPoints":[{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":42.0,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[2.6732088764802597e-15,3.8218284740811588,1.352034707786142,0.4253304394765631,1.0937048636780878e-15,0.1536120213245562,0.15119524913566243,0.07875264997066994,7.241115656602231e-16,0.047948399955734654,0.05513623287662578,0.03235759328251896,4.090564786249015e-16,0.023392967261331783,0.028681570420819397,0.017771008577657642,2.4672126293013086e-16,0.014016396066471176,0.017807158306087076,0.011386475028370351,1.714592334083309e-16,0.009474348163803117,0.012316262725232965,0.008041926341760328,1.1212044817922766e-16,0.006942375707716124,0.009173279000770123,0.006081161484994412,2.111598525624007e-17,0.00539520003414077,0.007217893995285857,0.004841101358177077,2.7775084398586677e-16,0.004388239086960148,0.005929248667932901,0.004014508967799147,9.486404471810101e-17,0.0037036469670079972,0.005045721135369522,0.0034434847131706546,7.352254058261595e-17,0.0032247860657257685,0.004424650744237962,0.0030404540187279,7.970041183113523e-17,0.002884979494813953,0.003983384372505874,0.0027540516952750064,7.285867231955118e-17,0.0026442785821523935,0.0036719561660329553,0.002552983887311817,1.6777299785974844e-16,0.002478058285030729,0.0034594522045109273,0.002417849611852019,9.246617731272175e-17,0.0023710817878942196,0.0033268345221317414,0.0023367953000158204,1.8921291035788944e-16,0.002314304311501306,0.0032630646073174008,0.002303167019442115],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":131904.75985510723,"average":2.609024107869118e-15,"dutyCycle":0.24999999999999975,"effectiveFrequency":131904.75985510723,"label":"Triangular","offset":-4.440892098500626e-16,"peak":4.999999999999999,"peakToPeak":10.0,"rms":2.8876908907560668,"thd":0.3765598442298381},"waveform":{"ancillaryLabel":null,"data":[-5.0,5.0,-5.0],"numberPeriods":null,"time":[0.0,2.4271844660194174e-06,9.70873786407767e-06]}},"frequency":103000.0,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":{"harmonics":{"amplitudes":[0.006798274035336369,0.6948836814766789,0.24493767435427322,0.07659321568695032,0.0011851486650484215,0.02713826890868042,0.026338258337802543,0.013524003490593392,0.001185148665048509,0.007934821893277648,0.00885178693142459,0.005121342057388876,0.0011851486650485484,0.003521567475518845,0.004036010998935581,0.0025392588578627244,0.001185148665048555,0.0019074529966875186,0.002056445686787696,0.0014922821772143813,0.0011851486650485722,0.00122003158014941,0.001056889452408924,0.0010459731468581636,0.0011851486650485724,0.0009399767435062919,0.0004847443569363152,0.0008801360125169896,0.0011851486650485668,0.0008503270906284321,0.00012878834865636887,0.0008391415445710964,0.0011851486650485696,0.0008389428243523104,0.0001057951454006002,0.0008448591470041428,0.001185148665048572,0.000853888382471694,0.000266631468075954,0.0008642251603076917,0.0011851486650485711,0.0008748056010587637,0.0003796903985003957,0.0008850169925296278,0.0011851486650485696,0.0008945181777242711,0.0004600180147137569,0.000903129917240823,0.0011851486650485716,0.000910768242017427,0.0005167100425910805,0.0009174039220630434,0.00118514866504857,0.0009230377659956728,0.0005553940161691154,0.000927685539989642,0.0011851486650485783,0.0009313687576957952,0.0005795355845633968,0.0009341090721450906,0.001185148665048575,0.0009359248919968734,0.0005911441855611932,0.0009368293872760598],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":133347.9143873257,"average":-0.006798274035336401,"dutyCycle":0.25,"effectiveFrequency":133342.3181567079,"label":"Custom","offset":0.0020903409525253602,"peak":0.9028033263894203,"peakToPeak":1.8014259708737899,"rms":0.5247313861450679,"thd":0.3745218047109},"waveform":{"ancillaryLabel":null,"data":[-0.8417355085620355,-0.8606978872028123,0.9028033263894203,-0.8986226444843696,-0.8417355085620355],"numberPeriods":null,"time":[0.0,7.58495145631068e-08,2.427184466019418e-06,9.632888349514574e-06,9.708737864077681e-06]}},"name":"Primary winding excitation","voltage":{"harmonics":{"amplitudes":[0.078125,4.390192203074813,3.1805418164042503,1.6082999269932894,0.15625,0.7853068032814723,1.0533519383460939,0.747228683429235,0.15625,0.38152854838936,0.6237849662140766,0.5097284946091541,0.15625,0.22395874930175205,0.4366894957016387,0.3969081532881607,0.15625,0.1387657119725448,0.3303628683669758,0.32998229455569084,0.15625,0.08454863298834173,0.2606873758724238,0.2849610019773425,0.15625,0.046391828940503435,0.21067873648230084,0.25205943526842706,0.15625,0.017597888807606342,0.1723953087083574,0.22653191551027138,0.15625,0.005294556301253328,0.14161674515924202,0.20579082584016142,0.15625,0.02426175514926864,0.11588289785500591,0.18829805684603593,0.15625,0.040521304735682134,0.09365264588780098,0.17307466014248524,0.15625,0.05487174035875039,0.07390074623301954,0.1594602594300079,0.15625,0.06786649949748992,0.05590714395539467,0.1469848919856514,0.15625,0.07990999700987583,0.03913858752989208,0.13529576856110245,0.15625,0.09131439305868466,0.023177498052867507,0.1241125131778793,0.15625,0.10233553928853857,0.007676070276481761,0.11319769974196572],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":711072.940860665,"average":-0.078125,"dutyCycle":0.24218749999999975,"effectiveFrequency":711013.8104343135,"label":"Rectangular","offset":0.0,"peak":7.5,"peakToPeak":10.0,"rms":4.284784125250653,"thd":0.9507084995254541},"waveform":{"ancillaryLabel":null,"data":[-2.5,7.5,7.5,-2.5,-2.5],"numberPeriods":null,"time":[0.0,0.0,2.4271844660194174e-06,2.4271844660194174e-06,9.70873786407767e-06]}}}],"name":"Operating Point No. 2"}]})";
        std::string weightsString = R"({"AREA_PRODUCT":1,"ENERGY_STORED":0,"COST":0,"EFFICIENCY":0,"DIMENSIONS":0})";
        std::string file_path = __FILE__;

        std::map<std::string, double> weightsKeysString = json::parse(weightsString);
        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        for (auto const& pair : weightsKeysString) {
            weights[magic_enum::enum_cast<OpenMagnetics::CoreAdviser::CoreAdviserFilters>(pair.first).value()] = pair.second;
        }

        auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/testData/cores_stock.ndjson");

        std::ifstream ndjsonFile(inventory_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            coresJson.push_back(jf);
        }

        std::string coresString = coresJson.dump(4);

        OpenMagnetics::InputsWrapper inputs(json::parse(inputsString));
        std::vector<OpenMagnetics::CoreWrapper> cores = json::parse(coresString);

        std::vector<OpenMagnetics::CoreWrapper> coresInInventory;
        for (auto& core : cores) {
            if (core.get_distributors_info()) {
                if (core.get_distributors_info().value().size() > 0) {
                    coresInInventory.push_back(core);
                }
            }
        }
        OpenMagnetics::CoreAdviser coreAdviser;
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &coresInInventory, 20);
        auto log = coreAdviser.read_log();
        std::vector<std::string> listOfNames;
        for (size_t i = 0; i < masMagnetics.size(); ++i ){
            listOfNames.push_back(masMagnetics[i].first.get_magnetic().get_manufacturer_info().value().get_reference().value());
        }

        auto scorings = coreAdviser.get_scorings();

        CHECK(masMagnetics[0].first.get_magnetic().get_manufacturer_info().value().get_reference().value() == "EQ 13 - 3C95 - Ungapped");
        CHECK(scorings[masMagnetics[0].first.get_magnetic().get_manufacturer_info().value().get_reference().value()][OpenMagnetics::CoreAdviser::CoreAdviserFilters::AREA_PRODUCT] == 1);

        settings->reset();
    }

    TEST(Test_CoreAdviser_Web_4) {
        OpenMagnetics::clear_databases();
        std::string inputsString = R"({"designRequirements":{"insulation":null,"leakageInductance":null,"magnetizingInductance":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1e-05},"market":null,"maximumDimensions":null,"maximumWeight":null,"name":null,"operatingTemperature":null,"strayCapacitance":null,"terminalType":null,"topology":null,"turnsRatios":[]},"operatingPoints":[{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":42.0,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[2.6732088764802597e-15,3.8218284740811588,1.352034707786142,0.4253304394765631,1.0937048636780878e-15,0.1536120213245562,0.15119524913566243,0.07875264997066994,7.241115656602231e-16,0.047948399955734654,0.05513623287662578,0.03235759328251896,4.090564786249015e-16,0.023392967261331783,0.028681570420819397,0.017771008577657642,2.4672126293013086e-16,0.014016396066471176,0.017807158306087076,0.011386475028370351,1.714592334083309e-16,0.009474348163803117,0.012316262725232965,0.008041926341760328,1.1212044817922766e-16,0.006942375707716124,0.009173279000770123,0.006081161484994412,2.111598525624007e-17,0.00539520003414077,0.007217893995285857,0.004841101358177077,2.7775084398586677e-16,0.004388239086960148,0.005929248667932901,0.004014508967799147,9.486404471810101e-17,0.0037036469670079972,0.005045721135369522,0.0034434847131706546,7.352254058261595e-17,0.0032247860657257685,0.004424650744237962,0.0030404540187279,7.970041183113523e-17,0.002884979494813953,0.003983384372505874,0.0027540516952750064,7.285867231955118e-17,0.0026442785821523935,0.0036719561660329553,0.002552983887311817,1.6777299785974844e-16,0.002478058285030729,0.0034594522045109273,0.002417849611852019,9.246617731272175e-17,0.0023710817878942196,0.0033268345221317414,0.0023367953000158204,1.8921291035788944e-16,0.002314304311501306,0.0032630646073174008,0.002303167019442115],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":131904.75985510723,"average":2.609024107869118e-15,"dutyCycle":0.24999999999999975,"effectiveFrequency":131904.75985510723,"label":"Triangular","offset":-4.440892098500626e-16,"peak":4.999999999999999,"peakToPeak":10.0,"rms":2.8876908907560668,"thd":0.3765598442298381},"waveform":{"ancillaryLabel":null,"data":[-5.0,5.0,-5.0],"numberPeriods":null,"time":[0.0,2.4271844660194174e-06,9.70873786407767e-06]}},"frequency":103000.0,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":{"harmonics":{"amplitudes":[0.006798274035336369,0.6948836814766789,0.24493767435427322,0.07659321568695032,0.0011851486650484215,0.02713826890868042,0.026338258337802543,0.013524003490593392,0.001185148665048509,0.007934821893277648,0.00885178693142459,0.005121342057388876,0.0011851486650485484,0.003521567475518845,0.004036010998935581,0.0025392588578627244,0.001185148665048555,0.0019074529966875186,0.002056445686787696,0.0014922821772143813,0.0011851486650485722,0.00122003158014941,0.001056889452408924,0.0010459731468581636,0.0011851486650485724,0.0009399767435062919,0.0004847443569363152,0.0008801360125169896,0.0011851486650485668,0.0008503270906284321,0.00012878834865636887,0.0008391415445710964,0.0011851486650485696,0.0008389428243523104,0.0001057951454006002,0.0008448591470041428,0.001185148665048572,0.000853888382471694,0.000266631468075954,0.0008642251603076917,0.0011851486650485711,0.0008748056010587637,0.0003796903985003957,0.0008850169925296278,0.0011851486650485696,0.0008945181777242711,0.0004600180147137569,0.000903129917240823,0.0011851486650485716,0.000910768242017427,0.0005167100425910805,0.0009174039220630434,0.00118514866504857,0.0009230377659956728,0.0005553940161691154,0.000927685539989642,0.0011851486650485783,0.0009313687576957952,0.0005795355845633968,0.0009341090721450906,0.001185148665048575,0.0009359248919968734,0.0005911441855611932,0.0009368293872760598],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":133347.9143873257,"average":-0.006798274035336401,"dutyCycle":0.25,"effectiveFrequency":133342.3181567079,"label":"Custom","offset":0.0020903409525253602,"peak":0.9028033263894203,"peakToPeak":1.8014259708737899,"rms":0.5247313861450679,"thd":0.3745218047109},"waveform":{"ancillaryLabel":null,"data":[-0.8417355085620355,-0.8606978872028123,0.9028033263894203,-0.8986226444843696,-0.8417355085620355],"numberPeriods":null,"time":[0.0,7.58495145631068e-08,2.427184466019418e-06,9.632888349514574e-06,9.708737864077681e-06]}},"name":"Primary winding excitation","voltage":{"harmonics":{"amplitudes":[0.078125,4.390192203074813,3.1805418164042503,1.6082999269932894,0.15625,0.7853068032814723,1.0533519383460939,0.747228683429235,0.15625,0.38152854838936,0.6237849662140766,0.5097284946091541,0.15625,0.22395874930175205,0.4366894957016387,0.3969081532881607,0.15625,0.1387657119725448,0.3303628683669758,0.32998229455569084,0.15625,0.08454863298834173,0.2606873758724238,0.2849610019773425,0.15625,0.046391828940503435,0.21067873648230084,0.25205943526842706,0.15625,0.017597888807606342,0.1723953087083574,0.22653191551027138,0.15625,0.005294556301253328,0.14161674515924202,0.20579082584016142,0.15625,0.02426175514926864,0.11588289785500591,0.18829805684603593,0.15625,0.040521304735682134,0.09365264588780098,0.17307466014248524,0.15625,0.05487174035875039,0.07390074623301954,0.1594602594300079,0.15625,0.06786649949748992,0.05590714395539467,0.1469848919856514,0.15625,0.07990999700987583,0.03913858752989208,0.13529576856110245,0.15625,0.09131439305868466,0.023177498052867507,0.1241125131778793,0.15625,0.10233553928853857,0.007676070276481761,0.11319769974196572],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":711072.940860665,"average":-0.078125,"dutyCycle":0.24218749999999975,"effectiveFrequency":711013.8104343135,"label":"Rectangular","offset":0.0,"peak":7.5,"peakToPeak":10.0,"rms":4.284784125250653,"thd":0.9507084995254541},"waveform":{"ancillaryLabel":null,"data":[-2.5,7.5,7.5,-2.5,-2.5],"numberPeriods":null,"time":[0.0,0.0,2.4271844660194174e-06,2.4271844660194174e-06,9.70873786407767e-06]}}}],"name":"Operating Point No. 2"}]})";
        std::string weightsString = R"({"AREA_PRODUCT":0.1,"ENERGY_STORED":0.1,"COST":0.1,"EFFICIENCY":1,"DIMENSIONS":0.1})";
        std::string file_path = __FILE__;

        std::map<std::string, double> weightsKeysString = json::parse(weightsString);
        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        for (auto const& pair : weightsKeysString) {
            weights[magic_enum::enum_cast<OpenMagnetics::CoreAdviser::CoreAdviserFilters>(pair.first).value()] = pair.second;
        }

        auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/testData/cores_stock.ndjson");

        std::ifstream ndjsonFile(inventory_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            coresJson.push_back(jf);
        }

        std::string coresString = coresJson.dump(4);

        OpenMagnetics::InputsWrapper inputs(json::parse(inputsString));
        std::vector<OpenMagnetics::CoreWrapper> cores = json::parse(coresString);

        std::vector<OpenMagnetics::CoreWrapper> coresInInventory;
        for (auto& core : cores) {
            if (core.get_distributors_info()) {
                if (core.get_distributors_info().value().size() > 0) {
                    coresInInventory.push_back(core);
                }
            }
        }
        OpenMagnetics::CoreAdviser coreAdviser;
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &coresInInventory, 500);
        auto log = coreAdviser.read_log();
        std::vector<std::string> listOfNames;
        for (size_t i = 0; i < masMagnetics.size(); ++i ){
            listOfNames.push_back(masMagnetics[i].first.get_magnetic().get_manufacturer_info().value().get_reference().value());
        }

        auto scorings = coreAdviser.get_scorings();

        settings->reset();
    }

    TEST(Test_CoreAdviser_Web_5) {
        OpenMagnetics::clear_databases();
        std::string inputsString = R"({"designRequirements":{"insulation":null,"leakageInductance":null,"magnetizingInductance":{"excludeMaximum":null,"excludeMinimum":null,"maximum":null,"minimum":null,"nominal":1e-05},"market":null,"maximumDimensions":null,"maximumWeight":null,"name":null,"operatingTemperature":null,"strayCapacitance":null,"terminalType":null,"topology":null,"turnsRatios":[]},"operatingPoints":[{"conditions":{"ambientRelativeHumidity":null,"ambientTemperature":42.0,"cooling":null,"name":null},"excitationsPerWinding":[{"current":{"harmonics":{"amplitudes":[2.6732088764802597e-15,3.8218284740811588,1.352034707786142,0.4253304394765631,1.0937048636780878e-15,0.1536120213245562,0.15119524913566243,0.07875264997066994,7.241115656602231e-16,0.047948399955734654,0.05513623287662578,0.03235759328251896,4.090564786249015e-16,0.023392967261331783,0.028681570420819397,0.017771008577657642,2.4672126293013086e-16,0.014016396066471176,0.017807158306087076,0.011386475028370351,1.714592334083309e-16,0.009474348163803117,0.012316262725232965,0.008041926341760328,1.1212044817922766e-16,0.006942375707716124,0.009173279000770123,0.006081161484994412,2.111598525624007e-17,0.00539520003414077,0.007217893995285857,0.004841101358177077,2.7775084398586677e-16,0.004388239086960148,0.005929248667932901,0.004014508967799147,9.486404471810101e-17,0.0037036469670079972,0.005045721135369522,0.0034434847131706546,7.352254058261595e-17,0.0032247860657257685,0.004424650744237962,0.0030404540187279,7.970041183113523e-17,0.002884979494813953,0.003983384372505874,0.0027540516952750064,7.285867231955118e-17,0.0026442785821523935,0.0036719561660329553,0.002552983887311817,1.6777299785974844e-16,0.002478058285030729,0.0034594522045109273,0.002417849611852019,9.246617731272175e-17,0.0023710817878942196,0.0033268345221317414,0.0023367953000158204,1.8921291035788944e-16,0.002314304311501306,0.0032630646073174008,0.002303167019442115],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":131904.75985510723,"average":2.609024107869118e-15,"dutyCycle":0.24999999999999975,"effectiveFrequency":131904.75985510723,"label":"Triangular","offset":-4.440892098500626e-16,"peak":4.999999999999999,"peakToPeak":10.0,"rms":2.8876908907560668,"thd":0.3765598442298381},"waveform":{"ancillaryLabel":null,"data":[-5.0,5.0,-5.0],"numberPeriods":null,"time":[0.0,2.4271844660194174e-06,9.70873786407767e-06]}},"frequency":103000.0,"magneticFieldStrength":null,"magneticFluxDensity":null,"magnetizingCurrent":{"harmonics":{"amplitudes":[0.006798274035336369,0.6948836814766789,0.24493767435427322,0.07659321568695032,0.0011851486650484215,0.02713826890868042,0.026338258337802543,0.013524003490593392,0.001185148665048509,0.007934821893277648,0.00885178693142459,0.005121342057388876,0.0011851486650485484,0.003521567475518845,0.004036010998935581,0.0025392588578627244,0.001185148665048555,0.0019074529966875186,0.002056445686787696,0.0014922821772143813,0.0011851486650485722,0.00122003158014941,0.001056889452408924,0.0010459731468581636,0.0011851486650485724,0.0009399767435062919,0.0004847443569363152,0.0008801360125169896,0.0011851486650485668,0.0008503270906284321,0.00012878834865636887,0.0008391415445710964,0.0011851486650485696,0.0008389428243523104,0.0001057951454006002,0.0008448591470041428,0.001185148665048572,0.000853888382471694,0.000266631468075954,0.0008642251603076917,0.0011851486650485711,0.0008748056010587637,0.0003796903985003957,0.0008850169925296278,0.0011851486650485696,0.0008945181777242711,0.0004600180147137569,0.000903129917240823,0.0011851486650485716,0.000910768242017427,0.0005167100425910805,0.0009174039220630434,0.00118514866504857,0.0009230377659956728,0.0005553940161691154,0.000927685539989642,0.0011851486650485783,0.0009313687576957952,0.0005795355845633968,0.0009341090721450906,0.001185148665048575,0.0009359248919968734,0.0005911441855611932,0.0009368293872760598],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":133347.9143873257,"average":-0.006798274035336401,"dutyCycle":0.25,"effectiveFrequency":133342.3181567079,"label":"Custom","offset":0.0020903409525253602,"peak":0.9028033263894203,"peakToPeak":1.8014259708737899,"rms":0.5247313861450679,"thd":0.3745218047109},"waveform":{"ancillaryLabel":null,"data":[-0.8417355085620355,-0.8606978872028123,0.9028033263894203,-0.8986226444843696,-0.8417355085620355],"numberPeriods":null,"time":[0.0,7.58495145631068e-08,2.427184466019418e-06,9.632888349514574e-06,9.708737864077681e-06]}},"name":"Primary winding excitation","voltage":{"harmonics":{"amplitudes":[0.078125,4.390192203074813,3.1805418164042503,1.6082999269932894,0.15625,0.7853068032814723,1.0533519383460939,0.747228683429235,0.15625,0.38152854838936,0.6237849662140766,0.5097284946091541,0.15625,0.22395874930175205,0.4366894957016387,0.3969081532881607,0.15625,0.1387657119725448,0.3303628683669758,0.32998229455569084,0.15625,0.08454863298834173,0.2606873758724238,0.2849610019773425,0.15625,0.046391828940503435,0.21067873648230084,0.25205943526842706,0.15625,0.017597888807606342,0.1723953087083574,0.22653191551027138,0.15625,0.005294556301253328,0.14161674515924202,0.20579082584016142,0.15625,0.02426175514926864,0.11588289785500591,0.18829805684603593,0.15625,0.040521304735682134,0.09365264588780098,0.17307466014248524,0.15625,0.05487174035875039,0.07390074623301954,0.1594602594300079,0.15625,0.06786649949748992,0.05590714395539467,0.1469848919856514,0.15625,0.07990999700987583,0.03913858752989208,0.13529576856110245,0.15625,0.09131439305868466,0.023177498052867507,0.1241125131778793,0.15625,0.10233553928853857,0.007676070276481761,0.11319769974196572],"frequencies":[0.0,103000.0,206000.0,309000.0,412000.0,515000.0,618000.0,721000.0,824000.0,927000.0,1030000.0,1133000.0,1236000.0,1339000.0,1442000.0,1545000.0,1648000.0,1751000.0,1854000.0,1957000.0,2060000.0,2163000.0,2266000.0,2369000.0,2472000.0,2575000.0,2678000.0,2781000.0,2884000.0,2987000.0,3090000.0,3193000.0,3296000.0,3399000.0,3502000.0,3605000.0,3708000.0,3811000.0,3914000.0,4017000.0,4120000.0,4223000.0,4326000.0,4429000.0,4532000.0,4635000.0,4738000.0,4841000.0,4944000.0,5047000.0,5150000.0,5253000.0,5356000.0,5459000.0,5562000.0,5665000.0,5768000.0,5871000.0,5974000.0,6077000.0,6180000.0,6283000.0,6386000.0,6489000.0]},"processed":{"acEffectiveFrequency":711072.940860665,"average":-0.078125,"dutyCycle":0.24218749999999975,"effectiveFrequency":711013.8104343135,"label":"Rectangular","offset":0.0,"peak":7.5,"peakToPeak":10.0,"rms":4.284784125250653,"thd":0.9507084995254541},"waveform":{"ancillaryLabel":null,"data":[-2.5,7.5,7.5,-2.5,-2.5],"numberPeriods":null,"time":[0.0,0.0,2.4271844660194174e-06,2.4271844660194174e-06,9.70873786407767e-06]}}}],"name":"Operating Point No. 2"}]})";
        std::string weightsString = R"({"AREA_PRODUCT":1,"ENERGY_STORED":1,"COST":1,"EFFICIENCY":1,"DIMENSIONS":1})";
        std::string file_path = __FILE__;

        std::map<std::string, double> weightsKeysString = json::parse(weightsString);
        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        for (auto const& pair : weightsKeysString) {
            weights[magic_enum::enum_cast<OpenMagnetics::CoreAdviser::CoreAdviserFilters>(pair.first).value()] = pair.second;
        }

        auto inventory_path = file_path.substr(0, file_path.rfind("/")).append("/testData/cores_stock.ndjson");

        std::ifstream ndjsonFile(inventory_path);
        std::string jsonLine;
        json coresJson = json::array();
        while (std::getline(ndjsonFile, jsonLine)) {
            json jf = json::parse(jsonLine);
            coresJson.push_back(jf);
        }

        std::string coresString = coresJson.dump(4);

        OpenMagnetics::InputsWrapper inputs(json::parse(inputsString));
        std::vector<OpenMagnetics::CoreWrapper> cores = json::parse(coresString);

        std::vector<OpenMagnetics::CoreWrapper> coresInInventory;
        for (auto& core : cores) {
            if (core.get_distributors_info()) {
                if (core.get_distributors_info().value().size() > 0) {
                    coresInInventory.push_back(core);
                }
            }
        }
        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &coresInInventory, 20);
        auto log = coreAdviser.read_log();
        std::vector<std::string> listOfNames;
        for (size_t i = 0; i < masMagnetics.size(); ++i ){
            listOfNames.push_back(masMagnetics[i].first.get_magnetic().get_manufacturer_info().value().get_reference().value());
        }

        auto scorings = coreAdviser.get_scorings();


        CHECK(bool(masMagnetics[0].first.get_outputs()[0].get_winding_losses().value().get_dc_resistance_per_turn()));
        CHECK(bool(masMagnetics[1].first.get_outputs()[0].get_winding_losses().value().get_dc_resistance_per_turn()));
        CHECK(bool(masMagnetics[2].first.get_outputs()[0].get_winding_losses().value().get_dc_resistance_per_turn()));
        CHECK(bool(masMagnetics[3].first.get_outputs()[0].get_winding_losses().value().get_dc_resistance_per_turn()));

        settings->reset();
    }

    TEST(Test_CoreAdviser_Web_6) {
        OpenMagnetics::clear_databases();
        std::string inputsString = R"({"designRequirements": {"magnetizingInductance": {"nominal": 0.0003 }, "name": "My Design Requirements", "turnsRatios": [] }, "operatingPoints": [{"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 65000, "current": {"waveform": {"ancillaryLabel": null, "data": [0, 0.5192640591887595, 1.0372573922108386, 1.5527123825794238, 2.0643676255575296, 2.5709710150267626, 3.0712828076008143, 3.5640786564852687, 4.048152607659378, 4.522320051047631, 4.985420619459117, 5.4363210282004895, 5.873917848413508, 6.297140207350368, 6.704952408978751, 7.096356468503582, 7.4703945546030495, 7.826151333402349, 8.162756208449071, 8.4793854512086, 8.775264216865875, 9.04966844050049, 9.301926608994888, 9.53142140433949, 9.737591214313307, 9.919931506843154, 10.077996064678178, 10.211398077358242, 10.31981108780395, 10.40296979121181, 10.460670684299556, 10.49277256331278, 10.499196869574178, 10.479927881729795, 10.435012754221821, 10.364561401893765, 10.268746231010402, 10.147801717350784, 10.002023832406682, 9.831769319090768, 9.637454818726928, 9.419555851459123, 9.178605652573903, 8.915193867584328, 8.629965109268564, 8.323617380194417, 7.996900364590075, 7.6506135937411655, 7.285604489403754, 6.902766290021262, 6.50303586482034, 6.08739142113488, 5.656850110568721, 5.212465539855277, 4.755325192505255, 4.286547767552297, 3.8072804419090045, 3.318696063032752, 2.8219902787714473, 2.3183786114129172, 1.8090934830981307, 1.2953811998777258, 0.7784989017922882, 0.2597114864399443, -0.25971148643994635, -0.778498901792281, -1.295381199877728, -1.8090934830981236, -2.3183786114129146, -2.8219902787714495, -3.3186960630327453, -3.8072804419090023, -4.286547767552299, -4.755325192505249, -5.212465539855274, -5.656850110568722, -6.087391421134874, -6.503035864820341, -6.902766290021259, -7.285604489403753, -7.650613593741167, -7.99690036459007, -8.323617380194417, -8.629965109268564, -8.915193867584328, -9.178605652573903, -9.419555851459119, -9.637454818726926, -9.831769319090766, -10.00202383240668, -10.147801717350784, -10.268746231010402, -10.364561401893761, -10.435012754221821, -10.479927881729795, -10.499196869574178, -10.49277256331278, -10.460670684299556, -10.40296979121181, -10.31981108780395, -10.211398077358242, -10.077996064678176, -9.919931506843156, -9.73759121431331, -9.531421404339488, -9.301926608994888, -9.04966844050049, -8.775264216865875, -8.479385451208605, -8.162756208449075, -7.826151333402345, -7.470394554603051, -7.0963564685035845, -6.704952408978752, -6.297140207350374, -5.873917848413506, -5.436321028200488, -4.985420619459123, -4.522320051047634, -4.048152607659381, -3.5640786564852793, -3.0712828076008147, -2.570971015026762, -2.0643676255575367, -1.55271238257943, -1.0372573922108435, -0.519264059188754, -2.5717582782094417e-15 ], "numberPeriods": null, "time": [0, 1.2113870381586917e-7, 2.4227740763173833e-7, 3.634161114476075e-7, 4.845548152634767e-7, 6.056935190793458e-7, 7.26832222895215e-7, 8.479709267110842e-7, 9.691096305269533e-7, 0.0000010902483343428223, 0.0000012113870381586916, 0.0000013325257419745608, 0.00000145366444579043, 0.000001574803149606299, 0.0000016959418534221685, 0.0000018170805572380375, 0.0000019382192610539067, 0.000002059357964869776, 0.0000021804966686856447, 0.0000023016353725015143, 0.000002422774076317383, 0.0000025439127801332524, 0.0000026650514839491216, 0.000002786190187764991, 0.00000290732889158086, 0.0000030284675953967292, 0.000003149606299212598, 0.0000032707450030284677, 0.000003391883706844337, 0.0000035130224106602057, 0.000003634161114476075, 0.0000037552998182919437, 0.000003876438522107813, 0.000003997577225923682, 0.000004118715929739552, 0.000004239854633555421, 0.000004360993337371289, 0.000004482132041187159, 0.000004603270745003029, 0.0000047244094488188975, 0.000004845548152634766, 0.000004966686856450636, 0.000005087825560266505, 0.0000052089642640823735, 0.000005330102967898243, 0.000005451241671714112, 0.000005572380375529982, 0.000005693519079345851, 0.00000581465778316172, 0.000005935796486977589, 0.0000060569351907934585, 0.000006178073894609327, 0.000006299212598425196, 0.000006420351302241066, 0.000006541490006056935, 0.000006662628709872804, 0.000006783767413688674, 0.000006904906117504543, 0.000007026044821320411, 0.00000714718352513628, 0.00000726832222895215, 0.000007389460932768019, 0.0000075105996365838874, 0.000007631738340399757, 0.000007752877044215627, 0.000007874015748031496, 0.000007995154451847364, 0.000008116293155663234, 0.000008237431859479104, 0.000008358570563294972, 0.000008479709267110841, 0.00000860084797092671, 0.000008721986674742579, 0.000008843125378558448, 0.000008964264082374318, 0.000009085402786190188, 0.000009206541490006057, 0.000009327680193821927, 0.000009448818897637795, 0.000009569957601453665, 0.000009691096305269533, 0.000009812235009085402, 0.000009933373712901272, 0.00001005451241671714, 0.00001017565112053301, 0.000010296789824348879, 0.000010417928528164747, 0.000010539067231980617, 0.000010660205935796486, 0.000010781344639612354, 0.000010902483343428224, 0.000011023622047244094, 0.000011144760751059963, 0.000011265899454875833, 0.000011387038158691702, 0.00001150817686250757, 0.00001162931556632344, 0.00001175045427013931, 0.000011871592973955178, 0.000011992731677771047, 0.000012113870381586917, 0.000012235009085402785, 0.000012356147789218655, 0.000012477286493034524, 0.000012598425196850392, 0.000012719563900666262, 0.000012840702604482131, 0.000012961841308298, 0.00001308298001211387, 0.00001320411871592974, 0.000013325257419745608, 0.000013446396123561478, 0.000013567534827377348, 0.000013688673531193216, 0.000013809812235009085, 0.000013930950938824953, 0.000014052089642640823, 0.000014173228346456692, 0.00001429436705027256, 0.00001441550575408843, 0.0000145366444579043, 0.000014657783161720168, 0.000014778921865536037, 0.000014900060569351907, 0.000015021199273167775, 0.000015142337976983646, 0.000015263476680799514, 0.000015384615384615384 ] }, "processed": {"dutyCycle": 0.25, "peakToPeak": 21, "offset": 0, "label": "Sinusoidal", "acEffectiveFrequency": 65000.000360815815, "effectiveFrequency": 65000.000360815815, "peak": 10.499196869574178, "rms": 7.395561908164924, "thd": 0.00005759172987244173 }, "harmonics": {"amplitudes": [1.851795626561792e-14, 10.497697022808778, 0.0005893000543027365, 0.0001234861844609525, 0.000046746660945999526, 0.00002280727141528787, 0.00001286335235090736, 0.000007976948556476238, 0.000005291078413078924, 0.0000036908418258942075, 0.000002677459292247501, 0.0000020042379606232204, 0.0000015393727539678285, 0.0000012079796198210968, 9.653252703145028e-7, 7.835471678152572e-7, 6.446657504660653e-7, 5.367210786865906e-7, 4.515448252188785e-7, 3.8342773867699213e-7, 3.28296907643635e-7, 2.831933622553178e-7, 2.4593201140179206e-7, 2.1487530029300243e-7, 1.887795081660895e-7, 1.6668846618496374e-7, 1.4785884454273265e-7, 1.3170683238510115e-7, 1.1776956206458994e-7, 1.056768351203901e-7, 9.513016376589268e-8, 8.588706124079704e-8, 7.774911777065636e-8, 7.055296075836333e-8, 6.41631899756511e-8, 5.84669460068633e-8, 5.336963842965032e-8, 4.879157084627916e-8, 4.466526231691933e-8, 4.093330047139088e-8, 3.7546619060158496e-8, 3.446310175801058e-8, 3.16464406754751e-8, 2.9065208779999603e-8, 2.6692094029812346e-8, 2.4503265897977743e-8, 2.24778489419825e-8, 2.0597488381568804e-8, 1.8845967199863474e-8, 1.7208923027832514e-8, 1.567355948695751e-8, 1.4228439782667332e-8, 1.286328939504102e-8, 1.1568836358136048e-8, 1.033666883870453e-8, 9.159114761840138e-9, 8.029130556351e-9, 6.9402166321552666e-9, 5.886324569094645e-9, 4.861790840873707e-9, 3.861265642846226e-9, 2.8796530980210227e-9, 1.912056254592697e-9, 9.537330925014646e-10 ], "frequencies": [0, 65000, 130000, 195000, 260000, 325000, 390000, 455000, 520000, 585000, 650000, 715000, 780000, 845000, 910000, 975000, 1040000, 1105000, 1170000, 1235000, 1300000, 1365000, 1430000, 1495000, 1560000, 1625000, 1690000, 1755000, 1820000, 1885000, 1950000, 2015000, 2080000, 2145000, 2210000, 2275000, 2340000, 2405000, 2470000, 2535000, 2600000, 2665000, 2730000, 2795000, 2860000, 2925000, 2990000, 3055000, 3120000, 3185000, 3250000, 3315000, 3380000, 3445000, 3510000, 3575000, 3640000, 3705000, 3770000, 3835000, 3900000, 3965000, 4030000, 4095000 ] } } } ] } ] })";
        std::string weightsString = R"({"AREA_PRODUCT":1,"ENERGY_STORED":1,"COST":1,"EFFICIENCY":1,"DIMENSIONS":1})";
        std::string file_path = __FILE__;

        std::map<std::string, double> weightsKeysString = json::parse(weightsString);
        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        for (auto const& pair : weightsKeysString) {
            weights[magic_enum::enum_cast<OpenMagnetics::CoreAdviser::CoreAdviserFilters>(pair.first).value()] = pair.second;
        }

        OpenMagnetics::InputsWrapper inputs(json::parse(inputsString));
        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;

        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, 2);
        auto log = coreAdviser.read_log();
        std::vector<std::string> listOfNames;
        for (size_t i = 0; i < masMagnetics.size(); ++i ){
            listOfNames.push_back(masMagnetics[i].first.get_magnetic().get_manufacturer_info().value().get_reference().value());
        }

        CHECK(masMagnetics.size() > 0);
        settings->reset();
    }

    TEST(Test_CoreAdviser_User_0) {
        OpenMagnetics::clear_databases();
        std::string inputsString = R"({"designRequirements": {"magnetizingInductance": {"nominal": 0.000006 }, "name": "My Design Requirements", "turnsRatios": [] }, "operatingPoints": [{"name": "Operating Point No. 1", "conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"name": "Primary winding excitation", "frequency": 200000, "current": {"waveform": {"ancillaryLabel": null, "data": [56.985, 67.415, 56.985 ], "numberPeriods": null, "time": [0, 0.0000042000000000000004, 0.000005 ] }, "processed": {"dutyCycle": 0.84, "peakToPeak": 10.43, "offset": 62.2, "label": "Triangular", "acEffectiveFrequency": 300258.7073651182, "effectiveFrequency": 20509.35261535522, "peak": 67.3645572916667, "rms": 62.272256562561566, "thd": 0.5135395888682869 }, "harmonics": {"amplitudes": [62.19940887451172, 3.788193673170306, 1.6600604274902384, 0.8723324955604976, 0.44502354588734777, 0.18510789567277652, 0.027454438456555246, 0.05922378818418705, 0.09496936245070094, 0.09574347587431536, 0.07516282442281835, 0.04477150722215702, 0.013737807079460167, 0.01171857202338386, 0.02774105038565934, 0.03361420481778365, 0.030566703298216306, 0.021295758815253026, 0.009168505436547483, 0.0030289326798604603, 0.011820093941227659, 0.01649770608447672, 0.016617062490040618, 0.012919847801793128, 0.006883846185546486, 0.0013008813280683239, 0.005890359494431914, 0.009462239200902774, 0.010413922344839593, 0.008859105864681853, 0.005510304072101526, 0.0017545672283819061, 0.003215298779201598, 0.00592474378040423, 0.007083641896483774, 0.00652395008903164, 0.004588941725857998, 0.002117512227581541, 0.0020138773922366935, 0.003943672284408362, 0.005077227310080852, 0.005024234574379809, 0.003921795582955405, 0.002305569711275414, 0.001669940867990107, 0.002807749879482047, 0.003783929724516936, 0.003983860219801076, 0.0034083765042437313, 0.0023989610715256663, 0.001768537738394057, 0.002235857850988618, 0.00293953534824593, 0.003224018804477077, 0.0029911395353770394, 0.0024433457016505566, 0.002018570807359561, 0.0020956739117972585, 0.002441654302943123, 0.002660136158579102, 0.0026339060269650392, 0.0024610090727368076, 0.0023132303486381877, 0.0022704164604682078 ], "frequencies": [0, 200000, 400000, 600000, 800000, 1000000, 1200000, 1400000, 1600000, 1800000, 2000000, 2200000, 2400000, 2600000, 2800000, 3000000, 3200000, 3400000, 3600000, 3800000, 4000000, 4200000, 4400000, 4600000, 4800000, 5000000, 5200000, 5400000, 5600000, 5800000, 6000000, 6200000, 6400000, 6600000, 6800000, 7000000, 7200000, 7400000, 7600000, 7800000, 8000000, 8200000, 8400000, 8600000, 8800000, 9000000, 9200000, 9400000, 9600000, 9800000, 10000000, 10200000, 10400000, 10600000, 10800000, 11000000, 11200000, 11400000, 11600000, 11800000, 12000000, 12200000, 12400000, 12600000 ] } }, "voltage": {"waveform": {"ancillaryLabel": null, "data": [-63.7224, 12.137600000000003, 12.137600000000003, -63.7224, -63.7224 ], "numberPeriods": null, "time": [0, 0, 0.0000042000000000000004, 0.0000042000000000000004, 0.000005 ] }, "processed": {"dutyCycle": 0.84, "peakToPeak": 75.86, "offset": 0, "label": "Rectangular", "acEffectiveFrequency": 1597694.1978822167, "effectiveFrequency": 1597646.1127299606, "peak": 12.137600000000003, "rms": 28.09509691004821, "thd": 1.3356374489325038 }, "harmonics": {"amplitudes": [0.30818124999999846, 23.80640349491498, 20.719882427214415, 16.10769435790146, 10.664994868582662, 5.180425361516884, 0.39637660780488254, 3.1172440421995313, 5.0517694299099345, 5.395226624849836, 4.40986507390017, 2.558980130353162, 0.4002310516205373, 1.53128962621708, 2.8260035663636964, 3.268726375557206, 2.861597513150363, 1.7995245536616769, 0.4067814140804575, 0.9487930071081266, 1.9437108520751825, 2.36937753861694, 2.1708188924640415, 1.4471352892692118, 0.41622639870013567, 0.6457200391696907, 1.4743321406437133, 1.8797993287695682, 1.7879655570329085, 1.251308954055062, 0.42886397748534805, 0.45847625314767393, 1.1853125000000002, 1.577328849337885, 1.5517759950109278, 1.1337151229433782, 0.44511436425132234, 0.3294280656605511, 0.9910355573978922, 1.376600449830656, 1.3981714250008515, 1.06234872857457, 0.46555489661099314, 0.23290000407082775, 0.8526319979265212, 1.2380854678815314, 1.297010142583099, 1.0221252498854263, 0.49097251315036106, 0.1555155662812027, 0.7499275477021183, 1.141223897515662, 1.2326838490208047, 1.0055172910305796, 0.522443369749807, 0.08936846001174778, 0.6714255963907871, 1.0745701353211696, 1.1968387110459648, 1.0090036687405854, 0.5614560087628049, 0.029168057047719874, 0.6101073099508907, 1.0316356960324993 ], "frequencies": [0, 200000, 400000, 600000, 800000, 1000000, 1200000, 1400000, 1600000, 1800000, 2000000, 2200000, 2400000, 2600000, 2800000, 3000000, 3200000, 3400000, 3600000, 3800000, 4000000, 4200000, 4400000, 4600000, 4800000, 5000000, 5200000, 5400000, 5600000, 5800000, 6000000, 6200000, 6400000, 6600000, 6800000, 7000000, 7200000, 7400000, 7600000, 7800000, 8000000, 8200000, 8400000, 8600000, 8800000, 9000000, 9200000, 9400000, 9600000, 9800000, 10000000, 10200000, 10400000, 10600000, 10800000, 11000000, 11200000, 11400000, 11600000, 11800000, 12000000, 12200000, 12400000, 12600000 ] } } } ] } ]})";
        std::string weightsString = R"({"AREA_PRODUCT":1,"ENERGY_STORED":1,"COST":1,"EFFICIENCY":0,"DIMENSIONS":1})";
        std::string file_path = __FILE__;

        std::map<std::string, double> weightsKeysString = json::parse(weightsString);
        std::map<OpenMagnetics::CoreAdviser::CoreAdviserFilters, double> weights;
        for (auto const& pair : weightsKeysString) {
            weights[magic_enum::enum_cast<OpenMagnetics::CoreAdviser::CoreAdviserFilters>(pair.first).value()] = pair.second;
        }

        OpenMagnetics::InputsWrapper inputs(json::parse(inputsString));

        settings->set_use_toroidal_cores(false);
        OpenMagnetics::CoreAdviser coreAdviser;
        auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, 20);
        // auto masMagnetics = coreAdviser.get_advised_core(inputs, weights, &coresInInventory, 20);
        auto log = coreAdviser.read_log();
        std::vector<std::string> listOfNames;
        for (size_t i = 0; i < masMagnetics.size(); ++i ){
            listOfNames.push_back(masMagnetics[i].first.get_magnetic().get_manufacturer_info().value().get_reference().value());
        }

        CHECK(masMagnetics.size() > 0);
        settings->reset();
    }

}