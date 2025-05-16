#include "support/Utils.h"
#include "support/Settings.h"
#include "json.hpp"

#include <UnitTest++.h>
#include <filesystem>
#include <cfloat>
#include <limits>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
using json = nlohmann::json;
#include <typeinfo>

SUITE(Utils) {
    TEST(CeilFloat) {
        double value = 1.263;
        double calculatedValue = OpenMagnetics::ceilFloat(value, 2);
        double expectedValue = 1.27;
        CHECK_EQUAL(expectedValue, calculatedValue);
    }

    TEST(FloorFloat) {
        double value = 1.263;
        double calculatedValue = OpenMagnetics::floorFloat(value, 2);
        double expectedValue = 1.26;
        CHECK_EQUAL(expectedValue, calculatedValue);
    }

    TEST(Modified_Bessel) {
        double calculatedValue = OpenMagnetics::modified_bessel_first_kind(0.0, std::complex<double>{1.0, 0.0}).real();
        double expectedValue = 1.2660658777520084;
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    }

    TEST(Bessel) {
        double calculatedValue = OpenMagnetics::bessel_first_kind(0.0, std::complex<double>{1.0, 0.0}).real();
        double expectedValue = 0.7651976865579666;
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);

        double calculatedBerValue = OpenMagnetics::kelvin_function_real(0.0, 1.0);
        double expectedBerValue = 0.98438178;
        CHECK_CLOSE(expectedBerValue, calculatedBerValue, expectedBerValue * 0.001);

        double calculatedBeiValue = OpenMagnetics::kelvin_function_imaginary(0.0, 1.0);
        double expectedBeiValue = 0.24956604;
        CHECK_CLOSE(expectedBeiValue, calculatedBeiValue, expectedBeiValue * 0.001);

        double calculatedBerpValue = OpenMagnetics::derivative_kelvin_function_real(0.0, 1.0);
        double expectedBerpValue = -0.06244575217903096;
        CHECK_CLOSE(expectedBerpValue, calculatedBerpValue, fabs(expectedBerpValue) * 0.001);

        double calculatedBeipValue = OpenMagnetics::derivative_kelvin_function_imaginary(0.0, 1.0);
        double expectedBeipValue = 0.49739651146809727;
        CHECK_CLOSE(expectedBeipValue, calculatedBeipValue, expectedBeipValue * 0.001);
    }

    TEST(Test_Complete_Ellipitical_1_0) {
        double calculatedValue = OpenMagnetics::comp_ellint_1(0);
        double expectedValue = std::comp_ellint_1(0);
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    }

    TEST(Test_Complete_Ellipitical_1_1) {
        double calculatedValue = OpenMagnetics::comp_ellint_1(1);
        double expectedValue = std::comp_ellint_1(1);
        CHECK(std::isnan(calculatedValue));
        CHECK(std::isnan(expectedValue));
    }

    TEST(Test_Complete_Ellipitical_1_2) {
        double calculatedValue = OpenMagnetics::comp_ellint_1(std::sin(std::numbers::pi / 18 / 2));
        double expectedValue = std::comp_ellint_1(std::sin(std::numbers::pi / 18 / 2));
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    }

    TEST(Test_Complete_Ellipitical_2_0) {
        double calculatedValue = OpenMagnetics::comp_ellint_2(0);
        double expectedValue = std::comp_ellint_2(0);
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    }

    TEST(Test_Complete_Ellipitical_2_1) {
        double calculatedValue = OpenMagnetics::comp_ellint_2(1);
        double expectedValue = std::comp_ellint_2(1);
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    }

    TEST(Test_Complete_Ellipitical_2_2) {
        double calculatedValue = OpenMagnetics::comp_ellint_2(std::sin(std::numbers::pi / 18 / 2));
        double expectedValue = std::comp_ellint_2(std::sin(std::numbers::pi / 18 / 2));
        CHECK_CLOSE(expectedValue, calculatedValue, expectedValue * 0.001);
    }

    TEST(Test_Find_By_Perimeter) {
        OpenMagnetics::clear_databases();
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_use_toroidal_cores(true);
        settings->set_use_concentric_cores(true);

        auto shape = OpenMagnetics::find_core_shape_by_winding_window_perimeter(0.03487);

        CHECK_EQUAL("UR 46/21/11", shape.get_name().value());
    }

    TEST(Test_Find_By_Perimeter_Only_Toroids) {
        OpenMagnetics::clear_databases();
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_use_toroidal_cores(true);
        settings->set_use_concentric_cores(false);

        auto shape = OpenMagnetics::find_core_shape_by_winding_window_perimeter(0.03487);

        CHECK_EQUAL("T 22/12.4/12.8", shape.get_name().value());
    }

    TEST(Test_Get_Shapes_By_Manufacturer) {
        OpenMagnetics::clear_databases();
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->reset();
        auto allShapeNames = OpenMagnetics::get_shape_names();
        auto magneticsShapeNames = OpenMagnetics::get_core_shapes_names("Magnetics");
        auto ferroxcubeShapeNames = OpenMagnetics::get_core_shapes_names("Ferroxcube");

        CHECK(allShapeNames.size() > magneticsShapeNames.size());
        CHECK(allShapeNames.size() > ferroxcubeShapeNames.size());
    }

    TEST(Test_Wire_Names_With_Types) {
        OpenMagnetics::clear_databases();
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->reset();
        auto allWireNames = OpenMagnetics::get_wire_names();
        for (auto wire : allWireNames) {
            CHECK(
                wire.starts_with("Round ") ||
                wire.starts_with("Litz ") ||
                wire.starts_with("Rectangular ") ||
                wire.starts_with("Planar ") ||
                wire.starts_with("Foil ")
                );
        }
    }

    TEST(Test_Wires_With_Type_And_Standard) {
        OpenMagnetics::clear_databases();
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->reset();
        auto allWires = OpenMagnetics::get_wires();
        auto someWires = OpenMagnetics::get_wires(OpenMagnetics::WireType::ROUND, OpenMagnetics::WireStandard::IEC_60317);
        CHECK(someWires.size() < allWires.size());
    }

    TEST(Test_Load_Toroidal_Cores) {
        OpenMagnetics::clear_databases();
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->reset();
        settings->set_use_toroidal_cores(true);
        settings->set_use_concentric_cores(false);
        OpenMagnetics::load_cores();

        auto allCores = coreDatabase;
        CHECK(allCores.size() > 0);
        for (auto core : allCores) {
            CHECK(core.get_type() == OpenMagnetics::CoreType::TOROIDAL);
        }
    }

    TEST(Test_Load_Two_Piece_Set_Cores) {
        OpenMagnetics::clear_databases();
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->reset();
        settings->set_use_toroidal_cores(false);
        settings->set_use_concentric_cores(true);
        OpenMagnetics::load_cores();

        auto allCores = coreDatabase;
        CHECK(allCores.size() > 0);
        for (auto core : allCores) {
            CHECK(core.get_type() == OpenMagnetics::CoreType::TWO_PIECE_SET);
        }
    }

    TEST(Test_Load_Cores_In_Stock) {
        OpenMagnetics::clear_databases();
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->reset();
        settings->set_use_only_cores_in_stock(false);
        OpenMagnetics::load_cores();

        auto allCores = coreDatabase;

        OpenMagnetics::clear_databases();
        settings->reset();
        settings->set_use_only_cores_in_stock(true);
        OpenMagnetics::load_cores();

        auto onlyCoresInStock = coreDatabase;
        CHECK(allCores.size() > onlyCoresInStock.size());
    }

    TEST(Test_Core_Materials_External) {
        std::string file_path = __FILE__;
        auto external_core_materials_path = file_path.substr(0, file_path.rfind("/")).append("/testData/external_core_materials.ndjson");

        std::ifstream file(external_core_materials_path, std::ios_base::binary | std::ios_base::in);
        if(!file.is_open())
            throw std::runtime_error("Failed to open " + external_core_materials_path);
        using Iterator = std::istreambuf_iterator<char>;
        std::string external_core_materials(Iterator{file}, Iterator{});

        OpenMagnetics::clear_databases();
        OpenMagnetics::load_core_materials(external_core_materials);

        auto allCoreMaterials = coreMaterialDatabase;
        CHECK(allCoreMaterials.size() == 4);
    }

    TEST(Test_Core_Materials_External_And_Internal) {
        std::string file_path = __FILE__;
        auto external_core_materials_path = file_path.substr(0, file_path.rfind("/")).append("/testData/external_core_materials.ndjson");

        std::ifstream file(external_core_materials_path, std::ios_base::binary | std::ios_base::in);
        if(!file.is_open())
            throw std::runtime_error("Failed to open " + external_core_materials_path);
        using Iterator = std::istreambuf_iterator<char>;
        std::string external_core_materials(Iterator{file}, Iterator{});

        OpenMagnetics::clear_databases();
        OpenMagnetics::load_core_materials();

        auto allCoreMaterials = coreMaterialDatabase;

        OpenMagnetics::load_core_materials(external_core_materials);

        auto allCoreMaterialsWithExternal = coreMaterialDatabase;

        CHECK(allCoreMaterialsWithExternal.size() > allCoreMaterials.size());
    }

    TEST(Test_Core_Shapes_Families) {
        CHECK(OpenMagnetics::get_shape_families().size() > 0);
    }

    TEST(Test_Core_Shapes_Dimensions_PQ) {
        auto dimensions = OpenMagnetics::get_shape_family_dimensions(OpenMagnetics::CoreShapeFamily::PQ);
        CHECK(dimensions.size() == 9);
        CHECK(dimensions[0] == "A");
        CHECK(dimensions[1] == "B");
        CHECK(dimensions[2] == "C");
        CHECK(dimensions[3] == "D");
        CHECK(dimensions[4] == "E");
        CHECK(dimensions[5] == "F");
        CHECK(dimensions[6] == "G");
        CHECK(dimensions[7] == "J");
        CHECK(dimensions[8] == "L");
    }

    TEST(Test_Core_Shapes_Dimensions_UR) {
        auto dimensions = OpenMagnetics::get_shape_family_dimensions(OpenMagnetics::CoreShapeFamily::UR, "2");
        CHECK(dimensions.size() == 7);
        CHECK(dimensions[0] == "A");
        CHECK(dimensions[1] == "B");
        CHECK(dimensions[2] == "C");
        CHECK(dimensions[3] == "D");
        CHECK(dimensions[4] == "E");
        CHECK(dimensions[5] == "H");
        CHECK(dimensions[6] == "S");
    }

    TEST(Test_Core_Shapes_Dimensions_PM) {
        auto dimensions = OpenMagnetics::get_shape_family_dimensions(OpenMagnetics::CoreShapeFamily::PM, "1");
        CHECK(dimensions.size() == 12);
        CHECK(dimensions[0] == "A");
        CHECK(dimensions[1] == "B");
        CHECK(dimensions[2] == "C");
        CHECK(dimensions[3] == "D");
        CHECK(dimensions[4] == "E");
        CHECK(dimensions[5] == "F");
        CHECK(dimensions[6] == "G");
        CHECK(dimensions[7] == "H");
        CHECK(dimensions[8] == "alpha");
        CHECK(dimensions[9] == "b");
        CHECK(dimensions[10] == "e");
        CHECK(dimensions[11] == "t");
    }

    TEST(Test_Core_Shapes_Dimensions_UR_No_Subtype) {
        auto dimensions = OpenMagnetics::get_shape_family_dimensions(OpenMagnetics::CoreShapeFamily::UR);
        CHECK(dimensions.size() == 9);
        CHECK(dimensions[0] == "A");
        CHECK(dimensions[1] == "B");
        CHECK(dimensions[2] == "C");
        CHECK(dimensions[3] == "D");
        CHECK(dimensions[4] == "E");
        CHECK(dimensions[5] == "F");
        CHECK(dimensions[6] == "G");
        CHECK(dimensions[7] == "H");
        CHECK(dimensions[8] == "S");
    }

    TEST(Test_Core_Shapes_Family_Subtypes_UR) {
        auto subtypes = OpenMagnetics::get_shape_family_subtypes(OpenMagnetics::CoreShapeFamily::UR);
        CHECK(subtypes.size() == 4);
        CHECK(subtypes[0] == "1");
        CHECK(subtypes[1] == "2");
        CHECK(subtypes[2] == "3");
        CHECK(subtypes[3] == "4");
    }

    TEST(Test_Core_Shapes_Family_Subtypes_E) {
        auto subtypes = OpenMagnetics::get_shape_family_subtypes(OpenMagnetics::CoreShapeFamily::E);
        CHECK(subtypes.size() == 0);
    }

    TEST(Test_Core_Shapes_External) {
        std::string file_path = __FILE__;
        auto external_core_shapes_path = file_path.substr(0, file_path.rfind("/")).append("/testData/external_core_shapes.ndjson");

        std::ifstream file(external_core_shapes_path, std::ios_base::binary | std::ios_base::in);
        if(!file.is_open())
            throw std::runtime_error("Failed to open " + external_core_shapes_path);
        using Iterator = std::istreambuf_iterator<char>;
        std::string external_core_shapes(Iterator{file}, Iterator{});
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->reset();
        settings->set_use_toroidal_cores(true);
        settings->set_use_concentric_cores(true);

        OpenMagnetics::clear_databases();
        OpenMagnetics::load_core_shapes(false, external_core_shapes);

        auto allCoreShapes = coreShapeDatabase;
        CHECK(allCoreShapes.size() == 10);
    }

    TEST(Test_Core_Shapes_External_And_Internal) {
        std::string file_path = __FILE__;
        auto external_core_shapes_path = file_path.substr(0, file_path.rfind("/")).append("/testData/external_core_shapes.ndjson");

        std::ifstream file(external_core_shapes_path, std::ios_base::binary | std::ios_base::in);
        if(!file.is_open())
            throw std::runtime_error("Failed to open " + external_core_shapes_path);
        using Iterator = std::istreambuf_iterator<char>;
        std::string external_core_shapes(Iterator{file}, Iterator{});
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->reset();
        settings->set_use_toroidal_cores(true);
        settings->set_use_concentric_cores(true);

        OpenMagnetics::clear_databases();
        OpenMagnetics::load_core_shapes(false);

        auto allCoreShapes = coreShapeDatabase;
        std::cout << "allCoreShapes.size(): " << allCoreShapes.size() << std::endl;

        OpenMagnetics::load_core_shapes(false, external_core_shapes);

        auto allCoreShapesWithExternal = coreShapeDatabase;
        std::cout << "allCoreShapesWithExternal.size(): " << allCoreShapesWithExternal.size() << std::endl;

        CHECK(allCoreShapesWithExternal.size() > allCoreShapes.size());
    }

    TEST(Test_Wires_External) {
        std::string file_path = __FILE__;
        auto external_wires_path = file_path.substr(0, file_path.rfind("/")).append("/testData/external_wires.ndjson");

        std::ifstream file(external_wires_path, std::ios_base::binary | std::ios_base::in);
        if(!file.is_open())
            throw std::runtime_error("Failed to open " + external_wires_path);
        using Iterator = std::istreambuf_iterator<char>;
        std::string external_wires(Iterator{file}, Iterator{});

        OpenMagnetics::clear_databases();
        OpenMagnetics::load_wires(external_wires);

        auto allCoreWires = wireDatabase;
        CHECK(allCoreWires.size() == 24);
    }

    TEST(Test_Wires_External_And_Internal) {
        std::string file_path = __FILE__;
        auto external_wires_path = file_path.substr(0, file_path.rfind("/")).append("/testData/external_wires.ndjson");

        std::ifstream file(external_wires_path, std::ios_base::binary | std::ios_base::in);
        if(!file.is_open())
            throw std::runtime_error("Failed to open " + external_wires_path);
        using Iterator = std::istreambuf_iterator<char>;
        std::string external_wires(Iterator{file}, Iterator{});

        OpenMagnetics::clear_databases();
        OpenMagnetics::load_wires();

        auto allCoreWires = wireDatabase;

        OpenMagnetics::load_wires(external_wires);

        auto allCoreWiresWithExternal = wireDatabase;

        CHECK(allCoreWiresWithExternal.size() > allCoreWires.size());
    }

    TEST(Test_Low_And_High_Harmonics) {
        OpenMagnetics::Harmonics harmonics = json::parse(R"({"amplitudes":[2.8679315866586563e-15,2.3315003998761155,1.1284261163900786],"frequencies":[0,50,400000000]})");
        auto mainHarmonicIndexes = OpenMagnetics::get_main_harmonic_indexes(harmonics, 0.05, 1);

        CHECK(mainHarmonicIndexes.size() == 2);
    }

    TEST(Test_Mas_Autocomplete) {
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_use_only_cores_in_stock(false);
        auto core = OpenMagnetics::find_core_by_name("PQ 32/30 - 3C90 - Gapped 0.492 mm");
        OpenMagnetics::CoilWrapper coil;
        OpenMagnetics::CoilFunctionalDescription dummyWinding;
        dummyWinding.set_name("");
        dummyWinding.set_number_turns(0);
        dummyWinding.set_number_parallels(0);
        dummyWinding.set_wire("");
        coil.get_mutable_functional_description().push_back(dummyWinding);
        OpenMagnetics::MagneticWrapper magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        OpenMagnetics::InputsWrapper inputs = OpenMagnetics::InputsWrapper::create_quick_operating_point(100000, 42e-6, 42, OpenMagnetics::WaveformLabel::SINUSOIDAL, 10, 0.5, 2);

        OpenMagnetics::MasWrapper mas;
        mas.set_inputs(inputs);
        mas.set_magnetic(magnetic);

        auto autocompletedMas = OpenMagnetics::mas_autocomplete(mas);

        CHECK(autocompletedMas.get_magnetic().get_core().get_geometrical_description());
        CHECK(autocompletedMas.get_magnetic().get_core().get_processed_description());
        CHECK(std::holds_alternative<OpenMagnetics::CoreShape>(autocompletedMas.get_magnetic().get_core().get_functional_description().get_shape()));
        CHECK(std::holds_alternative<OpenMagnetics::CoreMaterial>(autocompletedMas.get_magnetic().get_core().get_functional_description().get_material()));
        CHECK(std::holds_alternative<OpenMagnetics::Bobbin>(autocompletedMas.get_magnetic().get_coil().get_bobbin()));
        CHECK(autocompletedMas.get_mutable_magnetic().get_mutable_coil().resolve_bobbin().get_processed_description());
        CHECK(autocompletedMas.get_magnetic().get_coil().get_sections_description());
        CHECK(autocompletedMas.get_magnetic().get_coil().get_layers_description());
        CHECK(autocompletedMas.get_magnetic().get_coil().get_turns_description());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_harmonics());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_processed());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_waveform());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_harmonics());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_harmonics());
        CHECK_EQUAL(autocompletedMas.get_inputs().get_operating_points().size(), autocompletedMas.get_outputs().size());
        CHECK(autocompletedMas.get_outputs()[0].get_core_losses());
        CHECK(autocompletedMas.get_outputs()[0].get_winding_losses());
        CHECK(autocompletedMas.get_outputs()[0].get_magnetizing_inductance());

    }


    TEST(Test_Mas_Autocomplete_Json) {
        std::string masString = R"({"outputs": [], "inputs": {"designRequirements": {"isolationSides": ["primary" ], "magnetizingInductance": {"nominal": 0.00039999999999999996 }, "name": "My Design Requirements", "turnsRatios": [{"nominal": 1} ] }, "operatingPoints": [{"conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"frequency": 100000, "current": {"processed": {"label": "Triangular", "peakToPeak": 0.5, "offset": 0, "dutyCycle": 0.5 } }, "voltage": {"processed": {"label": "Rectangular", "peakToPeak": 20, "offset": 0, "dutyCycle": 0.5 } } } ], "name": "Operating Point No. 1" } ] }, "magnetic": {"coil": {"bobbin": "Basic", "functionalDescription":[{"name": "Primary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] }, "core": {"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } }, "manufacturerInfo": {"name": "", "reference": "Example" } } })";
        json masJson = json::parse(masString);
        OpenMagnetics::MasWrapper mas(masJson);

        auto autocompletedMas = OpenMagnetics::mas_autocomplete(mas);

        CHECK(autocompletedMas.get_magnetic().get_core().get_geometrical_description());
        CHECK(autocompletedMas.get_magnetic().get_core().get_processed_description());
        CHECK(std::holds_alternative<OpenMagnetics::CoreShape>(autocompletedMas.get_magnetic().get_core().get_functional_description().get_shape()));
        CHECK(std::holds_alternative<OpenMagnetics::CoreMaterial>(autocompletedMas.get_magnetic().get_core().get_functional_description().get_material()));
        CHECK(std::holds_alternative<OpenMagnetics::Bobbin>(autocompletedMas.get_magnetic().get_coil().get_bobbin()));
        CHECK(autocompletedMas.get_mutable_magnetic().get_mutable_coil().resolve_bobbin().get_processed_description());
        CHECK(autocompletedMas.get_magnetic().get_coil().get_sections_description());
        CHECK(autocompletedMas.get_magnetic().get_coil().get_layers_description());
        CHECK(autocompletedMas.get_magnetic().get_coil().get_turns_description());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_harmonics());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_processed());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_waveform());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_harmonics());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform());
        CHECK(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_harmonics());
        CHECK_EQUAL(autocompletedMas.get_inputs().get_operating_points().size(), autocompletedMas.get_outputs().size());
        CHECK(autocompletedMas.get_outputs()[0].get_core_losses());
        CHECK(autocompletedMas.get_outputs()[0].get_winding_losses());
        CHECK(autocompletedMas.get_outputs()[0].get_magnetizing_inductance());

    }

    TEST(Test_Mas_Autocomplete_Only_Magnetic) {
        auto settings = OpenMagnetics::Settings::GetInstance();
        settings->set_use_only_cores_in_stock(false);
        auto core = OpenMagnetics::find_core_by_name("PQ 32/30 - 3C90 - Gapped 0.492 mm");
        OpenMagnetics::CoilWrapper coil;
        OpenMagnetics::CoilFunctionalDescription dummyWinding;
        dummyWinding.set_name("");
        dummyWinding.set_number_turns(0);
        dummyWinding.set_number_parallels(0);
        dummyWinding.set_wire("");
        coil.get_mutable_functional_description().push_back(dummyWinding);
        OpenMagnetics::MagneticWrapper magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto autocompletedMagnetic = OpenMagnetics::magnetic_autocomplete(magnetic);

        CHECK(autocompletedMagnetic.get_core().get_geometrical_description());
        CHECK(autocompletedMagnetic.get_core().get_processed_description());
        CHECK(std::holds_alternative<OpenMagnetics::CoreShape>(autocompletedMagnetic.get_core().get_functional_description().get_shape()));
        CHECK(std::holds_alternative<OpenMagnetics::CoreMaterial>(autocompletedMagnetic.get_core().get_functional_description().get_material()));
        CHECK(std::holds_alternative<OpenMagnetics::Bobbin>(autocompletedMagnetic.get_coil().get_bobbin()));
        CHECK(autocompletedMagnetic.get_mutable_coil().resolve_bobbin().get_processed_description());
        CHECK(autocompletedMagnetic.get_coil().get_sections_description());
        CHECK(autocompletedMagnetic.get_coil().get_layers_description());
        CHECK(autocompletedMagnetic.get_coil().get_turns_description());

    }

}