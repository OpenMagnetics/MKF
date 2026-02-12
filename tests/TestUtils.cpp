#include <source_location>
#include "support/Painter.h"
#include "support/Utils.h"
#include "support/Settings.h"
#include "TestingUtils.h"
#include "json.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <cfloat>
#include <limits>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
using json = nlohmann::json;
#include <typeinfo>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    TEST_CASE("CeilFloat", "[support][utils][smoke-test]") {
        double value = 1.263;
        double calculatedValue = ceilFloat(value, 2);
        double expectedValue = 1.27;
        REQUIRE(expectedValue == calculatedValue);
    }

    TEST_CASE("FloorFloat", "[support][utils][smoke-test]") {
        double value = 1.263;
        double calculatedValue = floorFloat(value, 2);
        double expectedValue = 1.26;
        REQUIRE(expectedValue == calculatedValue);
    }

    TEST_CASE("Modified_Bessel", "[support][utils][smoke-test]") {
        double calculatedValue = modified_bessel_first_kind(0.0, std::complex<double>{1.0, 0.0}).real();
        double expectedValue = 1.2660658777520084;
        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(calculatedValue, expectedValue * 0.001));
    }

    TEST_CASE("Bessel", "[support][utils][smoke-test]") {
        double calculatedValue = bessel_first_kind(0.0, std::complex<double>{1.0, 0.0}).real();
        double expectedValue = 0.7651976865579666;
        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(calculatedValue, expectedValue * 0.001));

        double calculatedBerValue = kelvin_function_real(0.0, 1.0);
        double expectedBerValue = 0.98438178;
        REQUIRE_THAT(expectedBerValue, Catch::Matchers::WithinAbs(calculatedBerValue, expectedBerValue * 0.001));

        double calculatedBeiValue = kelvin_function_imaginary(0.0, 1.0);
        double expectedBeiValue = 0.24956604;
        REQUIRE_THAT(expectedBeiValue, Catch::Matchers::WithinAbs(calculatedBeiValue, expectedBeiValue * 0.001));

        double calculatedBerpValue = derivative_kelvin_function_real(0.0, 1.0);
        double expectedBerpValue = -0.06244575217903096;
        REQUIRE_THAT(expectedBerpValue, Catch::Matchers::WithinAbs(calculatedBerpValue, fabs(expectedBerpValue) * 0.001));

        double calculatedBeipValue = derivative_kelvin_function_imaginary(0.0, 1.0);
        double expectedBeipValue = 0.49739651146809727;
        REQUIRE_THAT(expectedBeipValue, Catch::Matchers::WithinAbs(calculatedBeipValue, expectedBeipValue * 0.001));
    }

    TEST_CASE("Test_Complete_Ellipitical_1_0", "[support][utils][smoke-test]") {
        double calculatedValue = comp_ellint_1(0);
        double expectedValue = M_PI / 2.0;
        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(calculatedValue, expectedValue * 0.001));
    }

    TEST_CASE("Test_Complete_Ellipitical_1_1", "[support][utils][smoke-test]") {
        double calculatedValue = comp_ellint_1(1);
        double expectedValue = std::numeric_limits<double>::infinity();
        // Our implementation returns NaN, standard library may return NaN or infinity
        REQUIRE(std::isnan(calculatedValue));
        // MSVC returns infinity, GCC/Clang return NaN - both are valid for the singularity at k=1
        REQUIRE((std::isnan(expectedValue) || std::isinf(expectedValue)));
    }

    TEST_CASE("Test_Complete_Ellipitical_1_2", "[support][utils][smoke-test]") {
        double calculatedValue = comp_ellint_1(std::sin(M_PI / 18 / 2));
        double expectedValue = 1.685058;
        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(calculatedValue, expectedValue * 0.001));
    }

    TEST_CASE("Test_Complete_Ellipitical_2_0", "[support][utils][smoke-test]") {
        double calculatedValue = comp_ellint_2(0);
        double expectedValue = M_PI / 2.0;
        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(calculatedValue, expectedValue * 0.001));
    }

    TEST_CASE("Test_Complete_Ellipitical_2_1", "[support][utils][smoke-test]") {
        double calculatedValue = comp_ellint_2(1);
        double expectedValue = 1.0;
        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(calculatedValue, expectedValue * 0.001));
    }

    TEST_CASE("Test_Complete_Ellipitical_2_2", "[support][utils][smoke-test]") {
        double calculatedValue = comp_ellint_2(std::sin(M_PI / 18 / 2));
        double expectedValue = 1.576;
        REQUIRE_THAT(expectedValue, Catch::Matchers::WithinAbs(calculatedValue, expectedValue * 0.001));
    }

    TEST_CASE("Test_Find_By_Area_Product", "[support][utils][smoke-test]") {
        clear_databases();
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);

        auto shape = find_core_shape_by_area_product(0.00000003568, MAS::CoreShapeFamily::PQ);

        REQUIRE("PQ 35/35" == shape.get_name().value());
    }

    TEST_CASE("Test_Find_By_Winding_Window_Area", "[support][utils][smoke-test]") {
        clear_databases();
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);

        auto shape = find_core_shape_by_winding_window_area(0.00022062, MAS::CoreShapeFamily::PQ);

        REQUIRE("PQ 35/35" == shape.get_name().value());
    }

    TEST_CASE("Test_Find_By_Winding_Window_Dimensions", "[support][utils][smoke-test]") {
        clear_databases();
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);
        auto shape = find_core_shape_by_winding_window_dimensions(0.008825, 0.025, MAS::CoreShapeFamily::PQ);

        REQUIRE("PQ 35/35" == shape.get_name().value());
    }

    TEST_CASE("Test_Find_By_Effective_Parameters", "[support][utils][smoke-test]") {
        clear_databases();
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);
        auto shape = find_core_shape_by_effective_parameters(0.079, 0.000171, 0.000014, MAS::CoreShapeFamily::PQ);

        REQUIRE("PQ 35/35" == shape.get_name().value());
    }

    TEST_CASE("Test_Find_By_Perimeter", "[support][utils][smoke-test]") {
        clear_databases();
        settings.set_use_toroidal_cores(true);
        settings.set_use_concentric_cores(true);

        auto shape = find_core_shape_by_winding_window_perimeter(0.03487);

        REQUIRE("UR 46/21/11" == shape.get_name().value());
    }

    TEST_CASE("Test_Find_By_Perimeter_Only_Toroids", "[support][utils][smoke-test]") {
        clear_databases();
        settings.set_use_toroidal_cores(true);
        settings.set_use_concentric_cores(false);

        auto shape = find_core_shape_by_winding_window_perimeter(0.03487);

        REQUIRE("T 22/12.4/12.8" == shape.get_name().value());
    }

    TEST_CASE("Test_Get_Shapes", "[support][utils][smoke-test]") {
        clear_databases();
        settings.reset();
        auto allShapeNames = get_core_shape_names();

        REQUIRE(allShapeNames.size() > 900);
    }

    TEST_CASE("Test_Get_Shapes_By_Manufacturer", "[support][utils][smoke-test]") {
        clear_databases();
        settings.reset();
        auto allShapeNames = get_core_shape_names();
        auto magneticsShapeNames = get_core_shape_names("Magnetics");
        auto ferroxcubeShapeNames = get_core_shape_names("Ferroxcube");

        REQUIRE(allShapeNames.size() > magneticsShapeNames.size());
        REQUIRE(allShapeNames.size() > ferroxcubeShapeNames.size());
    }

    TEST_CASE("Test_Wire_Names_With_Types", "[support][utils][smoke-test]") {
        clear_databases();
        settings.reset();
        auto allWireNames = get_wire_names();
        for (auto wire : allWireNames) {
            REQUIRE((
                wire.starts_with("Round ") ||
                wire.starts_with("Litz ") ||
                wire.starts_with("Rectangular ") ||
                wire.starts_with("Planar ") ||
                wire.starts_with("Foil ")
                ));
        }
    }

    TEST_CASE("Test_Wires_With_Type_And_Standard", "[support][utils][smoke-test]") {
        clear_databases();
        settings.reset();
        auto allWires = get_wires();
        auto someWires = get_wires(WireType::ROUND, WireStandard::IEC_60317);
        REQUIRE(someWires.size() < allWires.size());
    }

    TEST_CASE("Test_Load_Toroidal_Cores", "[support][utils][smoke-test]") {
        clear_databases();
        settings.reset();
        settings.set_use_toroidal_cores(true);
        settings.set_use_concentric_cores(false);
        load_cores();

        auto allCores = coreDatabase;
        REQUIRE(allCores.size() > 0);
        for (auto core : allCores) {
            REQUIRE(core.get_type() == CoreType::TOROIDAL);
        }
    }

    TEST_CASE("Test_Load_Two_Piece_Set_Cores", "[support][utils][smoke-test]") {
        clear_databases();
        settings.reset();
        settings.set_use_toroidal_cores(false);
        settings.set_use_concentric_cores(true);
        load_cores();

        auto allCores = coreDatabase;
        REQUIRE(allCores.size() > 0);
        for (auto core : allCores) {
            REQUIRE(core.get_type() == CoreType::TWO_PIECE_SET);
        }
    }

    TEST_CASE("Test_Load_Cores_In_Stock", "[support][utils][smoke-test]") {
        clear_databases();
        settings.reset();
        settings.set_use_only_cores_in_stock(false);
        load_cores();

        auto allCores = coreDatabase;

        clear_databases();
        settings.reset();
        settings.set_use_only_cores_in_stock(true);
        load_cores();

        auto onlyCoresInStock = coreDatabase;
        REQUIRE(allCores.size() > onlyCoresInStock.size());
    }

    TEST_CASE("Test_Core_Materials_External", "[support][utils][smoke-test]") {
        auto external_core_materials_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "external_core_materials.ndjson");

        std::ifstream file(external_core_materials_path, std::ios_base::binary | std::ios_base::in);
        if(!file.is_open())
            throw std::runtime_error("Failed to open " + external_core_materials_path.string());
        using Iterator = std::istreambuf_iterator<char>;
        std::string external_core_materials(Iterator{file}, Iterator{});

        clear_databases();
        load_core_materials(external_core_materials);

        auto allCoreMaterials = coreMaterialDatabase;
        REQUIRE(allCoreMaterials.size() == 4);
    }

    TEST_CASE("Test_Core_Materials_External_And_Internal", "[support][utils][smoke-test]") {
        auto external_core_materials_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "external_core_materials.ndjson");

        std::ifstream file(external_core_materials_path, std::ios_base::binary | std::ios_base::in);
        if(!file.is_open())
            throw std::runtime_error("Failed to open " + external_core_materials_path.string());
        using Iterator = std::istreambuf_iterator<char>;
        std::string external_core_materials(Iterator{file}, Iterator{});

        clear_databases();
        load_core_materials();

        auto allCoreMaterials = coreMaterialDatabase;

        load_core_materials(external_core_materials);

        auto allCoreMaterialsWithExternal = coreMaterialDatabase;

        REQUIRE(allCoreMaterialsWithExternal.size() > allCoreMaterials.size());
    }

    TEST_CASE("Test_Core_Shapes_Families", "[support][utils][smoke-test]") {
        REQUIRE(get_core_shape_families().size() > 0);
    }

    TEST_CASE("Test_Core_Shapes_Dimensions_PQ", "[support][utils][smoke-test]") {
        auto dimensions = get_shape_family_dimensions(CoreShapeFamily::PQ);
        REQUIRE(dimensions.size() == 9);
        REQUIRE(dimensions[0] == "A");
        REQUIRE(dimensions[1] == "B");
        REQUIRE(dimensions[2] == "C");
        REQUIRE(dimensions[3] == "D");
        REQUIRE(dimensions[4] == "E");
        REQUIRE(dimensions[5] == "F");
        REQUIRE(dimensions[6] == "G");
        REQUIRE(dimensions[7] == "J");
        REQUIRE(dimensions[8] == "L");
    }

    TEST_CASE("Test_Core_Shapes_Dimensions_UR", "[support][utils][smoke-test]") {
        auto dimensions = get_shape_family_dimensions(CoreShapeFamily::UR, "2");
        REQUIRE(dimensions.size() == 7);
        REQUIRE(dimensions[0] == "A");
        REQUIRE(dimensions[1] == "B");
        REQUIRE(dimensions[2] == "C");
        REQUIRE(dimensions[3] == "D");
        REQUIRE(dimensions[4] == "E");
        REQUIRE(dimensions[5] == "H");
        REQUIRE(dimensions[6] == "S");
    }

    TEST_CASE("Test_Core_Shapes_Dimensions_PM", "[support][utils][smoke-test]") {
        auto dimensions = get_shape_family_dimensions(CoreShapeFamily::PM, "1");
        REQUIRE(dimensions.size() == 12);
        REQUIRE(dimensions[0] == "A");
        REQUIRE(dimensions[1] == "B");
        REQUIRE(dimensions[2] == "C");
        REQUIRE(dimensions[3] == "D");
        REQUIRE(dimensions[4] == "E");
        REQUIRE(dimensions[5] == "F");
        REQUIRE(dimensions[6] == "G");
        REQUIRE(dimensions[7] == "H");
        REQUIRE(dimensions[8] == "alpha");
        REQUIRE(dimensions[9] == "b");
        REQUIRE(dimensions[10] == "e");
        REQUIRE(dimensions[11] == "t");
    }

    TEST_CASE("Test_Core_Shapes_Dimensions_UR_No_Subtype", "[support][utils][smoke-test]") {
        auto dimensions = get_shape_family_dimensions(CoreShapeFamily::UR);
        REQUIRE(dimensions.size() == 9);
        REQUIRE(dimensions[0] == "A");
        REQUIRE(dimensions[1] == "B");
        REQUIRE(dimensions[2] == "C");
        REQUIRE(dimensions[3] == "D");
        REQUIRE(dimensions[4] == "E");
        REQUIRE(dimensions[5] == "F");
        REQUIRE(dimensions[6] == "G");
        REQUIRE(dimensions[7] == "H");
        REQUIRE(dimensions[8] == "S");
    }

    TEST_CASE("Test_Core_Shapes_Family_Subtypes_UR", "[support][utils][smoke-test]") {
        auto subtypes = get_shape_family_subtypes(CoreShapeFamily::UR);
        REQUIRE(subtypes.size() == 4);
        REQUIRE(subtypes[0] == "1");
        REQUIRE(subtypes[1] == "2");
        REQUIRE(subtypes[2] == "3");
        REQUIRE(subtypes[3] == "4");
    }

    TEST_CASE("Test_Core_Shapes_Family_Subtypes_E", "[support][utils][smoke-test]") {
        auto subtypes = get_shape_family_subtypes(CoreShapeFamily::E);
        REQUIRE(subtypes.size() == 0);
    }

    TEST_CASE("Test_Core_Shapes_External", "[support][utils][smoke-test]") {
        auto external_core_shapes_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "external_core_shapes.ndjson");

        std::ifstream file(external_core_shapes_path, std::ios_base::binary | std::ios_base::in);
        if(!file.is_open())
            throw std::runtime_error("Failed to open " + external_core_shapes_path.string());
        using Iterator = std::istreambuf_iterator<char>;
        std::string external_core_shapes(Iterator{file}, Iterator{});
        settings.reset();
        settings.set_use_toroidal_cores(true);
        settings.set_use_concentric_cores(true);

        clear_databases();
        load_core_shapes(false, external_core_shapes);

        auto allCoreShapes = coreShapeDatabase;
        REQUIRE(allCoreShapes.size() == 10);
    }

    TEST_CASE("Test_Core_Shapes_External_And_Internal", "[support][utils][smoke-test]") {
        auto external_core_shapes_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "external_core_shapes.ndjson");

        std::ifstream file(external_core_shapes_path, std::ios_base::binary | std::ios_base::in);
        if(!file.is_open())
            throw std::runtime_error("Failed to open " + external_core_shapes_path.string());
        using Iterator = std::istreambuf_iterator<char>;
        std::string external_core_shapes(Iterator{file}, Iterator{});
        settings.reset();
        settings.set_use_toroidal_cores(true);
        settings.set_use_concentric_cores(true);

        clear_databases();
        load_core_shapes(false);

        auto allCoreShapes = coreShapeDatabase;

        load_core_shapes(false, external_core_shapes);

        auto allCoreShapesWithExternal = coreShapeDatabase;

        REQUIRE(allCoreShapesWithExternal.size() > allCoreShapes.size());
    }

    TEST_CASE("Test_Wires_External", "[support][utils][smoke-test]") {
        auto external_wires_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "external_wires.ndjson");

        std::ifstream file(external_wires_path, std::ios_base::binary | std::ios_base::in);
        if(!file.is_open())
            throw std::runtime_error("Failed to open " + external_wires_path.string());
        using Iterator = std::istreambuf_iterator<char>;
        std::string external_wires(Iterator{file}, Iterator{});

        clear_databases();
        load_wires(external_wires);

        auto allCoreWires = wireDatabase;
        REQUIRE(allCoreWires.size() == 24);
    }

    TEST_CASE("Test_Wires_External_And_Internal", "[support][utils][smoke-test]") {
        auto external_wires_path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "external_wires.ndjson");

        std::ifstream file(external_wires_path, std::ios_base::binary | std::ios_base::in);
        if(!file.is_open())
            throw std::runtime_error("Failed to open " + external_wires_path.string());
        using Iterator = std::istreambuf_iterator<char>;
        std::string external_wires(Iterator{file}, Iterator{});

        clear_databases();
        load_wires();

        auto allCoreWires = wireDatabase;

        load_wires(external_wires);

        auto allCoreWiresWithExternal = wireDatabase;

        REQUIRE(allCoreWiresWithExternal.size() > allCoreWires.size());
    }

    TEST_CASE("Test_Low_And_High_Harmonics", "[support][utils][smoke-test]") {
        Harmonics harmonics = json::parse(R"({"amplitudes":[2.8679315866586563e-15,2.3315003998761155,1.1284261163900786],"frequencies":[0,50,400000000]})");
        auto mainHarmonicIndexes = get_main_harmonic_indexes(harmonics, 0.05, 1);

        REQUIRE(mainHarmonicIndexes.size() == 2);
    }

    TEST_CASE("Test_Mas_Autocomplete", "[support][utils][smoke-test]") {
        settings.set_use_only_cores_in_stock(false);
        auto core = find_core_by_name("PQ 32/30 - 3C90 - Gapped 0.492 mm");
        OpenMagnetics::Coil coil;
        OpenMagnetics::Winding dummyWinding;
        dummyWinding.set_name("");
        dummyWinding.set_number_turns(0);
        dummyWinding.set_number_parallels(0);
        dummyWinding.set_wire("");
        coil.get_mutable_functional_description().push_back(dummyWinding);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        OpenMagnetics::Inputs inputs = OpenMagnetics::Inputs::create_quick_operating_point(100000, 42e-6, 42, WaveformLabel::SINUSOIDAL, 10, 0.5, 2);

        OpenMagnetics::Mas mas;
        mas.set_inputs(inputs);
        mas.set_magnetic(magnetic);

        auto autocompletedMas = mas_autocomplete(mas);

        REQUIRE(autocompletedMas.get_magnetic().get_core().get_geometrical_description());
        REQUIRE(autocompletedMas.get_magnetic().get_core().get_processed_description());
        REQUIRE(std::holds_alternative<CoreShape>(autocompletedMas.get_magnetic().get_core().get_functional_description().get_shape()));
        REQUIRE(std::holds_alternative<CoreMaterial>(autocompletedMas.get_magnetic().get_core().get_functional_description().get_material()));
        REQUIRE(std::holds_alternative<OpenMagnetics::Bobbin>(autocompletedMas.get_magnetic().get_coil().get_bobbin()));
        REQUIRE(autocompletedMas.get_mutable_magnetic().get_mutable_coil().resolve_bobbin().get_processed_description());
        REQUIRE(autocompletedMas.get_magnetic().get_coil().get_sections_description());
        REQUIRE(autocompletedMas.get_magnetic().get_coil().get_layers_description());
        REQUIRE(autocompletedMas.get_magnetic().get_coil().get_turns_description());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_rms());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_harmonics());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_processed());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()->get_rms());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_waveform());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_harmonics());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_rms());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_harmonics());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points().size() == autocompletedMas.get_outputs().size());
        REQUIRE(autocompletedMas.get_outputs()[0].get_core_losses());
        REQUIRE(autocompletedMas.get_outputs()[0].get_winding_losses());
        REQUIRE(autocompletedMas.get_outputs()[0].get_inductance());
    }

    TEST_CASE("Test_Mas_Autocomplete_Json", "[support][utils][bug][smoke-test]") {
        std::string masString = R"({"outputs": [], "inputs": {"designRequirements": {"isolationSides": ["primary" ], "magnetizingInductance": {"nominal": 0.00039999999999999996 }, "name": "My Design Requirements", "turnsRatios": [{"nominal": 1} ] }, "operatingPoints": [{"conditions": {"ambientTemperature": 42 }, "excitationsPerWinding": [{"frequency": 100000, "current": {"processed": {"label": "Triangular", "peakToPeak": 0.5, "offset": 0, "dutyCycle": 0.5 } }, "voltage": {"processed": {"label": "Rectangular", "peakToPeak": 20, "offset": 0, "dutyCycle": 0.5 } } } ], "name": "Operating Point No. 1" } ] }, "magnetic": {"coil": {"bobbin": "Basic", "functionalDescription":[{"name": "Primary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "primary", "wire": "Round 1.00 - Grade 1" }, {"name": "Secondary", "numberTurns": 4, "numberParallels": 1, "isolationSide": "secondary", "wire": "Round 1.00 - Grade 1" } ] }, "core": {"name": "core_E_19_8_5_N87_substractive", "functionalDescription": {"type": "two-piece set", "material": "N87", "shape": "PQ 32/20", "gapping": [{"type": "residual", "length": 0.000005 }], "numberStacks": 1 } }, "manufacturerInfo": {"name": "", "reference": "Example" } } })";
        json masJson = json::parse(masString);
        OpenMagnetics::Mas mas(masJson);

        auto autocompletedMas = mas_autocomplete(mas);

        REQUIRE(autocompletedMas.get_magnetic().get_core().get_geometrical_description());
        REQUIRE(autocompletedMas.get_magnetic().get_core().get_processed_description());
        REQUIRE(std::holds_alternative<CoreShape>(autocompletedMas.get_magnetic().get_core().get_functional_description().get_shape()));
        REQUIRE(std::holds_alternative<CoreMaterial>(autocompletedMas.get_magnetic().get_core().get_functional_description().get_material()));
        REQUIRE(std::holds_alternative<OpenMagnetics::Bobbin>(autocompletedMas.get_magnetic().get_coil().get_bobbin()));
        REQUIRE(autocompletedMas.get_mutable_magnetic().get_mutable_coil().resolve_bobbin().get_processed_description());
        REQUIRE(autocompletedMas.get_magnetic().get_coil().get_sections_description());
        REQUIRE(autocompletedMas.get_magnetic().get_coil().get_layers_description());
        REQUIRE(autocompletedMas.get_magnetic().get_coil().get_turns_description());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_rms());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_harmonics());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_processed());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()->get_rms());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_waveform());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_harmonics());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_rms());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_harmonics());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points().size() == autocompletedMas.get_outputs().size());
        REQUIRE(autocompletedMas.get_outputs()[0].get_core_losses());
        REQUIRE(autocompletedMas.get_outputs()[0].get_winding_losses());
        REQUIRE(autocompletedMas.get_outputs()[0].get_inductance());
    }

    TEST_CASE("Test_Mas_Autocomplete_Json_2", "[support][utils][bug][smoke-test]") {
        auto json_path_548 = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "test_mas_autocomplete_json_2_548.json");
        std::ifstream json_file_548(json_path_548);
        std::string masString((std::istreambuf_iterator<char>(json_file_548)), std::istreambuf_iterator<char>());
        json masJson = json::parse(masString);
        OpenMagnetics::Mas mas(masJson);

        auto autocompletedMas = mas_autocomplete(mas);

        REQUIRE(autocompletedMas.get_magnetic().get_core().get_geometrical_description());
        REQUIRE(autocompletedMas.get_magnetic().get_core().get_processed_description());
        REQUIRE(std::holds_alternative<CoreShape>(autocompletedMas.get_magnetic().get_core().get_functional_description().get_shape()));
        REQUIRE(std::holds_alternative<CoreMaterial>(autocompletedMas.get_magnetic().get_core().get_functional_description().get_material()));
        REQUIRE(std::holds_alternative<OpenMagnetics::Bobbin>(autocompletedMas.get_magnetic().get_coil().get_bobbin()));
        REQUIRE(autocompletedMas.get_mutable_magnetic().get_mutable_coil().resolve_bobbin().get_processed_description());
        REQUIRE(autocompletedMas.get_magnetic().get_coil().get_sections_description());
        REQUIRE(autocompletedMas.get_magnetic().get_coil().get_layers_description());
        REQUIRE(autocompletedMas.get_magnetic().get_coil().get_turns_description());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_processed()->get_rms());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_waveform());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_current()->get_harmonics());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_processed());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_processed()->get_rms());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_waveform());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_magnetizing_current()->get_harmonics());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_processed()->get_rms());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_waveform());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points()[0].get_excitations_per_winding()[0].get_voltage()->get_harmonics());
        REQUIRE(autocompletedMas.get_inputs().get_operating_points().size() == autocompletedMas.get_outputs().size());
        REQUIRE(autocompletedMas.get_outputs()[0].get_core_losses());
        REQUIRE(autocompletedMas.get_outputs()[0].get_winding_losses());
        REQUIRE(autocompletedMas.get_outputs()[0].get_inductance());
    }

    TEST_CASE("Test_Mas_Autocomplete_Only_Magnetic", "[support][utils][smoke-test]") {
        settings.set_use_only_cores_in_stock(false);
        auto core = find_core_by_name("PQ 32/30 - 3C90 - Gapped 0.492 mm");
        OpenMagnetics::Coil coil;
        OpenMagnetics::Winding dummyWinding;
        dummyWinding.set_name("");
        dummyWinding.set_number_turns(0);
        dummyWinding.set_number_parallels(0);
        dummyWinding.set_wire("");
        coil.get_mutable_functional_description().push_back(dummyWinding);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto autocompletedMagnetic = magnetic_autocomplete(magnetic);

        REQUIRE(autocompletedMagnetic.get_core().get_geometrical_description());
        REQUIRE(autocompletedMagnetic.get_core().get_processed_description());
        REQUIRE(std::holds_alternative<CoreShape>(autocompletedMagnetic.get_core().get_functional_description().get_shape()));
        REQUIRE(std::holds_alternative<CoreMaterial>(autocompletedMagnetic.get_core().get_functional_description().get_material()));
        REQUIRE(std::holds_alternative<OpenMagnetics::Bobbin>(autocompletedMagnetic.get_coil().get_bobbin()));
        REQUIRE(autocompletedMagnetic.get_mutable_coil().resolve_bobbin().get_processed_description());
        REQUIRE(autocompletedMagnetic.get_coil().get_sections_description());
        REQUIRE(autocompletedMagnetic.get_coil().get_layers_description());
        REQUIRE(autocompletedMagnetic.get_coil().get_turns_description());
    }

    TEST_CASE("Test_Mas_Autocomplete_Only_Magnetic_Edge_Wound", "[support][utils][smoke-test]") {
        settings.set_use_only_cores_in_stock(false);
        auto core = find_core_by_name("PQ 32/30 - 3C90 - Gapped 0.492 mm");
        OpenMagnetics::Coil coil;
        {
            OpenMagnetics::Winding winding;
            winding.set_name("");
            winding.set_number_turns(12);
            winding.set_number_parallels(1);
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0048);
            wire.set_nominal_value_conducting_height(0.00076);
            wire.set_nominal_value_outer_width(0.0049);
            wire.set_nominal_value_outer_height(0.0007676);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::RECTANGULAR);
            winding.set_wire(wire);
            coil.get_mutable_functional_description().push_back(winding);
        }
        {
            OpenMagnetics::Winding winding;
            winding.set_name("");
            winding.set_number_turns(1);
            winding.set_number_parallels(12);
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0038);
            wire.set_nominal_value_conducting_height(0.00076);
            wire.set_nominal_value_outer_width(0.0039);
            wire.set_nominal_value_outer_height(0.0007676);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::RECTANGULAR);
            winding.set_wire(wire);
            coil.get_mutable_functional_description().push_back(winding);
        }
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        auto autocompletedMagnetic = magnetic_autocomplete(magnetic);

        REQUIRE(autocompletedMagnetic.get_core().get_geometrical_description());
        REQUIRE(autocompletedMagnetic.get_core().get_processed_description());
        REQUIRE(std::holds_alternative<CoreShape>(autocompletedMagnetic.get_core().get_functional_description().get_shape()));
        REQUIRE(std::holds_alternative<CoreMaterial>(autocompletedMagnetic.get_core().get_functional_description().get_material()));
        REQUIRE(std::holds_alternative<OpenMagnetics::Bobbin>(autocompletedMagnetic.get_coil().get_bobbin()));
        REQUIRE(autocompletedMagnetic.get_mutable_coil().resolve_bobbin().get_processed_description());
        REQUIRE(autocompletedMagnetic.get_coil().get_sections_description());
        REQUIRE(autocompletedMagnetic.get_coil().get_layers_description());
        REQUIRE(autocompletedMagnetic.get_coil().get_turns_description());


        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Mas_Autocomplete_Only_Magnetic_Edge_Wound.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);

        painter.paint_core(autocompletedMagnetic);
        painter.paint_coil_turns(autocompletedMagnetic);
        painter.export_svg();
    }

    TEST_CASE("Test_Mas_Autocomplete_Only_Magnetic_Edge_Wound_Configured", "[support][utils][smoke-test]") {
        settings.set_use_only_cores_in_stock(false);
        auto core = find_core_by_name("PQ 32/30 - 3C90 - Gapped 0.492 mm");
        OpenMagnetics::Coil coil;
        {
            OpenMagnetics::Winding winding;
            winding.set_name("");
            winding.set_number_turns(12);
            winding.set_number_parallels(1);
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0048);
            wire.set_nominal_value_conducting_height(0.00056);
            wire.set_nominal_value_outer_width(0.0049);
            wire.set_nominal_value_outer_height(0.0005676);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::RECTANGULAR);
            winding.set_wire(wire);
            coil.get_mutable_functional_description().push_back(winding);
        }
        {
            OpenMagnetics::Winding winding;
            winding.set_name("");
            winding.set_number_turns(1);
            winding.set_number_parallels(12);
            OpenMagnetics::Wire wire;
            wire.set_nominal_value_conducting_width(0.0038);
            wire.set_nominal_value_conducting_height(0.00056);
            wire.set_nominal_value_outer_width(0.0039);
            wire.set_nominal_value_outer_height(0.0005676);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            wire.set_type(WireType::RECTANGULAR);
            winding.set_wire(wire);
            coil.get_mutable_functional_description().push_back(winding);
        }
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_coil(coil);

        json configuration;

        configuration["windingOrientation"] = "contiguous";
        configuration["sectionAlignment"] = "spread";
        configuration["interleavingPattern"] = std::vector<size_t>{0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
        configuration["turnsAlignment"] = "spread";

        auto autocompletedMagnetic = magnetic_autocomplete(magnetic, configuration);

        REQUIRE(autocompletedMagnetic.get_core().get_geometrical_description());
        REQUIRE(autocompletedMagnetic.get_core().get_processed_description());
        REQUIRE(std::holds_alternative<CoreShape>(autocompletedMagnetic.get_core().get_functional_description().get_shape()));
        REQUIRE(std::holds_alternative<CoreMaterial>(autocompletedMagnetic.get_core().get_functional_description().get_material()));
        REQUIRE(std::holds_alternative<OpenMagnetics::Bobbin>(autocompletedMagnetic.get_coil().get_bobbin()));
        REQUIRE(autocompletedMagnetic.get_mutable_coil().resolve_bobbin().get_processed_description());
        REQUIRE(autocompletedMagnetic.get_coil().get_sections_description());
        REQUIRE(autocompletedMagnetic.get_coil().get_layers_description());
        REQUIRE(autocompletedMagnetic.get_coil().get_turns_description());


        auto outputFilePath = std::filesystem::path{ std::source_location::current().file_name() }.parent_path().append("..").append("output");
        auto outFile = outputFilePath;
        outFile.append("Test_Mas_Autocomplete_Only_Magnetic_Edge_Wound_Configured.svg");
        std::filesystem::remove(outFile);
        Painter painter(outFile);

        painter.paint_core(autocompletedMagnetic);
        painter.paint_coil_turns(autocompletedMagnetic);
        painter.export_svg();
    }

    TEST_CASE("Test_Mas_Autocomplete_Web_0", "[support][utils][bug][smoke-test]") {
        auto path = OpenMagneticsTesting::get_test_data_path(std::source_location::current(), "geometricalDescriptionDebug.mas");
        auto mas = OpenMagneticsTesting::mas_loader(path.string());
        auto magnetic = mas.get_magnetic();

        auto autocompletedMagnetic = magnetic_autocomplete(magnetic);

        REQUIRE(autocompletedMagnetic.get_core().get_geometrical_description().value()[0].get_machining());
        REQUIRE(autocompletedMagnetic.get_core().get_geometrical_description().value()[0].get_machining().value()[0].get_length() > 0);
        REQUIRE(!autocompletedMagnetic.get_core().get_geometrical_description().value()[1].get_machining());
    }

}  // namespace
