#include "Utils.h"
#include "Settings.h"
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
        auto magneticsShapeNames = OpenMagnetics::get_shape_names("Magnetics");
        auto ferroxcubeShapeNames = OpenMagnetics::get_shape_names("Ferroxcube");

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

}