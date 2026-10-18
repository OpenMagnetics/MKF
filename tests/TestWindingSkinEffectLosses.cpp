#include "physical_models/WindingSkinEffectLosses.h"
#include "support/Utils.h"
#include "constructive_models/Core.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

using namespace MAS;
using namespace OpenMagnetics;

namespace {
    double maximumError = 0.01;

    TEST_CASE("Test_Skin_Depth_Wire_Material_Data_20C", "[physical-model][skin-losses][smoke-test]") {
        auto materialData = find_wire_material_by_name("copper");
        double frequency = 123000;
        double temperature = 20;

        auto windingSkinEffectLossesModel = WindingSkinEffectLosses();
        auto skinDepth = windingSkinEffectLossesModel.calculate_skin_depth(materialData, frequency, temperature);
        double expectedSkinDepth = 186.09e-6;
        REQUIRE_THAT(skinDepth, Catch::Matchers::WithinAbs(expectedSkinDepth, expectedSkinDepth * maximumError));
    }

    TEST_CASE("Test_Skin_Depth_Wire_Name_20C", "[physical-model][skin-losses][smoke-test]") {
        double frequency = 123000;
        double temperature = 20;

        auto windingSkinEffectLossesModel = WindingSkinEffectLosses();
        auto skinDepth = windingSkinEffectLossesModel.calculate_skin_depth("copper", frequency, temperature);
        double expectedSkinDepth = 186.09e-6;
        REQUIRE_THAT(skinDepth, Catch::Matchers::WithinAbs(expectedSkinDepth, expectedSkinDepth * maximumError));
    }

    TEST_CASE("Test_Skin_Depth_Wire_Material_Data_200C", "[physical-model][skin-losses][smoke-test]") {
        auto materialData = find_wire_material_by_name("copper");
        double frequency = 123000;
        double temperature = 120;

        auto windingSkinEffectLossesModel = WindingSkinEffectLosses();
        auto skinDepth = windingSkinEffectLossesModel.calculate_skin_depth(materialData, frequency, temperature);
        double expectedSkinDepth = 220e-6;
        REQUIRE_THAT(skinDepth, Catch::Matchers::WithinAbs(expectedSkinDepth, expectedSkinDepth * maximumError));
    }
    TEST_CASE("Test_Skin_Factor_Cache_Tells_Round_From_Rectangular", "[physical-model][skin-losses][abt-1498]") {
        // ABT #1498: the skin-factor cache keyed an unnamed wire by (conductors, outer width, outer height), so
        // a 1 mm round wire and a 1 x 1 mm rectangular wire were the same wire: whichever ran first set the
        // other's skin factor.
        auto make = [](WireType type) {
            OpenMagnetics::Wire wire;
            wire.set_type(type);
            wire.set_number_conductors(1);
            wire.set_material("copper");
            if (type == WireType::ROUND) {
                DimensionWithTolerance diameter;
                diameter.set_nominal(1e-3);
                wire.set_conducting_diameter(diameter);
                wire.set_outer_diameter(diameter);
            }
            else {
                // Same 1 x 1 mm outer envelope as the round wire (what the old key saw), but a 1 x 0.5 mm
                // copper section: ALBACH reads the thinner conducting side, so the factors really differ.
                DimensionWithTolerance side;
                side.set_nominal(1e-3);
                DimensionWithTolerance thin;
                thin.set_nominal(0.5e-3);
                wire.set_conducting_width(side);
                wire.set_conducting_height(thin);
                wire.set_outer_width(side);
                wire.set_outer_height(side);
            }
            return wire;
        };
        const double frequency = 1e6;  // skin depth ~66 um: every section below is several depths thick
        const double temperature = 25;
        const double dcLoss = 1;
        auto freshRectangular = WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::ALBACH);
        const double rectangularAlone = freshRectangular->calculate_turn_losses(make(WireType::RECTANGULAR), dcLoss, frequency, temperature);
        auto freshRound = WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::ALBACH);
        const double roundAlone = freshRound->calculate_turn_losses(make(WireType::ROUND), dcLoss, frequency, temperature);
        // the two wires must genuinely differ for this test to mean anything
        REQUIRE(std::abs(rectangularAlone / roundAlone - 1) > 0.01);

        auto shared = WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::ALBACH);
        CHECK(shared->calculate_turn_losses(make(WireType::ROUND), dcLoss, frequency, temperature) == Catch::Approx(roundAlone).epsilon(1e-12));
        CHECK(shared->calculate_turn_losses(make(WireType::RECTANGULAR), dcLoss, frequency, temperature) == Catch::Approx(rectangularAlone).epsilon(1e-12));

        // A named wire was keyed by its name alone: the same name re-cut to another thickness, or carrying
        // an inline material of another resistivity, read the first one's factor.
        auto named = [&](double width, double resistivity) {
            auto wire = make(WireType::RECTANGULAR);
            wire.set_name("Trace under test");
            DimensionWithTolerance side;
            side.set_nominal(width);
            wire.set_conducting_height(side);
            wire.set_outer_height(side);
            auto material = find_wire_material_by_name("copper");
            material.get_mutable_resistivity().set_reference_value(resistivity);
            wire.set_material(material);
            return wire;
        };
        const double copperResistivity = find_wire_material_by_name("copper").get_resistivity().get_reference_value();
        const double thinAlone = WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::ALBACH)->calculate_turn_losses(named(0.3e-3, copperResistivity), dcLoss, frequency, temperature);
        const double resistiveAlone = WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::ALBACH)->calculate_turn_losses(named(0.5e-3, 4 * copperResistivity), dcLoss, frequency, temperature);
        const double narrowAlone = WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::ALBACH)->calculate_turn_losses(named(0.5e-3, copperResistivity), dcLoss, frequency, temperature);
        REQUIRE(std::abs(thinAlone / narrowAlone - 1) > 0.01);
        REQUIRE(std::abs(resistiveAlone / narrowAlone - 1) > 0.01);
        auto sharedNamed = WindingSkinEffectLossesModel::factory(WindingSkinEffectLossesModels::ALBACH);
        CHECK(sharedNamed->calculate_turn_losses(named(0.5e-3, copperResistivity), dcLoss, frequency, temperature) == Catch::Approx(narrowAlone).epsilon(1e-12));
        CHECK(sharedNamed->calculate_turn_losses(named(0.3e-3, copperResistivity), dcLoss, frequency, temperature) == Catch::Approx(thinAlone).epsilon(1e-12));
        CHECK(sharedNamed->calculate_turn_losses(named(0.5e-3, 4 * copperResistivity), dcLoss, frequency, temperature) == Catch::Approx(resistiveAlone).epsilon(1e-12));
    }

// End of SUITE

}  // namespace
