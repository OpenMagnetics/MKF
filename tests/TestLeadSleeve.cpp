// ABT #1174 (WP5, MAS-RFC 0016): lead sleeving.
//
// The insulation coordinator decides whether a terminal lead needs a sleeve and which one, the coil
// records it on connections[].sleeve, reserves the sleeved lead in the margin band (the one lead the
// ABT #684 keep-out lets in), reports per-lead creepage, and the DFM rule R11 judges the chosen
// protection against the requirement.
#include "advisers/Manufacturability.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/Mas.h"
#include "support/Settings.h"
#include "support/Utils.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>
#include <source_location>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

const std::string enamelledWireName = "Round 0.3 - Grade 2";   // enamel, grade 2, OD max 0.352 mm
const std::string tripleInsulatedWireName = "Round TCA3 26 AWG";   // 3 extruded layers, 6 kV, 155 C

OpenMagnetics::Inputs make_reinforced_inputs(std::vector<InsulationStandards> standards, double ambientTemperature) {
    DimensionWithTolerance altitude;
    altitude.set_maximum(2000);
    DimensionWithTolerance mainSupplyVoltage;
    mainSupplyVoltage.set_maximum(250);
    auto inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, Cti::GROUP_II, IsolationClass::REINFORCED,
                                                                   mainSupplyVoltage, OvervoltageCategory::II,
                                                                   PollutionDegree::PD2, standards, 400, 600, 100000,
                                                                   WiringTechnology::WOUND);
    inputs.get_mutable_operating_points()[0].get_mutable_conditions().set_ambient_temperature(ambientTemperature);
    return inputs;
}

std::string sleeve_material_name(const ConnectionSleeve& sleeve) {
    REQUIRE(std::holds_alternative<std::string>(sleeve.get_material()));
    return std::get<std::string>(sleeve.get_material());
}

std::filesystem::path examples_dir() {
    return std::filesystem::path{std::source_location::current().file_name()}.parent_path().append("..").append("MAS").append("examples");
}

// The margin-wound flyback of the MAS examples: PQ 32/30, reinforced, IEC 60664-1 + IEC 62368-1,
// three enamelled windings. `edit` changes the document before it is autocompleted.
OpenMagnetics::Mas autocomplete_flyback(const std::function<void(json&)>& edit = [](json&) {}) {
    std::ifstream file((examples_dir() / "24_margin_interleaved_flyback_pq3230_3c94.json").string());
    json masJson = json::parse(file);
    edit(masJson);
    auto path = std::filesystem::temp_directory_path() / "mkf_abt1174_flyback.json";
    {
        std::ofstream out(path.string());
        out << masJson.dump();
    }
    auto mas = OpenMagneticsTesting::mas_loader(path.string());
    std::filesystem::remove(path);
    return OpenMagnetics::mas_autocomplete(mas, false);
}

LeadInsulation find_lead(const InsulationCoordinationResult& result, const std::string& winding, End end) {
    for (auto& lead : result.get_leads()) {
        if (lead.winding == winding && lead.end == end) {
            return lead;
        }
    }
    FAIL("no lead report for " << winding);
    throw std::runtime_error("unreachable");
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// The coordinator's rule
// ---------------------------------------------------------------------------------------------

TEST_CASE("Test_Lead_Sleeve_Enamelled_Lead_Crossing_A_Margin_Is_Sleeved", "[constructive-model][insulation][lead-sleeve]") {
    clear_databases();
    InsulationCoordinator coordinator;
    auto inputs = make_reinforced_inputs({InsulationStandards::IEC_623681}, 55);
    auto wire = find_wire_by_name(enamelledWireName);

    auto sleeve = coordinator.calculate_lead_sleeve_requirements(inputs, wire, true);
    REQUIRE(sleeve);
    // Thinnest wall >= max(0.4 mm, DTI) on offer is 0.41 mm, sold both in PTFE (260 C) and polyolefin
    // (135 C); the tie goes to the lowest class that still covers the 55 C ambient.
    REQUIRE(coordinator.calculate_distance_through_insulation(inputs) <= 0.41e-3);
    CHECK_THAT(sleeve->get_wall_thickness(), Catch::Matchers::WithinAbs(0.41e-3, 1e-12));
    CHECK(sleeve_material_name(sleeve.value()) == "Polyolefin heat-shrink tubing");
    CHECK_THAT(sleeve->get_inner_diameter(), Catch::Matchers::WithinAbs(wire.get_maximum_outer_width() + 0.1e-3, 1e-12));
    REQUIRE(sleeve->get_overlap_into_winding());
    CHECK(sleeve->get_overlap_into_winding().value() == 2e-3);
    REQUIRE(sleeve->get_number_layers());
    CHECK(sleeve->get_number_layers().value() == 1);
}

TEST_CASE("Test_Lead_Sleeve_Temperature_Class_Filters_The_Stock", "[constructive-model][insulation][lead-sleeve]") {
    clear_databases();
    InsulationCoordinator coordinator;
    auto inputs = make_reinforced_inputs({InsulationStandards::IEC_623681}, 150);
    auto wire = find_wire_by_name(enamelledWireName);

    auto sleeve = coordinator.calculate_lead_sleeve_requirements(inputs, wire, true);
    REQUIRE(sleeve);
    // 150 C rules out the 135 C polyolefin: PTFE at its 0.41 mm standard wall.
    CHECK(sleeve_material_name(sleeve.value()) == "PTFE extruded tubing");
    CHECK_THAT(sleeve->get_wall_thickness(), Catch::Matchers::WithinAbs(0.41e-3, 1e-12));
}

TEST_CASE("Test_Lead_Sleeve_Not_Needed_Without_A_Margin_Or_A_Safety_Class", "[constructive-model][insulation][lead-sleeve]") {
    clear_databases();
    InsulationCoordinator coordinator;
    auto wire = find_wire_by_name(enamelledWireName);
    auto inputs = make_reinforced_inputs({InsulationStandards::IEC_623681}, 55);
    CHECK_FALSE(coordinator.calculate_lead_sleeve_requirements(inputs, wire, false));

    auto insulation = inputs.get_design_requirements().get_insulation().value();
    insulation.set_insulation_type(IsolationClass::FUNCTIONAL);
    inputs.get_mutable_design_requirements().set_insulation(insulation);
    CHECK_FALSE(coordinator.calculate_lead_sleeve_requirements(inputs, wire, true));
}

TEST_CASE("Test_Lead_Sleeve_Triple_Insulated_Lead_Needs_None", "[constructive-model][insulation][lead-sleeve]") {
    clear_databases();
    InsulationCoordinator coordinator;
    auto inputs = make_reinforced_inputs({InsulationStandards::IEC_623681}, 55);
    auto wire = find_wire_by_name(tripleInsulatedWireName);
    CHECK_FALSE(coordinator.calculate_lead_sleeve_requirements(inputs, wire, true));
}

TEST_CASE("Test_Lead_Sleeve_IEC_60664_Disables_The_Insulated_Wire_Escape", "[constructive-model][insulation][lead-sleeve]") {
    clear_databases();
    InsulationCoordinator coordinator;
    // A low-voltage mains design, so the triple insulated wire's own layers each withstand the
    // voltage under both standard sets; ONLY the standard's refusal to accept a fully insulated wire
    // can then be what asks for the sleeve.
    auto make_inputs = [](std::vector<InsulationStandards> standards) {
        DimensionWithTolerance altitude;
        altitude.set_maximum(2000);
        DimensionWithTolerance mainSupplyVoltage;
        mainSupplyVoltage.set_maximum(120);
        auto inputs = OpenMagneticsTesting::get_quick_insulation_inputs(altitude, Cti::GROUP_II, IsolationClass::REINFORCED,
                                                                       mainSupplyVoltage, OvervoltageCategory::I,
                                                                       PollutionDegree::PD2, standards, 150, 212, 20000,
                                                                       WiringTechnology::WOUND);
        inputs.get_mutable_operating_points()[0].get_mutable_conditions().set_ambient_temperature(55);
        return inputs;
    };
    auto wire = find_wire_by_name(tripleInsulatedWireName);
    double breakdownVoltage = wire.resolve_coating()->get_breakdown_voltage().value();

    auto inputs62368 = make_inputs({InsulationStandards::IEC_623681});
    INFO("IEC 62368-1 withstand voltage " << coordinator.calculate_withstand_voltage(inputs62368));
    REQUIRE(coordinator.calculate_withstand_voltage(inputs62368) < breakdownVoltage);
    CHECK_FALSE(coordinator.calculate_lead_sleeve_requirements(inputs62368, wire, true));

    auto inputs60664 = make_inputs({InsulationStandards::IEC_623681, InsulationStandards::IEC_606641});
    REQUIRE_FALSE(InsulationCoordinator::can_fully_insulated_wire_be_used(inputs60664));
    INFO("IEC 60664-1 withstand voltage " << coordinator.calculate_withstand_voltage(inputs60664));
    REQUIRE(coordinator.calculate_withstand_voltage(inputs60664) < breakdownVoltage);
    auto sleeve = coordinator.calculate_lead_sleeve_requirements(inputs60664, wire, true);
    REQUIRE(sleeve);
    CHECK(sleeve->get_wall_thickness() >= 0.4e-3);
}

TEST_CASE("Test_Lead_Sleeve_Throws_When_No_Sleeve_Stock_Exists", "[constructive-model][insulation][lead-sleeve]") {
    clear_databases();
    load_insulation_materials();
    for (auto it = insulationMaterialDatabase.begin(); it != insulationMaterialDatabase.end();) {
        if (it->second.get_form() && it->second.get_form().value() == Form::SLEEVE) {
            it = insulationMaterialDatabase.erase(it);
        }
        else {
            ++it;
        }
    }
    InsulationCoordinator coordinator;
    auto inputs = make_reinforced_inputs({InsulationStandards::IEC_623681}, 55);
    auto wire = find_wire_by_name(enamelledWireName);
    CHECK_THROWS_WITH(coordinator.calculate_lead_sleeve_requirements(inputs, wire, true),
                      Catch::Matchers::ContainsSubstring("no insulation material with form 'sleeve'"));
    // A lead that needs no sleeve does not touch the stock at all.
    CHECK_FALSE(coordinator.calculate_lead_sleeve_requirements(inputs, wire, false));
    clear_databases();
}

TEST_CASE("Test_Lead_Sleeve_Stock_Is_Never_Wound_As_Tape", "[constructive-model][insulation][lead-sleeve]") {
    clear_databases();
    auto material = find_insulation_material_by_name("PTFE extruded tubing");
    REQUIRE(material.get_form());
    CHECK(material.get_form().value() == Form::SLEEVE);
    // The MAS -> MKF conversion keeps the form (it used to drop it).
    OpenMagnetics::InsulationMaterial converted(static_cast<MAS::InsulationMaterial>(material));
    REQUIRE(converted.get_form());
    CHECK(converted.get_form().value() == Form::SLEEVE);
}

// ---------------------------------------------------------------------------------------------
// The ABT #684 keep-out, conditional on the sleeve
// ---------------------------------------------------------------------------------------------

TEST_CASE("Test_Lead_Sleeve_Margin_Band_Refuses_An_Unsleeved_Lead", "[constructive-model][coil][lead-sleeve]") {
    // Top edge at y = 10 mm, 3 mm margin: the band is [7, 10] mm. A 1 mm run centred at 8 mm lies in it.
    const double edge = 10e-3;
    const double margin = 3e-3;
    CHECK_THROWS_WITH(OpenMagnetics::Coil::check_lead_margin_keep_out(edge, margin, true, 8e-3, 1e-3, false, "test lead"),
                      Catch::Matchers::ContainsSubstring("only a sleeved lead may cross a margin"));
    CHECK_NOTHROW(OpenMagnetics::Coil::check_lead_margin_keep_out(edge, margin, true, 8e-3, 1e-3, true, "test lead"));
    // Flush against the margin's inner face is outside the band.
    CHECK_NOTHROW(OpenMagnetics::Coil::check_lead_margin_keep_out(edge, margin, true, 6.5e-3, 1e-3, false, "test lead"));
    // Bottom edge mirrors it.
    CHECK_THROWS(OpenMagnetics::Coil::check_lead_margin_keep_out(-edge, margin, false, -8e-3, 1e-3, false, "test lead"));
    CHECK_NOTHROW(OpenMagnetics::Coil::check_lead_margin_keep_out(-edge, margin, false, -6.5e-3, 1e-3, false, "test lead"));
}

// ---------------------------------------------------------------------------------------------
// The coil: connections[].sleeve, reservation in the margin band, per-lead report
// ---------------------------------------------------------------------------------------------

TEST_CASE("Test_Lead_Sleeve_Margin_Wound_Flyback_Sleeves_Its_Enamelled_Leads", "[constructive-model][coil][masautocomplete][lead-sleeve]") {
    settings.reset();
    clear_databases();
    auto mas = autocomplete_flyback();
    auto& coil = mas.get_mutable_magnetic().get_mutable_coil();
    REQUIRE(coil.get_turns_description());
    REQUIRE(coil.winding_leads_cross_margin("Primary"));

    for (End end : {End::START, End::FINISH}) {
        auto sleeve = coil.get_recorded_lead_sleeve("Primary", end, 0);
        REQUIRE(sleeve);
        CHECK(sleeve->get_wall_thickness() >= 0.4e-3);
        CHECK_THAT(sleeve->get_wall_thickness(), Catch::Matchers::WithinAbs(0.41e-3, 1e-12));
        REQUIRE(sleeve->get_overlap_into_winding());
    }

    auto result = coil.calculate_insulation_coordination_result();
    REQUIRE(result.get_leads().size() == 2 * coil.get_functional_description().size());
    double margin = coil.get_winding_lead_margin("Primary");
    REQUIRE(margin > 0);
    for (End end : {End::START, End::FINISH}) {
        auto lead = find_lead(result, "Primary", end);
        CHECK(lead.crossesMargin);
        CHECK(lead.sleeved);
        REQUIRE(lead.creepageDistance);
        CHECK(lead.creepageDistance.value() == margin);
        CHECK(lead.requiredCreepageDistance == result.get_creepage_distance());
    }

    // The sleeved edge runs are reserved at the sleeve's diameter and may lie in the margin band.
    auto bobbin = coil.resolve_bobbin();
    auto window = bobbin.get_processed_description().value().get_winding_windows()[0];
    double top = window.get_coordinates().value()[1] + window.get_height().value() / 2;
    double bottom = window.get_coordinates().value()[1] - window.get_height().value() / 2;
    size_t sleevedRuns = 0;
    size_t sleevedRunsInBand = 0;
    for (auto& space : coil.get_connection_reserved_spaces()) {
        if (!space.isTerminal || !space.layer.empty() || !space.sleeveOuterDiameter) {
            continue;
        }
        if (space.kind != ConnectionKind::TERMINAL_ENTRANCE && space.kind != ConnectionKind::TERMINAL_EXIT) {
            continue;
        }
        if (std::abs(space.dimensions[1] - space.sleeveOuterDiameter.value()) > 1e-12) {
            continue;   // a stub or a radial exit, not an edge run
        }
        ++sleevedRuns;
        double faceDepth = std::min(top - (space.coordinates[1] + space.dimensions[1] / 2),
                                    (space.coordinates[1] - space.dimensions[1] / 2) - bottom);
        if (faceDepth < margin - 1e-9) {
            ++sleevedRunsInBand;
        }
    }
    CHECK(sleevedRuns > 0);
    CHECK(sleevedRunsInBand > 0);
    settings.reset();
}

TEST_CASE("Test_Lead_Sleeve_Triple_Insulated_Leads_Are_Not_Sleeved", "[constructive-model][coil][masautocomplete][lead-sleeve]") {
    settings.reset();
    clear_databases();
    auto mas = autocomplete_flyback([](json& masJson) {
        for (auto& winding : masJson["magnetic"]["coil"]["functionalDescription"]) {
            if (winding["isolationSide"] != "primary") {
                winding["wire"] = tripleInsulatedWireName;
            }
        }
        masJson["inputs"]["designRequirements"]["insulation"]["standards"] = json::array({"IEC 62368-1"});
    });
    auto& coil = mas.get_mutable_magnetic().get_mutable_coil();
    REQUIRE(coil.get_turns_description());
    auto result = coil.calculate_insulation_coordination_result();
    for (auto& winding : coil.get_functional_description()) {
        if (winding.get_isolation_side() == IsolationSide::PRIMARY) {
            continue;
        }
        for (End end : {End::START, End::FINISH}) {
            INFO(winding.get_name());
            CHECK_FALSE(coil.get_recorded_lead_sleeve(winding.get_name(), end, 0));
            CHECK_FALSE(find_lead(result, winding.get_name(), end).sleeved);
        }
    }
    for (auto& space : coil.get_connection_reserved_spaces()) {
        if (space.winding != "Primary") {
            CHECK_FALSE(space.sleeveOuterDiameter);
        }
    }
    settings.reset();
}

// ---------------------------------------------------------------------------------------------
// R11
// ---------------------------------------------------------------------------------------------

TEST_CASE("Test_Lead_Sleeve_R11_Passes_When_Sleeves_Follow_The_Requirement", "[adviser][manufacturability][lead-sleeve]") {
    settings.reset();
    clear_databases();
    auto mas = autocomplete_flyback();
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r11_termination_protection(mas.get_mutable_magnetic(), mas.get_mutable_inputs());
    CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
    REQUIRE(finding.get_measured_value());
    CHECK(finding.get_measured_value().value() > 0);
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("Manual labour"));
    settings.reset();
}

TEST_CASE("Test_Lead_Sleeve_R11_Fails_On_A_Missing_Or_Unneeded_Sleeve", "[adviser][manufacturability][lead-sleeve]") {
    settings.reset();
    clear_databases();
    auto mas = autocomplete_flyback();
    auto& coil = mas.get_mutable_magnetic().get_mutable_coil();
    Manufacturability manufacturability;

    SECTION("a required sleeve removed bridges the margin") {
        auto& winding = coil.get_mutable_functional_description()[coil.get_winding_index_by_name("Primary")];
        auto connections = winding.get_connections().value();
        for (auto& connection : connections) {
            if (connection.get_end() && connection.get_end().value() == End::START) {
                connection.set_sleeve(std::nullopt);
            }
        }
        winding.set_connections(connections);
        auto finding = manufacturability.evaluate_r11_termination_protection(mas.get_mutable_magnetic(), mas.get_mutable_inputs());
        CHECK(finding.get_status() == ManufacturabilityStatus::FAIL);
        CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("margin is bridged"));
        REQUIRE(finding.get_scope());
        CHECK_THAT(finding.get_scope().value(), Catch::Matchers::ContainsSubstring("Primary start"));

        auto result = coil.calculate_insulation_coordination_result();
        auto lead = find_lead(result, "Primary", End::START);
        CHECK_FALSE(lead.sleeved);
        REQUIRE(lead.creepageDistance);
        CHECK(lead.creepageDistance.value() == 0);
    }
    SECTION("a sleeve no requirement asks for is the costliest protection spent for nothing") {
        auto insulation = mas.get_inputs().get_design_requirements().get_insulation().value();
        insulation.set_insulation_type(IsolationClass::FUNCTIONAL);
        mas.get_mutable_inputs().get_mutable_design_requirements().set_insulation(insulation);
        auto finding = manufacturability.evaluate_r11_termination_protection(mas.get_mutable_magnetic(), mas.get_mutable_inputs());
        CHECK(finding.get_status() == ManufacturabilityStatus::FAIL);
        CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("costliest protection for nothing"));
    }
    settings.reset();
}
