#include "advisers/Manufacturability.h"
#include "advisers/MagneticFilter.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Wire.h"
#include "support/Settings.h"
#include "TestingUtils.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

using namespace MAS;
using namespace OpenMagnetics;

namespace {

// The wires the fixtures below build on. "Round TCA3 26 AWG" is a triple insulated wire
// (coating insulated, three layers); "Round 0.475 - Grade 1" and "Round 0.28 - Grade 1" are
// plain enamelled magnet wire, the second finer than the 28 AWG automatic-soldering limit.
const std::string tripleInsulatedWireName = "Round TCA3 26 AWG";
const std::string enamelledWireName = "Round 0.475 - Grade 1";
const std::string fineEnamelledWireName = "Round 0.28 - Grade 1";

OpenMagnetics::Inputs make_inputs(std::optional<std::vector<ConnectionType>> terminalTypes = std::nullopt,
                   std::optional<Topology> topology = std::nullopt,
                   bool withInsulationRequirements = false) {
    OpenMagnetics::Inputs inputs;
    DesignRequirements designRequirements;
    DimensionWithTolerance magnetizingInductance;
    magnetizingInductance.set_nominal(100e-6);
    designRequirements.set_magnetizing_inductance(magnetizingInductance);
    designRequirements.set_turns_ratios(std::vector<DimensionWithTolerance>{});
    if (terminalTypes) {
        designRequirements.set_terminal_type(terminalTypes);
    }
    if (topology) {
        designRequirements.set_topology(topology);
    }
    if (withInsulationRequirements) {
        InsulationRequirements insulationRequirements;
        insulationRequirements.set_insulation_type(IsolationClass::REINFORCED);
        designRequirements.set_insulation(insulationRequirements);
    }
    inputs.set_design_requirements(designRequirements);
    inputs.set_operating_points(std::vector<OperatingPoint>{});
    return inputs;
}

/// A wound two-winding magnetic on a named core shape, with a wire per winding.
OpenMagnetics::Magnetic make_wound_magnetic(std::string shapeName,
                             std::vector<int64_t> numberTurns,
                             std::vector<std::string> wireNames,
                             json gapping = json::parse("[]")) {
    std::vector<OpenMagnetics::Wire> wires;
    for (auto& wireName : wireNames) {
        wires.push_back(OpenMagnetics::find_wire_by_name(wireName));
    }
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns,
                                                     std::vector<int64_t>(numberTurns.size(), 1),
                                                     shapeName,
                                                     1,
                                                     WindingOrientation::OVERLAPPING,
                                                     WindingOrientation::OVERLAPPING,
                                                     CoilAlignment::CENTERED,
                                                     CoilAlignment::CENTERED,
                                                     wires);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, "N87");
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    return magnetic;
}

/// A wound single-winding coil on a bobbin of a given window, so the layer count is engineered.
OpenMagnetics::Magnetic make_magnetic_with_layer_count(int64_t numberTurns, double bobbinHeight, double bobbinWidth) {
    auto coil = OpenMagneticsTesting::get_quick_coil({numberTurns},
                                                     {1},
                                                     bobbinHeight,
                                                     bobbinWidth,
                                                     std::vector<double>({0.01, 0, 0}));
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagneticsTesting::get_quick_core("E 25/13/7", json::parse("[]"), 1, "N87"));
    magnetic.set_coil(coil);
    return magnetic;
}

size_t conduction_layers_in_first_section(OpenMagnetics::Magnetic& magnetic) {
    auto& coil = magnetic.get_mutable_coil();
    REQUIRE(coil.get_sections_description());
    size_t count = 0;
    auto sections = coil.get_sections_description().value();
    for (auto& section : sections) {
        if (section.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        for (auto& layer : coil.get_layers_by_section(section.get_name())) {
            if (layer.get_type() == ElectricalType::CONDUCTION) {
                ++count;
            }
        }
        break;
    }
    return count;
}

/// A quick bobbin is processed-only; the DFM fixtures give it the functional description a
/// catalogue bobbin carries, so the material and pinout have somewhere to live.
BobbinFunctionalDescription bobbin_functional_description(const OpenMagnetics::Bobbin& bobbin) {
    if (bobbin.get_functional_description()) {
        return bobbin.get_functional_description().value();
    }
    BobbinFunctionalDescription functionalDescription;
    functionalDescription.set_family(BobbinFamily::E);
    functionalDescription.set_shape("DFM fixture");
    functionalDescription.set_dimensions(std::map<std::string, Dimension>{});
    return functionalDescription;
}

void set_bobbin_material(OpenMagnetics::Magnetic& magnetic, std::optional<std::string> materialName) {
    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    auto functionalDescription = bobbin_functional_description(bobbin);
    if (materialName) {
        functionalDescription.set_material(InsulationMaterialDataOrNameUnion(materialName.value()));
    }
    else {
        functionalDescription.set_material(std::nullopt);
    }
    bobbin.set_functional_description(functionalDescription);
    magnetic.get_mutable_coil().set_bobbin(bobbin);
}

void set_bobbin_pin_pitch(OpenMagnetics::Magnetic& magnetic, double pitch) {
    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    auto functionalDescription = bobbin_functional_description(bobbin);
    Pinout pinout;
    pinout.set_number_pins(4);
    pinout.set_pitch(Pitch(pitch));
    functionalDescription.set_pinout(pinout);
    bobbin.set_functional_description(functionalDescription);
    magnetic.get_mutable_coil().set_bobbin(bobbin);
}

// -----------------------------------------------------------------------------------------
// The report as a whole
// -----------------------------------------------------------------------------------------

TEST_CASE("Test_Manufacturability_Report_Covers_Every_Rule", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    auto inputs = make_inputs();
    Manufacturability manufacturability;
    auto report = manufacturability.calculate_report(magnetic, inputs);

    REQUIRE(report.get_findings().size() == 17);
    for (int ruleNumber = 1; ruleNumber <= 17; ++ruleNumber) {
        auto ruleId = "R" + std::to_string(ruleNumber);
        INFO("rule " << ruleId);
        REQUIRE(report.has_finding(ruleId));
        auto& finding = report.get_finding(ruleId);
        CHECK(!finding.get_title().empty());
        CHECK(!finding.get_source().empty());
        CHECK(!finding.get_message().empty());
        if (finding.get_status() == ManufacturabilityStatus::NOT_EVALUATED) {
            CHECK(finding.get_reason());
        }
    }
}

TEST_CASE("Test_Manufacturability_Rules_Owned_Elsewhere_Are_Honest", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    auto inputs = make_inputs();
    Manufacturability manufacturability;
    auto report = manufacturability.calculate_report(magnetic, inputs);

    for (auto& ruleId : std::vector<std::string>{"R10", "R16"}) {
        INFO("rule " << ruleId);
        auto& finding = report.get_finding(ruleId);
        CHECK(finding.get_status() == ManufacturabilityStatus::NOT_EVALUATED);
        CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("not evaluated (owned by"));
    }
    // ABT #1172: WP3 owns R2, R3, R4 and R14 and evaluates them now. This fixture's quick bobbin
    // has no pins, so the pin rules say exactly that, and it has no shield.
    for (auto& ruleId : std::vector<std::string>{"R2", "R3", "R4"}) {
        INFO("rule " << ruleId);
        auto& finding = report.get_finding(ruleId);
        CHECK(finding.get_status() == ManufacturabilityStatus::NOT_EVALUATED);
        REQUIRE(finding.get_reason());
        CHECK_THAT(finding.get_reason().value(), Catch::Matchers::ContainsSubstring("pin"));
    }
    CHECK(report.get_finding("R14").get_status() == ManufacturabilityStatus::NOT_APPLICABLE);
}

TEST_CASE("Test_Manufacturability_Rule_Numbers_Come_From_The_Data_File", "[adviser][manufacturability]") {
    Manufacturability manufacturability;
    auto& rules = manufacturability.get_rules();
    // Every number the rules compare against is data, with the Wuerth source that states it.
    CHECK(rules.at("R8").at("outerDiameterGrowthTripleInsulated").get<double>() == 0.65);
    CHECK(rules.at("R9").at("typicalLayersPerSideMinimum").get<int>() == 10);
    CHECK(rules.at("R9").at("typicalLayersPerSideMaximum").get<int>() == 20);
    CHECK(rules.at("R13").at("finestAutomaticallySolderableAwg").get<int>() == 28);
    CHECK(rules.at("R17").at("partialDischargeTestPeakVoltage").get<double>() == 750);
    for (auto& ruleId : std::vector<std::string>{"R1", "R5", "R6", "R7", "R8", "R9", "R11", "R12", "R13", "R15", "R17"}) {
        INFO("rule " << ruleId);
        CHECK(!rules.at(ruleId).at("source").get<std::string>().empty());
    }
}

// -----------------------------------------------------------------------------------------
// R1 - even layer count
// -----------------------------------------------------------------------------------------

TEST_CASE("Test_Manufacturability_R1_Odd_Layer_Count_Warns", "[adviser][manufacturability]") {
    // 27 turns of 0.475 mm wire in a 5 mm high window land on three layers.
    auto magnetic = make_magnetic_with_layer_count(27, 0.005, 0.01);
    REQUIRE(conduction_layers_in_first_section(magnetic) == 3);

    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r1_layer_parity(magnetic);
    CHECK(finding.get_status() == ManufacturabilityStatus::WARNING);
    REQUIRE(finding.get_measured_value());
    CHECK(finding.get_measured_value().value() == 3);
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("drag-back"));
}

TEST_CASE("Test_Manufacturability_R1_Even_Layer_Count_Passes", "[adviser][manufacturability]") {
    auto magnetic = make_magnetic_with_layer_count(18, 0.005, 0.01);
    REQUIRE(conduction_layers_in_first_section(magnetic) == 2);

    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r1_layer_parity(magnetic);
    CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
}

TEST_CASE("Test_Manufacturability_R1_Single_Layer_Passes", "[adviser][manufacturability]") {
    auto magnetic = make_magnetic_with_layer_count(5, 0.005, 0.01);
    REQUIRE(conduction_layers_in_first_section(magnetic) == 1);

    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r1_layer_parity(magnetic);
    CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
}

TEST_CASE("Test_Manufacturability_R1_Unwound_Coil_Is_Not_Evaluated", "[adviser][manufacturability]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("E 25/13/7", OpenMagneticsTesting::get_residual_gap(), {20, 20});
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r1_layer_parity(magnetic);
    CHECK(finding.get_status() == ManufacturabilityStatus::NOT_EVALUATED);
    REQUIRE(finding.get_reason());
    CHECK_THAT(finding.get_reason().value(), Catch::Matchers::ContainsSubstring("has not been wound"));
}

TEST_CASE("Test_Manufacturability_R1_Layer_Parity_Filter", "[adviser][manufacturability]") {
    // The adviser-side half of R1: the filter scores the penalty, it never rejects.
    auto magneticOdd = make_magnetic_with_layer_count(27, 0.005, 0.01);
    auto sectionsOdd = magneticOdd.get_mutable_coil().get_sections_description().value();
    auto windingOdd = magneticOdd.get_coil().get_functional_description()[0];
    auto filter = MagneticFilterLayerParity();
    auto [validOdd, scoringOdd] = filter.evaluate_magnetic(windingOdd, sectionsOdd[0]);
    CHECK(validOdd);
    CHECK(scoringOdd == 1.0);

    auto magneticEven = make_magnetic_with_layer_count(18, 0.005, 0.01);
    auto sectionsEven = magneticEven.get_mutable_coil().get_sections_description().value();
    auto windingEven = magneticEven.get_coil().get_functional_description()[0];
    auto [validEven, scoringEven] = filter.evaluate_magnetic(windingEven, sectionsEven[0]);
    CHECK(validEven);
    CHECK(scoringEven == 0.0);
}

TEST_CASE("Test_Manufacturability_R1_Adviser_Penalty_Is_Opt_In", "[adviser][manufacturability]") {
    // The ranking of every existing adviser is unchanged until the caller opts in.
    CHECK(Settings::GetInstance().get_wire_adviser_penalize_odd_layer_count() == false);
}

// -----------------------------------------------------------------------------------------
// R5 - heavy wire on a surface-mount part
// -----------------------------------------------------------------------------------------

TEST_CASE("Test_Manufacturability_R5_Heavy_Wire_On_Smt_Warns", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    set_bobbin_pin_pitch(magnetic, 0.0005);  // 0.5 mm pitch against a ~0.5 mm wire
    auto inputs = make_inputs(std::vector<ConnectionType>{ConnectionType::SMT});

    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r5_smt_heavy_wire(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::WARNING);
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("Split the winding into parallels"));
}

TEST_CASE("Test_Manufacturability_R5_Fine_Wire_On_Smt_Passes", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    set_bobbin_pin_pitch(magnetic, 0.005);  // 5 mm pitch, ten times the wire
    auto inputs = make_inputs(std::vector<ConnectionType>{ConnectionType::SMT});

    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r5_smt_heavy_wire(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
}

TEST_CASE("Test_Manufacturability_R5_Through_Hole_Is_Not_Applicable", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    auto inputs = make_inputs(std::vector<ConnectionType>{ConnectionType::THT});
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r5_smt_heavy_wire(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::NOT_APPLICABLE);
}

TEST_CASE("Test_Manufacturability_R5_Missing_Terminal_Type_Is_Not_Evaluated", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    auto inputs = make_inputs();
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r5_smt_heavy_wire(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::NOT_EVALUATED);
    CHECK(finding.get_reason());
}

// -----------------------------------------------------------------------------------------
// R6 - thermoset bobbin
// -----------------------------------------------------------------------------------------

TEST_CASE("Test_Manufacturability_R6_Material_Classification", "[adviser][manufacturability]") {
    Manufacturability manufacturability;
    for (auto& name : std::vector<std::string>{"PA66", "PBT", "PET", "PPS", "LCP", "Nylon 6/6"}) {
        INFO(name);
        CHECK(manufacturability.classify_bobbin_material(name) == Manufacturability::BobbinMaterialClass::THERMOPLASTIC);
    }
    for (auto& name : std::vector<std::string>{"Phenolic", "DAP", "Diallyl phthalate", "Bakelite"}) {
        INFO(name);
        CHECK(manufacturability.classify_bobbin_material(name) == Manufacturability::BobbinMaterialClass::THERMOSET);
    }
    // A trade name nobody has classified is UNKNOWN, never one of the two by default.
    CHECK(manufacturability.classify_bobbin_material("SKYT.5220FR") == Manufacturability::BobbinMaterialClass::UNKNOWN);
}

TEST_CASE("Test_Manufacturability_R6_Thermoplastic_On_Smt_Fails", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    set_bobbin_material(magnetic, "PA66");
    auto inputs = make_inputs(std::vector<ConnectionType>{ConnectionType::SMT});

    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r6_thermoset_bobbin(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::FAIL);
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("thermoset bobbin"));
}

TEST_CASE("Test_Manufacturability_R6_Thermoset_On_Smt_Passes", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    set_bobbin_material(magnetic, "Phenolic");
    auto inputs = make_inputs(std::vector<ConnectionType>{ConnectionType::SMT});

    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r6_thermoset_bobbin(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
}

TEST_CASE("Test_Manufacturability_R6_Missing_Material_Is_Unknown_Not_A_Default", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    set_bobbin_material(magnetic, std::nullopt);
    auto inputs = make_inputs(std::vector<ConnectionType>{ConnectionType::SMT});

    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r6_thermoset_bobbin(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::NOT_EVALUATED);
    REQUIRE(finding.get_reason());
    CHECK_THAT(finding.get_reason().value(), Catch::Matchers::ContainsSubstring("unknown"));
}

TEST_CASE("Test_Manufacturability_R6_Unclassified_Material_Is_Unknown", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    set_bobbin_material(magnetic, "SKYT.5220FR");
    auto inputs = make_inputs(std::vector<ConnectionType>{ConnectionType::SMT});

    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r6_thermoset_bobbin(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::NOT_EVALUATED);
}

// -----------------------------------------------------------------------------------------
// R7 - size heuristic
// -----------------------------------------------------------------------------------------

TEST_CASE("Test_Manufacturability_R7_Insulated_Wire_On_A_Big_Core_Warns", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 55/28/21", {20, 20}, {enamelledWireName, tripleInsulatedWireName});
    auto inputs = make_inputs(std::nullopt, std::nullopt, true);

    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r7_size_heuristic(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::WARNING);
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("margin tape"));
}

TEST_CASE("Test_Manufacturability_R7_Insulated_Wire_On_A_Small_Core_Passes", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 16/8/5", {10, 10}, {enamelledWireName, tripleInsulatedWireName});
    auto inputs = make_inputs(std::nullopt, std::nullopt, true);

    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r7_size_heuristic(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
}

TEST_CASE("Test_Manufacturability_R7_Without_Insulation_Requirements_Is_Not_Applicable", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 55/28/21", {20, 20}, {enamelledWireName, tripleInsulatedWireName});
    auto inputs = make_inputs();
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r7_size_heuristic(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::NOT_APPLICABLE);
}

// -----------------------------------------------------------------------------------------
// R8 - insulated wire growth, placement and contact
// -----------------------------------------------------------------------------------------

TEST_CASE("Test_Manufacturability_R8_Wire_Classification", "[adviser][manufacturability]") {
    auto tripleInsulated = OpenMagnetics::find_wire_by_name(tripleInsulatedWireName);
    auto enamelled = OpenMagnetics::find_wire_by_name(enamelledWireName);
    CHECK(Manufacturability::is_insulated_wire(tripleInsulated));
    CHECK(Manufacturability::is_triple_insulated_wire(tripleInsulated));
    CHECK(!Manufacturability::is_enamelled_magnet_wire(tripleInsulated));
    CHECK(!Manufacturability::is_insulated_wire(enamelled));
    CHECK(Manufacturability::is_enamelled_magnet_wire(enamelled));
}

TEST_CASE("Test_Manufacturability_R8_Insulated_Next_To_Enamelled_Without_Tape_Fails", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {10, 40}, {enamelledWireName, tripleInsulatedWireName});
    // The winder puts an insulation section between windings; strip it so the triple
    // insulated winding sits directly on the enamelled one, which is the case the rule forbids.
    auto& coil = magnetic.get_mutable_coil();
    std::vector<Section> conductionOnly;
    auto woundSections = coil.get_sections_description().value();
    for (auto& section : woundSections) {
        if (section.get_type() == ElectricalType::CONDUCTION) {
            conductionOnly.push_back(section);
        }
    }
    REQUIRE(conductionOnly.size() == 2);
    coil.set_sections_description(conductionOnly);
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r8_insulated_wire(magnetic);
    CHECK(finding.get_status() == ManufacturabilityStatus::FAIL);
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("touches enamelled magnet wire"));
}

TEST_CASE("Test_Manufacturability_R8_Tape_Between_Windings_Only_Warns_About_Placement", "[adviser][manufacturability]") {
    // The same windings as wound, with the winder's insulation section between them: no
    // contact, but the insulated wire is on the 40-turn winding instead of the 10-turn one.
    auto magnetic = make_wound_magnetic("E 25/13/7", {10, 40}, {enamelledWireName, tripleInsulatedWireName});
    REQUIRE(!magnetic.get_mutable_coil().get_sections_description_insulation().empty());
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r8_insulated_wire(magnetic);
    CHECK(finding.get_status() == ManufacturabilityStatus::WARNING);
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("lowest-turn winding"));
}

TEST_CASE("Test_Manufacturability_R8_Insulated_Wire_On_The_Lowest_Turn_Winding_Passes", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {40, 10}, {enamelledWireName, tripleInsulatedWireName});
    REQUIRE(!magnetic.get_mutable_coil().get_sections_description_insulation().empty());
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r8_insulated_wire(magnetic);
    CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
}

TEST_CASE("Test_Manufacturability_R8_No_Insulated_Wire_Is_Not_Applicable", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {10, 40}, {enamelledWireName, enamelledWireName});
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r8_insulated_wire(magnetic);
    CHECK(finding.get_status() == ManufacturabilityStatus::NOT_APPLICABLE);
}

TEST_CASE("Test_Manufacturability_R8_Growth_Over_Enamelled_Is_Reported", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {10, 40}, {enamelledWireName, tripleInsulatedWireName});
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r8_insulated_wire(magnetic);
    REQUIRE(finding.get_measured_value());
    // A triple insulated wire is materially fatter than its enamelled equivalent.
    CHECK(finding.get_measured_value().value() > 0.1);
}

// -----------------------------------------------------------------------------------------
// R9 - margin tape
// -----------------------------------------------------------------------------------------

TEST_CASE("Test_Manufacturability_R9_No_Margins_Is_Not_Applicable", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r9_margin_tape(magnetic);
    CHECK(finding.get_status() == ManufacturabilityStatus::NOT_APPLICABLE);
}

TEST_CASE("Test_Manufacturability_R9_Margins_Report_The_Window_Loss", "[adviser][manufacturability]") {
    auto coil = OpenMagneticsTesting::get_quick_coil({20, 20}, {1, 1}, "E 25/13/7", 1);
    coil.add_margin_to_section_by_index(0, std::vector<double>{0.001, 0.001});
    coil.add_margin_to_section_by_index(1, std::vector<double>{0.001, 0.001});
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(OpenMagneticsTesting::get_quick_core("E 25/13/7", json::parse("[]"), 1, "N87"));
    magnetic.set_coil(coil);

    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r9_margin_tape(magnetic);
    REQUIRE(finding.get_measured_value());
    CHECK(finding.get_measured_value().value() > 0);
    CHECK(finding.get_status() != ManufacturabilityStatus::NOT_APPLICABLE);
    CHECK(finding.get_status() != ManufacturabilityStatus::NOT_EVALUATED);
    // A bare [top, bottom] margin knows no tape-layer count, and the rule does not invent one.
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("not recorded"));
    // The leakage delta is explicitly declared out of scope rather than silently omitted.
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("leakage"));
}

// -----------------------------------------------------------------------------------------
// R12 - flying leads
// -----------------------------------------------------------------------------------------

TEST_CASE("Test_Manufacturability_R12_Flying_Lead_Requirement_Warns", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    auto inputs = make_inputs(std::vector<ConnectionType>{ConnectionType::FLYING_LEAD});
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r12_flying_leads(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::WARNING);
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("least preferred"));
}

TEST_CASE("Test_Manufacturability_R12_Pin_Termination_Passes", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    auto inputs = make_inputs(std::vector<ConnectionType>{ConnectionType::PIN});
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r12_flying_leads(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
}

// -----------------------------------------------------------------------------------------
// R13 - manual termination
// -----------------------------------------------------------------------------------------

TEST_CASE("Test_Manufacturability_R13_Coarse_Enamel_Passes", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r13_manual_termination(magnetic);
    CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
}

TEST_CASE("Test_Manufacturability_R13_Fine_Wire_Needs_Manual_Termination", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, fineEnamelledWireName});
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r13_manual_termination(magnetic);
    CHECK(finding.get_status() == ManufacturabilityStatus::WARNING);
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("28 AWG"));
}

TEST_CASE("Test_Manufacturability_R13_Insulated_Wire_Needs_Stripping", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, tripleInsulatedWireName});
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r13_manual_termination(magnetic);
    CHECK(finding.get_status() == ManufacturabilityStatus::WARNING);
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("pre-stripped"));
}

// -----------------------------------------------------------------------------------------
// R15 - gap only in the centre leg
// -----------------------------------------------------------------------------------------

TEST_CASE("Test_Manufacturability_R15_Lateral_Additive_Gaps_On_A_Flyback_Warn", "[adviser][manufacturability]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("E 25/13/7", OpenMagneticsTesting::get_spacer_gap(0.001), {20, 20});
    auto inputs = make_inputs(std::nullopt, Topology::FLYBACK_CONVERTER);
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r15_lateral_gaps(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::WARNING);
    REQUIRE(finding.get_measured_value());
    CHECK(finding.get_measured_value().value() == 2);
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("centre leg"));
}

TEST_CASE("Test_Manufacturability_R15_Central_Gap_Only_Passes", "[adviser][manufacturability]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("E 25/13/7", OpenMagneticsTesting::get_ground_gap(0.001), {20, 20});
    auto inputs = make_inputs(std::nullopt, Topology::FLYBACK_CONVERTER);
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r15_lateral_gaps(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::PASS);
}

TEST_CASE("Test_Manufacturability_R15_Non_Isolated_Topology_Is_Not_Applicable", "[adviser][manufacturability]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("E 25/13/7", OpenMagneticsTesting::get_spacer_gap(0.001), {20, 20});
    auto inputs = make_inputs(std::nullopt, Topology::BUCK_CONVERTER);
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r15_lateral_gaps(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::NOT_APPLICABLE);
}

TEST_CASE("Test_Manufacturability_R15_Missing_Topology_Is_Not_Evaluated", "[adviser][manufacturability]") {
    auto magnetic = OpenMagneticsTesting::get_quick_magnetic("E 25/13/7", OpenMagneticsTesting::get_spacer_gap(0.001), {20, 20});
    auto inputs = make_inputs();
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r15_lateral_gaps(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::NOT_EVALUATED);
    CHECK(finding.get_reason());
}

// -----------------------------------------------------------------------------------------
// R17 - production test and tolerance fields
// -----------------------------------------------------------------------------------------

TEST_CASE("Test_Manufacturability_R17_Reports_The_Production_Test_Fields", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, enamelledWireName});
    auto inputs = make_inputs();
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r17_production_tests(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::INFORMATIONAL);
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("hipot"));
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("Partial-discharge test: not required"));
}

TEST_CASE("Test_Manufacturability_R17_Partial_Discharge_Depends_On_Insulated_Wire", "[adviser][manufacturability]") {
    auto magnetic = make_wound_magnetic("E 25/13/7", {20, 20}, {enamelledWireName, tripleInsulatedWireName});
    auto inputs = make_inputs();
    Manufacturability manufacturability;
    auto finding = manufacturability.evaluate_r17_production_tests(magnetic, inputs);
    CHECK(finding.get_status() == ManufacturabilityStatus::INFORMATIONAL);
    // No operating point, so the peak working voltage is not known and the rule says so
    // instead of assuming the threshold is or is not crossed.
    CHECK_THAT(finding.get_message(), Catch::Matchers::ContainsSubstring("not evaluated"));
}

} // namespace
