#include "advisers/Manufacturability.h"

#include "constructive_models/Bobbin.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Wire.h"
#include "support/Utils.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cmrc/cmrc.hpp>
#include <magic_enum.hpp>
#include <sstream>
#include <stdexcept>

CMRC_DECLARE(dfmData);

namespace OpenMagnetics {

namespace {

std::string to_upper(const std::string& text) {
    std::string result = text;
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c) { return std::toupper(c); });
    return result;
}

std::string join(const std::vector<std::string>& items, const std::string& separator) {
    std::string result;
    for (size_t index = 0; index < items.size(); ++index) {
        if (index > 0) {
            result += separator;
        }
        result += items[index];
    }
    return result;
}

/// The material name behind a MAS insulation-material union, or nullopt when the union
/// carries an object that does not name itself.
std::optional<std::string> insulation_material_name(const InsulationMaterialDataOrNameUnion& material) {
    if (std::holds_alternative<std::string>(material)) {
        return std::get<std::string>(material);
    }
    auto data = std::get<MAS::InsulationMaterial>(material);
    if (data.get_name().empty()) {
        return std::nullopt;
    }
    return data.get_name();
}

/// The name of the winding a section belongs to, or nullopt when the section names none.
std::optional<std::string> section_winding_name(const Section& section) {
    if (section.get_partial_windings().empty()) {
        return std::nullopt;
    }
    return section.get_partial_windings()[0].get_winding();
}

bool terminal_types_contain(const Inputs& inputs, ConnectionType wanted) {
    auto terminalTypes = inputs.get_design_requirements().get_terminal_type();
    if (!terminalTypes) {
        return false;
    }
    return std::find(terminalTypes->begin(), terminalTypes->end(), wanted) != terminalTypes->end();
}

} // namespace

void to_json(json& j, const ManufacturabilityFinding& x) {
    j = json::object();
    j["ruleId"] = x.get_rule_id();
    j["title"] = x.get_title();
    j["status"] = std::string(magic_enum::enum_name(x.get_status()));
    j["message"] = x.get_message();
    j["source"] = x.get_source();
    if (x.get_measured_value()) {
        j["measuredValue"] = x.get_measured_value().value();
    }
    if (x.get_threshold_value()) {
        j["thresholdValue"] = x.get_threshold_value().value();
    }
    if (x.get_unit()) {
        j["unit"] = x.get_unit().value();
    }
    if (x.get_scope()) {
        j["scope"] = x.get_scope().value();
    }
    if (x.get_reason()) {
        j["reason"] = x.get_reason().value();
    }
}

void to_json(json& j, const ManufacturabilityReport& x) {
    j = json::object();
    j["findings"] = json::array();
    for (auto& finding : x.get_findings()) {
        json findingJson;
        to_json(findingJson, finding);
        j["findings"].push_back(findingJson);
    }
    j["numberPass"] = x.count_by_status(ManufacturabilityStatus::PASS);
    j["numberWarnings"] = x.count_by_status(ManufacturabilityStatus::WARNING);
    j["numberFails"] = x.count_by_status(ManufacturabilityStatus::FAIL);
    j["numberNotEvaluated"] = x.count_by_status(ManufacturabilityStatus::NOT_EVALUATED);
    j["numberNotApplicable"] = x.count_by_status(ManufacturabilityStatus::NOT_APPLICABLE);
    j["numberInformational"] = x.count_by_status(ManufacturabilityStatus::INFORMATIONAL);
}

const ManufacturabilityFinding& ManufacturabilityReport::get_finding(const std::string& ruleId) const {
    for (auto& finding : findings) {
        if (finding.get_rule_id() == ruleId) {
            return finding;
        }
    }
    throw std::runtime_error("Manufacturability report carries no finding for rule " + ruleId);
}

bool ManufacturabilityReport::has_finding(const std::string& ruleId) const {
    for (auto& finding : findings) {
        if (finding.get_rule_id() == ruleId) {
            return true;
        }
    }
    return false;
}

size_t ManufacturabilityReport::count_by_status(ManufacturabilityStatus status) const {
    size_t count = 0;
    for (auto& finding : findings) {
        if (finding.get_status() == status) {
            ++count;
        }
    }
    return count;
}

Manufacturability::Manufacturability() {
    auto fs = cmrc::dfmData::get_filesystem();
    auto data = fs.open("src/data/dfm_rules.json");
    _rules = json::parse(std::string(data.begin(), data.end()));
}

const json& Manufacturability::rule(const std::string& ruleId) const {
    if (!_rules.contains(ruleId)) {
        throw std::runtime_error("src/data/dfm_rules.json has no rule block for " + ruleId);
    }
    return _rules.at(ruleId);
}

double Manufacturability::rule_number(const std::string& ruleId, const std::string& key) const {
    auto& block = rule(ruleId);
    if (!block.contains(key)) {
        throw std::runtime_error("src/data/dfm_rules.json rule " + ruleId + " has no number '" + key + "'");
    }
    return block.at(key).get<double>();
}

std::string Manufacturability::rule_string(const std::string& ruleId, const std::string& key) const {
    auto& block = rule(ruleId);
    if (!block.contains(key)) {
        throw std::runtime_error("src/data/dfm_rules.json rule " + ruleId + " has no string '" + key + "'");
    }
    return block.at(key).get<std::string>();
}

ManufacturabilityFinding Manufacturability::make_finding(const std::string& ruleId) const {
    return ManufacturabilityFinding(ruleId, rule_string(ruleId, "title"), rule_string(ruleId, "source"));
}

ManufacturabilityFinding Manufacturability::not_evaluated_rule(const std::string& ruleId) {
    if (!_rules.contains("notEvaluated") || !_rules.at("notEvaluated").contains(ruleId)) {
        throw std::runtime_error("src/data/dfm_rules.json has no notEvaluated entry for " + ruleId);
    }
    ManufacturabilityFinding finding(ruleId,
                                     _rules.at("notEvaluatedTitles").at(ruleId).get<std::string>(),
                                     _rules.at("notEvaluatedSources").at(ruleId).get<std::string>());
    finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
    finding.set_message(_rules.at("notEvaluated").at(ruleId).get<std::string>());
    finding.set_reason(_rules.at("notEvaluated").at(ruleId).get<std::string>());
    return finding;
}

Manufacturability::BobbinMaterialClass Manufacturability::classify_bobbin_material(const std::string& materialName) const {
    auto upperName = to_upper(materialName);
    // Thermoset first: an "epoxy-coated phenolic" must not be caught by a thermoplastic token.
    for (auto& token : rule("R6").at("thermosetMaterials")) {
        if (upperName.find(to_upper(token.get<std::string>())) != std::string::npos) {
            return BobbinMaterialClass::THERMOSET;
        }
    }
    for (auto& token : rule("R6").at("thermoplasticMaterials")) {
        if (upperName.find(to_upper(token.get<std::string>())) != std::string::npos) {
            return BobbinMaterialClass::THERMOPLASTIC;
        }
    }
    return BobbinMaterialClass::UNKNOWN;
}

bool Manufacturability::is_insulated_wire(Wire& wire) {
    auto coating = wire.resolve_coating();
    if (!coating || !coating->get_type()) {
        return false;
    }
    return coating->get_type().value() == InsulationWireCoatingType::INSULATED ||
           coating->get_type().value() == InsulationWireCoatingType::EXTRUDED;
}

bool Manufacturability::is_triple_insulated_wire(Wire& wire) {
    auto coating = wire.resolve_coating();
    if (!coating || !coating->get_type()) {
        return false;
    }
    if (coating->get_type().value() != InsulationWireCoatingType::INSULATED) {
        return false;
    }
    if (!coating->get_number_layers()) {
        return false;
    }
    return coating->get_number_layers().value() >= 3;
}

bool Manufacturability::is_fully_insulated_wire(Wire& wire) {
    // MKF's coating-label convention (Wire::encode_coating_label): an ENAMELLED wire whose
    // build grade is outside the standard 1..3 range is a fully insulated wire (FIW).
    auto coating = wire.resolve_coating();
    if (!coating || !coating->get_type()) {
        return false;
    }
    if (coating->get_type().value() != InsulationWireCoatingType::ENAMELLED) {
        return false;
    }
    if (!coating->get_grade()) {
        return false;
    }
    auto grade = coating->get_grade().value();
    return grade < 1 || grade > 3;
}

bool Manufacturability::is_enamelled_magnet_wire(Wire& wire) {
    auto coating = wire.resolve_coating();
    if (!coating || !coating->get_type()) {
        return false;
    }
    if (coating->get_type().value() != InsulationWireCoatingType::ENAMELLED) {
        return false;
    }
    return !is_fully_insulated_wire(wire);
}

// ---------------------------------------------------------------------------------------
// R1 - even layer count per winding
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r1_layer_parity(Magnetic& magnetic) {
    auto finding = make_finding("R1");
    finding.set_unit("layers");

    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil, so there are no layers to count");
        finding.set_message("Layer parity not evaluated: " + finding.get_reason().value());
        return finding;
    }

    auto& coil = magnetic.get_mutable_coil();
    if (!coil.get_sections_description() || !coil.get_layers_description()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the coil has not been wound (no sections or layers description), so no layer count exists");
        finding.set_message("Layer parity not evaluated: " + finding.get_reason().value());
        return finding;
    }

    auto singleLayerLimit = static_cast<size_t>(rule_number("R1", "maximumLayersConsideredSingleLayer"));

    std::vector<std::string> oddSections;
    size_t worstLayerCount = 0;
    size_t maximumLayerCount = 0;
    auto sectionsDescription = coil.get_sections_description().value();
    for (auto& section : sectionsDescription) {
        if (section.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        size_t layerCount = 0;
        for (auto& layer : coil.get_layers_by_section(section.get_name())) {
            if (layer.get_type() == ElectricalType::CONDUCTION) {
                ++layerCount;
            }
        }
        maximumLayerCount = std::max(maximumLayerCount, layerCount);
        if (layerCount > singleLayerLimit && (layerCount % 2) != 0) {
            oddSections.push_back(section.get_name() + " (" + std::to_string(layerCount) + " layers)");
            if (layerCount > worstLayerCount) {
                worstLayerCount = layerCount;
            }
        }
    }

    if (oddSections.empty()) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_measured_value(static_cast<double>(maximumLayerCount));
        finding.set_message("Every wound section ends on an even layer count (or a single layer); no drag-back is forced.");
        return finding;
    }

    finding.set_status(ManufacturabilityStatus::WARNING);
    finding.set_measured_value(static_cast<double>(worstLayerCount));
    finding.set_scope(join(oddSections, ", "));
    finding.set_message("An odd layer count forces a drag-back on: " + join(oddSections, ", ") +
                        ". Adjust the wire diameter or the strand count so the layer fills exactly and the section winds in an even number of layers.");
    return finding;
}

// ---------------------------------------------------------------------------------------
// R5 - heavy wire on a surface-mount part
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r5_smt_heavy_wire(Magnetic& magnetic, Inputs& inputs) {
    auto finding = make_finding("R5");
    finding.set_unit("outer diameter / pin pitch");

    if (!inputs.get_design_requirements().get_terminal_type()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("designRequirements has no terminalType, so it is unknown whether the part is surface mounted");
        finding.set_message("SMT heavy-wire check not evaluated: " + finding.get_reason().value());
        return finding;
    }
    if (!terminal_types_contain(inputs, ConnectionType::SMT)) {
        finding.set_status(ManufacturabilityStatus::NOT_APPLICABLE);
        finding.set_message("The part is not surface mounted, so terminal coplanarity is not at stake.");
        return finding;
    }
    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil, so there is no wire to measure");
        finding.set_message("SMT heavy-wire check not evaluated: " + finding.get_reason().value());
        return finding;
    }

    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_functional_description() || !bobbin.get_functional_description()->get_pinout() ||
        !bobbin.get_functional_description()->get_pinout()->get_pitch()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the bobbin record carries no pinout pitch, so the pin pitch the wire has to clear is unknown");
        finding.set_message("SMT heavy-wire check not evaluated: " + finding.get_reason().value());
        return finding;
    }

    auto pitchVariant = bobbin.get_functional_description()->get_pinout()->get_pitch().value();
    double pinPitch;
    if (std::holds_alternative<double>(pitchVariant)) {
        pinPitch = std::get<double>(pitchVariant);
    }
    else {
        auto pitches = std::get<std::vector<double>>(pitchVariant);
        if (pitches.empty()) {
            finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
            finding.set_reason("the bobbin pinout pitch list is empty");
            finding.set_message("SMT heavy-wire check not evaluated: " + finding.get_reason().value());
            return finding;
        }
        pinPitch = *std::min_element(pitches.begin(), pitches.end());
    }
    if (!(pinPitch > 0)) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the bobbin pinout pitch is not a positive length");
        finding.set_message("SMT heavy-wire check not evaluated: " + finding.get_reason().value());
        return finding;
    }

    double ratioLimit = rule_number("R5", "maximumWireOuterDiameterOverPinPitch");
    finding.set_threshold_value(ratioLimit);

    double worstRatio = 0;
    std::vector<std::string> heavyWindings;
    for (auto& winding : magnetic.get_coil().get_functional_description()) {
        auto wire = Coil::resolve_wire(winding);
        double outerDimension = wire.get_maximum_outer_dimension();
        double ratio = outerDimension / pinPitch;
        worstRatio = std::max(worstRatio, ratio);
        if (ratio > ratioLimit) {
            heavyWindings.push_back(winding.get_name());
        }
    }
    finding.set_measured_value(worstRatio);

    if (heavyWindings.empty()) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_message("No winding wire exceeds half the pin pitch on this surface-mount part.");
        return finding;
    }

    finding.set_status(ManufacturabilityStatus::WARNING);
    finding.set_scope(join(heavyWindings, ", "));
    finding.set_message("Heavy wire on a surface-mount part risks coplanarity and height failures at the terminals (" +
                        join(heavyWindings, ", ") + "). Split the winding into parallels and add terminals.");
    return finding;
}

// ---------------------------------------------------------------------------------------
// R6 - thermoset bobbin when the part sees reflow or a solder bath
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r6_thermoset_bobbin(Magnetic& magnetic, Inputs& inputs) {
    auto finding = make_finding("R6");

    if (!inputs.get_design_requirements().get_terminal_type()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("designRequirements has no terminalType, so it is unknown whether the part goes through a reflow oven");
        finding.set_message("Bobbin material class not evaluated: " + finding.get_reason().value());
        return finding;
    }
    if (!terminal_types_contain(inputs, ConnectionType::SMT)) {
        finding.set_status(ManufacturabilityStatus::NOT_APPLICABLE);
        finding.set_message("The part is not surface mounted, so it does not pass through a reflow oven and a thermoplastic bobbin is allowed.");
        return finding;
    }
    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil, so it names no bobbin");
        finding.set_message("Bobbin material class not evaluated: " + finding.get_reason().value());
        return finding;
    }

    auto bobbin = magnetic.get_mutable_coil().resolve_bobbin();
    if (!bobbin.get_functional_description() || !bobbin.get_functional_description()->get_material()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("unknown: the bobbin record carries no material, and this rule never assumes one");
        finding.set_message("Bobbin material class unknown: the bobbin record carries no material. A reflowed part must use a thermoset bobbin; record the bobbin material before this can be checked.");
        return finding;
    }

    auto materialName = insulation_material_name(bobbin.get_functional_description()->get_material().value());
    if (!materialName) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("unknown: the bobbin material object carries no name, and this rule never assumes one");
        finding.set_message("Bobbin material class unknown: the bobbin material object carries no name.");
        return finding;
    }
    finding.set_scope(materialName.value());

    auto materialClass = classify_bobbin_material(materialName.value());
    switch (materialClass) {
        case BobbinMaterialClass::THERMOSET:
            finding.set_status(ManufacturabilityStatus::PASS);
            finding.set_message("The bobbin material '" + materialName.value() + "' is a thermoset, which is what a reflowed part requires.");
            return finding;
        case BobbinMaterialClass::THERMOPLASTIC:
            finding.set_status(ManufacturabilityStatus::FAIL);
            finding.set_message("The bobbin material '" + materialName.value() +
                                "' is a thermoplastic. Its lower melting point lets the pins loosen and the bobbin deform in a reflow oven; a surface-mount part must use a thermoset bobbin.");
            return finding;
        case BobbinMaterialClass::UNKNOWN:
        default:
            finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
            finding.set_reason("unknown: the bobbin material '" + materialName.value() +
                               "' is in neither the thermoplastic nor the thermoset table of src/data/dfm_rules.json");
            finding.set_message("Bobbin material class unknown for '" + materialName.value() +
                                "'. Add it to the classification tables in src/data/dfm_rules.json rather than assuming a class.");
            return finding;
    }
}

// ---------------------------------------------------------------------------------------
// R7 - insulation-strategy size heuristic
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r7_size_heuristic(Magnetic& magnetic, Inputs& inputs) {
    auto finding = make_finding("R7");
    finding.set_unit("m (core effective length)");

    if (!inputs.get_design_requirements().get_insulation()) {
        finding.set_status(ManufacturabilityStatus::NOT_APPLICABLE);
        finding.set_message("The design requirements ask for no insulation coordination, so neither margin tape nor insulated wire is being chosen between.");
        return finding;
    }
    if (!magnetic.has_core() || !magnetic.get_core().get_processed_description()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the core has no processed description, so its effective length is unknown");
        finding.set_message("Size heuristic not evaluated: " + finding.get_reason().value());
        return finding;
    }
    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil, so the insulation strategy in use is unknown");
        finding.set_message("Size heuristic not evaluated: " + finding.get_reason().value());
        return finding;
    }

    auto& coil = magnetic.get_mutable_coil();
    if (!coil.get_sections_description()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the coil has not been wound, so it is unknown whether margin tape is in use");
        finding.set_message("Size heuristic not evaluated: " + finding.get_reason().value());
        return finding;
    }

    double effectiveLength = magnetic.get_mutable_core().get_effective_length();
    finding.set_measured_value(effectiveLength);

    bool usesInsulatedWire = false;
    for (auto& winding : coil.get_functional_description()) {
        auto wire = Coil::resolve_wire(winding);
        if (is_insulated_wire(wire) || is_fully_insulated_wire(wire)) {
            usesInsulatedWire = true;
        }
    }

    bool usesMarginTape = false;
    auto sectionsDescription = coil.get_sections_description().value();
    for (auto& section : sectionsDescription) {
        if (section.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        auto margins = Coil::resolve_margin(section);
        for (auto margin : margins) {
            if (margin > 0) {
                usesMarginTape = true;
            }
        }
    }

    double marginTapeTransition = rule_number("R7", "marginTapeTransitionEffectiveLength");
    double insulatedWireTransition = rule_number("R7", "insulatedWireTransitionEffectiveLength");

    if (effectiveLength >= marginTapeTransition && usesInsulatedWire && !usesMarginTape) {
        finding.set_status(ManufacturabilityStatus::WARNING);
        finding.set_threshold_value(marginTapeTransition);
        finding.set_message("This core is at or above the EE25/ER28 transition size (effective length " +
                            std::to_string(effectiveLength) + " m), where Wuerth achieve safety with magnet wire plus margin tape. The design uses insulated wire instead, which is billed by length and costs winding window.");
        return finding;
    }
    if (effectiveLength <= insulatedWireTransition && usesMarginTape && !usesInsulatedWire) {
        finding.set_status(ManufacturabilityStatus::WARNING);
        finding.set_threshold_value(insulatedWireTransition);
        finding.set_message("This core is at or below the EE20 size (effective length " + std::to_string(effectiveLength) +
                            " m), which is almost too small to achieve the safety distances with margin tape. Insulated wire is preferred there, to keep the copper area up.");
        return finding;
    }

    finding.set_status(ManufacturabilityStatus::PASS);
    finding.set_threshold_value(marginTapeTransition);
    finding.set_message("The insulation strategy in use matches Wuerth's size heuristic for a core of effective length " +
                        std::to_string(effectiveLength) + " m.");
    return finding;
}

// ---------------------------------------------------------------------------------------
// R8 - insulated-wire growth, placement and contact
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r8_insulated_wire(Magnetic& magnetic) {
    auto finding = make_finding("R8");
    finding.set_unit("outer-diameter growth over the enamelled equivalent, as a fraction");

    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil");
        finding.set_message("Insulated-wire rule not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto& coil = magnetic.get_mutable_coil();
    auto windings = coil.get_functional_description();

    std::vector<size_t> insulatedWindingIndexes;
    for (size_t index = 0; index < windings.size(); ++index) {
        auto wire = Coil::resolve_wire(windings[index]);
        if (is_insulated_wire(wire)) {
            insulatedWindingIndexes.push_back(index);
        }
    }
    if (insulatedWindingIndexes.empty()) {
        finding.set_status(ManufacturabilityStatus::NOT_APPLICABLE);
        finding.set_message("No winding uses insulated (extruded) wire, so neither the cost placement nor the contact rule applies.");
        return finding;
    }

    // The outer-diameter growth actually present, against the enamelled wire of the same conductor.
    double worstGrowth = 0;
    bool anyGrowthMeasured = false;
    for (auto index : insulatedWindingIndexes) {
        auto wire = Coil::resolve_wire(windings[index]);
        if (wire.get_type() != WireType::ROUND || !wire.get_conducting_diameter()) {
            continue;
        }
        double conductingDiameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
        if (!(conductingDiameter > 0)) {
            continue;
        }
        double enamelledEquivalent = Wire::get_outer_diameter_round(conductingDiameter, 1);
        if (!(enamelledEquivalent > 0)) {
            continue;
        }
        double growth = wire.get_maximum_outer_dimension() / enamelledEquivalent - 1.;
        worstGrowth = std::max(worstGrowth, growth);
        anyGrowthMeasured = true;
    }
    if (anyGrowthMeasured) {
        finding.set_measured_value(worstGrowth);
    }
    finding.set_threshold_value(rule_number("R8", "outerDiameterGrowthTripleInsulated"));

    std::vector<std::string> messages;
    if (anyGrowthMeasured) {
        messages.push_back("Insulated wire grows the outer diameter by up to " + std::to_string(worstGrowth * 100.) +
                           " % over the enamelled equivalent (Wuerth quote basic 15 %, supplementary 40 %, triple insulated 65 %).");
    }

    // Cost: insulated wire is billed by length, so it belongs on the lowest-turn winding.
    size_t lowestTurnIndex = 0;
    uint64_t lowestTurns = std::numeric_limits<uint64_t>::max();
    for (size_t index = 0; index < windings.size(); ++index) {
        auto turns = static_cast<uint64_t>(windings[index].get_number_turns());
        if (turns < lowestTurns) {
            lowestTurns = turns;
            lowestTurnIndex = index;
        }
    }
    bool misplaced = false;
    std::vector<std::string> misplacedWindings;
    for (auto index : insulatedWindingIndexes) {
        if (index != lowestTurnIndex) {
            misplaced = true;
            misplacedWindings.push_back(windings[index].get_name());
        }
    }
    if (misplaced) {
        messages.push_back("Insulated wire is billed by length, so it belongs on the lowest-turn winding ('" +
                           windings[lowestTurnIndex].get_name() + "', " + std::to_string(lowestTurns) +
                           " turns), not on " + join(misplacedWindings, ", ") + ".");
    }

    // Contact: insulated wire must not touch enamelled magnet wire without a tape layer.
    if (!coil.get_sections_description()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the coil has not been wound, so the insulated-to-enamelled contact check has no section order to walk");
        messages.push_back("Contact check not evaluated: " + finding.get_reason().value() + ".");
        finding.set_message(join(messages, " "));
        return finding;
    }

    std::vector<std::string> contactViolations;
    std::optional<std::string> previousConductionWinding;
    bool insulationSinceLastConduction = false;
    auto sectionsDescription = coil.get_sections_description().value();
    for (auto& section : sectionsDescription) {
        if (section.get_type() == ElectricalType::INSULATION) {
            insulationSinceLastConduction = true;
            continue;
        }
        if (section.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        auto windingName = section_winding_name(section);
        if (!windingName) {
            finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
            finding.set_reason("section '" + section.get_name() + "' names no winding, so its wire cannot be identified");
            messages.push_back("Contact check not evaluated: " + finding.get_reason().value() + ".");
            finding.set_message(join(messages, " "));
            return finding;
        }
        if (previousConductionWinding && !insulationSinceLastConduction &&
            previousConductionWinding.value() != windingName.value()) {
            auto previousWire = Coil::resolve_wire(coil.get_winding_by_name(previousConductionWinding.value()));
            auto currentWire = Coil::resolve_wire(coil.get_winding_by_name(windingName.value()));
            bool previousInsulated = is_insulated_wire(previousWire);
            bool currentInsulated = is_insulated_wire(currentWire);
            if ((previousInsulated && is_enamelled_magnet_wire(currentWire)) ||
                (currentInsulated && is_enamelled_magnet_wire(previousWire))) {
                contactViolations.push_back(previousConductionWinding.value() + " / " + windingName.value());
            }
        }
        previousConductionWinding = windingName;
        insulationSinceLastConduction = false;
    }

    if (!contactViolations.empty()) {
        finding.set_status(ManufacturabilityStatus::FAIL);
        finding.set_scope(join(contactViolations, ", "));
        messages.push_back("Insulated wire touches enamelled magnet wire with no tape layer between them (" +
                           join(contactViolations, ", ") + "); at least " +
                           std::to_string(static_cast<int>(rule_number("R8", "minimumTapeLayersBetweenInsulatedAndEnamelledWire"))) +
                           " tape layer is required.");
        finding.set_message(join(messages, " "));
        return finding;
    }
    if (misplaced) {
        finding.set_status(ManufacturabilityStatus::WARNING);
        finding.set_scope(join(misplacedWindings, ", "));
        finding.set_message(join(messages, " "));
        return finding;
    }

    finding.set_status(ManufacturabilityStatus::PASS);
    messages.push_back("Insulated wire sits on the lowest-turn winding and never touches enamelled magnet wire without tape.");
    finding.set_message(join(messages, " "));
    return finding;
}

// ---------------------------------------------------------------------------------------
// R9 - margin tape stack
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r9_margin_tape(Magnetic& magnetic) {
    auto finding = make_finding("R9");
    finding.set_unit("fraction of the winding window area lost to margin tape");

    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil");
        finding.set_message("Margin tape rule not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto& coil = magnetic.get_mutable_coil();
    if (!coil.get_sections_description()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the coil has not been wound, so it carries no sections and no margins");
        finding.set_message("Margin tape rule not evaluated: " + finding.get_reason().value());
        return finding;
    }

    auto sections = coil.get_sections_description().value();
    double marginArea = 0;
    size_t numberMarginedSections = 0;
    int64_t maximumLayersPerSide = 0;
    bool layerCountKnown = false;
    for (auto& section : sections) {
        if (section.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        if (section.get_coordinate_system() && section.get_coordinate_system().value() != CoordinateSystem::CARTESIAN) {
            finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
            finding.set_reason("section '" + section.get_name() +
                               "' is described in a polar coordinate system; this rule only measures the window loss of concentric (cartesian) coils");
            finding.set_message("Margin tape rule not evaluated: " + finding.get_reason().value());
            return finding;
        }
        auto margins = Coil::resolve_margin(section);
        double totalMargin = 0;
        for (auto margin : margins) {
            totalMargin += margin;
        }
        if (!(totalMargin > 0)) {
            continue;
        }
        ++numberMarginedSections;
        marginArea += totalMargin * section.get_dimensions()[0];
        if (section.get_margin() && std::holds_alternative<MarginInfo>(section.get_margin().value())) {
            // Only a MarginInfo margin knows how many tape layers it is; a bare [top, bottom]
            // pair does not, and this rule will not invent a number for it.
            layerCountKnown = true;
            maximumLayersPerSide = std::max(maximumLayersPerSide, std::get<MarginInfo>(section.get_margin().value()).get_number_layers());
        }
    }

    if (numberMarginedSections == 0) {
        finding.set_status(ManufacturabilityStatus::NOT_APPLICABLE);
        finding.set_message("No section carries margin tape.");
        return finding;
    }

    auto bobbin = coil.resolve_bobbin();
    if (!bobbin.get_processed_description() || bobbin.get_processed_description()->get_winding_windows().empty()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the bobbin has no processed winding window, so the fraction of window lost to margin tape cannot be computed");
        finding.set_message("Margin tape rule not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto windingWindow = bobbin.get_processed_description()->get_winding_windows()[0];
    if (!windingWindow.get_area()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the bobbin winding window carries no area");
        finding.set_message("Margin tape rule not evaluated: " + finding.get_reason().value());
        return finding;
    }
    double windowArea = windingWindow.get_area().value();
    if (!(windowArea > 0)) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the bobbin winding window area is not positive");
        finding.set_message("Margin tape rule not evaluated: " + finding.get_reason().value());
        return finding;
    }

    double lossRatio = marginArea / windowArea;
    finding.set_measured_value(lossRatio);

    double warnRatio = rule_number("R9", "windowWidthLossWarningRatio");
    double failRatio = rule_number("R9", "windowWidthLossFailRatio");

    std::vector<std::string> messages;
    messages.push_back("Margin tape on " + std::to_string(numberMarginedSections) + " section(s) consumes " +
                       std::to_string(lossRatio * 100.) + " % of the winding window area.");
    if (layerCountKnown) {
        messages.push_back("The thickest margin is " + std::to_string(maximumLayersPerSide) + " tape layers per side (Wuerth quote 10 to 20 per side per winding).");
        auto typicalMaximum = static_cast<int64_t>(rule_number("R9", "typicalLayersPerSideMaximum"));
        if (maximumLayersPerSide > typicalMaximum) {
            messages.push_back("That is above the 10 to 20 layers per side Wuerth quote as typical.");
        }
    }
    else {
        messages.push_back("The number of tape layers per side is not recorded: the sections carry bare margin widths rather than a marginInfo, and this rule does not infer a layer count.");
    }
    messages.push_back("The leakage-inductance delta caused by the margins is not reported here: MKF's leakage model would have to re-wind the coil without margins, which is an FEM-scale computation, not a report-time one.");

    if (lossRatio >= failRatio) {
        finding.set_status(ManufacturabilityStatus::FAIL);
        finding.set_threshold_value(failRatio);
    }
    else if (lossRatio >= warnRatio) {
        finding.set_status(ManufacturabilityStatus::WARNING);
        finding.set_threshold_value(warnRatio);
    }
    else {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_threshold_value(warnRatio);
    }
    finding.set_message(join(messages, " "));
    return finding;
}

// ---------------------------------------------------------------------------------------
// R11 - termination protection cost order
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r11_termination_protection(Magnetic& magnetic, Inputs& inputs) {
    auto finding = make_finding("R11");
    finding.set_unit("sleeved leads");
    std::vector<std::string> costOrder;
    for (auto& protection : rule("R11").at("protectionCostOrder")) {
        costOrder.push_back(protection.get<std::string>());
    }
    const std::string costOrderText = "Cost order, cheapest first: " + join(costOrder, " < ") + ".";

    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil, so there are no leads to protect");
        finding.set_message("Termination protection rule not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto& coil = magnetic.get_mutable_coil();
    if (!coil.get_sections_description()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the coil has not been wound, so it carries no margins for its leads to cross");
        finding.set_message("Termination protection rule not evaluated: " + finding.get_reason().value());
        return finding;
    }

    std::vector<std::string> sleevedLeads;
    for (auto& winding : coil.get_functional_description()) {
        for (End windingEnd : {End::START, End::FINISH}) {
            if (coil.get_recorded_lead_sleeve(winding.get_name(), windingEnd, 0)) {
                sleevedLeads.push_back(winding.get_name() + (windingEnd == End::START ? " start" : " finish"));
            }
        }
    }
    finding.set_measured_value(static_cast<double>(sleevedLeads.size()));

    if (!inputs.get_design_requirements().get_insulation()) {
        if (sleevedLeads.empty()) {
            finding.set_status(ManufacturabilityStatus::NOT_APPLICABLE);
            finding.set_message("The design has no insulation requirement and no lead is sleeved, so no termination protection is called for.");
            return finding;
        }
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_scope(join(sleevedLeads, ", "));
        finding.set_reason("leads are sleeved but the design states no insulation requirement to judge the sleeves against");
        finding.set_message("Termination protection rule not evaluated: " + finding.get_reason().value());
        return finding;
    }

    InsulationCoordinator coordinator;
    std::vector<std::string> unnecessarySleeves;
    std::vector<std::string> missingSleeves;
    std::vector<std::string> marginTapeWindings;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        auto windingName = coil.get_functional_description()[windingIndex].get_name();
        bool crossesMargin = coil.winding_leads_cross_margin(windingName);
        if (crossesMargin) {
            marginTapeWindings.push_back(windingName);
        }
        auto wire = coil.resolve_wire(windingIndex);
        bool sleeveRequired = coordinator.calculate_lead_sleeve_requirements(inputs, wire, crossesMargin).has_value();
        for (End windingEnd : {End::START, End::FINISH}) {
            const std::string lead = windingName + (windingEnd == End::START ? " start" : " finish");
            bool sleeved = coil.get_recorded_lead_sleeve(windingName, windingEnd, 0).has_value();
            if (sleeved && !sleeveRequired) {
                unnecessarySleeves.push_back(lead);
            }
            if (!sleeved && sleeveRequired) {
                missingSleeves.push_back(lead);
            }
        }
    }

    std::vector<std::string> messages;
    messages.push_back(costOrderText);
    if (!marginTapeWindings.empty()) {
        messages.push_back("Margin tape protects the windings " + join(marginTapeWindings, ", ") + ".");
    }
    if (!missingSleeves.empty() || !unnecessarySleeves.empty()) {
        finding.set_status(ManufacturabilityStatus::FAIL);
        std::vector<std::string> scope = missingSleeves;
        scope.insert(scope.end(), unnecessarySleeves.begin(), unnecessarySleeves.end());
        finding.set_scope(join(scope, ", "));
        if (!missingSleeves.empty()) {
            messages.push_back("These leads cross a margin unsleeved although their wire's own insulation does not cover the requirement, so the margin is bridged: " +
                               join(missingSleeves, ", ") + ".");
        }
        if (!unnecessarySleeves.empty()) {
            messages.push_back("These leads are sleeved although no requirement calls for it (no margin crossed, or the wire's own insulation covers it), spending the costliest protection for nothing: " +
                               join(unnecessarySleeves, ", ") + ".");
        }
        finding.set_message(join(messages, " "));
        return finding;
    }

    finding.set_status(ManufacturabilityStatus::PASS);
    if (sleevedLeads.empty()) {
        messages.push_back("No lead needs a sleeve.");
    }
    else {
        finding.set_scope(join(sleevedLeads, ", "));
        messages.push_back("Sleeving is used only where tape and the wire's own insulation cannot protect the lead. Manual labour: " +
                           std::to_string(sleevedLeads.size()) + " sleeved leads (" + join(sleevedLeads, ", ") + ").");
    }
    finding.set_message(join(messages, " "));
    return finding;
}

// ---------------------------------------------------------------------------------------
// R12 - flying leads
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r12_flying_leads(Magnetic& magnetic, Inputs& inputs) {
    auto finding = make_finding("R12");

    std::vector<std::string> flyingLeadWindings;
    if (magnetic.has_coil()) {
        for (auto& winding : magnetic.get_coil().get_functional_description()) {
            if (!winding.get_connections()) {
                continue;
            }
            auto connections = winding.get_connections().value();
            for (auto& connection : connections) {
                if (connection.get_type() && connection.get_type().value() == ConnectionType::FLYING_LEAD) {
                    flyingLeadWindings.push_back(winding.get_name());
                    break;
                }
            }
        }
    }

    bool requirementAsksForFlyingLeads = terminal_types_contain(inputs, ConnectionType::FLYING_LEAD);

    if (!requirementAsksForFlyingLeads && flyingLeadWindings.empty()) {
        if (!inputs.get_design_requirements().get_terminal_type() && !magnetic.has_coil()) {
            finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
            finding.set_reason("neither designRequirements.terminalType nor a coil connection says how the part is terminated");
            finding.set_message("Flying-lead check not evaluated: " + finding.get_reason().value());
            return finding;
        }
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_measured_value(0);
        finding.set_message("No termination is a flying lead.");
        return finding;
    }

    finding.set_status(ManufacturabilityStatus::WARNING);
    finding.set_measured_value(static_cast<double>(flyingLeadWindings.size()));
    if (!flyingLeadWindings.empty()) {
        finding.set_scope(join(flyingLeadWindings, ", "));
    }
    finding.set_message("Flying leads are Wuerth's least preferred termination: the leads are damaged by the remaining winding, soldering, taping, testing and varnishing steps and by PCB insertion, they must be inserted by hand, and being unshielded they emit and absorb noise." +
                        std::string(flyingLeadWindings.empty() ? " The design requirements ask for a flying-lead terminal."
                                                               : " Affected windings: " + join(flyingLeadWindings, ", ") + "."));
    return finding;
}

// ---------------------------------------------------------------------------------------
// R13 - manual termination
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r13_manual_termination(Magnetic& magnetic) {
    auto finding = make_finding("R13");
    finding.set_unit("m (conducting diameter)");

    if (!magnetic.has_coil()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no coil, so there is no wire to terminate");
        finding.set_message("Manual-termination check not evaluated: " + finding.get_reason().value());
        return finding;
    }

    double finestDiameterLimit = rule_number("R13", "finestAutomaticallySolderableConductingDiameter");
    finding.set_threshold_value(finestDiameterLimit);

    std::vector<std::string> reasons;
    double finestDiameter = std::numeric_limits<double>::max();
    bool anyDiameterKnown = false;
    for (auto& winding : magnetic.get_coil().get_functional_description()) {
        auto wire = Coil::resolve_wire(winding);
        if (wire.get_type() == WireType::LITZ) {
            reasons.push_back("'" + winding.get_name() + "' is litz, which needs crimp terminals or a special bobbin");
        }
        if (is_insulated_wire(wire)) {
            reasons.push_back("'" + winding.get_name() +
                              "' uses extruded insulation (TIW/TEX-E), which must be pre-stripped on a single-arbor machine");
        }
        // The conductor of a litz wire is its strand, whose fineness says nothing about how the
        // bundle is terminated; only solid wire is measured here.
        if (wire.get_type() == WireType::ROUND && wire.get_conducting_diameter()) {
            double conductingDiameter = resolve_dimensional_values(wire.get_conducting_diameter().value());
            anyDiameterKnown = true;
            finestDiameter = std::min(finestDiameter, conductingDiameter);
            if (conductingDiameter < finestDiameterLimit) {
                reasons.push_back("'" + winding.get_name() + "' is finer than " +
                                  std::to_string(static_cast<int>(rule_number("R13", "finestAutomaticallySolderableAwg"))) +
                                  " AWG (conducting diameter " + std::to_string(conductingDiameter) +
                                  " m), past the automatic-soldering limit");
            }
        }
    }
    if (anyDiameterKnown) {
        finding.set_measured_value(finestDiameter);
    }

    if (reasons.empty()) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_message("Every winding can be wound and terminated automatically on a multi-arbor machine.");
        return finding;
    }

    finding.set_status(ManufacturabilityStatus::WARNING);
    finding.set_message("Manual termination is needed: " + join(reasons, "; ") +
                        ". The solder bath runs at " + std::to_string(static_cast<int>(rule_number("R6", "solderBathTemperatureMinimum"))) +
                        " to " + std::to_string(static_cast<int>(rule_number("R6", "solderBathTemperatureMaximum"))) +
                        " degrees Celsius, so the dwell time a thermoplastic bobbin tolerates is limited.");
    return finding;
}

// ---------------------------------------------------------------------------------------
// R15 - additive gaps in the lateral columns
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r15_lateral_gaps(Magnetic& magnetic, Inputs& inputs) {
    auto finding = make_finding("R15");
    finding.set_unit("number of additive gaps outside the central column");
    finding.set_threshold_value(0);

    if (!inputs.get_design_requirements().get_topology()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("designRequirements has no topology, and this rule only applies to gapped isolated converters");
        finding.set_message("Lateral-gap EMI check not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto topologyName = std::string(magic_enum::enum_name(inputs.get_design_requirements().get_topology().value()));
    bool topologyChecked = false;
    for (auto& checked : rule("R15").at("topologiesChecked")) {
        if (checked.get<std::string>() == topologyName) {
            topologyChecked = true;
        }
    }
    if (!topologyChecked) {
        finding.set_status(ManufacturabilityStatus::NOT_APPLICABLE);
        finding.set_message("The topology " + topologyName + " is not one of the gapped isolated converters this rule covers.");
        return finding;
    }

    if (!magnetic.has_core()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the magnetic carries no core");
        finding.set_message("Lateral-gap EMI check not evaluated: " + finding.get_reason().value());
        return finding;
    }
    auto& core = magnetic.get_mutable_core();
    if (!core.get_processed_description()) {
        finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
        finding.set_reason("the core has no processed description, so its columns are unknown and a gap cannot be placed in one");
        finding.set_message("Lateral-gap EMI check not evaluated: " + finding.get_reason().value());
        return finding;
    }

    size_t numberLateralAdditiveGaps = 0;
    for (auto& column : core.get_columns()) {
        if (column.get_type() != ColumnType::LATERAL) {
            continue;
        }
        for (auto& gap : core.find_gaps_by_column(column)) {
            if (gap.get_type() == GapType::ADDITIVE) {
                ++numberLateralAdditiveGaps;
            }
        }
    }
    finding.set_measured_value(static_cast<double>(numberLateralAdditiveGaps));

    if (numberLateralAdditiveGaps == 0) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_message("Every additive gap is in the central column, where the winding shields it.");
        return finding;
    }

    finding.set_status(ManufacturabilityStatus::WARNING);
    finding.set_message("This " + topologyName + " has " + std::to_string(numberLateralAdditiveGaps) +
                        " additive gap(s) in a lateral column. An EI-style build with a spacer exposes three unshielded gaps and radiates; gap the centre leg only.");
    return finding;
}

// ---------------------------------------------------------------------------------------
// R17 - production test and tolerance fields (report only)
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r17_production_tests(Magnetic& magnetic, Inputs& inputs) {
    auto finding = make_finding("R17");
    finding.set_status(ManufacturabilityStatus::INFORMATIONAL);
    finding.set_unit("fraction (inductance tolerance)");

    std::vector<std::string> messages;
    messages.push_back("Production hipot test: " + std::to_string(static_cast<int>(rule_number("R17", "hipotDurationMinimum"))) +
                       " to " + std::to_string(static_cast<int>(rule_number("R17", "hipotDurationMaximum"))) +
                       " s, test voltage may be reduced by " + std::to_string(rule_number("R17", "hipotVoltageReduction") * 100.) + " %.");

    bool usesSpecialInsulatedWire = false;
    if (magnetic.has_coil()) {
        for (auto& winding : magnetic.get_coil().get_functional_description()) {
            auto wire = Coil::resolve_wire(winding);
            if (is_insulated_wire(wire) || is_fully_insulated_wire(wire)) {
                usesSpecialInsulatedWire = true;
            }
        }
    }
    double partialDischargeVoltage = rule_number("R17", "partialDischargeTestPeakVoltage");
    if (!usesSpecialInsulatedWire) {
        messages.push_back("Partial-discharge test: not required, the coil uses no FIW or TIW.");
    }
    else if (inputs.get_operating_points().empty()) {
        messages.push_back("Partial-discharge test: the coil uses FIW/TIW, but the peak working voltage is unknown (the inputs carry no operating point), so whether the " +
                           std::to_string(static_cast<int>(partialDischargeVoltage)) + " Vpk threshold is crossed is not evaluated.");
    }
    else {
        double peakVoltage = inputs.get_maximum_voltage_peak();
        finding.set_scope("peak working voltage " + std::to_string(peakVoltage) + " V");
        if (peakVoltage > partialDischargeVoltage) {
            messages.push_back("Partial-discharge test: REQUIRED. The coil uses FIW/TIW and the peak working voltage " +
                               std::to_string(peakVoltage) + " V is above " + std::to_string(static_cast<int>(partialDischargeVoltage)) + " Vpk.");
        }
        else {
            messages.push_back("Partial-discharge test: not required. The coil uses FIW/TIW but the peak working voltage " +
                               std::to_string(peakVoltage) + " V is at or below " + std::to_string(static_cast<int>(partialDischargeVoltage)) + " Vpk.");
        }
    }

    double typicalTolerance = rule_number("R17", "typicalInductanceTolerance");
    finding.set_threshold_value(typicalTolerance);
    auto magnetizingInductance = inputs.get_design_requirements().get_magnetizing_inductance();
    if (magnetizingInductance.get_nominal() && magnetizingInductance.get_maximum() &&
        magnetizingInductance.get_nominal().value() > 0) {
        double requestedTolerance = magnetizingInductance.get_maximum().value() / magnetizingInductance.get_nominal().value() - 1.;
        finding.set_measured_value(requestedTolerance);
        messages.push_back("Inductance tolerance requested: " + std::to_string(requestedTolerance * 100.) +
                           " % (Wuerth quote plus or minus " + std::to_string(typicalTolerance * 100.) + " % as typical).");
    }
    else {
        messages.push_back("Inductance tolerance requested: not specified as a nominal-and-maximum pair; Wuerth quote plus or minus " +
                           std::to_string(typicalTolerance * 100.) + " % as typical.");
    }

    finding.set_message(join(messages, " "));
    return finding;
}

// ---------------------------------------------------------------------------------------
// The whole report
// ---------------------------------------------------------------------------------------
ManufacturabilityReport Manufacturability::calculate_report(Magnetic& magnetic, Inputs& inputs) {
    ManufacturabilityReport report;
    report.add_finding(evaluate_r1_layer_parity(magnetic));
    report.add_finding(evaluate_r2_pin_order(magnetic));
    report.add_finding(evaluate_r3_rail_gauge_spread(magnetic));
    report.add_finding(evaluate_r4_wrap_height(magnetic));
    report.add_finding(evaluate_r5_smt_heavy_wire(magnetic, inputs));
    report.add_finding(evaluate_r6_thermoset_bobbin(magnetic, inputs));
    report.add_finding(evaluate_r7_size_heuristic(magnetic, inputs));
    report.add_finding(evaluate_r8_insulated_wire(magnetic));
    report.add_finding(evaluate_r9_margin_tape(magnetic));
    report.add_finding(not_evaluated_rule("R10"));
    report.add_finding(evaluate_r11_termination_protection(magnetic, inputs));
    report.add_finding(evaluate_r12_flying_leads(magnetic, inputs));
    report.add_finding(evaluate_r13_manual_termination(magnetic));
    report.add_finding(evaluate_r14_shield_terminations(magnetic));
    report.add_finding(evaluate_r15_lateral_gaps(magnetic, inputs));
    report.add_finding(not_evaluated_rule("R16"));
    report.add_finding(evaluate_r17_production_tests(magnetic, inputs));
    return report;
}

ManufacturabilityReport Manufacturability::calculate_report(Mas& mas) {
    return calculate_report(mas.get_mutable_magnetic(), mas.get_mutable_inputs());
}

} // namespace OpenMagnetics
