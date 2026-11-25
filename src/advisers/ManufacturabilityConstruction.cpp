// ABT #2548 (Alf, 2026-11-24, "MKF DFM rules + Henry playbook"): the construction rules R18 to R22.
//
// A design audit of 171 adviser and Henry designs (Alvier SMC study, 2026-11-24) found construction errors the
// rule pack R1..R17 let through: rectangular wire wound on edge in two or three layers, sections filled past what
// a winder can wind, toroid bores with no room left for the shuttle, overlapping turns, and solid wire many skin
// depths thick at a high ripple frequency. Every number these rules compare against is in src/data/dfm_rules.json
// with its source; the geometry and the physics are MKF's own (Coil::judge_copper_overlaps, the sections' stored
// filling factor, Wire::calculate_conducting_area, WindingSkinEffectLosses::calculate_skin_depth).

#include "advisers/Manufacturability.h"

#include "constructive_models/Coil.h"
#include "constructive_models/Core.h"
#include "constructive_models/Wire.h"
#include "physical_models/WindingSkinEffectLosses.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>

namespace OpenMagnetics {

namespace {

std::string join_items(const std::vector<std::string>& items) {
    std::string result;
    for (size_t index = 0; index < items.size(); ++index) {
        result += (index > 0 ? ", " : "") + items[index];
    }
    return result;
}

std::string format_number(double value, int precision = 3) {
    std::ostringstream stream;
    stream.precision(precision);
    stream << value;
    return stream.str();
}

ManufacturabilityFinding not_evaluated(ManufacturabilityFinding finding, const std::string& what, const std::string& reason) {
    finding.set_status(ManufacturabilityStatus::NOT_EVALUATED);
    finding.set_reason(reason);
    finding.set_message(what + " not evaluated: " + reason);
    return finding;
}

ManufacturabilityFinding not_applicable(ManufacturabilityFinding finding, const std::string& message) {
    finding.set_status(ManufacturabilityStatus::NOT_APPLICABLE);
    finding.set_message(message);
    return finding;
}

// The layers of a section that carry copper of the given winding.
std::vector<Layer> conduction_layers_of_winding(Coil& coil, const Section& section, const std::string& windingName) {
    std::vector<Layer> layers;
    for (auto& layer : coil.get_layers_by_section(section.get_name())) {
        if (layer.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        const auto& partials = layer.get_partial_windings();
        if (std::any_of(partials.begin(), partials.end(), [&](const PartialWinding& p) { return p.get_winding() == windingName; })) {
            layers.push_back(layer);
        }
    }
    return layers;
}

} // namespace

// ---------------------------------------------------------------------------------------
// R18 - flat wire in a single layer
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r18_flat_wire_single_layer(Magnetic& magnetic) {
    auto finding = make_finding("R18");
    const std::string what = "Flat-wire layer rule";
    if (!magnetic.has_coil()) {
        return not_evaluated(finding, what, "the magnetic carries no coil");
    }
    auto& coil = magnetic.get_mutable_coil();
    const auto rectangularLimit = static_cast<size_t>(rule_number("R18", "maximumLayersRectangular"));
    const auto foilTurnsPerLayerLimit = static_cast<size_t>(rule_number("R18", "maximumTurnsPerLayerFoil"));

    std::vector<std::string> flatWindings;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        auto wire = coil.resolve_wire(windingIndex);
        if (wire.get_type() == WireType::RECTANGULAR || wire.get_type() == WireType::FOIL) {
            flatWindings.push_back(coil.get_functional_description()[windingIndex].get_name());
        }
    }
    if (flatWindings.empty()) {
        return not_applicable(finding, "No winding uses rectangular or foil wire.");
    }
    if (!coil.get_sections_description() || !coil.get_layers_description()) {
        return not_evaluated(finding, what, "the coil has not been wound (no sections or layers description), so its layers are unknown");
    }
    if (!coil.get_turns_description() && std::any_of(flatWindings.begin(), flatWindings.end(), [&](const std::string& name) {
            return coil.resolve_wire(coil.get_winding_index_by_name(name)).get_type() == WireType::FOIL;
        })) {
        return not_evaluated(finding, what, "the coil carries no turns description, so the foil turns per layer are unknown");
    }

    std::vector<std::string> failures;
    double worst = 0;
    auto sections = coil.get_sections_description().value();
    const auto turns = coil.get_turns_description() ? coil.get_turns_description().value() : std::vector<Turn>{};
    for (const auto& windingName : flatWindings) {
        const auto windingIndex = coil.get_winding_index_by_name(windingName);
        const bool foil = coil.resolve_wire(windingIndex).get_type() == WireType::FOIL;
        for (const auto& section : sections) {
            if (section.get_type() != ElectricalType::CONDUCTION) {
                continue;
            }
            auto layers = conduction_layers_of_winding(coil, section, windingName);
            if (layers.empty()) {
                continue;
            }
            if (!foil) {
                if (layers.size() > rectangularLimit) {
                    failures.push_back(windingName + " in " + section.get_name() + ": " + std::to_string(layers.size()) + " layers of rectangular wire");
                    worst = std::max(worst, static_cast<double>(layers.size()));
                }
                continue;
            }
            // Foil: one turn of each parallel per layer is how foil is wound; more is a foil spiral in one layer.
            for (const auto& layer : layers) {
                std::map<int64_t, size_t> turnsPerParallel;
                for (const auto& turn : turns) {
                    if (turn.get_layer() && turn.get_layer().value() == layer.get_name() && turn.get_winding() == windingName) {
                        ++turnsPerParallel[turn.get_parallel()];
                    }
                }
                for (const auto& [parallel, count] : turnsPerParallel) {
                    if (count > foilTurnsPerLayerLimit) {
                        failures.push_back(windingName + " in " + layer.get_name() + ": " + std::to_string(count) +
                                           " foil turns of parallel " + std::to_string(parallel) + " side by side");
                        worst = std::max(worst, static_cast<double>(count));
                    }
                }
            }
        }
    }
    if (failures.empty()) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_message("Every rectangular-wire winding lies in a single layer, and every foil layer holds one turn per parallel.");
        return finding;
    }
    finding.set_status(ManufacturabilityStatus::FAIL);
    finding.set_measured_value(worst);
    finding.set_threshold_value(static_cast<double>(rectangularLimit));
    finding.set_scope(join_items(failures));
    finding.set_message("Rectangular wire must be wound in a single layer (either orientation), and foil one turn per layer: " +
                        join_items(failures) + ". Use fewer turns, a narrower strip, or round or litz wire for a multi-layer winding.");
    return finding;
}

// ---------------------------------------------------------------------------------------
// R19 - section filling factor
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r19_section_fill(Magnetic& magnetic) {
    auto finding = make_finding("R19");
    finding.set_unit("section filling factor");
    const std::string what = "Section fill rule";
    if (!magnetic.has_core()) {
        return not_evaluated(finding, what, "the magnetic carries no core, so it is unknown whether it is a toroid");
    }
    if (!magnetic.has_coil()) {
        return not_evaluated(finding, what, "the magnetic carries no coil");
    }
    auto& coil = magnetic.get_mutable_coil();
    if (!coil.get_sections_description()) {
        return not_evaluated(finding, what, "the coil has not been wound (no sections description)");
    }
    const bool toroid = magnetic.get_mutable_core().get_shape_family() == CoreShapeFamily::T;
    const double limit = rule_number("R19", toroid ? "maximumFillingFactorToroid" : "maximumFillingFactorOther");
    finding.set_threshold_value(limit);

    std::vector<std::string> over;
    std::vector<std::string> missing;
    double worst = 0;
    const auto sections = coil.get_sections_description().value();
    for (const auto& section : sections) {
        if (section.get_type() != ElectricalType::CONDUCTION) {
            continue;
        }
        if (!section.get_filling_factor()) {
            missing.push_back(section.get_name());
            continue;
        }
        const double fill = section.get_filling_factor().value();
        worst = std::max(worst, fill);
        if (fill > limit) {
            over.push_back(section.get_name() + " (" + format_number(fill * 100) + " %)");
        }
    }
    if (!missing.empty()) {
        return not_evaluated(finding, what, "the wound sections " + join_items(missing) + " store no filling factor");
    }
    finding.set_measured_value(worst);
    if (over.empty()) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_message("Every section is filled to at most " + format_number(limit * 100) + " % (" +
                            (toroid ? "toroid" : "bobbin core") + "); the fullest is at " + format_number(worst * 100) + " %.");
        return finding;
    }
    finding.set_status(ManufacturabilityStatus::FAIL);
    finding.set_scope(join_items(over));
    finding.set_message(std::string("Sections filled past the ") + format_number(limit * 100) + " % a " +
                        (toroid ? "toroid" : "bobbin core") + " can be wound to: " + join_items(over) +
                        ". Use fewer turns, a thinner wire or a core with a larger winding window.");
    return finding;
}

// ---------------------------------------------------------------------------------------
// R20 - toroid copper over the bore area
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r20_toroid_bore_copper(Magnetic& magnetic) {
    auto finding = make_finding("R20");
    finding.set_unit("copper area / bore area");
    const std::string what = "Toroid bore copper rule";
    if (!magnetic.has_core()) {
        return not_evaluated(finding, what, "the magnetic carries no core, so it is unknown whether it is a toroid");
    }
    auto& core = magnetic.get_mutable_core();
    if (core.get_shape_family() != CoreShapeFamily::T) {
        return not_applicable(finding, "The core is not a toroid.");
    }
    if (!magnetic.has_coil()) {
        return not_evaluated(finding, what, "the magnetic carries no coil");
    }
    const auto windingWindows = core.get_winding_windows();
    if (windingWindows.empty() || !windingWindows[0].get_area()) {
        return not_evaluated(finding, what, "the toroid's processed description has no winding window area (its bore)");
    }
    const double boreArea = windingWindows[0].get_area().value();
    auto& coil = magnetic.get_mutable_coil();
    double copperArea = 0;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        const auto& winding = coil.get_functional_description()[windingIndex];
        auto wire = coil.resolve_wire(windingIndex);
        // Every turn crosses the bore once: turns x parallels conductors of the wire's copper.
        copperArea += static_cast<double>(winding.get_number_turns() * winding.get_number_parallels()) * wire.calculate_conducting_area();
    }
    const double ratio = copperArea / boreArea;
    const double limit = rule_number("R20", "maximumCopperOverBoreArea");
    finding.set_measured_value(ratio);
    finding.set_threshold_value(limit);
    if (ratio <= limit) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_message("The copper through the bore takes " + format_number(ratio * 100) + " % of its area (at most " +
                            format_number(limit * 100) + " %), leaving room for the winding shuttle.");
        return finding;
    }
    finding.set_status(ManufacturabilityStatus::FAIL);
    finding.set_message("The copper through the bore takes " + format_number(ratio * 100) + " % of its area, more than the " +
                        format_number(limit * 100) + " % a toroid winder leaves for its shuttle. Use fewer turns, a thinner wire or a toroid with a larger bore.");
    return finding;
}

// ---------------------------------------------------------------------------------------
// R21 - no turn overlap on the real cross-sections
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r21_turn_overlap(Magnetic& magnetic) {
    auto finding = make_finding("R21");
    finding.set_unit("m");
    const std::string what = "Turn overlap rule";
    if (!magnetic.has_coil()) {
        return not_evaluated(finding, what, "the magnetic carries no coil");
    }
    auto& coil = magnetic.get_mutable_coil();
    if (!coil.get_turns_description()) {
        return not_evaluated(finding, what, "the coil carries no turns description, so there are no turns to place");
    }
    if (auto reason = coil.turn_outlines_unknown_reason()) {
        return not_evaluated(finding, what, reason.value());
    }
    const double limit = rule_number("R21", "maximumTurnOverlap");
    finding.set_threshold_value(limit);
    // The turns' real cross-sections at every station (Coil::judge_copper_overlaps, the overlaps the Painter draws in
    // red and ResultsAnalyzer reports): judged without the terminal leads, which this rule does not cover.
    const auto overlaps = coil.judge_copper_overlaps({}, {});
    double worst = 0;
    size_t count = 0;
    const CopperOverlap* worstOverlap = nullptr;
    for (const auto& overlap : overlaps) {
        if (overlap.kind != CopperOverlapKind::TURN_TURN || overlap.overlap <= limit) {
            continue;
        }
        ++count;
        if (overlap.overlap > worst) {
            worst = overlap.overlap;
            worstOverlap = &overlap;
        }
    }
    finding.set_measured_value(worst);
    if (count == 0) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_message("No two turns overlap by more than " + format_number(limit * 1e6) + " um on their real cross-sections.");
        return finding;
    }
    finding.set_status(ManufacturabilityStatus::FAIL);
    finding.set_scope(worstOverlap->first + " / " + worstOverlap->second);
    finding.set_message(std::to_string(count) + " pair" + (count == 1 ? "" : "s") + " of turn stations overlap by more than " +
                        format_number(limit * 1e6) + " um; the worst, " + worstOverlap->first + " and " + worstOverlap->second +
                        ", by " + format_number(worst * 1e3) + " mm. The turns do not fit where they are laid.");
    return finding;
}

// ---------------------------------------------------------------------------------------
// R22 - solid conductor much thicker than the skin depth carrying a large AC share
// ---------------------------------------------------------------------------------------
ManufacturabilityFinding Manufacturability::evaluate_r22_thick_solid_conductor(Magnetic& magnetic, Inputs& inputs, const std::vector<Outputs>& outputs) {
    auto finding = make_finding("R22");
    finding.set_unit("conductor thickness / skin depth at the AC loss's frequency");
    const std::string what = "Thick solid conductor rule";
    if (!magnetic.has_coil()) {
        return not_evaluated(finding, what, "the magnetic carries no coil");
    }
    auto& coil = magnetic.get_mutable_coil();
    std::vector<size_t> solidWindings;
    for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); ++windingIndex) {
        if (coil.resolve_wire(windingIndex).get_type() != WireType::LITZ) {
            solidWindings.push_back(windingIndex);
        }
    }
    if (solidWindings.empty()) {
        return not_applicable(finding, "Every winding is litz wire.");
    }
    const auto& operatingPoints = inputs.get_operating_points();
    if (outputs.empty()) {
        return not_evaluated(finding, what, "the design carries no outputs, so the winding-loss split is unknown");
    }
    if (outputs.size() != operatingPoints.size()) {
        return not_evaluated(finding, what, "the design carries " + std::to_string(outputs.size()) + " outputs for " +
                                                std::to_string(operatingPoints.size()) + " operating points");
    }
    const double thicknessLimit = rule_number("R22", "maximumConductorOverSkinDepth");
    const double shareLimit = rule_number("R22", "maximumAcShareOfWindingLoss");
    finding.set_threshold_value(thicknessLimit);

    std::vector<std::string> warnings;
    double worstRatio = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < outputs.size(); ++operatingPointIndex) {
        const auto& windingLosses = outputs[operatingPointIndex].get_winding_losses();
        if (!windingLosses || !windingLosses->get_winding_losses_per_winding() || !windingLosses->get_temperature()) {
            return not_evaluated(finding, what, "output " + std::to_string(operatingPointIndex) +
                                                    " has no winding losses per winding with their temperature");
        }
        const auto perWinding = windingLosses->get_winding_losses_per_winding().value();
        const double temperature = windingLosses->get_temperature().value();
        const auto& excitations = operatingPoints[operatingPointIndex].get_excitations_per_winding();
        for (auto windingIndex : solidWindings) {
            const auto& windingName = coil.get_functional_description()[windingIndex].get_name();
            auto losses = std::find_if(perWinding.begin(), perWinding.end(), [&](const WindingLossesPerElement& element) {
                return element.get_name() && element.get_name().value() == windingName;
            });
            if (losses == perWinding.end()) {
                return not_evaluated(finding, what, "output " + std::to_string(operatingPointIndex) + " has no winding losses named '" + windingName + "'");
            }
            if (!losses->get_ohmic_losses() || (!losses->get_skin_effect_losses() && !losses->get_proximity_effect_losses())) {
                return not_evaluated(finding, what, "output " + std::to_string(operatingPointIndex) + " splits no ohmic and AC losses for '" + windingName + "'");
            }
            if (windingIndex >= excitations.size()) {
                return not_evaluated(finding, what, "operating point " + std::to_string(operatingPointIndex) + " has no excitation for '" + windingName + "'");
            }
            auto wire = coil.resolve_wire(windingIndex);
            double thickness = 0;
            if (wire.get_type() == WireType::ROUND) {
                if (!wire.get_conducting_diameter()) {
                    return not_evaluated(finding, what, "the round wire of '" + windingName + "' states no conducting diameter");
                }
                thickness = resolve_dimensional_values(wire.get_conducting_diameter().value());
            }
            else {
                if (!wire.get_conducting_width() || !wire.get_conducting_height()) {
                    return not_evaluated(finding, what, "the wire of '" + windingName + "' states no conducting width and height");
                }
                thickness = std::min(resolve_dimensional_values(wire.get_conducting_width().value()),
                                     resolve_dimensional_values(wire.get_conducting_height().value()));
            }
            // The frequency the AC loss happens at: MKF's skin and proximity losses weighted over their harmonics. A PFC
            // choke's excitation frequency is the line's (50 Hz) while its AC loss sits at the switching ripple, so the
            // excitation frequency would never see the ripple (the Alvier audit judged at the switching frequency).
            const double ohmic = losses->get_ohmic_losses()->get_losses();
            double alternating = 0;
            double lossTimesFrequency = 0;
            for (const auto& element : {losses->get_skin_effect_losses(), losses->get_proximity_effect_losses()}) {
                if (!element) {
                    continue;
                }
                const auto& perHarmonic = element->get_losses_per_harmonic();
                const auto& frequencies = element->get_harmonic_frequencies();
                if (frequencies.size() != perHarmonic.size()) {
                    return not_evaluated(finding, what, "output " + std::to_string(operatingPointIndex) + " gives " +
                                                            std::to_string(perHarmonic.size()) + " AC losses of '" + windingName + "' for " +
                                                            std::to_string(frequencies.size()) + " harmonic frequencies");
                }
                for (size_t harmonic = 0; harmonic < perHarmonic.size(); ++harmonic) {
                    alternating += perHarmonic[harmonic];
                    lossTimesFrequency += perHarmonic[harmonic] * frequencies[harmonic];
                }
            }
            const double total = ohmic + alternating;
            if (total <= 0 || alternating <= 0) {
                continue;   // no AC loss: nothing for this rule to judge
            }
            const double frequency = lossTimesFrequency / alternating;
            const double skinDepth = WindingSkinEffectLosses::calculate_skin_depth(wire, frequency, temperature);
            const double ratio = thickness / skinDepth;
            const double share = alternating / total;
            worstRatio = std::max(worstRatio, ratio);
            if (ratio > thicknessLimit && share > shareLimit) {
                warnings.push_back(windingName + " at operating point " + std::to_string(operatingPointIndex) + ": " +
                                   format_number(thickness * 1e3) + " mm is " + format_number(ratio) + " skin depths at " +
                                   format_number(frequency / 1e3) + " kHz, AC share " + format_number(share * 100) + " %");
            }
        }
    }
    finding.set_measured_value(worstRatio);
    if (warnings.empty()) {
        finding.set_status(ManufacturabilityStatus::PASS);
        finding.set_message("No solid conductor is more than " + format_number(thicknessLimit) + " skin depths thick while skin and proximity carry more than " +
                            format_number(shareLimit * 100) + " % of its winding loss.");
        return finding;
    }
    finding.set_status(ManufacturabilityStatus::WARNING);
    finding.set_scope(join_items(warnings));
    finding.set_message("Solid conductor much thicker than the skin depth, with most of its loss from AC: " + join_items(warnings) +
                        ". Prefer litz wire (or a thinner conductor in more parallels).");
    return finding;
}

} // namespace OpenMagnetics
