#include "processors/ResultsAnalyzer.h"
#include "constructive_models/Core.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/Wire.h"
#include "physical_models/InitialPermeability.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/MagneticSimulator.h"
#include "support/Exceptions.h"
#include "support/Utils.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace OpenMagnetics {

std::string to_string(ResultFindingCode code) {
    switch (code) {
        case ResultFindingCode::CORE_SATURATES: return "CORE_SATURATES";
        case ResultFindingCode::VOLT_SECOND_IMBALANCE: return "VOLT_SECOND_IMBALANCE";
        case ResultFindingCode::CORE_TEMPERATURE_ABOVE_LIMIT: return "CORE_TEMPERATURE_ABOVE_LIMIT";
        case ResultFindingCode::CORE_TEMPERATURE_LIMIT_NOT_STATED: return "CORE_TEMPERATURE_LIMIT_NOT_STATED";
        case ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT: return "CORE_ABOVE_MATERIAL_DATA_LIMIT";
        case ResultFindingCode::WINDING_TEMPERATURE_ABOVE_LIMIT: return "WINDING_TEMPERATURE_ABOVE_LIMIT";
        case ResultFindingCode::WINDING_TEMPERATURE_LIMIT_NOT_STATED: return "WINDING_TEMPERATURE_LIMIT_NOT_STATED";
        case ResultFindingCode::NO_THERMAL_OPERATING_POINT: return "NO_THERMAL_OPERATING_POINT";
    }
    throw InvalidInputException(ErrorCode::INVALID_INPUT, "Unknown result finding code " + std::to_string(static_cast<int>(code)));
}

void to_json(nlohmann::json& j, const ResultFinding& finding) {
    j = nlohmann::json::object();
    j["operatingPointIndex"] = finding.operatingPointIndex;
    if (finding.operatingPointName) {
        j["operatingPointName"] = finding.operatingPointName.value();
    }
    j["code"] = to_string(finding.code);
    j["message"] = finding.message;
}

void to_json(nlohmann::json& j, const std::vector<ResultFinding>& findings) {
    j = nlohmann::json::array();
    for (const auto& finding : findings) {
        nlohmann::json entry;
        to_json(entry, finding);
        j.push_back(entry);
    }
}

namespace {

// ABT #1908: the ceiling the loss <-> temperature solve kept this operating point's core below
// (MagneticSimulator::get_core_data_ceiling), on the excitation it priced: the operating point's first
// winding with the flux density the output's core losses carry, from the point's ambient up, with the
// core-loss models the Settings select (as the simulator's).
std::optional<CoreDataCeiling> get_solved_core_data_ceiling(const Magnetic& magnetic, const OperatingPoint& operatingPoint,
                                                           const CoreLossesOutput& coreLosses) {
    if (operatingPoint.get_excitations_per_winding().empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "analyze_results: an operating point has no excitation");
    }
    if (!coreLosses.get_magnetic_flux_density()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "analyze_results: the core losses carry no magnetic flux density");
    }
    OperatingPointExcitation excitation = operatingPoint.get_excitations_per_winding()[0];
    excitation.set_magnetic_flux_density(coreLosses.get_magnetic_flux_density().value());
    CoreLosses coreLossesModel;
    return MagneticSimulator::get_core_data_ceiling(magnetic.get_core(), excitation,
                                                    operatingPoint.get_conditions().get_ambient_temperature(), coreLossesModel);
}

std::string describe_operating_point(const OperatingPoint& operatingPoint, size_t operatingPointIndex) {
    std::string description = "operating point " + std::to_string(operatingPointIndex);
    if (operatingPoint.get_name()) {
        description += " ('" + operatingPoint.get_name().value() + "')";
    }
    return description;
}

ResultFinding make_finding(const OperatingPoint& operatingPoint, size_t operatingPointIndex, ResultFindingCode code, const std::string& message) {
    ResultFinding finding;
    finding.operatingPointIndex = operatingPointIndex;
    finding.operatingPointName = operatingPoint.get_name();
    finding.code = code;
    finding.message = message;
    return finding;
}

std::string format_number(double value) {
    std::ostringstream stream;
    stream << value;
    return stream.str();
}

// A temperature limit a part states, and where it states it.
struct StatedTemperatureLimit {
    double temperature;
    std::string source;
};

} // namespace

SaturationVerdict ResultsAnalyzer::check_saturation(const CoreMaterial& coreMaterial, double peakFluxDensity, double temperature, double margin,
                                                    std::optional<DcBiasPoint> dcBiasPoint) {
    if (!std::isfinite(peakFluxDensity)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "Saturation check: the peak flux density must be finite, not " + format_number(peakFluxDensity) + " T");
    }
    if (!(margin >= 1)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "Saturation check: the margin multiplies the peak flux density and must be at least 1, not " + format_number(margin));
    }
    SaturationVerdict verdict;
    verdict.peakFluxDensity = peakFluxDensity;
    verdict.temperature = temperature;
    verdict.margin = margin;
    verdict.saturationFluxDensity = Core::get_magnetic_flux_density_saturation(coreMaterial, temperature, false);
    verdict.fluxDensityAboveSaturation = peakFluxDensity * margin > verdict.saturationFluxDensity;
    verdict.permeabilityCollapsed = false;

    // Criterion 2 (Alf, 2026-10-14): the material's permeability model, at the DC point the magnetizing
    // model converged on, has fallen to or below the vacuum's while the same material without the bias is
    // above it: the bias has driven the material past its knee and the magnetizing model took mu_r = 1
    // (MagnetizingInductance, InitialPermeability::get_initial_permeability). Until ABT #1855 this was the
    // throw condition (SaturationException). It is needed beside criterion 1 because a material's B_sat,
    // its magnetisation curve and its permeability model need not agree: 3C95 under DC bias goes air-core
    // from 44.6 A while its curve's B stays at 0.49 T against B_sat 0.53 T. The margin scales B and does
    // not apply here; the permeability is read at the bias point's own temperature and frequency.
    std::optional<double> unbiasedPermeability;
    if (dcBiasPoint && dcBiasPoint->magneticFieldDcBias > 0) {
        verdict.biasedPermeability = InitialPermeability::get_initial_permeability_model_value(
            coreMaterial, dcBiasPoint->temperature, dcBiasPoint->magneticFieldDcBias, dcBiasPoint->frequency);
        if (verdict.biasedPermeability.value() <= 1) {
            unbiasedPermeability = InitialPermeability::get_initial_permeability_model_value(
                coreMaterial, dcBiasPoint->temperature, std::nullopt, dcBiasPoint->frequency);
            verdict.permeabilityCollapsed = unbiasedPermeability.value() > 1;
        }
    }
    verdict.saturates = verdict.fluxDensityAboveSaturation || verdict.permeabilityCollapsed;

    std::ostringstream message;
    if (verdict.saturates) {
        message << "The core saturates:";
        if (verdict.fluxDensityAboveSaturation) {
            message << " B peak = " << peakFluxDensity << " T";
            if (margin != 1) {
                message << " (x " << margin << " margin = " << peakFluxDensity * margin << " T)";
            }
            message << " is above B_sat = " << verdict.saturationFluxDensity << " T of material " << coreMaterial.get_name()
                    << " at " << temperature << " C";
        }
        if (verdict.permeabilityCollapsed) {
            if (verdict.fluxDensityAboveSaturation) {
                message << ", and";
            }
            message << " at the DC bias H = " << dcBiasPoint->magneticFieldDcBias << " A/m (" << dcBiasPoint->temperature
                    << " C, " << dcBiasPoint->frequency << " Hz) the permeability model of material " << coreMaterial.get_name()
                    << " gives mu_r = " << verdict.biasedPermeability.value() << ", at or below the vacuum's 1 (unbiased "
                    << unbiasedPermeability.value() << "): the bias has driven it past its knee, and the model takes mu_r = 1 there";
            if (!verdict.fluxDensityAboveSaturation) {
                message << " (B peak " << peakFluxDensity << " T is still within B_sat = " << verdict.saturationFluxDensity
                        << " T at " << temperature << " C: the material's B_sat, magnetisation curve and permeability model disagree)";
            }
        }
        message << ". Numbers computed for this operating point are those of the saturated core;"
                << " more turns, a larger gap or core, or less volt-seconds are needed";
    }
    else {
        message << "B peak = " << peakFluxDensity << " T";
        if (margin != 1) {
            message << " (x " << margin << " margin = " << peakFluxDensity * margin << " T)";
        }
        message << " is within B_sat = " << verdict.saturationFluxDensity << " T of material " << coreMaterial.get_name()
                << " at " << temperature << " C";
    }
    verdict.message = message.str();
    return verdict;
}

double ResultsAnalyzer::get_peak_flux_density(const SignalDescriptor& magneticFluxDensity) {
    if (magneticFluxDensity.get_processed() && magneticFluxDensity.get_processed()->get_peak()) {
        return std::fabs(magneticFluxDensity.get_processed()->get_peak().value());
    }
    if (magneticFluxDensity.get_waveform() && !magneticFluxDensity.get_waveform()->get_data().empty()) {
        // get_waveform() returns the optional by value: bind it before reading its data.
        auto waveform = magneticFluxDensity.get_waveform().value();
        const auto& data = waveform.get_data();
        return std::max(std::fabs(*std::max_element(data.begin(), data.end())), std::fabs(*std::min_element(data.begin(), data.end())));
    }
    throw InvalidInputException(ErrorCode::MISSING_DATA, "The magnetic flux density has neither a processed peak nor a waveform");
}

std::vector<ResultFinding> ResultsAnalyzer::analyze_saturation(const Magnetic& magnetic, const OperatingPoint& operatingPoint,
                                                               const Outputs& output, size_t operatingPointIndex,
                                                               std::optional<Topology> topology) {
    std::string where = describe_operating_point(operatingPoint, operatingPointIndex);
    if (!output.get_core_losses()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "analyze_results: the output of " + where + " has no core losses, which carry its flux density");
    }
    auto coreLosses = output.get_core_losses().value();
    if (!coreLosses.get_magnetic_flux_density()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "analyze_results: the core losses of " + where + " carry no magnetic flux density");
    }
    double peakFluxDensity = get_peak_flux_density(coreLosses.get_magnetic_flux_density().value());
    const auto coreMaterial = magnetic.get_core().resolve_material();
    const double coreTemperature = coreLosses.get_temperature();
    // The DC point of criterion 2: the result does not carry it, so the magnetizing model is run again at the
    // temperature the core losses were evaluated at, as the loss solve ran it. A core the loss solve held at
    // its Curie temperature was magnetized there as a paramagnetic core (mu_r = 1, MagneticSimulator::
    // make_paramagnetic_magnetic), which has no knee to be driven past: only criterion 1 applies to it.
    std::optional<DcBiasPoint> dcBiasPoint;
    const auto dataCeiling = get_solved_core_data_ceiling(magnetic, operatingPoint, coreLosses);
    const bool heldAtCurie = dataCeiling && dataCeiling->source == CoreDataCeilingSource::CURIE_TEMPERATURE &&
                             coreTemperature > dataCeiling->temperature - MagneticSimulator::kOperatingTemperatureTolerance;
    if (!heldAtCurie) {
        OperatingPoint operatingPointAtCoreTemperature = operatingPoint;
        operatingPointAtCoreTemperature.get_mutable_conditions().set_ambient_temperature(coreTemperature);
        MagnetizingInductance magnetizingInductanceModel;
        magnetizingInductanceModel.set_topology(topology);
        magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(magnetic, &operatingPointAtCoreTemperature);
        dcBiasPoint = magnetizingInductanceModel.get_last_dc_bias_point();
    }
    auto verdict = check_saturation(coreMaterial, peakFluxDensity, coreTemperature, 1.0, dcBiasPoint);
    if (!verdict.saturates) {
        return {};
    }
    return {make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::CORE_SATURATES, "At " + where + ": " + verdict.message)};
}

std::vector<ResultFinding> ResultsAnalyzer::analyze_volt_second_balance(const OperatingPoint& operatingPoint, size_t operatingPointIndex) {
    std::vector<ResultFinding> findings;
    const auto& excitations = operatingPoint.get_excitations_per_winding();
    for (size_t windingIndex = 0; windingIndex < excitations.size(); ++windingIndex) {
        const auto& excitation = excitations[windingIndex];
        if (!excitation.get_voltage() || !excitation.get_voltage()->get_waveform()) {
            continue;
        }
        double frequency = excitation.get_frequency();
        auto voltage = Inputs::standardize_waveform(excitation.get_voltage().value(), frequency);
        auto sampledVoltage = Inputs::calculate_sampled_waveform(voltage.get_waveform().value(), frequency);
        std::string context = "At " + describe_operating_point(operatingPoint, operatingPointIndex) + ", winding " +
                              std::to_string(windingIndex) + " ('" + excitation.get_name().value_or("unnamed") + "')";
        auto imbalance = Inputs::describe_volt_second_imbalance(sampledVoltage, context);
        if (imbalance) {
            findings.push_back(make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::VOLT_SECOND_IMBALANCE,
                imbalance.value() + ". Where MKF derives a magnetizing current from this voltage it integrated its "
                "balanced part (the mean removed), so these numbers describe that part"));
        }
    }
    return findings;
}

std::vector<ResultFinding> ResultsAnalyzer::analyze_temperature(const Magnetic& magnetic, const OperatingPoint& operatingPoint,
                                                                const Outputs& output, size_t operatingPointIndex) {
    std::vector<ResultFinding> findings;
    std::string where = describe_operating_point(operatingPoint, operatingPointIndex);

    // The core: its temperature (the one its losses were evaluated at) against the limits its material states.
    if (!output.get_core_losses()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "analyze_results: the output of " + where + " has no core losses, which carry the core's temperature");
    }
    double coreTemperature = output.get_core_losses()->get_temperature();
    auto coreMaterial = magnetic.get_core().resolve_material();
    if (!output.get_temperature()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "analyze_results: the output of " + where + " has no temperature (simulated in fast mode?)");
    }
    double hottestTemperature = output.get_temperature()->get_maximum_temperature();

    // The core held at its material's data limit by the loss solve (ABT #1855, Alf 2026-10-14).
    // The ceiling covers the mu_i(T) / Curie data, the saturation field strength and the evaluated Steinmetz
    // ranges' ct(T) (ABT #1908).
    const auto dataCeiling = get_solved_core_data_ceiling(magnetic, operatingPoint, output.get_core_losses().value());
    if (MagneticSimulator::is_core_held_at_data_ceiling(dataCeiling, coreTemperature, hottestTemperature)) {
        double ceiling = dataCeiling->temperature;
        std::string limit;
        std::string magnetized;
        switch (dataCeiling->source) {
            case CoreDataCeilingSource::CURIE_TEMPERATURE:
                limit = "the Curie temperature of " + coreMaterial.get_name();
                magnetized = "as a paramagnetic core (mu_r = 1)";
                break;
            case CoreDataCeilingSource::PERMEABILITY_TABLE_END:
                limit = "the last temperature of " + coreMaterial.get_name() + "'s measured permeability data";
                magnetized = "with its last tabulated permeability";
                break;
            case CoreDataCeilingSource::SATURATION_FIELD_STRENGTH_END:
                limit = "the temperature where " + coreMaterial.get_name() + "'s saturation field strength data reach zero";
                magnetized = "with its own permeability at that temperature";
                break;
            case CoreDataCeilingSource::STEINMETZ_TEMPERATURE_FACTOR_END:
                limit = "the temperature where the Steinmetz temperature factor ct(T) of " + coreMaterial.get_name() +
                        "'s loss data at this frequency reaches zero";
                magnetized = "with its own permeability at that temperature";
                break;
        }
        findings.push_back(make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT,
            "At " + where + ": the losses drive the core above " + limit +
            " (" + format_number(ceiling) + " C), where its material data ends. Its core losses were evaluated at " +
            format_number(coreTemperature) + " C, the last temperature the data covers, with the core magnetized " +
            magnetized + "; the temperatures given for this operating point (hottest point " + format_number(hottestTemperature) +
            " C) follow from those losses, and the core's real temperature is higher"));
    }
    std::vector<StatedTemperatureLimit> coreLimits;
    if (coreMaterial.get_curie_temperature()) {
        coreLimits.push_back({coreMaterial.get_curie_temperature().value(), "its Curie temperature"});
    }
    if (coreMaterial.get_recommendations() && coreMaterial.get_recommendations()->get_maximum_operating_temperature()) {
        coreLimits.push_back({coreMaterial.get_recommendations()->get_maximum_operating_temperature().value(),
                              "its maximum recommended operating temperature"});
    }
    if (coreLimits.empty()) {
        findings.push_back(make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::CORE_TEMPERATURE_LIMIT_NOT_STATED,
            "At " + where + ": the core runs at " + format_number(coreTemperature) + " C, and material " + coreMaterial.get_name() +
            " states no temperature limit (no Curie temperature, no maximum operating temperature): its temperature is not checked"));
    }
    else {
        std::vector<std::string> exceeded;
        for (const auto& limit : coreLimits) {
            if (coreTemperature > limit.temperature) {
                exceeded.push_back(limit.source + " of " + format_number(limit.temperature) + " C");
            }
        }
        if (!exceeded.empty()) {
            std::string message = "At " + where + ": the core runs at " + format_number(coreTemperature) + " C, above material " +
                                  coreMaterial.get_name() + "'s ";
            for (size_t index = 0; index < exceeded.size(); ++index) {
                message += (index > 0 ? " and " : "") + exceeded[index];
            }
            findings.push_back(make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::CORE_TEMPERATURE_ABOVE_LIMIT, message));
        }
    }

    // The windings: the design's hottest point against the limits each winding's wire coating states.
    const auto& windings = magnetic.get_coil().get_functional_description();
    for (size_t windingIndex = 0; windingIndex < windings.size(); ++windingIndex) {
        auto wire = Coil::resolve_wire(windings[windingIndex]);
        std::string windingLabel = "winding " + std::to_string(windingIndex) + " ('" + windings[windingIndex].get_name() + "')";
        std::vector<StatedTemperatureLimit> windingLimits;
        auto coating = Wire::resolve_coating(wire);
        if (coating) {
            if (coating->get_temperature_rating()) {
                windingLimits.push_back({coating->get_temperature_rating().value(), "its wire coating's temperature rating"});
            }
            // Only a material the coating names: an unnamed one is not assumed (Wire::resolve_coating_insulation_material
            // would take the default enamel).
            if (coating->get_material()) {
                auto insulationMaterial = Wire::resolve_coating_insulation_material(coating.value());
                if (insulationMaterial.get_temperature_class()) {
                    windingLimits.push_back({InsulationCoordinator::temperature_class_to_celsius(insulationMaterial.get_temperature_class().value()),
                                             "the thermal class of its coating material " + insulationMaterial.get_name()});
                }
            }
        }
        if (windingLimits.empty()) {
            findings.push_back(make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::WINDING_TEMPERATURE_LIMIT_NOT_STATED,
                "At " + where + ": the design's hottest point is " + format_number(hottestTemperature) + " C, and the wire of " +
                windingLabel + " states no temperature limit (" + (coating ? "its coating has no temperature rating and no rated material"
                                                                            : "it has no coating") +
                "): its temperature is not checked"));
            continue;
        }
        std::vector<std::string> exceeded;
        for (const auto& limit : windingLimits) {
            if (hottestTemperature > limit.temperature) {
                exceeded.push_back(limit.source + " of " + format_number(limit.temperature) + " C");
            }
        }
        if (!exceeded.empty()) {
            std::string message = "At " + where + ": the design's hottest point is " + format_number(hottestTemperature) + " C, above " +
                                  windingLabel + "'s ";
            for (size_t index = 0; index < exceeded.size(); ++index) {
                message += (index > 0 ? " and " : "") + exceeded[index];
            }
            findings.push_back(make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::WINDING_TEMPERATURE_ABOVE_LIMIT, message));
        }
    }
    return findings;
}

std::vector<ResultFinding> ResultsAnalyzer::analyze_missing_temperature(const Magnetic& magnetic, const OperatingPoint& operatingPoint,
                                                                        size_t operatingPointIndex, std::optional<Topology> topology) {
    std::string where = describe_operating_point(operatingPoint, operatingPointIndex);
    // The same solve simulate ran for this point: wound as simulate winds it, the full winding-loss model.
    Magnetic woundMagnetic = magnetic;
    if (!woundMagnetic.get_coil().get_turns_description()) {
        woundMagnetic.get_mutable_coil().wind();
    }
    OperatingPoint solvedOperatingPoint = operatingPoint;
    MagneticSimulator simulator;
    simulator.set_topology(topology);
    std::optional<WindingLossesOutput> windingLossesOutput;
    std::optional<TemperatureOutput> temperatureOutput;
    try {
        simulator.calculate_core_losses(solvedOperatingPoint, woundMagnetic, &windingLossesOutput, &temperatureOutput,
                                        WindingLossesScope::FULL, NoThermalOperatingPointPolicy::THROW);
    }
    catch (const NoThermalOperatingPointException& noThermalOperatingPoint) {
        return {make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::NO_THERMAL_OPERATING_POINT,
            "At " + where + ": " + noThermalOperatingPoint.message() + " Its core and winding losses are those evaluated at the "
            "ambient temperature (" + format_number(operatingPoint.get_conditions().get_ambient_temperature()) +
            " C), and it has no temperature output: its temperatures are not checked")};
    }
    throw InvalidInputException(ErrorCode::MISSING_DATA,
        "analyze_results: the output of " + where + " has no temperature, yet its losses have a thermal operating point "
        "(hottest point " + format_number(temperatureOutput.value().get_maximum_temperature()) + " C): it was simulated in "
        "fast mode (or edited after simulate), which publishes no temperature; analyze_results needs simulate's full output");
}

std::vector<ResultFinding> ResultsAnalyzer::analyze_results(const Mas& mas) {
    const auto& operatingPoints = mas.get_inputs().get_operating_points();
    const auto& outputs = mas.get_outputs();
    if (outputs.size() != operatingPoints.size()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "analyze_results: the MAS has " + std::to_string(operatingPoints.size()) + " operating points and " +
            std::to_string(outputs.size()) + " outputs; it needs one simulated output per operating point");
    }
    std::vector<ResultFinding> findings;
    for (size_t index = 0; index < operatingPoints.size(); ++index) {
        for (auto& finding : analyze_saturation(mas.get_magnetic(), operatingPoints[index], outputs[index], index,
                                                     mas.get_inputs().get_design_requirements().get_topology())) {
            findings.push_back(finding);
        }
        for (auto& finding : analyze_volt_second_balance(operatingPoints[index], index)) {
            findings.push_back(finding);
        }
        if (!outputs[index].get_temperature()) {
            for (auto& finding : analyze_missing_temperature(mas.get_magnetic(), operatingPoints[index], index,
                                                             mas.get_inputs().get_design_requirements().get_topology())) {
                findings.push_back(finding);
            }
            continue;
        }
        for (auto& finding : analyze_temperature(mas.get_magnetic(), operatingPoints[index], outputs[index], index)) {
            findings.push_back(finding);
        }
    }
    return findings;
}

} // namespace OpenMagnetics
