#include "processors/ResultsAnalyzer.h"
#include "constructive_models/Core.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/Wire.h"
#include "physical_models/InitialPermeability.h"
#include "physical_models/MagnetizingInductance.h"
#include "processors/MagneticSimulator.h"
#include "support/Exceptions.h"
#include "support/Utils.h"
#include "support/NumberFormat.h"
#include "physical_models/WindingOhmicLosses.h"

#include <magic_enum.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
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
        case ResultFindingCode::LEAD_COLLIDES_WITH_TURNS: return "LEAD_COLLIDES_WITH_TURNS";
        case ResultFindingCode::LEAD_ROUTE_NOT_MODELLED: return "LEAD_ROUTE_NOT_MODELLED";
        case ResultFindingCode::EDGEWISE_BEND_BELOW_COATING_MANDREL: return "EDGEWISE_BEND_BELOW_COATING_MANDREL";
    }
    throw InvalidInputException(ErrorCode::INVALID_INPUT, "Unknown result finding code " + std::to_string(static_cast<int>(code)));
}

void to_json(nlohmann::json& j, const ResultFinding& finding) {
    j = nlohmann::json::object();
    if (finding.operatingPointIndex) {
        j["operatingPointIndex"] = finding.operatingPointIndex.value();
    }
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
    // Every finding reads as whole sentences (ABT #2001).
    finding.message = (!message.empty() && message.back() != '.') ? message + "." : message;
    return finding;
}

// A number in a message: three significant figures, all integer digits kept (2474.39 C -> "2474") (ABT #2001).
std::string format_number(double value) {
    return format_significant(value, 3);
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
    // Criterion 3 (ABT #1942, Alf 2026-10-18): where the magnetizing model read B on the magnetisation curve
    // sample by sample, the B peak is the curve's, which past the knee rises only on the vacuum slope and can stay
    // under B_sat while the core is driven deep into saturation. The field strength says it: its peak against the
    // material's H_sat at the check's temperature (the one B_sat is read at), times the same margin as B.
    verdict.fieldStrengthAboveSaturation = false;
    if (dcBiasPoint && dcBiasPoint->peakMagneticFieldStrength) {
        verdict.peakFieldStrength = dcBiasPoint->peakMagneticFieldStrength.value();
        verdict.saturationFieldStrength = Core::get_magnetic_field_strength_saturation(coreMaterial, temperature);
        verdict.fieldStrengthAboveSaturation = verdict.peakFieldStrength.value() * margin > verdict.saturationFieldStrength.value();
    }
    // Criterion 4 (Alf, 2026-10-24): the DC solve closed on a jump of the biased permeability, not on a root, and
    // took the saturated mu_r = 1 (MagnetizingInductance). The bias at mu_r = 1 can sit inside the material's
    // DC-bias table, where criterion 2 does not see it, so the solve's own mark says it.
    verdict.noDcOperatingPoint = dcBiasPoint && dcBiasPoint->noDcOperatingPoint;
    verdict.saturates = verdict.fluxDensityAboveSaturation || verdict.permeabilityCollapsed || verdict.fieldStrengthAboveSaturation ||
                        verdict.noDcOperatingPoint;

    // ABT #2001: the first sentence says the problem with rounded numbers; the explanation follows.
    std::ostringstream message;
    std::string marginNote = margin != 1 ? " (" + format_number(peakFluxDensity * margin) + " T with the " + format_number(margin) + " margin)" : "";
    if (verdict.saturates) {
        // One clause per criterion that fired, joined into the first sentence.
        std::vector<std::string> clauses;
        if (verdict.fluxDensityAboveSaturation) {
            std::ostringstream clause;
            clause << "B peak = " << format_number(peakFluxDensity) << " T" << marginNote << " is above B_sat = "
                   << format_number(verdict.saturationFluxDensity) << " T of " << coreMaterial.get_name() << " at "
                   << format_number(temperature) << " C";
            clauses.push_back(clause.str());
        }
        if (verdict.permeabilityCollapsed) {
            std::ostringstream clause;
            clause << "its DC bias of " << format_quantity(dcBiasPoint->magneticFieldDcBias, "A/m")
                   << " drives the permeability model of " << coreMaterial.get_name() << " to mu_r = " << format_number(verdict.biasedPermeability.value())
                   << " (unbiased " << format_number(unbiasedPermeability.value()) << ") at " << format_number(dcBiasPoint->temperature)
                   << " C and " << format_quantity(dcBiasPoint->frequency, "Hz") << ", past its knee, where the model takes mu_r = 1";
            clauses.push_back(clause.str());
        }
        if (verdict.fieldStrengthAboveSaturation) {
            std::ostringstream clause;
            clause << "its peak field strength on the magnetisation curve of " << coreMaterial.get_name() << ", H = "
                   << format_quantity(verdict.peakFieldStrength.value(), "A/m");
            if (margin != 1) {
                clause << " (" << format_quantity(verdict.peakFieldStrength.value() * margin, "A/m") << " with the " << format_number(margin) << " margin)";
            }
            clause << ", is above its H_sat = " << format_quantity(verdict.saturationFieldStrength.value(), "A/m") << " at "
                   << format_number(temperature) << " C, past its knee";
            clauses.push_back(clause.str());
        }
        if (verdict.noDcOperatingPoint) {
            std::ostringstream clause;
            clause << "under its DC bias at " << format_number(dcBiasPoint->temperature) << " C and " << format_quantity(dcBiasPoint->frequency, "Hz")
                   << " " << coreMaterial.get_name() << " has no operating point: the bias drives its permeability across the end of its"
                   << " DC-bias data, where it is taken as saturated (mu_r = 1), and no permeability reproduces the bias it implies";
            clauses.push_back(clause.str());
        }
        message << "The core saturates: ";
        for (size_t index = 0; index < clauses.size(); ++index) {
            if (index > 0) {
                message << (index + 1 == clauses.size() ? ", and " : ", ");
            }
            message << clauses[index];
        }
        message << ".";
        if (!verdict.fluxDensityAboveSaturation) {
            message << " B peak = " << format_number(peakFluxDensity) << " T is still below B_sat = " << format_number(verdict.saturationFluxDensity)
                    << " T at " << format_number(temperature) << " C";
            if (verdict.fieldStrengthAboveSaturation) {
                message << ": past the knee the curve's B rises only on the vacuum slope.";
            }
            else if (verdict.permeabilityCollapsed) {
                message << ": the material's B_sat, magnetisation curve and permeability model disagree.";
            }
            else {
                message << ": the DC bias, not the flux density, saturates it.";
            }
        }
        message << " The numbers for this operating point are those of the saturated core; it needs more turns, a larger gap or core,"
                << " or fewer volt-seconds";
    }
    else {
        message << "B peak = " << format_number(peakFluxDensity) << " T" << marginNote << " is within B_sat = "
                << format_number(verdict.saturationFluxDensity) << " T of " << coreMaterial.get_name() << " at " << format_number(temperature) << " C";
        if (verdict.peakFieldStrength) {
            message << ", and its peak field strength on the magnetisation curve, H = " << format_quantity(verdict.peakFieldStrength.value(), "A/m")
                    << ", within H_sat = " << format_quantity(verdict.saturationFieldStrength.value(), "A/m");
        }
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

std::vector<ResultFinding> ResultsAnalyzer::analyze_volt_second_balance(const OperatingPoint& operatingPoint, size_t operatingPointIndex,
                                                                        std::optional<std::vector<double>> windingDcResistances,
                                                                        std::optional<double> resistanceTemperature) {
    if (windingDcResistances.has_value() != resistanceTemperature.has_value()) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT,
            "analyze_volt_second_balance: the windings' DC resistances and the temperature they were computed at go together");
    }
    std::vector<ResultFinding> findings;
    const auto& excitations = operatingPoint.get_excitations_per_winding();
    for (size_t windingIndex = 0; windingIndex < excitations.size(); ++windingIndex) {
        auto excitation = excitations[windingIndex];
        if (!excitation.get_voltage() || !excitation.get_voltage()->get_waveform()) {
            continue;
        }
        // ABT #2002: the winding's own resistive drop i(t) R_dc is not flux drift; only what remains
        // once it is taken out is judged. Without a resistance (no wound coil) or a current, the voltage
        // is judged as given.
        std::optional<double> dcResistance;
        std::optional<double> meanResistiveDrop;
        if (windingDcResistances && windingIndex < windingDcResistances->size() && excitation.get_current()) {
            dcResistance = windingDcResistances->at(windingIndex);
            auto currentWaveform = Inputs::standardize_waveform(excitation.get_current().value(), excitation.get_frequency()).get_waveform().value();
            meanResistiveDrop = dcResistance.value() * Inputs::calculate_waveform_average(currentWaveform, excitation.get_frequency());
            excitation = Inputs::remove_resistive_drop(excitation, dcResistance.value());
        }
        double frequency = excitation.get_frequency();
        auto voltage = Inputs::standardize_waveform(excitation.get_voltage().value(), frequency);
        auto sampledVoltage = Inputs::calculate_sampled_waveform(voltage.get_waveform().value(), frequency);
        std::string context = "At " + describe_operating_point(operatingPoint, operatingPointIndex) + ", winding " +
                              std::to_string(windingIndex) + " ('" + excitation.get_name().value_or("unnamed") + "')";
        auto imbalance = Inputs::describe_volt_second_imbalance(sampledVoltage, context);
        if (imbalance) {
            std::string message = imbalance.value() + ". Where MKF derives a magnetizing current from this voltage, it integrates "
                                  "only its balanced part (the mean removed), so the numbers given describe that part";
            if (dcResistance) {
                message += ". The winding's own resistive drop (" + format_quantity(dcResistance.value(), "Ohm") + " at " +
                           format_number(resistanceTemperature.value()) + " C, " + format_quantity(meanResistiveDrop.value(), "V") +
                           " on average) was taken out of the voltage before this check";
            }
            findings.push_back(make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::VOLT_SECOND_IMBALANCE, message));
        }
    }
    return findings;
}

std::vector<ResultFinding> ResultsAnalyzer::analyze_volt_second_balance(const Magnetic& magnetic, const OperatingPoint& operatingPoint,
                                                                        const Outputs& output, size_t operatingPointIndex) {
    // The resistance at the temperature the output's winding losses were evaluated at (the solved one),
    // or, for an output without them, at the operating point's own temperature.
    double temperature = operatingPoint.get_conditions().get_ambient_temperature();
    if (output.get_winding_losses() && output.get_winding_losses()->get_temperature()) {
        temperature = output.get_winding_losses()->get_temperature().value();
    }
    auto windingDcResistances = Inputs::get_winding_dc_resistances(magnetic.get_coil(), temperature);
    return analyze_volt_second_balance(operatingPoint, operatingPointIndex, windingDcResistances,
                                       windingDcResistances ? std::optional<double>(temperature) : std::nullopt);
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
    // The ceiling covers the mu_i(T) / Curie data and the evaluated Steinmetz ranges' ct(T) (ABT #1908).
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
            case CoreDataCeilingSource::STEINMETZ_TEMPERATURE_FACTOR_END:
                limit = "the temperature where the Steinmetz temperature factor ct(T) of " + coreMaterial.get_name() +
                        "'s loss data at this frequency reaches zero";
                magnetized = "with its own permeability at that temperature";
                break;
        }
        findings.push_back(make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::CORE_ABOVE_MATERIAL_DATA_LIMIT,
            "At " + where + ": the losses drive the core above " + limit + " (" + format_number(ceiling) +
            " C), where its material data ends. MKF evaluated its core losses at " + format_number(coreTemperature) +
            " C, the last temperature the data covers, with the core magnetized " + magnetized + "; the temperatures given "
            "for this operating point (hottest point " + format_number(hottestTemperature) + " C) follow from those losses, "
            "and the core's real temperature is higher"));
    }
    std::vector<StatedTemperatureLimit> coreLimits;
    if (coreMaterial.get_curie_temperature()) {
        coreLimits.push_back({coreMaterial.get_curie_temperature().value(), "Curie temperature"});
    }
    if (coreMaterial.get_recommendations() && coreMaterial.get_recommendations()->get_maximum_operating_temperature()) {
        coreLimits.push_back({coreMaterial.get_recommendations()->get_maximum_operating_temperature().value(),
                              "maximum recommended operating temperature"});
    }
    if (coreLimits.empty()) {
        findings.push_back(make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::CORE_TEMPERATURE_LIMIT_NOT_STATED,
            "At " + where + ": the core runs at " + format_number(coreTemperature) + " C, and its material " + coreMaterial.get_name() +
            " states no temperature limit (no Curie temperature, no maximum operating temperature), so its temperature is not checked"));
    }
    else {
        std::vector<std::string> exceeded;
        for (const auto& limit : coreLimits) {
            if (coreTemperature > limit.temperature) {
                exceeded.push_back("the " + format_number(limit.temperature) + " C " + limit.source);
            }
        }
        if (!exceeded.empty()) {
            std::string message = "At " + where + ": the core runs at " + format_number(coreTemperature) + " C, above ";
            for (size_t index = 0; index < exceeded.size(); ++index) {
                message += (index > 0 ? " and " : "") + exceeded[index];
            }
            message += " of its material " + coreMaterial.get_name();
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
                windingLimits.push_back({coating->get_temperature_rating().value(), "temperature rating of the wire coating"});
            }
            // Only a material the coating names: an unnamed one is not assumed (Wire::resolve_coating_insulation_material
            // would take the default enamel).
            if (coating->get_material()) {
                auto insulationMaterial = Wire::resolve_coating_insulation_material(coating.value());
                if (insulationMaterial.get_temperature_class()) {
                    windingLimits.push_back({InsulationCoordinator::temperature_class_to_celsius(insulationMaterial.get_temperature_class().value()),
                                             "thermal class of the coating material (" + insulationMaterial.get_name() + ")"});
                }
            }
        }
        if (windingLimits.empty()) {
            findings.push_back(make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::WINDING_TEMPERATURE_LIMIT_NOT_STATED,
                "At " + where + ": the hottest point is " + format_number(hottestTemperature) + " C, and the wire of " +
                windingLabel + " states no temperature limit (" + (coating ? "its coating has no temperature rating and no rated material"
                                                                            : "it has no coating") +
                "), so its temperature is not checked"));
            continue;
        }
        std::vector<std::string> exceeded;
        for (const auto& limit : windingLimits) {
            if (hottestTemperature > limit.temperature) {
                exceeded.push_back("the " + format_number(limit.temperature) + " C " + limit.source);
            }
        }
        if (!exceeded.empty()) {
            std::string message = "At " + where + ": the hottest point, " + format_number(hottestTemperature) + " C, is above ";
            for (size_t index = 0; index < exceeded.size(); ++index) {
                message += (index > 0 ? " and " : "") + exceeded[index];
            }
            message += " of " + windingLabel;
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
            " C), and it has no temperature output, so its temperatures are not checked")};
    }
    throw InvalidInputException(ErrorCode::MISSING_DATA,
        "analyze_results: the output of " + where + " has no temperature, yet its losses have a thermal operating point "
        "(hottest point " + format_number(temperatureOutput.value().get_maximum_temperature()) + " C): it was simulated in "
        "fast mode (or edited after simulate), which publishes no temperature; analyze_results needs simulate's full output");
}

std::vector<ResultFinding> ResultsAnalyzer::analyze_flux_density_past_magnetisation_curve(const Magnetic& magnetic,
                                                                                          const OperatingPoint& operatingPoint,
                                                                                          size_t operatingPointIndex,
                                                                                          std::optional<Topology> topology) {
    std::string where = describe_operating_point(operatingPoint, operatingPointIndex);
    // The magnetizing model as simulate ran it on this point, at its ambient temperature.
    OperatingPoint operatingPointAtAmbient = operatingPoint;
    MagnetizingInductance magnetizingInductanceModel;
    magnetizingInductanceModel.set_topology(topology);
    std::string reason;
    try {
        magnetizingInductanceModel.calculate_inductance_and_magnetic_flux_density(magnetic, &operatingPointAtAmbient);
    }
    catch (const SaturationException& error) {
        reason = error.message();
    }
    catch (const MagnetisationCurveExhaustedException& error) {
        reason = error.message();
    }
    if (reason.empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA,
            "analyze_results: the output of " + where + " has no core losses, yet its flux density reads on the magnetisation "
            "curve at its ambient temperature: it was edited after simulate, which gives every other point its core losses");
    }
    return {make_finding(operatingPoint, operatingPointIndex, ResultFindingCode::CORE_SATURATES,
        "At " + where + ": the core saturates: " + reason + ". Its magnetizing inductance, flux density, core losses and "
        "temperatures are not computed for this operating point (its winding losses are those of its currents at the ambient "
        "temperature, " + format_number(operatingPoint.get_conditions().get_ambient_temperature()) + " C); more turns, a larger "
        "gap or core, or less volt-seconds are needed")};
}

std::vector<ResultFinding> ResultsAnalyzer::analyze_lead_collisions(const Magnetic& magnetic) {
    // A copy: judging resolves the coil's bobbin and column frames, and the magnetic is const here.
    OpenMagnetics::Coil coil = magnetic.get_coil();
    // Three significant figures in mm, as every finding's numbers (ABT #2001).
    auto millimetres = [](double metres) {
        return format_significant(metres * 1000, 3);
    };
    std::vector<ResultFinding> findings;
    for (const auto& lead : coil.judge_lead_collisions()) {
        const auto worst = std::max_element(lead.collisions.begin(), lead.collisions.end(),
                                            [](const LeadTurnCollision& a, const LeadTurnCollision& b) { return a.overlap < b.overlap; });
        std::string message = "Warning: winding '" + lead.windingName + "' has its " + lead.leadName +
                              " drawn across turn '" + worst->turnName + "': " + millimetres(worst->distanceOut) +
                              " mm out along the lead, that turn runs at height " + millimetres(worst->turnHeight) +
                              " mm and the lead at " + millimetres(worst->leadHeight) + " mm, so their copper overlaps by " +
                              millimetres(worst->overlap) + " mm";
        if (lead.collisions.size() > 1) {
            message += " (the worst of " + std::to_string(lead.collisions.size()) + " turns it crosses)";
        }
        message += ". The coil is drawn as routed; the lead collides with the winding there.";
        ResultFinding finding;
        finding.code = ResultFindingCode::LEAD_COLLIDES_WITH_TURNS;
        finding.message = message;
        findings.push_back(finding);
    }
    return findings;
}

std::vector<ResultFinding> ResultsAnalyzer::analyze_unmodelled_leads(const Magnetic& magnetic) {
    // A copy: judging resolves the coil's bobbin and column frames, and the magnetic is const here.
    OpenMagnetics::Coil coil = magnetic.get_coil();
    auto millimetres = [](double metres) {
        return format_significant(metres * 1000, 3);
    };
    std::vector<ResultFinding> findings;
    for (const auto& lead : coil.judge_unmodelled_tangent_leads()) {
        std::string message = "Warning: winding '" + lead.windingName + "' has its " + lead.leadName + " planned for a " +
                              millimetres(lead.plannedBendRadius) + " mm bend, which does not fit on the " +
                              millimetres(lead.connectionFaceHalfStraight) +
                              " mm connection-face half straight, so it leaves tangentially off the lateral face, on a route "
                              "MKF does not model in full: ";
        std::vector<std::string> conditions;
        for (const auto& condition : lead.conditions) {
            switch (condition.condition) {
                case TangentLeadCondition::DRAGBACK_RIDE:
                    conditions.push_back("a dragback lane rides " + millimetres(condition.rideHeight) +
                                         " mm over its connection face, which displaces that face and lengthens the lateral "
                                         "face by as much; the lead is drawn leaving where the unridden lateral face straight ends");
                    break;
                case TangentLeadCondition::LATERAL_COLUMN:
                    conditions.push_back("it is wound on a lateral column, " + millimetres(condition.columnAxisX) +
                                         " mm from the main one; the lead is drawn in that column's own frame with the main "
                                         "column's winding sense, which is not modelled for a lateral column");
                    break;
                case TangentLeadCondition::PIN_RUN:
                    conditions.push_back("it is run to " + (condition.pinName ? "pin '" + condition.pinName.value() + "'" : std::string("its pin")) +
                                         ", from the tangent lead's end on the window border as a level run; the lead's climb "
                                         "is not carried into it");
                    break;
                case TangentLeadCondition::LAYER_CHANGING_REVOLUTION:
                    conditions.push_back("the revolution it shortens, '" + condition.revolutionFrom + "' -> '" +
                                         condition.revolutionTo + "', changes layer and carries no pitch; the lead is drawn level "
                                         "at its station's height");
                    break;
                case TangentLeadCondition::SINGLE_REVOLUTION:
                    conditions.push_back("its conductor has a single turn and no pitch; the lead is drawn level at that turn's "
                                         "height, and both leads shorten that one turn");
                    break;
            }
        }
        for (size_t k = 0; k < conditions.size(); ++k) {
            message += (k == 0 ? "" : "; and ") + conditions[k];
        }
        message += ". The coil is drawn as routed.";
        ResultFinding finding;
        finding.code = ResultFindingCode::LEAD_ROUTE_NOT_MODELLED;
        finding.message = message;
        findings.push_back(finding);
    }
    return findings;
}

std::vector<ResultFinding> ResultsAnalyzer::analyze_winding_bends(const Magnetic& magnetic) {
    // A copy: judging resolves the coil's bobbin and column frames, and the magnetic is const here.
    OpenMagnetics::Coil coil = magnetic.get_coil();
    // A round or flatwise bend below its wire's limit is a layout real winding refuses (ABT #1290).
    // Reaching here with one means the coil was not wound by MKF under these settings; it is refused
    // the same way, never turned into a finding.
    coil.refuse_turns_bent_tighter_than_the_wire_allows();
    const auto judgement = coil.judge_corners_against_wire_flexibility();
    const auto leadJudgement = coil.judge_leads_against_wire_flexibility();

    // Three significant figures in mm, as every finding's numbers (ABT #2001).
    auto millimetres = [](double metres) {
        return format_significant(metres * 1000, 3);
    };
    auto describe = [&](const OpenMagnetics::Wire& wire) {
        std::string description = "rectangular";
        if (wire.get_name()) {
            description += " '" + wire.get_name().value() + "'";
        }
        if (wire.get_standard()) {
            description += ", " + std::string(magic_enum::enum_name(wire.get_standard().value()));
        }
        const double conductingWidth = resolve_dimensional_values(wire.get_conducting_width().value());
        const double conductingHeight = resolve_dimensional_values(wire.get_conducting_height().value());
        return description + ", conducting " + millimetres(conductingWidth) + " x " + millimetres(conductingHeight) + " mm";
    };
    // The limit a bend is held to, and how the mandrel's own centreline radius follows from it:
    // IEC 60851-3 5.1.1 winds the wire on the mandrel, so its centreline sits half the wire's outer
    // width beyond the mandrel's surface.
    auto limit = [&](const WoundCorner& asLaid) {
        return millimetres(asLaid.flexibilityRadius) + " mm its coating was qualified at on the " +
               millimetres(asLaid.mandrelDiameterUsed) + " mm mandrel of IEC 60317-0-2 Table 6";
    };
    auto derivation = [&](const WoundCorner& asLaid, const OpenMagnetics::Wire& wire) {
        const double outerWidth = OpenMagnetics::Wire::calculate_outer_width(wire);
        return "IEC 60851-3 5.1.1 winds the wire on that mandrel: centreline radius mandrel/2 + outer width/2 = " +
               millimetres(0.5 * asLaid.mandrelDiameterUsed) + " + " + millimetres(0.5 * outerWidth) + " = " +
               millimetres(asLaid.flexibilityRadius) + " mm";
    };
    const std::string trailer =
        " The coil is drawn as laid: edgewise winding is a forming process of its own and the mandrel is a "
        "coating test, not a forming limit, but the coating's integrity at this bend is outside the IEC qualification.";

    // One finding per winding, in the order the windings are met: its turns first, then its leads.
    std::vector<std::string> windingOrder;
    for (const auto& winding : judgement.windingsBelowFlexibility) {
        if (winding.axis != BendAxis::EDGEWISE) {
            // Unreachable after the refusal above, which throws for every non-edgewise winding here.
            throw InvalidInputException(ErrorCode::INVALID_INPUT,
                "analyze_winding_bends: winding '" + winding.windingName + "' is bent below its wire's flexibility radius "
                "on a non-edgewise axis, which real winding refuses");
        }
        windingOrder.push_back(winding.windingName);
    }
    for (const auto& winding : leadJudgement) {
        if (std::find(windingOrder.begin(), windingOrder.end(), winding.windingName) == windingOrder.end()) {
            windingOrder.push_back(winding.windingName);
        }
    }

    std::vector<ResultFinding> findings;
    for (const auto& windingName : windingOrder) {
        const auto turns = std::find_if(judgement.windingsBelowFlexibility.begin(), judgement.windingsBelowFlexibility.end(),
                                        [&](const WindingBendBelowFlexibility& winding) { return winding.windingName == windingName; });
        const auto leads = std::find_if(leadJudgement.begin(), leadJudgement.end(),
                                        [&](const WindingLeadsBelowFlexibility& winding) { return winding.windingName == windingName; });
        const bool hasTurns = turns != judgement.windingsBelowFlexibility.end();
        const bool hasLeads = leads != leadJudgement.end() && !leads->leads.empty();
        if (!hasTurns && !hasLeads) {
            throw std::logic_error("analyze_winding_bends: winding '" + windingName + "' was listed with no turn and no lead below its mandrel");
        }
        std::string message = "Warning: winding '" + windingName + "'";
        if (hasTurns) {
            message += " is bent edgewise to a " + millimetres(turns->bendRadius) + " mm centreline radius, tighter than the " +
                       limit(turns->asLaid) + ".";
        }
        if (hasLeads) {
            const auto describedWireName = hasTurns ? turns->wire.get_name() : leads->leads.front().wire.get_name();
            for (const auto& lead : leads->leads) {
                message += (hasTurns || &lead != &leads->leads.front()) ? " Its " : " has its ";
                message += lead.leadName;
                if (lead.wire.get_name() != describedWireName) {
                    message += " (" + describe(lead.wire) + ")";   // a lead conductor of its own
                }
                message += " planned for edgewise bends of " + millimetres(lead.plannedBendRadius) +
                           " mm centreline radius, tighter than the " + limit(lead.asLaid) + ".";
            }
        }
        const OpenMagnetics::Wire& wire = hasTurns ? turns->wire : leads->leads.front().wire;
        const WoundCorner& asLaid = hasTurns ? turns->asLaid : leads->leads.front().asLaid;
        message += " The wire is " + describe(wire);
        if (hasTurns) {
            message += "; the tightest bend is turn '" + turns->turnName + "', of " + std::to_string(turns->numberTightCorners) +
                       " such turns";
        }
        message += " (" + derivation(asLaid, wire) + ")." + trailer;

        ResultFinding finding;
        finding.code = ResultFindingCode::EDGEWISE_BEND_BELOW_COATING_MANDREL;
        finding.message = message;
        findings.push_back(finding);
    }
    return findings;
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
        if (!outputs[index].get_core_losses()) {
            // ABT #2052: simulate leaves the core losses (and the temperature) of a point whose flux density lies past
            // the end of the magnetisation curve at ambient unset; there is no flux density or temperature to check.
            for (auto& finding : analyze_flux_density_past_magnetisation_curve(mas.get_magnetic(), operatingPoints[index], index,
                                                                               mas.get_inputs().get_design_requirements().get_topology())) {
                findings.push_back(finding);
            }
            for (auto& finding : analyze_volt_second_balance(operatingPoints[index], index)) {
                findings.push_back(finding);
            }
            continue;
        }
        for (auto& finding : analyze_saturation(mas.get_magnetic(), operatingPoints[index], outputs[index], index,
                                                     mas.get_inputs().get_design_requirements().get_topology())) {
            findings.push_back(finding);
        }
        for (auto& finding : analyze_volt_second_balance(mas.get_magnetic(), operatingPoints[index], outputs[index], index)) {
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
    for (auto& finding : analyze_winding_bends(mas.get_magnetic())) {
        findings.push_back(finding);
    }
    for (auto& finding : analyze_lead_collisions(mas.get_magnetic())) {
        findings.push_back(finding);
    }
    for (auto& finding : analyze_unmodelled_leads(mas.get_magnetic())) {
        findings.push_back(finding);
    }
    return findings;
}

} // namespace OpenMagnetics
