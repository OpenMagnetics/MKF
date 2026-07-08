#include "physical_models/Impedance.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/InitialPermeability.h"
#include "physical_models/Reluctance.h"
#include "physical_models/StrayCapacitance.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/LeakageInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "physical_models/WindingLosses.h"
#include "constructive_models/NumberTurns.h"
#include "support/Utils.h"
#include "support/Settings.h"
#include <cmath>
#include <numbers>


namespace OpenMagnetics {

std::complex<double> Impedance::calculate_impedance(Magnetic magnetic, double frequency, double temperature) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_impedance(core, coil, frequency, temperature);
}

std::complex<double> Impedance::calculate_differential_mode_impedance(Magnetic magnetic, double frequency, double temperature) {
    return calculate_differential_mode_impedance(magnetic.get_core(), magnetic.get_coil(), frequency, temperature);
}

DifferentialModeParameters Impedance::calculate_differential_mode_parameters(Core core, Coil coil, double referenceFrequency, double temperature) {
    // Differential mode: the two windings carry opposing currents, so the flux
    // they drive into the core cancels. The inductance seen is therefore the
    // *leakage* inductance (the flux that does not couple through the core),
    // which is essentially air-cored and so frequency-flat. It parallel-
    // resonates with the *inter-winding* capacitance — the choke's second
    // self-resonance. All three terms are frequency-independent (the leakage is
    // taken at referenceFrequency), so a sweep computes them once.
    if (coil.get_functional_description().size() < 2) {
        throw std::runtime_error("Differential-mode impedance requires at least two windings");
    }

    Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);

    double leakageInductance = LeakageInductance()
        .calculate_leakage_inductance(magnetic, referenceFrequency, 0, 1)
        .get_leakage_inductance_per_winding()[0]
        .get_nominal()
        .value();

    // The leakage path is essentially air-cored (no core permeability), so the
    // inductive branch is purely reactive; the winding resistance is the
    // dominant loss term. Sign convention matches calculate_impedance().
    double windingResistance = WindingOhmicLosses::calculate_dc_resistance_per_winding(coil, temperature)[0];

    // Inter-winding capacitance: the off-diagonal term of the stray-capacitance
    // matrix (between the two windings). On a separated-winding CMC the two
    // windings have no adjacent turns, so this is the through-core path (passing
    // the core lets the model build it). The DM resonance is set by this, not by
    // the per-winding self-capacitance used in common mode.
    if (!coil.get_turns_description()) {
        coil.wind();
    }
    auto capacitanceMatrix = StrayCapacitance().calculate_capacitance(coil, core).get_capacitance_among_windings().value();
    auto primaryName = coil.get_functional_description()[0].get_name();
    auto secondaryName = coil.get_functional_description()[1].get_name();
    double interWindingCapacitance = capacitanceMatrix[primaryName][secondaryName];

    return {leakageInductance, windingResistance, interWindingCapacitance};
}

std::complex<double> Impedance::differential_mode_impedance_from_parameters(const DifferentialModeParameters& parameters, double frequency) {
    auto angularFrequency = 2 * std::numbers::pi * frequency;
    auto inductiveImpedance = std::complex<double>(parameters.windingResistance, angularFrequency * parameters.leakageInductance);
    if (parameters.interWindingCapacitance <= 0) {
        // No capacitive path: the DM impedance is purely the leakage branch.
        return inductiveImpedance;
    }
    auto capacitiveImpedance = std::complex<double>(0, -1.0 / (angularFrequency * parameters.interWindingCapacitance));
    return 1.0 / (1.0 / inductiveImpedance + 1.0 / capacitiveImpedance);
}

std::complex<double> Impedance::calculate_differential_mode_impedance(Core core, Coil coil, double frequency, double temperature) {
    auto parameters = calculate_differential_mode_parameters(core, coil, frequency, temperature);
    return differential_mode_impedance_from_parameters(parameters, frequency);
}

ImpedanceTank Impedance::build_magnetizing_tank(Core& core, Coil& coil) {
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory();
    double numberTurns = coil.get_functional_description()[0].get_number_turns();
    double reluctanceCoreUnityPermeability = reluctanceModel->get_core_reluctance(core, 1).get_core_reluctance();
    double airCoredInductance = numberTurns * numberTurns / reluctanceCoreUnityPermeability;

    double capacitance;
    auto& settings = Settings::GetInstance();
    if (_fastCapacitance) {
        capacitance = StrayCapacitanceOneLayer().calculate_capacitance(coil);
    }
    else {
        auto strayCapacitanceModel = settings.get_stray_capacitance_model();
        auto capacitanceMatrix = StrayCapacitance(strayCapacitanceModel).calculate_capacitance(coil).get_capacitance_among_windings().value();
        capacitance = capacitanceMatrix[coil.get_functional_description()[0].get_name()][coil.get_functional_description()[0].get_name()];
    }

    return ImpedanceTank{airCoredInductance, capacitance, true};
}

std::complex<double> Impedance::impedance_from_model(const WidebandImpedanceModel& model, double frequency) {
    auto angularFrequency = 2 * std::numbers::pi * frequency;

    // Core complex permeability for the magnetizing tank, evaluated once.
    // Standard e^{jωt} convention: Z_L = jωL_air·µ = jωL_air·(µ' − jµ'') = ωL_air·(µ'' + jµ').
    // So the real part (loss) is ωL_air·µ'' and the reactance is +ωL_air·µ' (inductive).
    double complexPermeabilityRealPart = 1.0;
    double complexPermeabilityImaginaryPart = 0.0;
    if (model.coreMaterial) {
        auto [muReal, muImag] = OpenMagnetics::ComplexPermeability().get_complex_permeability(model.coreMaterial.value(), frequency);
        complexPermeabilityRealPart = muReal * model.permeabilityScaling;
        complexPermeabilityImaginaryPart = muImag * model.permeabilityScaling;
    }

    // The leakage tanks are damped by their winding resistance, which is
    // frequency-dependent (skin/proximity). Evaluate one winding's resistance at
    // this frequency: fast => DC × analytic skin factor; slow => full field-based
    // effective resistance (adds proximity).
    auto windingResistanceAt = [&](size_t windingIndex) -> double {
        if (!model.fast && model.hasProximityModel && model.magnetic) {
            return WindingLosses::calculate_effective_resistance_of_winding(model.magnetic.value(), windingIndex, frequency, model.temperature);
        }
        if (model.hasWindingResistanceModel
                && windingIndex < model.windingDcResistancePerMeter.size()
                && model.windingDcResistancePerMeter[windingIndex] > 0) {
            double skinFactor = WindingOhmicLosses::calculate_effective_resistance_per_meter(model.windingWire[windingIndex], frequency, model.temperature) / model.windingDcResistancePerMeter[windingIndex];
            return model.windingResistanceDc[windingIndex] * skinFactor;
        }
        if (windingIndex < model.windingResistanceDc.size()) {
            return model.windingResistanceDc[windingIndex];  // no skin model for this wire -> DC
        }
        return 0.0;
    };

    // Winding 0's resistance is common to every leakage loop, so evaluate it once.
    bool hasLeakageTank = false;
    for (const auto& tank : model.tanks) {
        if (!tank.usesCorePermeability) {
            hasLeakageTank = true;
            break;
        }
    }
    double primaryResistance = hasLeakageTank ? windingResistanceAt(0) : 0.0;

    // Series cascade of parallel-RLC tanks (Foster ladder): each tank adds one
    // resonance. At low frequency the leakage tanks reduce to their series
    // R + jωL_leak arm (their shunt capacitance is negligible), so the dominant
    // magnetizing resonance is preserved; near each leakage resonance the shunt
    // inter-winding capacitance produces the corresponding higher self-resonance.
    std::complex<double> impedance = 0.0;
    for (const auto& tank : model.tanks) {
        std::complex<double> seriesArm;
        if (tank.usesCorePermeability) {
            seriesArm = angularFrequency * tank.inductance * std::complex<double>(complexPermeabilityImaginaryPart, complexPermeabilityRealPart);
        }
        else {
            // Leakage loop resistance referred to winding 0: R_0 + (N_0/N_j)²·R_j.
            double leakageResistance = primaryResistance + tank.turnsRatioSquared * windingResistanceAt(tank.secondaryWindingIndex);
            seriesArm = std::complex<double>(leakageResistance, angularFrequency * tank.inductance);
        }

        if (tank.capacitance > 0) {
            auto capacitiveArm = std::complex<double>(0, -1.0 / (angularFrequency * tank.capacitance));
            impedance += 1.0 / (1.0 / seriesArm + 1.0 / capacitiveArm);
        }
        else {
            impedance += seriesArm;
        }
    }

    return impedance;
}

WidebandImpedanceModel Impedance::build_wideband_impedance_model(Magnetic magnetic, double referenceFrequency, double temperature, bool fast) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();

    WidebandImpedanceModel model;
    model.coreMaterial = core.resolve_material();
    model.temperature = temperature;
    model.fast = fast;

    // Tank 0: the magnetizing resonance (core µ(f) ∥ winding self-capacitance),
    // present for every magnetic — the same first resonance as calculate_impedance.
    model.tanks.push_back(build_magnetizing_tank(core, coil));

    // For coupled magnetics each additional winding adds a leakage resonance: the
    // leakage inductance between winding 0 and winding j (the flux that does not
    // couple through the core, essentially air-cored and frequency-flat) resonates
    // with the inter-winding capacitance between them. This is the choke "second
    // resonance", generalized to one tank per extra winding. A single-winding
    // inductor adds no such tank and keeps its single resonance.
    size_t numberWindings = coil.get_functional_description().size();
    if (numberWindings >= 2) {
        if (!coil.get_turns_description()) {
            coil.wind();
            magnetic.set_coil(coil);
        }

        // Inter-winding capacitances: the off-diagonal terms of the stray-capacitance
        // matrix (the through-core path on a separated-winding choke). The whole
        // matrix is computed once here, not per frequency.
        auto capacitanceMatrix = StrayCapacitance().calculate_capacitance(coil, core).get_capacitance_among_windings().value();
        auto primaryName = coil.get_functional_description()[0].get_name();
        double primaryTurns = coil.get_functional_description()[0].get_number_turns();

        // Per-winding resistance data. The leakage path is air-cored (purely
        // reactive), so the winding resistance is its dominant loss term, matching
        // calculate_differential_mode_parameters(). Each leakage loop runs through
        // winding 0 and one secondary, so the evaluator damps it with both
        // windings' resistances (referred to winding 0).
        model.windingResistanceDc = WindingOhmicLosses::calculate_dc_resistance_per_winding(coil, temperature);
        if (fast) {
            // Precompute the DC per-meter resistance so the per-frequency skin factor
            // is just R_pm(f)/R_pm_dc — no coil field map needed.
            model.windingWire.reserve(numberWindings);
            model.windingDcResistancePerMeter.reserve(numberWindings);
            for (size_t windingIndex = 0; windingIndex < numberWindings; ++windingIndex) {
                auto wire = coil.resolve_wire(windingIndex);
                model.windingWire.push_back(wire);
                model.windingDcResistancePerMeter.push_back(WindingOhmicLosses::calculate_dc_resistance_per_meter(wire, temperature));
            }
            model.hasWindingResistanceModel = true;
        }
        else {
            // DC + skin + proximity: the per-point field-based effective resistance
            // needs the magnetic, so keep it for the evaluator.
            model.hasProximityModel = true;
            model.magnetic = magnetic;
        }

        for (size_t windingIndex = 1; windingIndex < numberWindings; ++windingIndex) {
            double leakageInductance = LeakageInductance()
                .calculate_leakage_inductance(magnetic, referenceFrequency, 0, windingIndex)
                .get_leakage_inductance_per_winding()[0]
                .get_nominal()
                .value();
            if (leakageInductance <= 0) {
                continue;  // no leakage path to this winding -> no extra resonance
            }
            auto secondaryName = coil.get_functional_description()[windingIndex].get_name();
            double interWindingCapacitance = capacitanceMatrix[primaryName][secondaryName];
            // Referral factor (N_0/N_j)² for the secondary resistance in this leakage loop.
            double secondaryTurns = coil.get_functional_description()[windingIndex].get_number_turns();
            double turnsRatioSquared = (secondaryTurns > 0) ? std::pow(primaryTurns / secondaryTurns, 2) : 1.0;
            model.tanks.push_back(ImpedanceTank{leakageInductance, interWindingCapacitance, false, windingIndex, turnsRatioSquared});
        }
    }

    return model;
}

WidebandImpedanceModel Impedance::build_common_mode_impedance_model(Magnetic magnetic, double temperature) {
    // All windings driven in parallel (the CM measurement): the leakage flux
    // between them is not excited, so the terminal impedance is the magnetizing
    // tank alone. The paralleled windings' series resistance (mΩ) is negligible
    // against the µ''-damped magnetizing branch, so no extra loss term is added.
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();

    WidebandImpedanceModel model;
    model.coreMaterial = core.resolve_material();
    model.temperature = temperature;
    model.tanks.push_back(build_magnetizing_tank(core, coil));
    return model;
}

std::complex<double> Impedance::calculate_impedance(Core core, Coil coil, double frequency, double temperature) {
    // Single-point common-mode/terminal impedance: just the magnetizing tank (the
    // first resonance). The wideband sweep adds the leakage tanks on top of this.
    WidebandImpedanceModel model;
    model.coreMaterial = core.resolve_material();
    model.temperature = temperature;
    model.tanks.push_back(build_magnetizing_tank(core, coil));
    return impedance_from_model(model, frequency);
}

std::complex<double> Impedance::calculate_impedance(Core core, Coil coil, double frequency, double magneticFieldDcBias, double temperature) {
    // DC-biased impedance: scale the small-signal complex permeability by
    // the ratio µ(H_dc) / µ(0) to account for permeability rolloff under
    // DC bias. This first-order correction is physically valid because the
    // DC bias shifts the operating point on the BH curve, reducing the
    // effective incremental permeability for both real and imaginary parts.
    auto coreMaterial = core.resolve_material();
    double muInitial = InitialPermeability::get_initial_permeability(coreMaterial, temperature, std::nullopt, frequency);
    double muBiased = InitialPermeability::get_initial_permeability(coreMaterial, temperature, magneticFieldDcBias, frequency);
    double biasRatio = (muInitial > 0) ? (muBiased / muInitial) : 1.0;

    // Same magnetizing tank as the unbiased path, with the complex permeability
    // scaled by µ(H_dc)/µ(0) to account for permeability rolloff under DC bias.
    WidebandImpedanceModel model;
    model.coreMaterial = coreMaterial;
    model.permeabilityScaling = biasRatio;
    model.temperature = temperature;
    model.tanks.push_back(build_magnetizing_tank(core, coil));
    return impedance_from_model(model, frequency);
}

int64_t Impedance::calculate_minimum_number_turns(Magnetic magnetic, Inputs inputs) {
    NumberTurns numberTurns(1, inputs.get_design_requirements());
    if (!inputs.get_design_requirements().get_minimum_impedance()) {
        throw std::invalid_argument("Missing impedance requirement");
    }
    auto minimumImpedanceRequirement = inputs.get_design_requirements().get_minimum_impedance().value();
    auto temperature = inputs.get_maximum_temperature();
    auto timeout = _maximumNumberTurns;
    int64_t calculatedNumberTurns = -1;
    while (true) {
        auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();
        calculatedNumberTurns = numberTurnsCombination[0];
        magnetic.get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(calculatedNumberTurns);
        auto selfResonantFrequency = calculate_self_resonant_frequency(magnetic);

        bool validDesign = true;
        for (auto impedanceAtFrequency : minimumImpedanceRequirement) {
            auto frequency = impedanceAtFrequency.get_frequency();
            auto requiredImpedanceMagnitude = impedanceAtFrequency.get_impedance().get_magnitude();
            auto impedanceMagnitude = abs(calculate_impedance(magnetic, frequency, temperature));
            if (impedanceMagnitude < requiredImpedanceMagnitude) {
                validDesign = false;
            }
            if (frequency > 0.25 * selfResonantFrequency) {  // hardcoded 20% of SRF
                validDesign = false;
                return -1;
            }
        }

        if (validDesign) {
            break;
        }

        if (timeout == 0) {
            return -1;
        }
        timeout--;
    }

    return calculatedNumberTurns;
}

double Impedance::calculate_q_factor(Magnetic magnetic, double frequency, double temperature) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_q_factor(core, coil, frequency, temperature);
}

double Impedance::calculate_q_factor(Core core, Coil coil, double frequency, double temperature) {
    auto impedance = calculate_impedance(core, coil, frequency, temperature);
    return impedance.imag() / impedance.real();
}

double Impedance::calculate_self_resonant_frequency(Magnetic magnetic, double temperature) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    return calculate_self_resonant_frequency(core, coil, temperature);
}

double Impedance::calculate_self_resonant_frequency(Core core, Coil coil, double temperature) {
    double capacitance;
    if (_fastCapacitance) {
        capacitance = StrayCapacitanceOneLayer().calculate_capacitance(coil);
    }
    else {
        if (!coil.get_turns_description()) {
            coil.wind();
        }
        auto capacitanceMatrix = StrayCapacitance().calculate_capacitance(coil).get_capacitance_among_windings().value();

        capacitance = capacitanceMatrix[coil.get_functional_description()[0].get_name()][coil.get_functional_description()[0].get_name()];
    }

    OperatingPoint operatingPoint;
    OperatingConditions conditions;
    conditions.set_cooling(std::nullopt);
    conditions.set_ambient_temperature(temperature);
    operatingPoint.set_conditions(conditions);
    MagnetizingInductance magnetizingInductanceModel("ZHANG");
    double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(core, coil, nullptr).get_magnetizing_inductance().get_nominal().value();

    double selfResonantFrequency = 1.0 / (2 * std::numbers::pi * sqrt(magnetizingInductance * capacitance));

    return selfResonantFrequency;
}

} // namespace OpenMagnetics
