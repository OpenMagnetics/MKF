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
#include <algorithm>
#include <complex>
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
    // ABT #383: a core carrying only its functionalDescription is a normal thing to hand
    // around — it is how MAS files are written when the constructive description is the
    // source of truth and the processed values are meant to be derived — but the reluctance
    // below reads the PROCESSED effective area and length. Without them it returned garbage
    // that looked like an answer: NaN for drumRing, exactly 0 for drumSemishielded, ~1e-10
    // ohm for molded. The NaN then surfaced far away, as "Waveform data contains NaN" out of
    // the waveform processor, which points at the wrong component entirely — the reporter
    // filed a model bug against impedance because of it. Refuse the input here, at the one
    // choke point all four sweep entry points share.
    if (!core.get_processed_description()) {
        throw CoreNotProcessedException(
            "Impedance needs the core's processed description (effective area and length) to "
            "build the magnetizing tank; the core carries only its functional description");
    }
    auto reluctanceModel = OpenMagnetics::ReluctanceModel::factory();
    double numberTurns = coil.get_functional_description()[0].get_number_turns();
    double reluctanceCoreUnityPermeability = reluctanceModel->get_core_reluctance(core, 1).get_core_reluctance();
    double airCoredInductance = numberTurns * numberTurns / reluctanceCoreUnityPermeability;

    double capacitance;
    auto& settings = Settings::GetInstance();
    if (_fastCapacitance) {
        capacitance = StrayCapacitanceOneLayer().calculate_capacitance(coil, core);
    }
    else {
        auto strayCapacitanceModel = settings.get_stray_capacitance_model();
        // The full model needs wound turns (OneLayer above does not), and the sweep
        // entry points hand over unwound coils; wind here or the branch throws
        // "Missing turns description" on 5 of the 107 WE catalogue chokes.
        if (!coil.get_turns_description()) {
            // A FAILED wind poisons the Coil: after one unsuccessful sectioned wind even a
            // fresh wind on the same object returns no turns, while the identical fresh
            // wind on an untouched copy succeeds (ABT #850). Keep a pristine copy so the
            // retry below starts clean.
            Coil pristineCoil = coil;
            coil.wind();
            if (!coil.get_turns_description()) {
                coil = pristineCoil;
            }
        }
        if (!coil.get_turns_description()) {
            // The winder can fail SILENTLY (no exception, no turns — ABT #850) when preset
            // sectional margins don't fit the bare-core window; several WE catalogue chokes
            // carry sheet-specified 3 mm spacers that fit the real part (wound on its case)
            // but not the modelled bare bore. The margins measurably do not change the full
            // model's capacitance (A/B over 107 measured chokes: bit-identical), so rather
            // than lose the part, retry on a fresh wind without the preset sections — and
            // only if THAT fails, throw.
            coil.set_sections_description(std::nullopt);
            coil.set_layers_description(std::nullopt);
            coil.wind();
            if (!coil.get_turns_description()) {
                throw InvalidInputException(ErrorCode::INVALID_INPUT,
                    "Impedance could not wind this coil: the winder returned no turns with the "
                    "preset sections and none on a fresh wind either (see ABT #850).");
            }
        }
        // Pass the core: the winding-to-core self term (ABT #848) needs it, and on
        // toroids it IS most of the self-capacitance.
        auto capacitanceMatrix = StrayCapacitance(strayCapacitanceModel).calculate_capacitance(coil, core).get_capacitance_among_windings().value();
        capacitance = capacitanceMatrix[coil.get_functional_description()[0].get_name()][coil.get_functional_description()[0].get_name()];
    }

    // The magnetizing tank models the COMMON-MODE measurement, which drives every winding
    // in parallel (that is how Zcm is measured: both winding pairs shorted together on each
    // side). Each winding then hangs its own self-capacitance across the same terminals, so
    // the tank capacitance is the SUM over windings — for a common-mode choke with two
    // equal windings, twice the single-winding value. The windings share one flux, so the
    // magnetizing inductance is NOT divided. The self-resonant-frequency full-model path
    // applies the same factor for a different physical reason: unity coupling mirrors the
    // driven winding's potential profile onto the open ones (see there). Verified on
    // 107 WE common-mode chokes: without this term the LC-governed families resonate a
    // consistent sqrt(2) high.
    size_t windingCount = coil.get_functional_description().size();
    if (windingCount > 1) {
        capacitance *= static_cast<double>(windingCount);
    }

    return ImpedanceTank{airCoredInductance, capacitance, true};
}


static double interpolate_permittivity_points(const MAS::Permittivity& data, double frequency) {
    // Log-log interpolation over the tabulated points; a single point is frequency-flat.
    // Clamps to the table ends: no extrapolation beyond the published data.
    if (std::holds_alternative<MAS::PermittivityPoint>(data)) {
        return std::get<MAS::PermittivityPoint>(data).get_value();
    }
    auto points = std::get<std::vector<MAS::PermittivityPoint>>(data);
    std::vector<std::pair<double, double>> table;
    for (auto& point : points) {
        if (point.get_frequency() && point.get_value() > 0) {
            table.push_back({point.get_frequency().value(), point.get_value()});
        }
    }
    if (table.empty()) {
        return points.empty() ? 0.0 : points[0].get_value();
    }
    std::sort(table.begin(), table.end());
    if (frequency <= table.front().first) return table.front().second;
    if (frequency >= table.back().first) return table.back().second;
    for (size_t i = 0; i + 1 < table.size(); ++i) {
        if (table[i].first <= frequency && frequency <= table[i + 1].first) {
            double t = std::log(frequency / table[i].first) / std::log(table[i + 1].first / table[i].first);
            return table[i].second * std::pow(table[i + 1].second / table[i].second, t);
        }
    }
    return table.back().second;
}

std::complex<double> Impedance::core_dimensional_attenuation(const CoreMaterial& material, double frequency, std::complex<double> complexPermeability, const std::vector<double>& crossSectionDimensions) {
    // ABT #848. A ferrite core is a lossy dielectric of enormous permittivity (MnZn: eps' ~ 1e5
    // at 1 MHz, Ferroxcube handbook Table 5) AND a conductor (MAS resistivity, eps'' = 1/(w eps0 rho)),
    // so the in-material wavelength is millimetres in the MHz band and the flux does not fill the
    // cross-section uniformly: the effective permeability of a slab of thickness d is
    // mu * tan(kd/2)/(kd/2) with k = w*sqrt(mu0 mu eps0 eps) (Snelling, Soft Ferrites, dimensional
    // resonance), complex because both mu and eps are. On a 15 mm MnZn cross-section this is
    // |factor| ~ 0.4 at 0.5 MHz and ~0.3 at 1 MHz — the simulated impedance of large MnZn chokes
    // was 2-4x high in exactly that band, which the stray-capacitance models had been blamed for.
    // The cross-section is solved as a 2D rectangle (exact double-series solution below). No
    // data -> factor 1: a material without permittivity gets no correction rather than an
    // invented one.
    if (!material.get_permittivity() || !material.get_permittivity()->get_complex() || crossSectionDimensions.empty()) {
        return 1.0;
    }
    auto complexPermittivity = material.get_permittivity()->get_complex().value();
    double epsReal = interpolate_permittivity_points(complexPermittivity.get_real(), frequency);
    double epsImag = interpolate_permittivity_points(complexPermittivity.get_imaginary(), frequency);
    if (epsReal <= 0) {
        return 1.0;
    }
    constexpr double vacuumPermeability = 4e-7 * std::numbers::pi;
    constexpr double vacuumPermittivity = 8.8541878128e-12;
    double angularFrequency = 2 * std::numbers::pi * frequency;
    // e^{jwt}: mu = mu' - j mu'', eps = eps' - j eps''
    std::complex<double> mu(complexPermeability.real(), -complexPermeability.imag());
    std::complex<double> eps(epsReal, -epsImag);
    std::complex<double> k = angularFrequency * std::sqrt(vacuumPermeability * mu * vacuumPermittivity * eps);
    // Average field over the cross-section with H = H0 on the boundary. One dimension (slab):
    // <H>/H0 = tan(kd/2)/(kd/2). Two dimensions a x b: the exact solution of d2H + k^2 H = 0
    // expands in a double sine series and gives
    //   <H>/H0 = 1 + sum_{m,n odd} 64 k^2 / (pi^4 m^2 n^2 (k_mn^2 - k^2)),  k_mn^2 = (m pi/a)^2 + (n pi/b)^2.
    // This is what a thin toroid needs: the field enters through the two faces a FEW mm apart
    // (the radial thickness), so the wave never has to cross the 15 mm height — the product of
    // two slab factors double-counted that and over-attenuated 2-5x (measured on the A05/A07
    // chokes before this form replaced it). The series converges as 1/(m^2 n^2); 41 odd terms
    // per axis is ample.
    std::vector<double> dims;
    for (double dimension : crossSectionDimensions) {
        if (dimension > 0) dims.push_back(dimension);
    }
    if (dims.empty()) {
        return 1.0;
    }
    std::complex<double> k2 = k * k;
    if (dims.size() == 1) {
        std::complex<double> x = k * dims[0] / 2.0;
        if (std::abs(x) < 1e-6) return 1.0;
        return std::tan(x) / x;
    }
    double a = dims[0], b = dims[1];
    std::complex<double> sum = 0.0;
    constexpr int maxOdd = 81;
    for (int m = 1; m <= maxOdd; m += 2) {
        for (int n = 1; n <= maxOdd; n += 2) {
            double kmn2 = std::pow(m * std::numbers::pi / a, 2) + std::pow(n * std::numbers::pi / b, 2);
            sum += 64.0 * k2 / (std::pow(std::numbers::pi, 4) * m * m * n * n * (kmn2 - k2));
        }
    }
    std::complex<double> factor = 1.0 + sum;
    return factor;
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
        // ABT #848: dimensional / eddy-dielectric attenuation across the core cross-section.
        // Complex multiply in e^{jwt}: (mu' - j mu'') * factor, then read back mu', mu''.
        auto factor = core_dimensional_attenuation(model.coreMaterial.value(), frequency, std::complex<double>(muReal, muImag), model.coreCrossSectionDimensions);
        std::complex<double> muEffective = std::complex<double>(muReal, -muImag) * factor;
        complexPermeabilityRealPart = muEffective.real() * model.permeabilityScaling;
        complexPermeabilityImaginaryPart = -muEffective.imag() * model.permeabilityScaling;
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
    if (core.get_processed_description() && !core.get_processed_description()->get_columns().empty()) {
        auto column = core.get_processed_description()->get_columns()[0];
        model.coreCrossSectionDimensions = {column.get_width(), column.get_depth()};
    }
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
    if (core.get_processed_description() && !core.get_processed_description()->get_columns().empty()) {
        auto column = core.get_processed_description()->get_columns()[0];
        model.coreCrossSectionDimensions = {column.get_width(), column.get_depth()};
    }
    model.temperature = temperature;
    model.tanks.push_back(build_magnetizing_tank(core, coil));
    return model;
}

std::complex<double> Impedance::calculate_impedance(Core core, Coil coil, double frequency, double temperature) {
    // Single-point common-mode/terminal impedance: just the magnetizing tank (the
    // first resonance). The wideband sweep adds the leakage tanks on top of this.
    WidebandImpedanceModel model;
    model.coreMaterial = core.resolve_material();
    if (core.get_processed_description() && !core.get_processed_description()->get_columns().empty()) {
        auto column = core.get_processed_description()->get_columns()[0];
        model.coreCrossSectionDimensions = {column.get_width(), column.get_depth()};
    }
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
    if (core.get_processed_description() && !core.get_processed_description()->get_columns().empty()) {
        auto column = core.get_processed_description()->get_columns()[0];
        model.coreCrossSectionDimensions = {column.get_width(), column.get_depth()};
    }
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
        capacitance = StrayCapacitanceOneLayer().calculate_capacitance(coil, core);
    }
    else {
        if (!coil.get_turns_description()) {
            coil.wind();
        }
        auto capacitanceMatrix = StrayCapacitance().calculate_capacitance(coil, core).get_capacitance_among_windings().value();

        capacitance = capacitanceMatrix[coil.get_functional_description()[0].get_name()][coil.get_functional_description()[0].get_name()];
        // An SRF measurement drives one winding, but the OTHERS ARE NOT ABSENT: unity
        // magnetic coupling forces every open winding to mirror the driven winding's
        // per-turn potential profile (a 1:1 transformer), so each one's turn-to-core
        // elements charge identically and add in parallel. Multiplying the single-winding
        // self term by the winding count is that mirror, and it lands both measured SRF
        // anchors (T12.5/7.5/5 110-turn: 180 kHz; Test_Impedance_0: 1.4 MHz) inside
        // tolerance where the bare single-winding term sat ~1.7x high (ABT #848).
        capacitance *= static_cast<double>(coil.get_functional_description().size());
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
