#include "converter_models/PushPull.h"
#include "converter_models/Topology.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>
#include <magic_enum.hpp>
#include "support/Exceptions.h"

namespace OpenMagnetics {

    // ----------------------------------------------------------------
    // SpiceSimulationConfig — central registry of per-topology defaults.
    //
    // Every topology that generates an ngspice netlist registers its
    // tuned values here. The struct's in-class defaults reflect a
    // generic "single-switch hard-switched DC-DC" baseline; topologies
    // diverge from those for documented physical reasons (snubber τ vs
    // Fs range, output cap vs power level, solver tolerances vs
    // switching events per period, etc.).
    //
    // Adding a new topology: insert a new entry below and override
    // `Topology::topology_kind()` in the concrete class.
    // ----------------------------------------------------------------
    const std::map<MAS::Topologies, SpiceSimulationConfig>&
    spice_simulation_defaults() {
        static const std::map<MAS::Topologies, SpiceSimulationConfig> defaults = [] {
            std::map<MAS::Topologies, SpiceSimulationConfig> m;

            // Boost — single-switch hard-switched, Fs 200 kHz – 2.2 MHz
            // (driven by reference-design EVMs TPS61089EVM-742,
            // TPS61178EVM-792, LM5122EVM-1PH). At 400 kHz a 1 µs snubber
            // τ would consume 40 % of the switching period and visibly
            // distort the inductor-current waveform, so the snubber is
            // sized for τ = 10 ns (negligible at all supported Fs).
            {
                SpiceSimulationConfig boost;
                boost.snubR = 100.0;       boost.snubC = 100e-12;
                boost.outputCapacitance = 100e-6;
                m[MAS::Topologies::BOOST_CONVERTER] = boost;
            }

            // Dab — dual active bridge, 8 switches in two H-bridges +
            // transformer, 4-corner phase-shift modulation. The original
            // hand-tuned values from Dab::generate_ngspice_circuit are
            // mirrored here verbatim so the registry is the single
            // source of truth. Snubber 1 kΩ/1 nF (τ = 1 µs) is fine at
            // the 50–250 kHz range Dab operates in. Solver tolerances
            // are looser than Boost's because the 4-position switching
            // produces 8 hard-switching events per period instead of 2.
            {
                SpiceSimulationConfig dab;
                dab.swModelVT = 2.5;       dab.swModelVH = 0.8;
                dab.swModelRON = 0.01;     dab.swModelROFF = 1e6;
                dab.snubR = 1e3;           dab.snubC = 1e-9;
                dab.diodeIS = 1e-12;       dab.diodeRS = 0.05;
                dab.outputCapacitance = 47e-6;
                dab.relTol = 0.01;         dab.absTol = 1e-7;
                dab.vnTol = 1e-4;
                dab.itl1 = 500;            dab.itl4 = 500;
                dab.method = "GEAR";       dab.trTol = 7.0;
                m[MAS::Topologies::DUAL_ACTIVE_BRIDGE_CONVERTER] = dab;
            }

            // Buck — single-switch hard-switched, sync or asynchronous
            // freewheeling diode. Reference EVMs span 400–500 kHz
            // (TPS54202EVM-716, LMR33630ADDAEVM, LM5146-Q1-EVM12V), so
            // the snubber τ = 10 ns mirrors Boost's choice. Output cap
            // 100 µF matches the integrator-style EVM output stage.
            {
                SpiceSimulationConfig buck;
                buck.snubR = 100.0;        buck.snubC = 100e-12;
                buck.outputCapacitance = 100e-6;
                m[MAS::Topologies::BUCK_CONVERTER] = buck;
            }

            // Flyback — single primary-switch isolated converter, energy
            // stored in coupled inductor airgap. Hard-switched at 50 kHz
            // – 200 kHz typical (current-mode controllers like UCC28C44,
            // LM5021, LM5180). Output capacitance 10 µF reflects the
            // smaller filter caps typical of low-power flybacks (vs the
            // 100 µF used in non-isolated buck/boost EVMs); the
            // demag-pulse-driven current ramp on the secondary doesn't
            // need as much bulk smoothing because peak current is
            // already metered by the controller. Snubber R = 100 kΩ
            // (× 100 pF, τ = 10 µs) intentionally large: with the
            // ideal-coupling K=1 transformer model leakage inductance is
            // zero so there is no spike to absorb, but a low snubR
            // (e.g. 1 kΩ) draws a continuous DC bias of tens of mA
            // from vin_dc into the primary winding sense during
            // off-time — the drain sits at -N·Vout (relative to
            // ground), so I_bias = (Vin + N·Vout) / Rsnub flows through
            // the primary current sense. On the low-power EVMs
            // (PMP30817 1.2 W, LM5180EVM 3 W) that is 30 % of peak
            // primary current and biases SPICE's primary-current
            // waveform away from the (correctly) zero-during-off
            // analytical model, polluting the NRMSE shape comparison.
            // 100 kΩ is high enough to drop the DC bias to <1 mA while
            // still presenting a small-signal damping path. GEAR +
            // TRTOL = 7 needed: the ideal-coupling Lp/Ls coupled-
            // inductor model produces near-step changes in iL at every
            // switching edge, which trapezoidal integration handles
            // poorly.
            {
                SpiceSimulationConfig flyback;
                flyback.swModelVT = 2.5;       flyback.swModelVH = 0.5;
                flyback.snubR = 1e5;           flyback.snubC = 100e-12;
                // Real-magnetic path retains the proven 1 kΩ damping for the
                // Lk·di/dt turn-off spike (Lk is non-zero with an actual
                // Magnetic component, vs the ideal-K=1 transformer used by
                // simulate_and_extract_topology_waveforms). 100 kΩ on the
                // real-magnetic path lets the spike run to ~12× Vin
                // (TestNgspiceRunner "Flyback ideal vs real magnetic
                // comparison" measures ~296 V at Vin = 24 V), which is the
                // *physically correct* uncladmped behaviour but well past
                // the 5× Vin envelope used to validate the design.
                flyback.snubRReal = 1e3;
                flyback.outputCapacitance = 10e-6;
                flyback.method = "GEAR";       flyback.trTol = 7.0;
                m[MAS::Topologies::FLYBACK_CONVERTER] = flyback;
            }

            // PushPull — center-tapped primary, two low-side switches 180°
            // out of phase, D ≤ 0.5. Reference designs span TI push-pull
            // transformer drivers SN6501 (410 kHz internal osc), SN6505B
            // (420 kHz), and SN6507 (programmable 100 kHz–2 MHz, wide-Vin
            // 3–36 V). Snubber R/C kept at 1 kΩ / 100 pF for symmetry
            // with Flyback/Dab but the *primary-switch* snubber is NOT
            // emitted by generate_ngspice_circuit() — see PushPull.cpp:
            // any cap from sw_node → 0 reflects through the K=0.9999
            // four-winding coupling into the secondary-diode
            // commutation path and breaks convergence at dsec_top. The
            // existing 1 MΩ "convergence helpers" on each winding
            // terminal plus the secondary RC diode snubbers (both
            // retained from the original hand-tuned netlist) provide
            // sufficient damping. snubR/snubC are still registered so
            // future MagneticAdviser passes that build a real magnetic
            // (with measurable leakage) can re-enable a primary snubber
            // without revisiting this comment. outputCapacitance =
            // 100 µF matches the LC-filter buck-derived secondary side
            // (forward-class topology, not energy-storage flyback).
            // Solver tolerances are looser than Boost/Buck (RELTOL 0.003,
            // ABSTOL 1e-7, VNTOL 1e-4, ITL1/ITL4 500/200) — matches
            // Dab's "many simultaneous switching events" tier: two
            // switches and two rectifier diodes commutate per Tsw on a
            // K=0.9999 multi-coupled transformer, which is too stiff for
            // Boost-tight tolerances. diodeRS=0.01 (vs Flyback's default
            // 1e-6) damps the rectifier reverse-recovery node enough to
            // clear "timestep too small" at commutation. METHOD=TRAP
            // (trapezoidal, ngspice default) mirrors the converter's
            // pre-registry hand-tuned netlist: GEAR was tried but the
            // simultaneous secondary-diode commutation against the
            // multi-coupled transformer pushes the GEAR step length
            // below 1 fs and aborts with "timestep too small"
            // (especially in step-up configurations like 24 V → 48 V).
            {
                SpiceSimulationConfig pushPull;
                pushPull.swModelVT = 2.5;      pushPull.swModelVH = 0.01;
                pushPull.snubR = 1e3;          pushPull.snubC = 100e-12;
                pushPull.diodeIS = 1e-14;      pushPull.diodeRS = 0.01;
                pushPull.outputCapacitance = 100e-6;
                pushPull.relTol = 0.003;       pushPull.absTol = 1e-7;
                pushPull.vnTol  = 1e-4;
                pushPull.itl1 = 500;           pushPull.itl4 = 200;
                pushPull.method = "TRAP";      pushPull.trTol = 7.0;
                m[MAS::Topologies::PUSH_PULL_CONVERTER] = pushPull;
            }

            // Cuk (V1 non-isolated, CCM) — fourth-order LC tank with hard
            // switching makes ngspice convergence stricter than Boost/Buck.
            // CUK_PLAN.md §6.2 mandates METHOD=GEAR + IC pre-charge for all
            // four reactive elements (L1, L2, C1, Co); the snubber τ is
            // sized for the 100 kHz–1.4 MHz fixture range (Erickson 100 kHz,
            // Bramble 300 kHz, LM2611 1.4 MHz). Snubber 100 Ω/100 pF
            // (τ = 10 ns) matches Boost; output cap is sized internally
            // per OP from coRipplePct so the registry default here is
            // unused (Cuk::generate_ngspice_circuit overrides it).
            {
                SpiceSimulationConfig cuk;
                cuk.swModelVT = 2.5;       cuk.swModelVH = 0.5;
                cuk.snubR = 100.0;         cuk.snubC = 100e-12;
                cuk.diodeIS = 1e-12;       cuk.diodeRS = 0.05;
                cuk.outputCapacitance = 10e-6;  // unused; Cuk overrides per-OP
                cuk.relTol = 0.01;         cuk.absTol = 1e-7;
                cuk.vnTol = 1e-4;
                cuk.itl1 = 500;            cuk.itl4 = 500;
                cuk.method = "GEAR";       cuk.trTol = 7.0;
                m[MAS::Topologies::CUK_CONVERTER] = cuk;
            }

            // SEPIC: same hard-switched, low-side-MOSFET shape as Cuk; Co
            // sized per OP from coRipplePct so the outputCapacitance default
            // is unused (Sepic::generate_ngspice_circuit overrides it).
            {
                SpiceSimulationConfig sepic;
                sepic.swModelVT = 2.5;       sepic.swModelVH = 0.5;
                sepic.snubR = 100.0;         sepic.snubC = 100e-12;
                sepic.diodeIS = 1e-12;       sepic.diodeRS = 0.05;
                sepic.outputCapacitance = 10e-6;
                sepic.relTol = 0.01;         sepic.absTol = 1e-7;
                sepic.vnTol = 1e-4;
                sepic.itl1 = 500;            sepic.itl4 = 500;
                sepic.method = "GEAR";       sepic.trTol = 7.0;
                m[MAS::Topologies::SEPIC_CONVERTER] = sepic;
            }

            // Zeta: same hard-switched 4th-order shape as SEPIC/Cuk but with a
            // high-side switch (modelled as an ideal SW between Vin and Vsw).
            // Co sized per OP from coRipplePct so the outputCapacitance default
            // is unused (Zeta::generate_ngspice_circuit overrides it).
            {
                SpiceSimulationConfig zeta;
                zeta.swModelVT = 2.5;        zeta.swModelVH = 0.5;
                zeta.snubR = 100.0;          zeta.snubC = 100e-12;
                zeta.diodeIS = 1e-12;        zeta.diodeRS = 0.05;
                zeta.outputCapacitance = 10e-6;
                zeta.relTol = 0.01;          zeta.absTol = 1e-7;
                zeta.vnTol = 1e-4;
                zeta.itl1 = 500;             zeta.itl4 = 500;
                zeta.method = "GEAR";        zeta.trTol = 7.0;
                m[MAS::Topologies::ZETA_CONVERTER] = zeta;
            }

            // Weinberg: isolated boost-derived topology with push-pull (V1) or
            // H-bridge (V2) primary, center-tapped forward rectifier on the
            // secondary, and (V1 only) D3a/D3b energy-recovery diodes from each
            // switch drain back to Vin to clamp leakage. SPICE convergence is
            // sensitive to the overlap (D > 0.5) interval where both primary
            // switches conduct simultaneously through the input inductor; GEAR
            // + ITL1=500 + UIC + IC pre-charge for L1a/L1b/Co are mandatory per
            // WEINBERG_PLAN.md §6.2. Co sized per OP from coRipplePct so the
            // outputCapacitance default here is unused (Weinberg overrides it).
            {
                SpiceSimulationConfig weinberg;
                weinberg.swModelVT = 2.5;       weinberg.swModelVH = 0.5;
                weinberg.snubR = 100.0;         weinberg.snubC = 100e-12;
                weinberg.diodeIS = 1e-12;       weinberg.diodeRS = 0.05;
                weinberg.outputCapacitance = 10e-6;
                weinberg.relTol = 0.01;         weinberg.absTol = 1e-7;
                weinberg.vnTol = 1e-4;
                weinberg.itl1 = 500;            weinberg.itl4 = 500;
                weinberg.method = "GEAR";       weinberg.trTol = 7.0;
                m[MAS::Topologies::WEINBERG_CONVERTER] = weinberg;
            }

            // 4-Switch Buck-Boost (FSBB): non-inverting H-bridge buck-boost with
            // single inductor and three operating regions (BUCK / BOOST /
            // BUCK_BOOST). Four MOSFETs (Q1/Q2 buck leg, Q3/Q4 boost leg) +
            // bootstrap drivers on the high sides + RC snubbers on both SW
            // nodes. The SPICE harness emits region-aware gate drives so the
            // netlist topology mirrors the analytical solver's region choice
            // (FOUR_SWITCH_BUCK_BOOST_PLAN.md §A.6).
            {
                SpiceSimulationConfig fsbb;
                fsbb.swModelVT = 2.5;            fsbb.swModelVH = 0.5;
                fsbb.snubR = 1000.0;             fsbb.snubC = 1e-9;
                fsbb.diodeIS = 1e-12;            fsbb.diodeRS = 0.05;
                fsbb.outputCapacitance = 47e-6;
                fsbb.relTol = 0.01;              fsbb.absTol = 1e-7;
                fsbb.vnTol = 1e-4;
                fsbb.itl1 = 500;                 fsbb.itl4 = 500;
                fsbb.method = "GEAR";            fsbb.trTol = 7.0;
                m[MAS::Topologies::FOUR_SWITCH_BUCK_BOOST_CONVERTER] = fsbb;
            }

            return m;
        }();
        return defaults;
    }

    SpiceSimulationConfig get_default_spice_config(MAS::Topologies t) {
        const auto& m = spice_simulation_defaults();
        auto it = m.find(t);
        if (it == m.end()) {
            throw std::invalid_argument(
                std::string("No SpiceSimulationConfig registered for topology ") +
                std::string(magic_enum::enum_name(t)) +
                ". Add an entry in Topology.cpp::spice_simulation_defaults().");
        }
        return it->second;
    }

    Topology::Topology(const json& j) {
        // Base class constructor - derived classes should handle their own JSON deserialization
    }

    OperatingPointExcitation Topology::complete_excitation(Waveform currentWaveform, Waveform voltageWaveform, double switchingFrequency, std::string name) {
        if (switchingFrequency <= 0 || !std::isfinite(switchingFrequency)) {
            throw std::invalid_argument("complete_excitation: Invalid switchingFrequency: " + std::to_string(switchingFrequency));
        }

        OperatingPointExcitation excitation;
        excitation.set_frequency(switchingFrequency);
        
        SignalDescriptor current;
        auto currentProcessed = Inputs::calculate_processed_data(currentWaveform, switchingFrequency, true);
        auto sampledCurrentWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(currentWaveform, switchingFrequency);
        auto currentHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency);
        // Store the resampled (power-of-2) waveform so downstream consumers
        // (MagnetizingInductance, harmonics derivation, FFT-based pipelines)
        // get the standardized waveform contract MKF expects. The raw
        // analytical waveform is not necessarily a power-of-2 length (e.g.
        // DAB/PSFB/PSHB use 2*N+1 = 513 samples, AHB has variable size due to
        // discontinuity duplicates), and the size-check gates in
        // MagnetizingInductance.cpp would otherwise throw.
        current.set_waveform(sampledCurrentWaveform);
        current.set_processed(currentProcessed);
        current.set_harmonics(currentHarmonics);
        excitation.set_current(current);
        SignalDescriptor voltage;
        auto voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, switchingFrequency, true);
        auto sampledVoltageWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(voltageWaveform, switchingFrequency);
        auto voltageHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency);
        voltage.set_waveform(sampledVoltageWaveform);
        voltage.set_processed(voltageProcessed);
        voltage.set_harmonics(voltageHarmonics);
        excitation.set_voltage(voltage);
        excitation.set_name(name);
        return excitation;
    }

    bool Topology::run_checks(bool assert) {
        throw NotImplementedException("Not implemented");
    }

    Inputs Topology::process() {
        run_checks(_assertErrors);

        Inputs inputs;
        auto designRequirements = process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto turnsRatio : designRequirements.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(turnsRatio));
        }
        auto desiredMagnetizingInductance = resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        auto operatingPoints = process_operating_points(turnsRatios, desiredMagnetizingInductance);
        for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
            operatingPoints[operatingPointIndex] = Inputs::prune_harmonics(operatingPoints[operatingPointIndex], Defaults().harmonicAmplitudeThreshold);
        }

        inputs.set_design_requirements(designRequirements);
        inputs.set_operating_points(operatingPoints);

        return inputs;
    }


    std::vector<std::variant<Inputs, CAS::Inputs>> Topology::get_extra_components_inputs(
        ExtraComponentsMode mode,
        std::optional<Magnetic> magnetic)
    {
        if (mode == ExtraComponentsMode::REAL && !magnetic.has_value()) {
            throw std::invalid_argument("get_extra_components_inputs: mode REAL requires a designed magnetic");
        }
        return {};
    }

} // namespace OpenMagnetics
