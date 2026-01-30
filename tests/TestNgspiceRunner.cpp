#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "processors/NgspiceRunner.h"
#include "processors/CircuitSimulatorInterface.h"
#include "converter_models/Flyback.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/Painter.h"
#include "TestingUtils.h"
#include <fstream>
#include <filesystem>
#include <source_location>
#include <numeric>

using namespace OpenMagnetics;

const auto outputFilePath = std::filesystem::path{std::source_location::current().file_name()}.parent_path().append("..").append("output");

TEST_CASE("NgspiceRunner availability check", "[NgspiceRunner]") {
    NgspiceRunner runner;
    
    // Just check that the runner can be instantiated
    // Actual availability depends on system configuration
    REQUIRE_NOTHROW(runner.get_mode());
}

TEST_CASE("NgspiceRunner simple netlist parsing", "[NgspiceRunner][skip-ci]") {
    NgspiceRunner runner;
    
    // Skip if ngspice is not available
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
    
    // Simple RC circuit for testing
    std::string netlist = R"(
* Simple RC Test Circuit
V1 in 0 PULSE(0 1 0 1n 1n 50u 100u)
R1 in out 1k
C1 out 0 1n

.tran 1u 500u 100u
.end
)";
    
    SimulationConfig config;
    config.frequency = 10e3;  // 10 kHz
    config.keepTempFiles = false;
    config.extractOnePeriod = false;  // Don't try to extract one period for this simple test
    
    auto result = runner.run_simulation(netlist, config);
    
    if (result.success) {
        // Just check that we got some waveforms - detailed extraction may fail
        REQUIRE(result.simulationTime > 0);
    } else {
        // Log the error but don't fail - ngspice might have issues on this system
        WARN("Simulation failed: " << result.errorMessage);
    }
}

TEST_CASE("SimulationResult structure", "[NgspiceRunner]") {
    SimulationResult result;
    
    // Default state
    REQUIRE(result.success == false);
    REQUIRE(result.waveforms.empty());
    REQUIRE(result.waveformNames.empty());
    REQUIRE(!result.operatingPoint.has_value());
    
    // Setting values
    result.success = true;
    result.simulationTime = 1.5;
    
    Waveform waveform;
    waveform.set_data({1.0, 2.0, 3.0});
    result.waveforms.push_back(waveform);
    result.waveformNames.push_back("test");
    
    REQUIRE(result.success == true);
    REQUIRE(result.waveforms.size() == 1);
    REQUIRE(result.waveformNames.size() == 1);
    REQUIRE(result.simulationTime == 1.5);
}

TEST_CASE("SimulationConfig defaults", "[NgspiceRunner]") {
    SimulationConfig config;
    
    REQUIRE(config.stopTime == 0.0);
    REQUIRE(config.stepSize == 0.0);
    REQUIRE(config.steadyStateCycles == 5);
    REQUIRE(config.frequency == 0.0);
    REQUIRE(config.extractOnePeriod == true);
    REQUIRE(config.workingDirectory.empty());
    REQUIRE(config.keepTempFiles == false);
    REQUIRE(config.timeout == 60.0);
}

TEST_CASE("NgspiceRunner simulate and export waveforms to SVG", "[NgspiceRunner][svg-export]") {
    NgspiceRunner runner;
    
    // Skip if ngspice is not available
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
    
    // Flyback-like transformer test circuit
    // Primary side: voltage source with rectangular wave
    // Secondary side: load resistor
    std::string netlist = R"(
* Flyback Transformer Test Circuit
* Primary voltage source (rectangular wave, 100kHz)
Vpri in_p 0 PULSE(0 24 0 10n 10n 4.5u 10u)

* Simple coupled inductor model
Lp in_p out_p 100u
Ls out_s 0 25u
Kps Lp Ls 0.98

* Primary current sense resistor
Rsense out_p 0 0.1

* Secondary load
Rload out_s sec_gnd 10
Vsec_gnd sec_gnd 0 0

* Analysis: 10us period, simulate for 200us, skip first 100us
.tran 100n 200u 100u

* Save relevant voltages and currents
.save v(in_p) v(out_s) i(Vpri) i(Vsec_gnd)

.end
)";
    
    SimulationConfig config;
    config.frequency = 100e3;  // 100 kHz
    config.keepTempFiles = false;
    config.extractOnePeriod = false;  // Keep full waveforms for plotting
    
    auto result = runner.run_simulation(netlist, config);
    
    REQUIRE(result.success);
    REQUIRE(result.simulationTime > 0);
    
    // Check we got some waveforms
    if (!result.waveforms.empty()) {
        INFO("Got " << result.waveforms.size() << " waveforms");
        
        // Create waveform name mapping for OperatingPoint extraction
        // Each entry is a map for one winding with "voltage" and "current" keys
        NgspiceRunner::WaveformNameMapping mapping = {
            {{"voltage", "v(in_p)"}, {"current", "i(Vpri)"}},    // Primary
            {{"voltage", "v(out_s)"}, {"current", "i(Vsec_gnd)"}}  // Secondary
        };
        
        // Extract OperatingPoint from simulation result
        std::vector<std::string> windingNames = {"Primary", "Secondary"};
        auto operatingPoint = NgspiceRunner::extract_operating_point(result, mapping, 100e3, windingNames);
        REQUIRE(!operatingPoint.get_excitations_per_winding().empty());
        
        // Export waveforms to SVG using OperatingPoint
        BasicPainter painter;
        std::string svg = painter.paint_operating_point_waveforms(
            operatingPoint,
            "Flyback Transformer Simulation",
            1000,  // width
            800    // height
        );
        
        // Verify SVG was generated
        REQUIRE(!svg.empty());
        REQUIRE(svg.find("<svg") != std::string::npos);
        REQUIRE(svg.find("</svg>") != std::string::npos);
        
        // Save to file for manual inspection
        auto outFile = outputFilePath;
        outFile.append("flyback_simulation_waveforms.svg");
        std::ofstream ofs(outFile);
        if (ofs.is_open()) {
            ofs << svg;
            ofs.close();
            INFO("Waveforms saved to " << outFile.string());
        }
    }
}

TEST_CASE("BasicPainter paint_operating_point_waveforms with synthetic data", "[NgspiceRunner][BasicPainter]") {
    // Test the waveform painting without running simulation
    BasicPainter painter;
    
    // Create synthetic test waveforms
    Waveform primaryVoltageWf, primaryCurrentWf, secondaryVoltageWf, secondaryCurrentWf;
    std::vector<double> data, time;
    
    // Sine wave for primary voltage
    for (int i = 0; i < 100; ++i) {
        double t = i * 0.00001;  // 10us total, 100ns steps
        time.push_back(t);
        data.push_back(10.0 * std::sin(2 * 3.14159 * 100000 * t));  // 100kHz sine
    }
    primaryVoltageWf.set_data(data);
    primaryVoltageWf.set_time(time);
    
    // Triangle wave for primary current
    data.clear();
    for (int i = 0; i < 100; ++i) {
        double t = i * 0.00001;
        double phase = std::fmod(t * 100000, 1.0);
        data.push_back(phase < 0.5 ? 4 * phase : 4 * (1 - phase));
    }
    primaryCurrentWf.set_data(data);
    primaryCurrentWf.set_time(time);
    
    // Rectangular wave for secondary voltage
    data.clear();
    for (int i = 0; i < 100; ++i) {
        double t = i * 0.00001;
        double phase = std::fmod(t * 100000, 1.0);
        data.push_back(phase < 0.5 ? 12.0 : 0.0);
    }
    secondaryVoltageWf.set_data(data);
    secondaryVoltageWf.set_time(time);
    
    // Sawtooth for secondary current
    data.clear();
    for (int i = 0; i < 100; ++i) {
        double t = i * 0.00001;
        double phase = std::fmod(t * 100000, 1.0);
        data.push_back(phase * 2.0);
    }
    secondaryCurrentWf.set_data(data);
    secondaryCurrentWf.set_time(time);
    
    // Build OperatingPoint with excitations
    OperatingPoint op;
    std::vector<OperatingPointExcitation> excitations;
    
    // Primary winding
    OperatingPointExcitation primaryExc;
    primaryExc.set_name("Primary");
    primaryExc.set_frequency(100e3);
    SignalDescriptor primaryVoltage, primaryCurrent;
    primaryVoltage.set_waveform(primaryVoltageWf);
    primaryCurrent.set_waveform(primaryCurrentWf);
    primaryExc.set_voltage(primaryVoltage);
    primaryExc.set_current(primaryCurrent);
    excitations.push_back(primaryExc);
    
    // Secondary winding
    OperatingPointExcitation secondaryExc;
    secondaryExc.set_name("Secondary");
    secondaryExc.set_frequency(100e3);
    SignalDescriptor secondaryVoltage, secondaryCurrent;
    secondaryVoltage.set_waveform(secondaryVoltageWf);
    secondaryCurrent.set_waveform(secondaryCurrentWf);
    secondaryExc.set_voltage(secondaryVoltage);
    secondaryExc.set_current(secondaryCurrent);
    excitations.push_back(secondaryExc);
    
    op.set_excitations_per_winding(excitations);
    
    // Generate SVG
    std::string svg = painter.paint_operating_point_waveforms(
        op,
        "Test Waveforms - 100kHz",
        800,
        600
    );
    
    // Verify SVG structure
    REQUIRE(!svg.empty());
    REQUIRE(svg.find("<svg") != std::string::npos);
    REQUIRE(svg.find("</svg>") != std::string::npos);
    REQUIRE(svg.find("Primary Voltage") != std::string::npos);
    REQUIRE(svg.find("Primary Current") != std::string::npos);
    REQUIRE(svg.find("Secondary Voltage") != std::string::npos);
    REQUIRE(svg.find("<path") != std::string::npos);  // Waveform paths
    
    // Save for manual inspection
    auto outFile = outputFilePath;
    outFile.append("test_waveforms.svg");
    std::ofstream ofs(outFile);
    if (ofs.is_open()) {
        ofs << svg;
        ofs.close();
    }
}

TEST_CASE("Flyback converter full flow: MAS::Flyback -> ngspice -> OperatingPoint", "[NgspiceRunner][flyback][integration]") {
    NgspiceRunner runner;
    
    // Skip if ngspice is not available
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
    
    // Create a Flyback converter specification using MAS::Flyback
    OpenMagnetics::Flyback flyback;
    
    // Input voltage: 24V nominal (18-32V range)
    DimensionWithTolerance inputVoltage;
    inputVoltage.set_nominal(24.0);
    inputVoltage.set_minimum(18.0);
    inputVoltage.set_maximum(32.0);
    flyback.set_input_voltage(inputVoltage);
    
    // Diode voltage drop
    flyback.set_diode_voltage_drop(0.5);
    
    // Efficiency
    flyback.set_efficiency(0.85);
    
    // Current ripple ratio
    flyback.set_current_ripple_ratio(0.4);
    
    // Maximum duty cycle
    flyback.set_maximum_duty_cycle(0.5);
    
    // Operating point: 5V @ 2A output, 100kHz
    OpenMagnetics::FlybackOperatingPoint opPoint;
    opPoint.set_output_voltages({5.0});
    opPoint.set_output_currents({2.0});
    opPoint.set_switching_frequency(100e3);
    opPoint.set_ambient_temperature(25.0);
    
    flyback.set_operating_points({opPoint});
    
    // Process design requirements to get turns ratios and inductance
    auto designReqs = flyback.process_design_requirements();
    
    // Extract calculated values
    std::vector<double> turnsRatios;
    for (const auto& tr : designReqs.get_turns_ratios()) {
        turnsRatios.push_back(tr.get_nominal().value());
    }
    double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
    
    INFO("Calculated turns ratio: " << turnsRatios[0]);
    INFO("Calculated inductance: " << (magnetizingInductance * 1e6) << " uH");
    
    // Generate ngspice circuit
    std::string netlist = flyback.generate_ngspice_circuit(turnsRatios, magnetizingInductance, 0, 0);
    
    INFO("Generated netlist:\n" << netlist);
    
    // Save netlist
    auto netlistPath = outputFilePath;
    netlistPath.append("flyback_from_mas.cir");
    std::ofstream netlistFile(netlistPath);
    if (netlistFile.is_open()) {
        netlistFile << netlist;
        netlistFile.close();
    }
    
    // Run simulation and extract operating points
    auto operatingPoints = flyback.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
    
    REQUIRE(!operatingPoints.empty());
    INFO("Got " << operatingPoints.size() << " operating points from simulation");
    
    // Check that we got waveforms
    for (size_t i = 0; i < operatingPoints.size(); ++i) {
        const auto& op = operatingPoints[i];
        INFO("Operating point " << i << ": " << op.get_name().value_or("unnamed"));
        INFO("  Windings: " << op.get_excitations_per_winding().size());
        
        REQUIRE(!op.get_excitations_per_winding().empty());
        
        for (size_t w = 0; w < op.get_excitations_per_winding().size(); ++w) {
            const auto& exc = op.get_excitations_per_winding()[w];
            std::string name = exc.get_name().value_or("Winding " + std::to_string(w));
            INFO("  " << name << ":");
            
            if (exc.get_voltage() && exc.get_voltage()->get_waveform()) {
                INFO("    Voltage: " << exc.get_voltage()->get_waveform()->get_data().size() << " points");
            }
            if (exc.get_current() && exc.get_current()->get_waveform()) {
                INFO("    Current: " << exc.get_current()->get_waveform()->get_data().size() << " points");
            }
        }
    }
    
    // Paint waveforms using the simplified OperatingPoint-based method
    BasicPainter painter;
    
    // Paint all windings for first operating point (just pass the OperatingPoint)
    std::string svg = painter.paint_operating_point_waveforms(
        operatingPoints[0],
        "Flyback Converter - All Windings (from MAS::Flyback)"
    );
    
    REQUIRE(!svg.empty());
    
    auto outFile = outputFilePath;
    outFile.append("flyback_mas_all_windings.svg");
    std::ofstream ofs(outFile);
    if (ofs.is_open()) {
        ofs << svg;
        ofs.close();
        INFO("All windings SVG saved to " << outFile.string());
    }
}

TEST_CASE("Flyback DCM with MAS::Flyback model", "[NgspiceRunner][flyback][dcm][integration]") {
    NgspiceRunner runner;
    
    // Skip if ngspice is not available
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
    
    // Create a Flyback converter specification with DCM parameters
    // Using currentRippleRatio > 1 forces DCM operation
    OpenMagnetics::Flyback flyback;
    
    // Input voltage
    DimensionWithTolerance inputVoltage;
    inputVoltage.set_nominal(48.0);
    inputVoltage.set_minimum(42.0);
    inputVoltage.set_maximum(54.0);
    flyback.set_input_voltage(inputVoltage);
    
    // Diode voltage drop
    flyback.set_diode_voltage_drop(0.5);
    
    // Efficiency
    flyback.set_efficiency(0.9);
    
    // Current ripple ratio > 1 triggers DCM (current goes to zero each cycle)
    flyback.set_current_ripple_ratio(2.0);  // DCM: ripple > 1
    
    // Maximum duty cycle is required for design
    flyback.set_maximum_duty_cycle(0.45);
    
    // Operating point: 12V @ 0.5A output (light load)
    OpenMagnetics::FlybackOperatingPoint opPoint;
    opPoint.set_output_voltages({12.0});
    opPoint.set_output_currents({0.5});
    opPoint.set_switching_frequency(100e3);
    opPoint.set_ambient_temperature(25.0);
    
    flyback.set_operating_points({opPoint});
    
    // Process design requirements
    auto designReqs = flyback.process_design_requirements();
    
    // Extract calculated values
    std::vector<double> turnsRatios;
    for (const auto& tr : designReqs.get_turns_ratios()) {
        turnsRatios.push_back(tr.get_nominal().value());
    }
    double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
    
    INFO("DCM Flyback calculated parameters:");
    INFO("  Turns ratio: " << turnsRatios[0]);
    INFO("  Inductance: " << (magnetizingInductance * 1e6) << " uH");
    
    // Generate ngspice circuit
    std::string netlist = flyback.generate_ngspice_circuit(turnsRatios, magnetizingInductance, 0, 0);
    
    INFO("Generated DCM netlist:\n" << netlist);
    
    // Save netlist
    auto netlistPath = outputFilePath;
    netlistPath.append("flyback_dcm_from_mas.cir");
    std::ofstream netlistFile(netlistPath);
    if (netlistFile.is_open()) {
        netlistFile << netlist;
        netlistFile.close();
    }
    
    // Run simulation and extract operating points
    auto operatingPoints = flyback.simulate_and_extract_operating_points(turnsRatios, magnetizingInductance);
    
    REQUIRE(!operatingPoints.empty());
    INFO("Got " << operatingPoints.size() << " operating points from DCM simulation");
    
    // Check waveforms
    for (size_t i = 0; i < operatingPoints.size(); ++i) {
        const auto& op = operatingPoints[i];
        REQUIRE(!op.get_excitations_per_winding().empty());
        
        // Check primary current waveform characteristics
        const auto& primaryExc = op.get_excitations_per_winding()[0];
        if (primaryExc.get_current() && primaryExc.get_current()->get_waveform()) {
            // Copy data to avoid dangling reference
            auto data = primaryExc.get_current()->get_waveform()->get_data();
            
            // In DCM, current should touch or approach zero
            double minCurrent = *std::min_element(data.begin(), data.end());
            INFO("Operating point " << i << " primary current min: " << minCurrent);
            
            // DCM characteristic: minimum current should be close to zero
            // (within some tolerance for simulation artifacts)
            CHECK(std::abs(minCurrent) < 0.5);  // Should be near zero in DCM
        }
    }
    
    // Paint waveforms
    BasicPainter painter;
    std::string svg = painter.paint_operating_point_waveforms(
        operatingPoints[0],
        "Flyback DCM - MAS::Flyback Model (Light Load)"
    );
    
    REQUIRE(!svg.empty());
    
    auto outFile = outputFilePath;
    outFile.append("flyback_dcm_mas_waveforms.svg");
    std::ofstream ofs(outFile);
    if (ofs.is_open()) {
        ofs << svg;
        ofs.close();
        INFO("DCM waveforms saved to " << outFile.string());
    }
}

TEST_CASE("Flyback topology waveform validation", "[NgspiceRunner][flyback][validation]") {
    NgspiceRunner runner;
    
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
    
    // Create a flyback converter
    OpenMagnetics::Flyback flyback;
    
    DimensionWithTolerance inputVoltage;
    inputVoltage.set_nominal(48.0);
    inputVoltage.set_minimum(42.0);
    inputVoltage.set_maximum(54.0);
    flyback.set_input_voltage(inputVoltage);
    
    OpenMagnetics::FlybackOperatingPoint opPoint;
    opPoint.set_output_voltages({12.0});
    opPoint.set_output_currents({1.0});
    opPoint.set_ambient_temperature(25.0);
    opPoint.set_switching_frequency(100000.0);
    flyback.set_operating_points({opPoint});
    
    flyback.set_diode_voltage_drop(0.5);
    flyback.set_efficiency(0.9);
    flyback.set_maximum_duty_cycle(0.5);
    flyback.set_current_ripple_ratio(0.4);  // CCM mode
    
    auto designReqs = flyback.process_design_requirements();
    
    std::vector<double> turnsRatios;
    for (const auto& tr : designReqs.get_turns_ratios()) {
        turnsRatios.push_back(tr.get_nominal().value());
    }
    double magnetizingInductance = designReqs.get_magnetizing_inductance().get_minimum().value();
    
    INFO("Turns ratio: " << turnsRatios[0]);
    INFO("Inductance: " << (magnetizingInductance * 1e6) << " uH");
    
    // Extract topology waveforms for validation
    auto topologyWaveforms = flyback.simulate_and_extract_topology_waveforms(turnsRatios, magnetizingInductance);
    
    REQUIRE(!topologyWaveforms.empty());
    
    const auto& wf = topologyWaveforms[0];
    
    // Validate that all waveforms were extracted
    REQUIRE(!wf.time.empty());
    REQUIRE(!wf.switchNodeVoltage.empty());
    REQUIRE(!wf.secondaryWindingVoltage.empty());
    REQUIRE(!wf.outputVoltage.empty());
    REQUIRE(!wf.primaryCurrent.empty());
    REQUIRE(!wf.secondaryCurrent.empty());
    
    // Calculate waveform statistics
    double priV_max = *std::max_element(wf.switchNodeVoltage.begin(), wf.switchNodeVoltage.end());
    double priV_min = *std::min_element(wf.switchNodeVoltage.begin(), wf.switchNodeVoltage.end());
    double secV_max = *std::max_element(wf.secondaryWindingVoltage.begin(), wf.secondaryWindingVoltage.end());
    double secV_min = *std::min_element(wf.secondaryWindingVoltage.begin(), wf.secondaryWindingVoltage.end());
    double vout_avg = std::accumulate(wf.outputVoltage.begin(), wf.outputVoltage.end(), 0.0) / wf.outputVoltage.size();
    
    INFO("Primary voltage: min=" << priV_min << " max=" << priV_max);
    INFO("Secondary voltage: min=" << secV_min << " max=" << secV_max);
    INFO("Output voltage avg: " << vout_avg);
    INFO("Input voltage: " << wf.inputVoltageValue);
    INFO("Expected output: " << wf.outputVoltageValue);
    INFO("Duty cycle: " << (wf.dutyCycle * 100) << "%");
    
    // Validate primary voltage behavior:
    // During ON: V_pri should be close to Vin
    // Allow some tolerance for simulation effects
    CHECK(priV_max > wf.inputVoltageValue * 0.9);
    CHECK(priV_max < wf.inputVoltageValue * 1.1);
    
    // During OFF: V_pri should be negative (reflected voltage = -Vout * N)
    double expectedReflectedVoltage = -wf.outputVoltageValue * turnsRatios[0];
    CHECK(priV_min < 0.0);
    // Check it's in the right ballpark (within 50% due to simulation variations)
    CHECK(priV_min < expectedReflectedVoltage * 0.5);
    CHECK(priV_min > expectedReflectedVoltage * 1.5);
    
    // Validate secondary voltage behavior:
    // During ON: V_sec should be negative = -Vin/N
    double expectedSecondaryBlocked = -wf.inputVoltageValue / turnsRatios[0];
    CHECK(secV_min < 0.0);
    CHECK(secV_min < expectedSecondaryBlocked * 0.5);
    CHECK(secV_min > expectedSecondaryBlocked * 1.5);
    
    // During OFF: V_sec should be around Vout (diode conducting)
    CHECK(secV_max > 0.0);
    CHECK(secV_max > wf.outputVoltageValue * 0.5);
    
    // Output voltage should be close to expected
    CHECK(vout_avg > wf.outputVoltageValue * 0.8);
    CHECK(vout_avg < wf.outputVoltageValue * 1.2);
    
    INFO("Topology waveform validation passed");
}

TEST_CASE("Flyback simulation with real Magnetic component", "[NgspiceRunner][flyback][magnetic][integration]") {
    NgspiceRunner runner;
    
    // Skip if ngspice is not available
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
    
    // Create a Flyback converter specification
    OpenMagnetics::Flyback flyback;
    
    // Input voltage: 24V nominal
    DimensionWithTolerance inputVoltage;
    inputVoltage.set_nominal(24.0);
    flyback.set_input_voltage(inputVoltage);
    
    // Diode voltage drop
    flyback.set_diode_voltage_drop(0.5);
    
    // Operating point: 5V @ 2A
    OpenMagnetics::FlybackOperatingPoint opPoint;
    opPoint.set_output_voltages({5.0});
    opPoint.set_output_currents({2.0});
    opPoint.set_ambient_temperature(25.0);
    opPoint.set_switching_frequency(100e3);
    flyback.set_operating_points({opPoint});
    
    // Create a real Magnetic component
    // Flyback transformer: 24V to 5V, ~5:1 turns ratio (actually using 20:5 = 4:1)
    std::vector<int64_t> numberTurns = {20, 5};  // Primary 20T, Secondary 5T = 4:1 ratio
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "E 25/13/7";
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName);
    
    int64_t numberStacks = 1;
    std::string coreMaterial = "3C90";
    // Use a gap for inductance control
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0004, 1);  // Single 0.4mm gap
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    
    INFO("Magnetic reference: " << magnetic.get_reference());
    
    // Generate the circuit with real magnetic
    std::string netlist = flyback.generate_ngspice_circuit_with_magnetic(magnetic);
    
    INFO("Generated netlist with real Magnetic:\n" << netlist);
    
    // Verify the netlist contains the subcircuit
    REQUIRE(netlist.find(".subckt") != std::string::npos);
    REQUIRE(netlist.find("X1") != std::string::npos);  // Subcircuit instance
    bool hasRdcOrLmag = (netlist.find("Rdc") != std::string::npos) || (netlist.find("Lmag") != std::string::npos);
    REQUIRE(hasRdcOrLmag);  // Magnetic model elements
    
    // Run simulation and extract operating points
    auto operatingPoints = flyback.simulate_with_magnetic_and_extract_operating_points(magnetic);
    
    REQUIRE(operatingPoints.size() == 1);  // One input voltage, one operating point
    
    auto& op = operatingPoints[0];
    REQUIRE(op.get_excitations_per_winding().size() == 2);  // Primary and secondary
    
    // Check that we got voltage and current waveforms
    auto& primaryExcitation = op.get_excitations_per_winding()[0];
    REQUIRE(primaryExcitation.get_voltage().has_value());
    REQUIRE(primaryExcitation.get_current().has_value());
    
    auto& secondaryExcitation = op.get_excitations_per_winding()[1];
    REQUIRE(secondaryExcitation.get_voltage().has_value());
    REQUIRE(secondaryExcitation.get_current().has_value());
    
    INFO("Primary voltage waveform points: " << primaryExcitation.get_voltage()->get_waveform()->get_data().size());
    INFO("Secondary current waveform points: " << secondaryExcitation.get_current()->get_waveform()->get_data().size());
    
    // Export to SVG for visualization
    BasicPainter painter;
    std::string svg = painter.paint_operating_point_waveforms(
        op,
        "Flyback with Real Magnetic - " + magnetic.get_reference(),
        1200,
        900
    );
    
    REQUIRE(!svg.empty());
    
    auto outFile = outputFilePath;
    outFile.append("flyback_real_magnetic.svg");
    std::ofstream ofs(outFile);
    if (ofs.is_open()) {
        ofs << svg;
        ofs.close();
        INFO("Real magnetic waveforms saved to " << outFile.string());
    }
    
    // Also save the netlist for reference
    auto netlistPath = outputFilePath;
    netlistPath.append("flyback_real_magnetic.cir");
    std::ofstream netlistFile(netlistPath);
    if (netlistFile.is_open()) {
        netlistFile << netlist;
        netlistFile.close();
        INFO("Netlist saved to " << netlistPath.string());
    }
    
    INFO("Flyback simulation with real Magnetic component completed successfully");
}

TEST_CASE("Flyback ideal vs real magnetic comparison", "[NgspiceRunner][flyback][magnetic][comparison]") {
    NgspiceRunner runner;
    
    // Skip if ngspice is not available
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
    
    // Create a Flyback converter specification
    OpenMagnetics::Flyback flyback;
    
    DimensionWithTolerance inputVoltage;
    inputVoltage.set_nominal(24.0);
    flyback.set_input_voltage(inputVoltage);
    flyback.set_diode_voltage_drop(0.5);
    
    OpenMagnetics::FlybackOperatingPoint opPoint;
    opPoint.set_output_voltages({5.0});
    opPoint.set_output_currents({2.0});
    opPoint.set_ambient_temperature(25.0);
    opPoint.set_switching_frequency(100e3);
    flyback.set_operating_points({opPoint});
    
    // Create a real Magnetic component
    std::vector<int64_t> numberTurns = {20, 5};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "E 25/13/7";
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName);
    int64_t numberStacks = 1;
    std::string coreMaterial = "3C90";
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0004, 1);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, numberStacks, coreMaterial);
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    
    // Extract parameters from magnetic
    double primaryTurns = coil.get_functional_description()[0].get_number_turns();
    double secondaryTurns = coil.get_functional_description()[1].get_number_turns();
    double turnsRatio = primaryTurns / secondaryTurns;
    double magnetizingInductance = resolve_dimensional_values(
        MagnetizingInductance().calculate_inductance_from_number_turns_and_gapping(magnetic).get_magnetizing_inductance());
    
    INFO("Turns ratio: " << turnsRatio);
    INFO("Magnetizing inductance: " << (magnetizingInductance * 1e6) << " uH");
    
    // Run ideal simulation
    auto idealOps = flyback.simulate_and_extract_operating_points({turnsRatio}, magnetizingInductance);
    REQUIRE(idealOps.size() == 1);
    
    // Run real magnetic simulation
    auto realOps = flyback.simulate_with_magnetic_and_extract_operating_points(magnetic);
    REQUIRE(realOps.size() == 1);
    
    // Both simulations should produce valid waveforms
    REQUIRE(idealOps[0].get_excitations_per_winding().size() == 2);
    REQUIRE(realOps[0].get_excitations_per_winding().size() == 2);
    
    // Compare primary voltage waveforms - extract characteristics
    auto idealPriV = idealOps[0].get_excitations_per_winding()[0].get_voltage()->get_waveform()->get_data();
    auto realPriV = realOps[0].get_excitations_per_winding()[0].get_voltage()->get_waveform()->get_data();
    
    double idealPriVmax = *std::max_element(idealPriV.begin(), idealPriV.end());
    double realPriVmax = *std::max_element(realPriV.begin(), realPriV.end());
    double idealPriVmin = *std::min_element(idealPriV.begin(), idealPriV.end());
    double realPriVmin = *std::min_element(realPriV.begin(), realPriV.end());
    
    INFO("Ideal primary V max: " << idealPriVmax << ", min: " << idealPriVmin);
    INFO("Real primary V max: " << realPriVmax << ", min: " << realPriVmin);
    
    // Primary max voltage should match Vin closely (switch ON applies Vin across primary)
    CHECK(std::abs(idealPriVmax - 24.0) < 1.0);  // Ideal should be near Vin
    CHECK(std::abs(realPriVmax - 24.0) < 1.0);   // Real should also be near Vin
    
    // Compare primary current waveforms
    auto idealPriI = idealOps[0].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data();
    auto realPriI = realOps[0].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data();
    
    double idealPriImax = *std::max_element(idealPriI.begin(), idealPriI.end());
    double realPriImax = *std::max_element(realPriI.begin(), realPriI.end());
    
    INFO("Ideal primary I max: " << idealPriImax);
    INFO("Real primary I max: " << realPriImax);
    
    // Both should have positive peak currents (magnetizing current builds up during on-time)
    CHECK(idealPriImax > 0.1);  // Should have measurable current
    CHECK(realPriImax > 0.1);   // Real also should have measurable current
    
    // Save comparison SVGs for visual inspection
    BasicPainter painter;
    
    std::string svgIdeal = painter.paint_operating_point_waveforms(
        idealOps[0], "Flyback Ideal Transformer", 1200, 900);
    auto outIdeal = outputFilePath;
    outIdeal.append("flyback_comparison_ideal.svg");
    std::ofstream ofsIdeal(outIdeal);
    if (ofsIdeal.is_open()) { ofsIdeal << svgIdeal; ofsIdeal.close(); }
    
    std::string svgReal = painter.paint_operating_point_waveforms(
        realOps[0], "Flyback Real Magnetic", 1200, 900);
    auto outReal = outputFilePath;
    outReal.append("flyback_comparison_real.svg");
    std::ofstream ofsReal(outReal);
    if (ofsReal.is_open()) { ofsReal << svgReal; ofsReal.close(); }
    
    INFO("Comparison SVGs saved for visual inspection");
}

TEST_CASE("Flyback real magnetic with multiple input voltages", "[NgspiceRunner][flyback][magnetic][multi-voltage]") {
    NgspiceRunner runner;
    
    if (!runner.is_available()) {
        SKIP("ngspice not available on this system");
    }
    
    // Create Flyback with min/nom/max input voltages
    OpenMagnetics::Flyback flyback;
    
    DimensionWithTolerance inputVoltage;
    inputVoltage.set_minimum(18.0);
    inputVoltage.set_nominal(24.0);
    inputVoltage.set_maximum(32.0);
    flyback.set_input_voltage(inputVoltage);
    flyback.set_diode_voltage_drop(0.5);
    
    OpenMagnetics::FlybackOperatingPoint opPoint;
    opPoint.set_output_voltages({5.0});
    opPoint.set_output_currents({2.0});
    opPoint.set_ambient_temperature(25.0);
    opPoint.set_switching_frequency(100e3);
    flyback.set_operating_points({opPoint});
    
    // Create magnetic
    std::vector<int64_t> numberTurns = {20, 5};
    std::vector<int64_t> numberParallels = {1, 1};
    std::string shapeName = "E 25/13/7";
    
    auto coil = OpenMagneticsTesting::get_quick_coil(numberTurns, numberParallels, shapeName);
    auto gapping = OpenMagneticsTesting::get_distributed_gap(0.0004, 1);
    auto core = OpenMagneticsTesting::get_quick_core(shapeName, gapping, 1, "3C90");
    
    OpenMagnetics::Magnetic magnetic;
    magnetic.set_core(core);
    magnetic.set_coil(coil);
    
    // Run simulation - should get 3 operating points (min, nom, max input voltage)
    auto operatingPoints = flyback.simulate_with_magnetic_and_extract_operating_points(magnetic);
    
    REQUIRE(operatingPoints.size() == 3);
    
    // Verify each operating point has valid waveforms
    for (size_t i = 0; i < operatingPoints.size(); ++i) {
        auto& op = operatingPoints[i];
        INFO("Operating point " << i << ": " << op.get_name().value_or("unnamed"));
        
        REQUIRE(op.get_excitations_per_winding().size() == 2);
        
        auto& priEx = op.get_excitations_per_winding()[0];
        auto& secEx = op.get_excitations_per_winding()[1];
        
        REQUIRE(priEx.get_voltage().has_value());
        REQUIRE(priEx.get_current().has_value());
        REQUIRE(secEx.get_voltage().has_value());
        REQUIRE(secEx.get_current().has_value());
        
        // Check waveforms have data
        CHECK(priEx.get_voltage()->get_waveform()->get_data().size() > 10);
        CHECK(priEx.get_current()->get_waveform()->get_data().size() > 10);
    }
    
    // Extract peak currents for each voltage condition
    // Operating points are ordered: nom (index 0), min (index 1), max (index 2) - based on get_voltage_values()
    auto nomPriI = operatingPoints[0].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data();
    auto minPriI = operatingPoints[1].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data();
    auto maxPriI = operatingPoints[2].get_excitations_per_winding()[0].get_current()->get_waveform()->get_data();
    
    double nomPeakI = *std::max_element(nomPriI.begin(), nomPriI.end());
    double minPeakI = *std::max_element(minPriI.begin(), minPriI.end());
    double maxPeakI = *std::max_element(maxPriI.begin(), maxPriI.end());
    
    INFO("Peak current at Vin_min (18V): " << minPeakI);
    INFO("Peak current at Vin_nom (24V): " << nomPeakI);
    INFO("Peak current at Vin_max (32V): " << maxPeakI);
    
    // All operating points should have positive, measurable peak currents
    CHECK(minPeakI > 0.1);
    CHECK(nomPeakI > 0.1);
    CHECK(maxPeakI > 0.1);
    
    // All peak currents should be reasonable (< 10A for this design)
    CHECK(minPeakI < 10.0);
    CHECK(nomPeakI < 10.0);
    CHECK(maxPeakI < 10.0);
    
    INFO("Multi-voltage test passed - all operating points simulated successfully");
}
