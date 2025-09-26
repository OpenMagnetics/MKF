#include "converter_models/CurrentTransformer.h"
#include "physical_models/MagnetizingInductance.h"
#include "physical_models/WindingOhmicLosses.h"
#include "support/Utils.h"
#include <cfloat>

namespace OpenMagnetics {
    CurrentTransformer::CurrentTransformer(const json& j) {
        from_json(j, *this);
    }

    DesignRequirements CurrentTransformer::process_design_requirements(Magnetic magnetic) {
        if (magnetic.get_coil().get_functional_description().size() != 2) {
            throw std::runtime_error("A current transformer must have exactly two windings");
        }
        double turnsRatio = double(magnetic.get_coil().get_functional_description()[0].get_number_turns()) / magnetic.get_coil().get_functional_description()[1].get_number_turns();
        return process_design_requirements(turnsRatio);
    }

    DesignRequirements CurrentTransformer::process_design_requirements(double turnsRatio) {
        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();

        DimensionWithTolerance turnsRatioWithTolerance;
        turnsRatioWithTolerance.set_nominal(roundFloat(turnsRatio, 2));
        designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);

        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(roundFloat(1e-6, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;

        isolationSides.push_back(IsolationSide::PRIMARY);
        isolationSides.push_back(IsolationSide::SECONDARY);

        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::CURRENT_TRANSFORMER);
        return designRequirements;
    }

    std::vector<OperatingPoint> CurrentTransformer::process_operating_points(double turnsRatio, double secondaryDcResistance) {
        std::vector<OperatingPoint> operatingPoints;

        double peakToPeak;
        switch (get_waveform_label()) {
            case WaveformLabel::SINUSOIDAL:
                peakToPeak = get_maximum_primary_current_peak() * 2;
                break;
            case WaveformLabel::UNIPOLAR_RECTANGULAR:
            case WaveformLabel::UNIPOLAR_TRIANGULAR:
                peakToPeak = get_maximum_primary_current_peak();
                break;
            default:
                throw std::runtime_error("Only SINUSOIDAL, UNIPOLAR_RECTANGULAR, UNIPOLAR_TRIANGULAR are allowed for current transformers");
        }


        auto primaryCurrentWaveform = Inputs::create_waveform(get_waveform_label(), peakToPeak, get_frequency(), get_maximum_duty_cycle());
        SignalDescriptor primaryCurrent;
        primaryCurrent.set_waveform(primaryCurrentWaveform);
        primaryCurrent.set_harmonics(Inputs::calculate_harmonics_data(primaryCurrentWaveform, get_frequency()));
        primaryCurrent.set_processed(Inputs::calculate_processed_data(primaryCurrent, primaryCurrentWaveform, true, primaryCurrent.get_processed()));

        SignalDescriptor secondaryVoltage;
        auto secondaryVoltageWaveform = Inputs::multiply_waveform(primaryCurrentWaveform, get_burden_resistor());
        secondaryVoltageWaveform = Inputs::sum_waveform(secondaryVoltageWaveform, get_diode_voltage_drop() + secondaryDcResistance);
        secondaryVoltage.set_waveform(secondaryVoltageWaveform);
        secondaryVoltage.set_harmonics(Inputs::calculate_harmonics_data(secondaryVoltageWaveform, get_frequency()));
        secondaryVoltage.set_processed(Inputs::calculate_processed_data(secondaryVoltage, secondaryVoltageWaveform, true, secondaryVoltage.get_processed()));

        auto primaryVoltage = Inputs::reflect_waveform(secondaryVoltage, turnsRatio);
        primaryVoltage.set_harmonics(Inputs::calculate_harmonics_data(primaryVoltage.get_waveform().value(), get_frequency()));
        primaryVoltage.set_processed(Inputs::calculate_processed_data(primaryVoltage, primaryVoltage.get_waveform().value(), true, primaryVoltage.get_processed()));
        auto secondaryCurrent = Inputs::reflect_waveform(primaryCurrent, turnsRatio);
        secondaryCurrent.set_harmonics(Inputs::calculate_harmonics_data(secondaryCurrent.get_waveform().value(), get_frequency()));
        secondaryCurrent.set_processed(Inputs::calculate_processed_data(secondaryCurrent, secondaryCurrent.get_waveform().value(), true, secondaryCurrent.get_processed()));
        OperatingPointExcitation primaryExcitation;
        primaryExcitation.set_current(primaryCurrent);
        primaryExcitation.set_frequency(get_frequency());
        primaryExcitation.set_voltage(primaryVoltage);

        OperatingPointExcitation secondaryExcitation;
        secondaryExcitation.set_current(secondaryCurrent);
        secondaryExcitation.set_frequency(get_frequency());
        secondaryExcitation.set_voltage(secondaryVoltage);

        OperatingPoint operatingPoint;
        operatingPoint.set_excitations_per_winding({primaryExcitation, secondaryExcitation});

        return {operatingPoint};
    }

    std::vector<OperatingPoint> CurrentTransformer::process_operating_points(Magnetic magnetic) {
        if (magnetic.get_coil().get_functional_description().size() != 2) {
            throw std::runtime_error("A current transformer must have exactly two windings");
        }
        double turnsRatio = double(magnetic.get_coil().get_functional_description()[0].get_number_turns()) / magnetic.get_coil().get_functional_description()[1].get_number_turns();
        double secondaryDcResistance = WindingOhmicLosses::calculate_dc_resistance_per_winding(magnetic.get_coil(), get_ambient_temperature())[1];
        return process_operating_points(turnsRatio, secondaryDcResistance);
    }

    Inputs CurrentTransformer::process(double turnsRatio, double secondaryDcResistance) {

        Inputs inputs;
        auto designRequirements = process_design_requirements(turnsRatio);
        auto operatingPoints = process_operating_points(turnsRatio, secondaryDcResistance);

        inputs.set_design_requirements(designRequirements);
        inputs.set_operating_points(operatingPoints);

        return inputs;
    }

    Inputs CurrentTransformer::process(Magnetic magnetic) {
        if (magnetic.get_coil().get_functional_description().size() != 2) {
            throw std::runtime_error("A current transformer must have exactly two windings");
        }
        double turnsRatio = double(magnetic.get_coil().get_functional_description()[0].get_number_turns()) / magnetic.get_coil().get_functional_description()[1].get_number_turns();
        double secondaryDcResistance = WindingOhmicLosses::calculate_dc_resistance_per_winding(magnetic.get_coil(), get_ambient_temperature())[1];
        return process(turnsRatio, secondaryDcResistance);
    }
} // namespace OpenMagnetics
