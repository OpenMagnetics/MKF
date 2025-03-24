#include <MAS.hpp>
#include "MasWrapper.h"
#include "MagnetizingInductance.h"
#include "Utils.h"

namespace OpenMagnetics {

MagneticWrapper MasWrapper::expand_magnetic(MagneticWrapper magnetic) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();
    auto coreMaterial = core.resolve_material();
    core.get_mutable_functional_description().set_material(coreMaterial);

    if (!core.get_processed_description()) {
        core.process_data();
    }

    if (core.get_functional_description().get_gapping().size() > 0 && !core.get_functional_description().get_gapping()[0].get_area()) {
        core.process_gap();
    }

    BobbinWrapper bobbin;

    if (std::holds_alternative<std::string>(coil.get_bobbin())) {
        auto bobbinName = std::get<std::string>(coil.get_bobbin());
        if (bobbinName == "Basic") {
            bobbin = BobbinWrapper::create_quick_bobbin(core, false);
            coil.set_bobbin(bobbin);
        }
        else if (bobbinName == "Dummy") {
            bobbin = BobbinWrapper::create_quick_bobbin(core, true);
            coil.set_bobbin(bobbin);
        }
    }
    else {
        bobbin = coil.resolve_bobbin();
    }


    if (!coil.get_turns_description()) {
        for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); windingIndex++) {
            coil.resolve_wire(windingIndex);
        }
        for (size_t windingIndex = 0; windingIndex < coil.get_functional_description().size(); windingIndex++)
        {
            auto wire = coil.get_wires()[windingIndex];
            if (wire.get_type() == WireType::FOIL) {
                if (!wire.get_conducting_height())
                {
                    auto bobbin = coil.resolve_bobbin();
                    DimensionWithTolerance aux;
                    aux.set_nominal(bobbin.get_processed_description()->get_winding_windows()[0].get_height().value() * 0.8);
                    wire.set_conducting_height(aux);
                }
                if (!wire.get_outer_height())
                {
                    wire.set_outer_height(wire.get_conducting_height().value());
                }
                if (!wire.get_outer_width())
                {
                    wire.set_outer_width(wire.get_conducting_width().value());
                }
            }
            if (wire.get_type() == WireType::RECTANGULAR) {
                if (!wire.get_outer_height())
                {
                    DimensionWithTolerance aux;
                    aux.set_nominal(WireWrapper::get_outer_height_rectangular(resolve_dimensional_values(wire.get_conducting_height().value())));
                    wire.set_outer_height(aux);
                }
                if (!wire.get_outer_width())
                {
                    DimensionWithTolerance aux;
                    aux.set_nominal(WireWrapper::get_outer_height_rectangular(resolve_dimensional_values(wire.get_conducting_width().value())));
                    wire.set_outer_width(aux);
                }
            }
            if (wire.get_type() == WireType::ROUND) {
                if (!wire.get_outer_diameter())
                {
                    auto coating = wire.resolve_coating();
                    if (coating->get_type() == InsulationWireCoatingType::ENAMELLED)
                    {
                        DimensionWithTolerance aux;
                        aux.set_nominal(WireWrapper::get_outer_diameter_round(resolve_dimensional_values(resolve_dimensional_values(wire.get_conducting_diameter().value()))));
                        wire.set_outer_diameter(aux);
                    }
                    
                    if (coating->get_type() == InsulationWireCoatingType::INSULATED)
                    {
                        int numberLayers = coating->get_number_layers().value();
                        int thicknessLayers = coating->get_thickness_layers().value();
                        DimensionWithTolerance aux;
                        aux.set_nominal(WireWrapper::get_outer_diameter_round(resolve_dimensional_values(wire.get_conducting_diameter().value()), numberLayers, thicknessLayers));
                        wire.set_outer_diameter(aux);
                    }
                }
            }
            if (wire.get_type() == WireType::LITZ) {
                if (!wire.get_outer_diameter())
                {
                    auto coating = wire.resolve_coating();
                    auto strand = wire.resolve_strand();
                    if (coating->get_type() == InsulationWireCoatingType::SERVED)
                    {
                        DimensionWithTolerance aux;
                        aux.set_nominal(WireWrapper::get_outer_diameter_served_litz(resolve_dimensional_values(strand.get_conducting_diameter()), wire.get_number_conductors().value()));
                        wire.set_outer_diameter(aux);
                    }
                    
                    if (coating->get_type() == InsulationWireCoatingType::INSULATED)
                    {
                        int numberLayers = coating->get_number_layers().value();
                        int thicknessLayers = coating->get_thickness_layers().value();
                        DimensionWithTolerance aux;
                        aux.set_nominal(WireWrapper::get_outer_diameter_insulated_litz(resolve_dimensional_values(strand.get_conducting_diameter()), wire.get_number_conductors().value(), numberLayers, thicknessLayers));
                    }
                }
            }
            coil.get_mutable_functional_description()[windingIndex].set_wire(wire);
        }
        if (!coil.get_sections_description())
        {
            coil.wind();
        }
        else {
            if (!coil.get_layers_description())
            {
                coil.wind_by_layers();
            }
            if (!coil.get_turns_description())
            {
                coil.wind_by_turns();
                coil.delimit_and_compact();
            }
        }
    }

    if (coil.get_layers_description()) {
        auto layers = coil.get_layers_description().value();
        for (size_t layerIndex = 0; layerIndex < layers.size(); layerIndex++) {
            auto insulationMaterial = coil.resolve_insulation_layer_insulation_material(layers[layerIndex]);
            layers[layerIndex].set_insulation_material(insulationMaterial);
        }
        coil.set_layers_description(layers);
    }

    magnetic.set_core(core);
    magnetic.set_coil(coil);

    return magnetic;
}

InputsWrapper MasWrapper::expand_inputs(MagneticWrapper magnetic, InputsWrapper inputs) {
    auto core = magnetic.get_core();
    auto coil = magnetic.get_coil();

    auto numberWindings = coil.get_functional_description().size();
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); operatingPointIndex++) {
        for (size_t excitationIndex = 0; excitationIndex < inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding().size(); excitationIndex++) {
            auto frequency = inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[excitationIndex].get_frequency();
            auto operatingPoint = inputs.get_operating_points()[operatingPointIndex];
            auto excitation = inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding()[excitationIndex];

            if (inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[excitationIndex].get_current()) {
                auto current = inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[excitationIndex].get_current().value();
                if (current.get_waveform()) {
                    if (current.get_processed()) {
                        auto waveform = InputsWrapper::create_waveform(current.get_processed().value(), frequency);
                        current.set_waveform(InputsWrapper::create_waveform(current.get_processed().value(), frequency));
                    }
                    else if (current.get_harmonics()) {
                        auto harmonics = current.get_harmonics().value();
                        auto waveform = InputsWrapper::reconstruct_signal(harmonics, frequency);
                        auto processed = InputsWrapper::calculate_processed_data(harmonics, waveform, true);
                        current.set_processed(processed);
                        current.set_waveform(waveform);
                    }
                }

                auto waveform = current.get_waveform().value();
                if (!current.get_harmonics()) {
                    auto sampledCurrentWaveform = InputsWrapper::calculate_sampled_waveform(waveform, frequency);
                    auto harmonics = InputsWrapper::calculate_harmonics_data(sampledCurrentWaveform, frequency);
                    current.set_harmonics(harmonics);
                }

                auto harmonics = current.get_harmonics().value();
                auto processed = InputsWrapper::calculate_processed_data(harmonics, waveform, true);
                current.set_processed(processed);
                inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[excitationIndex].set_current(current);
            }

            if (inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[excitationIndex].get_current()) {
                auto voltage = inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[excitationIndex].get_voltage().value();
                if (voltage.get_waveform()) {
                    if (voltage.get_processed()) {
                        auto waveform = InputsWrapper::create_waveform(voltage.get_processed().value(), frequency);
                        voltage.set_waveform(InputsWrapper::create_waveform(voltage.get_processed().value(), frequency));
                    }
                    else if (voltage.get_harmonics()) {
                        auto harmonics = voltage.get_harmonics().value();
                        auto waveform = InputsWrapper::reconstruct_signal(harmonics, frequency);
                        auto processed = InputsWrapper::calculate_processed_data(harmonics, waveform, true);
                        voltage.set_processed(processed);
                        voltage.set_waveform(waveform);
                    }
                }

                auto waveform = voltage.get_waveform().value();
                if (!voltage.get_harmonics()) {
                    auto sampledvoltageWaveform = InputsWrapper::calculate_sampled_waveform(waveform, frequency);
                    auto harmonics = InputsWrapper::calculate_harmonics_data(sampledvoltageWaveform, frequency);
                    voltage.set_harmonics(harmonics);
                }

                auto harmonics = voltage.get_harmonics().value();
                auto processed = InputsWrapper::calculate_processed_data(harmonics, waveform, true);
                voltage.set_processed(processed);
                inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[excitationIndex].set_voltage(voltage);
            }

            MagnetizingInductance magnetizingInductanceObj;
            double magnetizingInductance = magnetizingInductanceObj.calculate_inductance_from_number_turns_and_gapping(core, coil, &operatingPoint).get_magnetizing_inductance().get_nominal().value();

            auto magnetizingCurrent = InputsWrapper::calculate_magnetizing_current(excitation, magnetizingInductance, true, 0.0);

            // if (excitation.get_voltage()) {
            //     if (excitation.get_voltage().value().get_processed()) {
            //         if (excitation.get_voltage().value().get_processed().value().get_duty_cycle()) {
            //             auto processed = current.get_processed().value();
            //             processed.set_duty_cycle(excitation.get_voltage().value().get_processed().value().get_duty_cycle().value());
            //             current.set_processed(processed);
            //         }
            //     }
            // }

            if (excitationIndex == 0 && numberWindings == 2 && inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding().size() == 1) {

                auto primaryExcitation = inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding()[0];
                auto turnRatio = coil.get_number_turns(0) / coil.get_number_turns(1);

                OperatingPointExcitation secondaryExcitation(primaryExcitation);
                auto currentSignalDescriptorProcessed = InputsWrapper::calculate_basic_processed_data(primaryExcitation.get_current().value().get_waveform().value());
                auto voltageSignalDescriptorProcessed = InputsWrapper::calculate_basic_processed_data(primaryExcitation.get_voltage().value().get_waveform().value());

                auto voltageSignalDescriptor = InputsWrapper::reflect_waveform(primaryExcitation.get_voltage().value(), 1.0 / turnRatio, voltageSignalDescriptorProcessed.get_label());
                auto currentSignalDescriptor = InputsWrapper::reflect_waveform(primaryExcitation.get_current().value(), turnRatio, currentSignalDescriptorProcessed.get_label());

                auto voltageSampledWaveform = InputsWrapper::calculate_sampled_waveform(voltageSignalDescriptor.get_waveform().value(), frequency);
                voltageSignalDescriptor.set_harmonics(InputsWrapper::calculate_harmonics_data(voltageSampledWaveform, frequency));
                voltageSignalDescriptor.set_processed(InputsWrapper::calculate_processed_data(voltageSignalDescriptor, voltageSampledWaveform, true));

                auto currentSampledWaveform = InputsWrapper::calculate_sampled_waveform(currentSignalDescriptor.get_waveform().value(), frequency);
                currentSignalDescriptor.set_harmonics(InputsWrapper::calculate_harmonics_data(currentSampledWaveform, frequency));
                currentSignalDescriptor.set_processed(InputsWrapper::calculate_processed_data(currentSignalDescriptor, currentSampledWaveform, true));

                secondaryExcitation.set_voltage(voltageSignalDescriptor);
                secondaryExcitation.set_current(currentSignalDescriptor);

                inputs.get_mutable_operating_points()[operatingPointIndex].get_mutable_excitations_per_winding().push_back(secondaryExcitation);
            }
        }
    }

    return inputs;
}

} // namespace OpenMagnetics
