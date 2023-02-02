#pragma once

#include <fstream>
#include <numbers>
#include <iostream>
#include <cmath>
#include <vector>
#include <filesystem>
#include <streambuf>
#include <nlohmann/json-schema.hpp>
#include <MAS.hpp>
#include <magic_enum.hpp>
#include "json.hpp"
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;


namespace OpenMagnetics {
    using nlohmann::json;

    class InputsWrapper:public Inputs {
        public:
        InputsWrapper(const json & j) {
            from_json(j, *this);
            auto check_passed = check_integrity();
            if (!check_passed.first) {
                throw std::runtime_error("Missing inputs");
            }
            process_waveforms();
        }
        InputsWrapper() = default;
        virtual ~InputsWrapper() = default;

        static double get_requirement_value(NumericRequirement requirement);
        std::pair<bool, std::string> check_integrity();
        void process_waveforms();
        static OperationPoint process_operation_point(OperationPoint operationPoint, double magnetizingInductance);

        static Waveform get_sampled_waveform(Waveform waveform, double frequency);
        static Processed get_processed_data(ElectromagneticParameter excitation, Waveform sampledWaveform, bool force);
        static Harmonics get_harmonics_data(Waveform waveform, double frequency);
        ElectromagneticParameter reflect_waveform(ElectromagneticParameter excitation, double ratio);
        static ElectromagneticParameter standarize_waveform(ElectromagneticParameter parameter, double frequency);
        OperationPoint get_operation_point(size_t index);
        OperationPointExcitation get_winding_excitation(size_t operationPointIndex, size_t windingIndex);
        OperationPointExcitation get_primary_excitation(size_t operationPointIndex);

        static void get_magnetizing_current(OperationPointExcitation& excitation, Waveform sampledWaveform, double magnetizingInductance);
        static ElectromagneticParameter add_offset_to_excitation(ElectromagneticParameter electromagneticParameter, double offset, double frequency);

        static InputsWrapper create_quick_operation_point(double frequency, double magnetizingInductance, double temperature, 
                                                WaveformLabel waveShape, double peakToPeak, double dutyCycle, double dcCurrent);

        void set_operation_point_by_index(const OperationPoint & value, size_t index) { get_mutable_operation_points()[index] = value; }


        void from_json(const json & j, Inputs & x);
        void to_json(json & j, const Inputs & x);

        inline void from_json(const json & j, InputsWrapper& x) {
            x.set_design_requirements(j.at("designRequirements").get<DesignRequirements>());
            x.set_operation_points(j.at("operationPoints").get<std::vector<OperationPoint>>());
        }
        inline void to_json(json & j, const InputsWrapper & x) {
            j = json::object();
            j["designRequirements"] = x.get_design_requirements();
            j["operationPoints"] = x.get_operation_points();
        }

    };
}
