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
        virtual ~InputsWrapper() = default;

        double get_requirement_value(NumericRequirement requirement);
        std::pair<bool, std::string> check_integrity();
        void process_waveforms();
        Waveform get_sampled_waveform(Waveform waveform, double frequency);
        Processed get_processed_data(ElectromagneticParameter excitation, Waveform sampledWaveform);
        Harmonics get_harmonics_data(Waveform waveform, double frequency);
        std::shared_ptr<ElectromagneticParameter> reflect_waveform(ElectromagneticParameter excitation, double ratio);
        ElectromagneticParameter standarize_waveform(ElectromagneticParameter parameter, double frequency);
    };
}
