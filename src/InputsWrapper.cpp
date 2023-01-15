#include <fstream>
#include <numbers>
#include <iostream>
#include <cmath>
#include <vector>
#include <filesystem>
#include <streambuf>
#include <limits>
#include <algorithm>
#include <cctype>
#include <string>
#include <complex>
#include <valarray>

#include <nlohmann/json-schema.hpp>
#include <magic_enum.hpp>
#include "json.hpp"
#include "Utils.h"
#include "InputsWrapper.h"
#include <libInterpolate/Interpolate.hpp>
#include "Constants.h"

   
using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;


namespace OpenMagnetics {

    // Cooley-Tukey FFT (in-place, breadth-first, decimation-in-frequency)
    // Better optimized but less intuitive
    // !!! Warning : in some cases this code make result different from not optimased version above (need to fix bug)
    // The bug is now fixed @2017/05/30 
    // from https://rosettacode.org/wiki/Fast_Fourier_transform#C++
    void fft(std::vector<std::complex<double>> &x)
    {
        // DFT
        size_t N = x.size(), k = N, n;
        double thetaT = std::numbers::pi / N;
        std::complex<double> phiT = std::complex<double>(cos(thetaT), -sin(thetaT)), T;
        while (k > 1)
        {
            n = k;
            k >>= 1;
            phiT = phiT * phiT;
            T = 1.0L;
            for (size_t l = 0; l < k; l++)
            {
                for (size_t a = l; a < N; a += n)
                {
                    size_t b = a + k;
                    std::complex<double> t = x[a] - x[b];
                    x[a] += x[b];
                    x[b] = t * T;
                }
                T *= phiT;
            }
        }
        // Decimate
        size_t m = (size_t)log2(N);
        for (size_t a = 0; a < N; a++)
        {
            size_t b = a;
            // Reverse bits
            b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
            b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
            b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
            b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
            b = ((b >> 16) | (b << 16)) >> (32 - m);
            if (b > a)
            {
                std::complex<double> t = x[a];
                x[a] = x[b];
                x[b] = t;
            }
        }
    }
    void print(std::vector<double> data){
        for (auto i: data) {
            std::cout << i << ' ';
        }
        std::cout << std::endl;
    }

    void print(double data){
        std::cout << data << std::endl;
    }
    void print(std::string data){
        std::cout << data << std::endl;
    }
    void print_json(json data){
        std::cout << data << std::endl;
    }

    std::vector<double> linear_spaced_array(double startingValue, double endingValue, std::size_t numberPoints) {
        double h = (endingValue - startingValue) / static_cast<double>(numberPoints);
        std::vector<double> xs(numberPoints);
        std::vector<double>::iterator x;
        double val;
        for (x = xs.begin(), val = startingValue; x != xs.end(); ++x, val += h) {
            *x = val;
        }
        return xs;
    }

    // In case the waveform comes defined with processed data only, we create a MAS format waveform from it, as the rest of the code depends on it.
    ElectromagneticParameter InputsWrapper::standarize_waveform(ElectromagneticParameter parameter, double frequency) {
        ElectromagneticParameter standardized_parameter(parameter);
        if (parameter.get_waveform() == nullptr) {
            Waveform waveform;
            double period = 1 / frequency;
            auto processed = *parameter.get_processed();

            auto peakToPeak = processed.get_peak_to_peak();
            auto offset = processed.get_offset();
            double dutyCycle;
            if (processed.get_duty_cycle() != nullptr) {
                dutyCycle = *processed.get_duty_cycle();
            }
            std::vector<double> data;
            std::vector<double> time;
            switch (processed.get_label()) {
                case OpenMagnetics::WaveformLabel::TRIANGULAR:
                    data = {-peakToPeak / 2 + offset, peakToPeak / 2 + offset, -peakToPeak / 2 + offset};
                    time = {0., dutyCycle * period, period};
                    break;
                case OpenMagnetics::WaveformLabel::SQUARE: {
                    double max = peakToPeak * (1 - dutyCycle);
                    double min = -peakToPeak * dutyCycle;
                    double dc = dutyCycle * period;
                    data = {max, max, min, min, max};
                    time = {0, dc, dc, period, period};
                    break;
                }
                case OpenMagnetics::WaveformLabel::SQUARE_WITH_DEAD_TIME: {
                    double max = +peakToPeak / 2;
                    double min = -peakToPeak / 2;
                    double dc = dutyCycle * period;
                    data = {0, 0, max, max, 0, 0, min, min, 0, 0};
                    time = {0,
                            0.25 * period - dc / 2,
                            0.25 * period - dc / 2,
                            0.25 * period + dc / 2,
                            0.25 * period + dc / 2,
                            0.75 * period - dc / 2,
                            0.75 * period - dc / 2,
                            0.75 * period + dc / 2,
                            0.75 * period + dc / 2,
                            1};
                    break;
                }
                case OpenMagnetics::WaveformLabel::SINUSOIDAL: {
                    auto constants = Constants();
                    for (size_t i=0; i < constants.number_points_samples_waveforms; ++i ) {
                        double angle = i * 2 * std::numbers::pi / constants.number_points_samples_waveforms;
                        time.push_back(angle * period);
                        data.push_back((sin(angle) * peakToPeak / 2) + offset);
                    }
                    break;
                }
                default:
                    break;
            }

            waveform.set_data(data);
            if (processed.get_label() == OpenMagnetics::WaveformLabel::SINUSOIDAL) {
                waveform.set_time(nullptr);
            }
            else {
                waveform.set_time(std::make_shared<std::vector<double>>(time));
            }
            standardized_parameter.set_waveform(std::make_shared<Waveform>(waveform));
        }

        return standardized_parameter;
    }

    Waveform InputsWrapper::get_sampled_waveform(Waveform waveform, double frequency){
        auto constants = Constants();

        int n = waveform.get_data().size();
        _1D::LinearInterpolator<double> interp;
        _1D::LinearInterpolator<double>::VectorType xx(n), yy(n);

        std::vector<double> time;
        auto data = waveform.get_data();
        if (waveform.get_time() == nullptr) { // This means the waveform is equidistant
            time = linear_spaced_array(0, 1. / frequency, data.size());
        }
        else {
            time = *waveform.get_time();
        }
        for(int i = 0; i < n; i++) {
            xx(i) = time[i];
            yy(i) = data[i];
        }

        auto sampledTime = linear_spaced_array(0, 1. / frequency, constants.number_points_samples_waveforms);

        interp.setData( xx, yy );
        std::vector<double> sampledData;

        for(int i = 0; i < constants.number_points_samples_waveforms; i++) {
            sampledData.push_back(interp(sampledTime[i]));
        }

        Waveform sampledWaveform;
        sampledWaveform.set_data(sampledData);
        sampledWaveform.set_time(std::make_shared<std::vector<double>>(sampledTime));
        return sampledWaveform;
    }

    double InputsWrapper::get_requirement_value(NumericRequirement requirement){
        if (requirement.get_nominal() == nullptr) {
            if (requirement.get_minimum() != nullptr && requirement.get_maximum() != nullptr) {
                return ((*requirement.get_minimum()) + (*requirement.get_maximum())) / 2;
            }
            else if (requirement.get_minimum() != nullptr) {
                return (*requirement.get_minimum());
            }
            else {
                return (*requirement.get_maximum());
            }
        }
        else {
            return (*requirement.get_nominal());
        }
    }

    Waveform derivate_waveform(Waveform sourceWaveform) {
        std::vector<double> source = sourceWaveform.get_data();
        std::vector<double> derivative;
        Waveform resultWaveform(sourceWaveform);

        source.push_back(source[0]);

        std::adjacent_difference(source.begin(), source.end(), source.begin());
        derivative = std::vector<double>(source.begin() + 1, source.end());
        resultWaveform.set_data(derivative);
        return resultWaveform;

    }

    Waveform get_magnetizing_current_waveform(Waveform sourceWaveform, double frequency, double magnetizing_inductance) {
        std::vector<double> source = sourceWaveform.get_data();
        std::vector<double> integration;
        Waveform resultWaveform(sourceWaveform);

        double integral = 0;
        integration.push_back(integral);
        double time_per_point = 1 / frequency / source.size();
        for (auto& point : source) {
            integral += point / magnetizing_inductance * time_per_point;
            integration.push_back(integral);
        }

        integration = std::vector<double>(integration.begin(), integration.end() - 1);

        double integrationAverage = std::accumulate(std::begin(integration), std::end(integration), 0.0) / integration.size();
        for (size_t i=0; i<integration.size(); ++i) {
            integration[i] = integration[i] - integrationAverage;
        }

        resultWaveform.set_data(integration);
        return resultWaveform;
    }

    std::shared_ptr<ElectromagneticParameter> InputsWrapper::reflect_waveform(ElectromagneticParameter primaryWaveform, double ratio) {
        ElectromagneticParameter reflected_waveform;
        Waveform waveform(*primaryWaveform.get_waveform());
        json primaryWaveformJson;
        to_json(primaryWaveformJson, primaryWaveform);
 
        waveform.set_data({});
        for (auto& datum : (*primaryWaveform.get_waveform()).get_data()) {
            waveform.get_mutable_data().push_back(datum * ratio);
        }
        // waveform.set_time((*primaryWaveform.get_waveform()).get_time());
        reflected_waveform.set_waveform(std::make_shared<Waveform>(waveform));


        return std::make_shared<ElectromagneticParameter>(reflected_waveform);
    }

    std::pair<bool, std::string> InputsWrapper::check_integrity() {
        auto operation_points = get_mutable_operation_points();
        auto turnsRatios = get_design_requirements().get_turns_ratios();
        std::pair<bool, std::string> result;
        result.first = true;
        result.second = "";

        for (auto &operation_point : operation_points){
            if (operation_point.get_excitations_per_winding().size() == 0){
                result.first = false;
                throw std::invalid_argument( "Missing excitation for primary" );
            }
        }

        for (size_t i=0; i < operation_points.size(); ++i ) {
            if (turnsRatios.size() > operation_points[i].get_excitations_per_winding().size() - 1) {
                if (turnsRatios.size() == 1 && operation_points[i].get_excitations_per_winding().size() == 1) {
                    // We are missing excitation only for secondary
                    for (size_t i=0; i < turnsRatios.size(); ++i ) {
                        if (i >= operation_points[i].get_excitations_per_winding().size() - 1) {
                            double turnRatio = get_requirement_value(turnsRatios[i]);
                            auto excitationOfPrimaryWinding = operation_points[i].get_excitations_per_winding()[0];
                            OperationPointExcitationPerWinding excitationOfThisWinding(excitationOfPrimaryWinding);

                            excitationOfThisWinding.get_mutable_excitation().set_voltage(reflect_waveform(*excitationOfPrimaryWinding.get_excitation().get_voltage(), 1 / turnRatio));
                            excitationOfThisWinding.get_mutable_excitation().set_current(reflect_waveform(*excitationOfPrimaryWinding.get_excitation().get_current(), turnRatio));
                            operation_points[i].get_mutable_excitations_per_winding().push_back(excitationOfThisWinding);
                        }
                    }
                    result.second = "Had to create the excitations of some windings based on primary";
                }
                else {
                    throw std::invalid_argument( "Missing excitation for more than one secondary. Only one can be guessed" );
                }
            }
        }

        set_operation_points(operation_points);

        for (auto &operation_point : get_mutable_operation_points()){
            json operation_pointJson;
            to_json(operation_pointJson, operation_point);
        }
        return result;
    }

    Processed InputsWrapper::get_processed_data(ElectromagneticParameter excitation, Waveform sampledWaveform) {
        Processed processed;
        std::vector<double> dataToProcess;
        auto harmonics = *excitation.get_harmonics();
        auto sampledDataToProcess = sampledWaveform.get_data();

        if (excitation.get_processed() != nullptr){
            processed = *excitation.get_processed();
        }
        else {
            auto waveform = *excitation.get_waveform();
            dataToProcess = waveform.get_data();

            std::string labelString;
            if (waveform.get_ancillary_label() != nullptr) {
                labelString = *waveform.get_ancillary_label();
            }
            else {
                labelString = "custom";
            }
            std::transform(labelString.begin(), labelString.end(), labelString.begin(), [](unsigned char c){ return std::toupper(c); });
            std::optional<WaveformLabel> label = magic_enum::enum_cast<WaveformLabel>(labelString);
            processed.set_label(label.value());

            processed.set_offset(std::accumulate(std::begin(sampledDataToProcess), std::end(sampledDataToProcess), 0.0) / sampledDataToProcess.size());
            processed.set_peak_to_peak(*max_element(dataToProcess.begin(), dataToProcess.end()) - *min_element(dataToProcess.begin(), dataToProcess.end()));

        }
        if (processed.get_duty_cycle() == nullptr) {
            // TODO
        }
        if (processed.get_effective_frequency() == nullptr) {
            double effectiveFrequency = 0;
            std::vector<double> dividend;
            std::vector<double> divisor;


            for (size_t i=0; i < harmonics.get_amplitudes().size(); ++i ) {
                double dataSquared = pow(harmonics.get_amplitudes()[i], 2);
                double timeSquared = pow(harmonics.get_frequencies()[i], 2);
                dividend.push_back(dataSquared * timeSquared);
                divisor.push_back(dataSquared);
            }
            double sumDivisor = std::reduce(divisor.begin(), divisor.end());
            if (sumDivisor > 0)
                effectiveFrequency = sqrt(std::reduce(dividend.begin(), dividend.end()) / sumDivisor);


            processed.set_effective_frequency(std::make_shared<double>(effectiveFrequency));
        }

        if (processed.get_rms() == nullptr) {
            double rms = 0.0;
            for (int i=0; i<int(sampledDataToProcess.size()); ++i) {
                rms += sampledDataToProcess[i]*sampledDataToProcess[i];
            }
            rms /= sampledDataToProcess.size();
            rms = sqrt(rms);
            processed.set_rms(std::make_shared<double>(rms));
        }
        if (processed.get_thd() == nullptr) {
            std::vector<double> dividend;
            double divisor = harmonics.get_amplitudes()[1];
            double thd = 0;
            for (size_t i=2; i < harmonics.get_amplitudes().size(); ++i ) {
                dividend.push_back(pow(harmonics.get_amplitudes()[i], 2));
            }


            if (divisor > 0)
                thd = sqrt(std::reduce(dividend.begin(), dividend.end())) / divisor;
            processed.set_thd(std::make_shared<double>(thd));
        }

        return processed;

    }
    Harmonics InputsWrapper::get_harmonics_data(Waveform waveform, double frequency) {
        auto dataToProcess = waveform.get_data();
        Harmonics harmonics;

        std::vector<std::complex<double>> data;
        for (std::size_t i = 0;  i < dataToProcess.size();  ++i) {
            data.emplace_back(dataToProcess[i]);
        }

        fft(data);

        harmonics.get_mutable_amplitudes().push_back(abs(data[0] / static_cast<double>(data.size())));
        for (size_t i = 1; i < data.size() / 2; ++i) {
            harmonics.get_mutable_amplitudes().push_back(abs(2. * data[i] / static_cast<double>(data.size())));
        }
        for (size_t i = 0; i < data.size() / 2; ++i) {
            harmonics.get_mutable_frequencies().push_back(frequency * i);
        }

        return harmonics;
    }

    void InputsWrapper::process_waveforms() {
        auto operation_points = get_mutable_operation_points();
        std::vector<OperationPoint> processed_operation_points;
        for (auto &operation_point : operation_points){
            std::vector<OperationPointExcitationPerWinding> processed_excitations_per_winding;
            for (auto &excitation_per_windind : operation_point.get_mutable_excitations_per_winding()){
                auto excitation = excitation_per_windind.get_excitation();
                // Here we processed this excitation current
                if (excitation.get_current() != nullptr) {
                    auto current_excitation = *excitation.get_current();
                    current_excitation = standarize_waveform(current_excitation, excitation.get_frequency());
                    auto waveform = *current_excitation.get_waveform();
                    auto sampledWaveform = get_sampled_waveform(waveform, excitation.get_frequency());
                    current_excitation.set_harmonics(std::make_shared<Harmonics>(get_harmonics_data(sampledWaveform, excitation.get_frequency())));
                    current_excitation.set_processed(std::make_shared<Processed>(get_processed_data(current_excitation, sampledWaveform)));
                    excitation.set_current(std::make_shared<ElectromagneticParameter>(current_excitation));
                }
                // Here we processed this excitation voltage
                if (excitation.get_voltage() != nullptr) {
                    auto voltage_excitation = *excitation.get_voltage();
                    voltage_excitation = standarize_waveform(voltage_excitation, excitation.get_frequency());
                    auto waveform = *voltage_excitation.get_waveform();
                    auto sampledWaveform = get_sampled_waveform(waveform, excitation.get_frequency());
                    voltage_excitation.set_harmonics(std::make_shared<Harmonics>(get_harmonics_data(sampledWaveform, excitation.get_frequency())));
                    voltage_excitation.set_processed(std::make_shared<Processed>(get_processed_data(voltage_excitation, sampledWaveform)));
                    excitation.set_voltage(std::make_shared<ElectromagneticParameter>(voltage_excitation));

                    if (excitation.get_magnetizing_current() == nullptr) {
                        auto sampledMagnetizingCurrentWaveform = get_magnetizing_current_waveform(sampledWaveform, excitation.get_frequency(), get_requirement_value(get_design_requirements().get_magnetizing_inductance()));
                        ElectromagneticParameter magnetizing_current_excitation;

                        magnetizing_current_excitation.set_waveform(std::make_shared<Waveform>(sampledMagnetizingCurrentWaveform));
                        magnetizing_current_excitation.set_harmonics(std::make_shared<Harmonics>(get_harmonics_data(sampledMagnetizingCurrentWaveform, excitation.get_frequency())));
                        magnetizing_current_excitation.set_processed(std::make_shared<Processed>(get_processed_data(magnetizing_current_excitation, sampledMagnetizingCurrentWaveform)));
                        excitation.set_magnetizing_current(std::make_shared<ElectromagneticParameter>(magnetizing_current_excitation));
                    }

                }
                excitation_per_windind.set_excitation(excitation);
                processed_excitations_per_winding.push_back(excitation_per_windind);
            }
            operation_point.set_excitations_per_winding(processed_excitations_per_winding);
            processed_operation_points.push_back(operation_point);
        }
        set_operation_points(processed_operation_points);

    }
}
