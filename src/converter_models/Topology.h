#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include <string>
#include "constructive_models/Magnetic.h"


namespace OpenMagnetics {
using namespace MAS;


class Topology {
protected:
    std::string _magnetizingInductanceModel = "ZHANG";
    
public:
    bool _assertErrors = false;
    
    const std::string& get_magnetizing_inductance_model() const { return _magnetizingInductanceModel; }
    void set_magnetizing_inductance_model(const std::string& model) { _magnetizingInductanceModel = model; }

    Topology(const json& j);
    Topology() {
    };

    static OperatingPointExcitation complete_excitation(Waveform currentWaveform, Waveform voltageWaveform, double switchingFrequency, std::string name);
    virtual DesignRequirements process_design_requirements() = 0;
    virtual std::vector<OperatingPoint> process_operating_points(const std::vector<double>& turnsRatios, double magnetizingInductance) = 0;
    virtual bool run_checks(bool assert = false);
    Inputs process();

};

} // namespace OpenMagnetics