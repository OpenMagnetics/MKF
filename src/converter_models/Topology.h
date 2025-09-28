#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {


class Topology {
public:
    bool _assertErrors = false;

    Topology(const json& j);
    Topology() {
    };

    static OperatingPointExcitation complete_excitation(Waveform currentWaveform, Waveform voltageWaveform, double switchingFrequency, std::string name);
    virtual DesignRequirements process_design_requirements() = 0;
    virtual std::vector<OperatingPoint> process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) = 0;
    virtual bool run_checks(bool assert = false);
    Inputs process();

};

} // namespace OpenMagnetics