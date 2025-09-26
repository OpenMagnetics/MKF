#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {

class CurrentTransformer : public MAS::CurrentTransformer {
public:
    bool _assertErrors = false;

    CurrentTransformer(const json& j);
    CurrentTransformer() {};

    Inputs process(double turnsRatio, double secondaryDcResistance = 0);
    Inputs process(Magnetic magnetic);
    DesignRequirements process_design_requirements(double turnsRatio);
    DesignRequirements process_design_requirements(Magnetic magnetic);
    std::vector<OperatingPoint> process_operating_points(double turnsRatio, double secondaryDcResistance = 0);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);
};

} // namespace OpenMagnetics