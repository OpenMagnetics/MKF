#include "converter_models/PushPull.h"
#include <iostream>

int main() {
    using namespace OpenMagnetics;
    
    // Test push-pull with fixed duty cycle
    json pushPullJson = {
        {"currentRippleRatio", 0.3},
        {"diodeVoltageDrop", 0.7},
        {"efficiency", 0.95},
        {"inputVoltage", {
            {"minimum", 20.0},
            {"maximum", 30.0}
        }},
        {"maximumSwitchCurrent", 10.0},
        {"desiredTurnsRatios", {2.67}},
        {"desiredInductance", 10e-6},
        {"desiredDutyCycle", {0.45, 0.45}},
        {"operatingPoints", {{
            {"outputVoltages", {48.0}},
            {"outputCurrents", {0.7}},
            {"switchingFrequency", 100000.0},
            {"ambientTemperature", 25.0}
        }}}
    };
    
    try {
        AdvancedPushPull pushPull(pushPullJson);
        auto inputs = pushPull.process();
        
        std::cout << "Success! Generated " << inputs.get_operating_points().size() << " operating points" << std::endl;
        
        for (size_t i = 0; i < inputs.get_operating_points().size(); i++) {
            auto& op = inputs.get_operating_points()[i];
            std::cout << "\nOperating Point " << i << ": " << op.get_name() << std::endl;
            std::cout << "  Excitations: " << op.get_excitations_per_winding().size() << std::endl;
        }
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}
