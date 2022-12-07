#include <fstream>
#include <numbers>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <filesystem>
#include <streambuf>
#include <magic_enum.hpp>

#include "Constants.h"
#include <Core.h>

namespace OpenMagnetics {

    class Reluctance {
        private:
        public:
            virtual std::map<std::string, double> get_gap_reluctance(CoreGap gapInfo) = 0;
            virtual double get_core_reluctance() = 0;
            virtual double get_total_reluctance() = 0;
    };
}