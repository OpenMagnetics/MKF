#pragma once
#include "CoilWrapper.h"
#include "MasWrapper.h"
#include <MAS.hpp>


namespace OpenMagnetics {

class CoilAdviser {
    public:

        std::vector<std::pair<MasWrapper, double>> get_advised_coil(MasWrapper mas, size_t maximumNumberResults=1);
};


} // namespace OpenMagnetics