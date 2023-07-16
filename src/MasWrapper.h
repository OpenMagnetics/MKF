#pragma once

#include <MAS.hpp>
#include "InputsWrapper.h"
#include "MagneticWrapper.h"
#include "OutputsWrapper.h"

namespace OpenMagnetics {

class MasWrapper : public Mas {
    private:
        InputsWrapper inputs;
        MagneticWrapper magnetic;
        std::vector<OutputsWrapper> outputs;

    public:
        MasWrapper() = default;
        virtual ~MasWrapper() = default;
};
} // namespace OpenMagnetics
