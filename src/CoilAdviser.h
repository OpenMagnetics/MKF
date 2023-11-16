#pragma once
#include "WireAdviser.h"
#include "CoilWrapper.h"
#include "MasWrapper.h"
#include <MAS.hpp>


namespace OpenMagnetics {

class CoilAdviser : public WireAdviser {
    protected:
        int _interleavingLevel = 1;
    public:

        std::vector<std::pair<MasWrapper, double>> get_advised_coil(MasWrapper mas, size_t maximumNumberResults=1);
        std::vector<std::pair<MasWrapper, double>> get_advised_coil(std::vector<WireWrapper>* wires, MasWrapper mas, size_t maximumNumberResults=1);

        void set_interleaving_level(int interleavingLevel) {
            _interleavingLevel = interleavingLevel;
        }
};


} // namespace OpenMagnetics