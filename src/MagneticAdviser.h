#pragma once
#include "CoreAdviser.h"
#include "CoilAdviser.h"
#include "MasWrapper.h"
#include <MAS.hpp>


namespace OpenMagnetics {

class MagneticAdviser : public CoreAdviser, public CoilAdviser {
    private:
    public:

        std::vector<MasWrapper> get_advised_magnetic(InputsWrapper inputs, size_t maximumNumberResults=1);
        std::vector<MasWrapper> get_advised_magnetic(std::vector<CoreWrapper>* cores,
                                                     std::vector<WireWrapper>* wires,
                                                     InputsWrapper inputs,
                                                     size_t maximumNumberResults=1);

};


} // namespace OpenMagnetics