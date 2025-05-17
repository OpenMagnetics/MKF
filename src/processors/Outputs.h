#pragma once

#include <MAS.hpp>

using namespace MAS;

namespace OpenMagnetics {

class Outputs : public MAS::Outputs {
  public:
    Outputs() = default;
    virtual ~Outputs() = default;
};
} // namespace OpenMagnetics
