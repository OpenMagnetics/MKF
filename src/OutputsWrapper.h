#pragma once

#include <MAS.hpp>

namespace OpenMagnetics {

class OutputsWrapper : public Outputs {
  public:
    OutputsWrapper() = default;
    virtual ~OutputsWrapper() = default;
};
} // namespace OpenMagnetics
