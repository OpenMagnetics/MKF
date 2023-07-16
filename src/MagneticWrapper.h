#pragma once

#include <MAS.hpp>
#include "CoreWrapper.h"
#include "CoilWrapper.h"

namespace OpenMagnetics {

class MagneticWrapper : public Magnetic {
    private:
        CoreWrapper core;
        std::optional<std::vector<DistributorInfo>> distributors_info;
        std::optional<ManufacturerInfo> manufacturer_info;
        CoilWrapper winding;
    public:

        MagneticWrapper(const Magnetic magnetic) {
            set_core(magnetic.get_core());
            set_coil(magnetic.get_coil());

            if (magnetic.get_distributors_info()) {
                set_distributors_info(magnetic.get_distributors_info());
            }
            if (magnetic.get_manufacturer_info()) {
                set_manufacturer_info(magnetic.get_manufacturer_info());
            }
        }
        MagneticWrapper() = default;
        virtual ~MagneticWrapper() = default;
};
} // namespace OpenMagnetics
