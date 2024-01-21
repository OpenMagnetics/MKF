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

        public:
        /**
         * The description of the inputs that can be used to design a Magnetic
         */
        const InputsWrapper & get_inputs() const { return inputs; }
        InputsWrapper & get_mutable_inputs() { return inputs; }
        void set_inputs(const InputsWrapper & value) { this->inputs = value; }

        /**
         * The description of a magnetic
         */
        const MagneticWrapper & get_magnetic() const { return magnetic; }
        MagneticWrapper & get_mutable_magnetic() { return magnetic; }
        void set_magnetic(const MagneticWrapper & value) { this->magnetic = value; }

        /**
         * The description of the outputs that are produced after designing a Magnetic
         */
        const std::vector<OutputsWrapper> & get_outputs() const { return outputs; }
        std::vector<OutputsWrapper> & get_mutable_outputs() { return outputs; }
        void set_outputs(const std::vector<OutputsWrapper> & value) { this->outputs = value; }
};


void from_json(const json & j, MasWrapper & x);
void to_json(json & j, const MasWrapper & x);

inline void from_json(const json & j, MasWrapper& x) {
    x.set_inputs(j.at("inputs").get<InputsWrapper>());
    x.set_magnetic(j.at("magnetic").get<MagneticWrapper>());
    x.set_outputs(j.at("outputs").get<std::vector<OutputsWrapper>>());
}

inline void to_json(json & j, const MasWrapper & x) {
    j = json::object();
    j["inputs"] = x.get_inputs();
    j["magnetic"] = x.get_magnetic();
    j["outputs"] = x.get_outputs();
}
} // namespace OpenMagnetics
