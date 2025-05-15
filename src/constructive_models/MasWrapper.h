#pragma once

#include <MAS.hpp>
#include "processors/InputsWrapper.h"
#include "constructive_models/MagneticWrapper.h"
#include "processors/OutputsWrapper.h"

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

        static MagneticWrapper expand_magnetic(MagneticWrapper magnetic);
        static InputsWrapper expand_inputs(MagneticWrapper magnetic, InputsWrapper inputs);

};


void from_json(const json & j, MasWrapper & x);
void to_json(json & j, const MasWrapper & x);
void to_file(std::filesystem::path filepath, const MasWrapper & x);

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
inline void to_file(std::filesystem::path filepath, const MasWrapper & x) {
    json masJson;
    to_json(masJson, x);

    std::ofstream myfile;
    myfile.open(filepath);
    myfile << masJson;
}
} // namespace OpenMagnetics
