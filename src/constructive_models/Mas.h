#pragma once

#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "processors/Outputs.h"

using namespace MAS;

namespace OpenMagnetics {

class Mas : public MAS::Mas {
    private:
        Inputs inputs;
        Magnetic magnetic;
        std::vector<Outputs> outputs;

    public:
        Mas() = default;
        virtual ~Mas() = default;

        public:
        /**
         * The description of the inputs that can be used to design a Magnetic
         */
        const Inputs & get_inputs() const { return inputs; }
        Inputs & get_mutable_inputs() { return inputs; }
        void set_inputs(const Inputs & value) { this->inputs = value; }

        /**
         * The description of a magnetic
         */
        const Magnetic & get_magnetic() const { return magnetic; }
        Magnetic & get_mutable_magnetic() { return magnetic; }
        void set_magnetic(const Magnetic & value) { this->magnetic = value; }

        /**
         * The description of the outputs that are produced after designing a Magnetic
         */
        const std::vector<Outputs> & get_outputs() const { return outputs; }
        std::vector<Outputs> & get_mutable_outputs() { return outputs; }
        void set_outputs(const std::vector<Outputs> & value) { this->outputs = value; }

        static Magnetic expand_magnetic(Magnetic magnetic);
        static Inputs expand_inputs(Magnetic magnetic, Inputs inputs);

};


void from_json(const json & j, Mas & x);
void to_json(json & j, const Mas & x);
void to_file(std::filesystem::path filepath, const Mas & x);

inline void from_json(const json & j, Mas& x) {
    x.set_inputs(j.at("inputs").get<Inputs>());
    x.set_magnetic(j.at("magnetic").get<Magnetic>());
    x.set_outputs(j.at("outputs").get<std::vector<Outputs>>());
}

inline void to_json(json & j, const Mas & x) {
    j = json::object();
    j["inputs"] = x.get_inputs();
    j["magnetic"] = x.get_magnetic();
    j["outputs"] = x.get_outputs();
}
inline void to_file(std::filesystem::path filepath, const Mas & x) {
    json masJson;
    to_json(masJson, x);

    std::ofstream myfile;
    myfile.open(filepath);
    myfile << masJson;
}
} // namespace OpenMagnetics
