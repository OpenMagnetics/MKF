#pragma once
#include "constructive_models/Wire.h"
#include "constructive_models/Coil.h"
#include "constructive_models/Magnetic.h"
#include "constructive_models/Mas.h"
#include "processors/Inputs.h"
#include "support/Utils.h"
#include <MAS.hpp>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace OpenMagnetics {

// Enumerates series/parallel combinations of the windings of a real, multi-winding
// transformer to synthesize "virtual" magnetics that hit a set of required turns
// ratios (and a magnetizing-inductance requirement). Each virtual winding is a group
// of one or more physical windings connected in series or parallel; the combinator
// records how each virtual winding maps back to the physical ones so that a solution
// found against a virtual magnetic can be "devirtualized" onto the real part.
//
// Ported from Wuerth PyMMM's MagneticCombinator. The pure combinator logic lives here;
// the json / by-reference / magnetic-cache entry points stay in the consumers.
class MagneticCombinator
{
private:
    size_t maximumNumberCombinationLevels = 3;

public:
    class WindingCombination
    {
        private:
            std::string name;
            std::optional<std::vector<MAS::ConnectionElement>> connections;
            size_t numberTurns;
            size_t numberParallels;
            std::vector<size_t> windingIndexes;
            bool isParallel = false;
            MAS::IsolationSide isolationSide;
            std::vector<WindingCombination> windingCombinations;
        public:
            const std::string & get_name() const { return name; }
            std::string & get_mutable_name() { return name; }
            void set_name(const std::string & value) { this->name = value; }

            std::optional<std::vector<MAS::ConnectionElement>> get_connections() const { return connections; }
            void set_connections(std::optional<std::vector<MAS::ConnectionElement>> value) { this->connections = value; }

            const size_t & get_number_turns() const { return numberTurns; }
            size_t & get_mutable_number_turns() { return numberTurns; }
            void set_number_turns(const size_t & value) { this->numberTurns = value; }

            const size_t & get_number_parallels() const { return numberParallels; }
            size_t & get_mutable_number_parallels() { return numberParallels; }
            void set_number_parallels(const size_t & value) { this->numberParallels = value; }

            const std::vector<size_t> & get_winding_indexes() const { return windingIndexes; }
            std::vector<size_t> & get_mutable_winding_indexes() { return windingIndexes; }
            void set_winding_indexes(const std::vector<size_t> & value) { this->windingIndexes = value; }

            const bool & get_is_parallel() const { return isParallel; }
            bool & get_mutable_is_parallel() { return isParallel; }
            void set_is_parallel(const bool & value) { this->isParallel = value; }

            const MAS::IsolationSide & get_isolation_side() const { return isolationSide; }
            MAS::IsolationSide & get_mutable_isolation_side() { return isolationSide; }
            void set_isolation_side(const MAS::IsolationSide & value) { this->isolationSide = value; }

            const std::vector<WindingCombination> & get_winding_combinations() const { return windingCombinations; }
            std::vector<WindingCombination> & get_mutable_winding_combinations() { return windingCombinations; }
            void set_winding_combinations(const std::vector<WindingCombination> & value) { this->windingCombinations = value; }
    };

    OpenMagnetics::Wire get_most_restrictive_wire(OpenMagnetics::Coil coil, std::vector<size_t> windingIndexes);
    std::vector<size_t> get_winding_indexes(OpenMagnetics::Magnetic magnetic);
    std::pair<std::map<std::string, WindingCombination>, std::vector<OpenMagnetics::Magnetic>> calculate_virtual_magnetics(OpenMagnetics::Magnetic magnetic, OpenMagnetics::Inputs inputs);
    static OpenMagnetics::Mas devirtualize_mas(OpenMagnetics::Magnetic originalMagnetic, OpenMagnetics::Mas virtualMas, std::map<std::string, WindingCombination> devirtualizationData);
    std::vector<WindingCombination> calculate_coil_combinations(OpenMagnetics::Coil coil, int numberDesiredWindings);
    std::vector<WindingCombination> calculate_combinations(std::vector<WindingCombination> sourceCombinations, int maximumWindings);

};

} // namespace OpenMagnetics
