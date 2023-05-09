#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <nlohmann/json-schema.hpp>
#include <numbers>
#include <streambuf>
#include <vector>
#include "Constants.h"
#include "Utils.h"
#include "WindingWrapper.h"
#include "json.hpp"

#include <magic_enum.hpp>
#include "../tests/TestingUtils.h"

using nlohmann::json_uri;
using nlohmann::json_schema::json_validator;
using json = nlohmann::json;

namespace OpenMagnetics {
std::vector<bool> WindingWrapper::wind_by_consecutive_turns(std::vector<WindingFunctionalDescription> windings) {
    // std::cout << windings.size() << std::endl;
    std::vector<bool> windByConsecutiveTurns;
    for (auto & windingFunctionalDescription : windings) {
        uint64_t numberTurns = windingFunctionalDescription.get_number_turns();
        uint64_t numberParallels = windingFunctionalDescription.get_number_parallels();
        if (numberTurns == _interleavingLevel) {
            windByConsecutiveTurns.push_back(false);
            continue;
        }
        if (numberParallels == _interleavingLevel) {
            windByConsecutiveTurns.push_back(true);
            continue;
        }
        if (numberParallels % _interleavingLevel == 0) {
            windByConsecutiveTurns.push_back(true);
            continue;
        }
        if (numberTurns % _interleavingLevel == 0) {
            windByConsecutiveTurns.push_back(false);
            continue;
        }
        windByConsecutiveTurns.push_back(true);
    }
    return windByConsecutiveTurns;
}



void WindingWrapper::wind_by_sections() {
    std::vector<Section> sectionsDescription;
    // std::cout <<get_functional_description().size() << std::endl;
    auto windByConsecutiveTurns = wind_by_consecutive_turns(get_functional_description());
    auto windingWindows = std::get<Bobbin>(get_bobbin()).get_processed_description().value().get_winding_windows();
    if (windingWindows.size() > 1) {
        throw std::runtime_error("Bobbins with more than winding window not implemented yet");
    }
    double windingWindowHeight = windingWindows[0].get_height().value();
    double windingWindowWidth = windingWindows[0].get_width().value();
    std::vector<std::vector<double>> remainingParallelsProportion;
    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        remainingParallelsProportion.push_back(std::vector<double>(get_number_parallels(windingIndex), 1));
    }
    double interleavedHeight = 0;
    double interleavedWidth = 0;
    double currentSectionCenterWidth = 0;
    double currentSectionCenterHeight = 0;
    if (_windingOrientation == WindingOrientation::HORIZONTAL) {
        interleavedWidth = roundFloat(windingWindowWidth / _interleavingLevel / get_functional_description().size(), 9);
        interleavedHeight = windingWindowHeight;
        currentSectionCenterWidth = windingWindows[0].get_coordinates().value()[0] - windingWindowWidth / 2 + interleavedWidth / 2;
        currentSectionCenterHeight = windingWindows[0].get_coordinates().value()[1];
    }
    else if (_windingOrientation == WindingOrientation::VERTICAL) {
        interleavedWidth = windingWindowWidth;
        interleavedHeight = roundFloat(windingWindowHeight / _interleavingLevel / get_functional_description().size(), 9);
        currentSectionCenterWidth = windingWindows[0].get_coordinates().value()[0];
        currentSectionCenterHeight = windingWindows[0].get_coordinates().value()[1] + windingWindowHeight / 2 - interleavedHeight / 2;
    }
    else {
        throw std::runtime_error("Toroidal windings not implemented");
    }


    for (size_t sectionIndex = 0; sectionIndex < _interleavingLevel; ++sectionIndex) {

        for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
            PartialWinding partialWinding;
            Section section;
            std::vector<double> sectionParallelsProportion(get_number_parallels(windingIndex), 0);
            partialWinding.set_winding(get_name(windingIndex));
            if (windByConsecutiveTurns[windingIndex]) {
                double totalPhysicalTurns = get_number_turns(windingIndex) * get_number_parallels(windingIndex);

                size_t remainingPhysicalTurns = 0;
                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    remainingPhysicalTurns += round(remainingParallelsProportion[windingIndex][parallelIndex] * get_number_turns(windingIndex));
                }
                uint64_t maxPhysicalTurnsThisSection = std::min(remainingPhysicalTurns, uint64_t(ceil(totalPhysicalTurns / _interleavingLevel)));

                size_t currentParallel = 0;
                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    if (remainingParallelsProportion[windingIndex][parallelIndex] > 0) {
                        currentParallel = parallelIndex;
                        break;
                    }
                }

                while (maxPhysicalTurnsThisSection > 0) {
                    uint64_t numberTurnsToFittingInCurrentParallel = round(remainingParallelsProportion[windingIndex][currentParallel] * get_number_turns(windingIndex));
                    if (maxPhysicalTurnsThisSection >= numberTurnsToFittingInCurrentParallel) {
                        maxPhysicalTurnsThisSection -= numberTurnsToFittingInCurrentParallel;
                        sectionParallelsProportion[currentParallel] = double(numberTurnsToFittingInCurrentParallel) / get_number_turns(windingIndex);
                        currentParallel++;
                    }
                    else {
                        double proportionParallelsThisSection = double(maxPhysicalTurnsThisSection) / get_number_turns(windingIndex);
                        sectionParallelsProportion[currentParallel] += proportionParallelsThisSection;
                        maxPhysicalTurnsThisSection = 0;
                    }
                }
            }
            else {
                double numberTurnsToAddToCurrentParallel = ceil(double(get_number_turns(windingIndex)) / _interleavingLevel);
                double proportionParallelsThisSection = std::min(remainingParallelsProportion[windingIndex][0], numberTurnsToAddToCurrentParallel / get_number_turns(windingIndex));

                for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                    sectionParallelsProportion[parallelIndex] = proportionParallelsThisSection;
                }
            }
            partialWinding.set_parallels_proportion(sectionParallelsProportion);
            section.set_name(get_name(windingIndex) +  " section " + std::to_string(sectionIndex));
            section.set_partial_windings(std::vector<PartialWinding>{partialWinding});  // TODO: support more than one winding per section?
            section.set_layers_orientation(_layersOrientation);
            section.set_dimensions(std::vector<double>{interleavedWidth, interleavedHeight});
            section.set_coordinates(std::vector<double>{currentSectionCenterWidth, currentSectionCenterHeight, 0});
            sectionsDescription.push_back(section);

            for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
                remainingParallelsProportion[windingIndex][parallelIndex] -= sectionParallelsProportion[parallelIndex];
            }
            if (_windingOrientation == WindingOrientation::HORIZONTAL) {
                currentSectionCenterWidth += interleavedWidth;
            }
            else if (_windingOrientation == WindingOrientation::VERTICAL) {
                currentSectionCenterHeight -= interleavedHeight;
            }
            else {
                throw std::runtime_error("Toroidal windings not implemented");
            }
        }

    }

    for (size_t windingIndex = 0; windingIndex < get_functional_description().size(); ++windingIndex) {
        for (size_t parallelIndex = 0; parallelIndex < get_number_parallels(windingIndex); ++parallelIndex) {
            if (roundFloat(remainingParallelsProportion[windingIndex][parallelIndex], 9) != 0) {
                throw std::runtime_error("There are unassigned parallel proportion, something went wrong");
            }
        }
    }

    set_sections_description(sectionsDescription);

}

void WindingWrapper::wind_by_layers() {}

void WindingWrapper::wind_by_turns() {}

} // namespace OpenMagnetics
