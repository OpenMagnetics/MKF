#include "advisers/CoreAdviser.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"  // for is_energy_storing_topology()
#include "advisers/CoreAdviserInternal.h"      // log_stage/log_pruned, make_sinusoidal_excitation, CoreLossesModelPair
#include "constructive_models/Bobbin.h"
#include "physical_models/CoreTemperature.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/CoilMesher.h"
#include "constructive_models/NumberTurns.h"
#include "constructive_models/Wire.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/Bobbin.h"
#include <algorithm>
#include <cctype>
#include <set>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cfloat>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <execution>
#include <magic_enum_utility.hpp>
#include <list>
#include <cmrc/cmrc.hpp>
#include "support/Exceptions.h"
#include "support/Logger.h"
#include "support/Settings.h"

CMRC_DECLARE(data);


namespace OpenMagnetics {


void CoreAdviser::set_unique_core_shapes(bool value) {
    _uniqueCoreShapes = value;
}

bool CoreAdviser::get_unique_core_shapes() {
    return _uniqueCoreShapes;
}

void CoreAdviser::set_application(Application value) {
    _application = value;
}

Application CoreAdviser::get_application() {
    return _application;
}

void CoreAdviser::set_mode(CoreAdviser::CoreAdviserModes value) {
    _mode = value;
}

CoreAdviser::CoreAdviserModes CoreAdviser::get_mode() {
    return _mode;
}

void CoreAdviser::set_weights(std::map<CoreAdviserFilters, double> weights) {
    _weights = weights;
}

std::map<std::string, std::map<CoreAdviser::CoreAdviserFilters, double>> CoreAdviser::get_scorings(bool weighted){
    std::map<std::string, std::map<CoreAdviser::CoreAdviserFilters, double>> swappedScorings;
    for (auto& [filter, aux] : _scorings) {
        auto filterConfiguration = _filterConfiguration[filter];

        auto normalizedScorings = OpenMagnetics::normalize_scoring(aux, weighted? _weights[filter] : 1, false, false);

        for (auto& [name, scoring] : aux) {
            swappedScorings[name][filter] = normalizedScorings[name];
        }
    }
    return swappedScorings;
}

std::vector<double> normalize_scoring(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring, std::vector<double> newScoring, double weight, std::map<std::string, bool> filterConfiguration) {
    auto normalizedScorings = OpenMagnetics::normalize_scoring(newScoring, weight, filterConfiguration);

    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i) {
        (*magneticsWithScoring)[i].second += normalizedScorings[i];
    }

    // stable_sort((*magneticsWithScoring).begin(), (*magneticsWithScoring).end(), [](std::pair<Magnetic, double>& b1, std::pair<Magnetic, double>& b2) {
    //     return b1.second > b2.second;
    // });

    return normalizedScorings;
}

void sort_magnetics_by_scoring(std::vector<std::pair<Magnetic, double>>* magneticsWithScoring) {
    stable_sort((*magneticsWithScoring).begin(), (*magneticsWithScoring).end(), [](const std::pair<Magnetic, double>& b1, const std::pair<Magnetic, double>& b2) {
        return b1.second > b2.second;
    }); // F12 FIX: stable_sort for reproducible results
}



Coil get_dummy_coil(Inputs inputs) {
    double temperature = inputs.get_maximum_temperature();
    double frequency = 0;
    for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex) {
        frequency = std::max(frequency, Inputs::get_primary_excitation(inputs.get_operating_point(operatingPointIndex)).get_frequency());
    }

    // Set round wire with diameter to two times the skin depth 
    auto wire = Wire::get_wire_for_frequency(frequency, temperature, true);
    Winding primaryWinding;
    primaryWinding.set_isolation_side(IsolationSide::PRIMARY);
    primaryWinding.set_name("primary");
    primaryWinding.set_number_parallels(1);
    primaryWinding.set_number_turns(1);
    primaryWinding.set_wire(wire);

    Coil coil;
    coil.set_bobbin(DUMMY_SENTINEL_NAME);
    coil.set_functional_description({primaryWinding});
    return coil;
}

static CoreShape scale_core_shape(const CoreShape& ref, double scaleFactor) {
    CoreShape scaled = ref;
    if (!ref.get_dimensions()) return scaled;
    auto dims = ref.get_dimensions().value();
    for (auto& [key, dim] : dims) {
        if (std::holds_alternative<DimensionWithTolerance>(dim)) {
            auto& dwt = std::get<DimensionWithTolerance>(dim);
            if (dwt.get_nominal()) dwt.set_nominal(dwt.get_nominal().value() * scaleFactor);
            if (dwt.get_minimum()) dwt.set_minimum(dwt.get_minimum().value() * scaleFactor);
            if (dwt.get_maximum()) dwt.set_maximum(dwt.get_maximum().value() * scaleFactor);
        } else {
            dim = std::get<double>(dim) * scaleFactor;
        }
    }
    scaled.set_dimensions(dims);
    scaled.set_name("Custom " + ref.get_name().value_or("unknown"));
    return scaled;
}

std::vector<CoreShape> CoreAdviser::create_custom_core_shapes(Inputs inputs) {
    MagneticFilterAreaProduct apFilter(inputs);
    double apRequired = apFilter.get_estimated_area_product_required(inputs);

    if (coreShapeDatabase.empty()) {
        load_core_shapes();
    }

    std::vector<CoreShapeFamily> gappableFamilies = {
        CoreShapeFamily::E, CoreShapeFamily::ETD, CoreShapeFamily::EQ,
        CoreShapeFamily::RM, CoreShapeFamily::PQ, CoreShapeFamily::EP,
        CoreShapeFamily::EFD, CoreShapeFamily::EL, CoreShapeFamily::ER,
        CoreShapeFamily::EC, CoreShapeFamily::LP, CoreShapeFamily::EPX,
        CoreShapeFamily::U, CoreShapeFamily::T,
        CoreShapeFamily::PLANAR_E, CoreShapeFamily::PLANAR_EL, CoreShapeFamily::PLANAR_ER
    };

    std::vector<CoreShape> customShapes;
    for (auto family : gappableFamilies) {
        try {
            CoreShape refShape = find_core_shape_by_area_product(apRequired, family);
            Core refCore(refShape);
            refCore.process_data();
            auto columns = refCore.get_columns();
            if (columns.empty()) continue;
            double apRef = columns[0].get_area() * refCore.get_winding_window().get_area().value();
            if (apRef <= 0) continue;

            double s = std::pow(apRequired / apRef, 0.25);
            s = std::clamp(s, 0.5, 3.0);
            CoreShape scaled = scale_core_shape(refShape, s);
            customShapes.push_back(scaled);
        } catch (...) { continue; }
    }
    return customShapes;
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, size_t maximumNumberResults) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, size_t maximumNumberResults){
    if (get_mode() == CoreAdviserModes::AVAILABLE_CORES) {
        if (coreDatabase.empty()) {
            load_cores();
        }
        return get_advised_core(inputs, weights, &coreDatabase, maximumNumberResults);
    }
    else if (get_mode() == CoreAdviserModes::STANDARD_CORES) {
        if (coreMaterialDatabase.empty()) {
            load_core_materials();
        }
        if (coreShapeDatabase.empty()) {
            load_core_shapes();
        }
        std::vector<MAS::CoreShape> shapes;
        // Phase 1 fix: coreShapeDatabase is keyed by canonical-name AND by
        // each alias (see Utils.cpp:255-260), so a naive iteration yields the
        // same CoreShape (1 + #aliases) times. That was the root cause of
        // CoreAdviser STANDARD_CORES × POWER returning duplicate top-N
        // entries with identical names and scores (e.g. three copies of
        // "95 PQ 27/15 gapped 0.36 mm" in slots 0-2). Dedupe by canonical
        // shape.get_name() before feeding the dataset builder.
        std::set<std::string> seenShapeNames;
        for (auto& [key, shape] : coreShapeDatabase) {
            auto canonical = shape.get_name().value_or("");
            if (canonical.empty()) {
                throw std::runtime_error("CoreShape in coreShapeDatabase has no name");
            }
            if (!seenShapeNames.insert(canonical).second) {
                continue;
            }
            shapes.push_back(shape);
        }
        return get_advised_core(inputs, &shapes, maximumNumberResults);
    }
    else if (get_mode() == CoreAdviserModes::CUSTOM_CORES) {
        if (coreMaterialDatabase.empty()) {
            load_core_materials();
        }
        if (coreShapeDatabase.empty()) {
            load_core_shapes();
        }
        auto customShapes = create_custom_core_shapes(inputs);
        return get_advised_core(inputs, &customShapes, maximumNumberResults);
    }
    else {
        throw NotImplementedException("Unknown CoreAdviserMode");
    }
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, cores, maximumNumberResults);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::vector<Core>* cores, size_t maximumNumberResults, size_t maximumNumberCores) {
    std::map<CoreAdviserFilters, double> weights;
    magic_enum::enum_for_each<CoreAdviserFilters>([&] (auto val) {
        CoreAdviserFilters filter = val;
        weights[filter] = 1.0;
    });
    return get_advised_core(inputs, weights, cores, maximumNumberResults, maximumNumberCores);
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults, size_t maximumNumberCores) {

    std::vector<std::pair<Mas, double>> results;


    for(size_t i = 0; i < (*cores).size(); i += maximumNumberCores) {
        auto last = std::min((*cores).size(), i + maximumNumberCores);
        std::vector<Core> partialCores;
        partialCores.reserve(last - i);

        move((*cores).begin() + i, (*cores).begin() + last, std::back_inserter(partialCores));
        
        auto partialResult = get_advised_core(inputs, weights, &partialCores, maximumNumberResults);
        std::move(partialResult.begin(), partialResult.end(), std::back_inserter(results));
    }

    stable_sort(results.begin(), results.end(), [](const std::pair<Mas, double>& b1, const std::pair<Mas, double>& b2) {
        return b1.second > b2.second;
    });

    if (results.size() > maximumNumberResults) {
        results = std::vector<std::pair<Mas, double>>(results.begin(), results.end() - (results.size() - maximumNumberResults));
    }
    return results;
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::map<CoreAdviserFilters, double> weights, std::vector<Core>* cores, size_t maximumNumberResults){
    _weights = weights;
 
    size_t maximumMagneticsAfterFiltering = settings.get_core_adviser_maximum_magnetics_after_filtering();
    std::vector<std::pair<Magnetic, double>> magnetics;

    magnetics = create_magnetic_dataset(inputs, cores, false);

    // logEntry("We start the search with " + std::to_string(magnetics.size()) + " magnetics for the first filter, culling to " + std::to_string(maximumMagneticsAfterFiltering) + " for the remaining filters.", "CoreAdviser");
    // logEntry("We don't include stacks of cores in our search.", "CoreAdviser");
    std::vector<std::pair<Mas, double>> filteredMagnetics;
    
    if (get_application() == Application::POWER) {
        filteredMagnetics = filter_available_cores_power_application(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        if (filteredMagnetics.size() >= maximumNumberResults) {
            return filteredMagnetics;
        }

        auto globalIncludeStacks = settings.get_core_adviser_include_stacks();
        if (globalIncludeStacks) {
            expand_magnetic_dataset_with_stacks(inputs, cores, &magnetics);
        }

        logEntry("First attempt produced not enough results, so now we are searching again with " + std::to_string(magnetics.size()) + " magnetics, including up to " + std::to_string(defaults.coreAdviserMaximumNumberStacks) + " cores stacked when possible.", "CoreAdviser");
        maximumMagneticsAfterFiltering = magnetics.size();
        filteredMagnetics = filter_available_cores_power_application(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        return filteredMagnetics;
    }
    else {
        filteredMagnetics = filter_available_cores_suppression_application(&magnetics, inputs, weights, maximumMagneticsAfterFiltering, maximumNumberResults);
        return filteredMagnetics;
    }
}

std::vector<std::pair<Mas, double>> CoreAdviser::get_advised_core(Inputs inputs, std::vector<CoreShape>* shapes, size_t maximumNumberResults) {
    auto globalIncludeStacks = settings.get_core_adviser_include_stacks();
    auto magnetics = create_magnetic_dataset(inputs, shapes, globalIncludeStacks);

    size_t maximumMagneticsAfterFiltering = settings.get_core_adviser_maximum_magnetics_after_filtering();
    if (get_application() == Application::POWER) {
        return filter_standard_cores_power_application(&magnetics, inputs, _weights, maximumMagneticsAfterFiltering, maximumNumberResults);
    }
    else {
        return filter_standard_cores_interference_suppression_application(&magnetics, inputs, _weights, maximumMagneticsAfterFiltering, maximumNumberResults);
    }
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::create_magnetic_dataset(Inputs inputs, std::vector<Core>* cores, bool includeStacks) {
    std::vector<std::pair<Magnetic, double>> magnetics;
    Coil coil = get_dummy_coil(inputs);
    auto includeToroidalCores = settings.get_use_toroidal_cores();
    auto includeConcentricCores = settings.get_use_concentric_cores();
    auto globalIncludeStacks = settings.get_core_adviser_include_stacks();
    auto globalIncludeDistributedGaps = settings.get_core_adviser_include_distributed_gaps();
    double maximumHeight = std::numeric_limits<double>::infinity();
    if (inputs.get_design_requirements().get_maximum_dimensions()) {
        if (inputs.get_design_requirements().get_maximum_dimensions()->get_height()) {
            maximumHeight = inputs.get_design_requirements().get_maximum_dimensions()->get_height().value();
        }
    }

    Magnetic magnetic;

    magnetic.set_coil(std::move(coil));
    for (auto& core : *cores){
        auto coreMaterial = core.resolve_material();
        if (!Core::check_material_application(coreMaterial, get_application())) {
            continue;
        }

        if (get_application() == Application::INTERFERENCE_SUPPRESSION) {
            // Common-mode chokes MUST be wound on a single closed magnetic path so the
            // common-mode flux of both windings adds in the core: only toroidal cores
            // give the tight coupling required. Differential-mode chokes have no such
            // constraint — they're regular inductors and can use any core type
            // (powdered-iron toroids, gapped ferrite E/EE/PQ, etc.). Apply the
            // toroidal-only restriction only to CMC; for DMC and other suppression
            // topologies, honour the user's includeToroidalCores / includeConcentricCores
            // settings.
            auto topology = inputs.get_design_requirements().get_topology();
            bool isCommonModeChoke = topology.has_value() && topology.value() == Topologies::COMMON_MODE_CHOKE;
            if (isCommonModeChoke) {
                if (core.get_type() != CoreType::TOROIDAL) {
                    continue;
                }
            }
            else {
                if (!includeToroidalCores && core.get_type() == CoreType::TOROIDAL) {
                    continue;
                }
                if (!includeConcentricCores && (core.get_type() == CoreType::PIECE_AND_PLATE || core.get_type() == CoreType::TWO_PIECE_SET)) {
                    continue;
                }
            }
        }
        else {
            if (!includeToroidalCores && core.get_type() == CoreType::TOROIDAL) {
                continue;
            }

            if (!includeConcentricCores && (core.get_type() == CoreType::PIECE_AND_PLATE || core.get_type() == CoreType::TWO_PIECE_SET)) {
                continue;
            }
        }
        core.process_data();

        if (inputs.get_wiring_technology() == WiringTechnology::PRINTED) {
            if (core.get_type() == CoreType::TOROIDAL) {
                continue;
            }
            if (core.get_columns().size() == 2) {
                continue;
            }
            auto windingWindow = core.get_winding_window();
            if (windingWindow.get_height().value() > windingWindow.get_width().value()) {
                continue;
            }
        }

        if (!core.process_gap()) {
            continue;
        }

        if (core.get_type() == CoreType::TWO_PIECE_SET) {
            if (core.get_height() > maximumHeight) {
                continue;
            }
        }

        if (!globalIncludeDistributedGaps && core.get_gapping().size() > core.get_processed_description()->get_columns().size()) {
            continue;
        }

        if (includeStacks && globalIncludeStacks && (core.get_shape_family() == CoreShapeFamily::E || core.get_shape_family() == CoreShapeFamily::PLANAR_E || core.get_shape_family() == CoreShapeFamily::T || core.get_shape_family() == CoreShapeFamily::U || core.get_shape_family() == CoreShapeFamily::C)) {
            for (size_t i = 0; i < defaults.coreAdviserMaximumNumberStacks; ++i)
            {
                core.get_mutable_functional_description().set_number_stacks(1 + i);
                // process_data() resets processed description to base values, then calls scale_to_stacks internally
                core.process_data();
                core.process_gap(); // CA-OPT-2 FIX: reprocess gap data after stacking (was commented out)
                magnetic.set_core(core);
                MagneticManufacturerInfo magneticmanufacturerinfo;
                if (i != 0) {
                    magneticmanufacturerinfo.set_reference(core.get_name().value_or("unnamed") + " " + std::to_string(1 + i) + " stacks" );
                }
                else{
                    magneticmanufacturerinfo.set_reference(core.get_name().value_or("unnamed"));
                }
                magnetic.set_manufacturer_info(magneticmanufacturerinfo);
                magnetics.push_back(std::pair<Magnetic, double>{magnetic, 0});
            }
        }
        else {
            magnetic.set_core(core);
            MagneticManufacturerInfo magneticmanufacturerinfo;
            magneticmanufacturerinfo.set_reference(core.get_name().value_or("unnamed"));
            magnetic.set_manufacturer_info(magneticmanufacturerinfo);
            magnetics.push_back(std::pair<Magnetic, double>{magnetic, 0});
        }
    }

    return magnetics;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::create_magnetic_dataset(Inputs inputs, std::vector<CoreShape>* shapes, bool includeStacks) {
    std::vector<std::pair<Magnetic, double>> magnetics;
    Coil coil = get_dummy_coil(inputs);
    auto includeToroidalCores = settings.get_use_toroidal_cores();
    auto includeConcentricCores = settings.get_use_concentric_cores();
    auto globalIncludeStacks = settings.get_core_adviser_include_stacks();
    auto globalIncludeDistributedGaps = settings.get_core_adviser_include_distributed_gaps();
    double maximumHeight = std::numeric_limits<double>::infinity();
    if (inputs.get_design_requirements().get_maximum_dimensions()) {
        if (inputs.get_design_requirements().get_maximum_dimensions()->get_height()) {
            maximumHeight = inputs.get_design_requirements().get_maximum_dimensions()->get_height().value();
        }
    }

    Magnetic magnetic;

    magnetic.set_coil(std::move(coil));
    for (auto& shape : *shapes){
        if (shape.get_family() == CoreShapeFamily::PQI || shape.get_family() == CoreShapeFamily::UI || shape.get_family() == CoreShapeFamily::UT || shape.get_family() == CoreShapeFamily::C) {
            continue;
        }
        Core core(shape);

        if (get_application() == Application::INTERFERENCE_SUPPRESSION) {
            // See comment in create_magnetic_dataset(): only CMC requires toroidal-only;
            // DMC and other suppression topologies should honour user settings.
            auto topology = inputs.get_design_requirements().get_topology();
            bool isCommonModeChoke = topology.has_value() && topology.value() == Topologies::COMMON_MODE_CHOKE;
            if (isCommonModeChoke) {
                if (core.get_type() != CoreType::TOROIDAL) {
                    continue;
                }
            }
            else {
                if (!includeToroidalCores && core.get_type() == CoreType::TOROIDAL) {
                    continue;
                }
                if (!includeConcentricCores && (core.get_type() == CoreType::PIECE_AND_PLATE || core.get_type() == CoreType::TWO_PIECE_SET)) {
                    continue;
                }
            }
        }
        else {
            if (!includeToroidalCores && core.get_type() == CoreType::TOROIDAL) {
                continue;
            }

            if (!includeConcentricCores && (core.get_type() == CoreType::PIECE_AND_PLATE || core.get_type() == CoreType::TWO_PIECE_SET)) {
                continue;
            }
        }
        core.process_data();

        if (inputs.get_wiring_technology() == WiringTechnology::PRINTED) {
            if (core.get_type() == CoreType::TOROIDAL) {
                continue;
            }
            if (core.get_columns().size() == 2) {
                continue;
            }
            auto windingWindow = core.get_winding_window();
            if (windingWindow.get_height().value() > windingWindow.get_width().value()) {
                continue;
            }
        }


        if (!core.process_gap()) {
            continue;
        }

        if (core.get_type() == CoreType::TWO_PIECE_SET) {
            if (core.get_height() > maximumHeight) {
                continue;
            }
        }

        if (!globalIncludeDistributedGaps && core.get_gapping().size() > core.get_processed_description()->get_columns().size()) {
            continue;
        }

        if (includeStacks && globalIncludeStacks && (core.get_shape_family() == CoreShapeFamily::E || core.get_shape_family() == CoreShapeFamily::PLANAR_E || core.get_shape_family() == CoreShapeFamily::T || core.get_shape_family() == CoreShapeFamily::U || core.get_shape_family() == CoreShapeFamily::C)) {
            // Capture the un-stacked base name once so stack variants don't
            // accumulate suffixes across loop iterations.
            std::string baseName = core.get_name().value_or("unnamed");
            for (size_t i = 0; i < defaults.coreAdviserMaximumNumberStacks; ++i) {
                core.get_mutable_functional_description().set_number_stacks(1 + i);
                // process_data() resets processed description to base values, then calls scale_to_stacks internally
                core.process_data();
                core.process_gap(); // CA-OPT-2 FIX: reprocess gap data after stacking (was commented out)
                // Phase 1 fix: previously only manufacturer_info.reference got
                // the " N stacks" suffix; core.name stayed identical across
                // stack variants. add_gapping_standard_cores then produced
                // identical display names (e.g. three copies of
                // "95 PQ 27/15 gapped 0.36 mm") for 1/2/3-stack variants.
                // Put the stack suffix in the core name too so downstream
                // names stay unique.
                std::string variantName = baseName;
                if (i != 0) {
                    variantName += " " + std::to_string(1 + i) + " stacks";
                }
                core.set_name(variantName);
                magnetic.set_core(core);
                MagneticManufacturerInfo magneticManufacturerInfo;
                magneticManufacturerInfo.set_reference(variantName);
                magnetic.set_manufacturer_info(magneticManufacturerInfo);
                magnetics.push_back(std::pair<Magnetic, double>{magnetic, 0});
            }
        }
        else {
            magnetic.set_core(core);
            MagneticManufacturerInfo magneticManufacturerInfo;
            magneticManufacturerInfo.set_reference(core.get_name().value_or("unnamed"));
            magnetic.set_manufacturer_info(magneticManufacturerInfo);
            magnetics.push_back(std::pair<Magnetic, double>{magnetic, 0});
        }
    }

    return magnetics;
}

void CoreAdviser::expand_magnetic_dataset_with_stacks(Inputs inputs, std::vector<Core>* cores, std::vector<std::pair<Magnetic, double>>* magnetics) {
    Coil coil = get_dummy_coil(inputs);
    auto includeToroidalCores = settings.get_use_toroidal_cores();
    double maximumHeight = std::numeric_limits<double>::infinity();
    if (inputs.get_design_requirements().get_maximum_dimensions()) {
        if (inputs.get_design_requirements().get_maximum_dimensions()->get_height()) {
            maximumHeight = inputs.get_design_requirements().get_maximum_dimensions()->get_height().value();
        }
    }

    Magnetic magnetic;

    magnetic.set_coil(std::move(coil));
    for (auto& core : *cores){
        if (!includeToroidalCores && core.get_type() == CoreType::TOROIDAL) {
            continue;
        }

        if (core.get_type() == CoreType::TWO_PIECE_SET) {
            if (core.get_height() > maximumHeight) {
                continue;
            }
        }

        if (core.get_shape_family() == CoreShapeFamily::E || core.get_shape_family() == CoreShapeFamily::PLANAR_E || core.get_shape_family() == CoreShapeFamily::T || core.get_shape_family() == CoreShapeFamily::U || core.get_shape_family() == CoreShapeFamily::C) {

            core.process_data();
            if (!core.process_gap()) {
                continue;
            }

            for (size_t i = 1; i < defaults.coreAdviserMaximumNumberStacks; ++i)
            {
                core.get_mutable_functional_description().set_number_stacks(1 + i);
                // process_data() resets processed description to base values, then calls scale_to_stacks internally
                core.process_data();
                core.process_gap(); // CA-OPT-2 FIX: reprocess gap data after stacking (was commented out)
                MagneticManufacturerInfo magneticManufacturerInfo;
                if (i!=0) {
                    std::string name = core.get_name().value_or("unnamed");
                    name = std::regex_replace(std::string(name), std::regex(" [0-9] stacks"), "");
                    core.set_name(name + " " + std::to_string(1 + i) + " stacks" );
                    magneticManufacturerInfo.set_reference(core.get_name().value_or("unnamed"));
                }
                magnetic.set_manufacturer_info(magneticManufacturerInfo);
                magnetic.set_core(core);
                (*magnetics).push_back(std::pair<Magnetic, double>{magnetic, 0});
            }
        }

    }
}

void add_initial_turns_by_inductance(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    MagnetizingInductance magnetizingInductance;
    
    // Transformer vs Inductor/Energy-Storing Detection
    // =================================================
    // We need to determine if this application stores energy (needs gap) or transfers it (no gap).
    //
    // ENERGY-STORING (inductor/flyback-like):
    //   - Topologies: Flyback, Buck, Boost, Buck-Boost, SEPIC, Cuk, Zeta, etc.
    //   - Energy is stored in the magnetic field each switching cycle
    //   - Core typically needs a gap to store the required energy without saturating
    //   - Turns calculated from inductance requirement and gap
    //
    // TRANSFORMER (forward-like):
    //   - Topologies: Forward, PushPull, Half-Bridge, Full-Bridge, DAB, LLC, etc.
    //   - Energy is transferred directly, not stored in the magnetic field
    //   - Magnetizing inductance should be HIGH to minimize magnetizing current (parasitic)
    //   - Core should NOT be gapped (high permeability = high inductance)
    //   - Turns calculated from voltage/frequency to avoid saturation (Faraday's Law)
    //
    // Detection priority:
    // 1. Use topology if specified (most reliable)
    // 2. Fall back to inductance field heuristic (minimum-only = transformer)
    //
    auto topology = inputs.get_design_requirements().get_topology();
    bool isEnergyStoring = is_energy_storing_topology(topology);
    
    // If topology is not set, use the old heuristic based on inductance specification
    // minimum-only inductance suggests transformer (want high inductance, no specific value)
    // nominal or min+max suggests inductor (need specific value for energy storage)
    bool isTransformer;
    if (topology.has_value()) {
        isTransformer = !isEnergyStoring;
    } else {
        // Legacy heuristic: minimum-only = transformer
        isTransformer = inputs.get_design_requirements().get_magnetizing_inductance().get_minimum() &&
                        !inputs.get_design_requirements().get_magnetizing_inductance().get_nominal() &&
                        !inputs.get_design_requirements().get_magnetizing_inductance().get_maximum();
    }
    
    // Pre-compute shared transformer values once (not per-core)
    double transformerTemperature = 25.0;
    double maxVoltSeconds = 0.0;
    if (isTransformer) {
        for (size_t opIdx = 0; opIdx < inputs.get_operating_points().size(); ++opIdx) {
            auto op = inputs.get_operating_point(opIdx);
            transformerTemperature = std::max(transformerTemperature, op.get_conditions().get_ambient_temperature());
            auto excitation = Inputs::get_primary_excitation(op);

            if (excitation.get_voltage() && excitation.get_voltage()->get_waveform() &&
                excitation.get_voltage()->get_waveform()->get_time()) {
                auto voltageWaveform = excitation.get_voltage()->get_waveform().value();
                const auto& data = voltageWaveform.get_data();
                auto time = voltageWaveform.get_time().value();
                double integral = 0.0;
                for (size_t j = 0; j + 1 < std::min(data.size(), time.size()); ++j) {
                    integral += data[j] * (time[j + 1] - time[j]);
                    maxVoltSeconds = std::max(maxVoltSeconds, std::abs(integral));
                }
            } else if (excitation.get_voltage() && excitation.get_voltage()->get_processed() &&
                       excitation.get_voltage()->get_processed()->get_peak()) {
                double frequency = excitation.get_frequency() > 0 ? excitation.get_frequency() : 100000;
                double omega = 2 * std::numbers::pi * frequency;
                maxVoltSeconds = std::max(maxVoltSeconds,
                    excitation.get_voltage()->get_processed()->get_peak().value() / omega);
            }
        }
    }

    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i){

        Core core = (*magneticsWithScoring)[i].first.get_core();
        if (!core.get_processed_description()) {
            core.process_data();
            core.process_gap();
        }
        double initialNumberTurns = (*magneticsWithScoring)[i].first.get_coil().get_functional_description()[0].get_number_turns();

        if (initialNumberTurns == 1) {
            if (isTransformer) {
                // Feasibility seed: fewest turns that keep B_peak under the material's safe
                // operating flux. The loss filter refines N upward toward the loss optimum.
                double bMax = core.get_magnetic_flux_density_saturation(transformerTemperature, true);
                double effectiveArea = core.get_processed_description()->get_effective_parameters().get_effective_area();

                if (maxVoltSeconds > 0 && effectiveArea > 0 && bMax > 0) {
                    initialNumberTurns = std::max(1.0, std::ceil(maxVoltSeconds / (bMax * effectiveArea)));
                } else {
                    initialNumberTurns = 5;
                }
            } else {
                // For inductors/flybacks, calculate from inductance requirement
                initialNumberTurns = magnetizingInductance.calculate_number_turns_from_gapping_and_inductance(core, &inputs, DimensionalValues::MINIMUM);
            }
        }
        if (inputs.get_design_requirements().get_turns_ratios().size() > 0) {
            NumberTurns numberTurns(initialNumberTurns, inputs.get_design_requirements());
            auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            initialNumberTurns = numberTurnsCombination[0];
        }

        (*magneticsWithScoring)[i].first.get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(initialNumberTurns);
    }
}

std::vector<std::pair<Magnetic, double>> add_initial_turns_by_impedance(std::vector<std::pair<Magnetic, double>> magneticsWithScoring, Inputs inputs) {
    Impedance impedance;
    std::vector<std::pair<Magnetic, double>> magneticsWithScoringAndTurns;
    for (size_t i = 0; i < magneticsWithScoring.size(); ++i){
        auto [magnetic, scoring] = magneticsWithScoring[i];
        Core core = magnetic.get_core();
        if (!core.get_processed_description()) {
            core.process_data();
            core.process_gap();
        }
        Bobbin bobbin;
        if (inputs.get_wiring_technology() == WiringTechnology::PRINTED) {
            bobbin = Bobbin::create_quick_bobbin(core, true);
        }
        else {
            bobbin = Bobbin::create_quick_bobbin(core);
        }
        magnetic.get_mutable_coil().set_bobbin(bobbin);

        double initialNumberTurns = magnetic.get_coil().get_functional_description()[0].get_number_turns();

        try {
            initialNumberTurns = impedance.calculate_minimum_number_turns(magnetic, inputs);
            if (initialNumberTurns < 1) {
                continue;
            }
        }
        catch (...) { // XC-3 NOTE: consider catching std::exception for diagnostics
            continue;
        }
        if (inputs.get_design_requirements().get_turns_ratios().size() > 0) {
            NumberTurns numberTurns(initialNumberTurns, inputs.get_design_requirements());
            auto numberTurnsCombination = numberTurns.get_next_number_turns_combination();
            initialNumberTurns = numberTurnsCombination[0];
        }

        magnetic.get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(initialNumberTurns);
        magneticsWithScoringAndTurns.push_back({magnetic, scoring});
    }

    return magneticsWithScoringAndTurns;
}

void add_alternative_materials(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    CoreMaterialCrossReferencer coreMaterialCrossReferencer(std::map<std::string, std::string>{{"coreLosses", "Steinmetz"}});
    double temperature = inputs.get_maximum_temperature();
    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i){
        Core core = (*magneticsWithScoring)[i].first.get_core();
        auto coreMaterial = core.resolve_material();

        // Phase 1: per-candidate alternatives lookup is annotational metadata.
        // If the cross-reference can't proceed (e.g. this candidate's material
        // has no Steinmetz coefficients at this temperature, so it can't be
        // used as a reference for ranking others), record an empty
        // alternatives list and continue. Named, logged, bounded scope.
        std::vector<std::string> coreMaterialAlternatives;
        try {
            auto alternatives = coreMaterialCrossReferencer.get_cross_referenced_core_material(coreMaterial, temperature);
            for (auto [alternativeCoreMaterial, scoring] : alternatives) {
                coreMaterialAlternatives.push_back(alternativeCoreMaterial.get_name());
            }
        }
        catch (const std::exception& e) {
            logEntry(std::string("Skipping alternative-materials lookup for candidate with material '")
                         + coreMaterial.get_name() + "': " + e.what(),
                     "CoreAdviser", 2);
        }
        coreMaterial.set_alternatives(coreMaterialAlternatives);
        core.set_material(coreMaterial);

        (*magneticsWithScoring)[i].first.set_core(core);
    }
}



}