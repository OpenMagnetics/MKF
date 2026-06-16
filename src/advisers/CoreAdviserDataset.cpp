// Phase 5 (structural cleanup, dataset construction): extracted from
// CoreAdviser.cpp. This TU owns the helpers that build (and post-process)
// the pool of candidate Magnetics passed to the filter pipeline:
//
//   - get_dummy_coil           (single-winding stand-in for un-wound flows)
//   - scale_core_shape         (used by create_custom_core_shapes)
//   - CoreAdviser::create_custom_core_shapes
//   - CoreAdviser::create_magnetic_dataset (both overloads, from Cores and CoreShapes)
//   - CoreAdviser::expand_magnetic_dataset_with_stacks
//   - add_initial_turns_by_inductance / add_initial_turns_by_impedance
//   - add_alternative_materials
//
// Sibling TU of CoreAdviser.cpp; declarations live in advisers/CoreAdviser.h
// and shared helpers in advisers/CoreAdviserInternal.h.

#include "advisers/CoreAdviser.h"
#include "advisers/CoreAdviserInternal.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "advisers/MagneticFilter.h"
#include "advisers/MagneticFilterInternal.h"  // for is_energy_storing_topology()
#include "constructive_models/Bobbin.h"
#include "constructive_models/Insulation.h"
#include "constructive_models/NumberTurns.h"
#include "constructive_models/Wire.h"
#include "physical_models/ComplexPermeability.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/CoilMesher.h"
#include "support/Exceptions.h"
#include "support/Logger.h"
#include "support/Settings.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <execution>
#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <magic_enum_utility.hpp>

namespace OpenMagnetics {

Coil get_dummy_coil(const Inputs& inputs) {
    double temperature = inputs.get_maximum_temperature();
    double frequency = 0;
    // Phase 6 (perf): cache operating-points by const-ref.
    const auto& operatingPoints = inputs.get_operating_points();
    for (size_t operatingPointIndex = 0; operatingPointIndex < operatingPoints.size(); ++operatingPointIndex) {
        frequency = std::max(frequency, Inputs::get_primary_excitation(operatingPoints[operatingPointIndex]).get_frequency());
    }

    // Set round wire with diameter to two times the skin depth 
    auto wire = Wire::get_wire_for_frequency(frequency, temperature, true);

    // Determine the dummy-coil winding count from the operating point.
    //
    // Default = 1 winding (the historical behaviour). Char tests for
    // flyback/forward/push-pull/etc. were recorded against this shape;
    // widening the dummy coil for those topologies pushes new core/
    // material combos through MagnetizingInductance + CoreLosses and
    // exposes NaN-handling bugs (visible as "Too large losses" thrown
    // out of MagneticAdviser).
    //
    // Exception: for CMCs the single-winding dummy is wrong — Magnetizing
    // Inductance routes through the single-winding path, treats the full
    // line current as magnetising current, and rejects every candidate
    // for false saturation. CMCs need an N-winding dummy so the
    // can_be_common_mode_choke path is taken and the differential
    // magnetising current (tiny) is what feeds saturation.
    //
    // Use Inputs::can_be_common_mode_choke as the discriminator — it's
    // the same check MagnetizingInductance uses downstream, so the
    // two stay in sync by construction.
    size_t numberOfWindings = 1;
    if (!operatingPoints.empty() &&
        Inputs::can_be_common_mode_choke(operatingPoints[0])) {
        numberOfWindings = operatingPoints[0].get_excitations_per_winding().size();
        if (numberOfWindings < 1) numberOfWindings = 1;
    }

    // First winding is the primary; any additional windings are secondaries
    // with distinct names and SECONDARY isolation side. The original draft
    // of this loop made every winding identical ("primary"/PRIMARY), which
    // caused downstream code that looks up windings by name to either
    // collide on the duplicated key or fail to find the "secondary" the
    // CMC path expects — observed as a SEGV inside MagnetizingInductance
    // when called on a Mas built from this dummy coil.
    std::vector<Winding> windings;
    for (size_t w = 0; w < numberOfWindings; ++w) {
        Winding winding;
        if (w == 0) {
            winding.set_isolation_side(IsolationSide::PRIMARY);
            winding.set_name("primary");
        } else {
            winding.set_isolation_side(IsolationSide::SECONDARY);
            winding.set_name(numberOfWindings == 2 ? "secondary"
                                                   : "secondary_" + std::to_string(w));
        }
        winding.set_number_parallels(1);
        winding.set_number_turns(1);
        winding.set_wire(wire);
        windings.push_back(winding);
    }

    Coil coil;
    coil.set_bobbin(DUMMY_SENTINEL_NAME);
    coil.set_functional_description(windings);
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

void add_initial_turns_by_inductance(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, const Inputs& inputs) {
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
    // Inductor vs transformer drives both the turn/gap sizing below and the
    // saturation co-design. Use the single shared predicate (is_inductor, in
    // MagneticFilterInternal.h) so this seeder, the saturation filter, the
    // loss-sweep saturation cap, and the realism gate all agree on the
    // classification — a coupled inductor (single isolation side) or an
    // energy-storing / specific-Lm design is sized with the B-from-current path;
    // a transformer (high-Lm, minimum-only) with the B-from-voltage path.
    bool isTransformer = !is_inductor(inputs);
    
    // Pre-compute shared transformer values once (not per-core).
    // maxVoltSeconds is the peak V·s excursion over all OPs' primary
    // voltage waveforms — Faraday-integrated, non-sinusoidal-safe.
    double transformerTemperature = 25.0;
    double maxVoltSeconds = 0.0;
    if (isTransformer) {
        const auto& transformerOperatingPoints = inputs.get_operating_points();
        for (const auto& op : transformerOperatingPoints) {
            transformerTemperature = std::max(transformerTemperature, op.get_conditions().get_ambient_temperature());
            maxVoltSeconds = std::max(maxVoltSeconds,
                Inputs::calculate_max_volt_seconds(Inputs::get_primary_excitation(op)));
        }
    }

    // For suppression applications (CMC/DMC), the impedance requirement is the
    // primary constraint. Convert the most stringent |Z|≥Z_min@f into an
    // equivalent inductance L_eq = Z_min / (2πf) and use that to size turns.
    // Without this, a CMC with L=100µH and Z_req=100Ω@100kHz gets N turns
    // for 100µH (Z≈63Ω) and is then rejected by the impedance filter.
    double suppressionEquivalentInductance = 0.0;
    auto minimumImpedanceOpt = inputs.get_design_requirements().get_minimum_impedance();
    if (minimumImpedanceOpt && !minimumImpedanceOpt->empty()) {
        for (const auto& zSpec : minimumImpedanceOpt.value()) {
            double f = zSpec.get_frequency();
            double z = zSpec.get_impedance().get_magnitude();
            if (f > 0) {
                double l_eq = z / (2.0 * std::numbers::pi * f);
                suppressionEquivalentInductance = std::max(suppressionEquivalentInductance, l_eq);
            }
        }
    }
    auto reqAppOpt = inputs.get_design_requirements().get_application();
    bool isSuppression = reqAppOpt.has_value() && reqAppOpt.value() == Application::INTERFERENCE_SUPPRESSION;
    for (size_t i = 0; i < (*magneticsWithScoring).size(); ++i){

        Core core = (*magneticsWithScoring)[i].first.get_core();
        if (!core.get_processed_description()) {
            core.process_data();
            core.process_gap();
        }
        double initialNumberTurns = (*magneticsWithScoring)[i].first.get_coil().get_functional_description()[0].get_number_turns();

        if (initialNumberTurns == 1) {
            if (isTransformer) {
                // Feasibility seed: fewest turns that keep B_peak under the
                // material's safe operating flux. Use the shared single-source
                // ceiling maximum_allowed_magnetic_flux_density() — the smaller
                // of Bsat/margin and the maximumProportion cap — so the seed
                // sizes to exactly the same limit the gapping code targets (they
                // must not diverge). The loss filter refines N upward toward the
                // loss optimum later.
                // Evaluate B_sat at the hot junction corner (ABT #13), consistent
                // with the saturation filter and the gapping sizing.
                double bSatRaw = core.get_magnetic_flux_density_saturation(
                    saturation_derating_temperature(transformerTemperature), false);
                double bMax = maximum_allowed_magnetic_flux_density(bSatRaw);
                double nFromSaturation = magnetizingInductance
                    .calculate_turns_from_volt_seconds_and_max_flux_density(core, maxVoltSeconds, bMax);
                initialNumberTurns = nFromSaturation > 0 ? nFromSaturation : 5;

                // For resonant transformers with a specific Lm target
                // (e.g. CLLLC, CLLC) we must size BOTH N and the gap so
                // that L_actual(N, gap) == target Lm at the saturation-
                // safe N.
                //
                // Plan:
                //   1. Saturation floor for N (volt-seconds) is already
                //      in `initialNumberTurns` above.
                //   2. From that N, ask MagnetizingInductance for the
                //      gap length that produces the target Lm.
                //   3. Apply the gap to the core and persist N on the
                //      magnetic's coil so add_gapping_standard_cores
                //      (run earlier with the ungapped-transformer path)
                //      and filterMagneticInductance agree.
                //
                // Toroidal cores can't take a discrete gap — leave them
                // ungapped (they'll be rejected by the inductance filter
                // if AL_ungapped is too high for the target Lm, which
                // is the correct behaviour).
                auto lmDim = inputs.get_design_requirements().get_magnetizing_inductance();
                // Lm = 0 means "not specified" (same convention as
                // pre_process_inputs in CoreAdviserPipeline.cpp). Solving the
                // gap for a literal 0 H target makes the needed reluctance
                // infinite: the search grows the gap until it no longer fits
                // any column and every candidate dies with "Gap Area is not
                // set" (Test_FastAdviser_LLC_Lm_Zero).
                bool hasNominalLm = lmDim.get_nominal().has_value() && lmDim.get_nominal().value() > 0;
                bool hasMaxLm = lmDim.get_maximum().has_value() && lmDim.get_maximum().value() > 0;
                if ((hasNominalLm || hasMaxLm) && core.get_shape_family() != CoreShapeFamily::T) {
                    // Build a temporary coil with the saturation-safe N for the
                    // gap calculation. The wire is unimportant; only N matters
                    // for the reluctance-based gap solver.
                    Coil tempCoil = (*magneticsWithScoring)[i].first.get_coil();
                    if (!tempCoil.get_functional_description().empty()) {
                        tempCoil.get_mutable_functional_description()[0]
                                .set_number_turns(static_cast<int64_t>(initialNumberTurns));
                    }

                    // Pass an Inputs with only design requirements (no operating points)
                    // — the gap solver's reluctance model otherwise dives into the
                    // magnetizing-current pipeline and throws on resonant inputs
                    // (CLLLC/CLLC) that have voltage but no per-winding
                    // magnetizing_current waveform yet.
                    Inputs lmInputs;
                    lmInputs.set_design_requirements(inputs.get_design_requirements());

                    try {
                        auto gaps = magnetizingInductance
                            .calculate_gapping_from_number_turns_and_inductance(
                                core, tempCoil, &lmInputs, GappingType::GROUND);
                        if (!gaps.empty()) {
                            core.set_gapping(gaps);
                            core.process_gap();
                            (*magneticsWithScoring)[i].first.set_core(core);
                        }
                    } catch (const std::exception& e) {
                        // Solver failed (e.g. target Lm unreachable with this
                        // core/material combination). Leave the core ungapped;
                        // the inductance filter downstream rejects it cleanly.
                        logEntry(std::string("Transformer gap solver failed for core ")
                                 + core.get_name().value_or("?") + ": " + e.what()
                                 + " — leaving ungapped",
                                 "CoreAdviser", 2);
                    }
                }
            } else if (isSuppression && suppressionEquivalentInductance > 0.0) {
                // Use impedance-derived inductance for suppression apps
                auto reqCopy = inputs.get_design_requirements();
                DimensionWithTolerance lmSpec;
                lmSpec.set_minimum(suppressionEquivalentInductance);
                reqCopy.set_magnetizing_inductance(lmSpec);
                Inputs tempInputs(inputs);
                tempInputs.set_design_requirements(reqCopy);
                initialNumberTurns = magnetizingInductance.calculate_number_turns_from_gapping_and_inductance(core, (*magneticsWithScoring)[i].first.get_coil(), &tempInputs, DimensionalValues::MINIMUM);
            } else {
                // Inductors / flybacks. Turns for the target inductance come from
                // the single source of truth: the canonical inverse of the
                // inductance model (calculate_number_turns_from_gapping_and_inductance),
                // which is consistent with the downstream inductance filter.
                Inputs inputsCopy(inputs);
                initialNumberTurns = magnetizingInductance.calculate_number_turns_from_gapping_and_inductance(
                    core, (*magneticsWithScoring)[i].first.get_coil(), &inputsCopy, DimensionalValues::MINIMUM);

                // Energy-storing topologies (flyback / buck / boost / …): the
                // inductance-determined N minimises turns, which maximises peak flux
                // density and can leave the design saturating (isat < margin·I_pk).
                // The CoreAdviser — not the CoilAdviser, which only lays out turns and
                // never re-gaps — owns the final (N, gap), so size them here against the
                // GAP-AWARE saturation current (Magnetic::calculate_saturation_current,
                // the very identity the inductance filter's isat gate and the realism
                // gate use). Intervene only when the as-seeded design actually saturates
                // (zero perturbation otherwise); then raise N — re-solving the gap to
                // hold the target L, so at fixed L isat rises ∝ N — until it clears the
                // margin, and persist that exact (N, gap). Powder toroids (shape T) take
                // no discrete gap and keep the seed.
                if (is_inductor(inputs) && core.get_shape_family() != CoreShapeFamily::T && !inputs.get_operating_points().empty()) {
                    auto operatingPoint = inputsCopy.get_operating_point(0);
                    double peakCurrent = 0;
                    for (auto& op : inputs.get_operating_points()) {
                        auto exc = Inputs::get_primary_excitation(op);
                        if (exc.get_current() && exc.get_current()->get_processed() && exc.get_current()->get_processed()->get_peak()) {
                            peakCurrent = std::max(peakCurrent, std::abs(exc.get_current()->get_processed()->get_peak().value()));
                        }
                    }
                    // Derating (ABT #13): size against the gap-aware saturation
                    // current evaluated at the hot junction corner and with RAW
                    // B_sat (proportion=false), matching the saturation filter's
                    // inductor gate exactly — total derating is (hot temp × margin),
                    // with no 0.7 flux-proportion stacked on top.
                    double temperature = saturation_derating_temperature(operatingPoint.get_conditions().get_ambient_temperature());
                    double requiredIsat = settings.get_core_adviser_saturation_margin() * peakCurrent;

                    std::optional<double> seededIsat;
                    if (peakCurrent > 0 && requiredIsat > 0) {
                        Magnetic seeded = (*magneticsWithScoring)[i].first;
                        seeded.get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(static_cast<int64_t>(std::llround(initialNumberTurns)));
                        try { seededIsat = seeded.calculate_saturation_current(temperature, /*proportion=*/false); }
                        catch (const std::exception&) { seededIsat = std::nullopt; }
                    }

                    if (seededIsat && seededIsat.value() > 0 && seededIsat.value() < requiredIsat) {
                        Inputs lmInputs;
                        lmInputs.set_design_requirements(inputs.get_design_requirements());
                        auto regapAndIsat = [&](double n, Core& outCore) -> std::optional<double> {
                            Coil tempCoil = (*magneticsWithScoring)[i].first.get_coil();
                            tempCoil.get_mutable_functional_description()[0].set_number_turns(static_cast<int64_t>(std::llround(n)));
                            outCore = core;
                            try {
                                auto gaps = magnetizingInductance.calculate_gapping_from_number_turns_and_inductance(outCore, tempCoil, &lmInputs, GappingType::GROUND);
                                if (gaps.empty()) return std::nullopt;
                                outCore.set_gapping(gaps);
                                outCore.process_gap();
                            } catch (const std::exception&) { return std::nullopt; }
                            Magnetic tm = (*magneticsWithScoring)[i].first;
                            tm.set_core(outCore);
                            tm.get_mutable_coil().get_mutable_functional_description()[0].set_number_turns(static_cast<int64_t>(std::llround(n)));
                            try { return tm.calculate_saturation_current(temperature, /*proportion=*/false); }
                            catch (const std::exception&) { return std::nullopt; }
                        };

                        double n = initialNumberTurns;
                        double isatAtN = seededIsat.value();
                        Core acceptedCore;
                        bool accepted = false;
                        const size_t maxSatSteps = 30;
                        for (size_t step = 0; step < maxSatSteps; ++step) {
                            double nextN = std::ceil(n * (requiredIsat / isatAtN));
                            if (std::llround(nextN) <= std::llround(n)) nextN = n + 1;
                            n = nextN;
                            Core trialCore;
                            auto isat = regapAndIsat(n, trialCore);
                            if (!isat || isat.value() <= 0) break;
                            accepted = true;
                            acceptedCore = trialCore;
                            isatAtN = isat.value();
                            if (isatAtN >= requiredIsat) break;
                        }
                        if (accepted) {
                            initialNumberTurns = n;
                            core = acceptedCore;
                            (*magneticsWithScoring)[i].first.set_core(core);
                        }
                    }
                }
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

std::vector<std::pair<Magnetic, double>> add_initial_turns_by_impedance(std::vector<std::pair<Magnetic, double>> magneticsWithScoring, const Inputs& inputs) {
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
                logEntry("add_initial_turns_by_impedance: core " + core.get_name().value_or("?") + " returned turns=" + std::to_string(initialNumberTurns) + ", skipping", "CoreAdviser", 2);
                continue;
            }
        }
        catch (const std::exception& e) {
            logEntry(std::string("add_initial_turns_by_impedance: core ") + core.get_name().value_or("?") + " threw: " + e.what(), "CoreAdviser", 2);
            continue;
        }
        catch (...) {
            logEntry("add_initial_turns_by_impedance: core " + core.get_name().value_or("?") + " threw unknown exception", "CoreAdviser", 2);
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




} // namespace OpenMagnetics
