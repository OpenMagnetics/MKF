// Phase 5 (structural cleanup, material selection): extracted from
// CoreAdviser.cpp. This translation unit owns the helpers that build the
// pool of candidate core materials (powder + ferrite, by losses or by
// impedance) before the main scoring pipeline runs.
//
// Sibling TUs of CoreAdviser.cpp; declarations live in advisers/CoreAdviser.h
// and shared helpers live in advisers/CoreAdviserInternal.h.

#include "advisers/CoreAdviser.h"
#include "advisers/CoreAdviserInternal.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "physical_models/ComplexPermeability.h"
#include "support/Logger.h"
#include "support/Settings.h"
#include "support/Exceptions.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace OpenMagnetics {

bool CoreAdviser::should_include_powder(Inputs inputs) {
    // Check if powder cores are disabled in settings
    if (!settings.get_use_powder_cores()) {
        return false;
    }
    
    if (get_application() != MAS::MagneticApplication::POWER) {
        return false;
    }
    double maximumCurrentDcBias = inputs.get_maximum_current_dc_bias();
    if (maximumCurrentDcBias > 1e-3) {
        return true;
    }
    return false;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::add_powder_materials(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    size_t numberCoreMaterialsTouse = defaults.coreAdviserAlternativeMaterialsNumberToUse;
    double magneticFluxDensityReference = defaults.coreAdviserMagneticFluxDensityReferenceAlternative;
    std::vector<std::pair<Magnetic, double>> magneticsWithMaterials;
    std::vector<CoreMaterial> coreMaterialsToEvaluate;
    auto coreMaterials = get_core_material_names(settings.get_preferred_core_material_powder_manufacturer());
    for (auto coreMaterial : coreMaterials) {
        auto resolved = Core::resolve_material(coreMaterial);
        if (Core::check_material_application(resolved, _application)) {
            coreMaterialsToEvaluate.push_back(resolved);
        }
    }
    std::vector<CoreMaterial> coreMaterialsToUse;
    std::vector<std::pair<CoreMaterial, double>> evaluations;

    double temperature = inputs.get_maximum_temperature();
    double maximumCurrentDcBias = inputs.get_maximum_current_dc_bias();

    // Phase 3 (F8): sinusoidal excitation + loss-model pair are shared
    // with add_ferrite_materials_by_impedance.
    OperatingPointExcitation operatingPointExcitation = make_sinusoidal_excitation(
        magneticFluxDensityReference, maximumCurrentDcBias, 1);
    auto [coreLossesModelSteinmetz, coreLossesModelProprietary] = make_default_core_losses_model_pair();
    for (auto coreMaterial : coreMaterialsToEvaluate) {
        double averageVolumetricCoreLosses = 0;
        try {
            for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex){
                auto frequency = inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding()[0].get_frequency();
                Inputs::scale_time_to_frequency(operatingPointExcitation, frequency, false, false);
                auto coreLossesMethods = Core::get_available_core_losses_methods(coreMaterial);

                if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), VolumetricCoreLossesMethodType::STEINMETZ) != coreLossesMethods.end()) {
                    double volumetricCoreLosses = coreLossesModelSteinmetz->get_core_volumetric_losses(coreMaterial, operatingPointExcitation, temperature);
                    averageVolumetricCoreLosses += volumetricCoreLosses;
                }
                else {
                    double volumetricCoreLosses = coreLossesModelProprietary->get_core_volumetric_losses(coreMaterial, operatingPointExcitation, temperature);
                    averageVolumetricCoreLosses += volumetricCoreLosses;
                }
            }
        }
        catch (const std::exception& e) {
            // Candidate-pool construction: a material whose loss model cannot
            // be evaluated (e.g. no Steinmetz coefficients and no proprietary
            // method for its manufacturer) cannot compete in the ranking —
            // skip it loudly instead of letting one catalog gap abort the
            // whole adviser run. Named, logged, bounded scope.
            logEntry(std::string("Skipping core material '") + coreMaterial.get_name()
                     + "' in powder candidate ranking: " + e.what(), "CoreAdviser");
            continue;
        }
        averageVolumetricCoreLosses /= inputs.get_operating_points().size();
        evaluations.push_back({coreMaterial, averageVolumetricCoreLosses}); // NEW-BUG-FIX: was pushing 0, discarding computed losses
    }

    stable_sort((evaluations).begin(), (evaluations).end(), [](const std::pair<CoreMaterial, double>& b1, const std::pair<CoreMaterial, double>& b2) {
        return b1.second < b2.second;
    });

    numberCoreMaterialsTouse = std::min(numberCoreMaterialsTouse, evaluations.size());
    for (size_t magneticIndex = 0; magneticIndex < (*magneticsWithScoring).size(); ++magneticIndex) {
        for (size_t i = 0; i < numberCoreMaterialsTouse; ++i){
            auto [magnetic, scoring] = (*magneticsWithScoring)[magneticIndex];
            magnetic.get_mutable_core().set_material(evaluations[i].first);
            magnetic.get_mutable_core().set_name(evaluations[i].first.get_name() + " " + magnetic.get_core().get_name().value_or("unnamed"));
            if (magnetic.get_manufacturer_info()) {
                auto manufacturerInfo = magnetic.get_manufacturer_info().value();
                manufacturerInfo.set_reference(evaluations[i].first.get_name() + " " + magnetic.get_reference());
                magnetic.set_manufacturer_info(manufacturerInfo);
            }
            magneticsWithMaterials.push_back({magnetic, scoring});
        }
    }
    return magneticsWithMaterials;
}

// Phase 3 (F4): shared by add_ferrite_materials_by_losses /
// add_ferrite_materials_by_impedance — gather ferrite candidates filtered
// by current application; small but called twice.
std::vector<CoreMaterial> CoreAdviser::gather_ferrite_materials_for_application() const {
    std::vector<CoreMaterial> coreMaterialsToEvaluate;
    auto coreMaterials = get_core_material_names(settings.get_preferred_core_material_ferrite_manufacturer());
    for (auto coreMaterial : coreMaterials) {
        auto resolved = Core::resolve_material(coreMaterial);
        if (Core::check_material_application(resolved, _application)) {
            coreMaterialsToEvaluate.push_back(resolved);
        }
    }
    return coreMaterialsToEvaluate;
}

// Phase 3 (F4): the dummy-material fan-out block was duplicated verbatim
// between add_ferrite_materials_by_losses and add_ferrite_materials_by_impedance
// (24 lines × 2). Pass-through for non-Dummy candidates (the Phase 1 fix
// for the duplicate-entries bug) is preserved here.
std::vector<std::pair<Magnetic, double>> CoreAdviser::fan_out_dummy_into_top_materials(
    std::vector<std::pair<Magnetic, double>> *magneticsWithScoring,
    const std::vector<std::pair<CoreMaterial, double>>& sortedMaterialEvaluations,
    size_t numberCoreMaterialsTouse) {
    std::vector<std::pair<Magnetic, double>> magneticsWithMaterials;
    numberCoreMaterialsTouse = std::min(numberCoreMaterialsTouse, sortedMaterialEvaluations.size());
    for (size_t magneticIndex = 0; magneticIndex < (*magneticsWithScoring).size(); ++magneticIndex) {
        auto [magnetic, scoring] = (*magneticsWithScoring)[magneticIndex];
        if (magnetic.get_mutable_core().get_material_name() != DUMMY_SENTINEL_NAME) {
            // Already has a concrete material — pass through ONCE.
            // Phase 1 fix: previously the inner i-loop also ran in this branch
            // (with `continue`), so non-Dummy magnetics were pushed
            // numberCoreMaterialsTouse (=2) times, producing the duplicate
            // entries in CoreAdviser STANDARD_CORES × POWER top-N output.
            magneticsWithMaterials.push_back({magnetic, scoring});
            continue;
        }
        // Dummy material — fan out into the top-N candidate ferrite materials.
        for (size_t i = 0; i < numberCoreMaterialsTouse; ++i){
            auto magneticCopy = magnetic;
            magneticCopy.get_mutable_core().set_material(sortedMaterialEvaluations[i].first);
            magneticCopy.get_mutable_core().set_name(sortedMaterialEvaluations[i].first.get_name() + " " + magneticCopy.get_core().get_name().value_or("unnamed"));
            if (magneticCopy.get_manufacturer_info()) {
                auto manufacturerInfo = magneticCopy.get_manufacturer_info().value();
                manufacturerInfo.set_reference(sortedMaterialEvaluations[i].first.get_name() + " " + magneticCopy.get_reference());
                magneticCopy.set_manufacturer_info(manufacturerInfo);
            }
            magneticsWithMaterials.push_back({magneticCopy, scoring});
        }
    }
    return magneticsWithMaterials;
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::add_ferrite_materials_by_losses(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    size_t numberCoreMaterialsTouse = defaults.coreAdviserFanOutMaterialsNumberToUse;
    double magneticFluxDensityReference = defaults.coreAdviserMagneticFluxDensityReferenceAlternative;
    auto coreMaterialsToEvaluate = gather_ferrite_materials_for_application();
    std::vector<CoreMaterial> coreMaterialsToUse;
    std::vector<std::pair<CoreMaterial, double>> evaluations;

    double temperature = inputs.get_maximum_temperature();


    // Phase 3 (F8): sinusoidal excitation + loss-model pair are shared
    // with add_ferrite_materials_by_losses.
    OperatingPointExcitation operatingPointExcitation = make_sinusoidal_excitation(
        magneticFluxDensityReference, 0, 1);
    auto [coreLossesModelSteinmetz, coreLossesModelProprietary] = make_default_core_losses_model_pair();
    for (auto coreMaterial : coreMaterialsToEvaluate) {
        double averageVolumetricCoreLosses = 0;
        try {
            for (size_t operatingPointIndex = 0; operatingPointIndex < inputs.get_operating_points().size(); ++operatingPointIndex){
                auto frequency = inputs.get_operating_points()[operatingPointIndex].get_excitations_per_winding()[0].get_frequency();
                Inputs::scale_time_to_frequency(operatingPointExcitation, frequency, false, false);
                auto coreLossesMethods = Core::get_available_core_losses_methods(coreMaterial);

                if (std::find(coreLossesMethods.begin(), coreLossesMethods.end(), VolumetricCoreLossesMethodType::STEINMETZ) != coreLossesMethods.end()) {
                    double volumetricCoreLosses = coreLossesModelSteinmetz->get_core_volumetric_losses(coreMaterial, operatingPointExcitation, temperature);
                    averageVolumetricCoreLosses += volumetricCoreLosses;
                }
                else {
                    double volumetricCoreLosses = coreLossesModelProprietary->get_core_volumetric_losses(coreMaterial, operatingPointExcitation, temperature);
                    averageVolumetricCoreLosses += volumetricCoreLosses;
                }
            }
        }
        catch (const std::exception& e) {
            // Same skip-loudly policy as the powder ranking above: a material
            // with no usable loss model cannot compete; one catalog gap must
            // not abort the whole adviser run. Named, logged, bounded scope.
            logEntry(std::string("Skipping core material '") + coreMaterial.get_name()
                     + "' in ferrite candidate ranking: " + e.what(), "CoreAdviser");
            continue;
        }
        averageVolumetricCoreLosses /= inputs.get_operating_points().size();
        evaluations.push_back({coreMaterial, averageVolumetricCoreLosses});
    }

    stable_sort((evaluations).begin(), (evaluations).end(), [](const std::pair<CoreMaterial, double>& b1, const std::pair<CoreMaterial, double>& b2) {
        return b1.second < b2.second;
    });

    return fan_out_dummy_into_top_materials(magneticsWithScoring, evaluations, numberCoreMaterialsTouse);
}

std::vector<std::pair<Magnetic, double>> CoreAdviser::add_ferrite_materials_by_impedance(std::vector<std::pair<Magnetic, double>> *magneticsWithScoring, Inputs inputs) {
    size_t numberCoreMaterialsTouse = 2;
    auto coreMaterialsToEvaluate = gather_ferrite_materials_for_application();
    std::vector<CoreMaterial> coreMaterialsToUse;
    std::vector<std::pair<CoreMaterial, double>> evaluations;

    if (!inputs.get_design_requirements().get_minimum_impedance()) {
        throw std::invalid_argument("Missing impedance requirement");
    }

    auto minimumImpedanceRequirement = inputs.get_design_requirements().get_minimum_impedance().value();

    ComplexPermeability complexPermeabilityModel;
    for (auto coreMaterial : coreMaterialsToEvaluate) {
        double totalComplexPermeability = 0;
        try {
            for (auto impedanceAtFrequency : minimumImpedanceRequirement) {
                auto frequency = impedanceAtFrequency.get_frequency();
                auto [complexPermeabilityRealPart, complexPermeabilityImaginaryPart] = complexPermeabilityModel.get_complex_permeability(coreMaterial, frequency);
                totalComplexPermeability += hypot(complexPermeabilityRealPart, complexPermeabilityImaginaryPart);
            }
        }
        catch (const std::exception& e) {
            // Same skip-loudly policy as the loss-based rankings: a material
            // with no complex-permeability data cannot compete on impedance;
            // one catalog gap must not abort the whole adviser run.
            logEntry(std::string("Skipping core material '") + coreMaterial.get_name()
                     + "' in impedance candidate ranking: " + e.what(), "CoreAdviser");
            continue;
        }
        evaluations.push_back({coreMaterial, totalComplexPermeability});
    }

    // Impedance-driven (suppression) applications want the HIGHEST complex
    // permeability at the requirement frequencies, so sort descending (the
    // losses-based sibling sorts ascending because lower losses are better).
    stable_sort((evaluations).begin(), (evaluations).end(), [](const std::pair<CoreMaterial, double>& b1, const std::pair<CoreMaterial, double>& b2) {
        return b1.second > b2.second;
    });

    return fan_out_dummy_into_top_materials(magneticsWithScoring, evaluations, numberCoreMaterialsTouse);
}


} // namespace OpenMagnetics
