#include "constructive_models/Core.h"
#include "physical_models/InitialPermeability.h"
#include "advisers/CoreMaterialCrossReferencer.h"
#include "advisers/CrossReferencerCommon.h"
#include "Defaults.h"
#include "support/Settings.h"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <string>
#include <vector>
#include <execution>
#include <magic_enum_utility.hpp>
#include <list>
#include <cmrc/cmrc.hpp>
#include "support/Exceptions.h"
#include "support/Logger.h"

CMRC_DECLARE(data);


namespace OpenMagnetics {


std::map<std::string, std::map<CoreMaterialCrossReferencerFilters, double>> CoreMaterialCrossReferencer::get_scorings(bool weighted){
    // Phase 7 (Option D): same inline log-normalisation as CoreCross
    // Referencer::get_scorings was duplicated here. Both bodies now
    // forward to CrossReferencerCommon.h::compute_normalized_scorings.
    // The new common path includes this file's NaN-guard fixes plus
    // CCR's data-relative floor (B8/F4) and named-constant neutral
    // score (was 0.5 hardcoded here).
    return compute_normalized_scorings<CoreMaterialCrossReferencerFilters>(
        _scorings, _filterConfiguration, _weights, weighted);
}

void normalize_scoring(std::vector<std::pair<CoreMaterial, double>>* rankedCoreMaterials, std::vector<double>* newScoring, double weight, std::map<std::string, bool> filterConfiguration) {
    // Phase 3 (F2+F3): forwards to the shared templated implementation
    // in CrossReferencerCommon.h. CoreCrossReferencer carried a
    // byte-identical second copy; both fix histories (F5/B8/XC-6/F12)
    // now live in one place.
    normalize_scoring<CoreMaterial>(rankedCoreMaterials, newScoring, weight, filterConfiguration);
}

std::map<std::string, std::map<CoreMaterialCrossReferencerFilters, double>> CoreMaterialCrossReferencer::get_scored_values(){
    std::map<std::string, std::map<CoreMaterialCrossReferencerFilters, double>> swappedScoredValues;
    for (auto& [filter, aux] : _scoredValues) {
        auto filterConfiguration = _filterConfiguration[filter];

        for (auto& [name, scoredValue] : aux) {
            swappedScoredValues[name][filter] = scoredValue;
        }
    }
    return swappedScoredValues;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterInitialPermeability::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;
    OpenMagnetics::InitialPermeability initialPermeability;

    double referenceInitialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(referenceCoreMaterial, temperature, std::nullopt, std::nullopt);
    add_scored_value("Reference", CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY, referenceInitialPermeabilityValueWithTemperature);


    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double initialPermeabilityValueWithTemperature = initialPermeability.get_initial_permeability(coreMaterial, temperature, std::nullopt, std::nullopt);

        double scoring = fabs(referenceInitialPermeabilityValueWithTemperature - initialPermeabilityValueWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY, scoring);
        add_scored_value(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY, initialPermeabilityValueWithTemperature);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencerFilters::INITIAL_PERMEABILITY]);
    }
    return *unfilteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterRemanence::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    double referenceRemanenceWithTemperature = Core::get_remanence(referenceCoreMaterial, temperature, true);
    add_scored_value("Reference", CoreMaterialCrossReferencerFilters::REMANENCE, referenceRemanenceWithTemperature);

 
    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double remanenceWithTemperature = Core::get_remanence(coreMaterial, temperature, true);

        double scoring = fabs(referenceRemanenceWithTemperature - remanenceWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::REMANENCE, scoring);
        add_scored_value(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::REMANENCE, remanenceWithTemperature);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencerFilters::REMANENCE]);
    }
    return *unfilteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterCoerciveForce::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    double referenceCoerciveForceWithTemperature = Core::get_coercive_force(referenceCoreMaterial, temperature, true);
    add_scored_value("Reference", CoreMaterialCrossReferencerFilters::COERCIVE_FORCE, referenceCoerciveForceWithTemperature);


    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double coerciveForceWithTemperature = Core::get_coercive_force(coreMaterial, temperature, true);

        double scoring = fabs(referenceCoerciveForceWithTemperature - coerciveForceWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::COERCIVE_FORCE, scoring);
        add_scored_value(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::COERCIVE_FORCE, coerciveForceWithTemperature);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencerFilters::COERCIVE_FORCE]);
    }
    return *unfilteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterSaturation::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    double referenceSaturationWithTemperature = Core::get_magnetic_flux_density_saturation(referenceCoreMaterial, temperature);
    add_scored_value("Reference", CoreMaterialCrossReferencerFilters::SATURATION, referenceSaturationWithTemperature);


    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double saturationWithTemperature = Core::get_magnetic_flux_density_saturation(coreMaterial, temperature);

        double scoring = fabs(referenceSaturationWithTemperature - saturationWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::SATURATION, scoring);
        add_scored_value(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::SATURATION, saturationWithTemperature);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencerFilters::SATURATION]);
    }
    return *unfilteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterCurieTemperature::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    double referenceCurieTemperatureWithTemperature = Core::get_curie_temperature(referenceCoreMaterial);
    add_scored_value("Reference", CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE, referenceCurieTemperatureWithTemperature);


    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double curieTemperatureWithTemperature = Core::get_curie_temperature(coreMaterial);

        double scoring = fabs(referenceCurieTemperatureWithTemperature - curieTemperatureWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE, scoring);
        add_scored_value(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE, curieTemperatureWithTemperature);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencerFilters::CURIE_TEMPERATURE]);
    }
    return *unfilteredCoreMaterials;
}

// CMCR-OPT-1 NOTE: This computes losses at 25 (B,f) points per material.
// Consider pre-computing reference losses once, then using ratio-based comparison.
// Or use Steinmetz equation for initial screening before detailed model.
double CoreMaterialCrossReferencer::MagneticCoreFilterVolumetricLosses::calculate_average_volumetric_losses(CoreMaterial coreMaterial, double temperature, std::map<std::string, std::string> models) {
    if (models.find("coreLosses") == models.end()) {
        models["coreLosses"] = to_string(Defaults().coreLossesModelDefault);
    }

    OperatingPointExcitation excitation;
    SignalDescriptor magneticFluxDensity;
    ProcessedWaveform magneticFluxDensityProcessed;
    magneticFluxDensityProcessed.set_label(WaveformLabel::SINUSOIDAL);
    magneticFluxDensityProcessed.set_offset(0);
    magneticFluxDensityProcessed.set_duty_cycle(0.5);
    std::shared_ptr<CoreLossesModel> coreLossesModelForMaterial = nullptr;

    try {
        auto availableMethodsForMaterial = CoreLossesModel::get_methods(coreMaterial);
        for (auto& [modelName, coreLossesModel] : _coreLossesModels) {
            if (std::find(availableMethodsForMaterial.begin(), availableMethodsForMaterial.end(), modelName) != availableMethodsForMaterial.end()) {
                coreLossesModelForMaterial = coreLossesModel;
                break;
            }
        }

        if (coreLossesModelForMaterial == nullptr) {
            throw ModelNotAvailableException("No model found for material: " + coreMaterial.get_name());
        }

        double averageVolumetricLosses = 0;
        for (auto magneticFluxDensityPeak : _magneticFluxDensities) {
            magneticFluxDensityProcessed.set_peak(magneticFluxDensityPeak);
            magneticFluxDensityProcessed.set_peak_to_peak(magneticFluxDensityPeak * 2);
            magneticFluxDensity.set_processed(magneticFluxDensityProcessed);
            for (auto frequency : _frequencies) {
                magneticFluxDensity.set_waveform(Inputs::create_waveform(magneticFluxDensityProcessed, frequency));
                excitation.set_frequency(frequency);
                excitation.set_magnetic_flux_density(magneticFluxDensity);
                double coreVolumetricLosses = coreLossesModelForMaterial->get_core_volumetric_losses(coreMaterial, excitation, temperature);
                averageVolumetricLosses += coreVolumetricLosses;
            }
        }
        averageVolumetricLosses /= _magneticFluxDensities.size() * _frequencies.size();
        return averageVolumetricLosses;
    }
    catch(const ModelNotAvailableException& re)
    {
        // Phase 1: NaN sentinel for "uncomputable" — caller in filter_core_materials
        // (line ~403) explicitly std::isnan-checks and converts to DBL_MAX score.
        // Narrow exception type so unexpected throws still surface.
        return std::numeric_limits<double>::quiet_NaN();
    }
    catch(const InvalidInputException& re)
    {
        // Same as above — explicit NaN-as-uncomputable contract.
        return std::numeric_limits<double>::quiet_NaN();
    }
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterVolumetricLosses::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, std::map<std::string, std::string> models, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    try {
        double referenceVolumetricLossesWithTemperature = calculate_average_volumetric_losses(referenceCoreMaterial, temperature, models);
        if (std::isnan(referenceVolumetricLossesWithTemperature)) {
            // The reference material has no usable volumetric-loss model (e.g. an
            // interference-suppression ferrite characterised only by a loss factor
            // or complex permeability). Rather than abort the whole cross-reference,
            // skip the loss dimension and let the remaining dimensions (initial
            // permeability, etc.) rank the alternatives. Logged, scoped pass-through.
            logEntry("CoreMaterialCrossReferencer: reference material '"
                         + referenceCoreMaterial.get_name()
                         + "' has uncomputable volumetric losses at "
                         + std::to_string(temperature) + " degC; skipping the loss "
                         "dimension and cross-referencing on the remaining filters",
                     "CoreMaterialCrossReferencer");
            return *unfilteredCoreMaterials;
        }
        add_scored_value("Reference", CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES, referenceVolumetricLossesWithTemperature);


        OperatingPointExcitation excitation;
        SignalDescriptor magneticFluxDensity;
        ProcessedWaveform magneticFluxDensityProcessed;
        magneticFluxDensityProcessed.set_label(WaveformLabel::SINUSOIDAL);
        magneticFluxDensityProcessed.set_offset(0);
        magneticFluxDensityProcessed.set_duty_cycle(0.5);

        std::vector<std::pair<CoreMaterial, double>> filteredCoreMaterialsWithScoring;
        for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
            CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;
            double volumetricLossesWithTemperature = calculate_average_volumetric_losses(coreMaterial, temperature, models);
            if (std::isnan(volumetricLossesWithTemperature)) {
                // Uncomputable losses for this candidate: cull it (like the core
                // cross-referencer) rather than assigning DBL_MAX, which would blow up the
                // min-max/log normalization for every other candidate in the batch.
                continue;
            }

            // Score by absolute distance from the reference ALWAYS. The old branch hardcoded
            // scoring=0 ("perfect distance") for every better-than-reference material, tying
            // them all at the normalized maximum — the same score-ceiling bug the core
            // cross-referencer already fixed (CoreCrossReferencer Phase 1).
            double scoring = fabs(referenceVolumetricLossesWithTemperature - volumetricLossesWithTemperature);
            newScoring.push_back(scoring);
            add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES, scoring);

            add_scored_value(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES, volumetricLossesWithTemperature);
            filteredCoreMaterialsWithScoring.push_back((*unfilteredCoreMaterials)[coreMaterialIndex]);
        }

        if (filteredCoreMaterialsWithScoring.size() != newScoring.size()) {
            throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of filteredCoreMaterials: " + std::to_string(filteredCoreMaterialsWithScoring.size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
        }

        if (filteredCoreMaterialsWithScoring.size() > 0) {
            normalize_scoring(&filteredCoreMaterialsWithScoring, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencerFilters::VOLUMETRIC_LOSSES]);
        }
        return filteredCoreMaterialsWithScoring;
    }
    catch(const std::invalid_argument& re)
    {
        // Phase 1 fix: previously silently `return *unfilteredCoreMaterials;` —
        // an exception here meant the entire filter became a no-op pass-through,
        // disabling the volumetric-losses cross-reference ranking with no
        // signal to the caller. That's a textbook silent shortcut. Per
        // policy: propagate. If the reference material's losses can't be
        // computed, the whole cross-reference can't proceed honestly.
        throw;
    }
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::MagneticCoreFilterResistivity::filter_core_materials(std::vector<std::pair<CoreMaterial, double>>* unfilteredCoreMaterials, CoreMaterial referenceCoreMaterial, double temperature, double weight) {
    if (weight <= 0) {
        return *unfilteredCoreMaterials;
    }
    std::vector<double> newScoring;

    double referenceResistivityWithTemperature = Core::get_resistivity(referenceCoreMaterial, temperature);
    add_scored_value("Reference", CoreMaterialCrossReferencerFilters::RESISTIVITY, referenceResistivityWithTemperature);


    for (size_t coreMaterialIndex = 0; coreMaterialIndex < (*unfilteredCoreMaterials).size(); ++coreMaterialIndex){
        CoreMaterial coreMaterial = (*unfilteredCoreMaterials)[coreMaterialIndex].first;

        double resistivityWithTemperature = Core::get_resistivity(coreMaterial, temperature);

        double scoring = fabs(referenceResistivityWithTemperature - resistivityWithTemperature);
        newScoring.push_back(scoring);
        add_scoring(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::RESISTIVITY, scoring);
        add_scored_value(coreMaterial.get_name(), CoreMaterialCrossReferencerFilters::RESISTIVITY, resistivityWithTemperature);
    }

    if ((*unfilteredCoreMaterials).size() != newScoring.size()) {
        throw CalculationException(ErrorCode::CALCULATION_ERROR, "Something wrong happened while filtering, size of unfilteredCoreMaterials: " + std::to_string((*unfilteredCoreMaterials).size()) + ", size of newScoring: " + std::to_string(newScoring.size()));
    }

    if ((*unfilteredCoreMaterials).size() > 0) {
        normalize_scoring(unfilteredCoreMaterials, &newScoring, weight, (*_filterConfiguration)[CoreMaterialCrossReferencerFilters::RESISTIVITY]);
    }
    return *unfilteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::get_cross_referenced_core_material(CoreMaterial referenceCoreMaterial, double temperature, size_t maximumNumberResults) {
    return get_cross_referenced_core_material(referenceCoreMaterial, temperature, _weights, maximumNumberResults);
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::get_cross_referenced_core_material(CoreMaterial referenceCoreMaterial, double temperature, std::map<CoreMaterialCrossReferencerFilters, double> weights, size_t maximumNumberResults) {
    auto defaults = Defaults();
    _weights = weights;

    std::vector<std::pair<CoreMaterial, double>> coreMaterials;

    if (coreMaterialDatabase.empty()) {
        load_core_materials();
    }

    for (const auto& [name, coreMaterial] : coreMaterialDatabase){
        if (name != referenceCoreMaterial.get_name()) {
            if (!_onlyManufacturer || coreMaterial.get_manufacturer_info().get_name() == _onlyManufacturer.value()) {
                coreMaterials.push_back({coreMaterial, 0.0});
            }
        }
    }

    auto filteredCoreMaterials = apply_filters(&coreMaterials, referenceCoreMaterial, temperature, weights, maximumNumberResults);
    return filteredCoreMaterials;
}

std::vector<std::pair<CoreMaterial, double>> CoreMaterialCrossReferencer::apply_filters(std::vector<std::pair<CoreMaterial, double>>* coreMaterials, CoreMaterial referenceCoreMaterial, double temperature, std::map<CoreMaterialCrossReferencerFilters, double> weights, size_t maximumNumberResults) {
    MagneticCoreFilterInitialPermeability filterInitialPermeability;
    MagneticCoreFilterRemanence filterRemanence;
    MagneticCoreFilterCoerciveForce filterCoerciveForce;
    MagneticCoreFilterSaturation filterSaturation;
    MagneticCoreFilterCurieTemperature filterCurieTemperature;

    json modelJson = _models["coreLosses"];
    OpenMagnetics::CoreLossesModels coreLossesModel;
    from_json(modelJson, coreLossesModel);
    MagneticCoreFilterVolumetricLosses filterVolumetricLosses(coreLossesModel);
    MagneticCoreFilterResistivity filterResistivity;

    wire_cross_referencer_filters<CoreMaterialCrossReferencerFilters>(
        {&filterInitialPermeability, &filterRemanence, &filterCoerciveForce, &filterSaturation,
         &filterCurieTemperature, &filterVolumetricLosses, &filterResistivity},
        &_scorings, /*validScorings=*/nullptr, &_scoredValues, &_filterConfiguration);

    std::vector<std::pair<CoreMaterial, double>> rankedCoreMaterials;

    auto referenceCoreMaterialApplication = Core::guess_material_application(referenceCoreMaterial);
    auto referenceCoreMaterialType = referenceCoreMaterial.get_material();
    for (const auto& [coreMaterial, scoring] : (*coreMaterials)) {
        auto coreMaterialType = coreMaterial.get_material();
        auto coreMaterialApplication = Core::guess_material_application(coreMaterial);

        if (coreMaterialApplication == referenceCoreMaterialApplication && (settings.get_core_cross_referencer_allow_different_core_material_type() || coreMaterialType == referenceCoreMaterialType)) {
            rankedCoreMaterials.push_back({coreMaterial, scoring}); 
        }
    }

    using F = CoreMaterialCrossReferencerFilters;
    std::vector<CrossReferencerStep<CoreMaterial, F>> steps = {
        {F::INITIAL_PERMEABILITY,
         [&](auto* r){ return filterInitialPermeability.filter_core_materials(r, referenceCoreMaterial, temperature, weights[F::INITIAL_PERMEABILITY]); }},
        {F::REMANENCE,
         [&](auto* r){ return filterRemanence.filter_core_materials(r, referenceCoreMaterial, temperature, weights[F::REMANENCE]); }},
        {F::COERCIVE_FORCE,
         [&](auto* r){ return filterCoerciveForce.filter_core_materials(r, referenceCoreMaterial, temperature, weights[F::COERCIVE_FORCE]); }},
        {F::SATURATION,
         [&](auto* r){ return filterSaturation.filter_core_materials(r, referenceCoreMaterial, temperature, weights[F::SATURATION]); }},
        {F::CURIE_TEMPERATURE,
         [&](auto* r){ return filterCurieTemperature.filter_core_materials(r, referenceCoreMaterial, temperature, weights[F::CURIE_TEMPERATURE]); }},
        {F::VOLUMETRIC_LOSSES,
         [&](auto* r){ return filterVolumetricLosses.filter_core_materials(r, referenceCoreMaterial, temperature, _models, weights[F::VOLUMETRIC_LOSSES]); }},
        {F::RESISTIVITY,
         [&](auto* r){ return filterResistivity.filter_core_materials(r, referenceCoreMaterial, temperature, weights[F::RESISTIVITY]); }},
    };
    rankedCoreMaterials = run_cross_referencer_pipeline(rankedCoreMaterials, steps, "Core Material Cross Referencer");

    if (rankedCoreMaterials.size() > maximumNumberResults) {
        rankedCoreMaterials.resize(maximumNumberResults);
    }
    return rankedCoreMaterials;
}

} // namespace OpenMagnetics
