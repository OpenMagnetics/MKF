#include "support/Utils.h"
#include "support/Cache.h"
#include "physical_models/MagneticEnergy.h"
#include "Constants.h"


namespace OpenMagnetics {

void Cache::clear_cache() {
    _magneticsCache.clear();
    _magneticEnergyCache.clear();
}

size_t Cache::get_size(){
    return _magneticsCache.size();
}
size_t Cache::get_energy_cache_size(){
    return _magneticEnergyCache.size();
}

std::vector<std::string> Cache::get_references() {
    std::vector<std::string> filteredReferences;
    for (auto [reference, magnetic] : _magneticsCache) {
        filteredReferences.push_back(reference);
    }

    return filteredReferences;
}

std::vector<OpenMagnetics::Magnetic> Cache::get_magnetics(std::optional<std::vector<std::string>> references) {
    std::vector<OpenMagnetics::Magnetic> filteredMagnetics;
    if (references) {
        for (auto [reference, magnetic] : _magneticsCache) {
            auto referencesAux = references.value();
            if (std::find(referencesAux.begin(), referencesAux.end(), reference) != referencesAux.end()) {
                filteredMagnetics.push_back(magnetic);
            }
        }
    }
    else {
        for (auto [reference, magnetic] : _magneticsCache) {
            filteredMagnetics.push_back(magnetic);
        }
    }
    return filteredMagnetics;
}

void Cache::load_magnetic(std::string reference, OpenMagnetics::Magnetic magnetic) {
    _magneticsCache[reference] = magnetic;
}

OpenMagnetics::Magnetic Cache::read_magnetic(std::string reference) {
    if (_magneticsCache.contains(reference)) {
        return _magneticsCache[reference];
    } else {
        throw std::runtime_error("No magnetic found with reference: " + reference);
    }
}

OpenMagnetics::Magnetic Cache::evict_magnetic(std::string reference) {
    if (_magneticsCache.contains(reference)) {
        auto magnetic = std::move(_magneticsCache[reference]);
        _magneticsCache.erase(reference);
        return magnetic;
    } else {
        throw std::runtime_error("No magnetic found with reference: " + reference);
    }
}

void Cache::autocomplete_magnetics() {
    for (auto [reference, magnetic] : _magneticsCache) {
        _magneticsCache[reference] = magnetic_autocomplete(magnetic);
    }
}

void Cache::remove_magnetic(std::string reference) {
    _magneticsCache.erase(reference);
}

void Cache::compute_energy_cache(std::optional<OperatingPoint> operatingPoint) {
    _magneticEnergyCache.clear();
    auto constants = Constants();
    double temperature = Defaults().ambientTemperature;
    if (operatingPoint) {
        temperature = operatingPoint->get_conditions().get_ambient_temperature();
    }
    if (operatingPoint) {
        auto frequency = operatingPoint->get_excitations_per_winding()[0].get_frequency();
        compute_energy_cache(temperature, frequency);
    }
    else {
        compute_energy_cache(temperature);
    }
}

void Cache::compute_energy_cache(double temperature, std::optional<double> frequency) {
    _magneticEnergyCache.clear();
    for (auto [reference, magnetic] : _magneticsCache) {
        double coreMaximumMagneticEnergy = MagneticEnergy().calculate_core_maximum_magnetic_energy(magnetic.get_core(), temperature, frequency);
        _magneticEnergyCache[reference] = coreMaximumMagneticEnergy;
    }
}

std::vector<std::string> Cache::filter_magnetics_by_energy(double minimumEnergy, std::optional<double> maximumEnergy) {
    std::vector<std::string> filteredReferences;
    for (auto [reference, energy] : _magneticEnergyCache) {
        if (maximumEnergy? energy >= minimumEnergy && energy <= maximumEnergy.value() : energy >= minimumEnergy) {
            filteredReferences.push_back(reference);
        }
    }

    return filteredReferences;
}

} // namespace OpenMagnetics
