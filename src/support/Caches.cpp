#include "support/Utils.h"
#include "support/Caches.h"
#include "physical_models/MagneticEnergy.h"
#include "Constants.h"

std::map<std::string, OpenMagnetics::Magnetic> magneticsCache;
std::map<std::string, double> magneticEnergyCache;


namespace OpenMagnetics {

void clear_magnetic_cache() {
    magneticsCache.clear();
    magneticEnergyCache.clear();
}

std::vector<std::string> get_magnetic_cache_references() {
    std::vector<std::string> filteredReferences;
    for (auto [reference, magnetic] : magneticsCache) {
        filteredReferences.push_back(reference);
    }

    return filteredReferences;
}

std::vector<OpenMagnetics::Magnetic> get_magnetics_from_cache(std::optional<std::vector<std::string>> references) {
    std::vector<OpenMagnetics::Magnetic> filteredMagnetics;
    if (references) {
        for (auto [reference, magnetic] : magneticsCache) {
            auto referencesAux = references.value();
            if (std::find(referencesAux.begin(), referencesAux.end(), reference) != referencesAux.end()) {
                filteredMagnetics.push_back(magnetic);
            }
        }
    }
    else {
        for (auto [reference, magnetic] : magneticsCache) {
            filteredMagnetics.push_back(magnetic);
        }
    }
    return filteredMagnetics;
}

void load_magnetic_in_cache(std::string reference, OpenMagnetics::Magnetic magnetic) {
    magneticsCache[reference] = magnetic;
}

OpenMagnetics::Magnetic read_magnetic_from_cache(std::string reference) {
    if (magneticsCache.contains(reference)) {
        return magneticsCache[reference];
    } else {
        throw std::runtime_error("No magnetic found with reference: " + reference);
    }
}

OpenMagnetics::Magnetic evict_magnetic_from_cache(std::string reference) {
    if (magneticsCache.contains(reference)) {
        auto magnetic = std::move(magneticsCache[reference]);
        magneticsCache.erase(reference);
        return magnetic;
    } else {
        throw std::runtime_error("No magnetic found with reference: " + reference);
    }
}

void autocomplete_magnetics_in_cache() {
    for (auto [reference, magnetic] : magneticsCache) {
        magneticsCache[reference] = magnetic_autocomplete(magnetic);
    }
}

void remove_magnetic_from_cache(std::string reference) {
    magneticsCache.erase(reference);
}

void compute_energy_cache(std::optional<OperatingPoint> operatingPoint) {
    magneticEnergyCache.clear();
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

void compute_energy_cache(double temperature, std::optional<double> frequency) {
    magneticEnergyCache.clear();
    for (auto [reference, magnetic] : magneticsCache) {
        double coreMaximumMagneticEnergy = MagneticEnergy().calculate_core_maximum_magnetic_energy(magnetic.get_core(), temperature, frequency);
        magneticEnergyCache[reference] = coreMaximumMagneticEnergy;
    }
}

std::vector<std::string> filter_magnetics_by_energy(double minimumEnergy, std::optional<double> maximumEnergy) {
    std::vector<std::string> filteredReferences;
    for (auto [reference, energy] : magneticEnergyCache) {
        if (maximumEnergy? energy >= minimumEnergy && energy <= maximumEnergy.value() : energy >= minimumEnergy) {
            filteredReferences.push_back(reference);
        }
    }

    return filteredReferences;
}

} // namespace OpenMagnetics
