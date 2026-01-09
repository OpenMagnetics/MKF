#include "support/Utils.h"
#include "support/Cache.h"
#include "physical_models/MagneticEnergy.h"
#include "Constants.h"
#include "support/Exceptions.h"


namespace OpenMagnetics {


void MagneticsCache::autocomplete_magnetics() {
    for (auto [reference, magnetic] : _cache) {
        _cache[reference] = magnetic_autocomplete(magnetic);
    }
}

void MagneticsCache::clear() {
    _cache.clear();
    _magneticEnergyCache.clear();
}

size_t MagneticsCache::energy_cache_size(){
    return _magneticEnergyCache.size();
}

void MagneticsCache::compute_energy_cache(std::optional<OperatingPoint> operatingPoint, bool saturationProportion) {
    _magneticEnergyCache.clear();
    auto constants = Constants();
    double temperature = Defaults().ambientTemperature;
    if (operatingPoint) {
        temperature = operatingPoint->get_conditions().get_ambient_temperature();
    }
    if (operatingPoint) {
        auto frequency = operatingPoint->get_excitations_per_winding()[0].get_frequency();
        compute_energy_cache(temperature, frequency, saturationProportion);
    }
    else {
        compute_energy_cache(temperature, std::nullopt, saturationProportion);
    }
}

void MagneticsCache::compute_energy_cache(double temperature, std::optional<double> frequency, bool saturationProportion) {
    _magneticEnergyCache.clear();
    for (auto [reference, magnetic] : _cache) {
        double coreMaximumMagneticEnergy = MagneticEnergy().calculate_core_maximum_magnetic_energy(magnetic.get_core(), temperature, frequency, saturationProportion);
        _magneticEnergyCache[reference] = coreMaximumMagneticEnergy;
    }
}

std::pair<std::string, double> MagneticsCache::get_maximum_magnetic_energy_in_cache() {
    return *std::max_element(_magneticEnergyCache.begin(), _magneticEnergyCache.end(), [](const std::pair<std::string, double>& p1, const std::pair<std::string, double>& p2) {return p1.second < p2.second; });
}

std::vector<std::string> MagneticsCache::filter_magnetics_by_energy(double minimumEnergy, std::optional<double> maximumEnergy) {
    std::vector<std::string> filteredReferences;
    for (auto [reference, energy] : _magneticEnergyCache) {
        if (maximumEnergy? energy >= minimumEnergy && energy <= maximumEnergy.value() : energy >= minimumEnergy) {
            filteredReferences.push_back(reference);
        }
    }

    return filteredReferences;
}

std::map<std::string, double> MagneticsCache::read_magnetic_energy_cache() {
    return _magneticEnergyCache;
}

double MagneticsCache::read_magnetic_energy_cache(std::string reference) {
    if (_magneticEnergyCache.contains(reference)) {
        return _magneticEnergyCache[reference];
    } else {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "No value found with reference: " + reference);
    }
}

size_t MasCache::energy_cache_size(){
    return _magneticEnergyCache.size();
}

void MasCache::compute_energy_cache(std::optional<OperatingPoint> operatingPoint, bool saturationProportion) {
    _magneticEnergyCache.clear();
    auto constants = Constants();
    double temperature = Defaults().ambientTemperature;
    if (operatingPoint) {
        temperature = operatingPoint->get_conditions().get_ambient_temperature();
    }
    if (operatingPoint) {
        auto frequency = operatingPoint->get_excitations_per_winding()[0].get_frequency();
        compute_energy_cache(temperature, frequency, saturationProportion);
    }
    else {
        compute_energy_cache(temperature, std::nullopt, saturationProportion);
    }
}

void MasCache::compute_energy_cache(double temperature, std::optional<double> frequency, bool saturationProportion) {
    _magneticEnergyCache.clear();
    for (auto [reference, mas] : _cache) {
        double coreMaximumMagneticEnergy = MagneticEnergy().calculate_core_maximum_magnetic_energy(mas.get_magnetic().get_core(), temperature, frequency, saturationProportion);
        _magneticEnergyCache[reference] = coreMaximumMagneticEnergy;
    }
}

std::vector<std::string> MasCache::filter_magnetics_by_energy(double minimumEnergy, std::optional<double> maximumEnergy) {
    std::vector<std::string> filteredReferences;
    for (auto [reference, energy] : _magneticEnergyCache) {
        if (maximumEnergy? energy >= minimumEnergy && energy <= maximumEnergy.value() : energy >= minimumEnergy) {
            filteredReferences.push_back(reference);
        }
    }

    return filteredReferences;
}

std::pair<std::string, double> MasCache::get_maximum_magnetic_energy_in_cache() {
    return *std::max_element(_magneticEnergyCache.begin(), _magneticEnergyCache.end(), [](const std::pair<std::string, double>& p1, const std::pair<std::string, double>& p2) {return p1.second < p2.second; });
}

std::map<std::string, double> MasCache::read_magnetic_energy_cache() {
    return _magneticEnergyCache;
}
double MasCache::read_magnetic_energy_cache(std::string reference) {
    if (_magneticEnergyCache.contains(reference)) {
        return _magneticEnergyCache[reference];
    } else {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "No value found with reference: " + reference);
    }
}

} // namespace OpenMagnetics
