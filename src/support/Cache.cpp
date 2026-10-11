#include "support/Utils.h"
#include "support/Cache.h"
#include "physical_models/MagneticEnergy.h"
#include "Constants.h"
#include "support/Exceptions.h"


namespace OpenMagnetics {

// Per-thread derived memo over the shared catalogue storage (ABT #817).
thread_local std::map<std::string, double> MagneticsCache::_magneticEnergyCache;


void MagneticsCache::autocomplete_magnetics() {
    throw_if_frozen("autocomplete_magnetics");
    // ABT #1496: magnetic_autocomplete throws on a record whose winding MKF has to lay out (it states
    // no sections, layers or turns) and which does not fit its window.
    // A catalogue load must neither keep such a record half-wound nor drop it quietly, and one
    // exception per load would hide every unfit record after the first: autocomplete them all,
    // collect each failure, and fail the load once naming every one. The cache is only replaced
    // when every record autocompleted.
    std::map<std::string, Magnetic> autocompleted;
    std::vector<std::string> failures;
    std::optional<ErrorCode> commonCode;
    bool mixedCodes = false;
    for (const auto& [reference, magnetic] : _cache) {
        try {
            autocompleted.emplace(reference, magnetic_autocomplete(magnetic));
        }
        catch (const OpenMagneticsException& e) {
            failures.push_back("'" + reference + "': " + e.what());
            if (!commonCode) {
                commonCode = e.code();
            }
            else if (commonCode.value() != e.code()) {
                mixedCodes = true;
            }
        }
    }
    if (!failures.empty()) {
        std::string message = "magneticsCache::autocomplete_magnetics: " + std::to_string(failures.size()) + " of " +
                              std::to_string(_cache.size()) + " records could not be autocompleted:";
        for (const auto& failure : failures) {
            message += "\n  " + failure;
        }
        throw OpenMagneticsException(mixedCodes ? ErrorCode::INVALID_INPUT : commonCode.value(), message);
    }
    for (auto& [reference, magnetic] : autocompleted) {
        _cache[reference] = std::move(magnetic);
    }
}

void MagneticsCache::throw_if_frozen(const std::string& operation) {
    if (databases_frozen()) {
        throw InvalidInputException(
            ErrorCode::INVALID_INPUT,
            "magneticsCache::" + operation + " while the databases are frozen: the catalogue is "
            "shared between threads, so mutating it inside a parallel region is a data race. "
            "Load before the parallel region, or unfreeze first (ABT #817).");
    }
}

void MagneticsCache::clear() {
    throw_if_frozen("clear");
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
    if (_magneticEnergyCache.empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Energy cache is empty, call compute_energy_cache first");
    }
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

void MasCache::clear() {
    _cache.clear();
    _magneticEnergyCache.clear();
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
    if (_magneticEnergyCache.empty()) {
        throw InvalidInputException(ErrorCode::MISSING_DATA, "Energy cache is empty, call compute_energy_cache first");
    }
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
