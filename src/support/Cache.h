#pragma once
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {

template<class T> class Cache {
private:
protected:
    std::map<std::string, T> _cache;
public:
    Cache() = default;
    ~Cache() = default;

    void clear() {
        _cache.clear();
    }

    size_t size(){
        return _cache.size();
    }

    std::vector<std::string> references() {
        std::vector<std::string> filteredReferences;
        for (auto [reference, value] : _cache) {
            filteredReferences.push_back(reference);
        }

        return filteredReferences;
    }

    std::vector<T> read() {
        std::vector<T> filteredValues;
        for (auto [reference, value] : _cache) {
            filteredValues.push_back(value);
        }
        return filteredValues;
    }

    T read(std::string reference) {
        if (_cache.contains(reference)) {
            return _cache[reference];
        } else {
            throw std::runtime_error("No value found with reference: " + reference);
        }
    }

    std::vector<T> read(std::vector<std::string> references) {
        std::vector<T> filteredValues;
        for (auto reference : references) {
            if (_cache.contains(reference)) {
                filteredValues.push_back(_cache[reference]);
            }
        }
        return filteredValues;
    }

    void load(std::string reference, T value) {
        _cache[reference] = value;
    }

    std::vector<T> evict(std::vector<std::string> references) {
        std::vector<T> evictedValues;
        for (auto reference : references) {
            if (_cache.contains(reference)) {
                evictedValues.push_back(std::move(_cache[reference]));
                _cache.erase(reference);
            }
        }
        return evictedValues;
    }

    T evict(std::string reference) {
        if (_cache.contains(reference)) {
            auto value = std::move(_cache[reference]);
            _cache.erase(reference);
            return value;
        } else {
            throw std::runtime_error("No value found with reference: " + reference);
        }
    }

};

class MagneticsCache : public Cache<OpenMagnetics::Magnetic> {
private:
    std::map<std::string, double> _magneticEnergyCache;
public:
    void clear();
    void autocomplete_magnetics();
    size_t energy_cache_size();
    std::pair<std::string, double> get_maximum_magnetic_energy_in_cache();
    void compute_energy_cache(std::optional<OperatingPoint> operatingPoint = std::nullopt);
    void compute_energy_cache(double temperature, std::optional<double> frequency = std::nullopt);
    std::vector<std::string> filter_magnetics_by_energy(double minimumEnergy, std::optional<double> maximumEnergy = std::nullopt);
};

class MasCache : public Cache<OpenMagnetics::Mas> {
private:
    std::map<std::string, double> _magneticEnergyCache;
public:
    void autocomplete_mas();
    size_t energy_cache_size();
    std::pair<std::string, double> get_maximum_magnetic_energy_in_cache();
    void compute_energy_cache(std::optional<OperatingPoint> operatingPoint = std::nullopt);
    void compute_energy_cache(double temperature, std::optional<double> frequency = std::nullopt);
    std::vector<std::string> filter_magnetics_by_energy(double minimumEnergy, std::optional<double> maximumEnergy = std::nullopt);
};

} // namespace OpenMagnetics