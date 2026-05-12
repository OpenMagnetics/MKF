#pragma once
#include "constructive_models/Magnetic.h"
#include "support/Exceptions.h"

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

    // Returns a const reference to the underlying map. Previously this
    // returned the map by value, which copied every cached element on
    // every call — `magneticsCache.get()` is invoked in the catalog
    // adviser hot path with hundreds of fully-expanded Magnetics. Callers
    // that need a mutable copy can do `auto m = cache.get();` explicitly.
    const std::map<std::string, T>& get() const {
        return _cache;
    }

    std::vector<std::string> references() {
        std::vector<std::string> filteredReferences;
        filteredReferences.reserve(_cache.size());
        for (const auto& [reference, value] : _cache) {
            filteredReferences.push_back(reference);
        }

        return filteredReferences;
    }

    std::vector<T> read() {
        std::vector<T> filteredValues;
        filteredValues.reserve(_cache.size());
        for (const auto& [reference, value] : _cache) {
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
        filteredValues.reserve(references.size());
        for (auto reference : references) {
            if (_cache.contains(reference)) {
                filteredValues.push_back(_cache[reference]);
            }
        }
        return filteredValues;
    }

    // Move both the key and the value into the map. Previously this took
    // `T value` by value and then assigned by copy, costing two whole
    // Magnetic copies per insert.
    void load(std::string reference, T value) {
        _cache[std::move(reference)] = std::move(value);
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
    std::map<std::string, double> read_magnetic_energy_cache();
    double read_magnetic_energy_cache(std::string reference);
    std::pair<std::string, double> get_maximum_magnetic_energy_in_cache();
    void compute_energy_cache(std::optional<OperatingPoint> operatingPoint = std::nullopt, bool saturationProportion = true);
    void compute_energy_cache(double temperature, std::optional<double> frequency = std::nullopt, bool saturationProportion = true);
    std::vector<std::string> filter_magnetics_by_energy(double minimumEnergy, std::optional<double> maximumEnergy = std::nullopt);
};

class MasCache : public Cache<OpenMagnetics::Mas> {
private:
    std::map<std::string, double> _magneticEnergyCache;
public:
    void autocomplete_mas();
    size_t energy_cache_size();
    std::map<std::string, double> read_magnetic_energy_cache();
    double read_magnetic_energy_cache(std::string reference);
    std::pair<std::string, double> get_maximum_magnetic_energy_in_cache();
    void compute_energy_cache(std::optional<OperatingPoint> operatingPoint = std::nullopt, bool saturationProportion = true);
    void compute_energy_cache(double temperature, std::optional<double> frequency = std::nullopt, bool saturationProportion = true);
    std::vector<std::string> filter_magnetics_by_energy(double minimumEnergy, std::optional<double> maximumEnergy = std::nullopt);
};

} // namespace OpenMagnetics