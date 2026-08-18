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

    size_t size() const {
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

    std::vector<std::string> references() const {
        std::vector<std::string> filteredReferences;
        filteredReferences.reserve(_cache.size());
        for (const auto& [reference, value] : _cache) {
            filteredReferences.push_back(reference);
        }

        return filteredReferences;
    }

    std::vector<T> read() const {
        std::vector<T> filteredValues;
        filteredValues.reserve(_cache.size());
        for (const auto& [reference, value] : _cache) {
            filteredValues.push_back(value);
        }
        return filteredValues;
    }

    // Reads go through find(), never operator[]. operator[] is a MUTATING
    // member (it default-inserts a missing key), so concurrent readers of a
    // shared cache — which magneticsCache now is (ABT #817) — would be a data
    // race even when every key happens to be present.
    T read(std::string reference) const {
        auto entry = _cache.find(reference);
        if (entry == _cache.end()) {
            throw std::runtime_error("No value found with reference: " + reference);
        }
        return entry->second;
    }

    std::vector<T> read(std::vector<std::string> references) const {
        std::vector<T> filteredValues;
        filteredValues.reserve(references.size());
        for (const auto& reference : references) {
            auto entry = _cache.find(reference);
            if (entry != _cache.end()) {
                filteredValues.push_back(entry->second);
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

// The loaded part catalogue. Its storage (`_cache`, inherited) is SHARED
// between threads — see the magneticsCache declaration in Utils.h — because a
// catalogue of fully expanded Magnetics is the largest structure in the
// process (~1.9 GB for 5130 parts) and every thread wants the same read-only
// copy, exactly like coreDatabase and friends.
//
// The energy cache is different in kind: it is not the catalogue, it is a memo
// DERIVED from one operating point (compute_energy_cache clears and refills it
// per call). Two threads advising different designs would overwrite each
// other's energies, so it stays thread_local — shared storage, per-thread
// derived state (ABT #817).
class MagneticsCache : public Cache<OpenMagnetics::Magnetic> {
private:
    static thread_local std::map<std::string, double> _magneticEnergyCache;

    // Mutating the shared storage while other threads read it is undefined
    // behaviour, so the freeze flag that guards the other catalogues guards
    // this one too: a load or clear inside a frozen (parallel) region is a
    // loud failure rather than a race. Declared here rather than included from
    // Utils.h, which includes this header.
    static void throw_if_frozen(const std::string& operation);
public:
    void load(std::string reference, OpenMagnetics::Magnetic value) {
        throw_if_frozen("load");
        Cache<OpenMagnetics::Magnetic>::load(std::move(reference), std::move(value));
    }
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
    void clear();
    size_t energy_cache_size();
    std::map<std::string, double> read_magnetic_energy_cache();
    double read_magnetic_energy_cache(std::string reference);
    std::pair<std::string, double> get_maximum_magnetic_energy_in_cache();
    void compute_energy_cache(std::optional<OperatingPoint> operatingPoint = std::nullopt, bool saturationProportion = true);
    void compute_energy_cache(double temperature, std::optional<double> frequency = std::nullopt, bool saturationProportion = true);
    std::vector<std::string> filter_magnetics_by_energy(double minimumEnergy, std::optional<double> maximumEnergy = std::nullopt);
};

} // namespace OpenMagnetics