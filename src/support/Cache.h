#pragma once
#include <atomic>
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

    // The entries for exactly `references`, keyed as in the cache -- e.g. to run the catalogue
    // adviser over one family or one admin-chosen pool of a cache that holds several. A
    // reference the cache does not hold THROWS, listing every missing one: silently searching
    // fewer parts than asked for is how a failed or partial load goes unnoticed.
    std::map<std::string, T> subset(const std::vector<std::string>& references) const {
        std::map<std::string, T> entries;
        std::vector<std::string> missing;
        for (const auto& reference : references) {
            auto entry = _cache.find(reference);
            if (entry == _cache.end()) {
                missing.push_back(reference);
                continue;
            }
            entries.emplace(entry->first, entry->second);
        }
        if (!missing.empty()) {
            std::string list;
            for (size_t index = 0; index < missing.size() && index < 10; ++index) {
                list += (index ? ", " : "") + missing[index];
            }
            throw std::runtime_error(std::to_string(missing.size()) + " of " + std::to_string(references.size()) + " requested references are not in the cache (" + list + (missing.size() > 10 ? ", ..." : "") + ")");
        }
        return entries;
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
    // ABT #1626: the energy memo is per thread, the storage it was computed from is shared, so a
    // clear() or load() on one thread could not empty the memo of another, and a load() or evict()
    // never touched any memo: the memo then answered for a catalogue that no longer existed. Every
    // mutation of the storage bumps _contentVersion; the memo records the version (and the material
    // catalogue generation it read the cores' materials from) it was computed for, and a read of a
    // memo computed for another state throws.
    inline static std::atomic<uint64_t> _contentVersion{0};
    static thread_local uint64_t _energyCacheContentVersion;
    static thread_local uint64_t _energyCacheCatalogueGeneration;
    static void bump_content_version() { _contentVersion.fetch_add(1, std::memory_order_relaxed); }
    static void throw_if_energy_cache_stale();
    static void stamp_energy_cache();

    // Mutating the shared storage while other threads read it is undefined
    // behaviour, so the freeze flag that guards the other catalogues guards
    // this one too: a load or clear inside a frozen (parallel) region is a
    // loud failure rather than a race. Declared here rather than included from
    // Utils.h, which includes this header.
    static void throw_if_frozen(const std::string& operation);
public:
    void load(std::string reference, OpenMagnetics::Magnetic value) {
        throw_if_frozen("load");
        bump_content_version();
        Cache<OpenMagnetics::Magnetic>::load(std::move(reference), std::move(value));
    }
    std::vector<OpenMagnetics::Magnetic> evict(std::vector<std::string> references) {
        throw_if_frozen("evict");
        bump_content_version();
        return Cache<OpenMagnetics::Magnetic>::evict(std::move(references));
    }
    OpenMagnetics::Magnetic evict(std::string reference) {
        throw_if_frozen("evict");
        bump_content_version();
        return Cache<OpenMagnetics::Magnetic>::evict(std::move(reference));
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