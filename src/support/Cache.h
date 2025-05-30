#pragma once
#include "constructive_models/Magnetic.h"

using namespace MAS;

namespace OpenMagnetics {

class Cache {
private:
	std::map<std::string, OpenMagnetics::Magnetic> _magneticsCache;
	std::map<std::string, double> _magneticEnergyCache;
protected:
public:
    Cache() = default;
    ~Cache() = default;
	void clear_cache();
	size_t get_size();
	size_t get_energy_cache_size();
	std::vector<std::string> get_references();
	std::vector<OpenMagnetics::Magnetic> get_magnetics(std::optional<std::vector<std::string>> references);
	void autocomplete_magnetics();
	OpenMagnetics::Magnetic read_magnetic(std::string reference);
	OpenMagnetics::Magnetic evict_magnetic(std::string reference);
	void load_magnetic(std::string reference, OpenMagnetics::Magnetic magnetic);
	void remove_magnetic(std::string reference);
	void compute_energy_cache(std::optional<OperatingPoint> operatingPoint = std::nullopt);
	void compute_energy_cache(double temperature, std::optional<double> frequency = std::nullopt);
	std::vector<std::string> filter_magnetics_by_energy(double minimumEnergy, std::optional<double> maximumEnergy = std::nullopt);
};

} // namespace OpenMagnetics