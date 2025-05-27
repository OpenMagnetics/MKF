#pragma once
#include "constructive_models/Magnetic.h"


extern std::map<std::string, OpenMagnetics::Magnetic> magneticsCache;
extern std::map<std::string, double> magneticEnergyCache;

using namespace MAS;

namespace OpenMagnetics {

void clear_magnetic_cache();
std::vector<std::string> get_magnetic_cache_references();
std::vector<OpenMagnetics::Magnetic> get_magnetics_from_cache(std::optional<std::vector<std::string>> references);
void autocomplete_magnetics_in_cache();
OpenMagnetics::Magnetic read_magnetic_from_cache(std::string reference);
OpenMagnetics::Magnetic evict_magnetic_from_cache(std::string reference);
void load_magnetic_in_cache(std::string reference, OpenMagnetics::Magnetic magnetic);
void remove_magnetic_from_cache(std::string reference);
void compute_energy_cache(std::optional<OperatingPoint> operatingPoint = std::nullopt);
void compute_energy_cache(double temperature, std::optional<double> frequency = std::nullopt);
std::vector<std::string> filter_magnetics_by_energy(double minimumEnergy, std::optional<double> maximumEnergy = std::nullopt);

} // namespace OpenMagnetics