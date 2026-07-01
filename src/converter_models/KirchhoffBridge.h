#pragma once

// KirchhoffBridge — MKF's thin adapter over Kirchhoff (KH), the exportable converter-model library that
// replaces MKF's own converter_models. KH is linked as libKirchhoffApi.so: a shared library that exposes
// ONLY the Kirchhoff::api string/JSON facade with its internals (including KH's own `namespace MAS` types)
// hidden. That is deliberate — MKF ALSO generates its types into `namespace MAS`, so this TU includes only
// KirchhoffApi.hpp (plain std::string signatures) and never sees KH's MAS types: no symbol collision.
//
// KH also carries the ngspice simulator, which is why MKF no longer builds its own: any MKF path that needs
// a circuit solved goes through Kirchhoff::api::simulate_ngspice (see simulate_tas below).

#include <string>
#include <nlohmann/json.hpp>

#include "processors/Inputs.h"   // OpenMagnetics::Inputs (MKF's own MAS-derived type)

namespace OpenMagnetics {

// Design a converter from a high-level spec + topology name, returning the whole assembled TAS (JSON).
// `topology` is the lowercase KH name (flyback, buck, llc, …). Throws on a KH-side error.
nlohmann::json design_converter_tas(const std::string& topology, const nlohmann::json& spec);

// The main magnetic's Inputs — the DesignRequirements + operating points the MagneticAdviser designs
// around — extracted from a KH TAS and parsed into MKF's own Inputs (same MAS JSON schema).
Inputs design_requirements_from_tas(const nlohmann::json& tas);

// Run a TAS through KH's in-process ngspice, returning the per-vector summary JSON
// ({success,error,tStart,tEnd,points,vectors:{name:{average,min,max,last}}}). This is MKF's simulator now.
nlohmann::json simulate_tas(const nlohmann::json& tas, const nlohmann::json& fidelity);

}  // namespace OpenMagnetics
