// Test-only shim: compile Kirchhoff's in-process libngspice runner into MKF_tests.
//
// MKF consumes Kirchhoff as libKirchhoffApi.so, whose internals are hidden
// (-fvisibility=hidden + version script) — only the Kirchhoff::api string/JSON facade is
// exported, and that facade has no raw-deck entry point (api::simulate_ngspice takes a
// TAS, not a SPICE deck). Tests that must run an EMITTED MAGNETIC SUBCIRCUIT through the
// integrated engine (e.g. the mutual-resistance fit/emit verification, abt #76) therefore
// compile Kirchhoff::run_ngspice_in_process from source, against the same libngspice the
// Kirchhoff build links. This keeps every SPICE simulation on the in-process libngspice
// engine (never the external ngspice CLI, per the house rule).
//
// The include resolves through the Kirchhoff source include dir that the MKF target already
// propagates (KIRCHHOFF_SOURCE_DIR/src or the ExternalProject fetch location); the
// libngspice link for this TU is added to MKF_tests in CMakeLists.txt.
#define ENABLE_NGSPICE 1
#include "NgspiceRunner.cpp"
