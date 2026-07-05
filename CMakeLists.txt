cmake_minimum_required(VERSION 3.18.0)

if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("MKF_Project" VERSION 2.0.0 LANGUAGES CXX)

# ============================================================================
# C++ Standard and Compiler Settings
# ============================================================================
set(CMAKE_CXX_STANDARD 23) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Let user configure build type; default to Release if not set
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox")
    # Use /W3 instead of /W0 to catch potential issues; suppress only specific warnings if needed
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4996") # Suppress deprecation warnings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267") # Suppress size_t to int conversion warnings
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /FS")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /FS")
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-deprecated-declarations -Wno-unused-parameter -Wno-switch -Wno-format-truncation")
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    # Maximum-speed optimized build. NOTE: deliberately NO -ffast-math — this is a numeric/physics
    # codebase whose results (and pinned characterization tests) depend on IEEE FP semantics, NaN/Inf
    # handling and operation ordering, all of which -ffast-math breaks. -march=native tunes for the
    # build machine's CPU (binaries are not portable to older CPUs); -funroll-loops + -O3 + NDEBUG.
    # -march=native is x86/native-only; it is invalid for the wasm32 (Emscripten) target.
    if(EMSCRIPTEN)
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -funroll-loops")
        set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -funroll-loops")
    else()
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -funroll-loops")
        set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -funroll-loops")
    endif()
endif()


# ============================================================================
# Build Options
# ============================================================================
option(EMBED_MAS_DATA "Embed data in executable" ON)
option(EMBED_MAS_CORES "Embed cores in executable" ON)
option(EMBED_MAS_CORE_SHAPES "Embed core shapes in executable" ON)
option(EMBED_MAS_CORE_MATERIALS "Embed core materials in executable" ON)
option(EMBED_MAS_WIRES "Embed wires in executable" ON)
option(EMBED_MAS_WIRE_MATERIALS "Embed wire materials in executable" ON)
option(EMBED_MAS_BOBBINS "Embed bobbins in executable" ON)
option(EMBED_MAS_INSULATION_MATERIALS "Embed insulation materials in executable" ON)
option(INCLUDE_MKF_TESTS "Compile the test suite" ON)
option(BUILD_TESTS "Build tests (for dependencies)" OFF)
option(BUILD_EXAMPLES "Build examples (for dependencies)" OFF)
option(BUILD_DEMO "Build demo (for dependencies)" OFF)
# REMOVED: # REMOVED: option(HAVE_LAPACK "Use LAPACK" OFF) — no longer needed, levmar replaced by Eigen — no longer needed, levmar replaced by Eigen
option(ENABLE_NGSPICE "Enable ngspice integration for circuit simulation" OFF)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "ngspice integration: ${ENABLE_NGSPICE}")
message(STATUS "Embedding MAS Data: ${EMBED_MAS_DATA}")
message(STATUS "Embedding MAS Cores: ${EMBED_MAS_CORES}")
message(STATUS "Embedding MAS CoreShapes: ${EMBED_MAS_CORE_SHAPES}")
message(STATUS "Embedding MAS CoreMaterials: ${EMBED_MAS_CORE_MATERIALS}")
message(STATUS "Embedding MAS Wires: ${EMBED_MAS_WIRES}")
message(STATUS "Embedding MAS WireMaterials: ${EMBED_MAS_WIRE_MATERIALS}")
message(STATUS "Embedding MAS Bobbins: ${EMBED_MAS_BOBBINS}")
message(STATUS "Embedding MAS InsulationMaterials: ${EMBED_MAS_INSULATION_MATERIALS}")
message(STATUS "Including MKF Tests: ${INCLUDE_MKF_TESTS}")


# ============================================================================
# Dependencies via FetchContent
# ============================================================================
include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

# CMakeRC - Resource compiler (use FetchContent instead of file download)
FetchContent_Declare(cmrc
    GIT_REPOSITORY https://github.com/vector-of-bool/cmrc.git
    GIT_TAG 952ffddba731fc110bd50409e8d2b8a06abbd237)  # was 'master' (floating); pinned to known-good
FetchContent_MakeAvailable(cmrc)
include("${cmrc_SOURCE_DIR}/CMakeRC.cmake")

# JSON library
FetchContent_Declare(json
    GIT_REPOSITORY https://github.com/nlohmann/json/
    GIT_TAG  tags/v3.11.3)
FetchContent_MakeAvailable(json)
include_directories("${CMAKE_BINARY_DIR}/_deps/json-src/include/nlohmann/")
include_directories("${CMAKE_BINARY_DIR}/_deps/json-src/include/")


FetchContent_Declare(magic-enum
    GIT_REPOSITORY https://github.com/Neargye/magic_enum
    GIT_TAG  tags/v0.9.6)
FetchContent_MakeAvailable(magic-enum)
include_directories("${CMAKE_BINARY_DIR}/_deps/magic-enum-src/include/magic_enum")


# ============================================================================
# Kirchhoff — the exportable converter-model library that replaces MKF's converter_models (and carries the
# ngspice simulator, which is why MKF no longer builds its own: MKF's ENABLE_NGSPICE stays OFF by default
# and any circuit that must be solved goes through Kirchhoff::api::simulate_ngspice via KirchhoffBridge).
#
# Built ISOLATED with ExternalProject so KH's own dependency targets — critically its `namespace MAS`
# types, plus CAS/SAS/RAS/… — are configured in a SEPARATE CMake context and can never collide with MKF's
# identically-named ones. KH is consumed as libKirchhoffApi.so: a shared lib that exports ONLY the
# Kirchhoff::api string/JSON facade (KH's MAS:: symbols are hidden), so nothing but JSON strings crosses the
# boundary. No MAS material DB is pulled (KH has none; MKF's EMBED_MAS_* is independent and only for HS).
# ============================================================================
# Default ON natively, but OFF under Emscripten: the KirchhoffBridge links libKirchhoffApi.so (a native
# hidden-visibility shared lib) which has no wasm equivalent, and in the browser the converter surface is a
# SEPARATE module (webKirchhoff) — webMKF is magnetics-only. So the wasm MKF must not pull Kirchhoff.
if(EMSCRIPTEN)
    set(_mkf_use_kirchhoff_default OFF)
else()
    set(_mkf_use_kirchhoff_default ON)
endif()
option(MKF_USE_KIRCHHOFF "Link the Kirchhoff converter-model library (libKirchhoffApi.so)" ${_mkf_use_kirchhoff_default})
if(MKF_USE_KIRCHHOFF)
    include(ExternalProject)
    # Kirchhoff is configured in its OWN CMake context (isolated from MKF's targets). Two ways to provide it:
    #   * a GitHub clone (DEFAULT — GIT_REPOSITORY/GIT_TAG below; recursive, so KH's submodules come too), or
    #   * a local working tree by setting KIRCHHOFF_SOURCE_DIR to a path.
    # The github-fetch default became usable once Kirchhoff's `AAS` submodule was published on GitHub
    # (its siblings PEAS/RAS/CAS/SAS/MAS/CIAS/TAS/CONAS/CTAS and deps/MKF all resolve). To develop against
    # a local KH checkout instead, pass -DKIRCHHOFF_SOURCE_DIR=/path/to/Kirchhoff.
    set(KIRCHHOFF_GIT_TAG "main" CACHE STRING "Kirchhoff git tag/branch to build (when fetching from GitHub)")
    set(KIRCHHOFF_SOURCE_DIR "" CACHE PATH "Local Kirchhoff working tree (leave empty to fetch from GitHub)")
    set(KIRCHHOFF_BUILD_DIR "${CMAKE_BINARY_DIR}/kirchhoff-build")
    set(KIRCHHOFF_API_LIB "${KIRCHHOFF_BUILD_DIR}/libKirchhoffApi.so")
    if(KIRCHHOFF_SOURCE_DIR)
        set(_kirchhoff_src SOURCE_DIR "${KIRCHHOFF_SOURCE_DIR}")
    else()
        set(_kirchhoff_src GIT_REPOSITORY "git@github.com:OpenConverters/Kirchhoff.git"
                           GIT_TAG "${KIRCHHOFF_GIT_TAG}"
                           # KH's deps/MKF submodule is the circular back-reference to THIS repo, and its
                           # recorded pin can be a commit that never reached origin ("not our ref"). KH only
                           # needs MKF for the shared WaveformProcessor.cpp, which we hand it directly via
                           # -DMKF_DIR below — so init every KH submodule EXCEPT deps/MKF.
                           GIT_SUBMODULES "deps/PEAS;deps/RAS;deps/CAS;deps/SAS;deps/MAS;deps/CIAS;deps/TAS;deps/AAS;deps/CONAS;deps/CTAS")
    endif()

    # Kelvin — KH's shared TAS component selector (consumed via src/KelvinBridge.cpp). Since ~mid-2026 KH
    # requires it (add_subdirectory(KELVIN_DIR)), but it is NOT a KH submodule, so MKF must supply the
    # source. KH builds it in-tree reusing KH's own PEAS, so we only need the checkout, not a separate
    # build. MKF itself uses NOTHING from Kelvin (it is behind the libKirchhoffApi.so facade); this only
    # exists to let KH configure. Point -DKIRCHHOFF_KELVIN_DIR at a local checkout, else it is fetched.
    set(KIRCHHOFF_KELVIN_DIR "" CACHE PATH "Local Kelvin checkout for the Kirchhoff build (empty = fetch)")
    if(NOT KIRCHHOFF_KELVIN_DIR)
        set(KIRCHHOFF_KELVIN_GIT_TAG "main" CACHE STRING "Kelvin git tag/branch (when fetched for Kirchhoff)")
        include(FetchContent)
        FetchContent_Declare(kirchhoff_kelvin
            GIT_REPOSITORY "git@github.com:OpenConverters/Kelvin.git"
            GIT_TAG "${KIRCHHOFF_KELVIN_GIT_TAG}"
            GIT_SUBMODULES ""   # KH overrides KELVIN_PEAS_DIR to its own PEAS; Kelvin's own submodules are unused
            SOURCE_DIR "${CMAKE_BINARY_DIR}/kirchhoff-kelvin-src")
        # Download only — KH does the add_subdirectory(KELVIN_DIR) in its own context.
        FetchContent_GetProperties(kirchhoff_kelvin)
        if(NOT kirchhoff_kelvin_POPULATED)
            FetchContent_Populate(kirchhoff_kelvin)
        endif()
        set(KIRCHHOFF_KELVIN_DIR "${kirchhoff_kelvin_SOURCE_DIR}")
    endif()

    ExternalProject_Add(kirchhoff_ext
        ${_kirchhoff_src}
        BINARY_DIR "${KIRCHHOFF_BUILD_DIR}"
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=Release
            -DKIRCHHOFF_BUILD_SHARED_API=ON   # emit libKirchhoffApi.so (string facade, hidden internals)
            -DKIRCHHOFF_BUILD_PYBIND=OFF
            -DENABLE_NGSPICE=ON               # KH carries the simulator MKF stopped building
            -DMKF_DIR=${CMAKE_CURRENT_SOURCE_DIR}          # this MKF checkout supplies KH's shared WaveformProcessor
            -DKELVIN_DIR=${KIRCHHOFF_KELVIN_DIR}           # Kelvin source (fetched above or a local checkout)
        BUILD_ALWAYS OFF
        INSTALL_COMMAND ""
        BUILD_BYPRODUCTS "${KIRCHHOFF_API_LIB}")
    add_library(KirchhoffApi SHARED IMPORTED)
    set_target_properties(KirchhoffApi PROPERTIES IMPORTED_LOCATION "${KIRCHHOFF_API_LIB}")
    add_dependencies(KirchhoffApi kirchhoff_ext)
endif()


# Make profiling and ASan opt-in so users can build fast Release binaries by default
option(ENABLE_PROFILING "Enable gprof profiling (-pg)" OFF)
option(ENABLE_ASAN "Enable AddressSanitizer in Debug builds" OFF)

# Link-time optimization for the last few % of runtime speed. OFF by default because LTO linking of
# a library this size is memory-hungry and can OOM on machines that need a low -j (e.g. -j2). Turn
# it on (-DENABLE_LTO=ON) when the build host has the RAM to spare.
option(ENABLE_LTO "Enable link-time optimization (-flto) in Release builds" OFF)
if(ENABLE_LTO AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_msg)
    if(_ipo_ok)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    else()
        message(WARNING "LTO requested but not supported: ${_ipo_msg}")
    endif()
endif()

if(ENABLE_PROFILING)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()

# Set flags for Debug builds to include AddressSanitizer and debug symbols (-g)
# These are applied only when ENABLE_ASAN is enabled or when explicitly requested
if(ENABLE_ASAN)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -fsanitize=address")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -fsanitize=address")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address")
    set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -fsanitize=address")
endif()

option(ENABLE_COVERAGE "Build with gcov coverage instrumentation" OFF)
if(ENABLE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(--coverage -O0 -g)
    add_link_options(--coverage)
endif()

FetchContent_Declare(spline
    GIT_REPOSITORY https://github.com/AlfVII/spline.git
    GIT_TAG 10b4ca336b8ad3b3b842d4a2548a94ea00355a5f)  # pin default-branch HEAD to known-good
FetchContent_MakeAvailable(spline)
include_directories("${CMAKE_BINARY_DIR}/_deps/spline-src/src")


FetchContent_Declare(svg
    GIT_REPOSITORY https://github.com/AlfVII/svg
    GIT_TAG b92693f84ba6bbe79ea1965b68fffbcd1b991186)  # pin default-branch HEAD to known-good
FetchContent_MakeAvailable(svg)
include_directories("${CMAKE_BINARY_DIR}/_deps/svg-src/src")

FetchContent_Declare(rapidfuzz
  GIT_REPOSITORY https://github.com/rapidfuzz/rapidfuzz-cpp.git
  GIT_TAG b5830af53bd1b3c7460a8de1e9f7095df99b3470)  # was 'main' (floating); pinned to known-good
FetchContent_MakeAvailable(rapidfuzz)

# Eigen - Header-only linear algebra library (for Albach 2D boundary value solver)
FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE)
# Disable Eigen's BLAS/LAPACK and testing (avoids Fortran dependency for WASM builds)
set(EIGEN_BUILD_BLAS OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_LAPACK OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(eigen)
include_directories("${CMAKE_BINARY_DIR}/_deps/eigen-src")


# ============================================================================
# ngspice - Circuit Simulator (optional)
# ============================================================================
if(ENABLE_NGSPICE)
    message(STATUS "Configuring ngspice...")
    
    set(NGSPICE_PREFIX "${CMAKE_BINARY_DIR}/_deps/ngspice")
    set(NGSPICE_INSTALL_DIR "${NGSPICE_PREFIX}/install")
    
    # Check if ngspice is already installed (e.g. from a previous build)
    set(NGSPICE_ALREADY_INSTALLED FALSE)
    if(EXISTS "${NGSPICE_INSTALL_DIR}/lib/libngspice.so" AND EXISTS "${NGSPICE_INSTALL_DIR}/include/ngspice/sharedspice.h")
        set(NGSPICE_ALREADY_INSTALLED TRUE)
        message(STATUS "ngspice already installed at ${NGSPICE_INSTALL_DIR}, skipping build")
    endif()
    
    if(NOT NGSPICE_ALREADY_INSTALLED)
        # ngspice can be built as a shared library with the sharedspice API
        # We use ExternalProject for more control over the build process
        include(ExternalProject)
        
        set(NGSPICE_VERSION "45.2")
        set(NGSPICE_URL "file:///home/alf/OpenMagnetics/ngspice-45.2.tar.gz")
        
        # Configure options for ngspice shared library build
    # Basic configuration - NO CIDER, NO XSPICE for minimal build
    # Key options:
    # --with-ngshared: Build as shared library with C API
    # --disable-debug: Disable debug for performance
    # --without-x: No X11 dependencies
    
    if(EMSCRIPTEN)
        # Special configuration for WASM/Emscripten build
        # Based on ngspice patch #96: https://sourceforge.net/p/ngspice/patches/96/
        # Using --with-ngshared produces a static lib (.so extension but actually WASM static lib)
        # that provides the sharedspice.h API for WASM targets
        
        set(NGSPICE_CONFIGURE_FLAGS
            "--with-ngshared"
            "--disable-debug"
            "--without-x"
            "--disable-xspice"
            "--disable-cider"
            "--disable-openmp"
            "--with-readline=no"
            "CFLAGS=-O2 -fwasm-exceptions -sSUPPORT_LONGJMP=wasm"
            "CXXFLAGS=-O2 -fwasm-exceptions -sSUPPORT_LONGJMP=wasm"
            "LDFLAGS=-fwasm-exceptions -sSUPPORT_LONGJMP=wasm"
        )
        
        # Create a patch script to fix WASM compatibility issues
        # Based on ngspice patch #96 for WebAssembly support
        set(NGSPICE_PATCH_SCRIPT "${NGSPICE_PREFIX}/apply_wasm_patches.sh")
        file(WRITE ${NGSPICE_PATCH_SCRIPT} "#!/bin/bash
set -e
cd \$1

# Fix configure to recognize .wasm output files
if [ -f configure ]; then
    sed -i 's/ac_files=\\\"a.out/ac_files=\\\"a.out a.out.js a.out.wasm/' configure
fi

# Fix main.c - remove main function when building shared module
# The SHARED_MODULE macro should already handle this in newer ngspice versions
if [ -f src/main.c ]; then
    # Check if SHARED_MODULE guard already exists around main
    if ! grep -q '#ifndef SHARED_MODULE' src/main.c 2>/dev/null || ! grep -A1 '#ifndef SHARED_MODULE' src/main.c | grep -q 'int main'; then
        # Add guard around main function if not present
        sed -i 's/^int main(int argc, char \\*\\*argv)/\\n#ifndef SHARED_MODULE\\nint main(int argc, char **argv)/' src/main.c 2>/dev/null || true
        # Find the end of main and add endif (this is approximate - newer versions may already have this)
    fi
fi

# Fix misc_time.c - getrusage doesn't work properly in WASM
if [ -f src/misc/misc_time.c ]; then
    # Add EMSCRIPTEN guard around getrusage usage
    sed -i 's/#ifdef HAVE_GETRUSAGE/#if defined(HAVE_GETRUSAGE) \\&\\& !defined(__EMSCRIPTEN__)/' src/misc/misc_time.c 2>/dev/null || true
fi

# Fix sharedspice.c - skip user config file loading for WASM
if [ -f src/sharedspice.c ]; then
    # Add EMSCRIPTEN guard to skip config file access
    sed -i 's/read_initialisation_file()/#ifndef __EMSCRIPTEN__\\n    read_initialisation_file()\\n#endif/' src/sharedspice.c 2>/dev/null || true
fi

# Fix frontend/outitf.c - function signature mismatch
if [ -f src/frontend/outitf.c ]; then
    # The sh_vecinit signature fix - search for the function and ensure proper signature
    :
fi

echo 'WASM patches applied successfully'
")
        
        ExternalProject_Add(ngspice_external
            URL ${NGSPICE_URL}
            PREFIX ${NGSPICE_PREFIX}
            PATCH_COMMAND bash ${NGSPICE_PATCH_SCRIPT} <SOURCE_DIR>
            CONFIGURE_COMMAND emconfigure <SOURCE_DIR>/configure ${NGSPICE_CONFIGURE_FLAGS} --prefix=${NGSPICE_INSTALL_DIR}
            BUILD_COMMAND emmake make -j${CMAKE_BUILD_PARALLEL_LEVEL}
            INSTALL_COMMAND make install
            BUILD_IN_SOURCE 0
            # Note: Despite .so extension, this is actually a WASM static library
            # Version 45.2 produces libngspice.so.0.0.14
            BUILD_BYPRODUCTS ${NGSPICE_INSTALL_DIR}/lib/libngspice.so.0.0.14 ${NGSPICE_INSTALL_DIR}/lib/libngspice.so
        )
    else()
        # Native build configuration - basic ngspice only
        set(NGSPICE_CONFIGURE_FLAGS
            "--with-ngshared"
            "--disable-debug"
            "--without-x"
            "--disable-xspice"
            "--disable-cider"
            "--disable-openmp"
            "CFLAGS=-O2 -fPIC"
        )
        ExternalProject_Add(ngspice_external
            URL ${NGSPICE_URL}
            PREFIX ${NGSPICE_PREFIX}
            CONFIGURE_COMMAND <SOURCE_DIR>/configure ${NGSPICE_CONFIGURE_FLAGS} --prefix=${NGSPICE_INSTALL_DIR}
            BUILD_COMMAND make -j${CMAKE_BUILD_PARALLEL_LEVEL}
            INSTALL_COMMAND make install
            BUILD_IN_SOURCE 0
            BUILD_BYPRODUCTS ${NGSPICE_INSTALL_DIR}/lib/libngspice.so
        )
    endif()
    
    endif()
    
    # Create imported target for ngspice
    add_library(ngspice INTERFACE)
    if(NOT NGSPICE_ALREADY_INSTALLED)
        add_dependencies(ngspice ngspice_external)
    endif()
    
    # Set include and library paths
    target_include_directories(ngspice INTERFACE 
        ${NGSPICE_INSTALL_DIR}/include
        ${NGSPICE_INSTALL_DIR}/include/ngspice
    )
    
    if(WIN32)
        target_link_libraries(ngspice INTERFACE ${NGSPICE_INSTALL_DIR}/lib/ngspice.lib)
    elseif(EMSCRIPTEN)
        # Despite .so extension, this is a WASM static library when built with emscripten
        # Use the symlink for version independence
        target_link_libraries(ngspice INTERFACE ${NGSPICE_INSTALL_DIR}/lib/libngspice.so)
    else()
        target_link_libraries(ngspice INTERFACE ${NGSPICE_INSTALL_DIR}/lib/libngspice.so)
    endif()
    
    # Add compile definition
    add_compile_definitions(ENABLE_NGSPICE)
endif()


# Testing framework - Catch2
if(INCLUDE_MKF_TESTS)
    FetchContent_Declare(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.8.1)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    include(CTest)
    include(Catch)
endif()

# MAS
# NOTE: the previous form `find_program(NAMES quicktype REQUIRED)` was malformed — the
# result variable was literally named "NAMES", so REQUIRED never actually verified quicktype
# existed. The libgtest-dev probe was doubly wrong (gtest isn't used — Catch2 is — and it is
# a package, not a program) and has been dropped.
find_program(QUICKTYPE_EXECUTABLE NAMES quicktype REQUIRED)

if(NOT DEFINED MAS_DIRECTORY)
    set(MAS_DIRECTORY "${CMAKE_BINARY_DIR}/MAS/")
endif(NOT DEFINED MAS_DIRECTORY)

if(NOT DEFINED MAS_DIR)
    set(MAS_DIR "${PROJECT_SOURCE_DIR}/MAS")
endif(NOT DEFINED MAS_DIR)

if(NOT DEFINED PEAS_DIR)
    set(PEAS_DIR "${PROJECT_SOURCE_DIR}/PEAS")
endif()

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/MAS")

add_custom_command(
  OUTPUT "${MAS_DIRECTORY}/MAS.hpp"
  COMMAND quicktype -l c++ -s schema ${MAS_DIR}/schemas/MAS.json 
    -S ${MAS_DIR}/schemas/magnetic.json
    -S ${MAS_DIR}/schemas/magnetic/core.json
    -S ${MAS_DIR}/schemas/magnetic/coil.json
    -S ${MAS_DIR}/schemas/utils.json
    -S ${MAS_DIR}/schemas/magnetic/core/gap.json
    -S ${MAS_DIR}/schemas/magnetic/core/shape.json
    -S ${MAS_DIR}/schemas/magnetic/core/material.json
    -S ${MAS_DIR}/schemas/magnetic/core/coating.json
    -S ${MAS_DIR}/schemas/magnetic/insulation/material.json
    -S ${MAS_DIR}/schemas/magnetic/insulation/wireCoating.json
    -S ${MAS_DIR}/schemas/magnetic/bobbin.json
    -S ${MAS_DIR}/schemas/magnetic/core/piece.json
    -S ${MAS_DIR}/schemas/magnetic/core/spacer.json
    -S ${MAS_DIR}/schemas/magnetic/wire/basicWire.json
    -S ${MAS_DIR}/schemas/magnetic/wire/round.json
    -S ${MAS_DIR}/schemas/magnetic/wire/rectangular.json
    -S ${MAS_DIR}/schemas/magnetic/wire/foil.json
    -S ${MAS_DIR}/schemas/magnetic/wire/planar.json
    -S ${MAS_DIR}/schemas/magnetic/wire/litz.json
    -S ${MAS_DIR}/schemas/magnetic/wire/material.json
    -S ${MAS_DIR}/schemas/magnetic/wire.json
    -S ${MAS_DIR}/schemas/utils.json
    -S ${MAS_DIR}/schemas/magnetic/insulation/wireCoating.json
    -S ${MAS_DIR}/schemas/magnetic/insulation/material.json
    -S ${MAS_DIR}/schemas/inputs.json
    -S ${MAS_DIR}/schemas/outputs.json
    -S ${MAS_DIR}/schemas/outputs/coreLossesOutput.json
    -S ${MAS_DIR}/schemas/inputs/designRequirements.json
    -S ${MAS_DIR}/schemas/inputs/operatingConditions.json
    -S ${MAS_DIR}/schemas/inputs/operatingPoint.json
    -S ${MAS_DIR}/schemas/inputs/topologies/flyback.json
    -S ${MAS_DIR}/schemas/inputs/topologies/currentTransformer.json
    -S ${MAS_DIR}/schemas/inputs/topologies/boost.json
    -S ${MAS_DIR}/schemas/inputs/topologies/buck.json
    -S ${MAS_DIR}/schemas/inputs/topologies/forward.json
    -S ${MAS_DIR}/schemas/inputs/topologies/isolatedBuck.json
    -S ${MAS_DIR}/schemas/inputs/topologies/isolatedBuckBoost.json
    -S ${MAS_DIR}/schemas/inputs/topologies/pushPull.json
    -S ${MAS_DIR}/schemas/inputs/topologies/cllcResonant.json
    -S ${MAS_DIR}/schemas/inputs/topologies/dualActiveBridge.json
    -S ${MAS_DIR}/schemas/inputs/topologies/llcResonant.json
    -S ${MAS_DIR}/schemas/inputs/topologies/phaseShiftedFullBridge.json
    -S ${MAS_DIR}/schemas/inputs/topologies/phaseShiftedHalfBridge.json
    -S ${MAS_DIR}/schemas/inputs/topologies/asymmetricHalfBridge.json
    -S ${MAS_DIR}/schemas/inputs/topologies/cuk.json
    -S ${MAS_DIR}/schemas/inputs/topologies/clllcResonant.json
    -S ${MAS_DIR}/schemas/inputs/topologies/powerFactorCorrection.json
    -S ${MAS_DIR}/schemas/inputs/topologies/commonModeChoke.json
    -S ${MAS_DIR}/schemas/inputs/topologies/differentialModeChoke.json
    -S ${MAS_DIR}/schemas/inputs/topologies/sepic.json
    -S ${MAS_DIR}/schemas/inputs/topologies/zeta.json
    -S ${MAS_DIR}/schemas/inputs/topologies/weinberg.json
    -S ${MAS_DIR}/schemas/inputs/topologies/fourSwitchBuckBoost.json
    -S ${MAS_DIR}/schemas/inputs/topologies/seriesResonant.json
    -S ${MAS_DIR}/schemas/inputs/topologies/vienna.json
    -S ${PEAS_DIR}/schemas/utils.json
    -S ${PEAS_DIR}/schemas/inputs/operatingPointExcitation.json
    -S ${PEAS_DIR}/schemas/inputs/operatingConditions.json
    -S ${PEAS_DIR}/schemas/outputs/insulationCoordination.json
    -S ${PEAS_DIR}/schemas/outputs/lossElementPerHarmonic.json
    -o ${MAS_DIRECTORY}/MAS.hpp --namespace MAS --source-style single-source --type-style pascal-case --member-style underscore-case --enumerator-style upper-underscore-case --no-boost
  USES_TERMINAL)

add_custom_target(MAS ALL
    DEPENDS "${MAS_DIRECTORY}/MAS.hpp")

include_directories("${MAS_DIRECTORY}")

# CAS
if(NOT DEFINED CAS_DIRECTORY)
    set(CAS_DIRECTORY "${CMAKE_BINARY_DIR}/CAS/")
endif()

if(NOT DEFINED CAS_DIR)
    set(CAS_DIR "${PROJECT_SOURCE_DIR}/CAS")
endif()

if(NOT DEFINED PEAS_DIR)
    set(PEAS_DIR "${PROJECT_SOURCE_DIR}/PEAS")
endif()

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/CAS")

add_custom_command(
  OUTPUT "${CAS_DIRECTORY}/CAS.hpp"
  COMMAND quicktype -l c++ -s schema ${CAS_DIR}/schemas/inputs.json
    -S ${CAS_DIR}/schemas/inputs/designRequirements.json
    -S ${PEAS_DIR}/schemas/utils.json
    -S ${PEAS_DIR}/schemas/inputs/twoTerminalOperatingPoint.json
    -S ${PEAS_DIR}/schemas/inputs/operatingConditions.json
    -S ${PEAS_DIR}/schemas/inputs/operatingPointExcitation.json
    -o ${CAS_DIRECTORY}/CAS.hpp --namespace CAS --source-style single-source
    --type-style pascal-case --member-style underscore-case
    --enumerator-style upper-underscore-case --no-boost
    --top-level Inputs
  DEPENDS
    "${CAS_DIR}/schemas/inputs.json"
    "${CAS_DIR}/schemas/inputs/designRequirements.json"
    "${PEAS_DIR}/schemas/utils.json"
    "${PEAS_DIR}/schemas/inputs/twoTerminalOperatingPoint.json"
    "${PEAS_DIR}/schemas/inputs/operatingConditions.json"
    "${PEAS_DIR}/schemas/inputs/operatingPointExcitation.json"
  USES_TERMINAL)

add_custom_target(CAS ALL DEPENDS "${CAS_DIRECTORY}/CAS.hpp")
add_dependencies(CAS MAS)

include_directories("${CAS_DIRECTORY}")


# ============================================================================
# CCI strand coordinates — embed up to N=1000 into the binary at build time
# ============================================================================
set(CCI_COORDS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cci_coords"
    CACHE PATH "Path to cci_coords repository (submodule at MKF/cci_coords)")
set(CCI_GENERATED_CPP "${CMAKE_BINARY_DIR}/generated/CciCoordinatesData.cpp")

find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_custom_command(
    OUTPUT "${CCI_GENERATED_CPP}"
    COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/scripts/generate_cci_data.py"
            "${CCI_COORDS_DIR}"
            "${CCI_GENERATED_CPP}"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/scripts/generate_cci_data.py"
    COMMENT "Generating embedded CCI coordinates"
    VERBATIM)
add_custom_target(cci_data_gen DEPENDS "${CCI_GENERATED_CPP}")

file(GLOB SOURCES CONFIGURE_DEPENDS
    "src/*.cpp"
    "src/advisers/*.cpp"
    "src/constructive_models/*.cpp"
    "src/converter_models/*.cpp"
    "src/physical_models/*.cpp"
    "src/processors/*.cpp"
    "src/support/*.cpp")
list(APPEND SOURCES "${CCI_GENERATED_CPP}")

# The only surviving src/converter_models/*.cpp is KirchhoffBridge.cpp — it #includes KirchhoffApi.hpp and
# links libKirchhoffApi.so, so it can only be compiled when Kirchhoff is linked. Drop it otherwise (e.g. the
# Emscripten/webMKF build, MKF_USE_KIRCHHOFF=OFF).
if(NOT MKF_USE_KIRCHHOFF)
    list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/converter_models/KirchhoffBridge.cpp")
endif()

# The old converter_models were DELETED — replaced by Kirchhoff (linked via KirchhoffBridge, the only
# remaining src/converter_models/*.cpp). MKF's core never instantiated a converter class (the
# MagneticAdviser converter path is generic templates the CALLER feeds a converter to); the converter
# physics now lives in Kirchhoff and its reference designs are SPICE-gated there.

SET(LIBS data::insulation_standards)

# REMOVED: LINK_DIRECTORIES for LAPACKBLAS (no longer needed, levmar replaced by Eigen)

# libraries the demo depends on
IF(HAVE_PLASMA)
    LINK_DIRECTORIES(${PLASMA_DIR}/lib)
    SET(LIBS ${LIBS} ${PLASMA_LIB_NAMES})
ENDIF(HAVE_PLASMA)


# On Windows, build as STATIC to avoid DLL boundary issues with inline globals
# On other platforms, build as SHARED for better compatibility
if(WIN32)
    add_library(MKF STATIC ${SOURCES})
else()
    add_library(MKF SHARED ${SOURCES})
endif()
set_property(TARGET MKF PROPERTY POSITION_INDEPENDENT_CODE ON)
add_dependencies(MKF MAS CAS cci_data_gen)
target_include_directories(MKF PUBLIC "${CMAKE_BINARY_DIR}/generated")

cmrc_add_resource_library(insulation_standards ALIAS data::insulation_standards NAMESPACE insulationData src/data/insulation_standards/IEC_60664-1.json src/data/insulation_standards/IEC_60664-4.json src/data/insulation_standards/IEC_60664-5.json src/data/insulation_standards/IEC_62368-1.json src/data/insulation_standards/IEC_61558-1.json src/data/insulation_standards/IEC_61558-2-16.json src/data/insulation_standards/IEC_60335-1.json)
cmrc_add_resource_library(core_losses_data ALIAS data::core_losses_data NAMESPACE coreLossesData src/data/core_losses/ciGSE_coefficients.json)
target_link_libraries(MKF PUBLIC ${LIBS})

# Kirchhoff bridge: MKF links libKirchhoffApi.so and includes ONLY KirchhoffApi.hpp (the string facade),
# so KH's `namespace MAS` types never enter an MKF TU. src/converter_models/KirchhoffBridge.cpp is picked
# up by the SOURCES glob; here we add the include path + the imported .so + the build-order dependency.
if(MKF_USE_KIRCHHOFF)
    add_dependencies(MKF kirchhoff_ext)
    # The KH source tree lives at KIRCHHOFF_SOURCE_DIR when a local checkout is passed, or at the
    # ExternalProject fetch location when we clone from GitHub (the default). Pick the right one so the
    # KirchhoffApi.hpp include resolves in BOTH modes — an empty KIRCHHOFF_SOURCE_DIR must NOT collapse
    # to "/src" (which silently drops the include and breaks the KirchhoffBridge.cpp compile).
    if(KIRCHHOFF_SOURCE_DIR)
        set(_kirchhoff_include_dir "${KIRCHHOFF_SOURCE_DIR}/src")
    else()
        set(_kirchhoff_include_dir "${CMAKE_BINARY_DIR}/kirchhoff_ext-prefix/src/kirchhoff_ext/src")
    endif()
    # PUBLIC so consumers/tests that include KirchhoffBridge.h (-> KirchhoffApi.hpp) find it too.
    target_include_directories(MKF PUBLIC "${_kirchhoff_include_dir}")
    target_link_libraries(MKF PUBLIC KirchhoffApi)
endif()

set(DATA_LIST "")
if(EMBED_MAS_DATA)
    if(EMBED_MAS_CORES)
        list (APPEND DATA_LIST MAS/data/cores_stock.ndjson)
        list (APPEND DATA_LIST MAS/data/cores.ndjson)
    endif()
    if(EMBED_MAS_CORE_SHAPES)
        list (APPEND DATA_LIST MAS/data/core_shapes.ndjson)
    endif()
    if(EMBED_MAS_CORE_MATERIALS)
        list (APPEND DATA_LIST MAS/data/core_materials.ndjson)
    endif()
    if(EMBED_MAS_WIRES)
        list (APPEND DATA_LIST MAS/data/wires.ndjson)
    endif()
    if(EMBED_MAS_WIRE_MATERIALS)
        list (APPEND DATA_LIST MAS/data/wire_materials.ndjson)
    endif()
    if(EMBED_MAS_BOBBINS)
        list (APPEND DATA_LIST MAS/data/bobbins.ndjson)
    endif()
    if(EMBED_MAS_INSULATION_MATERIALS)
        list (APPEND DATA_LIST MAS/data/insulation_materials.ndjson)
    endif()
endif()

message(STATUS "MKF embedded MAS data files: ${DATA_LIST}")
cmrc_add_resource_library(data ALIAS data::data NAMESPACE data ${DATA_LIST})
target_link_libraries(MKF PUBLIC data::data)
target_link_libraries(MKF PUBLIC data::core_losses_data)

target_link_libraries(MKF PUBLIC rapidfuzz::rapidfuzz)
target_link_libraries(MKF PRIVATE nlohmann_json::nlohmann_json)

# Link ngspice if enabled
if(ENABLE_NGSPICE)
    target_link_libraries(MKF PUBLIC ngspice)
    if(NOT NGSPICE_ALREADY_INSTALLED)
        add_dependencies(MKF ngspice_external)
    endif()
endif()

target_include_directories(MKF PUBLIC
     "${PROJECT_SOURCE_DIR}/src/"
     "${PROJECT_SOURCE_DIR}/src/advisers/"
     "${PROJECT_SOURCE_DIR}/src/constructive_models/"
     "${PROJECT_SOURCE_DIR}/src/converter_models/"
     "${PROJECT_SOURCE_DIR}/src/physical_models/"
     "${PROJECT_SOURCE_DIR}/src/processors/"
     "${PROJECT_SOURCE_DIR}/src/support/"
     "${CMAKE_BINARY_DIR}/_deps/eigen-src"
     )

# add_executable(magnetic_adviser "src/advisers/MagneticAdviser.cpp")
# target_link_libraries(magnetic_adviser MKF)

if(INCLUDE_MKF_TESTS)
    file(GLOB ALL_TEST_SOURCES CONFIGURE_DEPENDS "tests/*.cpp")
    
    # TestingUtils.cpp holds shared fixtures; explicitly added below.
    list(FILTER ALL_TEST_SOURCES EXCLUDE REGEX "TestingUtils\\.cpp$")
    
    add_executable(MKF_tests ${ALL_TEST_SOURCES} "tests/TestingUtils.cpp")
    target_link_libraries(MKF_tests PRIVATE Catch2::Catch2WithMain)
    target_link_libraries(MKF_tests PUBLIC MKF)
    target_include_directories(MKF_tests PUBLIC "${PROJECT_SOURCE_DIR}/src/")

    # Precompiled header for the test target.
    #
    # The TestTopology*.cpp files (3k LOC each, 10+ of them) each parse the
    # generated MAS.hpp (~700KB of template-heavy quicktype output), Catch2's
    # header-only `catch_all.hpp`, nlohmann/json, MKF's own constructive +
    # converter model headers — every TU re-parses and re-instantiates the
    # same template machinery. Pre-compiling once amortises the cost across
    # the whole test build: empirically ~5-10x faster compile on the heavy
    # TestTopology* files with ZERO runtime impact (PCH is a parser-state
    # cache, not a code change — the generated object files are identical).
    target_precompile_headers(MKF_tests PRIVATE
        <catch2/catch_all.hpp>
        <nlohmann/json.hpp>
        <MAS.hpp>
        <vector>
        <string>
        <optional>
        <memory>
        <map>
        <variant>
    )

    # Discover and register tests with CTest
    catch_discover_tests(MKF_tests)


endif(INCLUDE_MKF_TESTS)