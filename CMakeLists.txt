cmake_minimum_required(VERSION 3.18.0)

project(MKF_Project)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gdwarf-3")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-deprecated-declarations -Wno-unused-parameter -Wno-switch")

option(BUILD_TESTS      "Build tests"    OFF)
option(BUILD_EXAMPLES   "Build examples" OFF)
option(UTPP_INCLUDE_TESTS_IN_BUILD   "Build tests" OFF)
option(INCLUDE_MKF_TESTS "Compiles also the tests" ON)
option(INCLUDE_PYMKF "Compiles also the Python interface" OFF)


if(NOT DEFINED MAS_DIRECTORY)
        set(MAS_DIRECTORY "${CMAKE_BINARY_DIR}/MAS/")
endif(NOT DEFINED MAS_DIRECTORY)

if(NOT DEFINED MAS_DIR)
        set(MAS_DIR "${PROJECT_SOURCE_DIR}/../MAS")
endif(NOT DEFINED MAS_DIR)

include(FetchContent)

if(INCLUDE_MKF_TESTS)
FetchContent_Declare(UnitTest++
        GIT_REPOSITORY https://github.com/unittest-cpp/unittest-cpp.git
        GIT_TAG  tags/v2.0.0)
endif(INCLUDE_MKF_TESTS)

FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/nlohmann/json/
        GIT_TAG  tags/v3.11.2)

FetchContent_Declare(magic-enum
        GIT_REPOSITORY https://github.com/Neargye/magic_enum
        GIT_TAG  tags/v0.9.3)

if(INCLUDE_PYMKF)
        FetchContent_Declare(pybind11
                GIT_REPOSITORY https://github.com/pybind/pybind11.git)

        FetchContent_Declare(pybind11_json
                GIT_REPOSITORY https://github.com/pybind/pybind11_json.git)
endif(INCLUDE_PYMKF)

FetchContent_Declare(spline
        GIT_REPOSITORY https://github.com/ttk592/spline.git)

FetchContent_Declare(svg
        GIT_REPOSITORY https://github.com/AlfVII/svg)

FetchContent_MakeAvailable(json magic-enum spline svg)

if(INCLUDE_MKF_TESTS)
        FetchContent_MakeAvailable(UnitTest++)
endif(INCLUDE_MKF_TESTS)

if(INCLUDE_PYMKF)
        FetchContent_MakeAvailable( pybind11 pybind11_json)
        set_target_properties(pybind11_json PROPERTIES EXCLUDE_FROM_ALL True)
endif(INCLUDE_PYMKF)

file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/_deps/svg-src/CMakeLists.txt")
file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/_deps/svg-src/tests")

set(CMAKE_CXX_FLAGS "-Wall -Wextra")

find_program(NAMES quicktype REQUIRED)
find_program(NAMES libgtest-dev REQUIRED)

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/MAS")



add_custom_command(
  OUTPUT "${MAS_DIRECTORY}/MAS.hpp"
  COMMAND quicktype -l c++ -s schema ${MAS_DIR}/schemas/MAS.json -S ${MAS_DIR}/schemas/magnetic.json -S ${MAS_DIR}/schemas/magnetic/core.json -S ${MAS_DIR}/schemas/magnetic/coil.json -S ${MAS_DIR}/schemas/utils.json -S ${MAS_DIR}/schemas/magnetic/core/gap.json -S ${MAS_DIR}/schemas/magnetic/core/shape.json -S ${MAS_DIR}/schemas/magnetic/core/material.json -S ${MAS_DIR}/schemas/magnetic/insulation/material.json -S ${MAS_DIR}/schemas/magnetic/insulation/wireCoating.json -S ${MAS_DIR}/schemas/magnetic/bobbin.json -S ${MAS_DIR}/schemas/magnetic/core/piece.json -S ${MAS_DIR}/schemas/magnetic/core/spacer.json -S ${MAS_DIR}/schemas/magnetic/wire/solid.json -S ${MAS_DIR}/schemas/magnetic/wire/stranded.json -S ${MAS_DIR}/schemas/magnetic/wire/material.json -S ${MAS_DIR}/schemas/magnetic/wire.json -S ${MAS_DIR}/schemas/utils.json -S ${MAS_DIR}/schemas/magnetic/insulation/wireCoating.json -S ${MAS_DIR}/schemas/magnetic/insulation/material.json -S ${MAS_DIR}/schemas/inputs.json -S ${MAS_DIR}/schemas/outputs.json -S ${MAS_DIR}/schemas/inputs/designRequirements.json -S ${MAS_DIR}/schemas/inputs/operatingConditions.json -S ${MAS_DIR}/schemas/inputs/operatingPointExcitation.json -o ${MAS_DIRECTORY}/MAS.hpp --namespace OpenMagnetics --source-style single-source --type-style pascal-case --member-style underscore-case --enumerator-style upper-underscore-case --no-boost
  USES_TERMINAL)

add_custom_command(
  OUTPUT "${MAS_DIRECTORY}/MAS.hpp"
  COMMAND quicktype -l c++ -s schema ${MAS_DIR}/schemas/MAS.json -S ${MAS_DIR}/schemas/magnetic.json -S ${MAS_DIR}/schemas/magnetic/core.json -S ${MAS_DIR}/schemas/magnetic/coil.json -S ${MAS_DIR}/schemas/utils.json -S ${MAS_DIR}/schemas/magnetic/core/gap.json -S ${MAS_DIR}/schemas/magnetic/core/shape.json -S ${MAS_DIR}/schemas/magnetic/core/material.json -S ${MAS_DIR}/schemas/magnetic/insulation/material.json -S ${MAS_DIR}/schemas/magnetic/insulation/wireCoating.json -S ${MAS_DIR}/schemas/magnetic/bobbin.json -S ${MAS_DIR}/schemas/magnetic/core/piece.json -S ${MAS_DIR}/schemas/magnetic/core/spacer.json -S ${MAS_DIR}/schemas/magnetic/wire/solid.json -S ${MAS_DIR}/schemas/magnetic/wire/stranded.json -S ${MAS_DIR}/schemas/magnetic/wire/material.json -S ${MAS_DIR}/schemas/magnetic/wire.json -S ${MAS_DIR}/schemas/utils.json -S ${MAS_DIR}/schemas/magnetic/insulation/wireCoating.json -S ${MAS_DIR}/schemas/magnetic/insulation/material.json -S ${MAS_DIR}/schemas/inputs.json -S ${MAS_DIR}/schemas/outputs.json -S ${MAS_DIR}/schemas/inputs/designRequirements.json -S ${MAS_DIR}/schemas/inputs/operatingConditions.json -S ${MAS_DIR}/schemas/inputs/operatingPointExcitation.json -o ${MAS_DIRECTORY}/MAS.hpp --namespace OpenMagnetics --source-style single-source --type-style pascal-case --member-style underscore-case --enumerator-style upper-underscore-case --no-boost
  USES_TERMINAL)

add_custom_target(MAS ALL
    DEPENDS "${MAS_DIRECTORY}/MAS.hpp")


file(GLOB SOURCES
    "src/*.cpp")

add_library(MKF SHARED ${SOURCES})


add_custom_target(MASGeneration
                  /bin/echo "RUNNING MASGeneration"
                  DEPENDS "${MAS_DIRECTORY}/MAS.hpp")
add_dependencies(MKF MASGeneration)

if(INCLUDE_PYMKF)
        add_custom_target(PyMASGeneration
                          /bin/echo "RUNNING PyMASGeneration"
                          DEPENDS "${MAS_DIRECTORY}/MAS.hpp")
        add_dependencies(MKF PyMASGeneration)
endif(INCLUDE_PYMKF)


if(INCLUDE_PYMKF)
        include_directories("${CMAKE_BINARY_DIR}/_deps/pybind11-src/include/")
        include_directories("${CMAKE_BINARY_DIR}/_deps/pybind11_json-src/include/pybind11_json/")
endif(INCLUDE_PYMKF)
if(INCLUDE_MKF_TESTS)
        include_directories("${CMAKE_BINARY_DIR}/_deps/unittest++-src/UnitTest++/")
endif(INCLUDE_MKF_TESTS)
include_directories("${CMAKE_BINARY_DIR}/_deps/json-src/include/nlohmann/")
include_directories("${CMAKE_BINARY_DIR}/_deps/magic-enum-src/include")
include_directories("${CMAKE_BINARY_DIR}/_deps/svg-src/src")
include_directories("${CMAKE_BINARY_DIR}/_deps/spline-src/src")
include_directories("${CMAKE_BINARY_DIR}/_deps/json-src/include/")
include_directories("${MAS_DIRECTORY}/")

target_link_libraries(MKF nlohmann_json::nlohmann_json)
if(INCLUDE_MKF_TESTS)
        target_link_libraries(MKF UnitTest++)
endif(INCLUDE_MKF_TESTS)

if(INCLUDE_PYMKF)
        target_link_libraries(MKF pybind11_json pybind11::embed)
endif(INCLUDE_PYMKF)

target_include_directories(MKF PUBLIC src/)

if(INCLUDE_MKF_TESTS)
        file(GLOB TESTS_SOURCES
            "tests/*.cpp")
        add_executable(MKF_tests ${TESTS_SOURCES})
        target_link_libraries(MKF_tests MKF)
endif(INCLUDE_MKF_TESTS)

if(INCLUDE_PYMKF)
        include_directories("${CMAKE_BINARY_DIR}/_deps/pybind11-src")
        pybind11_add_module(PyMKF src/PyMKF/PyMKFWrapper.cpp )
        target_link_libraries(PyMKF PRIVATE MKF pybind11_json)

        # EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
        # define (VERSION_INFO) here.

        target_compile_definitions(PyMKF
                                   PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
endif(INCLUDE_PYMKF)