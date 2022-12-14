cmake_minimum_required(VERSION 3.18.0)

project(MKF_Project)

SET(CMAKE_CXX_STANDARD 20)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

SET(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-deprecated-declarations -Wno-unused-parameter -Wno-switch")

option(BUILD_TESTS      "Build tests"    OFF)
option(BUILD_EXAMPLES   "Build examples" OFF)
option(UTPP_INCLUDE_TESTS_IN_BUILD   "Build tests" OFF)

include(FetchContent)

FetchContent_Declare(UnitTest++
        GIT_REPOSITORY https://github.com/unittest-cpp/unittest-cpp.git
        GIT_TAG  tags/v2.0.0)

FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/nlohmann/json/
        GIT_TAG  tags/v3.10.0)

FetchContent_Declare(json-schema-validator
        GIT_REPOSITORY https://github.com/pboettch/json-schema-validator.git
        GIT_TAG  tags/2.1.0)

FetchContent_Declare(magic-enum
        GIT_REPOSITORY https://github.com/Neargye/magic_enum
        GIT_TAG  tags/v0.8.1)

FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git)

FetchContent_Declare(pybind11_json
        GIT_REPOSITORY https://github.com/pybind/pybind11_json.git)

FetchContent_MakeAvailable(UnitTest++ json json-schema-validator magic-enum pybind11 pybind11_json)

set_target_properties(pybind11_json PROPERTIES EXCLUDE_FROM_ALL True)


SET(CMAKE_CXX_FLAGS "-Wall -Wextra")

find_program(NAMES quicktype REQUIRED)
find_program(NAMES libgtest-dev REQUIRED)

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/MAS")

set(MAS_DIR "${PROJECT_SOURCE_DIR}/../MAS")


add_custom_command(
  OUTPUT "${CMAKE_BINARY_DIR}/MAS/MAS.hpp"
  COMMAND quicktype -l c++ -s schema ${MAS_DIR}/schemas/MAS.json -S ${MAS_DIR}/schemas/magnetic.json -S ${MAS_DIR}/schemas/magnetic/core.json -S ${MAS_DIR}/schemas/magnetic/winding.json -S ${MAS_DIR}/schemas/utils.json -S ${MAS_DIR}/schemas/magnetic/core/gap.json -S ${MAS_DIR}/schemas/magnetic/core/shape.json -S ${MAS_DIR}/schemas/magnetic/core/material.json -S ${MAS_DIR}/schemas/magnetic/insulation/material.json -S ${MAS_DIR}/schemas/magnetic/insulation/wireCoating.json -S ${MAS_DIR}/schemas/magnetic/core/bobbin.json -S ${MAS_DIR}/schemas/magnetic/core/piece.json -S ${MAS_DIR}/schemas/magnetic/core/spacer.json -S ${MAS_DIR}/schemas/magnetic/wire/solid.json -S ${MAS_DIR}/schemas/magnetic/wire/stranded.json -S ${MAS_DIR}/schemas/magnetic/wire/material.json -S ${MAS_DIR}/schemas/utils.json -S ${MAS_DIR}/schemas/magnetic/insulation/wireCoating.json -S ${MAS_DIR}/schemas/magnetic/insulation/material.json -S ${MAS_DIR}/schemas/inputs.json -S ${MAS_DIR}/schemas/inputs/designRequirements.json -S ${MAS_DIR}/schemas/inputs/operationConditions.json -S ${MAS_DIR}/schemas/inputs/operationPointExcitation.json -o ${CMAKE_BINARY_DIR}/MAS/MAS.hpp --namespace OpenMagnetics --source-style single-source --type-style pascal-case --member-style underscore-case --enumerator-style upper-underscore-case --no-boost
  USES_TERMINAL)

add_custom_target(MAS ALL
    DEPENDS "${CMAKE_BINARY_DIR}/MAS/MAS.hpp")

file(GLOB SOURCES
    "src/*.cpp")

add_library(MKF SHARED ${SOURCES})

include_directories("${CMAKE_BINARY_DIR}/_deps/pybind11_json-src/include/pybind11_json/")
include_directories("${CMAKE_BINARY_DIR}/_deps/unittest++-src/UnitTest++/")
include_directories("${CMAKE_BINARY_DIR}/_deps/json-src/include/nlohmann/")
include_directories("${CMAKE_BINARY_DIR}/_deps/magic-enum-src/include")
include_directories("${CMAKE_BINARY_DIR}/_deps/json-src/include/")
include_directories("${CMAKE_BINARY_DIR}/MAS/")

target_link_libraries(MKF nlohmann_json::nlohmann_json nlohmann_json_schema_validator UnitTest++ pybind11_json)
target_include_directories(MKF PUBLIC src/)

file(GLOB TESTS_SOURCES
    "tests/*.cpp")

add_executable(MKF_tests ${TESTS_SOURCES})
target_link_libraries(MKF_tests MKF)

include_directories("${CMAKE_BINARY_DIR}/_deps/pybind11-src")
pybind11_add_module(PyMKF src/PyMKF/PyMKFWrapper.cpp )
target_link_libraries(PyMKF PRIVATE MKF pybind11_json)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.

target_compile_definitions(PyMKF
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})