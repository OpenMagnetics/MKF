cmake_minimum_required(VERSION 3.18.0)

project(MKF_tests)

SET(CMAKE_CXX_STANDARD 20)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

SET(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-deprecated-declarations -Wno-unused-parameter -Wno-switch")

option(BUILD_TESTS      "Build tests"    OFF)
option(BUILD_EXAMPLES   "Build examples" OFF)

include(FetchContent)

FetchContent_Declare(UnitTest++
        GIT_REPOSITORY https://github.com/unittest-cpp/unittest-cpp.git
        GIT_TAG  tags/v2.0.0)

FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/nlohmann/json/
        GIT_TAG  tags/v3.10.0)

FetchContent_Declare(magic-enum
        GIT_REPOSITORY https://github.com/Neargye/magic_enum
        GIT_TAG  tags/v0.8.1)

FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git)

FetchContent_Declare(pybind11_json
        GIT_REPOSITORY https://github.com/pybind/pybind11_json.git)

FetchContent_MakeAvailable(UnitTest++ json magic-enum pybind11 pybind11_json)

set_target_properties(pybind11_json PROPERTIES EXCLUDE_FROM_ALL True)


SET(CMAKE_CXX_FLAGS "-Wall -Wextra")

add_library(Constants STATIC src/Constants.h)

include_directories("${CMAKE_BINARY_DIR}/_deps/pybind11_json-src/include/pybind11_json/")
include_directories("${CMAKE_BINARY_DIR}/_deps/unittest++-src/UnitTest++/")
include_directories("${CMAKE_BINARY_DIR}/_deps/json-src/include/")
include_directories("${CMAKE_BINARY_DIR}/_deps/json-src/include/nlohmann/")
include_directories("${CMAKE_BINARY_DIR}/_deps/magic-enum-src/include")

target_link_libraries(Constants nlohmann_json::nlohmann_json UnitTest++ pybind11_json)
target_include_directories(Constants PUBLIC src/)

file(GLOB MKF_sources
    "tests/*.cpp")

add_executable(MKF_tests ${MKF_sources})
target_link_libraries(MKF_tests Constants)


add_subdirectory("${CMAKE_BINARY_DIR}/_deps/pybind11-src")
pybind11_add_module(PyMKF src/PyMKF/PyMKFWrapper.cpp )
target_link_libraries(PyMKF PRIVATE Constants pybind11_json)


# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.


target_compile_definitions(PyMKF
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})