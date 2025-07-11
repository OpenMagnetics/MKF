﻿cmake_minimum_required(VERSION 3.18.0)

if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("MKF_Project")


set(CMAKE_CXX_STANDARD 23) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /FS")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /FS")
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    set(CMAKE_BUILD_TYPE Release)
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-deprecated-declarations -Wno-unused-parameter -Wno-switch")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    # set(CMAKE_BUILD_TYPE RelWithDebInfo)
    # set(CMAKE_BUILD_TYPE MinSizeRel)
    set(CMAKE_BUILD_TYPE Release)
endif()


option(EMBED_MAS_DATA "Embed data in executable" ON)
option(EMBED_MAS_CORES "Embed cores in executable" ON)
option(EMBED_MAS_CORE_SHAPES "Embed cores in executable" ON)
option(EMBED_MAS_CORE_MATERIALS "Embed cores in executable" ON)
option(EMBED_MAS_WIRES "Embed cores in executable" ON)
option(EMBED_MAS_WIRE_MATERIALS "Embed cores in executable" ON)
option(EMBED_MAS_BOBBINS "Embed cores in executable" ON)
option(EMBED_MAS_INSULATION_MATERIALS "Embed cores in executable" ON)
option(INCLUDE_MKF_TESTS "Compiles also the tests" ON)
option(BUILD_TESTS      "Build tests"    OFF)
option(BUILD_EXAMPLES   "Build examples" OFF)
option(BUILD_DEMO   "Build examples" FALSE)
option(HAVE_LAPACK   "HAVE_LAPACK" 0)
option(UTPP_INCLUDE_TESTS_IN_BUILD   "Build tests" OFF)
option(INCLUDE_MKF_TESTS "Compiles also the tests" ON)
option(INCLUDE_PYMKF "Compiles also the Python interface" OFF)

message("Embeding MAS Data: ${EMBED_MAS_DATA}")
message("Embeding MAS Cores: ${EMBED_MAS_CORES}")
message("Embeding MAS CoreShapes: ${EMBED_MAS_CORE_SHAPES}")
message("Embeding MAS CoreMaterials: ${EMBED_MAS_CORE_MATERIALS}")
message("Embeding MAS Wires: ${EMBED_MAS_WIRES}")
message("Embeding MAS WireMaterials: ${EMBED_MAS_WIRE_MATERIALS}")
message("Embeding MAS Bobbins: ${EMBED_MAS_BOBBINS}")
message("Embeding MAS InsulationMaterials: ${EMBED_MAS_INSULATION_MATERIALS}")


file(DOWNLOAD "https://raw.githubusercontent.com/vector-of-bool/cmrc/master/CMakeRC.cmake"
                 "${CMAKE_BINARY_DIR}/CMakeRC.cmake")
include("${CMAKE_BINARY_DIR}/CMakeRC.cmake")


include(FetchContent)

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

FetchContent_Declare(matplotplusplus
    GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
    GIT_TAG tags/v1.2.1)
FetchContent_GetProperties(matplotplusplus)
if(NOT matplotplusplus_POPULATED)
    FetchContent_Populate(matplotplusplus)
    add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
    
FetchContent_Declare(spline
    GIT_REPOSITORY https://github.com/AlfVII/spline.git)
FetchContent_MakeAvailable(spline)
include_directories("${CMAKE_BINARY_DIR}/_deps/spline-src/src")

FetchContent_Declare(levmar
    GIT_REPOSITORY https://github.com/AlfVII/levmar.git
    GIT_TAG main)
FetchContent_MakeAvailable(levmar)
include_directories("${CMAKE_BINARY_DIR}/_deps/levmar-src")

FetchContent_Declare(svg
    GIT_REPOSITORY https://github.com/AlfVII/svg)
FetchContent_MakeAvailable(svg)
include_directories("${CMAKE_BINARY_DIR}/_deps/svg-src/src")


if(INCLUDE_MKF_TESTS)
FetchContent_Declare(UnitTest++
    GIT_REPOSITORY https://github.com/unittest-cpp/unittest-cpp.git
    GIT_TAG  tags/v2.0.0)
    FetchContent_MakeAvailable(UnitTest++)
    include_directories("${CMAKE_BINARY_DIR}/_deps/unittest++-src/UnitTest++/")
endif(INCLUDE_MKF_TESTS)

# MAS
find_program(NAMES quicktype REQUIRED)
find_program(NAMES libgtest-dev REQUIRED)

if(NOT DEFINED MAS_DIRECTORY)
    set(MAS_DIRECTORY "${CMAKE_BINARY_DIR}/MAS/")
endif(NOT DEFINED MAS_DIRECTORY)

if(NOT DEFINED MAS_DIR)
    set(MAS_DIR "${PROJECT_SOURCE_DIR}/MAS")
endif(NOT DEFINED MAS_DIR)

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/MAS")

add_custom_command(
  OUTPUT "${MAS_DIRECTORY}/MAS.hpp"
  COMMAND quicktype -l c++ -s schema ${MAS_DIR}/schemas/MAS.json -S ${MAS_DIR}/schemas/magnetic.json -S ${MAS_DIR}/schemas/magnetic/core.json -S ${MAS_DIR}/schemas/magnetic/coil.json -S ${MAS_DIR}/schemas/utils.json -S ${MAS_DIR}/schemas/magnetic/core/gap.json -S ${MAS_DIR}/schemas/magnetic/core/shape.json -S ${MAS_DIR}/schemas/magnetic/core/material.json -S ${MAS_DIR}/schemas/magnetic/insulation/material.json -S ${MAS_DIR}/schemas/magnetic/insulation/wireCoating.json -S ${MAS_DIR}/schemas/magnetic/bobbin.json -S ${MAS_DIR}/schemas/magnetic/core/piece.json -S ${MAS_DIR}/schemas/magnetic/core/spacer.json -S ${MAS_DIR}/schemas/magnetic/wire/basicWire.json -S ${MAS_DIR}/schemas/magnetic/wire/round.json -S ${MAS_DIR}/schemas/magnetic/wire/rectangular.json -S ${MAS_DIR}/schemas/magnetic/wire/foil.json -S ${MAS_DIR}/schemas/magnetic/wire/planar.json -S ${MAS_DIR}/schemas/magnetic/wire/litz.json -S ${MAS_DIR}/schemas/magnetic/wire/material.json -S ${MAS_DIR}/schemas/magnetic/wire.json -S ${MAS_DIR}/schemas/utils.json -S ${MAS_DIR}/schemas/magnetic/insulation/wireCoating.json -S ${MAS_DIR}/schemas/magnetic/insulation/material.json -S ${MAS_DIR}/schemas/inputs.json -S ${MAS_DIR}/schemas/outputs.json  -S ${MAS_DIR}/schemas/outputs/coreLossesOutput.json -S ${MAS_DIR}/schemas/inputs/designRequirements.json -S ${MAS_DIR}/schemas/inputs/operatingPoint.json -S ${MAS_DIR}/schemas/inputs/operatingConditions.json -S ${MAS_DIR}/schemas/inputs/operatingPointExcitation.json -S ${MAS_DIR}/schemas/inputs/topologies/flyback.json -o ${MAS_DIRECTORY}/MAS.hpp --namespace MAS --source-style single-source --type-style pascal-case --member-style underscore-case --enumerator-style upper-underscore-case --no-boost
  USES_TERMINAL)

add_custom_target(MAS ALL
    DEPENDS "${MAS_DIRECTORY}/MAS.hpp")

include_directories("${MAS_DIRECTORY}")

# message("Using for Levmar: ${CMAKE_CURRENT_BINARY_DIR}/_deps/levmar-2.6")
# if( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/_deps/levmar-2.6 )
#     if( NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/_deps/levmar-2.6.tgz )
#         file(DOWNLOAD https://users.ics.forth.gr/~lourakis/levmar/levmar-2.6.tgz ${CMAKE_CURRENT_BINARY_DIR}/_deps/levmar-2.6.tgz SHOW_PROGRESS)
#     endif()
    
#     execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xzvf ${CMAKE_CURRENT_BINARY_DIR}/_deps/levmar-2.6.tgz)
#     execute_process(COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_BINARY_DIR}/levmar-2.6 ${CMAKE_CURRENT_BINARY_DIR}/_deps/levmar-2.6)
# endif()
# add_subdirectory("${CMAKE_CURRENT_BINARY_DIR}/_deps/levmar-2.6" "${CMAKE_CURRENT_BINARY_DIR}/_deps/levmar-2.6")
# include_directories("${CMAKE_CURRENT_BINARY_DIR}/_deps/levmar-2.6")

file(GLOB SOURCES
    "src/*.cpp"
    "src/advisers/*.cpp"
    "src/constructive_models/*.cpp"
    "src/converter_models/*.cpp"
    "src/physical_models/*.cpp"
    "src/processors/*.cpp"
    "src/support/*.cpp")

SET(LIBS levmar data::insulation_standards)

LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) # location of the levmar library
LINK_DIRECTORIES(${LAPACKBLAS_DIR})

# libraries the demo depends on
IF(HAVE_PLASMA)
    LINK_DIRECTORIES(${PLASMA_DIR}/lib)
    SET(LIBS ${LIBS} ${PLASMA_LIB_NAMES})
ENDIF(HAVE_PLASMA)

IF(HAVE_LAPACK)
    IF(NEED_F2C)
        SET(LIBS ${LIBS} ${LAPACKBLAS_LIB_NAMES} ${F2C_LIB_NAME})
    ELSE(NEED_F2C)
        SET(LIBS ${LIBS} ${LAPACKBLAS_LIB_NAMES})
    ENDIF(NEED_F2C)
ENDIF(HAVE_LAPACK)

add_library(MKF SHARED ${SOURCES})
set_property(TARGET MKF PROPERTY POSITION_INDEPENDENT_CODE ON)
add_dependencies(MKF MAS)

cmrc_add_resource_library(insulation_standards ALIAS data::insulation_standards NAMESPACE insulationData src/data/insulation_standards/IEC_60664-1.json src/data/insulation_standards/IEC_60664-4.json src/data/insulation_standards/IEC_60664-5.json src/data/insulation_standards/IEC_62368-1.json src/data/insulation_standards/IEC_61558-1.json src/data/insulation_standards/IEC_61558-2-16.json src/data/insulation_standards/IEC_60335-1.json)
target_link_libraries(MKF ${LIBS})

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

message(${DATA_LIST})
cmrc_add_resource_library(data ALIAS data::data NAMESPACE data ${DATA_LIST})
target_link_libraries(MKF data::data)


target_link_libraries(MKF nlohmann_json::nlohmann_json matplot)

target_include_directories(MKF PUBLIC
     "${PROJECT_SOURCE_DIR}/src/"
     "${PROJECT_SOURCE_DIR}/src/advisers/"
     "${PROJECT_SOURCE_DIR}/src/constructive_models/"
     "${PROJECT_SOURCE_DIR}/src/converter_models/"
     "${PROJECT_SOURCE_DIR}/src/physical_models/"
     "${PROJECT_SOURCE_DIR}/src/processors/"
     "${PROJECT_SOURCE_DIR}/src/support/"
     )

add_executable(magnetic_adviser "src/advisers/MagneticAdviser.cpp")
target_link_libraries(magnetic_adviser MKF)

if(INCLUDE_MKF_TESTS)
    file(GLOB TESTS_SOURCES "tests/*.cpp" )
    add_executable(MKF_tests ${TESTS_SOURCES})
    target_link_libraries(MKF_tests MKF)
    target_link_libraries(MKF_tests UnitTest++)
endif(INCLUDE_MKF_TESTS)

if(INCLUDE_PYMKF)
        FetchContent_Declare(pybind11
                GIT_REPOSITORY https://github.com/pybind/pybind11.git)

        FetchContent_Declare(pybind11_json
                GIT_REPOSITORY https://github.com/pybind/pybind11_json.git)

        FetchContent_MakeAvailable( pybind11 pybind11_json)
        set_target_properties(pybind11_json PROPERTIES EXCLUDE_FROM_ALL True)
        add_custom_target(PyMASGeneration
                          /bin/echo "RUNNING PyMASGeneration"
                          DEPENDS "${MAS_DIRECTORY}/MAS.hpp")
        add_dependencies(MKF PyMASGeneration)
        include_directories("${CMAKE_BINARY_DIR}/_deps/pybind11-src/include/")
        include_directories("${CMAKE_BINARY_DIR}/_deps/pybind11_json-src/include/pybind11_json/")
        target_link_libraries(MKF pybind11_json pybind11::embed)
        include_directories("${CMAKE_BINARY_DIR}/_deps/pybind11-src")
        pybind11_add_module(PyMKF src/PyMKF/PyMKFWrapper.cpp )
        target_link_libraries(PyMKF PRIVATE MKF pybind11_json)

        # EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
        # define (VERSION_INFO) here.

        target_compile_definitions(PyMKF
                                   PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
endif(INCLUDE_PYMKF)