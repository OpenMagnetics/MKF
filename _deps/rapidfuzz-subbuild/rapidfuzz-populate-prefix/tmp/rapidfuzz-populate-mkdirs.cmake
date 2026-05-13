# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/alf/OpenMagnetics/MKF/_deps/rapidfuzz-src")
  file(MAKE_DIRECTORY "/home/alf/OpenMagnetics/MKF/_deps/rapidfuzz-src")
endif()
file(MAKE_DIRECTORY
  "/home/alf/OpenMagnetics/MKF/_deps/rapidfuzz-build"
  "/home/alf/OpenMagnetics/MKF/_deps/rapidfuzz-subbuild/rapidfuzz-populate-prefix"
  "/home/alf/OpenMagnetics/MKF/_deps/rapidfuzz-subbuild/rapidfuzz-populate-prefix/tmp"
  "/home/alf/OpenMagnetics/MKF/_deps/rapidfuzz-subbuild/rapidfuzz-populate-prefix/src/rapidfuzz-populate-stamp"
  "/home/alf/OpenMagnetics/MKF/_deps/rapidfuzz-subbuild/rapidfuzz-populate-prefix/src"
  "/home/alf/OpenMagnetics/MKF/_deps/rapidfuzz-subbuild/rapidfuzz-populate-prefix/src/rapidfuzz-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/alf/OpenMagnetics/MKF/_deps/rapidfuzz-subbuild/rapidfuzz-populate-prefix/src/rapidfuzz-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/alf/OpenMagnetics/MKF/_deps/rapidfuzz-subbuild/rapidfuzz-populate-prefix/src/rapidfuzz-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
