cmake_minimum_required(VERSION 3.18)
project(BRIXS LANGUAGES Fortran)

# BRIXS: BSE Calculations for RIXS spectra
#
# Build modes:
#   -DCMAKE_BUILD_TYPE=Release (default)  -> optimized build
#   -DCMAKE_BUILD_TYPE=Debug              -> -g plus compiler runtime checks
#   -DBRIXS_ENABLE_MPI=ON                 -> link MPI and an MPI-enabled (parallel) HDF5
#
# See README.md for the full build instructions.

option(BRIXS_ENABLE_MPI "Build with MPI / parallel HDF5 support" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Debug or Release)" FORCE)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # keep the historical install destination (BRIXS/bin) unless the user overrides it
  set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH "Install prefix" FORCE)
endif()

find_package(HDF5 REQUIRED COMPONENTS Fortran)

if(BRIXS_ENABLE_MPI)
  if(DEFINED HDF5_IS_PARALLEL AND NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR
      "BRIXS_ENABLE_MPI=ON requires an MPI-enabled (parallel) HDF5 installation, "
      "but the HDF5 found (include dir: ${HDF5_Fortran_INCLUDE_DIRS}) was built "
      "without MPI support. Point HDF5_ROOT at a parallel HDF5 install.")
  endif()
  find_package(MPI REQUIRED COMPONENTS Fortran)
endif()

find_package(LAPACK REQUIRED)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbacktrace -fcheck=all")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback -check all")
endif()

add_subdirectory(src)

enable_testing()
add_subdirectory(test)

find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
  add_custom_target(docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/src/doxygen.config
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM
  )
endif()
