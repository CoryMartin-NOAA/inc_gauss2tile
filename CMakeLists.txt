cmake_minimum_required(VERSION 3.19)

project(interp_inc
        VERSION 1.0.0
        LANGUAGES C Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_DIRECTORY_LABELS ${PROJECT_NAME})

include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT CMAKE_C_COMPILER_ID MATCHES "^(GNU|Intel|Clang|AppleClang)$")
  message(WARNING "${CMAKE_C_COMPILER_ID} is not supported.")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Intel)$")
  message(WARNING "${CMAKE_Fortran_COMPILER_ID} is not supported.")
endif()

# User options
option(OPENMP "Enable OpenMP Threading" OFF)
option(ENABLE_MKL "Use MKL for LAPACK implementation (if available)" ON)

# Dependencies
if(ENABLE_MKL)
  find_package(MKL QUIET)
endif()
if(MKL_FOUND)
  set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
else()
  set(ENABLE_MKL OFF CACHE INTERNAL "GSI Utils: Disable MKL since it was NOT FOUND")
  find_package(LAPACK REQUIRED)
endif()
find_package(MPI REQUIRED)
find_package(NetCDF REQUIRED Fortran)
if(OPENMP)
  find_package(OpenMP REQUIRED)
endif()

add_subdirectory(parm)
add_subdirectory(src)

