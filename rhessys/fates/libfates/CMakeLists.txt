# Standalone build of libfates for the RHESSys host (ED/FATES design doc 12,
# Phase 0). Compiles the FATES science dirs at the pinned tag
# (FATES_SRC = NGEET/fates sci.1.93.2_api.46.0.0 + rhessys_host_patches.patch)
# plus minimal CESM shr_* and CLM seed-dispersal stand-ins (shims/), so no
# CIME/ESMF/csm_share is needed. The FATES source and the build live outside
# this repo, by default in <parent of RHESSys>/fates_lib:
#
#   cmake -S <RHESSys>/rhessys/fates/libfates -B <parent>/fates_lib/build \
#         [-DFATES_SRC=<parent>/fates_lib/src] && cmake --build <parent>/fates_lib/build
#   -> build/lib/libfates.a, build/mod/*.mod (the RHESSys makefile default FATES_DIR)
cmake_minimum_required(VERSION 3.16)
project(rhessys_fates Fortran C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

get_filename_component(FATES_SRC_DEFAULT ${CMAKE_CURRENT_SOURCE_DIR}/../../../../fates_lib/src ABSOLUTE)
set(FATES_SRC ${FATES_SRC_DEFAULT} CACHE PATH "patched FATES source (NGEET/fates at the pinned tag)")
if(NOT EXISTS ${FATES_SRC}/main/EDTypesMod.F90)
  message(FATAL_ERROR "no FATES source in ${FATES_SRC}; clone NGEET/fates there, apply rhessys_host_patches.patch, or set -DFATES_SRC")
endif()
file(GLOB FATES_SOURCES
  ${FATES_SRC}/main/*.F90
  ${FATES_SRC}/biogeochem/*.F90
  ${FATES_SRC}/biogeophys/*.F90
  ${FATES_SRC}/parteh/*.F90
  ${FATES_SRC}/fire/*.F90
  ${FATES_SRC}/radiation/*.F90)
file(GLOB SHIM_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/shims/*.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/shims/*.c)

add_library(fates STATIC ${SHIM_SOURCES} ${FATES_SOURCES})
target_include_directories(fates PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/shims/include)
set_target_properties(fates PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  Fortran_PREPROCESS ON)
target_compile_options(fates PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none -fallow-argument-mismatch -fPIC -fopenmp>)
# -fopenmp: RHESSys runs FATES sites in parallel (rhessys_veg_do_all); it also
# makes FATES's local work arrays automatic (thread-private) instead of static
