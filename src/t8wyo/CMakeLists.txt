CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(t8wyo)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

macro(SET_O2_FLAGS)
  string(REPLACE "-O3" "-O2" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
  string(REPLACE "-O3" "-O2" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
  string(REPLACE "-O3" "-O2" CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}")
endmacro()

# AVX vector compiler flags
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx2 -mfma")

# override release -O3 flag
#SET_O2_FLAGS()

# set C++17 standard for OCCA compatibility
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# --------------------------------------------------------------------------- #
# force GNU linker to include all shared libraries
if (NOT APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-as-needed")
else ()
  set(CMAKE_MACOSX_RPATH ON)
endif ()

set(t8code_dir "" CACHE PATH  "t8code install directory location")
set(metis_dir  "" CACHE PATH  "metis  install directory location")
# --------------------------------------------------------------------------- #

# ========================================== #
# Intel Optimization Flag: AVX2 instructions
if (CMAKE_C_COMPILER_ID MATCHES "Intel")
  add_definitions( -xCore-AVX2 )
endif ()
# ========================================== #

# ======================================= #
# Library header files to be included
include_directories(${t8code_dir}/include
                    ${metis_dir}/include
                    ./include)

# Libraries to be linked
find_library(t8code_lb  NAMES t8    PATHS ${t8code_dir}/lib)
find_library(sc_lb      NAMES sc    PATHS ${t8code_dir}/lib)
find_library(metis_lb   NAMES metis PATHS ${metis_dir}/lib)
# ======================================= #

# ======================================= #
# Source directories for this project:
#   recursively build subdirectory
set(src_path src)
add_subdirectory(${src_path})
# ======================================= #