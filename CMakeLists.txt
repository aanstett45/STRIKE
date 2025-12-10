cmake_minimum_required(VERSION 3.15)
project(NZIM_Experiments CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options de compilation
option(WITH_OPENMP "Enable OpenMP" ON)
option(WITH_ROARING "Enable Roaring Bitmaps" ON)
option(WITH_TORCH "Enable PyTorch" OFF)

# Flags d'optimisation
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# OpenMP
if(WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# Roaring Bitmaps
if(WITH_ROARING)
    find_package(roaring REQUIRED)
    find_path(ROARING_INCLUDE roaring/roaring.h REQUIRED)
    message(STATUS "Roaring library: ${ROARING_LIB}")
    message(STATUS "Roaring include: ${ROARING_INCLUDE}")
    add_definitions(-DWITH_ROARING)
endif()

# PyTorch (optionnel)
if(WITH_TORCH)
    find_package(Torch REQUIRED)
    message(STATUS "Torch found: ${TORCH_VERSION}")
    add_definitions(-DWITH_TORCH)
endif()

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${ROARING_INCLUDE}
)

# Source files
set(MATRICE_SOURCES
    Matrice.cpp
    Papier.cpp
)

set(EXPERIMENT_SOURCES
    nzim_experiments.cpp
)

# Libraries
set(EXPERIMENT_LIBS
    ${ROARING_LIB}
)

if(WITH_OPENMP)
    list(APPEND EXPERIMENT_LIBS OpenMP::OpenMP_CXX)
endif()

if(WITH_TORCH)
    list(APPEND EXPERIMENT_LIBS "${TORCH_LIBRARIES}")
endif()

# Executable: Main experiments
add_executable(nzim_experiments ${EXPERIMENT_SOURCES} ${MATRICE_SOURCES})
target_link_libraries(nzim_experiments ${EXPERIMENT_LIBS})

# Executable: Beta calculator (standalone tool)
add_executable(compute_beta compute_beta_tool.cpp ${MATRICE_SOURCES})
target_link_libraries(compute_beta ${EXPERIMENT_LIBS})

# Installation
install(TARGETS nzim_experiments compute_beta
        RUNTIME DESTINATION bin)

# Tests (optionnel)
enable_testing()

add_test(NAME test_beta_zero
         COMMAND nzim_experiments
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Summary
message(STATUS "")
message(STATUS "========================================")
message(STATUS "NZIM SpGEMM Experimental Framework")
message(STATUS "========================================")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "OpenMP: ${WITH_OPENMP}")
message(STATUS "Roaring: ${WITH_ROARING}")
message(STATUS "PyTorch: ${WITH_TORCH}")
message(STATUS "========================================")
message(STATUS "")
