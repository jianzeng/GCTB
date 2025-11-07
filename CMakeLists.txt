cmake_minimum_required(VERSION 3.15)
project(gctb)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenMP REQUIRED)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/scr
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
)

# Core C++ source files
set(CORE_SOURCES
    scr/data.cpp
    scr/gadgets.cpp
    scr/gctb.cpp
    scr/hsq.cpp
    scr/mcmc.cpp
    scr/model.cpp
    scr/options.cpp
    scr/predict.cpp
    scr/stat.cpp
    scr/stratify.cpp
    scr/stratifyMixture.cpp
    scr/vgmaf.cpp
    scr/xci.cpp
    scr/eigen.cpp
)

# Create Python module
pybind11_add_module(_core 
    python/bindings/bindings.cpp
    ${CORE_SOURCES}
)

# Link libraries
target_link_libraries(_core PRIVATE 
    Eigen3::Eigen
    Boost::boost
    OpenMP::OpenMP_CXX
)

# Set output directory
set_target_properties(_core PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/gctb
)

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(_core PRIVATE 
        -O3 
        -Wall 
        -w
        -DUNIX
        -DNDEBUG
        -msse2
        -m64
    )
    if(APPLE)
        target_compile_definitions(_core PRIVATE Dfopen64=fopen)
    endif()
endif()

