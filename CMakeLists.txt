cmake_minimum_required(VERSION 3.15)

# Set macOS SDK before project() is called
if(APPLE)
    execute_process(
        COMMAND xcrun --show-sdk-path
        OUTPUT_VARIABLE MACOS_SDK_PATH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(CMAKE_OSX_SYSROOT "${MACOS_SDK_PATH}" CACHE PATH "macOS SDK path" FORCE)
endif()

project(gctb)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Ensure we use the C++ standard library
if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

# Find required packages
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Find pybind11 - try via Python first
execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE pybind11_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
find_package(pybind11 CONFIG REQUIRED PATHS "${pybind11_DIR}")

# Find Eigen3
find_package(Eigen3 QUIET)
if(NOT Eigen3_FOUND)
    if(DEFINED ENV{EIGEN3_INCLUDE_DIR})
        set(EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
    elseif(APPLE)
        set(EIGEN3_INCLUDE_DIR "/opt/homebrew/include/eigen3")
    endif()
endif()

# Find Boost
find_package(Boost QUIET)
if(NOT Boost_FOUND)
    if(DEFINED ENV{BOOST_LIB})
        set(Boost_INCLUDE_DIRS $ENV{BOOST_LIB})
    elseif(APPLE)
        set(Boost_INCLUDE_DIRS "/opt/homebrew/include")
    endif()
endif()

# Find OpenMP
if(APPLE)
    # On macOS, use libomp from Homebrew
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
endif()
find_package(OpenMP REQUIRED)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
)

# Core C++ source files
set(CORE_SOURCES
    src/data.cpp
    src/gadgets.cpp
    src/mcmc.cpp
    src/model.cpp
    src/stat.cpp
    src/stratify.cpp
    # src/stratifyMixture.cpp  # Temporarily removed - uses MPI
    src/vgmaf.cpp
    src/xci.cpp
    src/eigen.cpp
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

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(_core PRIVATE 
        -O3 
        -Wall 
        -w
        -DUNIX
        -DNDEBUG
    )
    
    # Add SSE2 only for x86_64 architecture
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
        target_compile_options(_core PRIVATE -msse2 -m64)
    endif()
    
    if(APPLE)
        target_compile_definitions(_core PRIVATE Dfopen64=fopen)
    endif()
endif()

