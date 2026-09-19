# -----------------------------------------------------------------------------
# Preliminary Settings
# -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.16)

# Set policies for modern CMake behavior
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)  # Use BoostConfig.cmake if available
endif()

project(
  pySMOKEPostProcessor
  VERSION 1.0.0
  DESCRIPTION "Python wrapper of the OpenSMOKEpp Graphical Post Processor"
  LANGUAGES CXX
)

# -----------------------------------------------------------------------------
# Set module path in order to use custom CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# -----------------------------------------------------------------------------
# C++ Standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# -----------------------------------------------------------------------------
# Standalone cxx_wrapper/testfiles/*.cpp drivers: bypass pybind entirely, fast
# iteration loop for the C++ backend. Off by default so it never affects the
# scikit-build-core wheel/pip path; enable with -DPYSMOKE_BUILD_TESTFILES=ON.
option(PYSMOKE_BUILD_TESTFILES "Build standalone testfiles/*.cpp drivers" OFF)

# -----------------------------------------------------------------------------
# Detect if building via scikit-build-core
if(DEFINED SKBUILD)
  message(STATUS "Building via scikit-build-core")
  set(SKBUILD_BUILD ON)
  # scikit-build sets installation directory
  set(CMAKE_INSTALL_PREFIX ${SKBUILD_PLATLIB_DIR}/pySMOKEPostProcessor)
else()
  message(STATUS "Building via standard CMake")
  set(SKBUILD_BUILD OFF)
  # For manual builds, install to project directory
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(
      CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/pySMOKEPostProcessor"
      CACHE PATH "Install prefix" FORCE
    )
  endif()
endif()

# -----------------------------------------------------------------------------
# Compulsory libraries
# -----------------------------------------------------------------------------
# Boost (Required by OpenSMOKEpp) Inside the installation to ease the life of
# everyone the boost package is installed using brew/apt/pip/conda
message(STATUS "Searching for Boost...")

# -----------------------------------------------------------------------------
# Try to find Boost with modern config first (Boost 1.70+)
# If not found, fall back to FindBoost module for older versions
# (e.g., manylinux boost169)
#find_package(Boost 1.70 QUIET CONFIG)
find_package(Boost 1.70 REQUIRED COMPONENTS filesystem system)

if(NOT Boost_FOUND)
  # Fallback to FindBoost module for older versions (Boost < 1.70)
  # This is needed for manylinux which uses boost169 (Boost 1.69)
  find_package(Boost REQUIRED
      COMPONENTS
          date_time
          filesystem
          program_options
          system
          regex
          timer
          chrono
  )
endif()

if(Boost_FOUND)
  message(STATUS "Boost found: ${Boost_VERSION}")
  if(Boost_INCLUDE_DIRS)
    message(STATUS "Boost include: ${Boost_INCLUDE_DIRS}") # This is the correct one
  elseif(Boost_INCLUDE_DIR)
    message(STATUS "Boost include: ${Boost_INCLUDE_DIR}")
  endif()
endif()

# -----------------------------------------------------------------------------
# External dependencies (FetchContent)
# -----------------------------------------------------------------------------
# Fetch Eigen3 (with system fallback)
include(cmake/external/Eigen3.cmake)

# -----------------------------------------------------------------------------
# Fetch pybind11 (with system fallback)
include(cmake/external/pybind11.cmake)

# -----------------------------------------------------------------------------
# Fetch OpenSMOKEpp v0.22
include(cmake/external/OpenSMOKEpp.cmake)

# -----------------------------------------------------------------------------
# Source code subdirectories
# -----------------------------------------------------------------------------
add_subdirectory(cxx_wrapper)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

set(LIBRARIES ${Boost_LIBRARIES}
              ${Python3_LIBRARIES}
              )
set(INCLUDE ${Boost_INCLUDE_DIRS}
            ${EIGEN3_INCLUDE_DIR}
            ${OPENSMOKE_INCLUDE}
            ${Python3_INCLUDE_DIRS}
            )

# testfiles/*.cpp executable targets (PYSMOKE_BUILD_TESTFILES) are defined in
# cxx_wrapper/CMakeLists.txt, next to the main module target whose Boost/Eigen
# linking logic they reuse.