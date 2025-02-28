CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

# Set DEAL_II_DIR explicitly to the location of deal.II
SET(DEAL_II_DIR "${CMAKE_HOME_DIRECTORY}/dealii")

# Print deal.II directory for debugging
MESSAGE(STATUS "Looking for deal.II in: ${DEAL_II_DIR}")

# Find deal.II package
FIND_PACKAGE(deal.II 9.0 REQUIRED HINTS ${DEAL_II_DIR})

# Check if deal.II is found
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** could not locate deal.II. ***\n\n"
    "you may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
)
ENDIF()

# Set the project name (not the target)
PROJECT(FEM_Practice)


# Set C++ standard to C++17 (or C++20 if you prefer)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # Optional: Disables compiler-specific extensions

# Add compiler flags for advanced C++ features
add_compile_options(-std=c++17)
add_compile_options(-fexceptions) # Allow exceptions for debugging


# Initialize deal.II cached variables
DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Add the practice_problems subdirectory
ADD_SUBDIRECTORY(practice_problems)

SET(CMAKE_EXPORT_COMPILE_COMMANDS OFF)

# Clean up CMake install files (only if they exist)

execute_process(
    COMMAND ${CMAKE_COMMAND} -E remove_if_exists ${CMAKE_BINARY_DIR}/cmake_install.cmake
)
execute_process(
    COMMAND ${CMAKE_COMMAND} -E remove_if_exists ${CMAKE_BINARY_DIR}/install_manifest.txt
)
execute_process(
    COMMAND ${CMAKE_COMMAND} -E remove_if_exists ${CMAKE_BINARY_DIR}/CMakeCache.txt
)


