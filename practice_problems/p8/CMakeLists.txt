# Minimum CMake version required
cmake_minimum_required(VERSION 3.14)

# Project name
project(p8 LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Specify build type if not set
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Find deal.II
find_package(deal.II 9.3.0 REQUIRED HINTS /path/to/dealii)

# Declare source files
set(TARGET_SRC
    ${PROJECT_NAME}.cc
)

# Create executable
add_executable(${PROJECT_NAME} ${TARGET_SRC})

# Set up deal.II target
deal_II_setup_target(${PROJECT_NAME})
