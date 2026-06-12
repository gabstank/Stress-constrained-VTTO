# 1. Initialization
cmake_minimum_required(VERSION 3.16)

# 2. Project Definition
project(topopt CXX)

# 3. Set Build Type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "No build type selected, defaulting to DebugRelease")
    set(CMAKE_BUILD_TYPE "DebugRelease" CACHE STRING "Build type (Debug, Release, DebugRelease)" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "DebugRelease")
endif()

# 4. C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 5. Find Dependencies
# A. Find deal.II
find_package(deal.II 9.4.0 REQUIRED
    HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
)

# B. Find MPI (Crucial: This fills the variables we use later)
find_package(MPI REQUIRED)

# 6. Check deal.II Configuration
if(NOT (DEAL_II_WITH_PETSC OR DEAL_II_WITH_TRILINOS))
    message(FATAL_ERROR "*** Error! This software requires deal.II with PETSc or Trilinos. ***")
endif()

# 7. Configure Project Options
option(WITH_PETSC "Use PETSc for Linear Algebra (if OFF, uses Trilinos)" OFF)

if(WITH_PETSC)
    if(NOT DEAL_II_WITH_PETSC)
        message(FATAL_ERROR "You requested PETSc, but deal.II was not configured with DEAL_II_WITH_PETSC=ON.")
    endif()
    message(STATUS "Configured for use of PETSC for LA")
else()
    if(NOT DEAL_II_WITH_TRILINOS)
        message(FATAL_ERROR "You requested Trilinos, but deal.II was not configured with DEAL_II_WITH_TRILINOS=ON.")
    endif()
    message(STATUS "Configured for use of Trilinos for LA")
endif()

# 8. Gather Source Files
file(GLOB_RECURSE TARGET_SRC "src/*.cc" "main.cc")

# 9. Create the Executable Target
add_executable(${PROJECT_NAME} ${TARGET_SRC})

# 10. Configure the Target

# A. Link deal.II
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
if (CMAKE_BUILD_TYPE_UPPER MATCHES "DEBUG" AND NOT CMAKE_BUILD_TYPE_UPPER MATCHES "DEBUGRELEASE")
    DEAL_II_SETUP_TARGET(${PROJECT_NAME} DEBUG)
else()
    DEAL_II_SETUP_TARGET(${PROJECT_NAME} RELEASE)
endif()

# B. Link MPI (The "Well-Designed" Fix)
target_link_libraries(${PROJECT_NAME} MPI::MPI_CXX)

# Force BFD linker (Solve "gold" linker issues)
target_link_options(${PROJECT_NAME} PRIVATE "-fuse-ld=bfd")

# Dynamic Path Extraction (Solve "cannot find -lopen-pal")
if(MPI_CXX_LIBRARIES)
    list(GET MPI_CXX_LIBRARIES 0 FIRST_MPI_LIB)
    if(EXISTS "${FIRST_MPI_LIB}")
        get_filename_component(MPI_LIB_DIR "${FIRST_MPI_LIB}" DIRECTORY)
        target_link_directories(${PROJECT_NAME} PRIVATE ${MPI_LIB_DIR})
    endif()
endif()

target_include_directories(${PROJECT_NAME} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/BVP
        ${CMAKE_SOURCE_DIR}/src/BVP/Linear
        ${CMAKE_SOURCE_DIR}/src/DesignField
        ${CMAKE_SOURCE_DIR}/src/Filters/Core
        ${CMAKE_SOURCE_DIR}/src/Filters/Projection
        ${CMAKE_SOURCE_DIR}/src/Filters/HeavisideKernels
        ${CMAKE_SOURCE_DIR}/src/Filters/Postprocessing
        ${CMAKE_SOURCE_DIR}/src/Materials
        ${CMAKE_SOURCE_DIR}/src/Mesh
        ${CMAKE_SOURCE_DIR}/src/Optimizer
        ${CMAKE_SOURCE_DIR}/src/Responses
        ${CMAKE_SOURCE_DIR}/src/Responses/Linear
        ${CMAKE_SOURCE_DIR}/src/Utilities
)

if(WITH_PETSC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE USE_PETSC)
else()
    target_compile_definitions(${PROJECT_NAME} PRIVATE USE_TRILINOS)
endif()

# 11. Cleanup
set(CLEAN_UP_FILES *.log *.gmv *.gnuplot *.gpl *.eps *.pov *.vtk *.ucd *.d2 *.vtu *.pvtu *.pvd)