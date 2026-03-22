cmake_minimum_required(VERSION 3.10)

project(MLCM_PROTON_BEAM)

# Use your own compiler only if absolutely necessary
#set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# Include pybind11 path
list(APPEND CMAKE_PREFIX_PATH "/home/alex_unix/.local/lib/python3.10/site-packages/pybind11/share/cmake/pybind11")

# Find pybind11
find_package(pybind11 REQUIRED)

# Ensure all code is position-independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Paths
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

# Common source files
set(COMMON_SOURCES
    "${SOURCE_DIR}/grids.cpp"
    "${SOURCE_DIR}/HU2material.cpp"
    "${SOURCE_DIR}/kinematics.cpp"
    "${SOURCE_DIR}/particles.cpp"
    "${SOURCE_DIR}/pencilbeam.cpp"
    "${SOURCE_DIR}/simulation.cpp"
    "${SOURCE_DIR}/voxeltraversal.cpp"
)

# Header files (not strictly needed in targets)
set(HEADERS
    "${INCLUDE_DIR}/physics_constants.hpp"
    "${INCLUDE_DIR}/simulation_constants.hpp"
    "${INCLUDE_DIR}/error_logging.hpp"
    "${INCLUDE_DIR}/grids.hpp"
    "${INCLUDE_DIR}/HU2material.hpp"
    "${INCLUDE_DIR}/kinematics.hpp"
    "${INCLUDE_DIR}/particles.hpp"
    "${INCLUDE_DIR}/pencilbeam.hpp"
    "${INCLUDE_DIR}/simulation.hpp"
    "${INCLUDE_DIR}/voxeltraversal.hpp"
)

# -------------------
# Shared C++ library
# -------------------
add_library(mlmc_core STATIC ${COMMON_SOURCES})
target_include_directories(mlmc_core PUBLIC ${INCLUDE_DIR})

# -------------------
# Pybind11 Python module
# -------------------
pybind11_add_module(mlmcprotons "${SOURCE_DIR}/bindings.cpp")
target_include_directories(mlmcprotons PRIVATE ${INCLUDE_DIR})
target_link_libraries(mlmcprotons PRIVATE mlmc_core pybind11::module)

set_target_properties(mlmcprotons PROPERTIES
    OUTPUT_NAME "mlmcprotons"
    PREFIX ""
    SUFFIX ".so"
)

# Optional install step
install(TARGETS mlmcprotons
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/site-packages
)

# -------------------
# Native C++ executable
# -------------------
add_executable(MLCM_PROTON_BEAM "${SOURCE_DIR}/main.cpp")
target_include_directories(MLCM_PROTON_BEAM PRIVATE ${INCLUDE_DIR})
target_link_libraries(MLCM_PROTON_BEAM PRIVATE mlmc_core)
