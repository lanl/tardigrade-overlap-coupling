#===================================================================================================== PROJECT SETUP ===
cmake_minimum_required(VERSION 3.14)
project(overlap_coupling)
set(PROJECT_LIBRARY_NAME "overlap")
set(PROJECT_SOURCE_NAME "overlapCoupling")
set(FILTER_NAME "filter")

# Set flag for building the python bindings
set(OVERLAP_COUPLING_BUILD_PYTHON_BINDINGS ON CACHE BOOL "Flag for whether to build the python bindings or not")

# Set common project paths relative to project root directory
set(CPP_SRC_PATH "src/cpp")
set(CPP_TEST_PATH "${CPP_SRC_PATH}/tests")
set(CMAKE_SRC_PATH "src/cmake")
set(DOXYGEN_SRC_PATH "docs/doxygen")
set(SPHINX_SRC_PATH "docs/sphinx")

# Add the cmake folder to locate project CMake module(s)
set(YAML-CPP_INSTALL_DIR "/path/to/yaml-cpp-config.cmake" CACHE PATH "The location of the yaml-cpp-config.cmake or yamlConfig.cmake file")
set(XDMF_INSTALL_DIR "/path/to/Xdmf-config.cmake" CACHE PATH "The location of the Xdmf-config.cmake or XdmfConfig.cmake file")
#set(VORO_INSTALL_DIR "/path/to/voro-config.cmake" CACHE PATH "The location of the voro-config.cmake or VOROConfig.cmake file")

set(ADDITIONAL_INSTALL_DIRECTORIES ${YAML-CPP_INSTALL_DIR} ${XDMF_INSTALL_DIR} ${VORO_INSTALL_DIR})
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${ADDITIONAL_INSTALL_DIRECTORIES})

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}" ${CMAKE_MODULE_PATH})

# Get version number from Git
set(VERSION_UPDATE_FROM_GIT True)
# FIXME: Figure out why include() doesn't pick up CMAKE_MODULE_PATH correctly
include(${CMAKE_SRC_PATH}/GetVersionFromGitTag.cmake)
project(${PROJECT_NAME} VERSION ${${PROJECT_NAME}_VERSION})

# Add installation directory variables
include(GNUInstallDirs)

# Make the code position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set the c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -ansi -pedantic -lrt")

# Enable CTest
enable_testing()

#================================================================================================= FIND DEPENDENCIES ===
# Find eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
if(EIGEN3_FOUND)
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
endif()

# Find yaml-cpp
find_package(YAML-CPP REQUIRED NO_MODULE)
include_directories(${YAML-CPP_INCLUDE_DIR})
if(YAML-CPP_FOUND)
    message(STATUS "Found yaml-cpp: ${YAML-CPP_INCLUDE_DIR}")
endif()

# Find Xdmf
find_package(xdmf CONFIG REQUIRED
                         NAMES "Xdmf" "xdmf")
include_directories(${XDMF_INCLUDE_DIRS})
set(XDMF_LIBRARY_DIRS ${XDMF_CMAKE_INSTALL_FULL_LIBDIR} "${XDMF_CMAKE_INSTALL_FULL_LIBDIR}/../lib")
link_directories(${XDMF_LIBRARY_DIRS})
if(xdmf_FOUND)
    message(STATUS "Found Xdmf: ${XDMF_INCLUDE_DIRS}")
endif()
set(XDMF_LIBRARIES "Xdmf" "XdmfCore" "XdmfUtils")

find_package(LibXml2 REQUIRED)
include_directories(${LIBXML2_INCLUDE_DIR})
if(LIBXML2_FOUND)
    message(STATUS "Found LibXml2: ${LIBXML2_INCLUDE_DIR}")
endif()


#find_package(VORO REQUIRED)
#include_directories(${VORO_INCLUDE_DIR})
#if(VORO_FOUND)
#    message(STATUS "Found voro++: ${VORO_INCLUDE_DIR}")
#endif()

# Set build type checks
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lower)
set(upstream_required "")
if(cmake_build_type_lower STREQUAL "release")
    set(upstream_required "REQUIRED")
endif()

# Set the related projects for the compile
set(LOCAL_SUPPORT_MODULES "dataFileInterface" "generateXDMFData" "element" "assembly" "DOFProjection" "geometry_decomposition" "inputFileProcessor" "volumeReconstruction")

# Set the required libraries for each of the support modules
set(dataFileInterface_SUPPORT_LIBS "yaml-cpp" ${XDMF_LIBRARIES} "xml2")
set(generateXDMFData_SUPPORT_LIBS "dataFileInterface")
set(element_SUPPORT_LIBS "error_tools")
set(assembly_SUPPORT_LIBS "")
set(DOFProjection_SUPPORT_LIBS "")
set(geometry_decomposition_SUPPORT_LIBS "")
set(inputFileProcessor_SUPPORT_LIBS "dataFileInterface")
set(volumeReconstruction_SUPPORT_LIBS "element" "yaml-cpp" ${XDMF_LIBRARIES}  "xml2")

set(test_dataFileInterface_SUPPORT_LIBS "")
set(test_generateXDMFData_SUPPORT_LIBS "")
set(test_element_SUPPORT_LIBS "")
set(test_assembly_SUPPORT_LIBS "")
set(test_DOFProjection_SUPPORT_LIBS "")
set(test_geometry_decomposition_SUPPORT_LIBS "")
set(test_inputFileProcessor_SUPPORT_LIBS "generateXDMFData")
set(test_volumeReconstruction_SUPPORT_LIBS "")

# Find related, but external, projects in installed environment
include(FetchContent)
set(upstream_packages "error_tools" "vector_tools" "solver_tools" "micromorphic_tools" "micromorphic_element")
foreach(package ${upstream_packages})
    string(TOUPPER "${package}" package_upper)
    set(${package_upper}_BUILD_PYTHON_BINDINGS ${OVERLAP_COUPLING_BUILD_PYTHON_BINDINGS})
    find_package(${package} ${upstream_required} CONFIG)
    if(${package}_FOUND)
        message(STATUS "Found ${package}: ${${package}_DIR}")
    else()
        # Find related, but external, projects using FetchContent and building locally
        message(WARNING "Did not find an installed ${package} package. Attempting local build with FetchContent.")
        # TODO: Decide how to handle main/dev build differences. bash env variables for GIT TAG?
        FetchContent_Declare(
            ${package}
            GIT_REPOSITORY ssh://git@re-git.lanl.gov:10022/aea/material-models/${package}.git
            GIT_TAG        origin/dev
        )
        FetchContent_MakeAvailable(${package})
    endif()
endforeach(package)

## Fetch additional header files
#FetchContent_Declare(
#    quickhull
#    GIT_REPOSITORY https://github.com/akuukka/quickhull.git
#    GIT_TAG        origin/master
#)
#FetchContent_MakeAvailable(quickhull)
#message("quickhull found: ${quickhull_SOURCE_DIR}")
#include_directories(${quickhull_SOURCE_DIR})

#=============================================================================================== ADD PROJECT TARGETS ===
# MUST COME AFTER DEPENDENCY LOCATING
# Add project source directories
include_directories(${CPP_SRC_PATH})
add_subdirectory(${CPP_SRC_PATH})

# Only add tests and documentation for current project builds. Protects downstream project builds.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Find Boost. Required for tests
    find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)
    # Add c++ tests and docs
    add_subdirectory(${CPP_TEST_PATH})
    add_subdirectory(${DOXYGEN_SRC_PATH})
    add_subdirectory(${SPHINX_SRC_PATH})
endif()

#==================================================================================== SETUP INSTALLATION CMAKE FILES ===
foreach(package ${PROJECT_LIBRARY_NAME} ${FILTER_NAME})
    include(CMakePackageConfigHelpers)
    write_basic_package_version_file("${package}ConfigVersion.cmake"
                                     VERSION ${PROJECT_VERSION}
                                     COMPATIBILITY SameMajorVersion)
    configure_package_config_file(
      "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}/Config.cmake.in"
      "${PROJECT_BINARY_DIR}/${package}Config.cmake"
      INSTALL_DESTINATION
      ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)

    # CMake won't build the targets for local builds of upstream projects
    if(cmake_build_type_lower STREQUAL release)
        install(EXPORT ${package}_Targets
                FILE ${package}Targets.cmake
                DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)
    endif()

    install(FILES "${PROJECT_BINARY_DIR}/${package}Config.cmake"
                  "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
            DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)
endforeach(package)
