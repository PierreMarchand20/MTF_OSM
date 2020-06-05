#=============================================================================#
#============================== Project ======================================#
#=============================================================================#
cmake_minimum_required(VERSION 3.0)
Project(MTF_OSM)

# To force c++11
if (${CMAKE_VERSION} VERSION_LESS 3.1)
	add_compile_options(-std=c++11)
elseif (${CMAKE_VERSION} VERSION_LESS 3.6.3 AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
	add_compile_options(-std=c++11)
else()
	set(CMAKE_CXX_STANDARD 11)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif ()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release)
endif()

# Files to do find_package for some module
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake_modules")

# Options
option(BEMTOOL_WITH_EXAMPLES "Build bemtool examples ?" ON)
option(BEMTOOL_WITH_TESTS    "Build bemtool tests ?" ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#=============================================================================#
#========================== External Libraries ===============================#
#=============================================================================#
#Bemtool
find_package(BemTool REQUIRED)

# Boost package - special functions
find_package(Boost REQUIRED)

# # EIGEN3
# find_package(Eigen)
# message("-- Eigen include files found in " "${Eigen_INCLUDE_DIRS}")

#=============================================================================#
#=========================== Libraries =======================================#
#=============================================================================#
add_library(bemtool INTERFACE)
target_include_directories(bemtool INTERFACE
					${BEMTOOL_INCLUDE_DIRS}
					${Boost_INCLUDE_DIRS})


# For headers to show in IDE
if (NOT "${CMAKE_VERSION}" VERSION_LESS 3.1)
	target_sources(bemtool INTERFACE 
		${BEMTOOL_INCLUDE_DIR}
	)
endif()

#=============================================================================#
#========================== Repertories ======================================#
#=============================================================================#
add_subdirectory(src)