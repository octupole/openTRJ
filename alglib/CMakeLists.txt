# Policy
cmake_minimum_required( VERSION 3.0) 

# Project name
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(SOURCE alglibinternal.cpp alglibmisc.cpp ap.cpp dataanalysis.cpp diffequations.cpp fasttransforms.cpp
	integration.cpp interpolation.cpp linalg.cpp optimization.cpp solvers.cpp specialfunctions.cpp
	statistics.cpp)
include_directories("include")
add_library(alglib STATIC ${SOURCE} )
