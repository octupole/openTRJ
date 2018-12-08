# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE Voronoi.cpp  VoronoiMicelles.cpp VoronoiSetter.cpp
	VoronoiBinary.cpp IteratorVoronoi.cpp VoronoiMicellesJSON.cpp)
add_library(Voronoi STATIC ${SOURCE})
target_include_directories(Voronoi PUBLIC "${PROJECT_SOURCE_DIR}/Molecules" "${PROJECT_SOURCE_DIR}/Properties" 
	"${PROJECT_SOURCE_DIR}/libtraj")
