# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE CenterMass.cpp CenterMassMicelles.cpp CenterMassBW3.cpp CenterMassWigner.cpp
	Quaternions.cpp QuaternionUtilities.cpp WignerDMatrices.cpp Combinatorics.cpp CMPick.cpp)
	
add_library(Properties STATIC ${SOURCE})
include_directories("./" "${PROJECT_SOURCE_DIR}/libtraj" "${PROJECT_SOURCE_DIR}/Molecules" "${PROJECT_SOURCE_DIR}/libtraj")
