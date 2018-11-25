# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(SOURCE CenterMass.cpp CenterMassMicelles.cpp CenterMassBW3.cpp CenterMassWigner.cpp
	Quaternions.cpp QuaternionUtilities.cpp WignerDMatrices.cpp Combinatorics.cpp 
	CMPick.cpp DensMode.cpp)
	
add_library(Properties STATIC ${SOURCE})
