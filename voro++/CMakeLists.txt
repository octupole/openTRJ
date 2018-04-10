# Policy
cmake_minimum_required( VERSION 3.0) 

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

include_directories("include")
add_library(voro++ STATIC cell.cc common.cc container.cc container_prd.cc c_loops.cc pre_container.cc 
	unitcell.cc v_compute.cc v_base.cc wall.cc)
