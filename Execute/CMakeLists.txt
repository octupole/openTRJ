# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
##set(SOURCE ExecbSaxs.cpp ExecbSaxsPost.cpp ExecbSaxsTraj.cpp ExecuteVoronoi.cpp TrjRead.cpp  trjInput.cpp ClearUsage.cpp 
##	ExecuteProp.cpp trjInput.cpp ClearUsage.cpp)
set(SOURCE ExecuteVoronoi.cpp TrjRead.cpp  trjInput.cpp ClearUsage.cpp 
	ExecuteProp.cpp trjInput.cpp ClearUsage.cpp)
set(SOURCE_SAXS Saxs/TrjReadSaxs.cpp Saxs/trjInputSaxs.cpp Saxs/ExecbSaxs.cpp Saxs/ExecbSaxsPost.cpp Saxs/ExecbSaxsTraj.cpp)
	
add_library(Execute STATIC ${SOURCE})
add_library(ExecuteSaxs STATIC ${SOURCE_SAXS})

include_directories("./" "${PROJECT_SOURCE_DIR}/libtraj" "${PROJECT_SOURCE_DIR}/Molecules" "${PROJECT_SOURCE_DIR}/Voronoi" "${PROJECT_SOURCE_DIR}/Properties" "${PROJECT_SOURCE_DIR}/Saxs" )

