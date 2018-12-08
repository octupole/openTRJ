# Policy
cmake_minimum_required( VERSION 3.0) 


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
##set(SOURCE ExecbSaxs.cpp ExecbSaxsPost.cpp ExecbSaxsTraj.cpp ExecuteVoronoi.cpp TrjRead.cpp  trjInput.cpp ClearUsage.cpp 
##	ExecuteProp.cpp trjInput.cpp ClearUsage.cpp)
set(SOURCE ExecuteVoronoi.cpp TrjRead.cpp  trjInput.cpp ClearUsage.cpp 
	ExecuteProp.cpp trjInput.cpp ClearUsage.cpp)
	
set(SOURCE_SAXS Saxs/TrjRead.cpp Saxs/trjInput.cpp Saxs/ExecbSaxs.cpp Saxs/ExecbSaxsPost.cpp 
	Saxs/ExecbSaxsTraj.cpp)
set(SOURCE_FIT FitSaxs/TrjRead.cpp FitSaxs/trjInput.cpp)
	
add_library(Execute STATIC ${SOURCE})
if(FFTW3_FOUND)
	add_library(ExecuteSaxs STATIC ${SOURCE_SAXS})
	add_library(ExecuteFit STATIC ${SOURCE_FIT})
endif()
target_include_directories(Execute PUBLIC "${PROJECT_SOURCE_DIR}/Properties" "${PROJECT_SOURCE_DIR}/Molecules" 
	"${PROJECT_SOURCE_DIR}/Voronoi" "${PROJECT_SOURCE_DIR}/libtraj")
target_include_directories(ExecuteSaxs PUBLIC "${PROJECT_SOURCE_DIR}/Properties" "${PROJECT_SOURCE_DIR}/Molecules" 
	"${PROJECT_SOURCE_DIR}/Saxs" "${PROJECT_SOURCE_DIR}/libtraj")
	