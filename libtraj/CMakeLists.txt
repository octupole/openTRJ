# Policy
cmake_minimum_required( VERSION 3.0) 


# Project name
project(xdr C)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(SOURCE FComms.cpp  Fstream.cpp  FstreamF.cpp  Metric.cpp  NewMPI.cpp  Timer.cpp
	Finalize.cpp  FstreamC.cpp  HeaderTrj.cpp  MyUtilClass.cpp  Split.cpp  histograms.cpp
	opgather.cpp BSpline.cpp BSpmod.cpp VecRotate.cpp DiffCoeffs.cpp Spline1DInterpolant.cpp)
add_library(traj STATIC ${SOURCE})
include_directories("./")
